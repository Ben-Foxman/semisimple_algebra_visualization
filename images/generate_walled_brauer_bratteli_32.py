"""Generate a JPG for the Bratteli diagram of the walled Brauer algebra B_{3,2}(d)."""

from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path


Partition = tuple[int, ...]
BiPartition = tuple[Partition, Partition]


CHAIN = [(0, 0), (1, 0), (1, 1), (2, 1), (2, 2), (3, 2)]
COL_X = [0.0, 8.0, 16.0, 24.0, 32.0, 40.0]
HEADER_Y = 18.2
ROW_GAP = 3.3

PART_BOX = 0.96
PAIR_H = 2.5
EMPTY_W = 0.60
EMPTY_H = 0.95
COMMA_GAP_LEFT = 0.12
COMMA_GAP_RIGHT = 0.26
OUTER_GAP = 0.10
COMMA_HALF = 0.10
PAREN_W = 0.34

LEVELS: list[list[BiPartition]] = [
    [(tuple(), tuple())],
    [((1,), tuple())],
    [((1,), (1,)), (tuple(), tuple())],
    [((2,), (1,)), ((1, 1), (1,)), ((1,), tuple())],
    [
        ((2,), (2,)),
        ((2,), (1, 1)),
        ((1, 1), (2,)),
        ((1, 1), (1, 1)),
        ((1,), (1,)),
        (tuple(), tuple()),
    ],
    [
        ((3,), (2,)),
        ((3,), (1, 1)),
        ((2, 1), (2,)),
        ((2, 1), (1, 1)),
        ((1, 1, 1), (2,)),
        ((1, 1, 1), (1, 1)),
        ((2,), (1,)),
        ((1, 1), (1,)),
        ((1,), tuple()),
    ],
]


def removable_partitions(part: Partition) -> list[Partition]:
    out: list[Partition] = []
    for i, row in enumerate(part):
        next_row = part[i + 1] if i + 1 < len(part) else 0
        if row > next_row:
            new = list(part)
            new[i] -= 1
            while new and new[-1] == 0:
                new.pop()
            out.append(tuple(new))
    return out


def addable_partitions(part: Partition) -> list[Partition]:
    if not part:
        return [(1,)]
    out: list[Partition] = []
    for i, row in enumerate(part):
        prev_row = part[i - 1] if i > 0 else None
        if i == 0 or row < prev_row:
            new = list(part)
            new[i] += 1
            out.append(tuple(new))
    out.append(tuple(list(part) + [1]))
    dedup: list[Partition] = []
    for cand in out:
        if cand not in dedup:
            dedup.append(cand)
    return dedup


def reachable(src: BiPartition, dst: BiPartition, col_idx: int) -> bool:
    src_left, src_right = src
    dst_left, dst_right = dst
    step = CHAIN[col_idx + 1]

    if step[0] > CHAIN[col_idx][0]:
        left_add = dst_right == src_right and dst_left in addable_partitions(src_left)
        right_rem = dst_left == src_left and dst_right in removable_partitions(src_right)
        return left_add or right_rem

    right_add = dst_left == src_left and dst_right in addable_partitions(src_right)
    left_rem = dst_right == src_right and dst_left in removable_partitions(src_left)
    return right_add or left_rem


def path_counts() -> list[dict[BiPartition, int]]:
    counts: list[dict[BiPartition, int]] = [{LEVELS[0][0]: 1}]
    for col_idx in range(len(LEVELS) - 1):
        current = counts[-1]
        nxt: dict[BiPartition, int] = {node: 0 for node in LEVELS[col_idx + 1]}
        for src in LEVELS[col_idx]:
            for dst in LEVELS[col_idx + 1]:
                if reachable(src, dst, col_idx):
                    nxt[dst] += current[src]
        counts.append(nxt)
    return counts


def y_positions(count: int) -> list[float]:
    top = (count - 1) * ROW_GAP / 2
    return [top - idx * ROW_GAP for idx in range(count)]


def draw_partition(part: Partition, center_x: float, center_y: float, lines: list[str]) -> None:
    if not part:
        lines.append(
            rf"\node[emptynode] at ({center_x:.2f}, {center_y:.2f}) "
            rf"{{\scalebox{{1.75}}{{$\emptyset$}}}};"
        )
        return

    total_h = len(part) * PART_BOX
    max_w = part[0] * PART_BOX
    top_y = center_y + total_h / 2
    left_x = center_x - max_w / 2

    for row_idx, row_len in enumerate(part):
        y_top = top_y - row_idx * PART_BOX
        for col_idx in range(row_len):
            x_left = left_x + col_idx * PART_BOX
            lines.append(
                rf"\filldraw[fill=white, draw=black, line width=0.85pt] "
                rf"({x_left:.2f}, {y_top:.2f}) rectangle "
                rf"({x_left + PART_BOX:.2f}, {y_top - PART_BOX:.2f});"
            )


def partition_dimensions(part: Partition) -> tuple[float, float]:
    if not part:
        return EMPTY_W, EMPTY_H
    return part[0] * PART_BOX, len(part) * PART_BOX


def pair_half_width(node: BiPartition) -> float:
    _, _, _, _, _, left_paren_x, right_paren_x, _, _ = pair_geometry(0.0, node)
    return max(-left_paren_x, right_paren_x) + 0.04


def pair_geometry(
    center_x: float, node: BiPartition
) -> tuple[float, float, float, float, float, float, float, float, float]:
    left, right = node
    left_w, left_h = partition_dimensions(left)
    right_w, right_h = partition_dimensions(right)
    total_w = (
        PAREN_W
        + OUTER_GAP
        + left_w
        + COMMA_GAP_LEFT
        + 2 * COMMA_HALF
        + COMMA_GAP_RIGHT
        + right_w
        + OUTER_GAP
        + PAREN_W
    )
    left_edge_x = center_x - total_w / 2
    left_paren_x = left_edge_x
    left_center_x = left_edge_x + PAREN_W + OUTER_GAP + left_w / 2
    comma_x = left_edge_x + PAREN_W + OUTER_GAP + left_w + COMMA_GAP_LEFT + COMMA_HALF
    right_center_x = comma_x + COMMA_HALF + COMMA_GAP_RIGHT + right_w / 2
    right_paren_x = left_edge_x + total_w
    return (
        left_w,
        right_w,
        left_h,
        right_h,
        left_center_x,
        left_paren_x,
        right_paren_x,
        comma_x,
        right_center_x,
    )


def parenthesis_lines(outer_x: float, center_y: float, content_h: float, side: str) -> list[str]:
    paren_h = max(content_h + 0.38, 1.20)
    top_y = center_y + paren_h / 2
    bot_y = center_y - paren_h / 2
    ctrl_y = center_y + paren_h * 0.20

    if side == "left":
        top_x = outer_x + PAREN_W
        mid_x = outer_x
        return [
            rf"\draw[paren] ({top_x:.2f}, {top_y:.2f}) "
            rf".. controls ({outer_x:.2f}, {top_y:.2f}) and ({outer_x:.2f}, {ctrl_y:.2f}) "
            rf".. ({mid_x:.2f}, {center_y:.2f}) "
            rf".. controls ({outer_x:.2f}, {-ctrl_y + 2*center_y:.2f}) and ({outer_x:.2f}, {bot_y:.2f}) "
            rf".. ({top_x:.2f}, {bot_y:.2f});"
        ]

    top_x = outer_x - PAREN_W
    mid_x = outer_x
    return [
        rf"\draw[paren] ({top_x:.2f}, {top_y:.2f}) "
        rf".. controls ({outer_x:.2f}, {top_y:.2f}) and ({outer_x:.2f}, {ctrl_y:.2f}) "
        rf".. ({mid_x:.2f}, {center_y:.2f}) "
        rf".. controls ({outer_x:.2f}, {-ctrl_y + 2*center_y:.2f}) and ({outer_x:.2f}, {bot_y:.2f}) "
        rf".. ({top_x:.2f}, {bot_y:.2f});"
    ]


def pair_lines(center_x: float, center_y: float, node: BiPartition) -> list[str]:
    left, right = node
    (
        _,
        _,
        left_h,
        right_h,
        left_center_x,
        left_paren_x,
        right_paren_x,
        comma_x,
        right_center_x,
    ) = pair_geometry(center_x, node)
    lines = [
        rf"\node[paircomma, anchor=mid] at ({comma_x:.2f}, {center_y - 0.06:.2f}) {{$,$}};",
    ]
    lines.extend(parenthesis_lines(left_paren_x, center_y, left_h, "left"))
    lines.extend(parenthesis_lines(right_paren_x, center_y, right_h, "right"))
    draw_partition(left, left_center_x, center_y, lines)
    draw_partition(right, right_center_x, center_y, lines)
    return lines


def edge_lines() -> list[str]:
    lines: list[str] = []
    for col_idx in range(len(LEVELS) - 1):
        left_x = COL_X[col_idx]
        right_x = COL_X[col_idx + 1]
        for src_y, src in zip(y_positions(len(LEVELS[col_idx])), LEVELS[col_idx]):
            for dst_y, dst in zip(y_positions(len(LEVELS[col_idx + 1])), LEVELS[col_idx + 1]):
                if not reachable(src, dst, col_idx):
                    continue
                _, _, _, _, _, _, src_right_paren_x, _, _ = pair_geometry(left_x, src)
                _, _, _, _, _, dst_left_paren_x, _, _, _ = pair_geometry(right_x, dst)
                lines.append(
                    rf"\draw[edge] ({src_right_paren_x:.2f}, {src_y:.2f}) -- "
                    rf"({dst_left_paren_x:.2f}, {dst_y:.2f});"
                )
    return lines


def separator_lines() -> list[str]:
    lines: list[str] = []
    for left_x, right_x in zip(COL_X, COL_X[1:]):
        sep_x = (left_x + right_x) / 2
        lines.append(
            rf"\draw[separator] ({sep_x:.2f}, -16.8) -- ({sep_x:.2f}, 18.9);"
        )
    return lines


def header_lines() -> list[str]:
    lines: list[str] = []
    for x, (r, s) in zip(COL_X, CHAIN):
        lines.append(
            rf"\node[colheader] at ({x:.2f}, {HEADER_Y:.2f}) "
            rf"{{\scalebox{{1.55}}{{$\displaystyle B_{{{r},{s}}}(d)$}}}};"
        )
    return lines


def final_dim_lines(counts: dict[BiPartition, int]) -> list[str]:
    lines: list[str] = []
    x = COL_X[-1]
    max_half = max(pair_half_width(node) for node in LEVELS[-1])
    dim_x = x + max_half + 0.90
    for y, node in zip(y_positions(len(LEVELS[-1])), LEVELS[-1]):
        lines.append(
            rf"\node[dimfont, anchor=west] at ({dim_x:.2f}, {y:.2f}) "
            rf"{{\scalebox{{1.25}}{{$\dim = {counts[node]}$}}}};"
        )
    return lines


def tikz_source() -> str:
    counts = path_counts()
    lines: list[str] = []
    lines.extend(separator_lines())
    lines.extend(edge_lines())
    for x, nodes in zip(COL_X, LEVELS):
        for y, node in zip(y_positions(len(nodes)), nodes):
            lines.extend(pair_lines(x, y, node))
    lines.extend(header_lines())
    lines.extend(final_dim_lines(counts[-1]))

    body = "\n".join(lines)
    return rf"""\documentclass{{article}}
\usepackage[paperwidth=52cm,paperheight=38cm,margin=0pt]{{geometry}}
\usepackage{{amsmath}}
\usepackage{{graphicx}}
\usepackage{{tikz}}
\pagestyle{{empty}}

\begin{{document}}
\thispagestyle{{empty}}
\vspace*{{\fill}}
\begin{{center}}
\begin{{tikzpicture}}[
    edge/.style={{black!58, line width=0.9pt}},
    separator/.style={{black!40, dashed, line width=0.8pt}},
    colheader/.style={{font=\fontsize{{42}}{{42}}\selectfont}},
    paren/.style={{black, line width=0.85pt, line cap=round}},
    paircomma/.style={{font=\fontsize{{34}}{{34}}\selectfont}},
    emptynode/.style={{font=\fontsize{{28}}{{28}}\selectfont}},
    dimfont/.style={{font=\fontsize{{24}}{{24}}\selectfont}}
]
{body}
\end{{tikzpicture}}
\end{{center}}
\vspace*{{\fill}}
\end{{document}}
"""


def run_command(command: list[str], cwd: Path) -> None:
    subprocess.run(
        command,
        cwd=cwd,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )


def build_pdf(workdir: Path, tex_path: Path) -> Path:
    run_command(
        ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", tex_path.name],
        cwd=workdir,
    )
    return tex_path.with_suffix(".pdf")


def convert_pdf_to_jpg(pdf_path: Path, jpg_path: Path) -> None:
    run_command(
        [
            "gs",
            "-dSAFER",
            "-dBATCH",
            "-dNOPAUSE",
            "-sDEVICE=jpeg",
            "-dJPEGQ=95",
            "-r300",
            f"-sOutputFile={jpg_path}",
            str(pdf_path),
        ],
        cwd=pdf_path.parent,
    )


def main() -> int:
    output = (
        Path(sys.argv[1]).resolve()
        if len(sys.argv) > 1
        else Path(__file__).with_name("walled_brauer_bratteli_32.jpg")
    )
    output.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as temp_dir:
        workdir = Path(temp_dir)
        tex_path = workdir / "walled_brauer_bratteli_32.tex"
        tex_path.write_text(tikz_source(), encoding="utf-8")
        pdf_path = build_pdf(workdir, tex_path)
        convert_pdf_to_jpg(pdf_path, output)

    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
