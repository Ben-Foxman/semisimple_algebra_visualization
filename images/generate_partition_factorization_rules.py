"""Generate a JPG illustrating three partition factorization rules."""

from __future__ import annotations

import subprocess
import sys
import tempfile
from itertools import permutations
from pathlib import Path


ROW_TOPS = [27.2, 15.2, 3.2]
LEFT_SHIFT = 1.8
RIGHT_SHIFT = 24.4
ARROW_LEFT_X = 15.60
ARROW_RIGHT_X = 21.60
ARROW_TEXT_X = 18.60
COL_STEP = 2.75
TOP_ROW_OFFSET = 3.2
BOTTOM_ROW_OFFSET = 0.0
TOP_LABEL_OFFSET = 4.15
BOTTOM_LABEL_OFFSET = -1.08
ROW_TITLE_OFFSET = 7.0
ROW_CAPTION_OFFSET = 5.95
BRACE_OFFSET = 5.30
ARROW_Y_OFFSET = 1.65
PERMUTATION_Y_OFFSET = 3.35


def x_positions(count: int) -> list[float]:
    return [idx * COL_STEP for idx in range(count)]


def top_label(col: int) -> str:
    return rf"$\displaystyle {col}$"


def bottom_label(col: int) -> str:
    return rf"$\displaystyle {col}'$"


def node_row(node: tuple[str, int]) -> str:
    return node[0]


def node_sort_key(node: tuple[str, int]) -> tuple[int, int]:
    row, idx = node
    return (0 if row == "top" else 1, idx)


def same_row_adjacent(node_a: tuple[str, int], node_b: tuple[str, int]) -> bool:
    row_a, idx_a = node_a
    row_b, idx_b = node_b
    return row_a == row_b and abs(idx_a - idx_b) == 1


def edge_priority(node_a: tuple[str, int], node_b: tuple[str, int]) -> int:
    if node_a[0] != node_b[0]:
        return 2
    if same_row_adjacent(node_a, node_b):
        return 0
    return 1


def edge_draw_color(color: str) -> str:
    return color


def curve_offsets(count: int) -> list[float]:
    if count <= 1:
        return [0.0] * count
    max_offset = 0.18 + 0.04 * max(0, count - 2)
    if count == 2:
        return [-max_offset, max_offset]
    step = (2 * max_offset) / (count - 1)
    return [-max_offset + i * step for i in range(count)]


def node_coords(
    x_shift: float,
    y_shift: float,
    node: tuple[str, int],
    xs: list[float],
) -> tuple[float, float]:
    row, idx = node
    y = y_shift + (TOP_ROW_OFFSET if row == "top" else BOTTOM_ROW_OFFSET)
    return x_shift + xs[idx], y


def component_cycle_edges(component: list[tuple[str, int]]) -> list[tuple[tuple[str, int], tuple[str, int]]]:
    if len(component) < 2:
        return []
    ordered = sorted(component, key=node_sort_key)
    if len(component) == 2:
        return [(ordered[0], ordered[1])]

    start = ordered[0]
    best_score: tuple[int, int, int, tuple[float, ...]] | None = None
    best_cycle: list[tuple[tuple[str, int], tuple[str, int]]] = []

    for perm in permutations(ordered[1:]):
        cycle = [start, *perm]
        edges = [(cycle[i], cycle[(i + 1) % len(cycle)]) for i in range(len(cycle))]
        counts = (
            sum(1 for a, b in edges if edge_priority(a, b) == 0),
            sum(1 for a, b in edges if edge_priority(a, b) == 1),
            sum(1 for a, b in edges if edge_priority(a, b) == 2),
        )
        lengths = tuple(
            sorted(
                (
                    (node_sort_key(a), node_sort_key(b))
                    if node_sort_key(a) <= node_sort_key(b)
                    else (node_sort_key(b), node_sort_key(a))
                )
                for a, b in edges
            )
        )
        score = (counts[0], counts[1], -counts[2], lengths)
        if best_score is None or score > best_score:
            best_score = score
            best_cycle = edges

    return best_cycle


def component_lines(
    x_shift: float,
    y_shift: float,
    top_colors: list[str],
    bottom_colors: list[str],
) -> list[str]:
    xs = x_positions(len(top_colors))
    components: dict[str, list[tuple[str, int]]] = {}

    for idx, color in enumerate(top_colors):
        components.setdefault(color, []).append(("top", idx))
    for idx, color in enumerate(bottom_colors):
        components.setdefault(color, []).append(("bottom", idx))

    lines: list[str] = []
    for color, component in components.items():
        edges = component_cycle_edges(component)
        curved_edges = [
            (node_a, node_b)
            for node_a, node_b in edges
            if node_row(node_a) == node_row(node_b) and edge_priority(node_a, node_b) != 0
        ]
        curved_offsets = curve_offsets(len(curved_edges))
        curve_idx = 0
        style = rf"componentedge, draw={edge_draw_color(color)}"
        for node_a, node_b in edges:
            x_a, y_a = node_coords(x_shift, y_shift, node_a, xs)
            x_b, y_b = node_coords(x_shift, y_shift, node_b, xs)
            if node_a[0] != node_b[0]:
                lines.append(rf"\draw[{style}] ({x_a:.2f}, {y_a:.2f}) -- ({x_b:.2f}, {y_b:.2f});")
            elif same_row_adjacent(node_a, node_b):
                lines.append(rf"\draw[{style}] ({x_a:.2f}, {y_a:.2f}) -- ({x_b:.2f}, {y_b:.2f});")
            else:
                center_y = y_shift + (TOP_ROW_OFFSET + BOTTOM_ROW_OFFSET) / 2 + curved_offsets[curve_idx]
                curve_idx += 1
                lines.append(
                    rf"\draw[{style}] ({x_a:.2f}, {y_a:.2f}) .. controls ({x_a:.2f}, {center_y:.2f}) "
                    rf"and ({x_b:.2f}, {center_y:.2f}) .. ({x_b:.2f}, {y_b:.2f});"
                )

    return lines


def draw_diagram(
    x_shift: float,
    y_shift: float,
    top_colors: list[str],
    bottom_colors: list[str],
    caption: str | None = None,
) -> list[str]:
    lines: list[str] = []
    xs = x_positions(len(top_colors))
    lines.extend(component_lines(x_shift, y_shift, top_colors, bottom_colors))

    for idx, (x, color) in enumerate(zip(xs, top_colors), start=1):
        lines.append(
            rf"\node[vertex, fill={color}] at ({x_shift + x:.2f}, {y_shift + TOP_ROW_OFFSET:.2f}) {{}};"
        )
        lines.append(
            rf"\node[labelfont, anchor=south] at ({x_shift + x:.2f}, {y_shift + TOP_LABEL_OFFSET:.2f}) "
            rf"{{{top_label(idx)}}};"
        )

    for idx, (x, color) in enumerate(zip(xs, bottom_colors), start=1):
        lines.append(
            rf"\node[vertex, fill={color}] at ({x_shift + x:.2f}, {y_shift + BOTTOM_ROW_OFFSET:.2f}) {{}};"
        )
        lines.append(
            rf"\node[labelfont, anchor=north] at ({x_shift + x:.2f}, {y_shift + BOTTOM_LABEL_OFFSET:.2f}) "
            rf"{{{bottom_label(idx)}}};"
        )

    if caption is not None:
        center_x = x_shift + xs[-1] / 2
        lines.append(
            rf"\node[captionfont] at ({center_x:.2f}, {y_shift + ROW_CAPTION_OFFSET:.2f}) "
            rf"{{\scalebox{{1.18}}{{{caption}}}}};"
        )

    return lines


def brace_lines(x_shift: float, y_shift: float, first_col: int, last_col: int, label: str) -> list[str]:
    xs = x_positions(last_col)
    left_x = x_shift + xs[first_col - 1] - 0.95
    right_x = x_shift + xs[last_col - 1] + 0.95
    y = y_shift + BRACE_OFFSET
    center_x = (left_x + right_x) / 2
    return [
        rf"\draw[brace] ({left_x:.2f}, {y:.2f}) -- ({right_x:.2f}, {y:.2f});",
        rf"\node[bracefont, anchor=south] at ({center_x:.2f}, {y + 0.40:.2f}) "
        rf"{{\scalebox{{1.18}}{{{label}}}}};",
    ]


def row_annotation(row_idx: int, title: str, arrow_formula: str) -> list[str]:
    y_shift = ROW_TOPS[row_idx]
    return [
        rf"\node[rowtitlefont] at ({ARROW_TEXT_X:.2f}, {y_shift + ROW_TITLE_OFFSET:.2f}) {{{title}}};",
        rf"\draw[rulearrow] ({ARROW_LEFT_X:.2f}, {y_shift + ARROW_Y_OFFSET:.2f}) -- "
        rf"({ARROW_RIGHT_X:.2f}, {y_shift + ARROW_Y_OFFSET:.2f});",
        rf"\node[permutationfont] at ({ARROW_TEXT_X:.2f}, {y_shift + PERMUTATION_Y_OFFSET:.2f}) "
        rf"{{\scalebox{{1.24}}{{{arrow_formula}}}}};",
    ]


def row_one_lines() -> list[str]:
    lines: list[str] = []
    y_shift = ROW_TOPS[0]

    left_top = ["blocka", "blockd", "blockb", "blockc", "blockb"]
    left_bottom = ["blocka", "blockb", "blockd", "blockc", "blockc"]
    right_top = ["blocka", "blockb", "blockb", "blockc", "blockd"]
    right_bottom = ["blocka", "blockb", "blockc", "blockc", "blockd"]

    lines.extend(
        row_annotation(
            0,
            r"{\scalebox{1.28}{$\textbf{Partition 1}$}}",
            r"$\begin{aligned}\pi_1&=(2\,5)\\[2pt]\pi_2&=(3'\,5')\end{aligned}$",
        )
    )
    lines.extend(draw_diagram(LEFT_SHIFT, y_shift, left_top, left_bottom, r"$D_a\in P_{5}(d)$"))
    lines.extend(draw_diagram(RIGHT_SHIFT, y_shift, right_top, right_bottom))
    lines.extend(
        brace_lines(
            RIGHT_SHIFT,
            y_shift,
            1,
            4,
            r"$D_b\in P_{4}(d)$",
        )
    )
    return lines


def row_two_lines() -> list[str]:
    lines: list[str] = []
    y_shift = ROW_TOPS[1]

    left_top = ["blockc", "blocka", "blockc", "blocka", "blockb"]
    left_bottom = ["blocka", "blockc", "blockb", "blockb", "blockb"]
    right_top = ["blocka", "blocka", "blockb", "blockc", "blockc"]
    right_bottom = ["blocka", "blockb", "blockb", "blockb", "blockc"]

    lines.extend(
        row_annotation(
            1,
            r"{\scalebox{1.28}{$\textbf{Partition 2}$}}",
            r"$\begin{aligned}\pi_1&=(1\,4)(3\,5)\\[2pt]\pi_2&=(2'\,5')\end{aligned}$",
        )
    )
    lines.extend(draw_diagram(LEFT_SHIFT, y_shift, left_top, left_bottom, r"$D_a\in P_{5}(d)$"))
    lines.extend(draw_diagram(RIGHT_SHIFT, y_shift, right_top, right_bottom))
    lines.extend(
        brace_lines(
            RIGHT_SHIFT,
            y_shift,
            1,
            4,
            r"$b_{4}D_b,\; D_b\in P_{4}(d)$",
        )
    )
    return lines


def row_three_lines() -> list[str]:
    lines: list[str] = []
    y_shift = ROW_TOPS[2]

    left_top = ["blocka", "blockb", "blockc", "blockc", "blockc"]
    left_bottom = ["blockb", "blocka", "blockb", "blocka", "blockd"]
    right_top = ["blocka", "blockc", "blockc", "blockc", "blockb"]
    right_bottom = ["blocka", "blocka", "blockd", "blockb", "blockb"]

    lines.extend(
        row_annotation(
            2,
            r"{\scalebox{1.28}{$\textbf{Partition 3}$}}",
            r"$\begin{aligned}\pi_1&=(2\,5)\\[2pt]\pi_2&=(1'\,4')(3'\,5')\end{aligned}$",
        )
    )
    lines.extend(draw_diagram(LEFT_SHIFT, y_shift, left_top, left_bottom, r"$D_a\in P_{5}(d)$"))
    lines.extend(draw_diagram(RIGHT_SHIFT, y_shift, right_top, right_bottom))
    lines.extend(
        brace_lines(
            RIGHT_SHIFT,
            y_shift,
            1,
            4,
            r"$D_bb_{4},\; D_b\in P_{4}(d)$",
        )
    )
    return lines


def tikz_source() -> str:
    lines: list[str] = []
    lines.extend(row_one_lines())
    lines.extend(row_two_lines())
    lines.extend(row_three_lines())
    body = "\n".join(lines)

    return rf"""\documentclass{{article}}
\usepackage[paperwidth=43cm,paperheight=39cm,margin=0pt]{{geometry}}
\usepackage{{amsmath}}
\usepackage{{graphicx}}
\usepackage{{tikz}}
\usetikzlibrary{{arrows.meta,decorations.pathreplacing}}
\definecolor{{blocka}}{{RGB}}{{88, 165, 242}}
\definecolor{{blockb}}{{RGB}}{{255, 155, 88}}
\definecolor{{blockc}}{{RGB}}{{122, 212, 147}}
\definecolor{{blockd}}{{RGB}}{{247, 206, 63}}
\definecolor{{blocke}}{{RGB}}{{186, 147, 245}}
\definecolor{{blockf}}{{RGB}}{{255, 98, 57}}
\pagestyle{{empty}}

\begin{{document}}
\thispagestyle{{empty}}
\vspace*{{\fill}}
\begin{{center}}
\begin{{tikzpicture}}[
    componentedge/.style={{
        line width=5.0pt,
        line cap=round
    }},
    vertex/.style={{
        circle,
        draw=black,
        line width=1.3pt,
        minimum size=1.56cm,
        inner sep=0pt
    }},
    labelfont/.style={{font=\fontsize{{25}}{{25}}\selectfont}},
    captionfont/.style={{font=\fontsize{{28}}{{28}}\selectfont}},
    rowtitlefont/.style={{font=\fontsize{{25}}{{25}}\selectfont}},
    permutationfont/.style={{font=\fontsize{{22}}{{22}}\selectfont}},
    rulearrow/.style={{->, >=Latex, line width=1.2pt}},
    bracefont/.style={{font=\fontsize{{23}}{{23}}\selectfont}},
    brace/.style={{decorate, decoration={{brace, amplitude=8pt}}, line width=1.25pt}}
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
        else Path(__file__).with_name("partition_factorization_rules.jpg")
    )
    output.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as temp_dir:
        workdir = Path(temp_dir)
        tex_path = workdir / "partition_factorization_rules.tex"
        tex_path.write_text(tikz_source(), encoding="utf-8")
        pdf_path = build_pdf(workdir, tex_path)
        convert_pdf_to_jpg(pdf_path, output)

    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
