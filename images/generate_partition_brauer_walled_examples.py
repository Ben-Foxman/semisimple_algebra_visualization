"""Generate a JPG with representative partition, Brauer, and walled Brauer diagrams."""

from __future__ import annotations

import subprocess
import sys
import tempfile
from itertools import permutations
from pathlib import Path


TOP_Y = 3.2
BOTTOM_Y = 0.0
TOP_LABEL_Y = 4.15
BOTTOM_LABEL_Y = -1.12
CAPTION_Y = -4.75
COL_STEP = 4.0
WALL_GAP = 2.20
PANEL_SHIFTS = [1.5, 19.0, 37.0]

PANELS = [
    {
        "caption": r"$D\in P_{4}(d)$",
        "wall_after": None,
        "top_colors": ["blocka", "blockb", "blockb", "blockc"],
        "bottom_colors": ["blocka", "blockb", "blockc", "blockc"],
        "components": {
            "blocka": [("top", 0), ("bottom", 0)],
            "blockb": [("top", 1), ("top", 2), ("bottom", 1)],
            "blockc": [("top", 3), ("bottom", 2), ("bottom", 3)],
        },
    },
    {
        "caption": r"$D\in B_{4}(d)$",
        "wall_after": None,
        "top_colors": ["blocka", "blocka", "blockb", "blockc"],
        "bottom_colors": ["blockd", "blockd", "blockc", "blockb"],
        "components": {
            "blocka": [("top", 0), ("top", 1)],
            "blockb": [("top", 2), ("bottom", 3)],
            "blockc": [("top", 3), ("bottom", 2)],
            "blockd": [("bottom", 0), ("bottom", 1)],
        },
    },
    {
        "caption": r"$D\in B_{2,2}(d)$",
        "wall_after": 2,
        "top_colors": ["blocka", "blockb", "blockb", "blockd"],
        "bottom_colors": ["blockc", "blocka", "blockc", "blockd"],
        "components": {
            "blocka": [("top", 0), ("bottom", 1)],
            "blockb": [("top", 1), ("top", 2)],
            "blockc": [("bottom", 0), ("bottom", 2)],
            "blockd": [("top", 3), ("bottom", 3)],
        },
    },
]


def x_positions(count: int, wall_after: int | None) -> list[float]:
    xs: list[float] = []
    for idx in range(count):
        x = idx * COL_STEP
        if wall_after is not None and idx >= wall_after:
            x += WALL_GAP
        xs.append(x)
    return xs


def top_label(col: int) -> str:
    return rf"$\displaystyle {col}$"


def bottom_label(col: int) -> str:
    return rf"$\displaystyle {col}'$"


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


def node_coords(
    x_shift: float,
    node: tuple[str, int],
    xs: list[float],
) -> tuple[float, float]:
    row, idx = node
    y = TOP_Y if row == "top" else BOTTOM_Y
    return x_shift + xs[idx], y


def component_cycle_edges(
    component: list[tuple[str, int]],
) -> list[tuple[tuple[str, int], tuple[str, int]]]:
    if len(component) < 2:
        return []
    ordered = sorted(component, key=node_sort_key)
    if len(component) == 2:
        return [(ordered[0], ordered[1])]

    start = ordered[0]
    best_score: tuple[int, int, int, tuple[tuple[tuple[int, int], tuple[int, int]], ...]] | None = None
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
    panel: dict[str, object],
) -> list[str]:
    xs = x_positions(len(panel["top_colors"]), panel["wall_after"])
    lines: list[str] = []
    curved_offset = 0
    for color, component in panel["components"].items():
        for node_a, node_b in component_cycle_edges(component):
            x_a, y_a = node_coords(x_shift, node_a, xs)
            x_b, y_b = node_coords(x_shift, node_b, xs)
            if node_a[0] != node_b[0]:
                lines.append(
                    rf"\draw[componentedge, draw={color}] ({x_a:.2f}, {y_a:.2f}) -- ({x_b:.2f}, {y_b:.2f});"
                )
            elif same_row_adjacent(node_a, node_b):
                lines.append(
                    rf"\draw[componentedge, draw={color}] ({x_a:.2f}, {y_a:.2f}) -- ({x_b:.2f}, {y_b:.2f});"
                )
            else:
                center_y = (TOP_Y + BOTTOM_Y) / 2
                bend = [0.0, 0.22, -0.22, 0.38, -0.38][curved_offset % 5]
                curved_offset += 1
                lines.append(
                    rf"\draw[componentedge, draw={color}] ({x_a:.2f}, {y_a:.2f}) .. controls ({x_a:.2f}, {center_y + bend:.2f}) "
                    rf"and ({x_b:.2f}, {center_y + bend:.2f}) .. ({x_b:.2f}, {y_b:.2f});"
                )
    return lines


def wall_x(x_shift: float, count: int, wall_after: int) -> float:
    xs = x_positions(count, wall_after)
    return x_shift + (xs[wall_after - 1] + xs[wall_after]) / 2


def panel_lines(x_shift: float, panel: dict[str, object]) -> list[str]:
    xs = x_positions(len(panel["top_colors"]), panel["wall_after"])
    lines = component_lines(x_shift, panel)

    if panel["wall_after"] is not None:
        lines.append(
            rf"\draw[walldivider] ({wall_x(x_shift, len(panel['top_colors']), panel['wall_after']):.2f}, -1.60) -- "
            rf"({wall_x(x_shift, len(panel['top_colors']), panel['wall_after']):.2f}, 4.90);"
        )

    for idx, (x, color) in enumerate(zip(xs, panel["top_colors"]), start=1):
        lines.append(rf"\node[vertexmask] at ({x_shift + x:.2f}, {TOP_Y:.2f}) {{}};")
        lines.append(rf"\node[vertex, fill={color}] at ({x_shift + x:.2f}, {TOP_Y:.2f}) {{}};")
        lines.append(
            rf"\node[labelfont, anchor=south] at ({x_shift + x:.2f}, {TOP_LABEL_Y:.2f}) {{{top_label(idx)}}};"
        )

    for idx, (x, color) in enumerate(zip(xs, panel["bottom_colors"]), start=1):
        lines.append(rf"\node[vertexmask] at ({x_shift + x:.2f}, {BOTTOM_Y:.2f}) {{}};")
        lines.append(rf"\node[vertex, fill={color}] at ({x_shift + x:.2f}, {BOTTOM_Y:.2f}) {{}};")
        lines.append(
            rf"\node[labelfont, anchor=north] at ({x_shift + x:.2f}, {BOTTOM_LABEL_Y:.2f}) {{{bottom_label(idx)}}};"
        )

    center_x = x_shift + (xs[0] + xs[-1]) / 2
    lines.append(
        rf"\node[captionfont] at ({center_x:.2f}, {CAPTION_Y:.2f}) {{\scalebox{{1.18}}{{{panel['caption']}}}}};"
    )
    return lines


def tikz_source() -> str:
    lines: list[str] = []
    for x_shift, panel in zip(PANEL_SHIFTS, PANELS):
        lines.extend(panel_lines(x_shift, panel))
    nodes = "\n".join(lines)
    return rf"""\documentclass{{article}}
\usepackage[paperwidth=54cm,paperheight=11.2cm,margin=0pt]{{geometry}}
\usepackage{{amsmath}}
\usepackage{{graphicx}}
\usepackage{{tikz}}
\definecolor{{blocka}}{{RGB}}{{59, 129, 239}}
\definecolor{{blockb}}{{RGB}}{{255, 128, 58}}
\definecolor{{blockc}}{{RGB}}{{123, 214, 146}}
\definecolor{{blockd}}{{RGB}}{{255, 192, 214}}
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
    walldivider/.style={{black!65, dashed, line width=0.9pt}},
    vertexmask/.style={{
        circle,
        draw=none,
        fill=white,
        minimum size=1.70cm,
        inner sep=0pt
    }},
    vertex/.style={{
        circle,
        draw=black,
        line width=1.3pt,
        minimum size=1.70cm,
        inner sep=0pt
    }},
    labelfont/.style={{
        font=\fontsize{{26}}{{26}}\selectfont
    }},
    captionfont/.style={{
        font=\fontsize{{32}}{{34}}\selectfont
    }}
]
{nodes}
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
        [
            "pdflatex",
            "-interaction=nonstopmode",
            "-halt-on-error",
            tex_path.name,
        ],
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
        else Path(__file__).with_name("partition_brauer_walled_examples.jpg")
    )
    output.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as temp_dir:
        workdir = Path(temp_dir)
        tex_path = workdir / "partition_brauer_walled_examples.tex"
        tex_path.write_text(tikz_source(), encoding="utf-8")
        pdf_path = build_pdf(workdir, tex_path)
        convert_pdf_to_jpg(pdf_path, output)

    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
