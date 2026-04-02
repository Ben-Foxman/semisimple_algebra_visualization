"""Generate a JPG showing s_i, b_i, and p_i for the partition algebra."""

from __future__ import annotations

import subprocess
import sys
import tempfile
from itertools import permutations
from pathlib import Path


PANEL_SHIFTS = [5.0, 23.1, 41.2]
TOP_Y = 3.2
BOTTOM_Y = 0.0
TITLE_Y = 6.15
TOP_LABEL_Y = 4.15
BOTTOM_LABEL_Y = -1.12

FIRST_COLOR = "firstcolor"
LAST_COLOR = "lastcolor"
BLOCK_A = "blocka"
BLOCK_B = "blockb"

PANELS = [
    {
        "title": r"$\displaystyle s_i$",
        "center_x": 5.9,
        "visible_x": [0.0, 4.2, 7.6, 11.8],
        "ellipsis_x": [2.1, 9.7],
        "top_labels": [r"$\displaystyle 1$", r"$\displaystyle i$", r"$\displaystyle i+1$", r"$\displaystyle n$"],
        "bottom_labels": [
            r"$\displaystyle 1'$",
            r"$\displaystyle i'$",
            r"$\displaystyle (i+1)'$",
            r"$\displaystyle n'$",
        ],
        "top_colors": [FIRST_COLOR, BLOCK_A, BLOCK_B, LAST_COLOR],
        "bottom_colors": [FIRST_COLOR, BLOCK_B, BLOCK_A, LAST_COLOR],
        "top_styles": ["lightvertex", "vertex", "vertex", "lightvertex"],
        "bottom_styles": ["lightvertex", "vertex", "vertex", "lightvertex"],
        "components": [
            [("top", 0), ("bottom", 0)],
            [("top", 1), ("bottom", 2)],
            [("top", 2), ("bottom", 1)],
            [("top", 3), ("bottom", 3)],
        ],
        "adjacent_pairs": {("top", 1, 2), ("bottom", 1, 2)},
    },
    {
        "title": r"$\displaystyle b_i$",
        "center_x": 5.9,
        "visible_x": [0.0, 4.2, 7.6, 11.8],
        "ellipsis_x": [2.1, 9.7],
        "top_labels": [r"$\displaystyle 1$", r"$\displaystyle i$", r"$\displaystyle i+1$", r"$\displaystyle n$"],
        "bottom_labels": [
            r"$\displaystyle 1'$",
            r"$\displaystyle i'$",
            r"$\displaystyle (i+1)'$",
            r"$\displaystyle n'$",
        ],
        "top_colors": [FIRST_COLOR, BLOCK_A, BLOCK_A, LAST_COLOR],
        "bottom_colors": [FIRST_COLOR, BLOCK_A, BLOCK_A, LAST_COLOR],
        "top_styles": ["lightvertex", "vertex", "vertex", "lightvertex"],
        "bottom_styles": ["lightvertex", "vertex", "vertex", "lightvertex"],
        "components": [
            [("top", 0), ("bottom", 0)],
            [("top", 1), ("top", 2), ("bottom", 1), ("bottom", 2)],
            [("top", 3), ("bottom", 3)],
        ],
        "adjacent_pairs": {("top", 1, 2), ("bottom", 1, 2)},
    },
    {
        "title": r"$\displaystyle p_i$",
        "center_x": 5.9,
        "visible_x": [0.0, 5.2, 11.8],
        "ellipsis_x": [2.5, 9.0],
        "top_labels": [r"$\displaystyle 1$", r"$\displaystyle i$", r"$\displaystyle n$"],
        "bottom_labels": [
            r"$\displaystyle 1'$",
            r"$\displaystyle i'$",
            r"$\displaystyle n'$",
        ],
        "top_colors": [FIRST_COLOR, BLOCK_A, LAST_COLOR],
        "bottom_colors": [FIRST_COLOR, BLOCK_B, LAST_COLOR],
        "top_styles": ["lightvertex", "vertex", "lightvertex"],
        "bottom_styles": ["lightvertex", "vertex", "lightvertex"],
        "components": [
            [("top", 0), ("bottom", 0)],
            [("top", 1)],
            [("bottom", 1)],
            [("top", 2), ("bottom", 2)],
        ],
        "adjacent_pairs": set(),
    },
]


def edge_priority(
    node_a: tuple[str, int], node_b: tuple[str, int], adjacent_pairs: set[tuple[str, int, int]]
) -> int:
    row_a, idx_a = node_a
    row_b, idx_b = node_b
    if row_a != row_b:
        return 2
    pair = (row_a, min(idx_a, idx_b), max(idx_a, idx_b))
    if pair in adjacent_pairs:
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


def edge_length(
    node_a: tuple[str, int],
    node_b: tuple[str, int],
    coords,
) -> float:
    x_a, y_a = coords(*node_a)
    x_b, y_b = coords(*node_b)
    return (x_a - x_b) ** 2 + (y_a - y_b) ** 2


def cycle_edges(
    component: list[tuple[str, int]],
    node_order,
    edge_priority_fn,
    edge_length_fn,
) -> list[tuple[tuple[str, int], tuple[str, int]]]:
    ordered = sorted(component, key=node_order)
    if len(ordered) < 2:
        return []
    if len(ordered) == 2:
        return [(ordered[0], ordered[1])]

    start = ordered[0]
    best_score: tuple[object, ...] | None = None
    best_edges: list[tuple[tuple[str, int], tuple[str, int]]] = []

    for perm in permutations(ordered[1:]):
        cycle = [start, *perm]
        edges = [
            (cycle[idx], cycle[(idx + 1) % len(cycle)])
            for idx in range(len(cycle))
        ]
        priorities = tuple(sorted(edge_priority_fn(a, b) for a, b in edges))
        lengths = tuple(sorted(edge_length_fn(a, b) for a, b in edges))
        score = (priorities, lengths, tuple(cycle))
        if best_score is None or score < best_score:
            best_score = score
            best_edges = edges

    return best_edges


def component_color(component: list[tuple[str, int]], panel: dict[str, object]) -> str:
    row, idx = component[0]
    if row == "top":
        return str(panel["top_colors"][idx])
    return str(panel["bottom_colors"][idx])


def edge_lines(shift: float, panel: dict[str, object]) -> list[str]:
    visible_x = panel["visible_x"]
    components = panel["components"]
    adjacent_pairs = panel["adjacent_pairs"]
    center_y = (TOP_Y + BOTTOM_Y) / 2

    def node_order(node: tuple[str, int]) -> tuple[int, int]:
        row, idx = node
        return (0 if row == "top" else 1, idx)

    def coords(row: str, idx: int) -> tuple[float, float]:
        y = TOP_Y if row == "top" else BOTTOM_Y
        return shift + visible_x[idx], y

    def edge_priority_fn(
        node_a: tuple[str, int], node_b: tuple[str, int]
    ) -> int:
        return edge_priority(node_a, node_b, adjacent_pairs)

    def edge_length_fn(
        node_a: tuple[str, int], node_b: tuple[str, int]
    ) -> float:
        return edge_length(node_a, node_b, coords)

    lines: list[str] = []
    for component in components:
        edges = cycle_edges(component, node_order, edge_priority_fn, edge_length_fn)
        curved_edges = [
            (node_a, node_b)
            for node_a, node_b in edges
            if node_a[0] == node_b[0] and edge_priority_fn(node_a, node_b) != 0
        ]
        curved_offsets = curve_offsets(len(curved_edges))
        curve_idx = 0
        color = component_color(component, panel)
        style = rf"componentedge, draw={edge_draw_color(color)}"
        for node_a, node_b in edges:
            x_a, y_a = coords(*node_a)
            x_b, y_b = coords(*node_b)

            if node_a[0] != node_b[0]:
                lines.append(rf"\draw[{style}] ({x_a:.2f}, {y_a:.2f}) -- ({x_b:.2f}, {y_b:.2f});")
                continue

            if edge_priority_fn(node_a, node_b) == 0:
                lines.append(rf"\draw[{style}] ({x_a:.2f}, {y_a:.2f}) -- ({x_b:.2f}, {y_b:.2f});")
            else:
                bend_y = center_y + curved_offsets[curve_idx]
                curve_idx += 1
                lines.append(
                    rf"\draw[{style}] ({x_a:.2f}, {y_a:.2f}) .. controls ({x_a:.2f}, {bend_y:.2f}) "
                    rf"and ({x_b:.2f}, {bend_y:.2f}) .. ({x_b:.2f}, {y_b:.2f});"
                )

    return lines


def vertex_lines(
    x: float,
    y: float,
    label_y: float,
    label: str,
    color: str,
    style: str,
    anchor: str,
) -> list[str]:
    return [
        rf"\node[vertexmask] at ({x:.2f}, {y:.2f}) {{}};",
        rf"\node[{style}, fill={color}] at ({x:.2f}, {y:.2f}) {{}};",
        rf"\node[labelfont, anchor={anchor}] at ({x:.2f}, {label_y:.2f}) {{{label}}};",
    ]


def panel_lines(shift: float, panel: dict[str, object]) -> list[str]:
    title = panel["title"]
    center_x = panel["center_x"]
    visible_x = panel["visible_x"]
    ellipsis_x = panel["ellipsis_x"]
    top_labels = panel["top_labels"]
    bottom_labels = panel["bottom_labels"]
    top_colors = panel["top_colors"]
    bottom_colors = panel["bottom_colors"]
    top_styles = panel["top_styles"]
    bottom_styles = panel["bottom_styles"]

    lines = edge_lines(shift, panel)
    lines.extend(
        [
        rf"\node[titlefont] at ({shift + center_x:.2f}, {TITLE_Y:.2f}) {{\scalebox{{2.00}}{{{title}}}}};",
        ]
    )

    for x in ellipsis_x:
        lines.append(
            rf"\node[ellipsisfont] at ({shift + x:.2f}, {TOP_Y:.2f}) {{$\cdots$}};"
        )
        lines.append(
            rf"\node[ellipsisfont] at ({shift + x:.2f}, {BOTTOM_Y:.2f}) {{$\cdots$}};"
        )

    for x, label, color, style in zip(visible_x, top_labels, top_colors, top_styles):
        lines.extend(
            vertex_lines(
                shift + x,
                TOP_Y,
                TOP_LABEL_Y,
                label,
                color,
                style,
                "south",
            )
        )

    for x, label, color, style in zip(visible_x, bottom_labels, bottom_colors, bottom_styles):
        lines.extend(
            vertex_lines(
                shift + x,
                BOTTOM_Y,
                BOTTOM_LABEL_Y,
                label,
                color,
                style,
                "north",
            )
        )

    return lines


def tikz_source() -> str:
    all_lines: list[str] = []
    for shift, panel in zip(PANEL_SHIFTS, PANELS):
        all_lines.extend(panel_lines(shift, panel))

    nodes = "\n".join(all_lines)
    return rf"""\documentclass{{article}}
\usepackage[paperwidth=58.0cm,paperheight=10.4cm,margin=0pt]{{geometry}}
\usepackage{{amsmath}}
\usepackage{{graphicx}}
\usepackage{{tikz}}
\definecolor{{firstcolor}}{{RGB}}{{141, 214, 168}}
\definecolor{{lastcolor}}{{RGB}}{{245, 213, 110}}
\definecolor{{blocka}}{{RGB}}{{59, 129, 239}}
\definecolor{{blockb}}{{RGB}}{{255, 128, 58}}
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
    lightvertex/.style={{
        vertex,
        fill opacity=0.45
    }},
    labelfont/.style={{
        font=\fontsize{{26}}{{26}}\selectfont
    }},
    titlefont/.style={{
        font=\fontsize{{28}}{{28}}\selectfont
    }},
    ellipsisfont/.style={{
        font=\fontsize{{28}}{{28}}\selectfont
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
        else Path(__file__).with_name("partition_generators.jpg")
    )
    output.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as temp_dir:
        workdir = Path(temp_dir)
        tex_path = workdir / "partition_generators.tex"
        tex_path.write_text(tikz_source(), encoding="utf-8")
        pdf_path = build_pdf(workdir, tex_path)
        convert_pdf_to_jpg(pdf_path, output)

    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
