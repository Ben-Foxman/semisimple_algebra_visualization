"""Generate a JPG partition diagram using LaTeX/TikZ."""

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

BLOCK_COLORS = {
    1: "blocka",
    2: "blocka",
    -4: "blocka",
    3: "blockb",
    -1: "blockb",
    -2: "blockb",
    4: "blockc",
    -3: "blockc",
}

POSITIONS = {
    1: (0.0, TOP_Y),
    2: (4.0, TOP_Y),
    3: (8.0, TOP_Y),
    4: (12.0, TOP_Y),
    -1: (0.0, BOTTOM_Y),
    -2: (4.0, BOTTOM_Y),
    -3: (8.0, BOTTOM_Y),
    -4: (12.0, BOTTOM_Y),
}

ORDER = [1, 2, 3, 4, -1, -2, -3, -4]


def latex_label(node: int) -> str:
    if node > 0:
        return rf"$\displaystyle {node}$"
    return rf"$\displaystyle {-node}'$"


def node_key(node: int) -> tuple[int, int]:
    return (0 if node > 0 else 1, abs(node) - 1)


def node_row(node: int) -> str:
    return "top" if node > 0 else "bottom"


def edge_priority(node_a: int, node_b: int) -> int:
    if node_row(node_a) != node_row(node_b):
        return 2
    if abs(abs(node_a) - abs(node_b)) == 1:
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


def component_groups() -> list[tuple[str, list[int]]]:
    groups: dict[str, list[int]] = {}
    for node in ORDER:
        groups.setdefault(BLOCK_COLORS[node], []).append(node)
    return [(color, sorted(nodes, key=node_key)) for color, nodes in groups.items()]


def cycle_edges(component: list[int]) -> list[tuple[int, int]]:
    if len(component) < 2:
        return []
    if len(component) == 2:
        return [(component[0], component[1])]

    start = component[0]
    best_score: tuple[object, ...] | None = None
    best_edges: list[tuple[int, int]] = []

    for perm in permutations(component[1:]):
        cycle = [start, *perm]
        edges = [(cycle[idx], cycle[(idx + 1) % len(cycle)]) for idx in range(len(cycle))]
        priorities = tuple(sorted(edge_priority(a, b) for a, b in edges))
        score = (priorities, tuple(sorted(tuple(sorted(edge)) for edge in edges)), tuple(cycle))
        if best_score is None or score < best_score:
            best_score = score
            best_edges = edges

    return best_edges


def edge_lines() -> list[str]:
    lines: list[str] = []
    for color, component in component_groups():
        curved_edges = [
            (node_a, node_b)
            for node_a, node_b in cycle_edges(component)
            if node_row(node_a) == node_row(node_b) and edge_priority(node_a, node_b) != 0
        ]
        curved_offsets = curve_offsets(len(curved_edges))
        curve_idx = 0
        for node_a, node_b in cycle_edges(component):
            x_a, y_a = POSITIONS[node_a]
            x_b, y_b = POSITIONS[node_b]
            style = rf"componentedge, draw={edge_draw_color(color)}"
            if node_row(node_a) != node_row(node_b):
                lines.append(rf"\draw[{style}] ({x_a:.1f}, {y_a:.1f}) -- ({x_b:.1f}, {y_b:.1f});")
                continue
            if edge_priority(node_a, node_b) == 0:
                lines.append(rf"\draw[{style}] ({x_a:.1f}, {y_a:.1f}) -- ({x_b:.1f}, {y_b:.1f});")
            else:
                center_y = (TOP_Y + BOTTOM_Y) / 2 + curved_offsets[curve_idx]
                curve_idx += 1
                lines.append(
                    rf"\draw[{style}] ({x_a:.1f}, {y_a:.1f}) .. controls ({x_a:.1f}, {center_y:.2f}) "
                    rf"and ({x_b:.1f}, {center_y:.1f}) .. ({x_b:.1f}, {y_b:.1f});"
                )
    return lines


def vertex_lines(node: int) -> list[str]:
    x, y = POSITIONS[node]
    anchor = "south" if node > 0 else "north"
    label_y = TOP_LABEL_Y if node > 0 else BOTTOM_LABEL_Y
    return [
        rf"\node[vertexmask] at ({x:.1f}, {y:.1f}) {{}};",
        rf"\node[vertex, fill={BLOCK_COLORS[node]}] at ({x:.1f}, {y:.1f}) {{}};",
        rf"\node[labelfont, anchor={anchor}] at ({x:.1f}, {label_y:.2f}) {{{latex_label(node)}}};",
    ]


def tikz_source() -> str:
    nodes = "\n".join(edge_lines() + [line for node in ORDER for line in vertex_lines(node)])
    return rf"""\documentclass{{article}}
\usepackage[paperwidth=17.0cm,paperheight=10.4cm,margin=0pt]{{geometry}}
\usepackage{{amsmath}}
\usepackage{{graphicx}}
\usepackage{{tikz}}
\definecolor{{blocka}}{{RGB}}{{59, 129, 239}}
\definecolor{{blockb}}{{RGB}}{{255, 128, 58}}
\definecolor{{blockc}}{{RGB}}{{123, 214, 146}}
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
    labelfont/.style={{
        font=\fontsize{{26}}{{26}}\selectfont
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
        else Path(__file__).with_name("partition_n4_three_blocks.jpg")
    )
    output.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as temp_dir:
        workdir = Path(temp_dir)
        tex_path = workdir / "partition_n4_three_blocks.tex"
        tex_path.write_text(tikz_source(), encoding="utf-8")
        pdf_path = build_pdf(workdir, tex_path)
        convert_pdf_to_jpg(pdf_path, output)

    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
