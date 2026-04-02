"""Generate a JPG with two n=3 partition diagrams and their outer-product formulas."""

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
FORMULA_Y = -4.25

X_POS = [0.0, 4.0, 8.0]
PANEL_SHIFTS = [2.5, 22.5]

PANELS = [
    {
        "top_colors": ["blocka", "blocka", "blockb"],
        "bottom_colors": ["blocka", "blockb", "blockc"],
        "formula": (
            r"$\displaystyle \mathbf{S}(b_{\scriptscriptstyle 1}\,s_{\scriptscriptstyle 2}\,p_{\scriptscriptstyle 3})"
            r"=\sum_{x,y,z}\,\ket{xxy}\!\bra{xyz}$"
        ),
    },
    {
        "top_colors": ["blockd", "blocke", "blockf"],
        "bottom_colors": ["blocke", "blockd", "blockg"],
        "formula": (
            r"$\displaystyle \mathbf{S}(s_{\scriptscriptstyle 1}\,p_{\scriptscriptstyle 3})"
            r"=\sum_{x,y,z,w}\,\ket{xyz}\!\bra{yxw}$"
        ),
    },
]


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


def component_groups(panel: dict[str, object]) -> list[tuple[str, list[int]]]:
    groups: dict[str, list[int]] = {}
    for node, color in zip([1, 2, 3, -1, -2, -3], panel["top_colors"] + panel["bottom_colors"]):
        groups.setdefault(color, []).append(node)
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


def edge_lines(shift: float, panel: dict[str, object]) -> list[str]:
    lines: list[str] = []
    center_y = (TOP_Y + BOTTOM_Y) / 2

    def coords(node: int) -> tuple[float, float]:
        x_idx = abs(node) - 1
        x = shift + X_POS[x_idx]
        y = TOP_Y if node > 0 else BOTTOM_Y
        return x, y

    for color, component in component_groups(panel):
        edges = cycle_edges(component)
        curved_edges = [
            (node_a, node_b)
            for node_a, node_b in edges
            if node_row(node_a) == node_row(node_b) and edge_priority(node_a, node_b) != 0
        ]
        curved_offsets = curve_offsets(len(curved_edges))
        curve_idx = 0
        style = rf"componentedge, draw={edge_draw_color(color)}"
        for node_a, node_b in edges:
            x_a, y_a = coords(node_a)
            x_b, y_b = coords(node_b)
            if node_row(node_a) != node_row(node_b):
                lines.append(rf"\draw[{style}] ({x_a:.2f}, {y_a:.2f}) -- ({x_b:.2f}, {y_b:.2f});")
                continue
            if edge_priority(node_a, node_b) == 0:
                lines.append(rf"\draw[{style}] ({x_a:.2f}, {y_a:.2f}) -- ({x_b:.2f}, {y_b:.2f});")
            else:
                bend_y = center_y + curved_offsets[curve_idx]
                curve_idx += 1
                lines.append(
                    rf"\draw[{style}] ({x_a:.2f}, {y_a:.2f}) .. controls ({x_a:.2f}, {bend_y:.2f}) "
                    rf"and ({x_b:.2f}, {bend_y:.2f}) .. ({x_b:.2f}, {y_b:.2f});"
                )
    return lines


def vertex_lines(shift: float, node: int, color: str) -> list[str]:
    x_idx = abs(node) - 1
    x = shift + X_POS[x_idx]
    y = TOP_Y if node > 0 else BOTTOM_Y
    anchor = "south" if node > 0 else "north"
    label_y = TOP_LABEL_Y if node > 0 else BOTTOM_LABEL_Y
    return [
        rf"\node[vertexmask] at ({x:.2f}, {y:.2f}) {{}};",
        rf"\node[vertex, fill={color}] at ({x:.2f}, {y:.2f}) {{}};",
        rf"\node[labelfont, anchor={anchor}] at ({x:.2f}, {label_y:.2f}) "
        rf"{{{latex_label(node)}}};",
    ]


def panel_lines(shift: float, panel: dict[str, object]) -> list[str]:
    lines: list[str] = edge_lines(shift, panel)
    top_colors = panel["top_colors"]
    bottom_colors = panel["bottom_colors"]

    for node, color in zip([1, 2, 3], top_colors):
        lines.extend(vertex_lines(shift, node, color))

    for node, color in zip([-1, -2, -3], bottom_colors):
        lines.extend(vertex_lines(shift, node, color))

    lines.append(
        rf"\node[formula, align=center] at ({shift + 4.0:.2f}, {FORMULA_Y:.2f}) "
        rf"{{{panel['formula']}}};"
    )
    return lines


def tikz_source() -> str:
    lines: list[str] = []
    for shift, panel in zip(PANEL_SHIFTS, PANELS):
        lines.extend(panel_lines(shift, panel))

    nodes = "\n".join(lines)
    return rf"""\documentclass{{article}}
\usepackage[paperwidth=36.0cm,paperheight=12.0cm,margin=0pt]{{geometry}}
\usepackage{{amsmath}}
\usepackage{{graphicx}}
\usepackage{{tikz}}
\newcommand{{\ket}}[1]{{\left|#1\right\rangle}}
\newcommand{{\bra}}[1]{{\left\langle#1\right|}}
\definecolor{{blocka}}{{RGB}}{{59, 129, 239}}
\definecolor{{blockb}}{{RGB}}{{255, 128, 58}}
\definecolor{{blockc}}{{RGB}}{{123, 214, 146}}
\definecolor{{blockd}}{{RGB}}{{255, 184, 205}}
\definecolor{{blocke}}{{RGB}}{{136, 202, 255}}
\definecolor{{blockf}}{{RGB}}{{255, 196, 122}}
\definecolor{{blockg}}{{RGB}}{{156, 231, 175}}
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
    }},
    formula/.style={{
        font=\fontsize{{34}}{{34}}\selectfont
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
        else Path(__file__).with_name("partition_schur_examples_n3.jpg")
    )
    output.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as temp_dir:
        workdir = Path(temp_dir)
        tex_path = workdir / "partition_schur_examples_n3.tex"
        tex_path.write_text(tikz_source(), encoding="utf-8")
        pdf_path = build_pdf(workdir, tex_path)
        convert_pdf_to_jpg(pdf_path, output)

    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
