"""Generate a JPG illustrating multiplication of two n=4 permutation diagrams."""

from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path


X_POS = [0.0, 4.0, 8.0, 12.0]
TOP_Y = 3.2
BOTTOM_Y = 0.0
TOP_LABEL_Y = 5.10
BOTTOM_LABEL_Y = -2.00
PANEL_LABEL_Y = -5.60

STACK_TOP_Y = 5.8
STACK_UPPER_Y = 3.4
STACK_LOWER_Y = -0.2
STACK_BOTTOM_Y = -2.6
STACK_TOP_LABEL_Y = 7.20
STACK_BOTTOM_LABEL_Y = -4.15

D1_SHIFT = 0.0
TIMES_X = 15.8
D2_SHIFT = 18.8
EQUAL1_X = 34.8
STACK_SHIFT = 37.8
EQUAL2_X = 55.8
RESULT_SHIFT = 58.8
OPERATOR_Y = 1.85

D1_PERM = [2, 1, 4, 3]
D2_PERM = [3, 4, 1, 2]

D1_TOP_COLORS = ["lightblue", "lightorange", "lightgreen", "lightyellow"]
D2_TOP_COLORS = ["darkorange", "darkblue", "darkyellow", "darkgreen"]
RESULT_TOP_COLORS = ["resulta", "resultb", "resultc", "resultd"]

COLOR_RGB = {
    "lightblue": (136, 202, 255),
    "lightorange": (255, 184, 122),
    "lightgreen": (156, 231, 175),
    "lightyellow": (255, 219, 122),
    "darkblue": (20, 92, 196),
    "darkorange": (201, 103, 28),
    "darkgreen": (18, 137, 84),
    "darkyellow": (214, 165, 34),
}

def label_text(node: int) -> str:
    if node > 0:
        return rf"$\displaystyle {node}$"
    return rf"$\displaystyle {-node}'$"


def bottom_colors(top_colors: list[str], perm: list[int]) -> list[str]:
    colors = [""] * len(top_colors)
    for idx, target in enumerate(perm):
        colors[target - 1] = top_colors[idx]
    return colors


def compose_perm(first: list[int], second: list[int]) -> list[int]:
    return [second[target - 1] for target in first]


def average_rgb(color_a: str, color_b: str) -> tuple[int, int, int]:
    rgb_a = COLOR_RGB[color_a]
    rgb_b = COLOR_RGB[color_b]
    return tuple((a + b) // 2 for a, b in zip(rgb_a, rgb_b))


def result_color_rgbs() -> dict[str, tuple[int, int, int]]:
    colors: dict[str, tuple[int, int, int]] = {}
    for idx, result_name in enumerate(RESULT_TOP_COLORS):
        joined_idx = D1_PERM[idx] - 1
        colors[result_name] = average_rgb(D1_TOP_COLORS[idx], D2_TOP_COLORS[joined_idx])
    return colors


def edge_color_name(color: str) -> str:
    return color


def edge_lines(
    shift: float,
    top_y: float,
    bottom_y: float,
    perm: list[int],
    edge_colors: list[str],
) -> list[str]:
    lines: list[str] = []
    for idx, target in enumerate(perm):
        lines.append(
            rf"\draw[componentedge, draw={edge_color_name(edge_colors[idx])}] "
            rf"({shift + X_POS[idx]:.2f}, {top_y:.2f}) -- "
            rf"({shift + X_POS[target - 1]:.2f}, {bottom_y:.2f});"
        )
    return lines


def vertex_lines(
    x: float,
    y: float,
    label_y: float,
    label: str,
    color: str,
    anchor: str,
) -> list[str]:
    return [
        rf"\node[vertexmask] at ({x:.2f}, {y:.2f}) {{}};",
        rf"\node[vertex, fill={color}] at ({x:.2f}, {y:.2f}) {{}};",
        rf"\node[labelfont, anchor={anchor}] at ({x:.2f}, {label_y:.2f}) "
        rf"{{\scalebox{{1.65}}{{{label}}}}};",
    ]


def solid_diagram_lines(
    shift: float,
    top_colors: list[str],
    bottom_colors_row: list[str],
    perm: list[int],
) -> list[str]:
    lines = edge_lines(shift, TOP_Y, BOTTOM_Y, perm, top_colors)
    for x, node, color in zip(X_POS, [1, 2, 3, 4], top_colors):
        lines.extend(
            vertex_lines(
                shift + x,
                TOP_Y,
                TOP_LABEL_Y,
                label_text(node),
                color,
                "south",
            )
        )
    for x, node, color in zip(X_POS, [-1, -2, -3, -4], bottom_colors_row):
        lines.extend(
            vertex_lines(
                shift + x,
                BOTTOM_Y,
                BOTTOM_LABEL_Y,
                label_text(node),
                color,
                "north",
            )
        )
    return lines


def split_node_lines(x: float, y: float, left: str, right: str) -> list[str]:
    return [
        rf"\begin{{scope}}[shift={{({x:.2f},{y:.2f})}}]",
        rf"\clip (0,0) circle (0.85cm);",
        rf"\fill[{left}] (-1cm,-1cm) rectangle (0,1cm);",
        rf"\fill[{right}] (0,-1cm) rectangle (1cm,1cm);",
        rf"\end{{scope}}",
        rf"\draw[black, line width=1.3pt] ({x:.2f},{y:.2f}) circle (0.85cm);",
    ]


def stacked_lines() -> list[str]:
    d1_bottom = bottom_colors(D1_TOP_COLORS, D1_PERM)
    d2_bottom = bottom_colors(D2_TOP_COLORS, D2_PERM)

    lines = edge_lines(STACK_SHIFT, STACK_TOP_Y, STACK_UPPER_Y, D1_PERM, D1_TOP_COLORS)
    lines.extend(edge_lines(STACK_SHIFT, STACK_LOWER_Y, STACK_BOTTOM_Y, D2_PERM, D2_TOP_COLORS))

    top_nodes = [(color, color) for color in D1_TOP_COLORS]
    middle_nodes = list(zip(d1_bottom, D2_TOP_COLORS))
    bottom_nodes = [(color, color) for color in d2_bottom]

    for x, node, colors in zip(X_POS, [1, 2, 3, 4], top_nodes):
        lines.extend(split_node_lines(STACK_SHIFT + x, STACK_TOP_Y, colors[0], colors[1]))
        lines.append(
            rf"\node[labelfont, anchor=south] at ({STACK_SHIFT + x:.2f}, {STACK_TOP_LABEL_Y:.2f}) "
            rf"{{\scalebox{{1.65}}{{{label_text(node)}}}}};"
        )

    for x, colors in zip(X_POS, middle_nodes):
        lines.extend(split_node_lines(STACK_SHIFT + x, STACK_UPPER_Y, colors[0], colors[1]))

    for x in X_POS:
        lines.append(
            rf"\draw[dashed, dash pattern=on 2.2pt off 4.0pt, line width=1.2pt, black!80] "
            rf"({STACK_SHIFT + x:.2f}, {STACK_UPPER_Y - 1.05:.2f}) -- "
            rf"({STACK_SHIFT + x:.2f}, {STACK_LOWER_Y + 1.05:.2f});"
        )

    for x, colors in zip(X_POS, middle_nodes):
        lines.extend(split_node_lines(STACK_SHIFT + x, STACK_LOWER_Y, colors[0], colors[1]))

    for x, node, colors in zip(X_POS, [-1, -2, -3, -4], bottom_nodes):
        lines.extend(split_node_lines(STACK_SHIFT + x, STACK_BOTTOM_Y, colors[0], colors[1]))
        lines.append(
            rf"\node[labelfont, anchor=north] at ({STACK_SHIFT + x:.2f}, {STACK_BOTTOM_LABEL_Y:.2f}) "
            rf"{{\scalebox{{1.65}}{{{label_text(node)}}}}};"
        )

    return lines


def tikz_source() -> str:
    d1_bottom = bottom_colors(D1_TOP_COLORS, D1_PERM)
    d2_bottom = bottom_colors(D2_TOP_COLORS, D2_PERM)
    result_perm = compose_perm(D1_PERM, D2_PERM)
    result_bottom = bottom_colors(RESULT_TOP_COLORS, result_perm)

    lines: list[str] = []
    lines.extend(solid_diagram_lines(D1_SHIFT, D1_TOP_COLORS, d1_bottom, D1_PERM))
    lines.extend(solid_diagram_lines(D2_SHIFT, D2_TOP_COLORS, d2_bottom, D2_PERM))
    lines.extend(stacked_lines())
    lines.extend(solid_diagram_lines(RESULT_SHIFT, RESULT_TOP_COLORS, result_bottom, result_perm))
    lines.extend(
        [
            rf"\node[opfont] at ({TIMES_X:.2f}, {OPERATOR_Y:.2f}) {{\scalebox{{2.0}}{{$\times$}}}};",
            rf"\node[opfont] at ({EQUAL1_X:.2f}, {OPERATOR_Y:.2f}) {{\scalebox{{2.0}}{{$=$}}}};",
            rf"\node[opfont] at ({EQUAL2_X:.2f}, {OPERATOR_Y:.2f}) {{\scalebox{{2.0}}{{$=$}}}};",
            rf"\node[panelnamefont] at ({D1_SHIFT + 6.0:.2f}, {PANEL_LABEL_Y:.2f}) {{\scalebox{{1.8}}{{$\pi_1$}}}};",
            rf"\node[panelnamefont] at ({D2_SHIFT + 6.0:.2f}, {PANEL_LABEL_Y:.2f}) {{\scalebox{{1.8}}{{$\pi_2$}}}};",
            rf"\node[panelnamefont] at ({RESULT_SHIFT + 6.0:.2f}, {PANEL_LABEL_Y:.2f}) {{\scalebox{{1.8}}{{$\pi_1\pi_2$}}}};",
        ]
    )

    nodes = "\n".join(lines)
    color_defs = "\n".join(
        rf"\definecolor{{{name}}}{{RGB}}{{{rgb[0]}, {rgb[1]}, {rgb[2]}}}"
        for name, rgb in COLOR_RGB.items()
    )
    edge_defs = ""
    result_rgbs = result_color_rgbs()
    result_defs = "\n".join(
        rf"\definecolor{{{name}}}{{RGB}}{{{rgb[0]}, {rgb[1]}, {rgb[2]}}}"
        for name, rgb in result_rgbs.items()
    )
    result_edge_defs = ""
    return rf"""\documentclass{{article}}
\usepackage[paperwidth=75cm,paperheight=18.5cm,margin=0pt]{{geometry}}
\usepackage{{amsmath}}
\usepackage{{graphicx}}
\usepackage{{tikz}}
\usepackage{{xcolor}}
{color_defs}
{result_defs}
{edge_defs}
{result_edge_defs}
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
        font=\fontsize{{62}}{{62}}\selectfont
    }},
    panelnamefont/.style={{
        font=\fontsize{{68}}{{68}}\selectfont
    }},
    opfont/.style={{
        font=\fontsize{{96}}{{96}}\selectfont
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
        else Path(__file__).with_name("partition_permutation_multiplication_example.jpg")
    )
    output.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as temp_dir:
        workdir = Path(temp_dir)
        tex_path = workdir / "partition_permutation_multiplication_example.tex"
        tex_path.write_text(tikz_source(), encoding="utf-8")
        pdf_path = build_pdf(workdir, tex_path)
        convert_pdf_to_jpg(pdf_path, output)

    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
