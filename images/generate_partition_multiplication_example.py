"""Generate a JPG illustrating multiplication of two n=4 partition diagrams."""

from __future__ import annotations

import subprocess
import sys
import tempfile
from itertools import permutations
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
RESULT_SHIFT = 61.0
SCALAR_X = 59.2
OPERATOR_Y = 1.85

COLOR_RGB = {
    "merged1": (136, 202, 255),
    "merged2": (0, 82, 179),
    "survive": (255, 196, 122),
    "internal1": (156, 231, 175),
    "internal2": (0, 113, 68),
    "loop1": (255, 184, 205),
    "loop2": (186, 88, 20),
    "merged_final": (0, 210, 255),
}
STACK_ROW_ORDER = {
    "top": 0,
    "upper": 1,
    "lower": 2,
    "bottom": 3,
}
STACK_ROW_Y = {
    "top": STACK_TOP_Y,
    "upper": STACK_UPPER_Y,
    "lower": STACK_LOWER_Y,
    "bottom": STACK_BOTTOM_Y,
}


def label_text(node: int) -> str:
    if node > 0:
        return rf"$\displaystyle {node}$"
    return rf"$\displaystyle {-node}'$"


def edge_color_name(color: str) -> str:
    return color


def result_color_rgbs() -> dict[str, tuple[int, int, int]]:
    return {
        "merged_final": COLOR_RGB["merged_final"],
        "survive": COLOR_RGB["survive"],
    }


def curve_bend(center_y: float, curve_idx: int) -> float:
    if curve_idx == 0:
        return center_y
    delta = 0.18 * ((curve_idx + 1) // 2)
    return center_y + delta if curve_idx % 2 else center_y - delta


def node_order(node: tuple[str, int]) -> tuple[int, int]:
    row, idx = node
    return (0 if row == "top" else 1, idx)


def node_xy(shift: float, node: tuple[str, int]) -> tuple[float, float]:
    row, idx = node
    y = TOP_Y if row == "top" else BOTTOM_Y
    return shift + X_POS[idx], y


def edge_priority(
    node_a: tuple[str, int],
    node_b: tuple[str, int],
) -> int:
    row_a, idx_a = node_a
    row_b, idx_b = node_b
    if row_a != row_b:
        return 2
    if abs(idx_a - idx_b) == 1:
        return 0
    return 1


def component_cycle_edges(
    shift: float,
    component: list[tuple[str, int]],
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
        priorities = tuple(sorted(edge_priority(a, b) for a, b in edges))
        lengths = tuple(
            sorted(
                (node_xy(shift, a)[0] - node_xy(shift, b)[0]) ** 2
                + (node_xy(shift, a)[1] - node_xy(shift, b)[1]) ** 2
                for a, b in edges
            )
        )
        score = (priorities, lengths, tuple(cycle))
        if best_score is None or score < best_score:
            best_score = score
            best_edges = edges

    return best_edges


def component_lines(shift: float, top_colors: list[str], bottom_colors: list[str]) -> list[str]:
    components: dict[str, list[tuple[str, int]]] = {}
    for idx, color in enumerate(top_colors):
        components.setdefault(color, []).append(("top", idx))
    for idx, color in enumerate(bottom_colors):
        components.setdefault(color, []).append(("bottom", idx))

    lines: list[str] = []
    for color_name, component in components.items():
        for node_a, node_b in component_cycle_edges(shift, component):
            x_a, y_a = node_xy(shift, node_a)
            x_b, y_b = node_xy(shift, node_b)
            lines.append(
                rf"\draw[componentedge, draw={edge_color_name(color_name)}] "
                rf"({x_a:.2f}, {y_a:.2f}) -- ({x_b:.2f}, {y_b:.2f});"
            )

    return lines


def solid_diagram_lines(
    shift: float,
    top_colors: list[str],
    bottom_colors: list[str],
) -> list[str]:
    lines: list[str] = component_lines(shift, top_colors, bottom_colors)
    for x, node, color in zip(X_POS, [1, 2, 3, 4], top_colors):
        lines.append(
            rf"\node[vertex, fill={color}] at ({shift + x:.2f}, {TOP_Y:.2f}) {{}};"
        )
        lines.append(
            rf"\node[labelfont, anchor=south] at ({shift + x:.2f}, {TOP_LABEL_Y:.2f}) {{\scalebox{{1.65}}{{{label_text(node)}}}}};"
        )
    for x, node, color in zip(X_POS, [-1, -2, -3, -4], bottom_colors):
        lines.append(
            rf"\node[vertex, fill={color}] at ({shift + x:.2f}, {BOTTOM_Y:.2f}) {{}};"
        )
        lines.append(
            rf"\node[labelfont, anchor=north] at ({shift + x:.2f}, {BOTTOM_LABEL_Y:.2f}) {{\scalebox{{1.65}}{{{label_text(node)}}}}};"
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


def stacked_node_order(node: tuple[str, int]) -> tuple[int, int]:
    row, idx = node
    return (STACK_ROW_ORDER[row], idx)


def stacked_node_xy(node: tuple[str, int]) -> tuple[float, float]:
    row, idx = node
    return STACK_SHIFT + X_POS[idx], STACK_ROW_Y[row]


def stacked_edge_priority(
    node_a: tuple[str, int],
    node_b: tuple[str, int],
) -> int:
    row_a, idx_a = node_a
    row_b, idx_b = node_b
    if row_a != row_b:
        return 2
    if abs(idx_a - idx_b) == 1:
        return 0
    return 1


def stacked_component_cycle_edges(
    component: list[tuple[str, int]],
) -> list[tuple[tuple[str, int], tuple[str, int]]]:
    ordered = sorted(component, key=stacked_node_order)
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
        priorities = tuple(sorted(stacked_edge_priority(a, b) for a, b in edges))
        lengths = tuple(
            sorted(
                (stacked_node_xy(a)[0] - stacked_node_xy(b)[0]) ** 2
                + (stacked_node_xy(a)[1] - stacked_node_xy(b)[1]) ** 2
                for a, b in edges
            )
        )
        score = (priorities, lengths, tuple(cycle))
        if best_score is None or score < best_score:
            best_score = score
            best_edges = edges

    return best_edges


def stacked_component_lines() -> list[str]:
    top_components: dict[str, list[tuple[str, int]]] = {}
    for idx, color in enumerate(["merged1", "merged1", "survive", "survive"]):
        top_components.setdefault(color, []).append(("top", idx))
    for idx, color in enumerate(["loop1", "merged1", "internal1", "internal1"]):
        top_components.setdefault(color, []).append(("upper", idx))

    bottom_components: dict[str, list[tuple[str, int]]] = {}
    for idx, color in enumerate(["loop2", "merged2", "internal2", "internal2"]):
        bottom_components.setdefault(color, []).append(("lower", idx))
    for idx, color in enumerate(["merged2", "merged2", "merged2", "merged2"]):
        bottom_components.setdefault(color, []).append(("bottom", idx))

    lines: list[str] = []
    for components in (top_components, bottom_components):
        for color_name, component in components.items():
            curve_idx = 0
            for node_a, node_b in stacked_component_cycle_edges(component):
                x_a, y_a = stacked_node_xy(node_a)
                x_b, y_b = stacked_node_xy(node_b)
                if node_a[0] != node_b[0] or abs(node_a[1] - node_b[1]) == 1:
                    lines.append(
                        rf"\draw[componentedge, draw={edge_color_name(color_name)}] "
                        rf"({x_a:.2f}, {y_a:.2f}) -- ({x_b:.2f}, {y_b:.2f});"
                    )
                else:
                    center_y = (
                        (STACK_TOP_Y + STACK_UPPER_Y) / 2
                        if node_a[0] in {"top", "upper"}
                        else (STACK_LOWER_Y + STACK_BOTTOM_Y) / 2
                    )
                    bend_y = curve_bend(center_y, curve_idx)
                    curve_idx += 1
                    lines.append(
                        rf"\draw[componentedge, draw={edge_color_name(color_name)}] "
                        rf"({x_a:.2f}, {y_a:.2f}) .. controls ({x_a:.2f}, {bend_y:.2f}) "
                        rf"and ({x_b:.2f}, {bend_y:.2f}) .. ({x_b:.2f}, {y_b:.2f});"
                    )

    return lines


def stacked_lines() -> list[str]:
    lines: list[str] = stacked_component_lines()

    top_nodes = [
        ("merged1", "merged2"),
        ("merged1", "merged2"),
        ("survive", "survive"),
        ("survive", "survive"),
    ]
    upper_middle = [
        ("loop1", "loop2"),
        ("merged1", "merged2"),
        ("internal1", "internal2"),
        ("internal1", "internal2"),
    ]
    lower_middle = [
        ("loop1", "loop2"),
        ("merged1", "merged2"),
        ("internal1", "internal2"),
        ("internal1", "internal2"),
    ]
    bottom_nodes = [
        ("merged1", "merged2"),
        ("merged1", "merged2"),
        ("merged1", "merged2"),
        ("merged1", "merged2"),
    ]

    for x, node, colors in zip(X_POS, [1, 2, 3, 4], top_nodes):
        lines.extend(split_node_lines(STACK_SHIFT + x, STACK_TOP_Y, colors[0], colors[1]))
        lines.append(
            rf"\node[labelfont, anchor=south] at ({STACK_SHIFT + x:.2f}, {STACK_TOP_LABEL_Y:.2f}) {{\scalebox{{1.65}}{{{label_text(node)}}}}};"
        )

    for x, colors in zip(X_POS, upper_middle):
        lines.extend(
            split_node_lines(STACK_SHIFT + x, STACK_UPPER_Y, colors[0], colors[1])
        )

    for x in X_POS:
        lines.append(
            rf"\draw[dashed, dash pattern=on 2.2pt off 4.0pt, line width=1.2pt, black!80] "
            rf"({STACK_SHIFT + x:.2f}, {STACK_UPPER_Y - 1.05:.2f}) -- "
            rf"({STACK_SHIFT + x:.2f}, {STACK_LOWER_Y + 1.05:.2f});"
        )

    for x, colors in zip(X_POS, lower_middle):
        lines.extend(
            split_node_lines(STACK_SHIFT + x, STACK_LOWER_Y, colors[0], colors[1])
        )

    for x, node, colors in zip(X_POS, [-1, -2, -3, -4], bottom_nodes):
        lines.extend(
            split_node_lines(STACK_SHIFT + x, STACK_BOTTOM_Y, colors[0], colors[1])
        )
        lines.append(
            rf"\node[labelfont, anchor=north] at ({STACK_SHIFT + x:.2f}, {STACK_BOTTOM_LABEL_Y:.2f}) {{\scalebox{{1.65}}{{{label_text(node)}}}}};"
        )

    return lines


def result_lines() -> list[str]:
    top_colors = ["merged_final", "merged_final", "survive", "survive"]
    bottom_colors = ["merged_final", "merged_final", "merged_final", "merged_final"]
    return solid_diagram_lines(RESULT_SHIFT, top_colors, bottom_colors)


def tikz_source() -> str:
    lines: list[str] = []
    lines.extend(
        solid_diagram_lines(
            D1_SHIFT,
            ["merged1", "merged1", "survive", "survive"],
            ["loop1", "merged1", "internal1", "internal1"],
        )
    )
    lines.extend(
        solid_diagram_lines(
            D2_SHIFT,
            ["loop2", "merged2", "internal2", "internal2"],
            ["merged2", "merged2", "merged2", "merged2"],
        )
    )
    lines.extend(stacked_lines())
    lines.extend(result_lines())

    lines.extend(
        [
            rf"\node[opfont] at ({TIMES_X:.2f}, {OPERATOR_Y:.2f}) {{\scalebox{{2.0}}{{$\times$}}}};",
            rf"\node[opfont] at ({EQUAL1_X:.2f}, {OPERATOR_Y:.2f}) {{\scalebox{{2.0}}{{$=$}}}};",
            rf"\node[opfont] at ({EQUAL2_X:.2f}, {OPERATOR_Y:.2f}) {{\scalebox{{2.0}}{{$=$}}}};",
            rf"\node[scalarfont, anchor=east] at ({SCALAR_X:.2f}, {OPERATOR_Y:.2f}) {{\scalebox{{1.8}}{{$d^2$}}}};",
            rf"\node[panelnamefont] at ({D1_SHIFT + 6.0:.2f}, {PANEL_LABEL_Y:.2f}) {{\scalebox{{1.8}}{{$D_1$}}}};",
            rf"\node[panelnamefont] at ({D2_SHIFT + 6.0:.2f}, {PANEL_LABEL_Y:.2f}) {{\scalebox{{1.8}}{{$D_2$}}}};",
            rf"\node[panelnamefont] at ({RESULT_SHIFT + 6.0:.2f}, {PANEL_LABEL_Y:.2f}) {{\scalebox{{1.8}}{{$D_1D_2$}}}};",
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
\usepackage[paperwidth=78cm,paperheight=18.5cm,margin=0pt]{{geometry}}
\usepackage{{amsmath}}
\usepackage{{graphicx}}
\usepackage{{tikz}}
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
    scalarfont/.style={{
        font=\fontsize{{86}}{{86}}\selectfont
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
        else Path(__file__).with_name("partition_multiplication_example.jpg")
    )
    output.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as temp_dir:
        workdir = Path(temp_dir)
        tex_path = workdir / "partition_multiplication_example.tex"
        tex_path.write_text(tikz_source(), encoding="utf-8")
        pdf_path = build_pdf(workdir, tex_path)
        convert_pdf_to_jpg(pdf_path, output)

    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
