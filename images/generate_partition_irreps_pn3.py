"""Generate a JPG showing a pn(D)=3 diagram and the partition-algebra irreps."""

from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path


TOP_Y = 7.4
BOTTOM_Y = 0.0
TOP_LABEL_Y = 10.35
BOTTOM_LABEL_Y = -2.95
PN_LABEL_Y = -11.60

LEFT_SHIFT = 0.0
X_POS = [0.0, 7.8, 15.6, 23.4]

ARROW_START_X = 32.0
ARROW_END_X = 44.0
ARROW_Y = 3.7
ARROW_LABEL_Y = 5.0

LEVEL_X = [50.6, 59.4, 69.8, 81.8, 95.4]
LEVEL_Y = {
    0: [3.7],
    1: [3.7],
    2: [6.7, 0.7],
    3: [9.7, 3.7, -2.3],
    4: [15.7, 9.7, 3.7, -2.3, -8.3],
}

PARTITIONS = {
    0: [()],
    1: [(1,)],
    2: [(2,), (1, 1)],
    3: [(3,), (2, 1), (1, 1, 1)],
    4: [(4,), (3, 1), (2, 2), (2, 1, 1), (1, 1, 1, 1)],
}

BLOCK_RGB = {
    "blocka": (59, 129, 239),
    "blockb": (255, 128, 58),
    "blockc": (123, 214, 146),
    "blockd": (255, 192, 214),
    "blocke": (245, 213, 110),
}

def latex_label(node: int) -> str:
    if node > 0:
        return rf"$\displaystyle {node}$"
    return rf"$\displaystyle {-node}'$"


def edge_color_name(color: str) -> str:
    return color


def curve_bend(center_y: float, curve_idx: int) -> float:
    if curve_idx == 0:
        return center_y
    delta = 0.24 * ((curve_idx + 1) // 2)
    return center_y + delta if curve_idx % 2 else center_y - delta


def left_component_edges() -> list[str]:
    top_positions = [("top", idx, x, TOP_Y) for idx, x in enumerate(X_POS)]
    bottom_positions = [("bottom", idx, x, BOTTOM_Y) for idx, x in enumerate(X_POS)]
    positions = top_positions + bottom_positions

    color_to_nodes: dict[str, list[tuple[str, int, float, float]]] = {
        "blocka": [positions[0], positions[5]],
        "blockb": [positions[2], positions[4]],
        "blockc": [positions[1], positions[7]],
        "blockd": [positions[3]],
        "blocke": [positions[6]],
    }

    lines: list[str] = []
    for color_name, nodes in color_to_nodes.items():
        if len(nodes) < 2:
            continue
        curve_idx = 0
        (row_a, idx_a, x_a, y_a), (row_b, idx_b, x_b, y_b) = nodes
        if row_a == row_b:
            if abs(idx_a - idx_b) == 1:
                lines.append(
                    rf"\draw[componentedge, draw={edge_color_name(color_name)}] "
                    rf"({x_a:.2f}, {y_a:.2f}) -- ({x_b:.2f}, {y_b:.2f});"
                )
            else:
                center_y = (TOP_Y + BOTTOM_Y) / 2
                bend_y = curve_bend(center_y, curve_idx)
                lines.append(
                    rf"\draw[componentedge, draw={edge_color_name(color_name)}] "
                    rf"({x_a:.2f}, {y_a:.2f}) .. controls ({x_a:.2f}, {bend_y:.2f}) and "
                    rf"({x_b:.2f}, {bend_y:.2f}) .. ({x_b:.2f}, {y_b:.2f});"
                )
        else:
            lines.append(
                rf"\draw[componentedge, draw={edge_color_name(color_name)}] "
                rf"({x_a:.2f}, {y_a:.2f}) -- ({x_b:.2f}, {y_b:.2f});"
            )

    return lines


def left_diagram_lines() -> list[str]:
    top_colors = ["blocka", "blockc", "blockb", "blockd"]
    bottom_colors = ["blockb", "blocka", "blocke", "blockc"]
    lines: list[str] = []
    ket_left_x = LEFT_SHIFT - 3.1
    ket_right_x = LEFT_SHIFT + 27.2
    ket_center_y = (TOP_Y + BOTTOM_Y) / 2
    ket_top_y = TOP_Y + 2.2
    ket_bottom_y = BOTTOM_Y - 2.2

    lines.extend(
        [
            rf"\draw[black, line width=1.5pt] ({ket_left_x:.2f}, {ket_top_y:.2f}) -- ({ket_left_x:.2f}, {ket_bottom_y:.2f});",
            rf"\draw[black, line width=1.5pt] ({ket_right_x:.2f}, {ket_center_y:.2f}) -- ({ket_right_x - 0.9:.2f}, {ket_top_y:.2f});",
            rf"\draw[black, line width=1.5pt] ({ket_right_x:.2f}, {ket_center_y:.2f}) -- ({ket_right_x - 0.9:.2f}, {ket_bottom_y:.2f});",
        ]
    )

    lines.extend(left_component_edges())

    for x, node, color in zip(X_POS, [1, 2, 3, 4], top_colors):
        lines.append(
            rf"\node[vertex, fill={color}] at ({LEFT_SHIFT + x:.2f}, {TOP_Y:.2f}) {{}};"
        )
        lines.append(
            rf"\node[labelfont, anchor=south] at ({LEFT_SHIFT + x:.2f}, {TOP_LABEL_Y:.2f}) "
            rf"{{{latex_label(node)}}};"
        )

    for x, node, color in zip(X_POS, [-1, -2, -3, -4], bottom_colors):
        lines.append(
            rf"\node[vertex, fill={color}] at ({LEFT_SHIFT + x:.2f}, {BOTTOM_Y:.2f}) {{}};"
        )
        lines.append(
            rf"\node[labelfont, anchor=north] at ({LEFT_SHIFT + x:.2f}, {BOTTOM_LABEL_Y:.2f}) "
            rf"{{{latex_label(node)}}};"
        )

    lines.append(
            rf"\node[pnfont] at ({LEFT_SHIFT + 11.7:.2f}, {PN_LABEL_Y:.2f}) "
            rf"{{\scalebox{{1.75}}{{$\displaystyle \mathrm{{pn}}(D)=3$}}}};"
    )
    return lines


def ket_dimensions(partition: tuple[int, ...]) -> tuple[float, float]:
    if not partition:
        return (7.2, 2.8)
    width = 1.1 * partition[0] + 6.2
    height = 1.1 * len(partition) + 1.9
    return (width, height)


def ket_fill(level: int) -> str | None:
    if level <= 2:
        return "lightred"
    if level == 3:
        return "strongred"
    return None


def young_ket_lines(center_x: float, center_y: float, partition: tuple[int, ...], level: int) -> list[str]:
    width, height = ket_dimensions(partition)
    left_x = center_x - width / 2
    right_x = center_x + width / 2
    top_y = center_y + height / 2
    bottom_y = center_y - height / 2

    lines: list[str] = []
    fill = ket_fill(level)
    if fill is not None:
        lines.append(
            rf"\fill[{fill}, rounded corners=0.14cm] "
            rf"({left_x - 0.35:.2f}, {bottom_y - 0.28:.2f}) rectangle "
            rf"({right_x + 0.60:.2f}, {top_y + 0.28:.2f});"
        )

    lines.extend(
        [
            rf"\draw[black, line width=1.2pt] ({left_x - 0.22:.2f}, {top_y + 0.12:.2f}) -- "
            rf"({left_x - 0.22:.2f}, {bottom_y - 0.12:.2f});",
            rf"\draw[black, line width=1.2pt] ({right_x + 0.30:.2f}, {center_y:.2f}) -- "
            rf"({right_x + 0.02:.2f}, {top_y + 0.12:.2f});",
            rf"\draw[black, line width=1.2pt] ({right_x + 0.30:.2f}, {center_y:.2f}) -- "
            rf"({right_x + 0.02:.2f}, {bottom_y - 0.12:.2f});",
        ]
    )

    box = 1.1
    left = center_x - width / 2 + 0.75

    if not partition:
        diagram_center_x = left + 0.70
        diagram_visual_right = diagram_center_x + 0.55
        lines.append(
            rf"\node[emptyfont] at ({diagram_center_x:.2f}, {center_y:.2f}) {{$\emptyset$}};"
        )
    else:
        total_h = len(partition) * box
        top = center_y + total_h / 2
        diagram_visual_right = left + partition[0] * box
        for row_index, row_len in enumerate(partition):
            y_top = top - row_index * box
            for col_index in range(row_len):
                x_left = left + col_index * box
                lines.append(
                    rf"\filldraw[fill=white, draw=black, line width=1.0pt] "
                    rf"({x_left:.2f}, {y_top:.2f}) rectangle ({x_left + box:.2f}, {y_top - box:.2f});"
                )

    item_gap = 1.95
    comma_gap = 0.62
    i_x = diagram_visual_right + item_gap
    j_x = i_x + item_gap
    comma1_x = diagram_visual_right + comma_gap
    comma2_x = i_x + comma_gap

    lines.append(
        rf"\node[ketcomma, anchor=base] at ({comma1_x:.2f}, {center_y - 0.12:.2f}) {{,}};"
    )
    lines.append(
        rf"\node[ketcontext, anchor=mid] at ({i_x:.2f}, {center_y - 0.28:.2f}) {{$i$}};"
    )
    lines.append(
        rf"\node[ketcomma, anchor=base] at ({comma2_x:.2f}, {center_y - 0.12:.2f}) {{,}};"
    )
    lines.append(
        rf"\node[ketcontext, anchor=mid] at ({j_x:.2f}, {center_y - 0.28:.2f}) {{$j$}};"
    )
    return lines


def irreps_lines() -> list[str]:
    lines: list[str] = []
    for level, partitions in PARTITIONS.items():
        for center_y, partition in zip(LEVEL_Y[level], partitions):
            lines.extend(young_ket_lines(LEVEL_X[level], center_y, partition, level))
    return lines


def tikz_source() -> str:
    lines = left_diagram_lines()
    lines.extend(
        [
            rf"\draw[-{{Latex[length=12pt,width=14pt]}}, line width=1.9pt] "
            rf"({ARROW_START_X:.2f}, {ARROW_Y:.2f}) -- ({ARROW_END_X:.2f}, {ARROW_Y:.2f});",
            rf"\node[arrowfont, anchor=south] at ({(ARROW_START_X + ARROW_END_X) / 2:.2f}, {ARROW_LABEL_Y:.2f}) "
            rf"{{\scalebox{{1.45}}{{$\displaystyle \widetilde{{\texttt{{FT}}_{{A}}}}$}}}};",
        ]
    )
    lines.extend(irreps_lines())

    nodes = "\n".join(lines)
    color_defs = "\n".join(
        rf"\definecolor{{{name}}}{{RGB}}{{{rgb[0]}, {rgb[1]}, {rgb[2]}}}"
        for name, rgb in BLOCK_RGB.items()
    )
    edge_defs = ""
    return rf"""\documentclass{{article}}
\usepackage[paperwidth=55cm,paperheight=21cm,margin=0pt]{{geometry}}
\usepackage{{amsmath}}
\usepackage{{graphicx}}
\usepackage{{tikz}}
\usetikzlibrary{{arrows.meta}}
{color_defs}
{edge_defs}
\definecolor{{lightred}}{{RGB}}{{255, 223, 223}}
\definecolor{{strongred}}{{RGB}}{{255, 122, 122}}
\pagestyle{{empty}}

\begin{{document}}
\thispagestyle{{empty}}
\vspace*{{\fill}}
\begin{{center}}
\begin{{tikzpicture}}[
    x=0.45cm,
    y=0.45cm,
    componentedge/.style={{
        line width=5.0pt,
        line cap=round
    }},
    vertex/.style={{
        circle,
        draw=black,
        line width=1.3pt,
        minimum size=2.10cm,
        inner sep=0pt
    }},
    labelfont/.style={{
        font=\fontsize{{54}}{{54}}\selectfont
    }},
    pnfont/.style={{
        font=\fontsize{{78}}{{78}}\selectfont
    }},
    arrowfont/.style={{
        font=\fontsize{{78}}{{78}}\selectfont
    }},
    emptyfont/.style={{
        font=\fontsize{{40}}{{40}}\selectfont
    }},
    ketcontext/.style={{
        font=\fontsize{{33}}{{33}}\selectfont
    }},
    ketcomma/.style={{
        font=\fontsize{{33}}{{33}}\selectfont
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
        else Path(__file__).with_name("partition_irreps_pn3.jpg")
    )
    output.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as temp_dir:
        workdir = Path(temp_dir)
        tex_path = workdir / "partition_irreps_pn3.tex"
        tex_path.write_text(tikz_source(), encoding="utf-8")
        pdf_path = build_pdf(workdir, tex_path)
        convert_pdf_to_jpg(pdf_path, output)

    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
