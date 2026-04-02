"""Generate a JPG illustrating Brauer and walled Brauer no-propagation rules."""

from __future__ import annotations

import subprocess
import sys
import tempfile
from itertools import permutations
from pathlib import Path


ROW_TOPS = [17.2, 3.4]
LEFT_SHIFT = 1.8
RIGHT_SHIFT = 24.4
ARROW_LEFT_X = 17.15
ARROW_RIGHT_X = 23.05
ARROW_TEXT_X = 20.10
COL_STEP = 2.75
TOP_ROW_OFFSET = 3.2
BOTTOM_ROW_OFFSET = 0.0
TOP_LABEL_OFFSET = 4.15
BOTTOM_LABEL_OFFSET = -1.08
ROW_TITLE_OFFSET = 8.45
ROW_CAPTION_OFFSET = 6.30
BRACE_OFFSET = 5.65
ARROW_Y_OFFSET = 1.65
PERMUTATION_Y_OFFSET = 3.35


def x_positions(count: int) -> list[float]:
    return [idx * COL_STEP for idx in range(count)]


def top_label(col: int) -> str:
    return rf"$\displaystyle {col}$"


def bottom_label(col: int) -> str:
    return rf"$\displaystyle {col}'$"


def wall_x(x_shift: float, count: int, r: int) -> float:
    xs = x_positions(count)
    return x_shift + (xs[r - 1] + xs[r]) / 2


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
    y_shift: float,
    node: tuple[str, int],
    xs: list[float],
) -> tuple[float, float]:
    row, idx = node
    y = y_shift + (TOP_ROW_OFFSET if row == "top" else BOTTOM_ROW_OFFSET)
    return x_shift + xs[idx], y


def component_cycle_edges(
    component: list[tuple[str, int]]
) -> list[tuple[tuple[str, int], tuple[str, int]]]:
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
    curved_offset = 0
    for color, component in components.items():
        draw_color = color
        for node_a, node_b in component_cycle_edges(component):
            x_a, y_a = node_coords(x_shift, y_shift, node_a, xs)
            x_b, y_b = node_coords(x_shift, y_shift, node_b, xs)
            if node_a[0] != node_b[0]:
                lines.append(
                    rf"\draw[componentedge, draw={draw_color}] ({x_a:.2f}, {y_a:.2f}) -- ({x_b:.2f}, {y_b:.2f});"
                )
            elif same_row_adjacent(node_a, node_b):
                lines.append(
                    rf"\draw[componentedge, draw={draw_color}] ({x_a:.2f}, {y_a:.2f}) -- ({x_b:.2f}, {y_b:.2f});"
                )
            else:
                center_y = y_shift + (TOP_ROW_OFFSET + BOTTOM_ROW_OFFSET) / 2
                bend = [0.0, 0.22, -0.22, 0.38, -0.38][curved_offset % 5]
                curved_offset += 1
                lines.append(
                    rf"\draw[componentedge, draw={draw_color}] ({x_a:.2f}, {y_a:.2f}) .. controls ({x_a:.2f}, {center_y + bend:.2f}) "
                    rf"and ({x_b:.2f}, {center_y + bend:.2f}) .. ({x_b:.2f}, {y_b:.2f});"
                )

    return lines


def draw_diagram(
    x_shift: float,
    y_shift: float,
    top_colors: list[str],
    bottom_colors: list[str],
    caption: str | None = None,
    wall_after: int | None = None,
) -> list[str]:
    lines: list[str] = []
    xs = x_positions(len(top_colors))
    lines.extend(component_lines(x_shift, y_shift, top_colors, bottom_colors))

    if wall_after is not None:
        lines.append(
            rf"\draw[walldivider] ({wall_x(x_shift, len(top_colors), wall_after):.2f}, {y_shift - 1.6:.2f}) -- "
            rf"({wall_x(x_shift, len(top_colors), wall_after):.2f}, {y_shift + 4.9:.2f});"
        )

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


def brace_segments(
    x_shift: float,
    y_shift: float,
    spans: list[tuple[int, int]],
    label: str,
) -> list[str]:
    xs = x_positions(max(b for _, b in spans))
    lines = []
    for first_col, last_col in spans:
        left_x = x_shift + xs[first_col - 1] - 0.95
        right_x = x_shift + xs[last_col - 1] + 0.95
        lines.append(rf"\draw[brace] ({left_x:.2f}, {y_shift + BRACE_OFFSET:.2f}) -- ({right_x:.2f}, {y_shift + BRACE_OFFSET:.2f});")
    label_x = x_shift + (xs[spans[0][0] - 1] + xs[spans[-1][1] - 1]) / 2
    lines.append(
        rf"\node[bracefont, anchor=south] at ({label_x:.2f}, {y_shift + BRACE_OFFSET + 0.40:.2f}) "
        rf"{{\scalebox{{1.18}}{{{label}}}}};"
    )
    return lines


def row_annotation(row_idx: int, title: str, arrow_formula: str) -> list[str]:
    y_shift = ROW_TOPS[row_idx]
    return [
        rf"\node[rowtitlefont] at ({ARROW_TEXT_X:.2f}, {y_shift + ROW_TITLE_OFFSET:.2f}) {{{title}}};",
        rf"\draw[rulearrow] ({ARROW_LEFT_X:.2f}, {y_shift + ARROW_Y_OFFSET:.2f}) -- "
        rf"({ARROW_RIGHT_X:.2f}, {y_shift + ARROW_Y_OFFSET:.2f});",
        rf"\node[permutationfont] at ({ARROW_TEXT_X:.2f}, {y_shift + PERMUTATION_Y_OFFSET:.2f}) "
        rf"{{\scalebox{{1.24}}{{{arrow_formula}}}}};",
    ]


def brauer_row_lines() -> list[str]:
    lines: list[str] = []
    y_shift = ROW_TOPS[0]

    right_top = ["blocka", "blocka", "blockb", "blockb", "blockc", "blockc"]
    right_bottom = ["blockd", "blocke", "blocke", "blockd", "blockf", "blockf"]
    left_top = ["blockc", "blocka", "blockb", "blockc", "blocka", "blockb"]
    left_bottom = ["blockd", "blockf", "blockf", "blocke", "blocke", "blockd"]

    lines.extend(
        row_annotation(
            0,
            r"{\scalebox{1.28}{$\textbf{Brauer No Propagation}$}}",
            r"$\begin{aligned}\pi_1&=(1\,5)(4\,6)\\[2pt]\pi_2&=(2'\,6')(3'\,5')\end{aligned}$",
        )
    )
    lines.extend(draw_diagram(LEFT_SHIFT, y_shift, left_top, left_bottom, r"$D_a\in B_{6}(d)$"))
    lines.extend(draw_diagram(RIGHT_SHIFT, y_shift, right_top, right_bottom))
    lines.extend(brace_lines(RIGHT_SHIFT, y_shift, 1, 4, r"$D_be_{5},\; D_b\in B_{4}(d)$"))
    return lines


def walled_row_lines() -> list[str]:
    lines: list[str] = []
    y_shift = ROW_TOPS[1]

    right_top = ["blocka", "blockb", "blocke", "blocka", "blockb", "blocke"]
    right_bottom = ["blockc", "blockd", "blockf", "blockd", "blockc", "blockf"]
    left_top = ["blocke", "blockb", "blocka", "blocke", "blockb", "blocka"]
    left_bottom = ["blockc", "blockf", "blockd", "blockd", "blockf", "blockc"]

    lines.extend(
        row_annotation(
            1,
            r"{\scalebox{1.28}{$\textbf{Walled Brauer No Propagation}$}}",
            r"$\begin{aligned}\pi_1&=(1\,3)(4\,6)\\[2pt]\pi_2&=(2'\,3')(5'\,6')\end{aligned}$",
        )
    )
    lines.extend(
        draw_diagram(
            LEFT_SHIFT,
            y_shift,
            left_top,
            left_bottom,
            r"$D_a\in B_{3,3}(d)$",
            wall_after=3,
        )
    )
    lines.extend(draw_diagram(RIGHT_SHIFT, y_shift, right_top, right_bottom, wall_after=3))
    lines.extend(
        brace_segments(
            RIGHT_SHIFT,
            y_shift,
            [(1, 2), (4, 5)],
            r"$D_bf_{3},\; D_b\in B_{2,2}(d)$",
        )
    )
    return lines


def tikz_source() -> str:
    lines: list[str] = []
    lines.extend(brauer_row_lines())
    lines.extend(walled_row_lines())
    body = "\n".join(lines)

    return rf"""\documentclass{{article}}
\usepackage[paperwidth=45cm,paperheight=28cm,margin=0pt]{{geometry}}
\usepackage{{amsmath}}
\usepackage{{graphicx}}
\usepackage{{tikz}}
\usetikzlibrary{{arrows.meta,decorations.pathreplacing}}
\definecolor{{blocka}}{{RGB}}{{88, 165, 242}}
\definecolor{{blockb}}{{RGB}}{{255, 155, 88}}
\definecolor{{blockc}}{{RGB}}{{122, 212, 147}}
\definecolor{{blockd}}{{RGB}}{{186, 147, 245}}
\definecolor{{blocke}}{{RGB}}{{247, 206, 63}}
\definecolor{{blockf}}{{RGB}}{{255, 98, 57}}
\pagestyle{{empty}}

\begin{{document}}
\thispagestyle{{empty}}
\vspace*{{\fill}}
\begin{{center}}
\begin{{tikzpicture}}[
    vertex/.style={{
        circle,
        draw=black,
        line width=1.3pt,
        minimum size=1.56cm,
        inner sep=0pt
    }},
    componentedge/.style={{line width=5.0pt, line cap=round}},
    walldivider/.style={{black!65, dashed, line width=0.9pt}},
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
        else Path(__file__).with_name("no_prop_factorization_rules.jpg")
    )
    output.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as temp_dir:
        workdir = Path(temp_dir)
        tex_path = workdir / "no_prop_factorization_rules.tex"
        tex_path.write_text(tikz_source(), encoding="utf-8")
        pdf_path = build_pdf(workdir, tex_path)
        convert_pdf_to_jpg(pdf_path, output)

    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
