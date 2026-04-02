"""Generate a JPG illustrating four partition no-propagation factorization rules."""

from __future__ import annotations

import math
import subprocess
import sys
import tempfile
from itertools import permutations
from pathlib import Path

from crop_centered_whitespace import image_size


PANEL_XS = [2.0, 35.0]
PANEL_YS = [17.6, 3.6]
PAGE_WIDTH_CM = 68.0
PAGE_HEIGHT_CM = 32.0
RENDER_DPI = 300
LEFT_LOCAL_X = 0.0
RIGHT_LOCAL_X = 19.0
ARROW_LEFT_LOCAL_X = 13.0
ARROW_RIGHT_LOCAL_X = 17.0
ARROW_TEXT_LOCAL_X = 15.0
PANEL_TITLE_LOCAL_X = 15.0
COL_STEP = 2.75
TOP_ROW_OFFSET = 3.2
BOTTOM_ROW_OFFSET = 0.0
TOP_LABEL_OFFSET = 4.15
BOTTOM_LABEL_OFFSET = -1.08
TITLE_OFFSET = 8.35
CAPTION_OFFSET = 5.95
BRACE_OFFSET = 5.30
ARROW_Y_OFFSET = 1.65
PERMUTATION_Y_OFFSET = 3.35
DIVIDER_X = 33.50
DIVIDER_BOTTOM_Y = 1.20
DIVIDER_TOP_Y = PANEL_YS[0] + TITLE_OFFSET


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
            rf"\node[captionfont] at ({center_x:.2f}, {y_shift + CAPTION_OFFSET:.2f}) "
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
        rf"{{\scalebox{{1.10}}{{{label}}}}};",
    ]


def panel_annotation(x_shift: float, y_shift: float, title: str, arrow_formula: str | None) -> list[str]:
    lines = [
        rf"\node[rowtitlefont] at ({x_shift + PANEL_TITLE_LOCAL_X:.2f}, {y_shift + TITLE_OFFSET:.2f}) {{{title}}};",
        rf"\draw[rulearrow] ({x_shift + ARROW_LEFT_LOCAL_X:.2f}, {y_shift + ARROW_Y_OFFSET:.2f}) -- "
        rf"({x_shift + ARROW_RIGHT_LOCAL_X:.2f}, {y_shift + ARROW_Y_OFFSET:.2f});",
    ]
    if arrow_formula:
        lines.append(
            rf"\node[permutationfont] at ({x_shift + ARROW_TEXT_LOCAL_X:.2f}, {y_shift + PERMUTATION_Y_OFFSET:.2f}) "
            rf"{{\scalebox{{1.20}}{{{arrow_formula}}}}};"
        )
    return lines


def panel_one_lines() -> list[str]:
    x_shift = PANEL_XS[0]
    y_shift = PANEL_YS[0]

    right_top = ["blocka", "blocka", "blockb", "blockb", "blockc"]
    right_bottom = ["blockd", "blocke", "blocke", "blockf", "blockg"]

    lines = panel_annotation(
        x_shift,
        y_shift,
        r"{\scalebox{1.18}{$\textbf{Partition No Propagation 1}$}}",
        None,
    )
    lines.extend(
        draw_diagram(
            x_shift + LEFT_LOCAL_X,
            y_shift,
            right_top,
            right_bottom,
            r"$D_a\in P_{5}(d)$",
        )
    )
    lines.extend(draw_diagram(x_shift + RIGHT_LOCAL_X, y_shift, right_top, right_bottom))
    lines.extend(
        brace_lines(
            x_shift + RIGHT_LOCAL_X,
            y_shift,
            1,
            4,
            r"$D_bp_{5},\; D_b\in P_{4}(d)$",
        )
    )
    return lines


def panel_two_lines() -> list[str]:
    x_shift = PANEL_XS[1]
    y_shift = PANEL_YS[0]

    right_top = ["blocka", "blocka", "blocka", "blockb", "blockc"]
    right_bottom = ["blockd", "blockd", "blockd", "blocke", "blocke"]
    left_top = right_top
    left_bottom = ["blockd", "blocke", "blockd", "blockd", "blocke"]

    lines = panel_annotation(
        x_shift,
        y_shift,
        r"{\scalebox{1.18}{$\textbf{Partition No Propagation 2}$}}",
        r"$\pi=(2'\,4')$",
    )
    lines.extend(
        draw_diagram(
            x_shift + LEFT_LOCAL_X,
            y_shift,
            left_top,
            left_bottom,
            r"$D_a\in P_{5}(d)$",
        )
    )
    lines.extend(draw_diagram(x_shift + RIGHT_LOCAL_X, y_shift, right_top, right_bottom))
    lines.extend(
        brace_lines(
            x_shift + RIGHT_LOCAL_X,
            y_shift,
            1,
            4,
            r"$D_bp_{5}b_{4},\; D_b\in P_{4}(d)$",
        )
    )
    return lines


def panel_three_lines() -> list[str]:
    x_shift = PANEL_XS[0]
    y_shift = PANEL_YS[1]

    right_top = ["blocka", "blocka", "blocka", "blockb", "blockb"]
    right_bottom = ["blockd", "blockd", "blockd", "blocke", "blockf"]
    left_top = ["blocka", "blockb", "blocka", "blocka", "blockb"]
    left_bottom = right_bottom

    lines = panel_annotation(
        x_shift,
        y_shift,
        r"{\scalebox{1.18}{$\textbf{Partition No Propagation 3}$}}",
        r"$\pi=(2\,4)$",
    )
    lines.extend(
        draw_diagram(
            x_shift + LEFT_LOCAL_X,
            y_shift,
            left_top,
            left_bottom,
            r"$D_a\in P_{5}(d)$",
        )
    )
    lines.extend(draw_diagram(x_shift + RIGHT_LOCAL_X, y_shift, right_top, right_bottom))
    lines.extend(
        brace_lines(
            x_shift + RIGHT_LOCAL_X,
            y_shift,
            1,
            4,
            r"$b_{4}D_bp_{5},\; D_b\in P_{4}(d)$",
        )
    )
    return lines


def panel_four_lines() -> list[str]:
    x_shift = PANEL_XS[1]
    y_shift = PANEL_YS[1]

    right_top = ["blocka", "blocka", "blockb", "blockb", "blockb"]
    right_bottom = ["blockd", "blockd", "blocke", "blocke", "blocke"]
    left_top = ["blocka", "blockb", "blockb", "blocka", "blockb"]
    left_bottom = ["blockd", "blocke", "blocke", "blockd", "blocke"]

    lines = panel_annotation(
        x_shift,
        y_shift,
        r"{\scalebox{1.18}{$\textbf{Partition No Propagation 4}$}}",
        r"$\begin{aligned}\pi_1&=(1\,4)\\[2pt]\pi_2&=(2'\,4')\end{aligned}$",
    )
    lines.extend(
        draw_diagram(
            x_shift + LEFT_LOCAL_X,
            y_shift,
            left_top,
            left_bottom,
            r"$D_a\in P_{5}(d)$",
        )
    )
    lines.extend(draw_diagram(x_shift + RIGHT_LOCAL_X, y_shift, right_top, right_bottom))
    lines.extend(
        brace_lines(
            x_shift + RIGHT_LOCAL_X,
            y_shift,
            1,
            4,
            r"$b_{4}D_bp_{5}b_{4},\; D_b\in P_{4}(d)$",
        )
    )
    return lines


def tikz_source() -> str:
    lines: list[str] = []
    lines.extend(panel_one_lines())
    lines.extend(panel_two_lines())
    lines.extend(panel_three_lines())
    lines.extend(panel_four_lines())
    lines.append(
        rf"\draw[dividerline] ({DIVIDER_X:.2f}, {DIVIDER_BOTTOM_Y:.2f}) -- ({DIVIDER_X:.2f}, {DIVIDER_TOP_Y:.2f});"
    )
    body = "\n".join(lines)

    return rf"""\documentclass{{article}}
\usepackage[paperwidth={PAGE_WIDTH_CM}cm,paperheight={PAGE_HEIGHT_CM}cm,margin=0pt]{{geometry}}
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
\definecolor{{blockg}}{{RGB}}{{116, 203, 232}}
\definecolor{{blockh}}{{RGB}}{{255, 214, 102}}
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
    rowtitlefont/.style={{font=\fontsize{{24}}{{24}}\selectfont}},
    permutationfont/.style={{font=\fontsize{{21}}{{21}}\selectfont}},
    rulearrow/.style={{->, >=Latex, line width=1.2pt}},
    bracefont/.style={{font=\fontsize{{21}}{{21}}\selectfont}},
    brace/.style={{decorate, decoration={{brace, amplitude=8pt}}, line width=1.25pt}},
    dividerline/.style={{black!60, densely dotted, line width=1.0pt}}
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
            f"-r{RENDER_DPI}",
            f"-sOutputFile={jpg_path}",
            str(pdf_path),
        ],
        cwd=pdf_path.parent,
    )


def pdf_content_bbox(pdf_path: Path) -> tuple[float, float, float, float]:
    result = subprocess.run(
        ["gs", "-q", "-dBATCH", "-dNOPAUSE", "-sDEVICE=bbox", str(pdf_path)],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    for line in result.stderr.splitlines():
        if line.startswith("%%HiResBoundingBox:"):
            _, left, bottom, right, top = line.split()
            return float(left), float(bottom), float(right), float(top)
    raise RuntimeError("Ghostscript did not report a bounding box.")


def bbox_to_pixel_bounds(
    bbox: tuple[float, float, float, float],
    image_width: int,
    image_height: int,
) -> tuple[int, int, int, int]:
    left_pt, bottom_pt, right_pt, top_pt = bbox
    scale = RENDER_DPI / 72.0
    page_height_pt = PAGE_HEIGHT_CM * 72.0 / 2.54

    left = max(0, math.ceil(left_pt * scale))
    right = min(image_width, math.floor(right_pt * scale))
    top = max(0, math.ceil((page_height_pt - top_pt) * scale))
    bottom = min(image_height, math.floor((page_height_pt - bottom_pt) * scale))
    return (left, top, right, bottom)


def crop_file_to_bounds(path: Path, bounds: tuple[int, int, int, int]) -> None:
    left, top, right, bottom = bounds
    width, height = image_size(path)
    if left == 0 and top == 0 and right == width and bottom == height:
        print(f"{path.name}: unchanged")
        return

    out_w = right - left
    out_h = bottom - top

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir) / path.name
        run_command(
            [
                "ffmpeg",
                "-v",
                "error",
                "-y",
                "-i",
                str(path),
                "-vf",
                f"crop={out_w}:{out_h}:{left}:{top}",
                str(temp_path),
            ],
            cwd=path.parent,
        )
        temp_path.replace(path)

    print(
        f"{path.name}: left={left}, top={top}, right={width - right}, bottom={height - bottom}"
    )


def main() -> int:
    output = (
        Path(sys.argv[1]).resolve()
        if len(sys.argv) > 1
        else Path(__file__).with_name("partition_no_prop_factorization_rules.jpg")
    )
    output.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as temp_dir:
        workdir = Path(temp_dir)
        tex_path = workdir / "partition_no_prop_factorization_rules.tex"
        tex_path.write_text(tikz_source(), encoding="utf-8")
        pdf_path = build_pdf(workdir, tex_path)
        convert_pdf_to_jpg(pdf_path, output)
        crop_file_to_bounds(
            output,
            bbox_to_pixel_bounds(pdf_content_bbox(pdf_path), *image_size(output)),
        )

    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
