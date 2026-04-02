"""Generate a JPG illustrating three walled Brauer factorization rules."""

from __future__ import annotations

import subprocess
import sys
import tempfile
from itertools import permutations
from pathlib import Path

from crop_centered_whitespace import crop_bounds, image_size

LEFT_SHIFT = 1.0
RIGHT_SHIFT = 28.0
ARROW_LEFT_X = 19.6
ARROW_RIGHT_X = 25.4
ARROW_TEXT_X = 22.5
COL_STEP = 2.75
WALL_GAP = 2.20
TOP_ROW_OFFSET = 3.2
BOTTOM_ROW_OFFSET = 0.0
TOP_LABEL_OFFSET = 4.15
BOTTOM_LABEL_OFFSET = -1.08
ROW_TITLE_OFFSET = 7.0
ROW_CAPTION_OFFSET = 5.95
BRACE_OFFSET = 5.30
ARROW_Y_OFFSET = 1.65
PERMUTATION_Y_OFFSET = 3.35


def x_positions(count: int, r: int) -> list[float]:
    xs: list[float] = []
    for idx in range(1, count + 1):
        x = (idx - 1) * COL_STEP
        if idx > r:
            x += WALL_GAP
        xs.append(x)
    return xs


def top_label(col: int) -> str:
    return rf"$\displaystyle {col}$"


def bottom_label(col: int) -> str:
    return rf"$\displaystyle {col}'$"


def wall_x(x_shift: float, count: int, r: int) -> float:
    xs = x_positions(count, r)
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
    r: int,
    top_colors: list[str],
    bottom_colors: list[str],
) -> list[str]:
    xs = x_positions(len(top_colors), r)
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
    r: int,
    top_colors: list[str],
    bottom_colors: list[str],
    caption: str | None = None,
) -> list[str]:
    lines: list[str] = []
    xs = x_positions(len(top_colors), r)
    lines.extend(component_lines(x_shift, y_shift, r, top_colors, bottom_colors))

    lines.append(
        rf"\draw[walldivider] ({wall_x(x_shift, len(top_colors), r):.2f}, {y_shift - 1.6:.2f}) -- "
        rf"({wall_x(x_shift, len(top_colors), r):.2f}, {y_shift + 4.9:.2f});"
    )

    for idx, (x, color) in enumerate(zip(xs, top_colors), start=1):
        lines.append(
            rf"\node[vertexmask] at ({x_shift + x:.2f}, {y_shift + TOP_ROW_OFFSET:.2f}) {{}};"
        )
        lines.append(
            rf"\node[vertex, fill={color}] at ({x_shift + x:.2f}, {y_shift + TOP_ROW_OFFSET:.2f}) {{}};"
        )
        lines.append(
            rf"\node[labelfont, anchor=south] at ({x_shift + x:.2f}, {y_shift + TOP_LABEL_OFFSET:.2f}) "
            rf"{{{top_label(idx)}}};"
        )

    for idx, (x, color) in enumerate(zip(xs, bottom_colors), start=1):
        lines.append(
            rf"\node[vertexmask] at ({x_shift + x:.2f}, {y_shift + BOTTOM_ROW_OFFSET:.2f}) {{}};"
        )
        lines.append(
            rf"\node[vertex, fill={color}] at ({x_shift + x:.2f}, {y_shift + BOTTOM_ROW_OFFSET:.2f}) {{}};"
        )
        lines.append(
            rf"\node[labelfont, anchor=north] at ({x_shift + x:.2f}, {y_shift + BOTTOM_LABEL_OFFSET:.2f}) "
            rf"{{{bottom_label(idx)}}};"
        )

    if caption is not None:
        center_x = x_shift + (xs[0] + xs[-1]) / 2
        lines.append(
            rf"\node[captionfont] at ({center_x:.2f}, {y_shift + ROW_CAPTION_OFFSET:.2f}) "
            rf"{{\scalebox{{1.18}}{{{caption}}}}};"
        )

    return lines


def brace_segment(
    x_shift: float, y_shift: float, r: int, first_col: int, last_col: int
) -> str:
    xs = x_positions(last_col, r)
    left_x = x_shift + xs[first_col - 1] - 0.95
    right_x = x_shift + xs[last_col - 1] + 0.95
    y = y_shift + BRACE_OFFSET
    return rf"\draw[brace] ({left_x:.2f}, {y:.2f}) -- ({right_x:.2f}, {y:.2f});"


def brace_lines(
    x_shift: float,
    y_shift: float,
    r: int,
    spans: list[tuple[int, int]],
    label: str,
) -> list[str]:
    xs = x_positions(max(b for _, b in spans), r)
    label_x = x_shift + (xs[spans[0][0] - 1] + xs[spans[-1][1] - 1]) / 2
    lines = [brace_segment(x_shift, y_shift, r, a, b) for a, b in spans]
    lines.append(
        rf"\node[bracefont, anchor=south] at ({label_x:.2f}, {y_shift + BRACE_OFFSET + 0.40:.2f}) "
        rf"{{\scalebox{{1.18}}{{{label}}}}};"
    )
    return lines


def row_annotation(x_shift: float, y_shift: float, title: str, arrow_formula: str) -> list[str]:
    return [
        rf"\node[rowtitlefont] at ({x_shift + ARROW_TEXT_X:.2f}, {y_shift + ROW_TITLE_OFFSET:.2f}) {{{title}}};",
        rf"\draw[rulearrow] ({x_shift + ARROW_LEFT_X:.2f}, {y_shift + ARROW_Y_OFFSET:.2f}) -- "
        rf"({x_shift + ARROW_RIGHT_X:.2f}, {y_shift + ARROW_Y_OFFSET:.2f});",
        rf"\node[permutationfont] at ({x_shift + ARROW_TEXT_X:.2f}, {y_shift + PERMUTATION_Y_OFFSET:.2f}) "
        rf"{{\scalebox{{1.24}}{{{arrow_formula}}}}};",
    ]


def row_one_lines(x_shift: float, y_shift: float) -> list[str]:
    lines: list[str] = []
    r = 2

    left_top = ["propa", "propb", "prope", "propf", "propb", "propa"]
    left_bottom = ["propc", "propd", "propd", "propc", "propf", "prope"]
    right_top = ["propa", "propb", "propa", "propb", "prope", "propf"]
    right_bottom = ["propc", "propd", "propc", "propd", "prope", "propf"]

    lines.extend(
        row_annotation(
            x_shift,
            y_shift,
            r"{\scalebox{1.28}{$\textbf{Walled Brauer 1}$}}",
            r"$\begin{aligned}\pi_1&=(3\,6)(4\,5)\\[2pt]\pi_2&=(3'\,4')(5'\,6')\end{aligned}$",
        )
    )
    lines.extend(draw_diagram(x_shift + LEFT_SHIFT, y_shift, r, left_top, left_bottom, r"$D_a\in B_{2,4}(d)$"))
    lines.extend(draw_diagram(x_shift + RIGHT_SHIFT, y_shift, r, right_top, right_bottom))
    lines.extend(brace_lines(x_shift + RIGHT_SHIFT, y_shift, r, [(1, 4)], r"$D_b\in B_{2,2}(d)$"))
    return lines


def row_two_lines(x_shift: float, y_shift: float) -> list[str]:
    lines: list[str] = []
    r = 2
    centered_left_shift = x_shift + LEFT_SHIFT + COL_STEP / 2
    centered_right_shift = x_shift + RIGHT_SHIFT + COL_STEP / 2

    left_top = ["propa", "propb", "prope", "propd", "propb"]
    left_bottom = ["propa", "propc", "prope", "propd", "propc"]
    right_top = ["propa", "propb", "propb", "propd", "prope"]
    right_bottom = ["propa", "propc", "propc", "propd", "prope"]

    lines.extend(
        row_annotation(
            x_shift,
            y_shift,
            r"{\scalebox{1.28}{$\textbf{Walled Brauer 2}$}}",
            r"$\begin{aligned}\pi_1&=(3\,5)\\[2pt]\pi_2&=(3'\,5')\end{aligned}$",
        )
    )
    lines.extend(
        draw_diagram(
            centered_left_shift,
            y_shift,
            r,
            left_top,
            left_bottom,
            r"$D_a\in B_{2,3}(d)$",
        )
    )
    lines.extend(draw_diagram(centered_right_shift, y_shift, r, right_top, right_bottom))
    lines.extend(brace_lines(centered_right_shift, y_shift, r, [(1, 4)], r"$D_b\in B_{2,2}(d)$"))
    return lines


def row_three_lines(x_shift: float, y_shift: float) -> list[str]:
    lines: list[str] = []
    r = 3

    left_top = ["prope", "propb", "propa", "propf", "propd", "propb"]
    left_bottom = ["propa", "prope", "propc", "propf", "propd", "propc"]
    right_top = ["propa", "propb", "prope", "propb", "propd", "propf"]
    right_bottom = ["propa", "propc", "prope", "propc", "propd", "propf"]

    lines.extend(
        row_annotation(
            x_shift,
            y_shift,
            r"{\scalebox{1.28}{$\textbf{Walled Brauer 3}$}}",
            r"$\begin{aligned}\pi_1&=(1\,3)(4\,6)\\[2pt]\pi_2&=(2'\,3')(4'\,6')\end{aligned}$",
        )
    )
    lines.extend(draw_diagram(x_shift + LEFT_SHIFT, y_shift, r, left_top, left_bottom, r"$D_a\in B_{3,3}(d)$"))
    lines.extend(draw_diagram(x_shift + RIGHT_SHIFT, y_shift, r, right_top, right_bottom))
    lines.extend(brace_lines(x_shift + RIGHT_SHIFT, y_shift, r, [(1, 2), (4, 5)], r"$D_b\in B_{2,2}(d)$"))
    return lines


def tikz_source(body_lines: list[str], page_width_cm: float, page_height_cm: float) -> str:
    body = "\n".join(body_lines)
    return rf"""\documentclass{{article}}
\usepackage[paperwidth={page_width_cm}cm,paperheight={page_height_cm}cm,margin=0pt]{{geometry}}
\usepackage{{amsmath}}
\usepackage{{graphicx}}
\usepackage{{tikz}}
\usetikzlibrary{{arrows.meta,decorations.pathreplacing}}
\definecolor{{propa}}{{RGB}}{{88, 165, 242}}
\definecolor{{propb}}{{RGB}}{{255, 155, 88}}
\definecolor{{propc}}{{RGB}}{{122, 212, 147}}
\definecolor{{propd}}{{RGB}}{{186, 147, 245}}
\definecolor{{prope}}{{RGB}}{{255, 211, 77}}
\definecolor{{propf}}{{RGB}}{{255, 98, 57}}
\definecolor{{focusa}}{{RGB}}{{255, 98, 57}}
\definecolor{{focusb}}{{RGB}}{{247, 206, 63}}
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
        minimum size=1.56cm,
        inner sep=0pt
    }},
    vertex/.style={{
        circle,
        draw=black,
        line width=1.3pt,
        minimum size=1.56cm,
        inner sep=0pt
    }},
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


def rules_12_source() -> str:
    lines: list[str] = []
    lines.extend(row_one_lines(0.0, 14.0))
    lines.extend(row_two_lines(0.0, 0.0))
    return tikz_source(lines, page_width_cm=48.0, page_height_cm=28.0)


def rule_3_source() -> str:
    return tikz_source(row_three_lines(0.0, 0.0), page_width_cm=48.0, page_height_cm=18.0)


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


def convert_pdf_to_png(pdf_path: Path, png_path: Path) -> None:
    run_command(
        [
            "gs",
            "-dSAFER",
            "-dBATCH",
            "-dNOPAUSE",
            "-sDEVICE=png16m",
            "-r300",
            f"-sOutputFile={png_path}",
            str(pdf_path),
        ],
        cwd=pdf_path.parent,
    )


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


def render_image(tex_source: str, output: Path, stem: str) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as temp_dir:
        workdir = Path(temp_dir)
        tex_path = workdir / f"{stem}.tex"
        tex_path.write_text(tex_source, encoding="utf-8")
        pdf_path = build_pdf(workdir, tex_path)
        png_path = workdir / f"{stem}.png"
        convert_pdf_to_jpg(pdf_path, output)
        convert_pdf_to_png(pdf_path, png_path)
        crop_file_to_bounds(output, crop_bounds(png_path))


def default_outputs() -> tuple[Path, Path]:
    root = Path(__file__).with_name("walled_brauer_factorization_rules")
    return (
        root.with_name("walled_brauer_factorization_rules_12.jpg"),
        root.with_name("walled_brauer_factorization_rules_3.jpg"),
    )


def resolve_outputs(argv: list[str]) -> tuple[Path, Path]:
    if len(argv) == 1:
        return default_outputs()
    if len(argv) == 3:
        return Path(argv[1]).resolve(), Path(argv[2]).resolve()
    raise SystemExit("Usage: generate_walled_brauer_factorization_rules.py [output_12.jpg output_3.jpg]")


def main() -> int:
    output_12, output_3 = resolve_outputs(sys.argv)
    render_image(
        rules_12_source(),
        output_12,
        "walled_brauer_factorization_rules_12",
    )
    render_image(
        rule_3_source(),
        output_3,
        "walled_brauer_factorization_rules_3",
    )

    print(output_12)
    print(output_3)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
