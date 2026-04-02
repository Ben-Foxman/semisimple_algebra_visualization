"""Generate a JPG showing G, s_2(G), b_2(G), and p_2(G) on a five-vertex graph."""

from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path


PANEL_SHIFTS = [4.0, 20.5, 37.0, 53.5]
TITLE_Y = -4.85
TITLE_SCALE = 1.75

COLORS = {
    "1": "blocka",
    "2": "blockb",
    "3": "blockc",
    "4": "blockd",
    "5": "blocke",
    "{2,3}": "mergedcolor",
}

BASE_POSITIONS = {
    "1": (4.0, 3.3),
    "2": (1.15, 1.23),
    "3": (6.85, 1.23),
    "4": (2.24, -2.13),
    "5": (5.76, -2.13),
}

ORIGINAL_EDGES = [("1", "2"), ("2", "3"), ("2", "4"), ("3", "4"), ("3", "5")]

PANELS = [
    {
        "title": r"$\displaystyle G$",
        "positions": BASE_POSITIONS,
        "edges": ORIGINAL_EDGES,
        "labels": ["1", "2", "3", "4", "5"],
        "display_labels": {"1": "1", "2": "2", "3": "3", "4": "4", "5": "5"},
        "merged": None,
    },
    {
        "title": r"$\displaystyle s_2(G)$",
        "positions": BASE_POSITIONS,
        "edges": ORIGINAL_EDGES,
        "labels": ["1", "2", "3", "4", "5"],
        "display_labels": {"1": "1", "2": "3", "3": "2", "4": "4", "5": "5"},
        "merged": None,
    },
    {
        "title": r"$\displaystyle b_2(G)$",
        "positions": {
            "1": (4.0, 3.55),
            "{2,3}": (4.0, 0.15),
            "4": (1.80, -2.35),
            "5": (6.20, -2.35),
        },
        "edges": [("1", "{2,3}"), ("{2,3}", "4"), ("{2,3}", "5")],
        "labels": ["1", "{2,3}", "4", "5"],
        "display_labels": {"1": "1", "{2,3}": r"\{2,3\}", "4": "4", "5": "5"},
        "merged": "{2,3}",
    },
    {
        "title": r"$\displaystyle p_2(G)$",
        "positions": {
            "1": BASE_POSITIONS["1"],
            "2": BASE_POSITIONS["2"],
            "3": BASE_POSITIONS["3"],
            "4": BASE_POSITIONS["4"],
            "5": BASE_POSITIONS["5"],
        },
        "edges": [("3", "4"), ("3", "5")],
        "labels": ["1", "2", "3", "4", "5"],
        "display_labels": {"1": "1", "2": "2", "3": "3", "4": "4", "5": "5"},
        "merged": None,
    },
]


def edge_lines(shift: float, panel: dict[str, object]) -> list[str]:
    positions = panel["positions"]
    edges = panel["edges"]
    lines: list[str] = []
    for left, right in edges:
        x_a, y_a = positions[left]
        x_b, y_b = positions[right]
        lines.append(
            rf"\draw[graphedge] ({shift + x_a:.2f}, {y_a:.2f}) -- ({shift + x_b:.2f}, {y_b:.2f});"
        )
    return lines


def vertex_lines(
    shift: float,
    label: str,
    display_label: str,
    x: float,
    y: float,
    merged: bool,
) -> list[str]:
    color = COLORS[label]
    if merged:
        return [
            rf"\node[mergedmask] at ({shift + x:.2f}, {y:.2f}) {{}};",
            rf"\node[mergedvertex, fill={color}] at ({shift + x:.2f}, {y:.2f}) {{}};",
            rf"\node[mergedlabel] at ({shift + x:.2f}, {y:.2f}) {{$\displaystyle {display_label}$}};",
        ]

    return [
        rf"\node[vertexmask] at ({shift + x:.2f}, {y:.2f}) {{}};",
        rf"\node[vertex, fill={color}] at ({shift + x:.2f}, {y:.2f}) {{}};",
        rf"\node[vertexlabel] at ({shift + x:.2f}, {y:.2f}) {{$\displaystyle {display_label}$}};",
    ]


def panel_lines(shift: float, panel: dict[str, object]) -> list[str]:
    lines = edge_lines(shift, panel)
    lines.append(
        rf"\node[titlefont] at ({shift + 4.0:.2f}, {TITLE_Y:.2f}) {{\scalebox{{{TITLE_SCALE}}}{{{panel['title']}}}}};"
    )
    merged_label = panel["merged"]
    display_labels = panel["display_labels"]
    for label in panel["labels"]:
        x, y = panel["positions"][label]
        lines.extend(
            vertex_lines(
                shift,
                label,
                display_labels[label],
                x,
                y,
                label == merged_label,
            )
        )
    return lines


def tikz_source() -> str:
    lines: list[str] = []
    for shift, panel in zip(PANEL_SHIFTS, PANELS):
        lines.extend(panel_lines(shift, panel))
    nodes = "\n".join(lines)

    return rf"""\documentclass{{article}}
\usepackage[paperwidth=67cm,paperheight=15cm,margin=0pt]{{geometry}}
\usepackage{{amsmath}}
\usepackage{{graphicx}}
\usepackage{{tikz}}
\definecolor{{blocka}}{{RGB}}{{59, 129, 239}}
\definecolor{{blockb}}{{RGB}}{{255, 128, 58}}
\definecolor{{blockc}}{{RGB}}{{123, 214, 146}}
\definecolor{{blockd}}{{RGB}}{{255, 192, 214}}
\definecolor{{blocke}}{{RGB}}{{245, 213, 110}}
\definecolor{{mergedcolor}}{{RGB}}{{187, 171, 102}}
\pagestyle{{empty}}

\begin{{document}}
\thispagestyle{{empty}}
\vspace*{{\fill}}
\begin{{center}}
\begin{{tikzpicture}}[
    graphedge/.style={{
        line width=4.4pt,
        line cap=round,
        draw=black!78
    }},
    vertexmask/.style={{
        circle,
        draw=none,
        fill=white,
        minimum size=1.75cm,
        inner sep=0pt
    }},
    vertex/.style={{
        circle,
        draw=black,
        line width=1.25pt,
        minimum size=1.75cm,
        inner sep=0pt
    }},
    mergedmask/.style={{
        circle,
        draw=none,
        fill=white,
        minimum size=2.30cm,
        inner sep=0pt
    }},
    mergedvertex/.style={{
        circle,
        draw=black,
        line width=1.25pt,
        minimum size=2.30cm,
        inner sep=0pt
    }},
    titlefont/.style={{
        font=\fontsize{{48}}{{52}}\selectfont
    }},
    vertexlabel/.style={{
        font=\fontsize{{34}}{{36}}\selectfont
    }},
    mergedlabel/.style={{
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
        else Path(__file__).with_name("graph_s2_b2_p2_example.jpg")
    )
    output.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as temp_dir:
        workdir = Path(temp_dir)
        tex_path = workdir / "graph_s2_b2_p2_example.tex"
        tex_path.write_text(tikz_source(), encoding="utf-8")
        pdf_path = build_pdf(workdir, tex_path)
        convert_pdf_to_jpg(pdf_path, output)

    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
