"""Generate a JPG illustrating block diagonalization by the group Fourier transform."""

from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path


MATRIX_DIM = 20
MATRIX_ANCHOR_X = 31.0
MATRIX_ANCHOR_Y = 17.2
CROP_RIGHT_PX = 17175
CROP_BOTTOM_PX = 906

BLOCKS = [
    (1, 1, 3, 3, r"$\displaystyle \rho_1(g)$", "blocka", 0.92),
    (4, 4, 7, 7, r"$\displaystyle \rho_2(g)$", "blockb", 1.00),
    (8, 8, 11, 11, r"$\displaystyle \rho_2(g)$", "blockb", 1.00),
    (17, 17, 20, 20, r"$\displaystyle \rho_n(g)$", "blockd", 1.00),
]


def matrix_rows() -> str:
    rows: list[str] = []
    for _ in range(MATRIX_DIM):
        rows.append(" & ".join("" for _ in range(MATRIX_DIM)) + r" \\")
    return "\n".join(rows)


def block_lines() -> list[str]:
    lines: list[str] = []
    for idx, (r1, c1, r2, c2, label, fill, label_scale) in enumerate(BLOCKS, start=1):
        lines.append(
            rf"\node[blockbox, fill={fill}, fit=(M-{r1}-{c1})(M-{r2}-{c2})] (B{idx}) {{}};"
        )
        lines.append(
            rf"\node[blockfont, align=center, inner sep=0pt] at (B{idx}.center) {{\scalebox{{{label_scale:.2f}}}{{{label}}}}};"
        )
    return lines


def bracket_lines() -> list[str]:
    return [
        r"\node[fit=(M-1-1)(M-20-20), inner sep=0pt] (MatrixBox) {};",
        r"\draw[black, line width=1.8pt, line cap=round] ($(MatrixBox.north west)+(-0.34,0.10)$) .. controls ($(MatrixBox.north west)+(-1.18,-3.10)$) and ($(MatrixBox.south west)+(-1.18,3.10)$) .. ($(MatrixBox.south west)+(-0.34,-0.10)$);",
        r"\draw[black, line width=1.8pt, line cap=round] ($(MatrixBox.north east)+(0.34,0.10)$) .. controls ($(MatrixBox.north east)+(1.18,-3.10)$) and ($(MatrixBox.south east)+(1.18,3.10)$) .. ($(MatrixBox.south east)+(0.34,-0.10)$);",
    ]


def annotation_lines() -> list[str]:
    lines = [
        r"\node[eqfont, anchor=east] (EqNode) at ($(MatrixBox.west)+(-3.2,0)$) {\scalebox{1.8}{$\displaystyle =$}};",
        r"\node[lhsfont, anchor=east] (LHSNode) at ($(EqNode.west)+(-3.6,0)$) {\scalebox{2.10}{$\displaystyle \mathtt{FT}_G\,U_g\,\mathtt{FT}_G^\dagger$}};",
        r"\coordinate (CaptionBase) at ($(MatrixBox.south)+(0,-2.9)$);",
        r"\node[ketfont, anchor=north] at ($(LHSNode.center |- CaptionBase)$) {\scalebox{1.9}{$\displaystyle \left\{\ket{g}\right\}$}};",
        r"\node[ketfont, anchor=north] at ($(MatrixBox.center |- CaptionBase)$) {\scalebox{1.9}{$\displaystyle \left\{\ket{\rho,i,j}\right\}$}};",
    ]
    for k in range(12, 17):
        lines.append(rf"\node[dotfont] at (M-{k}-{k}.center) {{$\cdot$}};")
    return lines


def tikz_source() -> str:
    rows = matrix_rows()
    lines = [
        rf"\matrix (M) [matrix of nodes, nodes in empty cells, anchor=north west, column sep=0pt, row sep=0pt, nodes={{matrixcell}}] at ({MATRIX_ANCHOR_X:.2f}, {MATRIX_ANCHOR_Y:.2f}) {{",
        rows,
        r"};",
    ]
    lines.extend(block_lines())
    lines.extend(bracket_lines())
    lines.extend(annotation_lines())
    nodes = "\n".join(lines)

    return rf"""\documentclass{{article}}
\usepackage[paperwidth=76cm,paperheight=25cm,margin=0.55cm]{{geometry}}
\usepackage{{amsmath}}
\usepackage{{graphicx}}
\usepackage{{tikz}}
\usetikzlibrary{{matrix,positioning,fit,calc,backgrounds}}
\definecolor{{blocka}}{{RGB}}{{59, 129, 239}}
\definecolor{{blockb}}{{RGB}}{{255, 128, 58}}
\definecolor{{blockd}}{{RGB}}{{245, 213, 110}}
\pagestyle{{empty}}

\newcommand{{\ket}}[1]{{\lvert #1 \rangle}}

\begin{{document}}
\thispagestyle{{empty}}
\noindent
\begin{{tikzpicture}}[
    scale=1.16,
    transform shape,
    matrixcell/.style={{
        minimum width=0.82cm,
        minimum height=0.82cm,
        inner sep=0pt,
        outer sep=0pt
    }},
    blockbox/.style={{
        draw=black,
        line width=0.9pt,
        rounded corners=0pt,
        fill opacity=0.22,
        inner sep=0pt
    }},
    lhsfont/.style={{
        font=\fontsize{{260}}{{268}}\selectfont
    }},
    eqfont/.style={{
        font=\fontsize{{300}}{{308}}\selectfont
    }},
    blockfont/.style={{
        font=\fontsize{{48}}{{52}}\selectfont
    }},
    ketfont/.style={{
        font=\fontsize{{230}}{{236}}\selectfont
    }},
    dotfont/.style={{
        font=\fontsize{{74}}{{74}}\selectfont
    }}
]
{nodes}
\end{{tikzpicture}}
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
            "-dJPEGQ=100",
            "-r1800",
            f"-sOutputFile={jpg_path}",
            str(pdf_path),
        ],
        cwd=pdf_path.parent,
    )


def jpg_size(jpg_path: Path) -> tuple[int, int]:
    result = subprocess.run(
        ["sips", "-g", "pixelWidth", "-g", "pixelHeight", str(jpg_path)],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    width: int | None = None
    height: int | None = None
    for line in result.stdout.splitlines():
        stripped = line.strip()
        if stripped.startswith("pixelWidth:"):
            width = int(stripped.split(":", 1)[1].strip())
        elif stripped.startswith("pixelHeight:"):
            height = int(stripped.split(":", 1)[1].strip())
    if width is None or height is None:
        raise RuntimeError("Could not determine JPG size with sips")
    return width, height


def crop_jpg(jpg_path: Path) -> None:
    width, height = jpg_size(jpg_path)
    out_width = width - CROP_RIGHT_PX
    out_height = height - CROP_BOTTOM_PX
    if out_width <= 0 or out_height <= 0:
        raise RuntimeError("Crop would produce a non-positive image size")

    offset_y = -(CROP_BOTTOM_PX // 2)
    offset_x = -(CROP_RIGHT_PX // 2)

    with tempfile.TemporaryDirectory() as temp_dir:
        cropped = Path(temp_dir) / jpg_path.name
        subprocess.run(
            [
                "sips",
                "-c",
                str(out_height),
                str(out_width),
                "--cropOffset",
                str(offset_y),
                str(offset_x),
                str(jpg_path),
                "--out",
                str(cropped),
            ],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        cropped.replace(jpg_path)


def main() -> int:
    output = (
        Path(sys.argv[1]).resolve()
        if len(sys.argv) > 1
        else Path(__file__).with_name("group_ft_block_diagonalization.jpg")
    )
    output.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as temp_dir:
        workdir = Path(temp_dir)
        tex_path = workdir / "group_ft_block_diagonalization.tex"
        tex_path.write_text(tikz_source(), encoding="utf-8")
        pdf_path = build_pdf(workdir, tex_path)
        convert_pdf_to_jpg(pdf_path, output)
        crop_jpg(output)

    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
