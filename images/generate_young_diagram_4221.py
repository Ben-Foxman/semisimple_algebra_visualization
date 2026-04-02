"""Generate a JPG showing the Young diagram for the partition (4, 2, 2, 1)."""

from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path


PARTITION = [4, 2, 2, 1]
BOX_SIZE = 1.8
LINE_WIDTH = 1.4
FILL_COLOR = "boxfill"


def diagram_lines() -> list[str]:
    lines: list[str] = []
    for row_index, row_length in enumerate(PARTITION):
        y_top = -row_index * BOX_SIZE
        for col_index in range(row_length):
            x_left = col_index * BOX_SIZE
            lines.append(
                rf"\filldraw[fill={FILL_COLOR}, draw=black, line width={LINE_WIDTH}pt] "
                rf"({x_left:.2f}, {y_top:.2f}) rectangle "
                rf"({x_left + BOX_SIZE:.2f}, {y_top - BOX_SIZE:.2f});"
            )
    return lines


def tikz_source() -> str:
    nodes = "\n".join(diagram_lines())
    return rf"""\documentclass{{article}}
\usepackage[paperwidth=12.8cm,paperheight=11.6cm,margin=0pt]{{geometry}}
\usepackage{{tikz}}
\definecolor{{boxfill}}{{RGB}}{{173, 214, 255}}
\pagestyle{{empty}}

\begin{{document}}
\thispagestyle{{empty}}
\vspace*{{\fill}}
\begin{{center}}
\begin{{tikzpicture}}[
    x=1cm,
    y=1cm
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
        else Path(__file__).with_name("young_diagram_4_2_2_1.jpg")
    )
    output.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as temp_dir:
        workdir = Path(temp_dir)
        tex_path = workdir / "young_diagram_4_2_2_1.tex"
        tex_path.write_text(tikz_source(), encoding="utf-8")
        pdf_path = build_pdf(workdir, tex_path)
        convert_pdf_to_jpg(pdf_path, output)

    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
