"""Crop outer whitespace tightly from JPG assets in the images directory."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np


THRESHOLD = 245


def image_size(path: Path) -> tuple[int, int]:
    result = subprocess.run(
        [
            "ffprobe",
            "-v",
            "error",
            "-select_streams",
            "v:0",
            "-show_entries",
            "stream=width,height",
            "-of",
            "json",
            str(path),
        ],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    payload = json.loads(result.stdout)
    stream = payload["streams"][0]
    return int(stream["width"]), int(stream["height"])


def read_rgb(path: Path, width: int, height: int) -> np.ndarray:
    result = subprocess.run(
        [
            "ffmpeg",
            "-v",
            "error",
            "-i",
            str(path),
            "-frames:v",
            "1",
            "-f",
            "rawvideo",
            "-pix_fmt",
            "rgb24",
            "-",
        ],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    data = np.frombuffer(result.stdout, dtype=np.uint8)
    return data.reshape((height, width, 3))


def read_rgba(path: Path, width: int, height: int) -> np.ndarray:
    result = subprocess.run(
        [
            "ffmpeg",
            "-v",
            "error",
            "-i",
            str(path),
            "-frames:v",
            "1",
            "-f",
            "rawvideo",
            "-pix_fmt",
            "rgba",
            "-",
        ],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    data = np.frombuffer(result.stdout, dtype=np.uint8)
    return data.reshape((height, width, 4))


def crop_bounds(path: Path) -> tuple[int, int, int, int]:
    width, height = image_size(path)
    image = read_rgb(path, width, height)
    mask = np.any(image < THRESHOLD, axis=2)

    if not np.any(mask):
        return (0, 0, width, height)

    rows = np.any(mask, axis=1)
    cols = np.any(mask, axis=0)

    top = int(np.argmax(rows))
    bottom = int(height - 1 - np.argmax(rows[::-1]))
    left = int(np.argmax(cols))
    right = int(width - 1 - np.argmax(cols[::-1]))

    return (left, top, right + 1, bottom + 1)


def alpha_crop_bounds(path: Path, threshold: int = 0) -> tuple[int, int, int, int]:
    width, height = image_size(path)
    image = read_rgba(path, width, height)
    mask = image[:, :, 3] > threshold

    if not np.any(mask):
        return (0, 0, width, height)

    rows = np.any(mask, axis=1)
    cols = np.any(mask, axis=0)

    top = int(np.argmax(rows))
    bottom = int(height - 1 - np.argmax(rows[::-1]))
    left = int(np.argmax(cols))
    right = int(width - 1 - np.argmax(cols[::-1]))

    return (left, top, right + 1, bottom + 1)


def crop_file(path: Path) -> None:
    left, top, right, bottom = crop_bounds(path)
    width, height = image_size(path)
    if left == 0 and top == 0 and right == width and bottom == height:
        print(f"{path.name}: unchanged")
        return

    out_w = right - left
    out_h = bottom - top

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir) / path.name
        subprocess.run(
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
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        temp_path.replace(path)

    print(
        f"{path.name}: left={left}, top={top}, right={width - right}, bottom={height - bottom}"
    )


def main() -> int:
    root = Path(__file__).resolve().parent
    paths = sorted(root.glob("*.jpg"))
    if not paths:
        return 0
    for path in paths:
        crop_file(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
