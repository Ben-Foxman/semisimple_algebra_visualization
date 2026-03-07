# Diagram Algebra Viewer

This repository provides a desktop GUI for exploring several diagram algebras:

- partition algebra
- half-partition algebra
- Brauer algebra
- walled Brauer algebra
- symmetric group algebra

The application is built around a single workflow:

1. choose an algebra and its parameters,
2. load the basis and structural data,
3. inspect matrices, irreducible representations, and Fourier-state data in the GUI.

## Run

```bash
python gui.py
```

## Requirements

- Python 3
- SymPy
- PyQt6 or PyQt5

Example:

```bash
pip install sympy PyQt6
```

## How The GUI Works

After you press `Load`, the GUI constructs an algebra object from `algebras/`, then lazily computes and caches the expensive data attached to that algebra.

The main pieces shown in the interface are:

- `Basis`: the diagram basis of the algebra, with labels generated from the algebra generators.
- `Gram matrix`: the trace-form Gram matrix of the basis.
- `Dual basis`: the basis dual to the trace form.
- `Irrep matrices`: matrix realizations of the generators in each irreducible representation.
- `Matrix units / Fourier states`: matrix-unit coordinates and the normalized Fourier-state vectors derived from them.

Most large computations are cached under `cache/` locally so that reloading the same algebra/parameter choice is much faster on later runs.

## Parameters

- Partition algebra: choose `k` and numeric `d`.
- Half-partition algebra: choose the integer part `k` and numeric `d`.
- Brauer algebra: choose `k` and numeric `d`.
- Walled Brauer algebra: choose `r`, `s`, and numeric `d`.
- Symmetric group algebra: choose `n`.

`d` is numeric only in the current GUI.

## Notes

- The matrix and Fourier-state views can become expensive for larger parameters.
- The GUI displays exact SymPy expressions when practical and numerical approximations when needed for readability.
- Cache files, `__pycache__`, and one-off verification scripts are local artifacts and are not intended to be committed.
