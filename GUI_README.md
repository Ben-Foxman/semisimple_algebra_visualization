# Diagram Algebra Viewer

GUI for visualizing diagram algebras (Partition, Brauer, Walled Brauer, Symmetric group): basis diagrams, Gram matrix, dual basis, and irreducible representations.

## Run

```bash
python gui.py
```

## Requirements

- **Python 3**
- **SymPy:** `pip install sympy`
- **PyQt5 or PyQt6:** `pip install PyQt6` (or `PyQt5`)

## Usage

1. Select an **algebra** (Partition, Brauer, Walled Brauer, Symmetric).
2. Enter **parameters** (e.g., `k` for partition/brauer; `r`, `s` for walled Brauer).
3. Click **Load** to compute and display:
   - Basis diagrams (click to inspect)
   - Gram matrix
   - Dual basis elements
   - Irrep decompositions and matrices
4. Use **symbolic d** for formulas in `d`; uncheck for numeric evaluation.
