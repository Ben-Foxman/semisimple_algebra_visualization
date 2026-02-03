# Partition Algebra / Brauer Algebra

Diagram algebras with symbolic parameter `d`, implemented in Python using SymPy.

## Requirements

- **Python 3**
- **SymPy:** `pip install sympy`
- **Optional (faster dual basis):** SageMath with `sagelib` available to Python.
  Set `REPTHEORY_USE_SAGE=1` to enable Sage-backed matrix inversion.

## Command line

```text
python main.py <algebra> [args] [--print-matrix]
```

- **algebra:** `partition`, `brauer`, `walled_brauer`, or `symmetric`
- **partition/brauer args:** `k` (positive integer; 2k points total)
- **walled_brauer args:** `r s` with `r + s = k`
- **symmetric args:** `k` (permutation degree)

Examples:

```bash
python main.py partition 1
python main.py partition 2 --print-matrix
python main.py brauer 2
python main.py brauer 3 --print-matrix
python main.py walled_brauer 2 1
python main.py symmetric 3 --print-matrix
```

Output: basis elements and dual basis elements, formatted.

## Project layout

- `algebras.py` — `DiagramAlgebra` (abstract), `PartitionAlgebra`, `BrauerAlgebra`; `d` is a SymPy symbol.
- `main.py` — CLI: algebra name + k, prints basis and dual basis.
- `irreps.py` — (for future use)
- `cache/` — saved computed data (basis, multiplication table, Gram matrix, dual basis)
