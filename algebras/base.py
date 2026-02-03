"""
Diagram algebras (Partition algebra, Brauer algebra, etc.) using SymPy.

Requires: pip install sympy

- d is a symbolic variable (SymPy Symbol), not a numeric parameter.
- DiagramAlgebra: abstract base; subclasses implement _generate_basis() (and optionally
  _mult_table(), _is_zero_result() for Brauer).
- Multiplication table: (i,j) -> (basis_index, coefficient in Q[d]); loops contribute d^n.
- Dual basis: Gram matrix G[i,j] = tr(R(e_i * e_j)); dual e_i^* has coefficients (G^{-1})_ji in Q(d).
"""

from abc import ABC, abstractmethod
from collections import defaultdict, deque
import itertools
import json
import os

import sympy

try:
    from sage.all import SR as _sage_SR
    from sage.all import matrix as _sage_matrix

    _SAGE_AVAILABLE = True
except Exception:
    _SAGE_AVAILABLE = False


def _sage_expr_from_sympy(expr):
    return _sage_SR(str(expr))


def _sympy_expr_from_sage(expr, d_symbol):
    text = str(expr).replace("^", "**")
    return sympy.sympify(text, locals={"d": d_symbol})


def _sage_inverse_matrix_sympy(sympy_matrix, d_symbol):
    n = sympy_matrix.rows
    rows = [
        [_sage_expr_from_sympy(sympy_matrix[i, j]) for j in range(n)]
        for i in range(n)
    ]
    sage_matrix = _sage_matrix(_sage_SR, rows)
    sage_inv = sage_matrix.inverse()
    return [
        [_sympy_expr_from_sage(sage_inv[i, j], d_symbol) for j in range(n)]
        for i in range(n)
    ]


def _sage_enabled():
    if not _SAGE_AVAILABLE:
        return False
    flag = os.getenv("REPTHEORY_USE_SAGE", "").strip().lower()
    return flag in {"1", "true", "yes", "y", "on"}


def _set_partition(blocks):
    """Return a normalized tuple of frozensets for use as partition key."""
    return tuple(sorted(frozenset(sorted(b)) for b in blocks if b))


def _partition_size(part):
    return sum(part)


def _addable_partitions(part):
    if not part:
        return [((1,))]
    out = []
    for i in range(len(part)):
        if i == 0 or part[i - 1] > part[i]:
            new = list(part)
            new[i] += 1
            out.append(tuple(new))
    out.append(tuple(list(part) + [1]))
    return out


def _removable_partitions(part):
    out = []
    for i in range(len(part)):
        if i == len(part) - 1 or part[i] > part[i + 1]:
            new = list(part)
            new[i] -= 1
            if new[i] == 0:
                new = new[:i]
            out.append(tuple(new))
    return out


def _partitions_of(n):
    if n == 0:
        return [tuple()]
    out = []

    def rec(remaining, max_part, current):
        if remaining == 0:
            out.append(tuple(current))
            return
        for k in range(min(remaining, max_part), 0, -1):
            rec(remaining - k, k, current + [k])

    rec(n, n, [])
    return out


def _addable_nodes(part):
    nodes = []
    if not part:
        return [(1, 1)]
    for i in range(len(part)):
        if i == 0 or part[i - 1] > part[i]:
            nodes.append((i + 1, part[i] + 1))
    nodes.append((len(part) + 1, 1))
    return nodes


def _removable_nodes(part):
    nodes = []
    for i in range(len(part)):
        if i == len(part) - 1 or part[i] > part[i + 1]:
            nodes.append((i + 1, part[i]))
    return nodes


def _content(cell):
    r, c = cell
    return c - r


def _add_cell(part, cell):
    r, c = cell
    rows = list(part)
    while len(rows) < r:
        rows.append(0)
    if rows[r - 1] + 1 != c:
        return None
    rows[r - 1] += 1
    return tuple(x for x in rows if x > 0)


def _remove_cell(part, cell):
    r, c = cell
    rows = list(part)
    if r > len(rows) or rows[r - 1] != c:
        return None
    rows[r - 1] -= 1
    if rows[r - 1] == 0:
        rows = rows[: r - 1]
    return tuple(rows)


def _diff_cell(prev, curr):
    if _partition_size(curr) == _partition_size(prev) + 1:
        # added
        for i in range(max(len(curr), len(prev))):
            prev_len = prev[i] if i < len(prev) else 0
            curr_len = curr[i] if i < len(curr) else 0
            if curr_len == prev_len + 1:
                return ("add", (i + 1, curr_len))
        return None
    if _partition_size(curr) + 1 == _partition_size(prev):
        for i in range(max(len(curr), len(prev))):
            prev_len = prev[i] if i < len(prev) else 0
            curr_len = curr[i] if i < len(curr) else 0
            if prev_len == curr_len + 1:
                return ("rem", (i + 1, prev_len))
        return None
    return None


def _dominance_geq(a, b):
    max_len = max(len(a), len(b))
    sa = sb = 0
    for i in range(max_len):
        sa += a[i] if i < len(a) else 0
        sb += b[i] if i < len(b) else 0
        if sa < sb:
            return False
    return True


def _path_compare(path_a, path_b):
    # reverse lexicographic on (lambda,l)
    for idx in range(len(path_a) - 1, -1, -1):
        la = path_a[idx]
        lb = path_b[idx]
        if la != lb:
            return la, lb
    return None


def _partition_label(part):
    if part is None:
        return "()"
    if isinstance(part, tuple):
        return "(" + ",".join(str(x) for x in part) + ")"
    return str(part)


class SetPartition:
    """Set partition of {1,...,k} ∪ {-1,...,-k} (bottom row)."""

    def __init__(self, blocks):
        blocks = [frozenset(sorted(b)) for b in blocks if b]
        self._blocks = tuple(sorted(blocks, key=lambda b: (len(b), min(b))))

    @property
    def blocks(self):
        return list(self._blocks)

    def __eq__(self, other):
        return isinstance(other, SetPartition) and self._blocks == other._blocks

    def __hash__(self):
        return hash(self._blocks)

    def __repr__(self):
        return f"SetPartition({list(self._blocks)})"

    def __str__(self):
        return "|".join("{" + ",".join(str(x) for x in sorted(b)) + "}" for b in self._blocks)


class DiagramAlgebra(ABC):
    """
    Abstract base class for a diagram algebra with symbolic parameter d.
    Basis: diagram basis. Multiplication: diagram composition; loops contribute factor d.
    Dual basis: defined by inverting the Gram matrix of the trace form tr(e_i * e_j).
    """

    algebra_id = "diagram"
    cache_version = "v19"

    def __init__(self, k, d=5, symbolic_d=False):
        self.k = k
        self.d = self._normalize_d(d, symbolic_d)

        self._multiplication_table = None
        self._gram_matrix_cache = None
        self._dual_basis_cache = None
        self._basis_traces = None
        self._label_list_cache = None
        self._irreps_cache = None
        self._bratteli_cache = None
        self._irrep_matrices_cache = None
        self._matrix_units_cache = None
        self._matrix_unit_labels_cache = None

        self.basis = []
        self.dim = 0
        self._basis_lookup = {}

        if not self._load_cache():
            self.basis = self._generate_basis()
            self.dim = len(self.basis)
            self._basis_lookup = {self._key(p): i for i, p in enumerate(self.basis)}
            self._save_cache()

    @abstractmethod
    def _generate_basis(self):
        """Return list of SetPartition forming the basis."""
        pass

    def _key(self, elem):
        """Hashable key for basis lookup."""
        return _set_partition(elem.blocks)

    def _serialize_basis_element(self, elem):
        return [sorted(list(b)) for b in elem.blocks]

    def _deserialize_basis_element(self, data):
        return SetPartition([set(b) for b in data])

    def _serialize_coeff(self, coeff):
        return sympy.srepr(coeff)

    def _deserialize_coeff(self, text):
        return sympy.sympify(text)

    def _serialize_matrix(self, mat):
        return [[self._serialize_coeff(c) for c in row] for row in mat.tolist()]

    def _deserialize_matrix(self, data):
        rows = [[self._deserialize_coeff(c) for c in row] for row in data]
        return sympy.Matrix(rows)

    def _serialize_shape(self, shape):
        if (
            isinstance(shape, tuple)
            and len(shape) == 2
            and all(isinstance(x, tuple) for x in shape)
        ):
            return [list(shape[0]), list(shape[1])]
        return list(shape)

    def _deserialize_shape(self, data):
        if (
            isinstance(data, list)
            and len(data) == 2
            and all(isinstance(x, list) for x in data)
        ):
            return (tuple(data[0]), tuple(data[1]))
        return tuple(data)

    def _cache_params(self):
        return {"k": self.k, "d": self._d_cache_token()}

    def _normalize_d(self, d, symbolic_d):
        if symbolic_d or d is None:
            self.is_symbolic_d = True
            return sympy.Symbol("d")
        if isinstance(d, sympy.Symbol):
            self.is_symbolic_d = True
            return d
        self.is_symbolic_d = False
        return sympy.sympify(d)

    def _d_cache_token(self):
        if getattr(self, "is_symbolic_d", False):
            return "sym"
        text = str(self.d)
        return text.replace("/", "div").replace(" ", "")

    def _identity_element(self):
        return SetPartition([{i, -i} for i in range(1, self.k + 1)])

    def _star(self, elem):
        blocks = [{-x for x in b} for b in elem.blocks]
        return SetPartition(blocks)

    def _diagram_s(self, i):
        blocks = [{i, -(i + 1)}, {i + 1, -i}]
        for j in range(1, self.k + 1):
            if j in (i, i + 1):
                continue
            blocks.append({j, -j})
        return SetPartition(blocks)

    def _diagram_e(self, i):
        blocks = [{i, i + 1}, {-i, -(i + 1)}]
        for j in range(1, self.k + 1):
            if j in (i, i + 1):
                continue
            blocks.append({j, -j})
        return SetPartition(blocks)

    def _diagram_p(self, i):
        blocks = [{i}, {-i}]
        for j in range(1, self.k + 1):
            if j == i:
                continue
            blocks.append({j, -j})
        return SetPartition(blocks)

    def _generator_diagrams(self):
        raise NotImplementedError

    def _cache_path(self):
        base_dir = os.path.dirname(os.path.dirname(__file__))
        cache_dir = os.path.join(base_dir, "cache")
        params = self._cache_params()
        parts = [f"{k}{params[k]}" for k in sorted(params.keys())]
        filename = f"{self.algebra_id}_{'_'.join(parts)}.json"
        return cache_dir, os.path.join(cache_dir, filename)

    def _load_cache(self):
        cache_dir, path = self._cache_path()
        if not os.path.exists(path):
            return False
        try:
            with open(path, "r", encoding="utf-8") as fh:
                data = json.load(fh)
        except (OSError, json.JSONDecodeError):
            return False

        if data.get("version") != self.cache_version:
            return False
        if data.get("algebra") != self.algebra_id:
            return False
        if data.get("params") != self._cache_params():
            return False

        basis_data = data.get("basis")
        if not isinstance(basis_data, list):
            return False
        self.basis = [self._deserialize_basis_element(b) for b in basis_data]
        self.dim = len(self.basis)
        self._basis_lookup = {self._key(p): i for i, p in enumerate(self.basis)}

        mult_data = data.get("multiplication_table")
        if isinstance(mult_data, dict):
            table = {}
            for key, value in mult_data.items():
                i_str, j_str = key.split(",")
                idx, coeff = value
                table[(int(i_str), int(j_str))] = (int(idx), self._deserialize_coeff(coeff))
            self._multiplication_table = table

        gram_data = data.get("gram_matrix")
        if isinstance(gram_data, list):
            rows = [[self._deserialize_coeff(c) for c in row] for row in gram_data]
            self._gram_matrix_cache = sympy.Matrix(rows)

        dual_data = data.get("dual_basis")
        if isinstance(dual_data, list):
            dual = []
            for row in dual_data:
                entry = {int(j): self._deserialize_coeff(c) for j, c in row.items()}
                dual.append(entry)
            self._dual_basis_cache = dual

        units_data = data.get("matrix_units")
        if isinstance(units_data, list):
            units = []
            for row in units_data:
                entry = {int(j): self._deserialize_coeff(c) for j, c in row.items()}
                units.append(entry)
            self._matrix_units_cache = units

        unit_labels = data.get("matrix_unit_labels")
        if isinstance(unit_labels, list):
            self._matrix_unit_labels_cache = unit_labels

        labels_data = data.get("labels")
        if isinstance(labels_data, list):
            if self.dim == len(labels_data):
                self._label_list_cache = labels_data
            else:
                self._label_list_cache = None

        irreps_data = data.get("irreps")
        paths_data = data.get("bratteli_paths")
        if isinstance(irreps_data, list) and isinstance(paths_data, list):
            self._irreps_cache = [self._deserialize_shape(x) for x in irreps_data]
            self._bratteli_cache = [
                [[self._deserialize_shape(step) for step in path] for path in paths]
                for paths in paths_data
            ]

        irrep_mats = data.get("irrep_matrices")
        if isinstance(irrep_mats, list):
            self._irrep_matrices_cache = []
            for entry in irrep_mats:
                if isinstance(entry, dict):
                    mats = {
                        k: self._deserialize_matrix(v) for k, v in entry.items()
                    }
                    self._irrep_matrices_cache.append(mats)

        return True

    def _save_cache(self):
        cache_dir, path = self._cache_path()
        os.makedirs(cache_dir, exist_ok=True)

        payload = {
            "version": self.cache_version,
            "algebra": self.algebra_id,
            "params": self._cache_params(),
            "basis": [self._serialize_basis_element(b) for b in self.basis],
        }

        if self._multiplication_table is not None:
            mult = {}
            for (i, j), (idx, coeff) in self._multiplication_table.items():
                mult[f"{i},{j}"] = [idx, self._serialize_coeff(coeff)]
            payload["multiplication_table"] = mult

        if self._gram_matrix_cache is not None:
            payload["gram_matrix"] = [
                [self._serialize_coeff(c) for c in row]
                for row in self._gram_matrix_cache.tolist()
            ]

        if self._dual_basis_cache is not None:
            dual = []
            for row in self._dual_basis_cache:
                dual.append({str(j): self._serialize_coeff(c) for j, c in row.items()})
            payload["dual_basis"] = dual

        if self._matrix_units_cache is not None:
            units = []
            for row in self._matrix_units_cache:
                units.append({str(j): self._serialize_coeff(c) for j, c in row.items()})
            payload["matrix_units"] = units

        if self._matrix_unit_labels_cache is not None:
            payload["matrix_unit_labels"] = self._matrix_unit_labels_cache

        if self._label_list_cache is not None:
            payload["labels"] = self._label_list_cache

        if self._irreps_cache is not None and self._bratteli_cache is not None:
            payload["irreps"] = [self._serialize_shape(x) for x in self._irreps_cache]
            payload["bratteli_paths"] = [
                [[self._serialize_shape(step) for step in path] for path in paths]
                for paths in self._bratteli_cache
            ]

        if getattr(self, "_irrep_matrices_cache", None) is not None:
            payload["irrep_matrices"] = [
                {k: self._serialize_matrix(v) for k, v in entry.items()}
                for entry in self._irrep_matrices_cache
            ]

        try:
            with open(path, "w", encoding="utf-8") as fh:
                json.dump(payload, fh, indent=2, sort_keys=True)
        except OSError:
            pass

    def _compose(self, p1, p2):
        """
        Compose two diagrams: p1 * p2 (bottom of p1 with top of p2).
        Returns (result_partition, number_of_loops).
        """
        adj = defaultdict(list)

        def add_edge(u, v):
            adj[u].append(v)
            adj[v].append(u)

        def connect_block(layer, block):
            nodes = [(layer, x) for x in block]
            for n in nodes:
                adj.setdefault(n, [])
            if len(nodes) <= 1:
                return
            head = nodes[0]
            for n in nodes[1:]:
                add_edge(head, n)

        for b in p1.blocks:
            connect_block("p1", b)
        for b in p2.blocks:
            connect_block("p2", b)
        for i in range(1, self.k + 1):
            add_edge(("p1", -i), ("p2", i))

        visited = set()
        blocks = []
        loops = 0

        def is_boundary(node):
            side, val = node
            return (side == "p1" and val > 0) or (side == "p2" and val < 0)

        def boundary_label(node):
            return node[1]

        for start in adj:
            if start in visited:
                continue
            stack = [start]
            comp = set()
            while stack:
                v = stack.pop()
                if v in visited:
                    continue
                visited.add(v)
                comp.add(v)
                for w in adj[v]:
                    if w not in visited:
                        stack.append(w)

            boundary = [node for node in comp if is_boundary(node)]
            if len(boundary) == 0:
                loops += 1
            else:
                labels = {boundary_label(n) for n in boundary}
                blocks.append(labels)

        return SetPartition(blocks), loops

    def _trace_basis_index(self, j):
        """
        Right regular trace of basis element e_j:
        tr(R(e_j)) = sum_i (coefficient of e_i in e_i * e_j).
        """
        total = sympy.S.Zero
        for i in range(self.dim):
            idx, coeff = self.multiplication_table[(i, j)]
            if idx == i:
                total += coeff
        return total

    def _mult_table(self):
        """Multiplication table: (i,j) -> (basis_index, coefficient in R)."""
        table = {}
        for i, p1 in enumerate(self.basis):
            for j, p2 in enumerate(self.basis):
                res, nloops = self._compose(p1, p2)
                coeff = self.d ** nloops
                key = self._key(res)
                idx = self._basis_lookup[key]
                table[(i, j)] = (idx, coeff)
        return table

    @property
    def multiplication_table(self):
        if self._multiplication_table is None:
            self._multiplication_table = self._mult_table()
            self._save_cache()
        return self._multiplication_table

    def _gram_matrix(self):
        """Gram matrix G[i,j] = tr(R(e_i * e_j)); entries are SymPy expressions in d."""
        if self._basis_traces is None:
            self._basis_traces = [self._trace_basis_index(j) for j in range(self.dim)]
        n = self.dim
        rows = []
        for i in range(n):
            row = []
            for j in range(n):
                idx, coeff = self.multiplication_table.get((i, j), (0, sympy.S.Zero))
                # e_i * e_j = coeff * e_idx, so tr(e_i * e_j) = coeff * tr(e_idx)
                row.append(coeff * self._basis_traces[idx])
            rows.append(row)
        return sympy.Matrix(rows)

    @property
    def gram_matrix(self):
        if self._gram_matrix_cache is None:
            self._gram_matrix_cache = self._gram_matrix()
            self._save_cache()
        return self._gram_matrix_cache

    def _dual_basis(self):
        """
        Dual basis with respect to trace form: tr(e_i * e_j^*) = δ_ij.
        e_i^* = sum_j (G^{-1})_ji e_j; coefficients are SymPy expressions in d.
        """
        G = self.gram_matrix
        Gi_entries = None
        if _sage_enabled():
            try:
                Gi_entries = _sage_inverse_matrix_sympy(G, self.d)
            except Exception:
                Gi_entries = None
        if Gi_entries is None:
            try:
                Gi = G.inv()
                Gi_entries = [[Gi[i, j] for j in range(G.rows)] for i in range(G.rows)]
            except (ValueError, ZeroDivisionError):
                return [{j: sympy.S.One if j == i else sympy.S.Zero for j in range(self.dim)}
                        for i in range(self.dim)]
        n = self.dim
        dual = []
        for i in range(n):
            row = {}
            for j in range(n):
                c = Gi_entries[j][i]
                if c != 0:
                    row[j] = c
            dual.append(row)
        return dual

    @property
    def dual_basis(self):
        if self._dual_basis_cache is None:
            self._dual_basis_cache = self._dual_basis()
            self._save_cache()
        return self._dual_basis_cache

    def _compute_matrix_units(self):
        # E^rho_{ij} = d_rho * sum_a rho(a^*)_{ij} * a, where a^* is dual basis.
        basis = self.basis
        basis_labels = [self.label_of(a) for a in basis]
        dual = self.dual_basis
        units = []
        for ir_idx, paths in enumerate(self.bratteli_paths):
            d_rho = len(paths)
            rho_basis = [
                self.irrep_matrix_for_label(ir_idx, lbl) for lbl in basis_labels
            ]
            irrep_units = [dict() for _ in range(d_rho * d_rho)]
            for a_idx, coeffs in enumerate(dual):
                rho_a_star = sympy.zeros(d_rho, d_rho)
                for j_idx, coeff in coeffs.items():
                    rho_a_star += coeff * rho_basis[j_idx]
                for i in range(d_rho):
                    for j in range(d_rho):
                        coeff = sympy.S(d_rho) * rho_a_star[j, i]
                        if coeff != 0:
                            irrep_units[i * d_rho + j][a_idx] = coeff
            units.extend(irrep_units)
        return units

    def _compute_matrix_unit_labels(self):
        labels = []
        for irrep, paths in zip(self.irreps, self.bratteli_paths):
            if isinstance(irrep, tuple) and len(irrep) == 2 and all(
                isinstance(x, tuple) for x in irrep
            ):
                ir_label = f"{_partition_label(irrep[0])} | {_partition_label(irrep[1])}"
            else:
                ir_label = _partition_label(irrep)
            d_rho = len(paths)
            for i in range(d_rho):
                for j in range(d_rho):
                    labels.append(f"E^{ir_label}_{{{i+1},{j+1}}}")
        return labels

    @property
    def matrix_units(self):
        if self._matrix_units_cache is None:
            self._matrix_units_cache = self._compute_matrix_units()
            self._save_cache()
        return self._matrix_units_cache

    @property
    def matrix_unit_labels(self):
        if self._matrix_unit_labels_cache is None:
            self._matrix_unit_labels_cache = self._compute_matrix_unit_labels()
            self._save_cache()
        return self._matrix_unit_labels_cache

    def _compute_irreps_and_paths(self):
        raise NotImplementedError

    def _compute_irrep_matrices(self):
        raise NotImplementedError

    @property
    def irreps(self):
        if self._irreps_cache is None or self._bratteli_cache is None:
            self._irreps_cache, self._bratteli_cache = self._compute_irreps_and_paths()
            self._save_cache()
        return self._irreps_cache

    @property
    def bratteli_paths(self):
        if self._irreps_cache is None or self._bratteli_cache is None:
            self._irreps_cache, self._bratteli_cache = self._compute_irreps_and_paths()
            self._save_cache()
        return self._bratteli_cache

    @property
    def irrep_matrices(self):
        if getattr(self, "_irrep_matrices_cache", None) is None:
            self._irrep_matrices_cache = self._compute_irrep_matrices()
            self._save_cache()
        return self._irrep_matrices_cache

    def irrep_matrix_for_label(self, irrep_idx, label):
        raise NotImplementedError

    def _compute_labels(self):
        generators = self._generator_diagrams()
        if not generators:
            if self.dim == 1:
                return ["1"]
            return ["" for _ in range(self.dim)]
        start = self._identity_element()
        labels = ["" for _ in range(self.dim)]
        start_idx = self._basis_lookup[self._key(start)]
        labels[start_idx] = "1"
        queue = deque([start])

        while queue:
            current = queue.popleft()
            current_idx = self._basis_lookup[self._key(current)]
            current_label = labels[current_idx]
            for name, gen in generators:
                res, _ = self._compose(current, gen)
                res_idx = self._basis_lookup[self._key(res)]
                if labels[res_idx]:
                    continue
                labels[res_idx] = name if current_label == "1" else f"{current_label} {name}"
                queue.append(res)
                if all(labels):
                    break
            if all(labels):
                break
        return labels

    def label_of(self, elem):
        if self._label_list_cache is None or len(self._label_list_cache) != self.dim:
            self._label_list_cache = self._compute_labels()
            self._save_cache()
        idx = self._basis_lookup[self._key(elem)]
        return self._label_list_cache[idx] or str(elem)

    def multiply(self, coeffs1, coeffs2):
        """Multiply two elements given as dicts: index -> coefficient (SymPy expression)."""
        result = defaultdict(lambda: sympy.S.Zero)
        for i, c1 in coeffs1.items():
            for j, c2 in coeffs2.items():
                idx, coeff = self.multiplication_table.get((i, j), (0, sympy.S.Zero))
                result[idx] += c1 * c2 * coeff
        return dict(result)

    def __repr__(self):
        return f"{self.__class__.__name__}(k={self.k}, dim={self.dim})"

__all__ = [name for name in globals() if not name.startswith('__')]
