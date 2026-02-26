from .base import *

class BrauerAlgebra(DiagramAlgebra):
    """Brauer algebra B_k(d): pair partitions (all blocks of size 2)."""

    algebra_id = "brauer"

    def _generator_diagrams(self):
        gens = []
        for i in range(1, self.k):
            gens.append((f"s_{i}", self._diagram_s(i)))
        for i in range(1, self.k):
            gens.append((f"e_{i}", self._diagram_e(i)))
        return gens

    def _compute_irreps_and_paths(self):
        k = self.k
        all_paths = defaultdict(list)

        def rec(level, part, path):
            if level == k:
                all_paths[part].append(path[:])
                return
            for nxt in _addable_partitions(part):
                rec(level + 1, nxt, path + [nxt])
            for nxt in _removable_partitions(part):
                rec(level + 1, nxt, path + [nxt])

        rec(0, tuple(), [tuple()])
        irreps = []
        paths_by_irrep = []
        for r in range(0, k // 2 + 1):
            for part in _partitions_of(k - 2 * r):
                if part in all_paths:
                    irreps.append(part)
                    paths_by_irrep.append(all_paths[part])
        for paths in paths_by_irrep:
            paths.sort(key=_path_lex_key)
        return irreps, paths_by_irrep

    def _brauer_contents(self, path):
        contents = {}
        for i in range(1, len(path)):
            prev = path[i - 1]
            curr = path[i]
            diff = _diff_cell(prev, curr)
            if diff is None:
                raise ValueError("Invalid Brauer path step.")
            kind, cell = diff
            base = (self.d - 1) / 2 + _content(cell)
            contents[i] = base if kind == "add" else -base
        return contents

    def _brauer_i_equivalent(self, s, t, i):
        for j in range(len(s)):
            if j != i and s[j] != t[j]:
                return False
        return True

    def _brauer_mid_candidates(self, prev, nxt):
        mids = []
        for add in _addable_nodes(prev):
            mid = _add_cell(prev, add)
            if mid is None:
                continue
            if _diff_cell(mid, nxt) is not None:
                mids.append(mid)
        for rem in _removable_nodes(prev):
            mid = _remove_cell(prev, rem)
            if mid is None:
                continue
            if _diff_cell(mid, nxt) is not None:
                mids.append(mid)
        # unique
        out = []
        for m in mids:
            if m not in out:
                out.append(m)
        return out

    def _brauer_xk(self, path, k):
        prev = path[k - 1]
        curr = path[k]
        diff = _diff_cell(prev, curr)
        if diff is None:
            raise ValueError("Invalid Brauer path step.")
        kind, cell = diff
        base = (self.d - 1) / 2 + _content(cell)
        return base if kind == "add" else -base

    def _brauer_b_set(self, mu):
        B_set = []
        for add in _addable_nodes(mu):
            B_set.append((self.d - 1) / 2 + _content(add))
        for rem in _removable_nodes(mu):
            B_set.append(-((self.d - 1) / 2 + _content(rem)))
        return B_set

    def _brauer_same(self, a, b):
        try:
            return sympy.simplify(a - b) == 0
        except Exception:
            return a == b

    def _brauer_ek_diag(self, mu, b):
        B_set = self._brauer_b_set(mu)
        prod = sympy.S.One
        for bj in B_set:
            if self._brauer_same(bj, b):
                continue
            prod *= (b + bj) / (b - bj)
        diag_general = (2 * b + 1) * prod
        diag_special = -prod
        return sympy.Piecewise(
            (diag_special, sympy.simplify(b + sympy.S.Half) == 0),
            (diag_general, True),
        )

    def _brauer_swap_path(self, t, i):
        prev = t[i - 1]
        curr = t[i]
        nxt = t[i + 1]
        candidates = self._brauer_mid_candidates(prev, nxt)
        if len(candidates) != 2:
            return None
        alt = candidates[0] if candidates[1] == curr else candidates[1]
        t2 = list(t)
        t2[i] = alt
        return tuple(t2)

    def _brauer_ek(self, s_path, t_path, i):
        if not self._brauer_i_equivalent(s_path, t_path, i):
            return sympy.S.Zero
        if t_path[i - 1] != t_path[i + 1]:
            return sympy.S.Zero
        mu = t_path[i - 1]
        if s_path == t_path:
            b = self._brauer_xk(t_path, i)
            return self._brauer_ek_diag(mu, b)
        b_s = self._brauer_xk(s_path, i)
        b_t = self._brauer_xk(t_path, i)
        e_ss = self._brauer_ek_diag(mu, b_s)
        e_tt = self._brauer_ek_diag(mu, b_t)
        return sympy.sqrt(e_ss * e_tt)

    def _brauer_sk(self, s_path, t_path, i):
        if not self._brauer_i_equivalent(s_path, t_path, i):
            return sympy.S.Zero
        xk_t = self._brauer_xk(t_path, i)
        xk1_t = self._brauer_xk(t_path, i + 1)
        if t_path[i - 1] != t_path[i + 1]:
            if s_path == t_path:
                return sympy.S.One / (xk1_t - xk_t)
            t_swap = self._brauer_swap_path(t_path, i)
            if t_swap is None or s_path != t_swap:
                return sympy.S.Zero
            return sympy.sqrt(sympy.S.One - sympy.S.One / ((xk1_t - xk_t) ** 2))
        e_st = self._brauer_ek(s_path, t_path, i)
        delta = sympy.S.One if s_path == t_path else sympy.S.Zero
        denom = xk_t + self._brauer_xk(s_path, i)
        if s_path == t_path and sympy.simplify(denom) == 0:
            return sympy.S.One
        if sympy.simplify(denom) == 0:
            return sympy.S.Zero
        return (e_st - delta) / denom

    def _compute_irrep_matrices(self):
        irreps = self.irreps
        paths = self.bratteli_paths
        matrices = []
        for ir_idx, irrep in enumerate(irreps):
            basis_paths = [tuple(p) for p in paths[ir_idx]]
            n = len(basis_paths)
            gen_mats = {}
            for i in range(1, self.k):
                M_s = sympy.zeros(n, n)
                M_e = sympy.zeros(n, n)
                for a, s_path in enumerate(basis_paths):
                    for b, t_path in enumerate(basis_paths):
                        M_s[a, b] = self._brauer_sk(s_path, t_path, i)
                        M_e[a, b] = self._brauer_ek(s_path, t_path, i)
                gen_mats[f"s_{i}"] = M_s
                gen_mats[f"e_{i}"] = M_e
            matrices.append(gen_mats)
        return matrices

    def irrep_matrix_for_label(self, irrep_idx, label):
        if label == "1" or not label:
            n = len(self.bratteli_paths[irrep_idx])
            return sympy.eye(n)
        mats = self.irrep_matrices[irrep_idx]
        tokens = label.split()
        result = None
        for tok in tokens:
            if tok not in mats:
                raise ValueError(f"Unknown generator {tok} for Brauer algebra.")
            if result is None:
                result = mats[tok]
            else:
                result = result * mats[tok]
        return result
    def _generate_basis(self):
        elts = list(range(1, self.k + 1)) + [-i for i in range(1, self.k + 1)]
        return self._pair_partitions(elts)

    def _pair_partitions(self, elts):
        if not elts:
            return [SetPartition([])]
        if len(elts) < 2:
            return []
        first, rest = elts[0], elts[1:]
        out = []
        for i, other in enumerate(rest):
            pair = {first, other}
            rem = rest[:i] + rest[i + 1:]
            for p in self._pair_partitions(rem):
                out.append(SetPartition(p.blocks + [pair]))
        return out
