from .base import *

class WalledBrauerAlgebra(DiagramAlgebra):
    """
    Walled Brauer algebra with parameters r, s (k = r + s).
    Pairs must respect the wall: through-strings stay on each side;
    horizontal pairs cross the wall within top or within bottom.
    """

    algebra_id = "walled_brauer"

    def __init__(self, r, s, d=5, symbolic_d=False):
        self.r = r
        self.s = s
        super().__init__(r + s, d=d, symbolic_d=symbolic_d)

    def _cache_params(self):
        return {"r": self.r, "s": self.s, "d": self._d_cache_token()}

    def _allowed_pair(self, a, b):
        side_a = "L" if abs(a) <= self.r else "R"
        side_b = "L" if abs(b) <= self.r else "R"
        if a > 0 and b > 0:
            return side_a != side_b
        if a < 0 and b < 0:
            return side_a != side_b
        return side_a == side_b

    def _generate_basis(self):
        elts = list(range(1, self.k + 1)) + [-i for i in range(1, self.k + 1)]
        return self._pair_partitions(elts)

    def _generator_diagrams(self):
        gens = []
        for i in range(1, self.k):
            if i == self.r:
                continue
            gens.append((f"s_{i}", self._diagram_s(i)))
        if 1 <= self.r < self.k:
            gens.append((f"e_{self.r}", self._diagram_e(self.r)))
        return gens

    def _compute_irreps_and_paths(self):
        r = self.r
        s = self.s

        all_paths = defaultdict(list)

        # Walled Brauer tableaux:
        # first r steps: add a box to left partition
        # next s steps: add a box to right partition OR remove a box from left partition
        def rec_left(level, left_part, path):
            if level == r:
                rec_right(0, left_part, tuple(), path)
                return
            for nxt in _addable_partitions(left_part):
                rec_left(level + 1, nxt, path + [(nxt, tuple())])

        def rec_right(level, left_part, right_part, path):
            if level == s:
                all_paths[(left_part, right_part)].append(path)
                return
            for nxt in _addable_partitions(right_part):
                rec_right(level + 1, left_part, nxt, path + [(left_part, nxt)])
            for nxt in _removable_partitions(left_part):
                rec_right(level + 1, nxt, right_part, path + [(nxt, right_part)])

        rec_left(0, tuple(), [(tuple(), tuple())])

        irreps = []
        paths_by_irrep = []
        for t in range(0, min(r, s) + 1):
            for left in _partitions_of(r - t):
                for right in _partitions_of(s - t):
                    key = (left, right)
                    if key in all_paths:
                        irreps.append(key)
                        paths_by_irrep.append(all_paths[key])
        return irreps, paths_by_irrep

    def _wb_step_info(self, path):
        info = {}
        for i in range(1, len(path)):
            prev_l, prev_r = path[i - 1]
            curr_l, curr_r = path[i]
            if i <= self.r:
                diff = _diff_cell(prev_l, curr_l)
                if diff is None or diff[0] != "add":
                    raise ValueError("Invalid left-add step in walled path.")
                info[i] = ("L_add", diff[1])
            else:
                if prev_l == curr_l:
                    diff = _diff_cell(prev_r, curr_r)
                    if diff is None or diff[0] != "add":
                        raise ValueError("Invalid right-add step in walled path.")
                    info[i] = ("R_add", diff[1])
                elif prev_r == curr_r:
                    diff = _diff_cell(prev_l, curr_l)
                    if diff is None or diff[0] != "rem":
                        raise ValueError("Invalid left-rem step in walled path.")
                    info[i] = ("L_rem", diff[1])
                else:
                    raise ValueError("Invalid walled step.")
        return info

    def _wb_wcont(self, path):
        info = self._wb_step_info(path)
        wcont = {}
        for i, (kind, cell) in info.items():
            if i <= self.r:
                wcont[i] = _content(cell)
            else:
                if kind == "R_add":
                    wcont[i] = _content(cell) + self.d
                else:
                    wcont[i] = -_content(cell)
        return wcont

    def _wb_swap_path(self, path, i):
        info = self._wb_step_info(path)
        if i not in info or i + 1 not in info:
            return None
        info2 = dict(info)
        info2[i], info2[i + 1] = info2[i + 1], info2[i]
        left = tuple()
        right = tuple()
        new_path = [(left, right)]
        for step in range(1, len(path)):
            kind, cell = info2[step]
            if step <= self.r:
                if kind != "L_add":
                    return None
                new_left = _add_cell(left, cell)
                if new_left is None:
                    return None
                left = new_left
            else:
                if kind == "R_add":
                    new_right = _add_cell(right, cell)
                    if new_right is None:
                        return None
                    right = new_right
                elif kind == "L_rem":
                    new_left = _remove_cell(left, cell)
                    if new_left is None:
                        return None
                    left = new_left
                else:
                    return None
            new_path.append((left, right))
        return tuple(new_path)

    def _wb_mobile_paths(self, path):
        if self.r <= 0 or self.r >= len(path):
            return []
        if path[self.r - 1] != path[self.r + 1]:
            return []
        prev_l, prev_r = path[self.r - 1]
        if prev_r != tuple():
            return []
        paths = []
        for cell in _addable_nodes(prev_l):
            new_l = _add_cell(prev_l, cell)
            if new_l is None:
                continue
            new_path = list(path)
            new_path[self.r] = (new_l, tuple())
            paths.append(tuple(new_path))
        return paths

    def _wb_c_coeff(self, path):
        if path[self.r - 1] != path[self.r + 1]:
            return sympy.S.Zero
        prev_l, _ = path[self.r - 1]
        cur_l, _ = path[self.r]
        diff = _diff_cell(prev_l, cur_l)
        if diff is None or diff[0] != "add":
            return sympy.S.Zero
        a = diff[1]
        num = sympy.S.One
        den = sympy.S.One
        for c in _removable_nodes(prev_l):
            num *= (_content(a) - _content(c))
        for b in _addable_nodes(prev_l):
            if b != a:
                den *= (_content(a) - _content(b))
        return sympy.sqrt((self.d + _content(a)) * num / den)

    def _compute_irrep_matrices(self):
        irreps = self.irreps
        paths = self.bratteli_paths
        matrices = []
        for ir_idx, irrep in enumerate(irreps):
            basis_paths = [tuple(p) for p in paths[ir_idx]]
            index = {p: idx for idx, p in enumerate(basis_paths)}
            n = len(basis_paths)
            gen_mats = {}
            for i in range(1, self.k):
                if i == self.r:
                    M = sympy.zeros(n, n)
                    for a, T in enumerate(basis_paths):
                        cT = self._wb_c_coeff(T)
                        if cT == 0:
                            continue
                        mobile = self._wb_mobile_paths(T)
                        for T2 in mobile:
                            cT2 = self._wb_c_coeff(T2)
                            if cT2 == 0:
                                continue
                            b = index[T2]
                            M[b, a] += cT * cT2
                    gen_mats[f"e_{i}"] = M
                else:
                    M = sympy.zeros(n, n)
                    for a, T in enumerate(basis_paths):
                        wcont = self._wb_wcont(T)
                        r_i = wcont[i + 1] - wcont[i]
                        if r_i == 0:
                            continue
                        M[a, a] += sympy.S.One / r_i
                        T2 = self._wb_swap_path(T, i)
                        if T2 is not None:
                            b = index[T2]
                            M[b, a] += sympy.sqrt(1 - sympy.S.One / (r_i ** 2))
                    gen_mats[f"s_{i}"] = M
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
                raise ValueError(f"Unknown generator {tok} for walled Brauer.")
            if result is None:
                result = mats[tok]
            else:
                result = result * mats[tok]
        return result

    def _pair_partitions(self, elts):
        if not elts:
            return [SetPartition([])]
        if len(elts) < 2:
            return []
        first, rest = elts[0], elts[1:]
        out = []
        for i, other in enumerate(rest):
            if not self._allowed_pair(first, other):
                continue
            pair = {first, other}
            rem = rest[:i] + rest[i + 1:]
            for p in self._pair_partitions(rem):
                out.append(SetPartition(p.blocks + [pair]))
        return out
