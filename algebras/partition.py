from .base import *

class PartitionAlgebra(DiagramAlgebra):
    """Partition algebra P_k(d): all set partitions of {1,...,k} ∪ {-1,...,-k}."""

    algebra_id = "partition"

    def _diagram_e(self, i):
        blocks = [{i, i + 1, -i, -(i + 1)}]
        for j in range(1, self.k + 1):
            if j in (i, i + 1):
                continue
            blocks.append({j, -j})
        return SetPartition(blocks)

    def _generator_diagrams(self):
        gens = []
        for i in range(1, self.k):
            gens.append((f"s_{i}", self._diagram_s(i)))
        for i in range(1, self.k):
            gens.append((f"b_{i}", self._diagram_e(i)))
        for i in range(1, self.k + 1):
            gens.append((f"p_{i}", self._diagram_p(i)))
        return gens

    def _compute_irreps_and_paths(self):
        k = self.k
        irreps = []
        paths_by_irrep = []

        # Paths in hat{A} to level 2k, with nodes (lambda, l).
        all_paths = defaultdict(list)

        def node_for(level, part):
            l_val = (level // 2) - _partition_size(part)
            return (part, l_val)

        def valid_node(level, part):
            l_val = (level // 2) - _partition_size(part)
            return l_val >= 0

        def rec(level, part, path):
            if level == 2 * k:
                all_paths[path[-1]].append(path[:])
                return
            if level % 2 == 0:
                # level 2i -> 2i+1: stay or remove
                if valid_node(level + 1, part):
                    rec(level + 1, part, path + [node_for(level + 1, part)])
                for nxt in _removable_partitions(part):
                    if valid_node(level + 1, nxt):
                        rec(level + 1, nxt, path + [node_for(level + 1, nxt)])
            else:
                # level 2i+1 -> 2i+2: stay or add
                if valid_node(level + 1, part):
                    rec(level + 1, part, path + [node_for(level + 1, part)])
                for nxt in _addable_partitions(part):
                    if valid_node(level + 1, nxt):
                        rec(level + 1, nxt, path + [node_for(level + 1, nxt)])

        rec(0, tuple(), [node_for(0, tuple())])
        for l in range(k + 1):
            size = k - l
            for part in _partitions_of(size):
                key = (part, l)
                if key in all_paths:
                    irreps.append(key)
                    paths_by_irrep.append(all_paths[key])
        return irreps, paths_by_irrep

    def _pa_l(self, node):
        return node[1]

    def _pa_content(self, path, i):
        prev = path[i - 1][0]
        curr = path[i][0]
        diff = _diff_cell(prev, curr)
        if i % 2 == 0:
            if diff is None:
                return self.d - _partition_size(curr)
            if diff[0] == "add":
                return _content(diff[1])
        else:
            if diff is None:
                return _partition_size(curr)
            if diff[0] == "rem":
                return self.d - _content(diff[1])
        raise ValueError("Invalid partition path step for content.")

    def _pa_contents(self, path):
        return {i: self._pa_content(path, i) for i in range(1, len(path))}

    def _pa_norm(self, path):
        # Recursive norm from branching factors (Proposition 5.13).
        val = sympy.S.One
        for level in range(1, len(path)):
            prev = path[level - 1]
            curr = path[level]
            l_prev = self._pa_l(prev)
            l_curr = self._pa_l(curr)
            if l_curr == l_prev:
                gamma = self._pa_gamma(prev[0], curr[0], level % 2 == 0)
            else:
                # l increases by 1; use e_{level-1}(t,t) * gamma_{curr -> prev}
                e_diag = self._pa_e_diag(path, level - 1)
                gamma = e_diag * self._pa_gamma(curr[0], prev[0], level % 2 == 0)
            val *= sympy.simplify(gamma)
        return sympy.simplify(val)

    def _pa_i_equiv(self, s, t, i):
        for j in range(len(s)):
            if j != i and s[j] != t[j]:
                return False
        return True

    def _pa_i_approx(self, s, t, i):
        for j in range(len(s)):
            if j != i and j != i - 1 and s[j] != t[j]:
                return False
        return True

    def _pa_next_nodes(self, node, level):
        part, _ = node
        out = []
        if level % 2 == 0:
            # stay or remove
            l_next = (level + 1 - _partition_size(part)) // 2
            if l_next >= 0:
                out.append((part, l_next))
            for r in _removable_nodes(part):
                nxt = _remove_cell(part, r)
                if nxt is not None:
                    l_next = (level + 1 - _partition_size(nxt)) // 2
                    if l_next >= 0:
                        out.append((nxt, l_next))
        else:
            # stay or add
            l_next = (level + 1 - _partition_size(part)) // 2
            if l_next >= 0:
                out.append((part, l_next))
            for a in _addable_nodes(part):
                nxt = _add_cell(part, a)
                if nxt is not None:
                    l_next = (level + 1 - _partition_size(nxt)) // 2
                    if l_next >= 0:
                        out.append((nxt, l_next))
        # unique
        uniq = []
        for n in out:
            if n not in uniq:
                uniq.append(n)
        return uniq

    def _pa_approx_candidates(self, t, k):
        candidates = []
        prev2 = t[k - 2]
        next1 = t[k + 1]
        for mid1 in self._pa_next_nodes(prev2, k - 2):
            for mid2 in self._pa_next_nodes(mid1, k - 1):
                if next1 in self._pa_next_nodes(mid2, k):
                    cand = list(t)
                    cand[k - 1] = mid1
                    cand[k] = mid2
                    candidates.append(tuple(cand))
        # unique
        out = []
        for c in candidates:
            if c not in out:
                out.append(c)
        return out

    def _pa_compare_paths(self, s, t):
        for idx in range(len(s) - 1, -1, -1):
            if s[idx] != t[idx]:
                l_s = self._pa_l(s[idx])
                l_t = self._pa_l(t[idx])
                if l_s != l_t:
                    return 1 if l_s > l_t else -1
                if _dominance_geq(s[idx][0], t[idx][0]) and s[idx] != t[idx]:
                    return 1
                return -1
        return 0

    def _pa_gamma_add(self, lam, add_cell):
        num = sympy.S.One
        den = sympy.S.One
        for beta in _addable_nodes(lam):
            if beta != add_cell and beta[0] > add_cell[0]:
                num *= (_content(add_cell) - _content(beta))
        for beta in _removable_nodes(lam):
            if beta[0] > add_cell[0]:
                den *= (_content(add_cell) - _content(beta))
        return sympy.simplify(num / den)

    def _pa_gamma(self, lam, mu, step_even):
        if lam == mu:
            return sympy.S.One
        diff = _diff_cell(lam, mu)
        if diff:
            if diff[0] == "add":
                return self._pa_gamma_add(lam, diff[1])
            if diff[0] == "rem":
                return self._pa_gamma_add(mu, diff[1])
        return sympy.S.One

    def _pa_e_diag(self, path, k):
        lam = path[k + 1][0]
        if path[k - 1][0] != lam:
            return sympy.S.Zero
        curr = path[k][0]
        if k % 2 == 0:
            if curr == lam:
                num = sympy.S.One
                den = sympy.S.One
                for beta in _removable_nodes(lam):
                    num *= (self.d - _content(beta) - _partition_size(lam))
                for beta in _addable_nodes(lam):
                    den *= (self.d - _content(beta) - _partition_size(lam))
                return sympy.simplify(num / den)
            diff = _diff_cell(lam, curr)
            if diff and diff[0] == "add":
                alpha = diff[1]
                num = (self.d - _content(alpha) - _partition_size(lam) - 1)
                den = (self.d - _content(alpha) - _partition_size(lam))
                num *= sympy.prod([_content(alpha) - _content(b) for b in _removable_nodes(lam)])
                den *= sympy.prod(
                    [_content(alpha) - _content(b) for b in _addable_nodes(lam) if b != alpha]
                )
                return sympy.simplify(num / den)
        else:
            if curr == lam:
                num = sympy.S.One
                den = sympy.S.One
                for beta in _addable_nodes(lam):
                    num *= (self.d - _content(beta) - _partition_size(lam))
                for beta in _removable_nodes(lam):
                    den *= (self.d - _content(beta) - _partition_size(lam))
                return sympy.simplify(num / den)
            diff = _diff_cell(lam, curr)
            if diff and diff[0] == "rem":
                alpha = diff[1]
                num = -(self.d - _content(alpha) - _partition_size(lam) + 1)
                den = (self.d - _content(alpha) - _partition_size(lam))
                num *= sympy.prod([_content(b) - _content(alpha) for b in _addable_nodes(lam)])
                den *= sympy.prod(
                    [_content(b) - _content(alpha) for b in _removable_nodes(lam) if b != alpha]
                )
                return sympy.simplify(num / den)
        return sympy.S.Zero

    def _pa_e(self, s, t, k):
        if not self._pa_i_equiv(s, t, k):
            return sympy.S.Zero
        lam_node = t[k + 1]
        lam = lam_node[0]
        if t[k - 1] != (lam, lam_node[1] - 1):
            return sympy.S.Zero
        e_tt = self._pa_e_diag(t, k)
        e_ss = self._pa_e_diag(s, k)

        if k % 2 == 0:
            if t[k][1] == lam_node[1] and s[k][1] == lam_node[1] - 1:
                mu = s[k][0]
                gamma = self._pa_gamma(lam, mu, True)
                return gamma ** -1
            if t[k][1] == lam_node[1] - 1 and s[k][1] == lam_node[1]:
                mu = t[k][0]
                gamma = self._pa_gamma(lam, mu, True)
                return e_ss * e_tt * gamma
            if t[k][1] == lam_node[1] - 1 and s[k][1] == lam_node[1] - 1:
                mu = t[k][0]
                rho = s[k][0]
                gamma_mu = self._pa_gamma(lam, mu, True)
                gamma_rho = self._pa_gamma(lam, rho, True)
                return (gamma_mu / gamma_rho) * e_tt
        else:
            if t[k][1] == lam_node[1] and s[k][1] == lam_node[1] - 1:
                mu = t[k][0]
                gamma = self._pa_gamma(mu, lam, False)
                return gamma
            if t[k][1] == lam_node[1] - 1 and s[k][1] == lam_node[1]:
                mu = s[k][0]
                gamma = self._pa_gamma(mu, lam, False)
                return e_ss * e_tt / gamma
            if t[k][1] == lam_node[1] and s[k][1] == lam_node[1]:
                mu = t[k][0]
                rho = s[k][0]
                gamma_mu = self._pa_gamma(mu, lam, False)
                gamma_rho = self._pa_gamma(rho, lam, False)
                return (gamma_mu / gamma_rho) * e_ss
            if t[k][1] == lam_node[1] - 1 and s[k][1] == lam_node[1] - 1:
                mu = t[k][0]
                rho = s[k][0]
                gamma_mu = self._pa_gamma(mu, lam, False)
                gamma_rho = self._pa_gamma(rho, lam, False)
                return (gamma_mu / gamma_rho) * e_tt
        return sympy.S.Zero

    def _pa_sigma(self, s, t, k, basis_paths):
        if not self._pa_i_approx(s, t, k):
            return sympy.S.Zero

        lam_kp1 = t[k + 1][0]
        lam_km1 = t[k - 1][0]
        lam_km2 = t[k - 2][0]
        lam_k = t[k][0]

        c_t_km1 = self._pa_content(t, k - 1)
        c_t_k = self._pa_content(t, k)
        c_t_kp1 = self._pa_content(t, k + 1)
        c_s_k = self._pa_content(s, k)
        c_s_kp1 = self._pa_content(s, k + 1)

        e_km1_tt = self._pa_e_diag(t, k - 1)

        t_swap = None
        for cand in self._pa_approx_candidates(t, k):
            if cand != t and self._pa_i_approx(cand, t, k):
                t_swap = cand
                break

        def delta_st():
            return sympy.S.One if s == t else sympy.S.Zero

        if k % 2 == 0:
            # Theorem 5.8 (sigmaeven)
            if lam_km1 == lam_kp1 and lam_km2 == lam_k:
                if s == t:
                    return sympy.simplify(c_t_k / e_km1_tt)
                num = self.d - c_s_k - c_t_km1 - e_km1_tt
                den = c_s_kp1 - c_t_km1
                return sympy.simplify((num / den) * self._pa_e(s, t, k))

            if lam_km1 != lam_kp1 and lam_km2 == lam_k:
                v = None
                for cand in basis_paths:
                    if self._pa_i_equiv(cand, t, k - 1) and cand[k - 1][0] == lam_kp1:
                        v = cand
                        break
                if v is None:
                    return sympy.S.Zero
                if s == v:
                    return sympy.simplify(c_t_k * self._pa_e(v, t, k))
                num = delta_st() - self._pa_e(v, t, k) * self._pa_e(s, v, k)
                den = c_s_kp1 - c_t_km1
                return sympy.simplify(num / den)

            if lam_km1 == lam_kp1 and lam_km2 != lam_k:
                if s[k][0] != lam_km2:
                    num = delta_st() + (self.d - c_t_km1 - c_s_k) * self._pa_e(s, t, k)
                    den = c_s_kp1 - c_t_km1
                    return sympy.simplify(num / den)
                # use symmetry with norms
                sigma_ts = self._pa_sigma(t, s, k, basis_paths)
                return sympy.simplify(sigma_ts * self._pa_norm(t) / self._pa_norm(s))

            if lam_km1 != lam_kp1 and lam_km2 != lam_k and t_swap is None:
                return sympy.simplify(delta_st() / (c_t_kp1 - c_t_km1))

            if lam_km1 != lam_kp1 and lam_km2 != lam_k and t_swap is not None:
                if s == t:
                    return sympy.simplify(sympy.S.One / (c_t_kp1 - c_t_km1))
                if s == t_swap:
                    comp = self._pa_compare_paths(s, t)
                    if comp == 1:
                        return sympy.simplify(sympy.S.One - sympy.S.One / ((c_t_kp1 - c_t_km1) ** 2))
                    return sympy.S.One
                return sympy.S.Zero

        else:
            # Theorem 5.9 (sigmaodd)
            if lam_km1 == lam_kp1 and lam_km2 == lam_k:
                if s == t:
                    return sympy.simplify(c_t_km1 / e_km1_tt)
                num = -e_km1_tt * self._pa_e(s, t, k - 1)
                den = c_s_kp1 - c_t_km1
                return sympy.simplify(num / den)

            if lam_km1 != lam_kp1 and lam_km2 == lam_k:
                v = None
                for cand in basis_paths:
                    if self._pa_i_equiv(cand, t, k - 1) and cand[k - 1][0] == lam_kp1:
                        v = cand
                        break
                if v is None:
                    return sympy.S.Zero
                if s == v:
                    return sympy.simplify(c_t_km1 * self._pa_e(v, t, k))
                num = (delta_st()
                       - self._pa_e(v, t, k) * self._pa_e(s, v, k)
                       + (self._pa_content(v, k - 1) - c_t_km1) * self._pa_e(s, t, k - 1))
                den = c_s_kp1 - c_t_km1
                return sympy.simplify(num / den)

            if lam_km1 == lam_kp1 and lam_km2 != lam_k:
                if s[k][0] != lam_km2:
                    return sympy.simplify(delta_st() / (c_t_kp1 - c_t_km1))
                sigma_ts = self._pa_sigma(t, s, k, basis_paths)
                return sympy.simplify(sigma_ts * self._pa_norm(t) / self._pa_norm(s))

            if lam_km1 != lam_kp1 and lam_km2 != lam_k and t_swap is None:
                return sympy.simplify(delta_st() / (c_t_kp1 - c_t_km1))

            if lam_km1 != lam_kp1 and lam_km2 != lam_k and t_swap is not None:
                if s == t:
                    return sympy.simplify(sympy.S.One / (c_t_kp1 - c_t_km1))
                if s == t_swap:
                    comp = self._pa_compare_paths(s, t)
                    if comp == 1:
                        return sympy.simplify(sympy.S.One - sympy.S.One / ((c_t_kp1 - c_t_km1) ** 2))
                    return sympy.S.One
                return sympy.S.Zero

        return sympy.S.Zero

    def _compute_irrep_matrices(self):
        irreps = self.irreps
        paths = self.bratteli_paths
        matrices = []
        for ir_idx, irrep in enumerate(irreps):
            basis_paths = paths[ir_idx]
            n = len(basis_paths)
            gen_mats = {}
            # compute e_i and sigma_i matrices for i=1..2k-1
            max_i = 2 * self.k - 1
            e_mats = {}
            sigma_mats = {}
            for i in range(1, max_i + 1):
                M_e = sympy.zeros(n, n)
                for a, s_path in enumerate(basis_paths):
                    for b, t_path in enumerate(basis_paths):
                        M_e[a, b] = self._pa_e(s_path, t_path, i)
                e_mats[i] = M_e
                if i >= 2:
                    M_sg = sympy.zeros(n, n)
                    for a, s_path in enumerate(basis_paths):
                        for b, t_path in enumerate(basis_paths):
                            M_sg[a, b] = self._pa_sigma(s_path, t_path, i, basis_paths)
                    sigma_mats[i] = M_sg

            for i in range(1, self.k + 1):
                gen_mats[f"p_{i}"] = e_mats[2 * i - 1]
            for i in range(1, self.k):
                gen_mats[f"b_{i}"] = e_mats[2 * i]
            for i in range(1, self.k):
                gen_mats[f"s_{i}"] = sigma_mats[2 * i]
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
                raise ValueError(f"Unknown generator {tok} for partition algebra.")
            if result is None:
                result = mats[tok]
            else:
                result = result * mats[tok]
        return result

    def _generate_basis(self):
        elts = list(range(1, self.k + 1)) + [-i for i in range(1, self.k + 1)]
        return self._all_partitions(elts)

    def _all_partitions(self, elts):
        if not elts:
            return [SetPartition([])]
        if len(elts) == 1:
            return [SetPartition([{elts[0]}])]
        first, rest = elts[0], elts[1:]
        out = []
        for p in self._all_partitions(rest):
            out.append(SetPartition(p.blocks + [{first}]))
            for idx, b in enumerate(p.blocks):
                new_blocks = list(p.blocks)
                new_blocks[idx] = set(b) | {first}
                out.append(SetPartition(new_blocks))
        return out
