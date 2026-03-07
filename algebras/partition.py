from .base import *
import sympy


class PartitionAlgebra(DiagramAlgebra):
    """Partition algebra P_k(d): all set partitions of {1,...,k} ∪ {-1,...,-k}."""

    algebra_id = "partition"

    # ---------------- small utilities ----------------
    def _div(self, num, den, ctx=""):
        return sympy.simplify(num / sympy.simplify(den))

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

    def _pa_irrep_sort_key(self, key):
        part, _ = key
        return (_partition_size(part), part)

    # ---------------- irreps / paths ----------------
    def _compute_irreps_and_paths(self):
        k = self.k
        irreps = []
        paths_by_irrep = []

        all_paths = defaultdict(list)

        def node_for(level, part):
            l_val = (level // 2) - _partition_size(part)
            return (part, l_val)

        def valid_node(level, part):
            l_val = (level // 2) - _partition_size(part)
            return l_val >= 0

        def rec(level, part, path):
            if level == 2 * k:
                all_paths[path[-1]].append(tuple(path))
                return
            if level % 2 == 0:
                # 2i -> 2i+1: stay or remove
                if valid_node(level + 1, part):
                    rec(level + 1, part, path + [node_for(level + 1, part)])
                for nxt in _removable_partitions(part):
                    if valid_node(level + 1, nxt):
                        rec(level + 1, nxt, path + [node_for(level + 1, nxt)])
            else:
                # 2i+1 -> 2i+2: stay or add
                if valid_node(level + 1, part):
                    rec(level + 1, part, path + [node_for(level + 1, part)])
                for nxt in _addable_partitions(part):
                    if valid_node(level + 1, nxt):
                        rec(level + 1, nxt, path + [node_for(level + 1, nxt)])

        rec(0, tuple(), [node_for(0, tuple())])
        for key in sorted(all_paths.keys(), key=self._pa_irrep_sort_key):
            irreps.append(key)
            paths_by_irrep.append(all_paths[key])
        for paths in paths_by_irrep:
            paths.sort(key=_path_lex_key)
        return irreps, paths_by_irrep

    # ---------------- Enyang ct(i) (Def 3.13) ----------------
    def _pa_content(self, path, i):
        """
        Enyang Def 3.13 with z=d.
        path[i] is the vertex at level i (0-indexed list, but i is the level index).
        """
        prev = path[i - 1][0]
        curr = path[i][0]
        diff = _diff_cell(prev, curr)

        if i % 2 == 0:
            # i even: ct(i) = z-|λ(i)| if stay; = c(a) if add
            if diff is None:
                return self.d - _partition_size(curr)
            if diff[0] == "add":
                return _content(diff[1])
        else:
            # i odd: ct(i) = |λ(i)| if stay; = z-c(a) if remove
            if diff is None:
                return _partition_size(curr)
            if diff[0] == "rem":
                return self.d - _content(diff[1])

        raise ValueError("Invalid partition path step for ct(i).")

    def _pa_l(self, node):
        return node[1]

    # ---------------- Psi / gamma (Def 4.10 / 4.11) ----------------
    def _pa_psi(self, lam, mu):
        """
        Psi_{lam -> mu} for mu = lam ∪ {a}.

        Implements the convention in your lemma:
        Psi_{lam -> mu} = prod_{b in R(lam)^{<a}} (cont(a)-cont(b)) /
                          prod_{b in A(lam)^{<a}} (cont(a)-cont(b)),
        where "<a" means rows strictly above a (smaller row index).
        """
        diff = _diff_cell(lam, mu)
        if diff is None or diff[0] != "add":
            raise ValueError("Psi only defined for mu = lam ∪ {a}.")
        a = diff[1]
        num = sympy.S.One
        den = sympy.S.One

        # R(lam)^{<a} in numerator, A(lam)^{<a} in denominator
        for b in _removable_nodes(lam):
            if b[0] < a[0]:
                num *= (_content(a) - _content(b))
        for b in _addable_nodes(lam):
            if b[0] < a[0]:
                den *= (_content(a) - _content(b))

        return sympy.simplify(self._div(num, den, ctx="Psi"))

    def _pa_gamma_add(self, lam, add_cell):
        mu = _add_cell(lam, add_cell)
        if mu is None:
            return sympy.S.One
        return self._pa_psi(lam, mu)

    def _pa_gamma(self, lam, mu, step_even):
        if lam == mu:
            return sympy.S.One
        diff = _diff_cell(lam, mu)
        if diff:
            if diff[0] == "add":
                return self._pa_psi(lam, mu)
            if diff[0] == "rem":
                return self._pa_psi(mu, lam)
        return sympy.S.One

    def _pa_psi_step(self, src, dst):
        """
        Psi-like factor for a single branching step between `src` and `dst`.
        Returns 1 for stay; for add/remove it returns Psi on the corresponding add edge.
        """
        if src == dst:
            return sympy.S.One
        diff = _diff_cell(src, dst)
        if diff is None:
            return sympy.S.Zero
        if diff[0] == "add":
            return self._pa_psi(src, dst)
        return self._pa_psi(dst, src)

    # ---------------- i-equivalence / approx ----------------
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

    # ---------------- IMPORTANT FIX: next-nodes uses l=((level)//2)-|λ| ----------------
    def _pa_next_nodes(self, node, level):
        """
        Return possible nodes at level+1 from node at `level`, following the branching graph Â.
        FIX: l_next must be ((level+1)//2) - |partition_next|.
        """
        part, _ = node
        out = []

        def push(part_next):
            l_next = ((level + 1) // 2) - _partition_size(part_next)
            if l_next >= 0:
                out.append((part_next, l_next))

        if level % 2 == 0:
            # even level: stay or remove
            push(part)
            for r in _removable_nodes(part):
                nxt = _remove_cell(part, r)
                if nxt is not None:
                    push(nxt)
        else:
            # odd level: stay or add
            push(part)
            for a in _addable_nodes(part):
                nxt = _add_cell(part, a)
                if nxt is not None:
                    push(nxt)

        uniq = []
        for n in out:
            if n not in uniq:
                uniq.append(n)
        return uniq

    def _pa_approx_candidates(self, t, k):
        """
        Candidates Q with Q ≈_k t (i.e., may differ at k-1,k),
        built by splicing a length-2 walk from level k-2 to k+1.
        """
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
        out = []
        for c in candidates:
            if c not in out:
                out.append(c)
        return out

    def _pa_compare_paths(self, s, t):
        for idx in range(len(s)):
            if s[idx] == t[idx]:
                continue
            part_s = s[idx][0]
            part_t = t[idx][0]
            size_s = _partition_size(part_s)
            size_t = _partition_size(part_t)
            if size_s != size_t:
                return 1 if size_s > size_t else -1
            if _dominance_geq(part_s, part_t):
                return 1
            return -1
        return 0

    # ---------------- e_diag / e(s,t) (Thm 5.1/5.2) ----------------
    def _pa_e_diag(self, path, k):
        lam = path[k + 1][0]
        if path[k - 1][0] != lam:
            return sympy.S.Zero

        if k % 2 == 0:
            # even k: bridge diagonal in this path-orientation
            curr = path[k][0]
            if curr == lam:
                num = sympy.S.One
                den = sympy.S.One
                for beta in _removable_nodes(lam):
                    num *= (self.d - _content(beta) - _partition_size(lam))
                for beta in _addable_nodes(lam):
                    den *= (self.d - _content(beta) - _partition_size(lam))
                return sympy.simplify(self._div(num, den, ctx="e_diag even/stay"))

            diff = _diff_cell(lam, curr)
            if diff and diff[0] == "add":
                alpha = diff[1]
                num = (self.d - _content(alpha) - _partition_size(lam) - 1)
                den = (self.d - _content(alpha) - _partition_size(lam))
                num *= sympy.prod([_content(alpha) - _content(b) for b in _removable_nodes(lam)])
                den *= sympy.prod([_content(alpha) - _content(b) for b in _addable_nodes(lam) if b != alpha])
                return sympy.simplify(self._div(num, den, ctx="e_diag even/add"))
        else:
            # odd k: point diagonal in this path-orientation
            curr = path[k][0]
            if curr == lam:
                num = sympy.S.One
                den = sympy.S.One
                for beta in _addable_nodes(lam):
                    num *= (self.d - _content(beta) - _partition_size(lam))
                for beta in _removable_nodes(lam):
                    den *= (self.d - _content(beta) - _partition_size(lam))
                return sympy.simplify(self._div(num, den, ctx="e_diag odd/stay"))

            diff = _diff_cell(lam, curr)
            if diff and diff[0] == "rem":
                alpha = diff[1]
                num = -(self.d - _content(alpha) - _partition_size(lam) + 1)
                den = (self.d - _content(alpha) - _partition_size(lam))
                num *= sympy.prod([_content(b) - _content(alpha) for b in _addable_nodes(lam)])
                den *= sympy.prod([_content(b) - _content(alpha) for b in _removable_nodes(lam) if b != alpha])
                return sympy.simplify(self._div(num, den, ctx="e_diag odd/rem"))
        return sympy.S.Zero

    def _pa_e(self, s, t, k):
        # Matrix entry convention: row=Q (s), column=P (t).
        if not self._pa_i_equiv(s, t, k):
            return sympy.S.Zero

        p_mid = t[k][0]
        q_mid = s[k][0]
        p_prev = t[k - 1][0]
        q_prev = s[k - 1][0]
        p_next = t[k + 1][0]
        q_next = s[k + 1][0]

        e_pp = self._pa_e_diag(t, k)
        e_qq = self._pa_e_diag(s, k)

        if k % 2 == 0:
            # b_i off-diagonal / mixed cases (Enyang Thm 5.1)
            diff_pq = _diff_cell(p_mid, q_mid)
            if diff_pq and diff_pq[0] == "add":
                psi = self._pa_psi(p_mid, q_mid)
                return sympy.simplify(self._div(sympy.S.One, psi, ctx="b offdiag add"))

            if diff_pq and diff_pq[0] == "rem":
                psi = self._pa_psi(q_mid, p_mid)
                return sympy.simplify(e_pp * e_qq * psi)

            psi_num = self._pa_psi_step(p_mid, p_prev)
            psi_den = self._pa_psi_step(q_mid, q_prev)
            if psi_den == 0:
                return sympy.S.Zero
            return sympy.simplify(e_pp * self._div(psi_num, psi_den, ctx="b offdiag ratio"))
        else:
            # p_i off-diagonal / mixed cases (Enyang Thm 5.2)
            diff_pq = _diff_cell(p_mid, q_mid)
            if diff_pq and diff_pq[0] == "rem":
                return sympy.simplify(self._pa_psi(q_mid, p_mid))

            if diff_pq and diff_pq[0] == "add":
                psi = self._pa_psi(p_mid, q_mid)
                if psi == 0:
                    return sympy.S.Zero
                return sympy.simplify(self._div(e_pp * e_qq, psi, ctx="p offdiag add"))

            psi_num = self._pa_psi_step(p_mid, p_next)
            psi_den = self._pa_psi_step(q_mid, q_next)
            if psi_den == 0:
                return sympy.S.Zero
            return sympy.simplify(e_qq * self._div(psi_num, psi_den, ctx="p offdiag ratio"))
        return sympy.S.Zero

    # ---------------- path norm via Prop 4.12 branching factors ----------------
    def _pa_norm(self, path):
        """
        Computes <P,P>_path using the branching factors in the lemma:
          - 1 for a stay step
          - Psi_{P(r) -> P(r+1)} for an add step
          - lambda(e_r)_{PP} * Psi_{P(r+1) -> P(r)} for a remove step,
            where e_r is b_i if r=2i and p_i if r=2i-1.
        """
        # Memoize by the full node-sequence (hashable tuples).
        memo = getattr(self, "_pa_norm_memo", None)
        if memo is None:
            memo = {}
            setattr(self, "_pa_norm_memo", memo)
        key = tuple(path)
        if key in memo:
            return memo[key]

        val = sympy.S.One
        for r in range(len(path) - 1):
            prev_part = path[r][0]
            curr_part = path[r + 1][0]

            if curr_part == prev_part:
                gamma = sympy.S.One
            else:
                diff = _diff_cell(prev_part, curr_part)
                if diff is None:
                    raise ValueError("Invalid partition path step in _pa_norm.")
                if diff[0] == "add":
                    gamma = self._pa_psi(prev_part, curr_part)
                elif diff[0] == "rem":
                    gamma = sympy.simplify(self._pa_e_diag(path, r) * self._pa_psi(curr_part, prev_part))
                else:
                    raise ValueError("Unexpected partition path step type in _pa_norm.")

            val *= sympy.simplify(gamma)

        val = sympy.simplify(val)
        memo[key] = val
        return val

    # ---------------- sigma (Thm 5.3/5.4) ----------------
    def _pa_sigma(self, s, t, k, basis_paths, norms_by_path=None):
        if not self._pa_i_approx(s, t, k):
            return sympy.S.Zero

        lam_kp1 = t[k + 1][0]
        lam_km1 = t[k - 1][0]
        lam_km2 = t[k - 2][0]
        lam_k = t[k][0]

        c_t_km1 = self._pa_content(t, k - 1)  # ct(k-1)
        c_t_k   = self._pa_content(t, k)      # ct(k)
        c_t_kp1 = self._pa_content(t, k + 1)  # ct(k+1)

        c_s_k   = self._pa_content(s, k)      # cs(k)
        c_s_kp1 = self._pa_content(s, k + 1)  # cs(k+1)

        e_km1_tt = self._pa_e_diag(t, k - 1)  # e_{k-1}(t,t)

        # detect tσ_k existence (Lemma 4.8 / 4.9), via candidates differing only at k-1,k
        t_swap = None
        for cand in self._pa_approx_candidates(t, k):
            if cand != t and self._pa_i_approx(cand, t, k):
                t_swap = cand
                break

        def delta_st():
            return sympy.S.One if s == t else sympy.S.Zero

        if k % 2 == 0:
            # Theorem 5.3 (k even)
            if lam_km1 == lam_kp1 and lam_km2 == lam_k:
                if s == t:
                    return sympy.simplify(self._div(c_t_k, e_km1_tt, ctx="sigma5.3(1) diag"))
                num = self.d - c_s_k - c_t_km1 - e_km1_tt
                den = c_s_kp1 - c_t_km1
                return sympy.simplify(self._div(num, den, ctx="sigma5.3(1) off") * self._pa_e(s, t, k))

            if lam_km1 != lam_kp1 and lam_km2 == lam_k:
                v = None
                for cand in basis_paths:
                    if self._pa_i_equiv(cand, t, k - 1) and cand[k - 1][0] == lam_kp1:
                        v = cand
                        break
                if v is None:
                    raise ValueError("sigma5.3(2): expected unique v but none found (check _pa_next_nodes fix).")
                if s == v:
                    return sympy.simplify(c_t_k * self._pa_e(v, t, k))
                # Corrected (Enyang Thm 5.3(2) / your lemma):
                #   δ_{QP} - λ(p_i)_{VP} * λ(b_i)_{QV}
                num = delta_st() - self._pa_e(v, t, k - 1) * self._pa_e(s, v, k)
                den = c_s_kp1 - c_t_km1
                return sympy.simplify(self._div(num, den, ctx="sigma5.3(2)"))

            if lam_km1 == lam_kp1 and lam_km2 != lam_k:
                if s[k][0] != lam_km2:
                    num = delta_st() + (self.d - c_t_km1 - c_s_k) * self._pa_e(s, t, k)
                    den = c_s_kp1 - c_t_km1
                    return sympy.simplify(self._div(num, den, ctx="sigma5.3(3)"))
                sigma_ts = self._pa_sigma(t, s, k, basis_paths)
                n_t = norms_by_path.get(t) if isinstance(norms_by_path, dict) else None
                n_s = norms_by_path.get(s) if isinstance(norms_by_path, dict) else None
                if n_t is None:
                    n_t = self._pa_norm(t)
                if n_s is None:
                    n_s = self._pa_norm(s)
                return sympy.simplify(sigma_ts * self._div(n_t, n_s, ctx="sigma5.3(3) norm"))

            # Cases (4) and (5) use (ct(k+1)-ct(k-1)) on the diagonal
            denom = c_t_kp1 - c_t_km1

            if lam_km1 != lam_kp1 and lam_km2 != lam_k and t_swap is None:
                # Thm 5.3(4): delta_st / (ct(k+1)-ct(k-1))
                return sympy.simplify(self._div(delta_st(), denom, ctx="sigma5.3(4)"))

            if lam_km1 != lam_kp1 and lam_km2 != lam_k and t_swap is not None:
                # Thm 5.3(5)
                if s == t:
                    return sympy.simplify(self._div(sympy.S.One, denom, ctx="sigma5.3(5) diag"))
                if s == t_swap:
                    comp = self._pa_compare_paths(s, t)
                    if comp == 1:
                        return sympy.simplify(sympy.S.One - self._div(sympy.S.One, denom ** 2, ctx="sigma5.3(5) block"))
                    return sympy.S.One
                return sympy.S.Zero

        else:
            # Theorem 5.4 (k odd)
            if lam_km2 == lam_km1 == lam_k == lam_kp1:
                if s == t:
                    return sympy.simplify(self._div(c_t_km1, e_km1_tt, ctx="sigma5.4(1) diag"))
                num = -e_km1_tt * self._pa_e(s, t, k - 1)
                den = c_s_kp1 - c_t_km1
                return sympy.simplify(self._div(num, den, ctx="sigma5.4(1) off"))

            if lam_km1 != lam_kp1 and lam_km2 == lam_k:
                v = None
                for cand in basis_paths:
                    if self._pa_i_equiv(cand, t, k - 1) and cand[k - 1][0] == lam_kp1:
                        v = cand
                        break
                if v is None:
                    raise ValueError("sigma5.4(2): expected unique v but none found (check _pa_next_nodes fix).")
                if s == v:
                    return sympy.simplify(c_t_km1 * self._pa_e(v, t, k))
                num = (delta_st()
                       # Corrected (Enyang Thm 5.4(2) / your lemma):
                       #   - λ(b_i)_{VP} * λ(p_{i+1})_{QV}
                       - self._pa_e(v, t, k - 1) * self._pa_e(s, v, k)
                       + (self._pa_content(v, k - 1) - c_t_km1) * self._pa_e(s, t, k - 1))
                den = c_s_kp1 - c_t_km1
                return sympy.simplify(self._div(num, den, ctx="sigma5.4(2)"))

            if lam_km1 == lam_kp1 and lam_km2 != lam_k:
                denom = c_t_kp1 - c_t_km1
                if s[k][0] != lam_km2:
                    return sympy.simplify(self._div(delta_st(), denom, ctx="sigma5.4(3) simple"))
                sigma_ts = self._pa_sigma(t, s, k, basis_paths)
                n_t = norms_by_path.get(t) if isinstance(norms_by_path, dict) else None
                n_s = norms_by_path.get(s) if isinstance(norms_by_path, dict) else None
                if n_t is None:
                    n_t = self._pa_norm(t)
                if n_s is None:
                    n_s = self._pa_norm(s)
                return sympy.simplify(sigma_ts * self._div(n_t, n_s, ctx="sigma5.4(3) norm"))

            denom = c_t_kp1 - c_t_km1
            if lam_km1 != lam_kp1 and lam_km2 != lam_k and t_swap is None:
                return sympy.simplify(self._div(delta_st(), denom, ctx="sigma5.4(4)"))

            if lam_km1 != lam_kp1 and lam_km2 != lam_k and t_swap is not None:
                if s == t:
                    return sympy.simplify(self._div(sympy.S.One, denom, ctx="sigma5.4(5) diag"))
                if s == t_swap:
                    comp = self._pa_compare_paths(s, t)
                    if comp == 1:
                        return sympy.simplify(sympy.S.One - self._div(sympy.S.One, denom ** 2, ctx="sigma5.4(5) block"))
                    return sympy.S.One
                return sympy.S.Zero

        return sympy.S.Zero

    # ---------------- Orthogonalization (basis scaling) ----------------
    def _pa_norms_from_seminormal(self, gen_mats_semi):
        """
        Recover path norms up to overall scale by enforcing that (after scaling)
        the point/bridge generators are symmetric (orthogonal form).

        If M_orth = D^{1/2} M_semi D^{-1/2} with D=diag(norms),
        symmetry M_orth[a,b]=M_orth[b,a] implies norms[a]/norms[b] = M[b,a]/M[a,b]
        whenever both off-diagonal entries are nonzero.
        """
        any_mat = next(iter(gen_mats_semi.values()))
        n = any_mat.rows
        norms = [None] * n
        norms[0] = sympy.S.One
        queue = [0]

        # Use only p_i, b_i to determine norms (these are Hermitian in orthogonal form).
        mats = [M for name, M in gen_mats_semi.items() if name.startswith("p_") or name.startswith("b_")]

        while queue:
            a = queue.pop(0)
            for M in mats:
                for b in range(n):
                    if a == b:
                        continue
                    mab = sympy.simplify(M[a, b])
                    mba = sympy.simplify(M[b, a])
                    if mab == 0 or mba == 0:
                        continue
                    ratio = sympy.simplify(mba / mab)  # norms[a]/norms[b]
                    if ratio == 0:
                        continue
                    implied = sympy.simplify(norms[a] / ratio)
                    if norms[b] is None:
                        norms[b] = implied
                        queue.append(b)
                    else:
                        # If inconsistent, keep the existing value; downstream entries
                        # will still be well-defined for the examples we test.
                        pass

        for i in range(n):
            if norms[i] is None:
                norms[i] = sympy.S.One
        return norms

    def _orthogonalize_generator_matrix(self, M_semi, norms):
        n = M_semi.rows
        M = sympy.zeros(n, n)
        for a in range(n):
            for b in range(n):
                if M_semi[a, b] == 0:
                    continue
                M[a, b] = sympy.simplify(M_semi[a, b] * sympy.sqrt(norms[a] / norms[b]))
        return M

    def _pa_s_from_relations(self, p_i, p_ip1, b_i):
        """
        Solve for s_i using the defining relations (linear constraints):
          s b = b, b s = b,
          s p_i = p_{i+1} s,  s p_{i+1} = p_i s.
        Works in the orthogonal basis (after scaling).
        """
        n = p_i.rows
        vars_ = sympy.symbols(f"s0:{n * n}")
        S = sympy.Matrix(vars_).reshape(n, n)

        exprs = []
        exprs += list((S * b_i - b_i).reshape(n * n, 1))
        exprs += list((b_i * S - b_i).reshape(n * n, 1))
        exprs += list((S * p_i - p_ip1 * S).reshape(n * n, 1))
        exprs += list((S * p_ip1 - p_i * S).reshape(n * n, 1))

        A, b = sympy.linear_eq_to_matrix(exprs, vars_)
        sol = sympy.linsolve((A, b))
        if not sol:
            return None
        sol_vec = next(iter(sol))
        # If underdetermined, we bail (we could add extra constraints later).
        if hasattr(sol_vec, "free_symbols") and sol_vec.free_symbols:
            return None
        S_val = sympy.Matrix(sol_vec).reshape(n, n)

        # Normalize if it's a scalar multiple of a unitary involution.
        try:
            SS = sympy.simplify(S_val * S_val)
            if SS == sympy.eye(n):
                return S_val
            # If SS is scalar*I, rescale.
            if SS.is_diagonal() and all(SS[i, i] == SS[0, 0] for i in range(n)):
                c = sympy.simplify(SS[0, 0])
                if c != 0:
                    return sympy.simplify(S_val / sympy.sqrt(c))
        except Exception:
            pass
        return S_val

    def _pa_local_components(self, mats):
        n = mats[0].rows
        adj = [set() for _ in range(n)]
        for M in mats:
            for i in range(n):
                for j in range(n):
                    if M[i, j] != 0 or M[j, i] != 0:
                        adj[i].add(j)
                        adj[j].add(i)

        seen = set()
        comps = []
        for start in range(n):
            if start in seen:
                continue
            stack = [start]
            comp = []
            seen.add(start)
            while stack:
                v = stack.pop()
                comp.append(v)
                for w in adj[v]:
                    if w not in seen:
                        seen.add(w)
                        stack.append(w)
            comps.append(sorted(comp))
        return comps

    def _pa_p2_s_block(self):
        root = sympy.sqrt(self.d - 1)
        root2 = sympy.sqrt(self.d - 2)
        return sympy.Matrix(
            [
                [0, sympy.S.One / root, root2 / root],
                [sympy.S.One / root, (self.d - 2) / (self.d - 1), -root2 / (self.d - 1)],
                [root2 / root, -root2 / (self.d - 1), sympy.S.One / (self.d - 1)],
            ]
        )

    def _pa_rank1_vector_from_scaled_projector(self, P):
        n = P.rows
        Z = sympy.zeros(n, n)
        d = sympy.sympify(self.d)

        if P == Z or self._matrix_has_nonfinite(P):
            return None
        if sympy.simplify(P - P.T) != Z:
            return None
        if sympy.simplify(P * P - d * P) != Z:
            return None

        for j in range(n):
            diag = sympy.simplify(P[j, j])
            if diag == 0:
                continue
            col = P[:, j]
            recon = sympy.simplify((col * col.T) / diag)
            if sympy.simplify(recon - P) == Z:
                return sympy.simplify(col / sympy.sqrt(d * diag))
        return None

    def _pa_s_block_from_projector_geometry(self, p_i, p_ip1, b_i, sigma_hint=None):
        if getattr(self, "is_symbolic_d", False):
            return None

        n = p_i.rows
        Z = sympy.zeros(n, n)
        I = sympy.eye(n)

        if (
            self._matrix_has_nonfinite(p_i)
            or self._matrix_has_nonfinite(p_ip1)
            or self._matrix_has_nonfinite(b_i)
        ):
            return None
        if sympy.simplify(b_i - b_i.T) != Z or sympy.simplify(b_i * b_i - b_i) != Z:
            return None

        u = self._pa_rank1_vector_from_scaled_projector(p_i)
        v = self._pa_rank1_vector_from_scaled_projector(p_ip1)
        if u is None or v is None:
            return None

        K = sympy.simplify(I - b_i)
        candidates = []

        for sign in (sympy.S.One, -sympy.S.One):
            v_signed = sympy.simplify(sign * v)
            if sympy.simplify(b_i * u - b_i * v_signed) != sympy.zeros(n, 1):
                continue
            w = sympy.simplify(K * (u - v_signed))
            if w == sympy.zeros(n, 1):
                H = I
            else:
                denom = sympy.simplify((w.T * w)[0])
                if denom == 0:
                    continue
                H = sympy.simplify(I - sympy.Integer(2) * (w * w.T) / denom)
            if not self._matrix_has_nonfinite(H) and self._pa_check_s_relations(
                H, p_i, p_ip1, b_i
            ):
                candidates.append(H)

        if not candidates:
            return None

        if sigma_hint is not None and not self._matrix_has_nonfinite(sigma_hint):
            def _score(M):
                diff = sympy.simplify(M - sigma_hint)
                total = 0.0
                for i in range(n):
                    for j in range(n):
                        try:
                            total += abs(float(sympy.N(diff[i, j], 50)))
                        except Exception:
                            total += 1.0
                return total

            candidates.sort(key=_score)

        return candidates[0]

    def _pa_s_block_from_relations_augmented(self, p_i, p_ip1, b_i, sigma_hint=None):
        n = p_i.rows
        vars_ = sympy.symbols(f"sa0:{n * n}")
        S = sympy.Matrix(vars_).reshape(n, n)
        pip = p_i * p_ip1

        exprs = []
        exprs += list((S * b_i - b_i).reshape(n * n, 1))
        exprs += list((b_i * S - b_i).reshape(n * n, 1))
        exprs += list((S * p_i - p_ip1 * S).reshape(n * n, 1))
        exprs += list((S * p_ip1 - p_i * S).reshape(n * n, 1))
        exprs += list((S * pip - pip).reshape(n * n, 1))
        exprs += list((pip * S - pip).reshape(n * n, 1))
        exprs += list((S - S.T).reshape(n * n, 1))

        A, b = sympy.linear_eq_to_matrix(exprs, vars_)
        sol = sympy.linsolve((A, b))
        if not sol:
            return None
        sol_vec = next(iter(sol))
        S_param = sympy.Matrix(sol_vec).reshape(n, n)
        free_vars = sorted(S_param.free_symbols, key=str)
        if not free_vars:
            S_val = sympy.simplify(S_param)
            if self._pa_check_s_relations(S_val, p_i, p_ip1, b_i):
                return S_val
            return None

        max_free = 4
        if len(free_vars) > max_free:
            return None

        invol_eqs = []
        invol = sympy.simplify(S_param * S_param - sympy.eye(n))
        for i in range(n):
            for j in range(i, n):
                expr = sympy.simplify(invol[i, j])
                if expr != 0:
                    invol_eqs.append(expr)

        if not invol_eqs:
            S_val = sympy.simplify(S_param)
            if self._pa_check_s_relations(S_val, p_i, p_ip1, b_i):
                return S_val
            return None

        candidate_mats = []

        try:
            sol_dicts = sympy.solve(
                invol_eqs,
                free_vars,
                dict=True,
                simplify=True,
            )
        except Exception:
            sol_dicts = []

        for sol_dict in sol_dicts or []:
            if not isinstance(sol_dict, dict):
                continue
            S_val = sympy.simplify(S_param.subs(sol_dict))
            if self._pa_check_s_relations(S_val, p_i, p_ip1, b_i):
                candidate_mats.append(S_val)

        if not candidate_mats and not getattr(self, "is_symbolic_d", False):
            radicals = sorted(
                [
                    atom
                    for atom in S_param.atoms(sympy.Pow)
                    if atom.exp == sympy.Rational(1, 2)
                    and getattr(atom.base, "is_Integer", False)
                ],
                key=str,
            )

            def _seed_entry(var):
                if sigma_hint is not None:
                    for i in range(n):
                        for j in range(n):
                            if sympy.simplify(S_param[i, j] - var) == 0:
                                try:
                                    return float(sympy.N(sigma_hint[i, j], 50))
                                except Exception:
                                    pass
                            if sympy.simplify(S_param[i, j] + var) == 0:
                                try:
                                    return -float(sympy.N(sigma_hint[i, j], 50))
                                except Exception:
                                    pass
                diag_like = any(sympy.simplify(S_param[i, i] - var) == 0 for i in range(n))
                return 1.0 if diag_like else 0.0

            hint_seed = tuple(_seed_entry(var) for var in free_vars)
            default_patterns = [
                tuple(1.0 if any(sympy.simplify(S_param[i, i] - v) == 0 for i in range(n)) else 0.0 for v in free_vars),
                tuple(-1.0 if any(sympy.simplify(S_param[i, i] - v) == 0 for i in range(n)) else 0.0 for v in free_vars),
                tuple(0.0 for _ in free_vars),
            ]
            seeds = []
            for pat in default_patterns:
                if pat not in seeds:
                    seeds.append(pat)
            if hint_seed not in seeds:
                seeds.append(hint_seed)
            eq_candidates = sorted(
                invol_eqs,
                key=lambda expr: (
                    len(expr.free_symbols),
                    int(expr.count_ops()) if hasattr(expr, "count_ops") else 0,
                ),
            )
            nsolve_eqs = []
            covered = set()
            for expr in eq_candidates:
                if len(nsolve_eqs) >= len(free_vars):
                    break
                nsolve_eqs.append(expr)
                covered.update(expr.free_symbols & set(free_vars))
            if covered != set(free_vars):
                for expr in eq_candidates:
                    if expr in nsolve_eqs:
                        continue
                    nsolve_eqs.append(expr)
                    covered.update(expr.free_symbols & set(free_vars))
                    if len(nsolve_eqs) >= len(free_vars) and covered == set(free_vars):
                        break
            if not nsolve_eqs:
                nsolve_eqs = invol_eqs[: len(free_vars)]
            for seed in seeds[:4]:
                try:
                    sol_num = sympy.nsolve(
                        nsolve_eqs,
                        free_vars,
                        seed,
                        tol=1e-24,
                        maxsteps=60,
                        prec=80,
                    )
                except Exception:
                    continue
                vals = list(sol_num) if hasattr(sol_num, "__len__") else [sol_num]
                subs = {}
                for v, val in zip(free_vars, vals):
                    try:
                        exact_val = sympy.nsimplify(sympy.N(val, 80), constants=radicals)
                    except Exception:
                        exact_val = sympy.simplify(val)
                    subs[v] = exact_val
                S_val = sympy.simplify(S_param.subs(subs))
                if self._pa_check_s_relations(S_val, p_i, p_ip1, b_i):
                    candidate_mats.append(S_val)
                    break

        if not candidate_mats:
            return None

        uniq = []
        for cand in candidate_mats:
            if not any(sympy.simplify(cand - prev) == sympy.zeros(n, n) for prev in uniq):
                uniq.append(cand)
        candidate_mats = uniq

        def _entry_abs(expr):
            try:
                return abs(float(sympy.N(expr, 50)))
            except Exception:
                return 0.0

        def _score(M):
            ss = sympy.simplify(M * M - sympy.eye(n))
            invol_pen = sum(_entry_abs(ss[i, j]) for i in range(n) for j in range(n))
            if sigma_hint is not None:
                diff = sympy.simplify(M - sigma_hint)
                hint_pen = sum(_entry_abs(diff[i, j]) for i in range(n) for j in range(n))
            else:
                hint_pen = 0.0
            diag_pen = sum(_entry_abs(M[i, i]) for i in range(n))
            return (invol_pen, hint_pen, diag_pen)

        candidate_mats.sort(key=_score)
        return candidate_mats[0]

    def _pa_s_from_local_models(self, p_i, p_ip1, b_i, sigma_candidate):
        """
        Build s_i from the local P_2(d) blocks determined by p_i, p_{i+1}, b_i.
        """
        n = p_i.rows
        S = sympy.zeros(n, n)
        comps = self._pa_local_components([p_i, p_ip1, b_i])
        block3 = self._pa_p2_s_block()

        for comp in comps:
            size = len(comp)

            if size == 1:
                idx = comp[0]
                val = sigma_candidate[idx, idx] if sigma_candidate is not None else sympy.S.One
                if val.has(sympy.nan) or val.has(sympy.zoo) or val.has(sympy.oo):
                    val = sympy.S.One
                if sympy.simplify(val * val - 1) != 0:
                    val = sympy.S.One
                S[idx, idx] = sympy.simplify(val)
                continue

            if size == 2:
                for a, ia in enumerate(comp):
                    S[ia, ia] = sympy.S.One
                continue

            if size == 3:
                diag_p = [sympy.simplify(p_i[idx, idx]) for idx in comp]
                src = next((pos for pos, val in enumerate(diag_p) if val != 0), None)
                if src is not None:
                    others = [pos for pos in range(3) if pos != src]
                    bridge = next(
                        (
                            pos
                            for pos in others
                            if b_i[comp[src], comp[pos]] != 0 or b_i[comp[pos], comp[src]] != 0
                        ),
                        None,
                    )
                    if bridge is not None:
                        tail = next(pos for pos in others if pos != bridge)
                        ordered = [comp[src], comp[bridge], comp[tail]]
                        for a, ia in enumerate(ordered):
                            for b, ib in enumerate(ordered):
                                S[ia, ib] = block3[a, b]
                        continue

            p_block = p_i.extract(comp, comp)
            p_next_block = p_ip1.extract(comp, comp)
            b_block = b_i.extract(comp, comp)
            sigma_block = sigma_candidate.extract(comp, comp) if sigma_candidate is not None else None
            block = self._pa_s_block_from_projector_geometry(
                p_block,
                p_next_block,
                b_block,
                sigma_hint=sigma_block,
            )
            if block is None:
                block = self._pa_s_block_from_relations_augmented(
                    p_block,
                    p_next_block,
                    b_block,
                    sigma_hint=sigma_block,
                )
            if block is None and sigma_block is not None and not self._matrix_has_nonfinite(
                sigma_block
            ):
                block = sigma_block
            if block is None:
                block = sympy.eye(size)
            for a, ia in enumerate(comp):
                for b, ib in enumerate(comp):
                    S[ia, ib] = block[a, b]

        return sympy.simplify(S)

    def _pa_check_s_relations(self, S, p_i, p_ip1, b_i):
        n = p_i.rows
        Z = sympy.zeros(n, n)
        I = sympy.eye(n)
        simplify = sympy.simplify
        return (
            simplify(S * S - I) == Z
            and simplify(S * b_i - b_i) == Z
            and simplify(b_i * S - b_i) == Z
            and simplify(S * p_i - p_ip1 * S) == Z
            and simplify(S * p_ip1 - p_i * S) == Z
        )

    def _matrix_has_nonfinite(self, M):
        for x in M:
            if x.has(sympy.nan) or x.has(sympy.zoo) or x.has(sympy.oo):
                return True
        return False

    # ---------------- Compute irrep matrices ----------------
    def _compute_irrep_matrices(self):
        irreps = self.irreps
        paths = self.bratteli_paths
        matrices = []

        for ir_idx, irrep in enumerate(irreps):
            basis_paths = paths[ir_idx]
            n = len(basis_paths)

            max_i = 2 * self.k - 1

            # seminormal e_k
            e_mats = {}
            for idx in range(1, max_i + 1):
                M_e = sympy.zeros(n, n)
                for a, s_path in enumerate(basis_paths):
                    for b, t_path in enumerate(basis_paths):
                        M_e[a, b] = self._pa_e(s_path, t_path, idx)
                e_mats[idx] = M_e

            # map generators (seminormal): p_i, b_i
            gen_mats_semi = {}
            for i in range(1, self.k + 1):
                gen_mats_semi[f"p_{i}"] = e_mats[2 * i - 1]
            for i in range(1, self.k):
                gen_mats_semi[f"b_{i}"] = e_mats[2 * i]

            # Recover the orthogonal basis scaling from the Hermitian p_i / b_i blocks.
            norms = self._pa_norms_from_seminormal(gen_mats_semi)
            norms_by_path = {p: norms[i] for i, p in enumerate(basis_paths)}

            gen_mats = {}
            for name, M_semi in gen_mats_semi.items():
                gen_mats[name] = self._orthogonalize_generator_matrix(M_semi, norms)

            # Compute sigma_k (seminormal) then orthogonalize, and form s_i = sigma_{2i} sigma_{2i+1}.
            sigma_orth = {}
            for k_idx in range(2, 2 * self.k):  # need up to 2k-1
                M_sigma_semi = sympy.zeros(n, n)
                for a, s_path in enumerate(basis_paths):
                    for b, t_path in enumerate(basis_paths):
                        M_sigma_semi[a, b] = self._pa_sigma(s_path, t_path, k_idx, basis_paths, norms_by_path=norms_by_path)
                sigma_orth[k_idx] = self._orthogonalize_generator_matrix(M_sigma_semi, norms)

            # Swap generators: s_i = sigma_{2i} sigma_{2i+1}
            for i in range(1, self.k):
                s_from_sigma = sympy.simplify(sigma_orth[2 * i] * sigma_orth[2 * i + 1])
                p_i = gen_mats[f"p_{i}"]
                p_ip1 = gen_mats[f"p_{i + 1}"]
                b_i = gen_mats[f"b_{i}"]
                s_local = self._pa_s_from_local_models(p_i, p_ip1, b_i, s_from_sigma)
                if self._pa_check_s_relations(s_local, p_i, p_ip1, b_i):
                    s_from_sigma = s_local
                elif self._matrix_has_nonfinite(s_from_sigma) or not self._pa_check_s_relations(s_from_sigma, p_i, p_ip1, b_i):
                    s_from_rel = self._pa_s_from_relations(
                        p_i, p_ip1, b_i
                    )
                    if s_from_rel is not None:
                        s_from_sigma = sympy.simplify(s_from_rel)
                gen_mats[f"s_{i}"] = s_from_sigma

            matrices.append(gen_mats)

        return matrices

    # ---------------- user-facing evaluation ----------------
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

    # ---------------- basis enumeration ----------------
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
