from .base import *
import sympy
from collections import defaultdict

from .partition import PartitionAlgebra  # reuse all the Enyang machinery


class HalfPartitionAlgebra(PartitionAlgebra):
    """
    Half-partition algebra P_{k+1/2}(d).

    Implementation model (IMPORTANT):
      - This algebra has a diagram basis indexed by set partitions of 2k+1 points.
      - We model it using integers {1,...,k+1} (top) and {-1,...,-k} (bottom),
        so the distinguished point (k+1) appears ONLY ONCE (it is "shared").
      - Composition glues (-i) from the left factor to (+i) of the right factor
        for i=1..k, and also identifies the distinguished points (+k+1) across factors.
      - Bratteli tower stops at level 2k+1 (so e-indices are 1..2k)
      - generators: s_1..s_{k-1}, b_1..b_k, p_1..p_k
    """

    algebra_id = "half_partition"

    def __init__(self, k: int, d=None, symbolic_d=False):
        # Here k means the integer part in P_{k+1/2}.
        # We still store self.k = k+1 so the renderer spacing stays reasonable,
        # but the *basis elements* live on {1..k+1} ∪ {-1..-k}.
        self.k_int = int(k)
        super().__init__(self.k_int + 1, d=d, symbolic_d=symbolic_d)

    # ---------------- identity / involution ----------------
    def _identity_element(self):
        # Identity: i is paired to -i for i=1..k, and the distinguished top point (k+1)
        # is a singleton.
        blocks = [{i, -i} for i in range(1, self.k_int + 1)]
        blocks.append({self.k_int + 1})
        return SetPartition(blocks)

    # ---------------- generators ----------------
    def _diagram_s(self, i):
        # swap generator for i=1..k-1
        blocks = [{i, -(i + 1)}, {i + 1, -i}]
        for j in range(1, self.k_int + 1):
            if j in (i, i + 1):
                continue
            blocks.append({j, -j})
        blocks.append({self.k_int + 1})
        return SetPartition(blocks)

    def _diagram_p(self, i):
        # point generator for i=1..k
        blocks = [{i}, {-i}]
        for j in range(1, self.k_int + 1):
            if j == i:
                continue
            blocks.append({j, -j})
        blocks.append({self.k_int + 1})
        return SetPartition(blocks)

    def _diagram_e(self, i):
        # bridge generators b_i for i=1..k.
        #
        # For i<k, this is the usual 4-block {i,i+1,-i,-(i+1)}.
        # For i=k, there is no bottom -(k+1), so we use the 3-block {k, k+1, -k}
        # connecting the last ordinary strand to the distinguished point.
        blocks = []
        if i < self.k_int:
            blocks.append({i, i + 1, -i, -(i + 1)})
        else:
            blocks.append({self.k_int, self.k_int + 1, -self.k_int})
        for j in range(1, self.k_int + 1):
            if j in (i, i + 1) or (i == self.k_int and j == self.k_int):
                continue
            blocks.append({j, -j})
        # distinguished point is already in the big block if i==k; otherwise singleton
        if i < self.k_int:
            blocks.append({self.k_int + 1})
        return SetPartition(blocks)

    def _generator_diagrams(self):
        gens = []
        # s_i for i=1..k-1 (do not swap with last fixed column)
        for i in range(1, self.k_int):
            gens.append((f"s_{i}", self._diagram_s(i)))
        # b_i for i=1..k
        for i in range(1, self.k_int + 1):
            gens.append((f"b_{i}", self._diagram_e(i)))
        # p_i for i=1..k  (no p_{k+1})
        for i in range(1, self.k_int + 1):
            gens.append((f"p_{i}", self._diagram_p(i)))
        return gens

    # ---------------- basis enumeration ----------------
    def _generate_basis(self):
        # basis elements are set partitions of {1..k+1} ∪ {-1..-k}
        top = list(range(1, self.k_int + 2))          # 1..k+1
        bottom = [-i for i in range(1, self.k_int + 1)]  # -1..-k
        elts = top + bottom
        return self._all_partitions(elts)

    # ---------------- composition (diagram multiplication) ----------------
    def _compose(self, p1, p2):
        """
        Compose two half-partition diagrams.

        Boundary of result:
          - top:   1..k+1 from the LEFT factor (p1)
          - bottom: -1..-k from the RIGHT factor (p2)

        Gluing:
          - glue (-i) from p1 to (+i) from p2 for i=1..k
          - identify the distinguished points (+k+1) across p1 and p2
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

        # glue the ordinary bottom nodes
        for i in range(1, self.k_int + 1):
            add_edge(("p1", -i), ("p2", i))

        # identify the distinguished point across the two factors
        add_edge(("p1", self.k_int + 1), ("p2", self.k_int + 1))

        visited = set()
        blocks = []
        loops = 0

        def is_boundary(node):
            side, val = node
            if side == "p1" and val > 0:
                return True
            if side == "p2" and val < 0:
                return True
            return False

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

    # ---------------- irreps / paths ----------------
    def _compute_irreps_and_paths(self):
        """
        Same recursion as PartitionAlgebra, but stop at level 2k+1 (k=self.k_int),
        and only keep irreps that appear at that top level.

        Node convention unchanged: node=(partition, l) with l = floor(level/2) - |partition|.
        """
        k = self.k_int
        top_level = 2 * k + 1

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
            if level == top_level:
                all_paths[path[-1]].append(tuple(path))
                return
            if level % 2 == 0:
                # even -> odd: stay or remove
                if valid_node(level + 1, part):
                    rec(level + 1, part, path + [node_for(level + 1, part)])
                for nxt in _removable_partitions(part):
                    if valid_node(level + 1, nxt):
                        rec(level + 1, nxt, path + [node_for(level + 1, nxt)])
            else:
                # odd -> even: stay or add
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

    # ---------------- Compute irrep matrices ----------------
    def _compute_irrep_matrices(self):
        """
        Same as PartitionAlgebra, but max e-index is 2k (not 2(k+1)-1).
        That matches generators:
           p_i = e_{2i-1} for i<=k
           b_i = e_{2i}   for i<=k
           s_i = sigma_{2i} sigma_{2i+1} for i<=k-1
        """
        irreps = self.irreps
        paths = self.bratteli_paths
        matrices = []

        max_i = 2 * self.k_int  # NOTE: differs from PartitionAlgebra

        for ir_idx, irrep in enumerate(irreps):
            basis_paths = paths[ir_idx]
            n = len(basis_paths)

            # seminormal e_k for k=1..2k
            e_mats = {}
            for idx in range(1, max_i + 1):
                M_e = sympy.zeros(n, n)
                for a, s_path in enumerate(basis_paths):
                    for b, t_path in enumerate(basis_paths):
                        M_e[a, b] = self._pa_e(s_path, t_path, idx)
                e_mats[idx] = M_e

            # generator mats (seminormal): p_i, b_i
            gen_mats_semi = {}
            for i in range(1, self.k_int + 1):
                gen_mats_semi[f"p_{i}"] = e_mats[2 * i - 1]
                gen_mats_semi[f"b_{i}"] = e_mats[2 * i]

            # Recover the orthogonal basis scaling from the Hermitian p_i / b_i blocks.
            norms = self._pa_norms_from_seminormal(gen_mats_semi)
            norms_by_path = {p: norms[i] for i, p in enumerate(basis_paths)}

            gen_mats = {}
            for name, M_semi in gen_mats_semi.items():
                gen_mats[name] = self._orthogonalize_generator_matrix(M_semi, norms)

            # Compute sigma_k (seminormal) then orthogonalize, and form s_i = sigma_{2i} sigma_{2i+1}.
            # Need sigma indices up to 2k-1 (since s_{k-1} uses sigma_{2k-2}, sigma_{2k-1}).
            sigma_orth = {}
            for k_idx in range(2, 2 * self.k_int):  # up to 2k-1
                M_sigma_semi = sympy.zeros(n, n)
                for a, s_path in enumerate(basis_paths):
                    for b, t_path in enumerate(basis_paths):
                        M_sigma_semi[a, b] = self._pa_sigma(
                            s_path, t_path, k_idx, basis_paths, norms_by_path=norms_by_path
                        )
                sigma_orth[k_idx] = self._orthogonalize_generator_matrix(M_sigma_semi, norms)

            for i in range(1, self.k_int):
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
                raise ValueError(f"Unknown generator {tok} for half-partition algebra.")
            result = mats[tok] if result is None else (result * mats[tok])
        return result
