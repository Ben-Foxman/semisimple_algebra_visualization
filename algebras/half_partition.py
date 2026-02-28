from .base import *
import sympy
from collections import defaultdict

from .partition import PartitionAlgebra  # reuse all the Enyang machinery


class HalfPartitionAlgebra(PartitionAlgebra):
    """
    Half-partition algebra P_{k+1/2}(d).

    Implementation model:
      - diagrams live on columns 1..(k+1) (so renderer works with self.k = k+1)
      - basis is partitions of ±[k+1] with the constraint {k+1, -(k+1)} in the same block
      - Bratteli tower stops at level 2k+1 (so e-indices are 1..2k)
      - generators: s_1..s_{k-1}, b_1..b_k, p_1..p_k
    """

    algebra_id = "half_partition"

    def __init__(self, k: int, d=None, symbolic_d=False):
        # Here k means the integer part in P_{k+1/2}.
        # We store (k+1) columns so the renderer/table logic stays unchanged.
        self.k_int = int(k)
        super().__init__(self.k_int + 1, d=d, symbolic_d=symbolic_d)

    # ---------------- generators ----------------
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
        # basis elements are set partitions of ±[k+1] with {k+1, -(k+1)} in same block
        ncols = self.k_int + 1
        elts = list(range(1, ncols + 1)) + [-i for i in range(1, ncols + 1)]
        all_parts = self._all_partitions(elts)

        fixed_a = ncols
        fixed_b = -ncols
        out = []
        for P in all_parts:
            ok = False
            for blk in P.blocks:
                if fixed_a in blk and fixed_b in blk:
                    ok = True
                    break
            if ok:
                out.append(P)
        return out

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
                all_paths[path[-1]].append(path[:])
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

        # At level 2k+1, floor(level/2)=k, so l = k - |λ|, hence |λ| <= k.
        for l in range(k + 1):
            size = k - l
            for part in _partitions_of(size):
                key = (part, l)
                if key in all_paths:
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
           s_i computed for i<=k-1
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

            # generator mats (seminormal)
            gen_mats_semi = {}
            for i in range(1, self.k_int + 1):
                gen_mats_semi[f"p_{i}"] = e_mats[2 * i - 1]
                gen_mats_semi[f"b_{i}"] = e_mats[2 * i]

            # orthogonalize from p_i,b_i
            norms = self._pa_norms_from_seminormal(gen_mats_semi)
            gen_mats = {}
            for name, M_semi in gen_mats_semi.items():
                gen_mats[name] = self._orthogonalize_generator_matrix(M_semi, norms)

            # Solve for s_i (i=1..k-1) using relations with p_i,p_{i+1},b_i
            for i in range(1, self.k_int):
                S = self._pa_s_from_relations(
                    gen_mats[f"p_{i}"], gen_mats[f"p_{i + 1}"], gen_mats[f"b_{i}"]
                )
                if S is None:
                    S = sympy.eye(n)
                gen_mats[f"s_{i}"] = sympy.simplify(S)

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