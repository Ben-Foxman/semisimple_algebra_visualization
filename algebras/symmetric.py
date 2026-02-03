from .base import *

class SymmetricGroupAlgebra(DiagramAlgebra):
    """Group algebra of S_k with basis indexed by permutations."""

    algebra_id = "symmetric_group"

    def _generate_basis(self):
        return [tuple(p) for p in itertools.permutations(range(1, self.k + 1))]

    def _key(self, elem):
        return elem

    def _identity_element(self):
        return tuple(range(1, self.k + 1))

    def _star(self, elem):
        inv = [0] * self.k
        for i, val in enumerate(elem, start=1):
            inv[val - 1] = i
        return tuple(inv)

    def _compose(self, p1, p2):
        # permutation composition: p1 after p2
        comp = tuple(p1[i - 1] for i in p2)
        return comp, 0

    def _serialize_basis_element(self, elem):
        return list(elem)

    def _deserialize_basis_element(self, data):
        return tuple(data)

    def _generator_diagrams(self):
        gens = []
        for i in range(1, self.k):
            perm = list(range(1, self.k + 1))
            perm[i - 1], perm[i] = perm[i], perm[i - 1]
            gens.append((f"s_{i}", tuple(perm)))
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

        rec(0, tuple(), [tuple()])
        irreps = []
        paths_by_irrep = []
        for part in _partitions_of(k):
            if part in all_paths:
                irreps.append(part)
                paths_by_irrep.append(all_paths[part])
        return irreps, paths_by_irrep

    def _tableau_positions(self, path):
        positions = {}
        for step in range(1, len(path)):
            prev = path[step - 1]
            curr = path[step]
            if _partition_size(curr) != _partition_size(prev) + 1:
                raise ValueError("Non-additive step in symmetric path.")
            # find added box
            row = len(curr)
            col = 1
            for i in range(len(curr)):
                prev_len = prev[i] if i < len(prev) else 0
                if curr[i] == prev_len + 1:
                    row = i + 1
                    col = curr[i]
                    break
            positions[step] = (row, col)
        return positions

    def _path_from_positions(self, positions):
        max_row = max((r for r, _ in positions.values()), default=0)
        counts = [0] * max_row
        path = [tuple()]
        for t in range(1, self.k + 1):
            r, _ = positions[t]
            if r > len(counts):
                counts.extend([0] * (r - len(counts)))
            counts[r - 1] += 1
            shape = tuple(c for c in counts if c > 0)
            path.append(shape)
        return tuple(path)

    def _compute_irrep_matrices(self):
        irreps = self.irreps
        paths = self.bratteli_paths
        matrices = []
        for ir_idx, irrep in enumerate(irreps):
            basis_paths = paths[ir_idx]
            index = {tuple(p): i for i, p in enumerate(basis_paths)}
            n = len(basis_paths)
            gen_mats = {}
            for i in range(1, self.k):
                M = sympy.zeros(n, n)
                for p_idx, p in enumerate(basis_paths):
                    pos = self._tableau_positions(p)
                    r1, c1 = pos[i]
                    r2, c2 = pos[i + 1]
                    if r1 == r2:
                        M[p_idx, p_idx] = 1
                        continue
                    if c1 == c2:
                        M[p_idx, p_idx] = -1
                        continue
                    # swap i and i+1 labels in tableau to get p'
                    pos2 = dict(pos)
                    pos2[i], pos2[i + 1] = pos2[i + 1], pos2[i]
                    p2 = self._path_from_positions(pos2)
                    if p2 not in index:
                        raise ValueError("Invalid swapped tableau.")
                    q_idx = index[p2]
                    c_i = c1 - r1
                    c_j = c2 - r2
                    a = sympy.Rational(1, c_j - c_i)
                    b = sympy.sqrt(1 - a * a)
                    M[p_idx, p_idx] = a
                    M[p_idx, q_idx] = b
                    M[q_idx, p_idx] = b
                    M[q_idx, q_idx] = -a
                gen_mats[f"s_{i}"] = M
            matrices.append(gen_mats)
        return matrices

    def irrep_matrix_for_label(self, irrep_idx, label):
        if label == "1":
            n = len(self.bratteli_paths[irrep_idx])
            return sympy.eye(n)
        if not label:
            n = len(self.bratteli_paths[irrep_idx])
            return sympy.eye(n)
        mats = self.irrep_matrices[irrep_idx]
        tokens = label.split()
        result = None
        for tok in tokens:
            if tok not in mats:
                raise ValueError(f"Unknown generator {tok} for symmetric group.")
            if result is None:
                result = mats[tok]
            else:
                result = result * mats[tok]
        return result

    def _compose(self, p1, p2):
        # permutation composition: p1 ∘ p2
        return tuple(p1[i - 1] for i in p2), 0
