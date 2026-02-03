import os
import sys
import sympy

sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from algebras import PartitionAlgebra


def _element_add(a, b):
    out = dict(a)
    for idx, coeff in b.items():
        out[idx] = out.get(idx, sympy.S.Zero) + coeff
        if out[idx] == 0:
            out.pop(idx)
    return out


def _element_sub(a, b):
    return _element_add(a, {idx: -coeff for idx, coeff in b.items()})


def _element_scale(a, scalar):
    if scalar == 0:
        return {}
    return {idx: scalar * coeff for idx, coeff in a.items() if coeff != 0}


def _element_mul(algebra, a, b):
    return algebra.multiply(a, b)


def _element_mul_many(algebra, elems):
    if not elems:
        return {}
    out = elems[0]
    for elem in elems[1:]:
        out = _element_mul(algebra, out, elem)
    return out


def _element_identity(algebra):
    idx = algebra._basis_lookup[algebra._key(algebra._identity_element())]
    return {idx: sympy.S.One}


def _element_from_diagram(algebra, diagram):
    idx = algebra._basis_lookup[algebra._key(diagram)]
    return {idx: sympy.S.One}


def _element_matrix(algebra, elem):
    n = algebra.dim
    mat = sympy.zeros(n, n)
    for j in range(n):
        prod = algebra.multiply(elem, {j: sympy.S.One})
        for idx, coeff in prod.items():
            mat[idx, j] = coeff
    return mat


def _element_e(algebra, m):
    # e_{2i-1} -> p_i, e_{2i} -> b_i
    if m % 2 == 1:
        return _element_from_diagram(algebra, algebra._diagram_p((m + 1) // 2))
    return _element_from_diagram(algebra, algebra._diagram_e(m // 2))


def _element_s(algebra, i):
    return _element_from_diagram(algebra, algebra._diagram_s(i))


def build_jm_elements(algebra):
    # Based on Enyang's recursive definitions in Section 2.3.
    L = {}
    sigma = {}

    L[1] = {}
    L[2] = _element_e(algebra, 1)
    sigma[2] = _element_identity(algebra)
    sigma[3] = _element_s(algebra, 1)

    max_i = 2 * algebra.k
    if max_i >= 3:
        e_2 = _element_e(algebra, 2)
        term1 = _element_mul_many(algebra, [L[2], e_2])
        term2 = _element_mul_many(algebra, [e_2, L[2]])
        term3 = _element_mul_many(algebra, [_element_scale(_element_identity(algebra), algebra.d), e_2])
        L[3] = _element_add(
            _element_add(_element_scale(term1, -1), _element_scale(term2, -1)),
            _element_add(term3, sigma[2]),
        )
    if max_i >= 4:
        e_2 = _element_e(algebra, 2)
        e_3 = _element_e(algebra, 3)
        s_1 = _element_s(algebra, 1)
        term1 = _element_mul_many(algebra, [s_1, L[2], e_2])
        term2 = _element_mul_many(algebra, [e_2, L[2], s_1])
        term3 = _element_mul_many(algebra, [e_2, L[2], e_3, e_2])
        term4 = _element_mul_many(algebra, [s_1, L[2], s_1])
        L[4] = _element_add(
            _element_add(_element_add(_element_scale(term1, -1), _element_scale(term2, -1)), term3),
            _element_add(term4, sigma[3]),
        )
    for i in range(2, algebra.k):
        s_i = _element_s(algebra, i)
        e_2i = _element_e(algebra, 2 * i)
        e_2i_minus_1 = _element_e(algebra, 2 * i - 1)
        e_2i_minus_2 = _element_e(algebra, 2 * i - 2)
        e_2i_plus_1 = _element_e(algebra, 2 * i + 1)

        # sigma_{2i+1}
        if 2 * i + 1 not in sigma:
            term1 = _element_mul_many(algebra, [s_i, _element_s(algebra, i - 1), sigma[2 * i - 1],
                                                _element_s(algebra, i), _element_s(algebra, i - 1)])
            term2 = _element_mul_many(algebra, [s_i, e_2i_minus_2, L[2 * i - 2], s_i, e_2i_minus_2, s_i])
            term3 = _element_mul_many(algebra, [e_2i_minus_2, L[2 * i - 2], s_i, e_2i_minus_2])
            term4 = _element_mul_many(algebra, [s_i, e_2i_minus_2, L[2 * i - 2], _element_s(algebra, i - 1),
                                                e_2i, e_2i_minus_1, e_2i_minus_2])
            term5 = _element_mul_many(algebra, [e_2i_minus_2, e_2i_minus_1, e_2i, _element_s(algebra, i - 1),
                                                L[2 * i - 2], e_2i_minus_2, s_i])
            sigma[2 * i + 1] = _element_sub(
                _element_add(_element_add(term1, term2), term3),
                _element_add(term4, term5),
            )

        # sigma_{2i}
        if 2 * i not in sigma:
            term1 = _element_mul_many(algebra, [_element_s(algebra, i - 1), s_i, sigma[2 * i - 2],
                                                s_i, _element_s(algebra, i - 1)])
            term2 = _element_mul_many(algebra, [e_2i_minus_2, L[2 * i - 2], s_i, e_2i_minus_2, s_i])
            term3 = _element_mul_many(algebra, [s_i, e_2i_minus_2, L[2 * i - 2], s_i, e_2i_minus_2])
            term4 = _element_mul_many(algebra, [e_2i_minus_2, L[2 * i - 2], _element_s(algebra, i - 1),
                                                e_2i, e_2i_minus_1, e_2i_minus_2])
            term5 = _element_mul_many(algebra, [s_i, e_2i_minus_2, e_2i_minus_1, e_2i,
                                                _element_s(algebra, i - 1), L[2 * i - 2],
                                                e_2i_minus_2, s_i])
            sigma[2 * i] = _element_sub(
                _element_add(_element_add(term1, term2), term3),
                _element_add(term4, term5),
            )

        # L_{2i+1}
        if 2 * i + 1 <= max_i:
            term1 = _element_mul_many(algebra, [L[2 * i], e_2i])
            term2 = _element_mul_many(algebra, [e_2i, L[2 * i]])
            term3 = _element_mul_many(algebra, [_element_sub(_element_scale(_element_identity(algebra), algebra.d),
                                                            L[2 * i - 1]), e_2i])
            term4 = _element_mul_many(algebra, [s_i, L[2 * i - 1], s_i])
            L[2 * i + 1] = _element_add(
                _element_add(_element_add(_element_scale(term1, -1), _element_scale(term2, -1)), term3),
                _element_add(term4, sigma[2 * i]),
            )

        # L_{2i+2}
        if 2 * i + 2 <= max_i:
            term1 = _element_mul_many(algebra, [s_i, L[2 * i], e_2i])
            term2 = _element_mul_many(algebra, [e_2i, L[2 * i], s_i])
            term3 = _element_mul_many(algebra, [e_2i, L[2 * i], e_2i_plus_1, e_2i])
            term4 = _element_mul_many(algebra, [s_i, L[2 * i], s_i])
            L[2 * i + 2] = _element_add(
                _element_add(_element_add(_element_scale(term1, -1), _element_scale(term2, -1)), term3),
                _element_add(term4, sigma[2 * i + 1]),
            )

    return L, sigma


def joint_eigenvector_for_path(algebra, L_mats, path):
    # Solve (L_i - c_i I) v = 0 for all i.
    n = algebra.dim
    mats = []
    for i in range(1, 2 * algebra.k + 1):
        if i not in L_mats:
            continue
        c_i = algebra._pa_content(path, i) if i < len(path) else sympy.S.Zero
        mats.append(L_mats[i] - c_i * sympy.eye(n))
    if not mats:
        return None
    stacked = sympy.Matrix.vstack(*mats)
    null = stacked.nullspace()
    if not null:
        return None
    # Use the first basis vector.
    v = null[0]
    return v


def eigenbasis_by_irrep(algebra, L_mats):
    basis_by_irrep = []
    for paths in algebra.bratteli_paths:
        eigenvectors = []
        for path in paths:
            v = joint_eigenvector_for_path(algebra, L_mats, path)
            if v is None:
                eigenvectors.append(None)
            else:
                eigenvectors.append(v)
        basis_by_irrep.append(eigenvectors)
    return basis_by_irrep


def generator_matrices_in_eigenbasis(algebra, eigenvectors):
    gen_mats = []
    for paths, vecs in zip(algebra.bratteli_paths, eigenvectors):
        if any(v is None for v in vecs):
            gen_mats.append(None)
            continue
        B = sympy.Matrix.hstack(*vecs)
        if B.rank() != B.cols:
            gen_mats.append(None)
            continue
        mats = {}
        gen_elems = {}
        for i in range(1, algebra.k + 1):
            gen_elems[f"p_{i}"] = _element_matrix(algebra, _element_e(algebra, 2 * i - 1))
        for i in range(1, algebra.k):
            gen_elems[f"b_{i}"] = _element_matrix(algebra, _element_e(algebra, 2 * i))
            gen_elems[f"s_{i}"] = _element_matrix(algebra, _element_s(algebra, i))
        for label, G in gen_elems.items():
            cols = []
            for j in range(B.cols):
                rhs = G * B[:, j]
                coeffs, _ = B.gauss_jordan_solve(rhs)
                cols.append(coeffs)
            mats[label] = sympy.Matrix.hstack(*cols)
        gen_mats.append(mats)
    return gen_mats


def compare_with_formula(algebra, eigen_mats):
    results = []
    for idx, mats in enumerate(eigen_mats):
        if mats is None:
            results.append((algebra.irreps[idx], "missing_eigenbasis"))
            continue
        formula = algebra.irrep_matrices[idx]
        for key, mat in mats.items():
            if key not in formula:
                results.append((algebra.irreps[idx], key, "missing_formula"))
                continue
            diff = (mat - formula[key]).applyfunc(sympy.simplify)
            if diff != sympy.zeros(diff.rows, diff.cols):
                results.append((algebra.irreps[idx], key, "mismatch"))
    return results


def run(k=2, d=5):
    alg = PartitionAlgebra(k, d=d, symbolic_d=False)
    L, _ = build_jm_elements(alg)
    L_mats = {i: _element_matrix(alg, elem) for i, elem in L.items()}
    eigenvectors = eigenbasis_by_irrep(alg, L_mats)
    eigen_mats = generator_matrices_in_eigenbasis(alg, eigenvectors)
    mismatches = compare_with_formula(alg, eigen_mats)
    return mismatches


if __name__ == "__main__":
    issues = run(k=2, d=5)
    print("mismatches:", issues if issues else "none")
