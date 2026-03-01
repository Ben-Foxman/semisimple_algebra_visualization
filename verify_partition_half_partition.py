"""
Verification checks for PartitionAlgebra and HalfPartitionAlgebra.

What this script checks (numerically, at user-chosen d):
  1) P_2(d) orthogonal-form generator matrices match the printed example (up to
     harmless square-root simplification choices).
  2) Diagram-level relations: b_i^2=b_i, p_i^2=d p_i (exact in the diagram basis).
  3) Branching “approximate norm preservation” identity for matrix units:
        ||E_{ij}^σ||^2 ≈ Σ_{ρ: σ ∈ Res(ρ)} ||E_{i+ρ, j+ρ}^ρ||^2
     where σ is an irrep of P_k(d) and ρ an irrep of P_{k+1/2}(d),
     with the identification i ↦ (i+ρ) determined by Bratteli path truncation.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import sympy as sp

from algebras.partition import PartitionAlgebra
from algebras.half_partition import HalfPartitionAlgebra


def _as_float(x) -> float:
    return float(sp.N(x, 40))


def _dict_star(A, coeffs: dict[int, sp.Expr]) -> dict[int, sp.Expr]:
    out: dict[int, sp.Expr] = {}
    for idx, c in coeffs.items():
        b = A.basis[idx]
        b_star = A._star(b)
        idx_star = A._basis_lookup[A._key(b_star)]
        out[idx_star] = out.get(idx_star, sp.S.Zero) + c
    return out


def _ensure_basis_traces(A) -> None:
    if getattr(A, "_basis_traces", None) is None:
        A._basis_traces = [A._trace_basis_index(j) for j in range(A.dim)]


def _trace(A, coeffs: dict[int, sp.Expr]) -> sp.Expr:
    _ensure_basis_traces(A)
    total = sp.S.Zero
    for idx, c in coeffs.items():
        total += c * A._basis_traces[idx]
    return sp.simplify(total)


def _inner(A, x: dict[int, sp.Expr], y: dict[int, sp.Expr]) -> sp.Expr:
    # <x,y> := tr(x * y^*)
    return _trace(A, A.multiply(x, _dict_star(A, y)))


def _norm2(A, x: dict[int, sp.Expr]) -> sp.Expr:
    return sp.simplify(_inner(A, x, x))


def _unit_offset(A, ir_idx: int) -> int:
    off = 0
    for r in range(ir_idx):
        d_r = len(A.bratteli_paths[r])
        off += d_r * d_r
    return off


def _irrep_basis_matrices(A, ir_idx: int):
    """
    Cache-friendly list of rho(e_a) for each diagram basis element e_a.
    Uses the cached word labels for diagrams.
    """
    labels = [A.label_of(a) for a in A.basis]
    return [A.irrep_matrix_for_label(ir_idx, lbl) for lbl in labels]


def _matrix_unit_coeffs_via_gram_solve(A, ir_idx: int, i: int, j: int) -> dict[int, sp.Expr]:
    """
    Compute E^rho_{ij} without building dual basis / full inverse:

      E^rho_{ij} = d_rho * sum_a rho(e_a^*)_{ij} e_a
    and rho(e_a^*)_{ij} is the (a)-entry of G^{-T} v, where v_a=rho(e_a)_{ij}.

    So solve (G^T) w = v, then coeffs are d_rho * w_a.
    """
    d_r = len(A.bratteli_paths[ir_idx])
    G = A.gram_matrix
    rho_basis = _irrep_basis_matrices(A, ir_idx)
    # Match the convention used in DiagramAlgebra._compute_matrix_units:
    # it uses rho(a^*)[j,i] as the coefficient contributing to E_{ij}.
    v = sp.Matrix([rho_basis[a][j, i] for a in range(A.dim)])
    w = G.T.LUsolve(v)
    out: dict[int, sp.Expr] = {}
    for a in range(A.dim):
        c = sp.simplify(sp.Integer(d_r) * w[a])
        if c != 0:
            out[a] = c
    return out


def _E(A, ir_idx: int, i: int, j: int) -> dict[int, sp.Expr]:
    return _matrix_unit_coeffs_via_gram_solve(A, ir_idx, i, j)


def check_diagram_relations_P2() -> None:
    d = sp.Symbol("d")
    A = PartitionAlgebra(2, d=d, symbolic_d=True)

    b1 = A._diagram_e(1)
    p1 = A._diagram_p(1)

    def as_dict(elem):
        return {A._basis_lookup[A._key(elem)]: sp.S.One}

    b1b1 = A.multiply(as_dict(b1), as_dict(b1))
    p1p1 = A.multiply(as_dict(p1), as_dict(p1))

    b1_idx = next(iter(as_dict(b1).keys()))
    p1_idx = next(iter(as_dict(p1).keys()))
    assert b1b1 == {b1_idx: sp.S.One}, "Expected b_1^2=b_1 in diagram basis."
    assert p1p1 == {p1_idx: d}, "Expected p_1^2=d p_1 in diagram basis."


def check_P2_example_numeric(d_val: float = 10.0, tol: float = 1e-10) -> None:
    d = sp.Symbol("d")
    A = PartitionAlgebra(2, d=d, symbolic_d=True)

    # Expected matrices from your lemma (same conventions).
    E = {}
    E[((), 2)] = {
        "p_1": sp.Matrix([[d, 0], [0, 0]]),
        "b_1": sp.Matrix([[sp.Rational(1, 1) / d, sp.sqrt(d - 1) / d], [sp.sqrt(d - 1) / d, (d - 1) / d]]),
        "p_2": sp.Matrix([[d, 0], [0, 0]]),
        "s_1": sp.eye(2),
    }
    E[((1,), 1)] = {
        "p_1": sp.Matrix([[d, 0, 0], [0, 0, 0], [0, 0, 0]]),
        "b_1": sp.Matrix([[sp.Rational(1, 1) / d, sp.sqrt(d - 1) / d, 0], [sp.sqrt(d - 1) / d, (d - 1) / d, 0], [0, 0, 0]]),
        "p_2": sp.Matrix([[0, 0, 0], [0, d / (d - 1), d * sp.sqrt(d - 2) / (d - 1)], [0, d * sp.sqrt(d - 2) / (d - 1), d * (d - 2) / (d - 1)]]),
        "s_1": sp.Matrix(
            [
                [0, 1 / sp.sqrt(d - 1), sp.sqrt(d - 2) / sp.sqrt(d - 1)],
                [1 / sp.sqrt(d - 1), (d - 2) / (d - 1), -sp.sqrt(d - 2) / (d - 1)],
                [sp.sqrt(d - 2) / sp.sqrt(d - 1), -sp.sqrt(d - 2) / (d - 1), 1 / (d - 1)],
            ]
        ),
    }
    E[((1, 1), 0)] = {"p_1": sp.Matrix([[0]]), "b_1": sp.Matrix([[0]]), "p_2": sp.Matrix([[0]]), "s_1": sp.Matrix([[-1]])}
    E[((2,), 0)] = {"p_1": sp.Matrix([[0]]), "b_1": sp.Matrix([[0]]), "p_2": sp.Matrix([[0]]), "s_1": sp.Matrix([[1]])}

    subs = {d: sp.nsimplify(d_val)}
    for ir_idx, ir in enumerate(A.irreps):
        if ir not in E:
            continue
        mats = A.irrep_matrices[ir_idx]
        for name, Mexp in E[ir].items():
            M = mats[name]
            diff = (M - Mexp).subs(subs).evalf(50)
            max_abs = max(abs(complex(diff[i, j])) for i in range(diff.rows) for j in range(diff.cols))
            assert max_abs <= tol, f"P2 example mismatch at irrep={ir} gen={name}: max|diff|={max_abs}"


@dataclass(frozen=True)
class BranchEmbedding:
    # For a fixed (sigma irrep, rho irrep) pair, map sigma basis index -> rho basis index.
    sigma_ir_idx: int
    rho_ir_idx: int
    sigma_to_rho: dict[int, int]


def _build_branch_embeddings(Pk: PartitionAlgebra, H: HalfPartitionAlgebra) -> dict[tuple[int, int], BranchEmbedding]:
    """
    Build embeddings σ ↪ ρ on bases via Bratteli path truncation:
      - σ is an irrep of P_k (top level 2k)
      - ρ is an irrep of P_{k+1/2} (top level 2k+1)
    For each ρ, each basis path T (length 2k+2) truncates to a P_k path t (length 2k+1).
    The predecessor node at level 2k determines which σ appears.
    """

    # Map from trunc_path (tuple of nodes) to (sigma_ir_idx, sigma_basis_idx)
    trunc_to_sigma = {}
    for s_ir_idx, paths in enumerate(Pk.bratteli_paths):
        for s_basis_idx, p in enumerate(paths):
            trunc_to_sigma[tuple(p)] = (s_ir_idx, s_basis_idx)

    out: dict[tuple[int, int], BranchEmbedding] = {}

    for r_ir_idx, rho_paths in enumerate(H.bratteli_paths):
        # For each rho basis vector, compute its truncation to length 2k+1
        for r_basis_idx, T in enumerate(rho_paths):
            trunc = tuple(T[:-1])  # drop the last node (level 2k+1)
            if trunc not in trunc_to_sigma:
                # If towers are inconsistent, we can't compare.
                raise RuntimeError("Half-partition path truncation did not match any P_k path.")
            s_ir_idx, s_basis_idx = trunc_to_sigma[trunc]
            key = (s_ir_idx, r_ir_idx)
            if key not in out:
                out[key] = BranchEmbedding(sigma_ir_idx=s_ir_idx, rho_ir_idx=r_ir_idx, sigma_to_rho={})
            out[key].sigma_to_rho[s_basis_idx] = r_basis_idx

    return out


def check_branching_norm_identity(k: int = 2, d_val: float = 200.0, tol_rel: float = 5e-3, max_checks_per_sigma: int = 6) -> None:
    """
    Numerical check of
      ||E^σ_{ij}||^2 ≈ Σ_ρ ||E^ρ_{I,J}||^2
    where (I,J) is the embedding of (i,j) inside ρ determined by Bratteli truncation.
    """

    # Use exact integer d to avoid numerical instability in Gram solves.
    d_int = int(d_val)
    Pk = PartitionAlgebra(k, d=sp.Integer(d_int), symbolic_d=False)
    H = HalfPartitionAlgebra(k, d=sp.Integer(d_int), symbolic_d=False)  # P_{k+1/2}

    embeddings = _build_branch_embeddings(Pk, H)

    # Group rho's by sigma
    rhos_by_sigma: dict[int, list[BranchEmbedding]] = {}
    for (s_ir_idx, _r_ir_idx), emb in embeddings.items():
        rhos_by_sigma.setdefault(s_ir_idx, []).append(emb)

    # For each sigma, check a handful of (i,j) pairs.
    for s_ir_idx, sigma_paths in enumerate(Pk.bratteli_paths):
        d_sigma = len(sigma_paths)
        if d_sigma == 0:
            continue
        rho_embs = rhos_by_sigma.get(s_ir_idx, [])
        if not rho_embs:
            raise AssertionError(f"No half-partition irreps contained sigma irrep index {s_ir_idx}.")

        # Choose some test pairs (i,j). Prefer diagonal and a couple off-diagonal if possible.
        pairs = []
        pairs.append((0, 0))
        if d_sigma > 1:
            pairs.append((1, 1))
            pairs.append((0, 1))
            pairs.append((1, 0))
        if d_sigma > 2:
            pairs.append((2, 2))
            pairs.append((0, 2))
        pairs = pairs[:max_checks_per_sigma]

        for (i, j) in pairs:
            lhs = _as_float(_norm2(Pk, _E(Pk, s_ir_idx, i, j)))
            rhs = 0.0
            for emb in rho_embs:
                I = emb.sigma_to_rho.get(i)
                J = emb.sigma_to_rho.get(j)
                if I is None or J is None:
                    # This rho doesn't contain the (i,j) basis vectors from sigma.
                    continue
                rhs += _as_float(_norm2(H, _E(H, emb.rho_ir_idx, I, J)))

            if lhs == 0.0 and rhs == 0.0:
                continue
            denom = max(1e-12, abs(lhs))
            rel = abs(lhs - rhs) / denom
            assert rel <= tol_rel, (
                f"Branching norm check failed at k={k}, d={d_val}, "
                f"sigma_ir={s_ir_idx}, (i,j)=({i},{j}): lhs={lhs}, rhs={rhs}, rel={rel}"
            )


def main() -> None:
    # 1) Exact diagram relations at k=2 (symbolic d).
    check_diagram_relations_P2()

    # 2) P2 printed example, numeric spot check.
    check_P2_example_numeric(d_val=10.0, tol=1e-10)

    # 3) Branching norm identity (numeric, large d).
    # Start with k=2; k=3 is substantially heavier.
    check_branching_norm_identity(k=2, d_val=200.0, tol_rel=7e-3)

    print("All checks passed.")


if __name__ == "__main__":
    main()

