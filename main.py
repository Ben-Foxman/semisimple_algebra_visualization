"""
Command-line script for diagram algebras (Partition, Brauer).

Usage: python main.py <algebra> [args] [--print-matrix]

  algebra: "partition", "brauer", "walled_brauer", "symmetric"
  args:
    partition/brauer: k (positive integer)
    walled_brauer:    r s (non-negative integers, r + s = k)
    symmetric:        k (permutation degree)
"""

import argparse
from algebras import (
    PartitionAlgebra,
    BrauerAlgebra,
    WalledBrauerAlgebra,
    SymmetricGroupAlgebra,
)


ALGEBRAS = {
    "partition": PartitionAlgebra,
    "brauer": BrauerAlgebra,
    "walled_brauer": WalledBrauerAlgebra,
    "symmetric": SymmetricGroupAlgebra,
}


def format_coeff(coeff):
    """Format a symbolic coefficient for display."""
    s = str(coeff)
    if s == "1":
        return ""
    if s == "-1":
        return "-"
    return f"({s}) * "


def format_dual_element(alg, i):
    """Format dual basis element e_i* as a linear combination."""
    dual = alg.dual_basis[i]
    terms = []
    for j, c in sorted(dual.items()):
        try:
            if hasattr(c, "is_zero") and c.is_zero():
                continue
        except Exception:
            if c == 0:
                continue
        basis_label = alg.label_of(alg.basis[j])
        terms.append(f"{format_coeff(c)}{basis_label}")
    if not terms:
        return "0"
    return " + ".join(terms).replace(" + -", " - ")


def main():
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument(
        "--print-matrix",
        action="store_true",
        help="print the Gram matrix",
    )
    common.add_argument(
        "--d",
        default="5",
        help="fixed value for d (ignored if --symbolic-d)",
    )
    common.add_argument(
        "--symbolic-d",
        action="store_true",
        help="use symbolic d instead of a fixed numeric value",
    )
    parser = argparse.ArgumentParser(
        description="Diagram algebras with numeric d (default 5) or symbolic d.",
        parents=[common],
    )
    subparsers = parser.add_subparsers(dest="algebra", required=True)

    part = subparsers.add_parser("partition", parents=[common])
    part.add_argument("k", type=int, help="positive integer (points per side)")
    part.set_defaults(
        factory=lambda args: PartitionAlgebra(
            args.k, d=args.d, symbolic_d=args.symbolic_d
        ),
        label="Partition",
    )

    brauer = subparsers.add_parser("brauer", parents=[common])
    brauer.add_argument("k", type=int, help="positive integer (points per side)")
    brauer.set_defaults(
        factory=lambda args: BrauerAlgebra(
            args.k, d=args.d, symbolic_d=args.symbolic_d
        ),
        label="Brauer",
    )

    walled = subparsers.add_parser("walled_brauer", parents=[common])
    walled.add_argument("r", type=int, help="left points")
    walled.add_argument("s", type=int, help="right points")
    walled.set_defaults(
        factory=lambda args: WalledBrauerAlgebra(
            args.r, args.s, d=args.d, symbolic_d=args.symbolic_d
        ),
        label="Walled Brauer",
    )

    symm = subparsers.add_parser("symmetric", parents=[common])
    symm.add_argument("k", type=int, help="permutation degree")
    symm.set_defaults(
        factory=lambda args: SymmetricGroupAlgebra(
            args.k, d=args.d, symbolic_d=args.symbolic_d
        ),
        label="Symmetric",
    )

    args = parser.parse_args()

    if getattr(args, "k", None) is not None and args.k < 1:
        parser.error("k must be at least 1.")
    if getattr(args, "r", None) is not None and args.r < 0:
        parser.error("r must be non-negative.")
    if getattr(args, "s", None) is not None and args.s < 0:
        parser.error("s must be non-negative.")
    if getattr(args, "r", None) is not None and getattr(args, "s", None) is not None:
        if args.r + args.s < 1:
            parser.error("r + s must be at least 1.")

    if args.algebra == "walled_brauer":
        print(f"Initializing {args.label} algebra (r={args.r}, s={args.s})...")
    else:
        print(f"Initializing {args.label} algebra (k={getattr(args, 'k', None)})...")
    alg = args.factory(args)
    print(f"Dimension: {alg.dim}\n")

    # Basis elements
    print("=" * 70)
    print("BASIS ELEMENTS")
    print("=" * 70)
    for i in range(alg.dim):
        label = alg.label_of(alg.basis[i])
        print(f"  e_{i}:  {label}")
    print()

    # Dual basis elements
    print("=" * 70)
    print("DUAL BASIS ELEMENTS")
    print("=" * 70)
    for i in range(alg.dim):
        combo = format_dual_element(alg, i)
        label = alg.label_of(alg.basis[i])
        print(f"  {label}* = {combo}")
    print("=" * 70)

    if args.print_matrix:
        print("\n" + "=" * 70)
        print("GRAM MATRIX")
        print("=" * 70)
        print(alg.gram_matrix)
        print("=" * 70)


if __name__ == "__main__":
    main()
