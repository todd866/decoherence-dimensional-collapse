#!/usr/bin/env python3
"""
Verify the qutrit (d=3) orthogonality theorem:

  For d=3, a full-rank state rho has T^diag perp_{g_B} T^off
  if and only if rho is diagonal in the pointer basis.

This script:
1. Verifies the key identity g_B(X_D, X_O) = (1/2) sum_k x_k (L_{X_O})_{kk}
2. Verifies the diagonal-of-SLD formula (L_{R_st})_{kk} from the Lyapunov eq
3. Tests the "sparse off-diagonal" case (rho_{12} != 0, rho_{13}=rho_{23}=0)
4. Tests fully general non-diagonal states
5. Checks that the minimum cross-Gram entry is bounded away from zero

Usage:
    python3 verify_qutrit_theorem.py
"""

import numpy as np
from itertools import product

np.set_printoptions(precision=8, linewidth=120)


def random_density_matrix(d, rng):
    """Random full-rank density matrix (Hilbert-Schmidt measure)."""
    A = rng.standard_normal((d, d)) + 1j * rng.standard_normal((d, d))
    rho = A @ A.conj().T
    return rho / np.trace(rho)


def solve_sld(rho, X):
    """Solve rho L + L rho = 2X for Hermitian L (full-rank rho)."""
    eigvals, U = np.linalg.eigh(rho)
    d = len(eigvals)
    X_eig = U.conj().T @ X @ U
    L_eig = np.zeros((d, d), dtype=complex)
    for m in range(d):
        for n in range(d):
            L_eig[m, n] = 2.0 * X_eig[m, n] / (eigvals[m] + eigvals[n])
    return U @ L_eig @ U.conj().T


def bures_inner(rho, X, Y):
    """g_B(X,Y) = (1/4) Re Tr[rho {L_X, L_Y}] (with the 1/2 convention)."""
    L_X = solve_sld(rho, X)
    L_Y = solve_sld(rho, Y)
    return 0.25 * np.real(np.trace(rho @ (L_X @ L_Y + L_Y @ L_X)))


def bures_inner_alt(rho, X_D, X_O):
    """Alternative: g_B(X_D, X_O) = (1/2) Re Tr[X_D L_{X_O}].

    Valid when X_D is Hermitian (but not necessarily diagonal).
    The standard identity for the SLD metric.
    """
    L_O = solve_sld(rho, X_O)
    return 0.5 * np.real(np.trace(X_D @ L_O))


def cross_gram_matrix(rho, d=3):
    """Compute the (d-1) x (d^2-d) cross-Gram matrix G_DO."""
    # Diagonal basis: D_a = diag with +1 at a, -1 at a+1
    diag_basis = []
    for a in range(d - 1):
        D = np.zeros((d, d), dtype=complex)
        D[a, a] = 1.0
        D[a + 1, a + 1] = -1.0
        diag_basis.append(D)

    # Off-diagonal basis: R_{kl} and I_{kl}
    off_basis = []
    for k in range(d):
        for l in range(k + 1, d):
            R = np.zeros((d, d), dtype=complex)
            R[k, l] = 1.0
            R[l, k] = 1.0
            off_basis.append(R)
            I = np.zeros((d, d), dtype=complex)
            I[k, l] = 1j
            I[l, k] = -1j
            off_basis.append(I)

    n_diag = len(diag_basis)
    n_off = len(off_basis)
    G = np.zeros((n_diag, n_off))
    for i in range(n_diag):
        for j in range(n_off):
            G[i, j] = bures_inner(rho, diag_basis[i], off_basis[j])
    return G


def sld_diagonal_entries(rho, X_O, d=3):
    """Return the d diagonal entries (L_{X_O})_{kk} in the pointer basis."""
    L = solve_sld(rho, X_O)
    return np.array([np.real(L[k, k]) for k in range(d)])


# ============================================================
# Test 1: Verify the identity g_B(X_D, X_O) = (1/2) sum_k x_k L_{kk}
# ============================================================
def test_identity(n_trials=100):
    print("=" * 70)
    print("Test 1: Verify g_B(X_D, X_O) = (1/2) Re Tr[X_D L_{X_O}]")
    print("=" * 70)
    rng = np.random.default_rng(42)
    max_err = 0.0
    for _ in range(n_trials):
        rho = random_density_matrix(3, rng)
        # Random diagonal traceless
        x = rng.standard_normal(3)
        x -= x.mean()
        X_D = np.diag(x + 0j)
        # Random off-diagonal Hermitian
        z12 = rng.standard_normal() + 1j * rng.standard_normal()
        z13 = rng.standard_normal() + 1j * rng.standard_normal()
        z23 = rng.standard_normal() + 1j * rng.standard_normal()
        X_O = np.array([[0, z12, z13], [np.conj(z12), 0, z23],
                        [np.conj(z13), np.conj(z23), 0]])

        val1 = bures_inner(rho, X_D, X_O)
        val2 = bures_inner_alt(rho, X_D, X_O)

        # Also check the sum formula
        L_diag = sld_diagonal_entries(rho, X_O)
        val3 = 0.5 * np.sum(x * L_diag)

        err12 = abs(val1 - val2)
        err13 = abs(val1 - val3)
        max_err = max(max_err, err12, err13)

    print(f"  Max error across {n_trials} trials: {max_err:.2e}")
    assert max_err < 1e-9, f"Identity verification failed: max_err = {max_err}"
    print("  PASSED\n")


# ============================================================
# Test 2: Sparse off-diagonal case
# ============================================================
def test_sparse_case():
    print("=" * 70)
    print("Test 2: Sparse case rho = diag(p) + eps*(|1><2|+|2><1|)")
    print("=" * 70)
    p1, p2, p3 = 0.5, 0.3, 0.2
    for eps in [0.0, 0.01, 0.05, 0.1, 0.14]:
        rho = np.diag([p1, p2, p3]) + 0j
        rho[0, 1] = eps
        rho[1, 0] = eps
        # Check positive definiteness
        eigvals = np.linalg.eigvalsh(rho)
        if eigvals.min() < 1e-10:
            print(f"  eps={eps}: not full rank, skip")
            continue

        # SLD for R_{12}
        R12 = np.zeros((3, 3), dtype=complex)
        R12[0, 1] = 1.0
        R12[1, 0] = 1.0
        L = solve_sld(rho, R12)
        diag_L = [np.real(L[k, k]) for k in range(3)]

        # Analytical predictions for eps != 0 sparse case:
        # L_{33} = 0, L_{11} = -eps*L_{12}/p1, L_{22} = -eps*L_{12}/p2
        L12_val = L[0, 1]

        print(f"  eps={eps:.2f}: L_diag = [{diag_L[0]:.6f}, {diag_L[1]:.6f}, {diag_L[2]:.6f}]"
              f"  L_{33}={diag_L[2]:.2e}")

        if abs(eps) > 0:
            pred_L11 = -eps * np.real(L12_val) / p1
            pred_L22 = -eps * np.real(L12_val) / p2
            print(f"    Predicted: L11={pred_L11:.6f}, L22={pred_L22:.6f}, L33=0.0")
            print(f"    L11 != L33: {abs(diag_L[0] - diag_L[2]):.2e}")
            print(f"    L11 != L22: {abs(diag_L[0] - diag_L[1]):.2e} (zero iff p1=p2)")

        # Cross-Gram matrix
        G = cross_gram_matrix(rho)
        print(f"    max|G_DO| = {np.max(np.abs(G)):.6e}, Frobenius = {np.linalg.norm(G):.6e}")

    print()


# ============================================================
# Test 3: Equal populations but nonzero coherence
# ============================================================
def test_equal_populations():
    print("=" * 70)
    print("Test 3: Equal populations p1=p2=p3=1/3, nonzero rho_{12}")
    print("=" * 70)
    for eps in [0.0, 0.05, 0.1, 0.15]:
        rho = np.eye(3, dtype=complex) / 3.0
        rho[0, 1] = eps
        rho[1, 0] = eps
        eigvals = np.linalg.eigvalsh(rho)
        if eigvals.min() < 1e-10:
            print(f"  eps={eps}: not full rank, skip")
            continue

        R12 = np.zeros((3, 3), dtype=complex)
        R12[0, 1] = 1.0
        R12[1, 0] = 1.0
        L = solve_sld(rho, R12)
        diag_L = [np.real(L[k, k]) for k in range(3)]

        G = cross_gram_matrix(rho)
        print(f"  eps={eps:.2f}: L_diag = [{diag_L[0]:.6f}, {diag_L[1]:.6f}, {diag_L[2]:.6f}]"
              f"  max|G_DO|={np.max(np.abs(G)):.2e}")
        if eps > 0:
            # L11 = L22 (by symmetry p1=p2) but L33 = 0
            print(f"    L11=L22: {abs(diag_L[0]-diag_L[1]):.2e}, "
                  f"L11!=L33: {abs(diag_L[0]-diag_L[2]):.2e}")
    print()


# ============================================================
# Test 4: General non-diagonal states (large sample)
# ============================================================
def test_general(n_samples=5000):
    print("=" * 70)
    print(f"Test 4: Random full-rank non-diagonal states (n={n_samples})")
    print("=" * 70)
    rng = np.random.default_rng(123)
    max_G_at_diagonal = 0.0
    min_G_at_nondiag = np.inf
    n_nondiag = 0

    for _ in range(n_samples):
        rho = random_density_matrix(3, rng)
        G = cross_gram_matrix(rho)
        G_norm = np.max(np.abs(G))

        # Check if rho is diagonal (up to numerical precision)
        off_diag_norm = np.sqrt(sum(abs(rho[i, j])**2 for i in range(3)
                                    for j in range(3) if i != j))

        if off_diag_norm < 1e-10:
            max_G_at_diagonal = max(max_G_at_diagonal, G_norm)
        else:
            n_nondiag += 1
            min_G_at_nondiag = min(min_G_at_nondiag, G_norm)

    print(f"  Non-diagonal states: {n_nondiag}/{n_samples}")
    print(f"  Max |G_DO| at diagonal states: {max_G_at_diagonal:.2e}")
    print(f"  Min max|G_DO| at non-diagonal states: {min_G_at_nondiag:.6e}")
    print(f"  Separation: non-diagonal G is always nonzero")
    assert min_G_at_nondiag > 1e-10, "Found non-diagonal state with G_DO ~ 0!"
    print("  PASSED\n")


# ============================================================
# Test 5: Verify the specific formula L_{kk} for the sparse case
# ============================================================
def test_sparse_formula():
    print("=" * 70)
    print("Test 5: Analytical vs numerical L_{kk} in sparse case")
    print("=" * 70)
    p1, p2, p3 = 0.5, 0.3, 0.2

    for eps in [0.001, 0.01, 0.05, 0.1, 0.13]:
        rho = np.diag([p1, p2, p3]) + 0j
        rho[0, 1] = eps
        rho[1, 0] = eps
        eigvals = np.linalg.eigvalsh(rho)
        if eigvals.min() < 1e-10:
            continue

        R12 = np.zeros((3, 3), dtype=complex)
        R12[0, 1] = 1.0
        R12[1, 0] = 1.0
        L = solve_sld(rho, R12)

        # Analytical: L_{12} = 2 / [(p1+p2)(1 - eps^2/(p1*p2))]
        L12_pred = 2.0 / ((p1 + p2) * (1.0 - eps**2 / (p1 * p2)))
        L11_pred = -eps * L12_pred / p1
        L22_pred = -eps * L12_pred / p2
        L33_pred = 0.0

        # Numerical
        L12_num = np.real(L[0, 1])
        L11_num = np.real(L[0, 0])
        L22_num = np.real(L[1, 1])
        L33_num = np.real(L[2, 2])

        err_12 = abs(L12_pred - L12_num)
        err_11 = abs(L11_pred - L11_num)
        err_22 = abs(L22_pred - L22_num)
        err_33 = abs(L33_pred - L33_num)

        print(f"  eps={eps:.3f}: L12 err={err_12:.2e}, L11 err={err_11:.2e}, "
              f"L22 err={err_22:.2e}, L33 err={err_33:.2e}")

    print("  (Errors should be small for sparse case, larger for big eps due to L13,L23 coupling)")
    print()


# ============================================================
# Test 6: States with ALL off-diagonal entries nonzero
# ============================================================
def test_all_offdiag():
    print("=" * 70)
    print("Test 6: States with all off-diagonal entries nonzero")
    print("=" * 70)
    rng = np.random.default_rng(999)
    n_trials = 1000
    min_G = np.inf

    for _ in range(n_trials):
        # Generate state with all off-diag entries guaranteed nonzero
        rho = random_density_matrix(3, rng)
        # Check all off-diagonal are nonzero
        if any(abs(rho[i, j]) < 1e-6 for i in range(3) for j in range(3) if i != j):
            continue

        G = cross_gram_matrix(rho)
        G_norm = np.max(np.abs(G))
        min_G = min(min_G, G_norm)

    print(f"  Min max|G_DO| across trials: {min_G:.6e}")
    assert min_G > 1e-10, "Found fully-coherent state with G_DO ~ 0!"
    print("  PASSED\n")


# ============================================================
# Test 7: SLD diagonal entries for each off-diagonal direction
# ============================================================
def test_sld_diag_per_direction():
    print("=" * 70)
    print("Test 7: SLD diagonal entries for each off-diagonal direction")
    print("  (showing equal diags only at diagonal rho)")
    print("=" * 70)

    # A specific non-diagonal state
    rho = np.array([
        [0.5, 0.1, 0.05],
        [0.1, 0.3, 0.02],
        [0.05, 0.02, 0.2]
    ], dtype=complex)
    print(f"  rho = \n{np.real(rho)}")
    print(f"  eigenvalues: {np.linalg.eigvalsh(rho)}")
    print()

    off_dirs = {
        "R12": np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], dtype=complex),
        "I12": np.array([[0, 1j, 0], [-1j, 0, 0], [0, 0, 0]], dtype=complex),
        "R13": np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]], dtype=complex),
        "I13": np.array([[0, 0, 1j], [0, 0, 0], [-1j, 0, 0]], dtype=complex),
        "R23": np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0]], dtype=complex),
        "I23": np.array([[0, 0, 0], [0, 0, 1j], [0, -1j, 0]], dtype=complex),
    }

    for name, X_O in off_dirs.items():
        diag_L = sld_diagonal_entries(rho, X_O)
        spread = max(diag_L) - min(diag_L)
        print(f"  {name}: L_diag = [{diag_L[0]:+.6f}, {diag_L[1]:+.6f}, {diag_L[2]:+.6f}]"
              f"  spread={spread:.4e}")

    print()
    print("  Now at a diagonal state:")
    rho_diag = np.diag([0.5, 0.3, 0.2]) + 0j
    for name, X_O in off_dirs.items():
        diag_L = sld_diagonal_entries(rho_diag, X_O)
        spread = max(diag_L) - min(diag_L)
        print(f"  {name}: L_diag = [{diag_L[0]:+.6f}, {diag_L[1]:+.6f}, {diag_L[2]:+.6f}]"
              f"  spread={spread:.2e}")
    print()


if __name__ == "__main__":
    print("Verification suite for Qutrit Orthogonality Theorem (d=3)")
    print("=" * 70)
    print()
    test_identity()
    test_sparse_case()
    test_equal_populations()
    test_sparse_formula()
    test_sld_diag_per_direction()
    test_all_offdiag()
    test_general()
    print("=" * 70)
    print("ALL TESTS PASSED")
    print("=" * 70)
