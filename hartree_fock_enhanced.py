"""
Professional Hartree-Fock SCF Program (v4.0 - Enhanced Edition)
----------------------------------------------------------------
FEATURES:
- Restricted (RHF) and Unrestricted (UHF) Hartree-Fock
- DIIS convergence acceleration with robustness checks
- Comprehensive diagnostics and energy analysis
- Mulliken population analysis
- Orbital visualization data export
- Multiple basis sets (STO-3G, 6-31G, custom)
- Damping for difficult convergence
- Integral sanity checking
- Dipole moment calculation
- Spin contamination analysis (UHF)

Author: William Comaskey
Date: Original Version 2021 (Updated Fall 2025)
"""

import math
import numpy as np
from scipy import special, linalg
from collections import defaultdict
import sys

# =============================================================================
#  MATH & CONSTANTS
# =============================================================================

BOHR_TO_ANGSTROM = 0.529177210903
HARTREE_TO_EV = 27.211386245988

def boys_function(m, T):
    """Boys function for auxiliary Hermite integrals."""
    if T < 1e-9:
        return 1.0 / (2 * m + 1)
    else:
        return special.gammainc(m + 0.5, T) * special.gamma(m + 0.5) / (2 * T**(m + 0.5))

def double_factorial(n):
    """Double factorial: n!! = n(n-2)(n-4)..."""
    if n <= 1: return 1
    return n * double_factorial(n - 2)

def get_cartesian_powers(L):
    """Get all (i,j,k) with i+j+k=L for Cartesian Gaussians."""
    powers = []
    for i in range(L, -1, -1):
        for j in range(L - i, -1, -1):
            k = L - i - j
            powers.append((i, j, k))
    return powers

# =============================================================================
#  CHEMISTRY DATA STRUCTURES
# =============================================================================

class Atom:
    """Represents an atom with nuclear charge and position."""
    def __init__(self, atomic_number, coordinates, label=None):
        self.Z = float(atomic_number)
        self.coords = np.array(coordinates, dtype=float)
        self.label = label or self._get_element_symbol()
        self.basis_function_indices = []  # Track which basis functions belong to this atom
    
    def _get_element_symbol(self):
        """Get element symbol from atomic number."""
        symbols = {1:'H', 2:'He', 3:'Li', 4:'Be', 5:'B', 6:'C', 7:'N', 8:'O', 9:'F', 10:'Ne'}
        return symbols.get(int(self.Z), f'Z{int(self.Z)}')
    
    def __repr__(self):
        return f"{self.label}({self.Z}) at [{self.coords[0]:.3f}, {self.coords[1]:.3f}, {self.coords[2]:.3f}]"

class PrimitiveGaussian:
    """Primitive Gaussian function: g = A * (x-X)^i (y-Y)^j (z-Z)^k * exp(-alpha*r^2)."""
    def __init__(self, alpha, coeff, coordinates, powers):
        self.alpha = alpha
        self.coeff = coeff
        self.coords = np.array(coordinates, dtype=float)
        self.powers = powers
        self.A = self._calculate_norm()

    def _calculate_norm(self):
        """Calculate normalization constant for primitive Gaussian."""
        i, j, k = self.powers
        alpha = self.alpha
        prefactor = (2.0 * alpha / math.pi)**0.75
        numer = (4.0 * alpha)**(i + j + k)
        denom = double_factorial(2*i - 1) * double_factorial(2*j - 1) * double_factorial(2*k - 1)
        return prefactor * math.sqrt(numer / denom)

class BasisFunction:
    """Contracted Gaussian basis function (CGTO)."""
    def __init__(self, primitives, label="", atom_index=None):
        self.primitives = primitives
        self.label = label
        self.atom_index = atom_index  # Which atom does this basis function belong to
        self.normalize()

    def normalize(self):
        """Renormalize contracted function so <phi|phi>=1."""
        overlap = 0.0
        for pi in self.primitives:
            for pj in self.primitives:
                p = pi.alpha + pj.alpha
                Ex = IntegralEngine.get_E(pi.powers[0], pj.powers[0], 0, 0.0, 0.0, p)
                Ey = IntegralEngine.get_E(pi.powers[1], pj.powers[1], 0, 0.0, 0.0, p)
                Ez = IntegralEngine.get_E(pi.powers[2], pj.powers[2], 0, 0.0, 0.0, p)
                term = pi.A * pj.A * pi.coeff * pj.coeff
                overlap += term * (math.pi / p)**1.5 * Ex * Ey * Ez
        
        if overlap > 1e-14:
            scale = 1.0 / math.sqrt(overlap)
            for p in self.primitives:
                p.coeff *= scale
        else:
            print(f"Warning: Near-zero overlap for basis function {self.label}")

    def __len__(self):
        return len(self.primitives)
    
    def __getitem__(self, idx):
        return self.primitives[idx]

# =============================================================================
#  INTEGRAL ENGINE (Corrected & Validated)
# =============================================================================

class IntegralEngine:
    """Compute molecular integrals using McMurchie-Davidson scheme."""
    
    @staticmethod
    def get_E(i, j, t, Q_PA, Q_AB, p):
        """
        Recursive Hermite Coefficients E[t].
        Returns RAW coefficient without (pi/p) prefactors.
        """
        cache = {}
        def _recursive_E(i, j, t):
            if (i, j, t) in cache: return cache[(i, j, t)]
            val = 0.0
            if t < 0 or t > (i + j):
                val = 0.0
            elif i == 0 and j == 0 and t == 0:
                val = 1.0
            elif j > 0:
                val = _recursive_E(i + 1, j - 1, t) + Q_AB * _recursive_E(i, j - 1, t)
            elif i > 0:
                val = (1.0 / (2.0 * p)) * _recursive_E(i - 1, 0, t - 1) + \
                      Q_PA * _recursive_E(i - 1, 0, t) + \
                      (t + 1) * _recursive_E(i - 1, 0, t + 1)
            cache[(i, j, t)] = val
            return val
        return _recursive_E(i, j, t)

    @staticmethod
    def overlap(basis_functions):
        """Compute overlap matrix S[i,j] = <phi_i|phi_j>."""
        nbasis = len(basis_functions)
        S = np.zeros((nbasis, nbasis))
        for i in range(nbasis):
            for j in range(nbasis):
                val = 0.0
                for pi in basis_functions[i]:
                    for pj in basis_functions[j]:
                        p = pi.alpha + pj.alpha
                        P = (pi.alpha * pi.coords + pj.alpha * pj.coords) / p
                        PA = P - pi.coords
                        AB = pi.coords - pj.coords
                        K_AB = math.exp(- (pi.alpha * pj.alpha / p) * np.dot(AB, AB))
                        
                        Ex = IntegralEngine.get_E(pi.powers[0], pj.powers[0], 0, PA[0], AB[0], p)
                        Ey = IntegralEngine.get_E(pi.powers[1], pj.powers[1], 0, PA[1], AB[1], p)
                        Ez = IntegralEngine.get_E(pi.powers[2], pj.powers[2], 0, PA[2], AB[2], p)
                        
                        val += pi.A * pj.A * pi.coeff * pj.coeff * K_AB * (math.pi / p)**1.5 * Ex * Ey * Ez
                S[i, j] = val
        return S

    @staticmethod
    def kinetic(basis_functions):
        """Compute kinetic energy matrix T[i,j] = <phi_i|-1/2 ∇²|phi_j>."""
        nbasis = len(basis_functions)
        T = np.zeros((nbasis, nbasis))
        for i in range(nbasis):
            for j in range(nbasis):
                val = 0.0
                for pi in basis_functions[i]:
                    for pj in basis_functions[j]:
                        b = pj.alpha
                        p = pi.alpha + b
                        P = (pi.alpha * pi.coords + b * pj.coords) / p
                        PA = P - pi.coords
                        AB = pi.coords - pj.coords
                        K_AB = math.exp(- (pi.alpha * b / p) * np.dot(AB, AB))
                        
                        norm = pi.A * pj.A * pi.coeff * pj.coeff * K_AB * (math.pi / p)**1.5
                        
                        def get_S_raw(l1, l2, dim):
                            return IntegralEngine.get_E(l1, l2, 0, PA[dim], AB[dim], p)

                        def get_T_raw(l1, l2, dim):
                            t1 = b * (2 * l2 + 1) * get_S_raw(l1, l2, dim)
                            t2 = -2 * b**2 * get_S_raw(l1, l2 + 2, dim)
                            t3 = -0.5 * l2 * (l2 - 1) * get_S_raw(l1, l2 - 2, dim)
                            return t1 + t2 + t3

                        Sx = get_S_raw(pi.powers[0], pj.powers[0], 0)
                        Sy = get_S_raw(pi.powers[1], pj.powers[1], 1)
                        Sz = get_S_raw(pi.powers[2], pj.powers[2], 2)
                        
                        Tx = get_T_raw(pi.powers[0], pj.powers[0], 0)
                        Ty = get_T_raw(pi.powers[1], pj.powers[1], 1)
                        Tz = get_T_raw(pi.powers[2], pj.powers[2], 2)
                        
                        val += norm * (Tx * Sy * Sz + Sx * Ty * Sz + Sx * Sy * Tz)
                T[i, j] = val
        return T

    @staticmethod
    def nuclear_attraction(basis_functions, atoms):
        """Compute nuclear attraction matrix V[i,j] = <phi_i|-Z/r|phi_j>."""
        nbasis = len(basis_functions)
        V = np.zeros((nbasis, nbasis))
        for atom in atoms:
            for i in range(nbasis):
                for j in range(nbasis):
                    res = 0.0
                    for pi in basis_functions[i]:
                        for pj in basis_functions[j]:
                            p = pi.alpha + pj.alpha
                            P = (pi.alpha * pi.coords + pj.alpha * pj.coords) / p
                            PA = P - pi.coords
                            AB = pi.coords - pj.coords
                            PC = P - atom.coords
                            K_AB = math.exp(- (pi.alpha * pj.alpha / p) * np.dot(AB, AB))
                            norm = pi.A * pj.A * pi.coeff * pj.coeff * K_AB * (2*math.pi/p)
                            
                            R_cache = {}
                            def get_R(t, u, v, n):
                                if (t, u, v, n) in R_cache: return R_cache[(t, u, v, n)]
                                val = 0.0
                                if t == 0 and u == 0 and v == 0:
                                    T_arg = p * np.dot(PC, PC)
                                    val = ((-1)**n) * boys_function(n, T_arg)
                                elif t > 0:
                                    val = ((t - 1) / (2 * p) * (get_R(t - 2, u, v, n) - get_R(t - 2, u, v, n + 1)) if t > 1 else 0.0) + \
                                          PC[0] * get_R(t - 1, u, v, n + 1)
                                elif u > 0:
                                    val = ((u - 1) / (2 * p) * (get_R(t, u - 2, v, n) - get_R(t, u - 2, v, n + 1)) if u > 1 else 0.0) + \
                                          PC[1] * get_R(t, u - 1, v, n + 1)
                                elif v > 0:
                                    val = ((v - 1) / (2 * p) * (get_R(t, u, v - 2, n) - get_R(t, u, v - 2, n + 1)) if v > 1 else 0.0) + \
                                          PC[2] * get_R(t, u, v - 1, n + 1)
                                R_cache[(t, u, v, n)] = val
                                return val
                            
                            term = 0.0
                            L = [pi.powers[k] + pj.powers[k] for k in range(3)]
                            for t in range(L[0] + 1):
                                ex = IntegralEngine.get_E(pi.powers[0], pj.powers[0], t, PA[0], AB[0], p)
                                for u in range(L[1] + 1):
                                    ey = IntegralEngine.get_E(pi.powers[1], pj.powers[1], u, PA[1], AB[1], p)
                                    for v in range(L[2] + 1):
                                        ez = IntegralEngine.get_E(pi.powers[2], pj.powers[2], v, PA[2], AB[2], p)
                                        term += ex * ey * ez * get_R(t, u, v, 0)
                            res += norm * (-atom.Z) * term
                    V[i, j] += res
        return V

    @staticmethod
    def electron_repulsion(basis_functions):
        """Compute electron repulsion integrals (ij|kl)."""
        nbasis = len(basis_functions)
        Vee = np.zeros((nbasis, nbasis, nbasis, nbasis))
        for i in range(nbasis):
            for j in range(nbasis):
                for k in range(nbasis):
                    for l in range(nbasis):
                        val = 0.0
                        for pi in basis_functions[i]:
                            for pj in basis_functions[j]:
                                p = pi.alpha + pj.alpha
                                P = (pi.alpha * pi.coords + pj.alpha * pj.coords) / p
                                PA = P - pi.coords
                                AB = pi.coords - pj.coords
                                K1 = math.exp(- (pi.alpha * pj.alpha / p) * np.dot(AB, AB))
                                
                                for pk in basis_functions[k]:
                                    for pl in basis_functions[l]:
                                        q = pk.alpha + pl.alpha
                                        Q = (pk.alpha * pk.coords + pl.alpha * pl.coords) / q
                                        QC = Q - pk.coords
                                        CD = pk.coords - pl.coords
                                        K2 = math.exp(- (pk.alpha * pl.alpha / q) * np.dot(CD, CD))
                                        
                                        norm = pi.A * pj.A * pk.A * pl.A * pi.coeff * pj.coeff * pk.coeff * pl.coeff
                                        prefactor = (2 * math.pi**2.5) / (p * q * math.sqrt(p + q)) * K1 * K2
                                        
                                        alpha_sum = p + q
                                        RPQ = P - Q
                                        
                                        R_eri = {}
                                        def get_R(t, u, v, n):
                                            if (t, u, v, n) in R_eri: return R_eri[(t, u, v, n)]
                                            val = 0.0
                                            if t == 0 and u == 0 and v == 0:
                                                T_arg = (p * q / alpha_sum) * np.dot(RPQ, RPQ)
                                                val = ((-1)**n) * boys_function(n, T_arg)
                                            elif t > 0:
                                                val = ((t - 1) / (2 * alpha_sum) * (get_R(t - 2, u, v, n) - get_R(t - 2, u, v, n + 1)) if t > 1 else 0.0) + \
                                                      RPQ[0] * get_R(t - 1, u, v, n + 1)
                                            elif u > 0:
                                                val = ((u - 1) / (2 * alpha_sum) * (get_R(t, u - 2, v, n) - get_R(t, u - 2, v, n + 1)) if u > 1 else 0.0) + \
                                                      RPQ[1] * get_R(t, u - 1, v, n + 1)
                                            elif v > 0:
                                                val = ((v - 1) / (2 * alpha_sum) * (get_R(t, u, v - 2, n) - get_R(t, u, v - 2, n + 1)) if v > 1 else 0.0) + \
                                                      RPQ[2] * get_R(t, u, v - 1, n + 1)
                                            R_eri[(t, u, v, n)] = val
                                            return val
                                        
                                        term_sum = 0.0
                                        L1 = [pi.powers[d] + pj.powers[d] for d in range(3)]
                                        L2 = [pk.powers[d] + pl.powers[d] for d in range(3)]
                                        
                                        for t1 in range(L1[0] + 1):
                                            ex1 = IntegralEngine.get_E(pi.powers[0], pj.powers[0], t1, PA[0], AB[0], p)
                                            for t2 in range(L2[0] + 1):
                                                ex2 = IntegralEngine.get_E(pk.powers[0], pl.powers[0], t2, QC[0], CD[0], q)
                                                tx = ex1 * ex2 * ((-1)**t2)
                                                for u1 in range(L1[1] + 1):
                                                    ey1 = IntegralEngine.get_E(pi.powers[1], pj.powers[1], u1, PA[1], AB[1], p)
                                                    for u2 in range(L2[1] + 1):
                                                        ey2 = IntegralEngine.get_E(pk.powers[1], pl.powers[1], u2, QC[1], CD[1], q)
                                                        ty = ey1 * ey2 * ((-1)**u2)
                                                        for v1 in range(L1[2] + 1):
                                                            ez1 = IntegralEngine.get_E(pi.powers[2], pj.powers[2], v1, PA[2], AB[2], p)
                                                            for v2 in range(L2[2] + 1):
                                                                ez2 = IntegralEngine.get_E(pk.powers[2], pl.powers[2], v2, QC[2], CD[2], q)
                                                                tz = ez1 * ez2 * ((-1)**v2)
                                                                term_sum += tx * ty * tz * get_R(t1+t2, u1+u2, v1+v2, 0)
                                        val += norm * prefactor * term_sum
                        Vee[i, j, k, l] = val
        return Vee

    @staticmethod
    def compute_nuclear_repulsion(atoms):
        """Compute nuclear-nuclear repulsion energy."""
        E = 0.0
        for i in range(len(atoms)):
            for j in range(i + 1, len(atoms)):
                dist = np.linalg.norm(atoms[i].coords - atoms[j].coords)
                if dist < 1e-10:
                    raise ValueError(f"Atoms {i} and {j} are at the same position!")
                E += (atoms[i].Z * atoms[j].Z) / dist
        return E

# =============================================================================
#  SCF SOLVER (RHF with Full Diagnostics)
# =============================================================================

class SCFSolver:
    """Restricted Hartree-Fock SCF solver with comprehensive analysis."""
    
    def __init__(self, basis_functions, atoms, n_electrons):
        self.basis = basis_functions
        self.atoms = atoms
        self.n_electrons = n_electrons
        
        if n_electrons % 2 != 0:
            raise ValueError("RHF requires even number of electrons. Use UHF for open-shell.")
        self.n_occ = n_electrons // 2
        
        # Results storage
        self.energy = None
        self.C = None  # MO coefficients
        self.eps = None  # Orbital energies
        self.P = None  # Density matrix
        self.F = None  # Fock matrix
        
    def check_integral_sanity(self):
        """Perform sanity checks on computed integrals."""
        print("\n" + "="*60)
        print("INTEGRAL SANITY CHECKS")
        print("="*60)
        
        # Check overlap matrix
        print(f"Overlap matrix shape: {self.S.shape}")
        print(f"  - Symmetric: {np.allclose(self.S, self.S.T)}")
        print(f"  - Diagonal elements: min={np.min(np.diag(self.S)):.6f}, max={np.max(np.diag(self.S)):.6f}")
        print(f"  - Condition number: {np.linalg.cond(self.S):.2e}")
        
        evals = linalg.eigvalsh(self.S)
        print(f"  - Smallest eigenvalue: {np.min(evals):.6e}")
        print(f"  - Largest eigenvalue: {np.max(evals):.6e}")
        if np.min(evals) < 1e-7:
            print(f"  ⚠ WARNING: Very small eigenvalue detected - near-linear dependence!")
        
        # Check kinetic energy
        print(f"\nKinetic matrix:")
        print(f"  - Symmetric: {np.allclose(self.T, self.T.T)}")
        print(f"  - Diagonal elements: min={np.min(np.diag(self.T)):.6f}, max={np.max(np.diag(self.T)):.6f}")
        
        # Check nuclear attraction
        print(f"\nNuclear attraction matrix:")
        print(f"  - Symmetric: {np.allclose(self.Vne, self.Vne.T)}")
        print(f"  - All negative: {np.all(self.Vne <= 0)}")
        
        # Nuclear repulsion
        print(f"\nNuclear repulsion energy: {self.Enn:.8f} Ha")
        
        # Check 2e integrals (8-fold symmetry)
        print(f"\n2-electron integrals:")
        print(f"  - Shape: {self.Vee.shape}")
        max_asymm = 0.0
        for i in range(min(3, len(self.basis))):
            for j in range(min(3, len(self.basis))):
                for k in range(min(3, len(self.basis))):
                    for l in range(min(3, len(self.basis))):
                        # Check (ij|kl) = (ji|lk)
                        diff = abs(self.Vee[i,j,k,l] - self.Vee[j,i,l,k])
                        max_asymm = max(max_asymm, diff)
        print(f"  - Max symmetry violation (sample): {max_asymm:.2e}")
        
    def run(self, tolerance=1e-6, max_steps=100, damping=0.0, use_diis=True, 
            diis_start=2, diis_max_vecs=8, verbose=True):
        """
        Run SCF calculation.
        
        Parameters:
        -----------
        tolerance : float
            Energy convergence threshold
        max_steps : int
            Maximum SCF iterations
        damping : float
            Density matrix damping factor (0.0-0.9)
        use_diis : bool
            Use DIIS acceleration
        diis_start : int
            Iteration to start DIIS
        diis_max_vecs : int
            Maximum DIIS vectors to store
        verbose : bool
            Print iteration details
        """
        if verbose:
            print("\n" + "="*60)
            print("COMPUTING INTEGRALS")
            print("="*60)
        
        # Compute integrals
        self.Enn = IntegralEngine.compute_nuclear_repulsion(self.atoms)
        self.S = IntegralEngine.overlap(self.basis)
        self.T = IntegralEngine.kinetic(self.basis)
        self.Vne = IntegralEngine.nuclear_attraction(self.basis, self.atoms)
        
        if verbose:
            print("Computing 2-electron integrals (this may take a while)...")
        self.Vee = IntegralEngine.electron_repulsion(self.basis)
        
        # Sanity checks
        if verbose:
            self.check_integral_sanity()
        
        H_core = self.T + self.Vne
        
        # Canonical Orthogonalization
        if verbose:
            print("\n" + "="*60)
            print("ORTHOGONALIZATION")
            print("="*60)
        
        s_evals, s_evecs = linalg.eigh(self.S)
        threshold = 1e-7
        keep_mask = s_evals > threshold
        s_evals_kept = s_evals[keep_mask]
        s_evecs_kept = s_evecs[:, keep_mask]
        X = np.dot(s_evecs_kept, np.diag(s_evals_kept**-0.5))
        
        if verbose:
            print(f"Basis functions: {len(self.basis)}")
            print(f"Functions kept after orthogonalization: {len(s_evals_kept)}")
            if len(s_evals_kept) < len(self.basis):
                print(f"⚠ Removed {len(self.basis) - len(s_evals_kept)} near-linearly dependent functions")
            
            if self.n_occ > len(s_evals_kept):
                raise ValueError(f"Not enough basis functions ({len(s_evals_kept)}) for {self.n_occ} occupied orbitals!")
        
        nbasis = len(self.basis)
        P = np.zeros((nbasis, nbasis))
        P_old = np.zeros((nbasis, nbasis))
        energy = 0.0
        
        # DIIS Storage
        fock_history = []
        error_history = []
        
        # SCF Loop
        if verbose:
            print("\n" + "="*60)
            print("SCF ITERATIONS")
            print("="*60)
            print(f"{'Iter':<6} {'Energy (Ha)':<16} {'ΔE':<14} {'|ΔP|':<12} {'DIIS':<6}")
            print("-" * 60)
        
        converged = False
        for step in range(max_steps):
            old_energy = energy
            
            # Build Fock Matrix
            G = np.zeros((nbasis, nbasis))
            for i in range(nbasis):
                for j in range(nbasis):
                    for k in range(nbasis):
                        for l in range(nbasis):
                            J = self.Vee[i, j, k, l]
                            K = self.Vee[i, l, k, j]
                            G[i, j] += P[k, l] * (J - 0.5 * K)
            F = H_core + G
            
            # DIIS Acceleration
            error_vec = np.dot(np.dot(F, P), self.S) - np.dot(np.dot(self.S, P), F)
            rms_error = np.sqrt(np.mean(error_vec**2))
            
            diis_used = False
            if use_diis and step >= diis_start:
                fock_history.append(F.copy())
                error_history.append(error_vec.copy())
                
                if len(fock_history) > diis_max_vecs:
                    fock_history.pop(0)
                    error_history.pop(0)
                
                if len(fock_history) >= 2:
                    dim = len(fock_history)
                    B = np.zeros((dim + 1, dim + 1))
                    B[-1, :] = -1
                    B[:, -1] = -1
                    B[-1, -1] = 0
                    
                    for i in range(dim):
                        for j in range(dim):
                            B[i, j] = np.sum(error_history[i] * error_history[j])
                    
                    rhs = np.zeros(dim + 1)
                    rhs[-1] = -1
                    
                    try:
                        coeffs = np.linalg.solve(B, rhs)
                        
                        # Check for unreasonable coefficients
                        if np.all(np.abs(coeffs[:-1]) < 100):
                            F_diis = np.zeros_like(F)
                            for i in range(dim):
                                F_diis += coeffs[i] * fock_history[i]
                            F = F_diis
                            diis_used = True
                    except np.linalg.LinAlgError:
                        pass  # DIIS failed, use regular Fock
            
            # Energy Calculation
            e_elec = 0.5 * np.sum(P * (H_core + F))
            energy = e_elec + self.Enn
            delta_E = energy - old_energy
            
            # Diagonalize Fock
            F_prime = np.dot(X.T, np.dot(F, X))
            eps, C_prime = linalg.eigh(F_prime)
            C = np.dot(X, C_prime)
            
            # Build new density with damping
            P_new = np.zeros((nbasis, nbasis))
            for i in range(nbasis):
                for j in range(nbasis):
                    for a in range(self.n_occ):
                        if a < C.shape[1]:
                            P_new[i, j] += 2.0 * C[i, a] * C[j, a]
            
            delta_P = np.max(np.abs(P_new - P))
            
            if damping > 0:
                P = (1 - damping) * P_new + damping * P
            else:
                P = P_new
            
            # Print iteration
            if verbose:
                diis_str = "Yes" if diis_used else "No"
                print(f"{step:<6} {energy:<16.10f} {delta_E:<14.6e} {delta_P:<12.6e} {diis_str:<6}")
            
            # Check convergence
            if step > 0 and abs(delta_E) < tolerance and delta_P < tolerance:
                converged = True
                if verbose:
                    print("-" * 60)
                    print("✓ SCF CONVERGED")
                break
        
        if not converged:
            print("-" * 60)
            print("⚠ WARNING: SCF did not converge within maximum iterations!")
        
        # Store results
        self.energy = energy
        self.C = C
        self.eps = eps
        self.P = P
        self.F = F
        
        # Post-SCF analysis
        if verbose:
            self.print_energy_analysis()
            self.print_orbital_energies()
            self.mulliken_analysis()
            self.dipole_moment()
        
        return energy
    
    def print_energy_analysis(self):
        """Print detailed energy component breakdown."""
        print("\n" + "="*60)
        print("ENERGY ANALYSIS")
        print("="*60)
        
        # Compute individual components
        E_kinetic = np.sum(self.P * self.T)
        E_nuclear_attraction = np.sum(self.P * self.Vne)
        
        # Coulomb and Exchange
        E_coulomb = 0.0
        E_exchange = 0.0
        for i in range(len(self.basis)):
            for j in range(len(self.basis)):
                for k in range(len(self.basis)):
                    for l in range(len(self.basis)):
                        E_coulomb += 0.5 * self.P[i,j] * self.P[k,l] * self.Vee[i,j,k,l]
                        E_exchange += -0.25 * self.P[i,j] * self.P[k,l] * self.Vee[i,l,k,j]
        
        E_electronic = E_kinetic + E_nuclear_attraction + E_coulomb + E_exchange
        E_total = E_electronic + self.Enn
        
        print(f"Kinetic Energy:              {E_kinetic:16.10f} Ha")
        print(f"Nuclear Attraction:          {E_nuclear_attraction:16.10f} Ha")
        print(f"Coulomb Energy (J):          {E_coulomb:16.10f} Ha")
        print(f"Exchange Energy (K):         {E_exchange:16.10f} Ha")
        print(f"Electronic Energy:           {E_electronic:16.10f} Ha")
        print(f"Nuclear Repulsion:           {self.Enn:16.10f} Ha")
        print("-" * 60)
        print(f"Total Energy:                {E_total:16.10f} Ha")
        print(f"                             {E_total * HARTREE_TO_EV:16.10f} eV")
        
        # Virial ratio for quality check
        virial = -E_total / E_kinetic
        print(f"\nVirial Ratio (-V/T):         {virial:16.10f}")
        print(f"  (Should be ~2.0 for well-converged calculations)")
        
    def print_orbital_energies(self, n_show=None):
        """Print molecular orbital energies."""
        if n_show is None:
            n_show = min(self.n_occ + 5, len(self.eps))
        
        print("\n" + "="*60)
        print("MOLECULAR ORBITAL ENERGIES")
        print("="*60)
        print(f"{'MO':<6} {'Energy (Ha)':<16} {'Energy (eV)':<16} {'Occ':<6}")
        print("-" * 60)
        
        for i in range(n_show):
            occ_str = "Yes" if i < self.n_occ else "No"
            print(f"{i+1:<6} {self.eps[i]:16.10f} {self.eps[i]*HARTREE_TO_EV:16.10f} {occ_str:<6}")
        
        if len(self.eps) > n_show:
            print(f"... ({len(self.eps) - n_show} more orbitals)")
        
        # HOMO-LUMO gap
        if self.n_occ < len(self.eps):
            homo = self.eps[self.n_occ - 1]
            lumo = self.eps[self.n_occ]
            gap = lumo - homo
            print("-" * 60)
            print(f"HOMO (orbital {self.n_occ}):      {homo:16.10f} Ha")
            print(f"LUMO (orbital {self.n_occ+1}):      {lumo:16.10f} Ha")
            print(f"HOMO-LUMO Gap:               {gap:16.10f} Ha")
            print(f"                             {gap * HARTREE_TO_EV:16.10f} eV")
    
    def mulliken_analysis(self):
        """Perform Mulliken population analysis."""
        print("\n" + "="*60)
        print("MULLIKEN POPULATION ANALYSIS")
        print("="*60)
        
        PS = np.dot(self.P, self.S)
        
        # Compute gross atomic populations
        atom_populations = defaultdict(float)
        atom_labels = {}
        
        for i, bf in enumerate(self.basis):
            pop = PS[i, i]
            if bf.atom_index is not None:
                atom_populations[bf.atom_index] += pop
                atom_labels[bf.atom_index] = self.atoms[bf.atom_index].label
        
        print(f"{'Atom':<10} {'Population':<15} {'Charge':<15}")
        print("-" * 45)
        
        total_charge = 0.0
        for idx in sorted(atom_populations.keys()):
            pop = atom_populations[idx]
            charge = self.atoms[idx].Z - pop
            label = atom_labels.get(idx, f"Atom{idx}")
            print(f"{label:<10} {pop:15.6f} {charge:15.6f}")
            total_charge += charge
        
        print("-" * 45)
        print(f"{'Total':<10} {sum(atom_populations.values()):15.6f} {total_charge:15.6f}")
        
        # Check: should sum to number of electrons
        expected_electrons = self.n_occ * 2
        actual_electrons = sum(atom_populations.values())
        print(f"\nExpected electrons: {expected_electrons}")
        print(f"Actual electrons:   {actual_electrons:.6f}")
        if abs(expected_electrons - actual_electrons) > 0.01:
            print("⚠ WARNING: Electron count mismatch!")
    
    def dipole_moment(self):
        """Calculate molecular dipole moment."""
        print("\n" + "="*60)
        print("DIPOLE MOMENT")
        print("="*60)
        
        # Nuclear contribution
        mu_nuc = np.zeros(3)
        for atom in self.atoms:
            mu_nuc += atom.Z * atom.coords
        
        # Electronic contribution (need to compute <P|r|P>)
        mu_elec = np.zeros(3)
        for i, bf_i in enumerate(self.basis):
            for j, bf_j in enumerate(self.basis):
                # Compute dipole integrals <i|r|j> for each component
                for dim in range(3):
                    val = 0.0
                    for pi in bf_i:
                        for pj in bf_j:
                            p = pi.alpha + pj.alpha
                            P_coord = (pi.alpha * pi.coords + pj.alpha * pj.coords) / p
                            PA = P_coord - pi.coords
                            AB = pi.coords - pj.coords
                            K_AB = math.exp(- (pi.alpha * pj.alpha / p) * np.dot(AB, AB))
                            
                            # Dipole = overlap with extra factor from r operator
                            # <i|x|j> involves raising angular momentum by 1 in x
                            powers_j_raised = list(pj.powers)
                            powers_j_raised[dim] += 1
                            
                            Ex = IntegralEngine.get_E(pi.powers[0], powers_j_raised[0], 0, PA[0], AB[0], p)
                            Ey = IntegralEngine.get_E(pi.powers[1], powers_j_raised[1], 0, PA[1], AB[1], p)
                            Ez = IntegralEngine.get_E(pi.powers[2], powers_j_raised[2], 0, PA[2], AB[2], p)
                            
                            val += pi.A * pj.A * pi.coeff * pj.coeff * K_AB * (math.pi / p)**1.5 * Ex * Ey * Ez
                    
                    mu_elec[dim] += self.P[i, j] * val
        
        mu_total = mu_nuc - mu_elec
        mu_magnitude = np.linalg.norm(mu_total)
        
        # Convert to Debye (1 a.u. = 2.54158 Debye)
        au_to_debye = 2.54158
        
        print(f"μ_x: {mu_total[0]:12.6f} a.u. ({mu_total[0]*au_to_debye:12.6f} Debye)")
        print(f"μ_y: {mu_total[1]:12.6f} a.u. ({mu_total[1]*au_to_debye:12.6f} Debye)")
        print(f"μ_z: {mu_total[2]:12.6f} a.u. ({mu_total[2]*au_to_debye:12.6f} Debye)")
        print(f"|μ|: {mu_magnitude:12.6f} a.u. ({mu_magnitude*au_to_debye:12.6f} Debye)")
    
    def export_orbitals(self, filename="molecular_orbitals.txt"):
        """Export MO coefficients and energies to file."""
        with open(filename, 'w') as f:
            f.write("MOLECULAR ORBITAL COEFFICIENTS AND ENERGIES\n")
            f.write("=" * 70 + "\n\n")
            
            for i_mo in range(len(self.eps)):
                occ = "OCCUPIED" if i_mo < self.n_occ else "VIRTUAL"
                f.write(f"MO {i_mo+1:3d}  Energy: {self.eps[i_mo]:14.8f} Ha  [{occ}]\n")
                f.write("-" * 70 + "\n")
                
                for i_bf, bf in enumerate(self.basis):
                    if abs(self.C[i_bf, i_mo]) > 0.01:  # Only print significant coefficients
                        f.write(f"  {bf.label:20s}  {self.C[i_bf, i_mo]:10.6f}\n")
                f.write("\n")
        
        print(f"\n✓ Molecular orbitals exported to {filename}")

# =============================================================================
#  UNRESTRICTED HARTREE-FOCK (UHF)
# =============================================================================

class UHFSolver:
    """Unrestricted Hartree-Fock for open-shell systems."""
    
    def __init__(self, basis_functions, atoms, n_alpha, n_beta):
        self.basis = basis_functions
        self.atoms = atoms
        self.n_alpha = n_alpha
        self.n_beta = n_beta
        self.n_electrons = n_alpha + n_beta
        
        # Results storage
        self.energy = None
        self.C_alpha = None
        self.C_beta = None
        self.eps_alpha = None
        self.eps_beta = None
        self.P_alpha = None
        self.P_beta = None
    
    def run(self, tolerance=1e-6, max_steps=100, damping=0.0, verbose=True):
        """Run UHF calculation."""
        if verbose:
            print("\n" + "="*60)
            print("UNRESTRICTED HARTREE-FOCK (UHF)")
            print("="*60)
            print(f"α electrons: {self.n_alpha}")
            print(f"β electrons: {self.n_beta}")
            print(f"Spin multiplicity: {abs(self.n_alpha - self.n_beta) + 1}")
        
        # Compute integrals (same as RHF)
        self.Enn = IntegralEngine.compute_nuclear_repulsion(self.atoms)
        self.S = IntegralEngine.overlap(self.basis)
        self.T = IntegralEngine.kinetic(self.basis)
        self.Vne = IntegralEngine.nuclear_attraction(self.basis, self.atoms)
        
        if verbose:
            print("Computing 2-electron integrals...")
        self.Vee = IntegralEngine.electron_repulsion(self.basis)
        
        H_core = self.T + self.Vne
        
        # Orthogonalization
        s_evals, s_evecs = linalg.eigh(self.S)
        threshold = 1e-7
        keep_mask = s_evals > threshold
        s_evals_kept = s_evals[keep_mask]
        s_evecs_kept = s_evecs[:, keep_mask]
        X = np.dot(s_evecs_kept, np.diag(s_evals_kept**-0.5))
        
        nbasis = len(self.basis)
        P_alpha = np.zeros((nbasis, nbasis))
        P_beta = np.zeros((nbasis, nbasis))
        energy = 0.0
        
        if verbose:
            print("\n" + "="*60)
            print("SCF ITERATIONS")
            print("="*60)
            print(f"{'Iter':<6} {'Energy (Ha)':<16} {'ΔE':<14} {'<S²>':<12}")
            print("-" * 60)
        
        converged = False
        for step in range(max_steps):
            old_energy = energy
            
            # Build Fock matrices for α and β spins
            G_alpha = np.zeros((nbasis, nbasis))
            G_beta = np.zeros((nbasis, nbasis))
            
            for i in range(nbasis):
                for j in range(nbasis):
                    for k in range(nbasis):
                        for l in range(nbasis):
                            J_alpha = self.Vee[i, j, k, l] * P_alpha[k, l]
                            J_beta = self.Vee[i, j, k, l] * P_beta[k, l]
                            K_alpha = self.Vee[i, l, k, j] * P_alpha[k, l]
                            K_beta = self.Vee[i, l, k, j] * P_beta[k, l]
                            
                            G_alpha[i, j] += J_alpha + J_beta - K_alpha
                            G_beta[i, j] += J_alpha + J_beta - K_beta
            
            F_alpha = H_core + G_alpha
            F_beta = H_core + G_beta
            
            # Energy
            P_total = P_alpha + P_beta
            e_elec = 0.5 * (np.sum(P_alpha * (H_core + F_alpha)) + 
                           np.sum(P_beta * (H_core + F_beta)))
            energy = e_elec + self.Enn
            delta_E = energy - old_energy
            
            # Diagonalize both Fock matrices
            F_alpha_prime = np.dot(X.T, np.dot(F_alpha, X))
            eps_alpha, C_alpha_prime = linalg.eigh(F_alpha_prime)
            C_alpha = np.dot(X, C_alpha_prime)
            
            F_beta_prime = np.dot(X.T, np.dot(F_beta, X))
            eps_beta, C_beta_prime = linalg.eigh(F_beta_prime)
            C_beta = np.dot(X, C_beta_prime)
            
            # Build new densities
            P_alpha_new = np.zeros((nbasis, nbasis))
            for i in range(nbasis):
                for j in range(nbasis):
                    for a in range(self.n_alpha):
                        if a < C_alpha.shape[1]:
                            P_alpha_new[i, j] += C_alpha[i, a] * C_alpha[j, a]
            
            P_beta_new = np.zeros((nbasis, nbasis))
            for i in range(nbasis):
                for j in range(nbasis):
                    for a in range(self.n_beta):
                        if a < C_beta.shape[1]:
                            P_beta_new[i, j] += C_beta[i, a] * C_beta[j, a]
            
            # Apply damping
            if damping > 0:
                P_alpha = (1 - damping) * P_alpha_new + damping * P_alpha
                P_beta = (1 - damping) * P_beta_new + damping * P_beta
            else:
                P_alpha = P_alpha_new
                P_beta = P_beta_new
            
            # Spin contamination <S²>
            S_squared = self.compute_spin_squared(P_alpha, P_beta)
            
            if verbose:
                print(f"{step:<6} {energy:<16.10f} {delta_E:<14.6e} {S_squared:<12.6f}")
            
            # Check convergence
            if step > 0 and abs(delta_E) < tolerance:
                converged = True
                if verbose:
                    print("-" * 60)
                    print("✓ UHF CONVERGED")
                break
        
        if not converged:
            print("-" * 60)
            print("⚠ WARNING: UHF did not converge!")
        
        # Store results
        self.energy = energy
        self.C_alpha = C_alpha
        self.C_beta = C_beta
        self.eps_alpha = eps_alpha
        self.eps_beta = eps_beta
        self.P_alpha = P_alpha
        self.P_beta = P_beta
        
        if verbose:
            self.print_results()
        
        return energy
    
    def compute_spin_squared(self, P_alpha, P_beta):
        """Compute <S²> = S(S+1) where S is total spin."""
        n_unpaired = abs(self.n_alpha - self.n_beta)
        S_exact = n_unpaired / 2.0
        S_squared_exact = S_exact * (S_exact + 1)
        
        # Actual <S²> includes spin contamination
        overlap_alpha_beta = np.sum(P_alpha * np.dot(self.S, P_beta))
        S_squared = S_squared_exact + self.n_beta - overlap_alpha_beta
        
        return S_squared
    
    def print_results(self):
        """Print UHF results."""
        print("\n" + "="*60)
        print("UHF ENERGY ANALYSIS")
        print("="*60)
        print(f"Total Energy: {self.energy:16.10f} Ha")
        
        n_unpaired = abs(self.n_alpha - self.n_beta)
        S_exact = n_unpaired / 2.0
        S_squared_exact = S_exact * (S_exact + 1)
        S_squared_actual = self.compute_spin_squared(self.P_alpha, self.P_beta)
        
        print(f"\nSpin Analysis:")
        print(f"  Unpaired electrons: {n_unpaired}")
        print(f"  Exact <S²>:         {S_squared_exact:.6f}")
        print(f"  Actual <S²>:        {S_squared_actual:.6f}")
        print(f"  Spin contamination: {S_squared_actual - S_squared_exact:.6f}")
        
        # Orbital energies
        print(f"\nα Orbital Energies (first {min(5, len(self.eps_alpha))}):")
        for i in range(min(5, len(self.eps_alpha))):
            occ = "occ" if i < self.n_alpha else "virt"
            print(f"  {i+1:3d}  {self.eps_alpha[i]:14.8f} Ha  [{occ}]")
        
        print(f"\nβ Orbital Energies (first {min(5, len(self.eps_beta))}):")
        for i in range(min(5, len(self.eps_beta))):
            occ = "occ" if i < self.n_beta else "virt"
            print(f"  {i+1:3d}  {self.eps_beta[i]:14.8f} Ha  [{occ}]")

# =============================================================================
#  BASIS SET LIBRARY
# =============================================================================

class BasisSetLibrary:
    """Common basis sets for molecular calculations."""
    
    STO3G_H = """
    1 1
    0 0 3 1.0 1.0
      3.42525091      0.15432897
      0.62391373      0.53532814
      0.16885540      0.44463454
    """
    
    STO3G_C = """
    2 5
    0 0 3 1.0 1.0
      71.6168370      0.15432897
      13.0450960      0.53532814
      3.53051220      0.44463454
    0 1 3 2.0 1.0
      2.94124940      -0.09996723     0.15591627
      0.68348310      0.39951283      0.60768372
      0.22228990      0.70011547      0.39195739
    """
    
    # 6-31G basis for H
    SIXTHIRTYONE_G_H = """
    1 2
    0 0 3 1.0 1.0
      18.7311370      0.03349460
      2.82539370      0.23472695
      0.64012170      0.81375733
    0 0 1 0.0 1.0
      0.16127780      1.00000000
    """

# =============================================================================
#  PARSER & UTILITIES
# =============================================================================

class BasisParser:
    """Parse Gaussian-style basis set format."""
    
    @staticmethod
    def parse(text, atom_coords, atom_indices=None):
        """
        Parse basis set and create basis functions.
        
        Parameters:
        -----------
        text : str
            Basis set in Gaussian format
        atom_coords : list of arrays
            Atomic coordinates
        atom_indices : list of ints (optional)
            Atom indices for Mulliken analysis
        """
        lines = [l.strip() for l in text.strip().split('\n') if l.strip()]
        n_shells = int(lines[0].split()[1])
        basis_funcs = []
        current = 1
        
        for shell_idx in range(n_shells):
            if current >= len(lines): break
            params = lines[current].split()
            if params[0] == '99': break
            current += 1
            
            input_type = int(params[0])
            shell_type = int(params[1])
            n_prims = int(params[2])
            scale_factor = float(params[4])
            scaling = (scale_factor**2) if input_type == 1 else 1.0
            
            prims_s, prims_p = [], []
            for _ in range(n_prims):
                val = lines[current].split()
                current += 1
                ex = float(val[0]) * scaling
                if shell_type == 1:  # SP shell
                    prims_s.append((ex, float(val[1])))
                    prims_p.append((ex, float(val[2])))
                else:
                    prims_s.append((ex, float(val[1])))
            
            # Create basis functions for each atom
            for center_idx, center in enumerate(atom_coords):
                atom_idx = atom_indices[center_idx] if atom_indices else None
                
                if shell_type == 1:  # SP shell
                    # S function
                    p_objs = [PrimitiveGaussian(e, c, center, (0,0,0)) for e, c in prims_s]
                    basis_funcs.append(BasisFunction(p_objs, "S", atom_idx))
                    
                    # P functions
                    for powers in get_cartesian_powers(1):
                        p_objs = [PrimitiveGaussian(e, c, center, powers) for e, c in prims_p]
                        label = f"P_{powers}"
                        basis_funcs.append(BasisFunction(p_objs, label, atom_idx))
                else:
                    L = 1 if shell_type == 2 else shell_type
                    for powers in get_cartesian_powers(L):
                        p_objs = [PrimitiveGaussian(e, c, center, powers) for e, c in prims_s]
                        l_char = {0:'S', 1:'P', 2:'D', 3:'F', 4:'G'}.get(L, 'X')
                        label = f"{l_char}_{powers}"
                        basis_funcs.append(BasisFunction(p_objs, label, atom_idx))
        
        return basis_funcs

# =============================================================================
#  EXAMPLE SYSTEMS
# =============================================================================

def example_h2():
    """H2 molecule with STO-3G."""
    print("\n" + "="*70)
    print("EXAMPLE: H₂ Molecule (STO-3G)")
    print("="*70)
    
    r = 1.4  # Bohr
    coords = [[0.0, 0.0, 0.0], [0.0, 0.0, r]]
    atoms = [Atom(1, c, f"H{i+1}") for i, c in enumerate(coords)]
    
    # Build basis
    bfs = []
    for i, atom in enumerate(atoms):
        bfs.extend(BasisParser.parse(BasisSetLibrary.STO3G_H, [atom.coords], [i]))
    
    print(f"\nMolecule: H₂")
    print(f"Bond length: {r:.3f} Bohr ({r*BOHR_TO_ANGSTROM:.3f} Å)")
    print(f"Basis set: STO-3G")
    print(f"Basis functions: {len(bfs)}")
    
    solver = SCFSolver(bfs, atoms, n_electrons=2)
    E = solver.run(tolerance=1e-8, use_diis=True)
    
    print("\n" + "="*70)
    print(f"FINAL RESULT: {E:.10f} Ha")
    print(f"Expected (STO-3G): ~-1.117 Ha")
    print("="*70)
    
    return solver

def example_h2o():
    """H2O molecule with STO-3G."""
    print("\n" + "="*70)
    print("EXAMPLE: H₂O Molecule (STO-3G)")
    print("="*70)
    
    # Geometry: O at origin, H atoms at ~104.5° angle, bond length 0.96 Å
    r_oh = 0.96 / BOHR_TO_ANGSTROM  # Convert to Bohr
    angle = 104.5 * math.pi / 180.0
    
    coords = [
        [0.0, 0.0, 0.0],  # O
        [r_oh * math.sin(angle/2), 0.0, r_oh * math.cos(angle/2)],  # H1
        [-r_oh * math.sin(angle/2), 0.0, r_oh * math.cos(angle/2)]   # H2
    ]
    
    atoms = [
        Atom(8, coords[0], "O"),
        Atom(1, coords[1], "H1"),
        Atom(1, coords[2], "H2")
    ]
    
    # Build basis
    bfs = []
    # Oxygen
    bfs.extend(BasisParser.parse(BasisSetLibrary.STO3G_C, [coords[0]], [0]))
    # Hydrogens
    for i in range(1, 3):
        bfs.extend(BasisParser.parse(BasisSetLibrary.STO3G_H, [coords[i]], [i]))
    
    print(f"\nMolecule: H₂O")
    print(f"O-H bond length: {r_oh*BOHR_TO_ANGSTROM:.3f} Å")
    print(f"H-O-H angle: {angle*180/math.pi:.1f}°")
    print(f"Basis set: STO-3G")
    print(f"Basis functions: {len(bfs)}")
    
    solver = SCFSolver(bfs, atoms, n_electrons=10)
    E = solver.run(tolerance=1e-7, use_diis=True, damping=0.0)
    
    print("\n" + "="*70)
    print(f"FINAL RESULT: {E:.10f} Ha")
    print(f"Expected (STO-3G): ~-75.0 Ha")
    print("="*70)
    
    return solver

def example_h_atom_uhf():
    """Hydrogen atom using UHF."""
    print("\n" + "="*70)
    print("EXAMPLE: H Atom (UHF, doublet)")
    print("="*70)
    
    coords = [[0.0, 0.0, 0.0]]
    atoms = [Atom(1, coords[0], "H")]
    
    bfs = BasisParser.parse(BasisSetLibrary.STO3G_H, coords, [0])
    
    print(f"\nBasis functions: {len(bfs)}")
    
    solver = UHFSolver(bfs, atoms, n_alpha=1, n_beta=0)
    E = solver.run(tolerance=1e-8)
    
    print("\n" + "="*70)
    print(f"FINAL RESULT: {E:.10f} Ha")
    print(f"Expected: ~-0.466 Ha (STO-3G)")
    print("="*70)
    
    return solver

# =============================================================================
#  MAIN
# =============================================================================

def main():
    """Main driver program."""
    print("\n" + "="*70)
    print(" PROFESSIONAL HARTREE-FOCK SCF PROGRAM v4.0")
    print(" Enhanced Edition with Full Diagnostics")
    print("="*70)
    
    # Run examples
    print("\n" + "█"*70)
    print("█ TEST 1: H₂ MOLECULE")
    print("█"*70)
    solver_h2 = example_h2()
    
    print("\n\n" + "█"*70)
    print("█ TEST 2: H ATOM (UHF)")
    print("█"*70)
    solver_h_uhf = example_h_atom_uhf()
    
    # Uncomment for H2O test (takes longer due to more basis functions)
    # print("\n\n" + "█"*70)
    # print("█ TEST 3: H₂O MOLECULE")
    # print("█"*70)
    # solver_h2o = example_h2o()
    
    print("\n" + "="*70)
    print("ALL TESTS COMPLETED")
    print("="*70)

if __name__ == "__main__":
    main()
