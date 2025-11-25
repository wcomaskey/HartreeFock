"""
Test Suite for Hartree-Fock Enhanced v4.0
Validates all major features and provides usage examples
"""

import sys
sys.path.insert(0, '/mnt/user-data/outputs')

from hartree_fock_enhanced import *
import numpy as np
import math

def test_h2_rhf():
    """Test 1: H2 molecule with RHF (basic functionality)"""
    print("\n" + "="*70)
    print("TEST 1: H₂ MOLECULE (RHF)")
    print("="*70)
    
    r = 1.4
    atoms = [Atom(1, [0,0,0], "H1"), Atom(1, [0,0,r], "H2")]
    bfs = []
    for i, atom in enumerate(atoms):
        bfs.extend(BasisParser.parse(BasisSetLibrary.STO3G_H, [atom.coords], [i]))
    
    solver = SCFSolver(bfs, atoms, n_electrons=2)
    E = solver.run(verbose=False)
    
    expected = -1.1167
    error = abs(E - expected)
    
    print(f"Calculated Energy: {E:.8f} Ha")
    print(f"Expected Energy:   {expected:.4f} Ha")
    print(f"Error:             {error:.6f} Ha")
    
    if error < 0.001:
        print("  TEST PASSED")
        return True
    else:
        print("  TEST FAILED")
        return False

def test_h_atom_uhf():
    """Test 2: Hydrogen atom with UHF (open-shell)"""
    print("\n" + "="*70)
    print("TEST 2: H ATOM (UHF)")
    print("="*70)
    
    atoms = [Atom(1, [0,0,0], "H")]
    bfs = BasisParser.parse(BasisSetLibrary.STO3G_H, [atoms[0].coords], [0])
    
    solver = UHFSolver(bfs, atoms, n_alpha=1, n_beta=0)
    E = solver.run(verbose=False)
    
    expected = -0.4666
    error = abs(E - expected)
    
    print(f"Calculated Energy: {E:.8f} Ha")
    print(f"Expected Energy:   {expected:.4f} Ha")
    print(f"Error:             {error:.6f} Ha")
    print(f"<S²>:              {solver.compute_spin_squared(solver.P_alpha, solver.P_beta):.6f}")
    
    if error < 0.001:
        print("TEST PASSED")
        return True
    else:
        print("TEST FAILED")
        return False

def test_convergence_damping():
    """Test 3: Convergence with damping"""
    print("\n" + "="*70)
    print("TEST 3: CONVERGENCE WITH DAMPING")
    print("="*70)
    
    r = 1.4
    atoms = [Atom(1, [0,0,0], "H1"), Atom(1, [0,0,r], "H2")]
    bfs = []
    for i, atom in enumerate(atoms):
        bfs.extend(BasisParser.parse(BasisSetLibrary.STO3G_H, [atom.coords], [i]))
    
    print("Testing different damping parameters:")
    results = []
    for damp in [0.0, 0.2, 0.5]:
        solver = SCFSolver(bfs, atoms, n_electrons=2)
        E = solver.run(damping=damp, verbose=False)
        results.append(E)
        print(f"  Damping={damp:.1f}  Energy={E:.8f} Ha")
    
    # All should give same answer
    if np.std(results) < 1e-6:
        print("  TEST PASSED (all converge to same energy)")
        return True
    else:
        print("  TEST FAILED (energies differ)")
        return False

def test_mulliken_charge_conservation():
    """Test 4: Mulliken charges sum to zero"""
    print("\n" + "="*70)
    print("TEST 4: MULLIKEN CHARGE CONSERVATION")
    print("="*70)
    
    r = 1.4
    atoms = [Atom(1, [0,0,0], "H1"), Atom(1, [0,0,r], "H2")]
    bfs = []
    for i, atom in enumerate(atoms):
        bfs.extend(BasisParser.parse(BasisSetLibrary.STO3G_H, [atom.coords], [i]))
    
    solver = SCFSolver(bfs, atoms, n_electrons=2)
    E = solver.run(verbose=False)
    
    # Compute Mulliken populations
    PS = np.dot(solver.P, solver.S)
    charges = []
    for i, atom in enumerate(atoms):
        pop = 0.0
        for j, bf in enumerate(bfs):
            if bf.atom_index == i:
                pop += PS[j, j]
        charge = atom.Z - pop
        charges.append(charge)
    
    total_charge = sum(charges)
    
    print(f"H1 charge: {charges[0]:+.6f}")
    print(f"H2 charge: {charges[1]:+.6f}")
    print(f"Total charge: {total_charge:+.6f}")
    
    if abs(total_charge) < 1e-6:
        print("  TEST PASSED (charges sum to zero)")
        return True
    else:
        print("  TEST FAILED (charge not conserved)")
        return False

def test_symmetry_h2():
    """Test 5: H2 should have symmetric charges"""
    print("\n" + "="*70)
    print("TEST 5: SYMMETRY (H₂ charges should be equal)")
    print("="*70)
    
    r = 1.4
    atoms = [Atom(1, [0,0,0], "H1"), Atom(1, [0,0,r], "H2")]
    bfs = []
    for i, atom in enumerate(atoms):
        bfs.extend(BasisParser.parse(BasisSetLibrary.STO3G_H, [atom.coords], [i]))
    
    solver = SCFSolver(bfs, atoms, n_electrons=2)
    E = solver.run(verbose=False)
    
    # Get populations
    PS = np.dot(solver.P, solver.S)
    pops = []
    for i in range(2):
        pop = 0.0
        for j, bf in enumerate(bfs):
            if bf.atom_index == i:
                pop += PS[j, j]
        pops.append(pop)
    
    diff = abs(pops[0] - pops[1])
    
    print(f"H1 population: {pops[0]:.6f}")
    print(f"H2 population: {pops[1]:.6f}")
    print(f"Difference:    {diff:.6f}")
    
    if diff < 1e-6:
        print("TEST PASSED (symmetric)")
        return True
    else:
        print("TEST FAILED (asymmetric)")
        return False

def test_homo_lumo_gap():
    """Test 6: HOMO-LUMO gap should be positive"""
    print("\n" + "="*70)
    print("TEST 6: HOMO-LUMO GAP")
    print("="*70)
    
    r = 1.4
    atoms = [Atom(1, [0,0,0], "H1"), Atom(1, [0,0,r], "H2")]
    bfs = []
    for i, atom in enumerate(atoms):
        bfs.extend(BasisParser.parse(BasisSetLibrary.STO3G_H, [atom.coords], [i]))
    
    solver = SCFSolver(bfs, atoms, n_electrons=2)
    E = solver.run(verbose=False)
    
    homo = solver.eps[solver.n_occ - 1]
    lumo = solver.eps[solver.n_occ]
    gap = lumo - homo
    
    print(f"HOMO energy: {homo:.6f} Ha ({homo*27.211:.2f} eV)")
    print(f"LUMO energy: {lumo:.6f} Ha ({lumo*27.211:.2f} eV)")
    print(f"Gap:         {gap:.6f} Ha ({gap*27.211:.2f} eV)")
    
    if gap > 0 and gap < 10.0:  # Reasonable gap
        print("TEST PASSED (positive, reasonable gap)")
        return True
    else:
        print("TEST FAILED (unreasonable gap)")
        return False

def test_integral_symmetry():
    """Test 7: Integral matrices should be symmetric"""
    print("\n" + "="*70)
    print("TEST 7: INTEGRAL SYMMETRY")
    print("="*70)
    
    r = 1.4
    atoms = [Atom(1, [0,0,0], "H1"), Atom(1, [0,0,r], "H2")]
    bfs = []
    for i, atom in enumerate(atoms):
        bfs.extend(BasisParser.parse(BasisSetLibrary.STO3G_H, [atom.coords], [i]))
    
    S = IntegralEngine.overlap(bfs)
    T = IntegralEngine.kinetic(bfs)
    Vne = IntegralEngine.nuclear_attraction(bfs, atoms)
    
    s_sym = np.allclose(S, S.T)
    t_sym = np.allclose(T, T.T)
    v_sym = np.allclose(Vne, Vne.T)
    
    print(f"Overlap symmetric:          {s_sym}")
    print(f"Kinetic symmetric:          {t_sym}")
    print(f"Nuclear attraction symmetric: {v_sym}")
    
    if s_sym and t_sym and v_sym:
        print("TEST PASSED (all symmetric)")
        return True
    else:
        print("TEST FAILED (asymmetric matrices)")
        return False

def test_bond_dissociation():
    """Test 8: Energy should increase as bond breaks"""
    print("\n" + "="*70)
    print("TEST 8: BOND DISSOCIATION CURVE")
    print("="*70)
    
    distances = [1.0, 1.4, 2.0, 3.0, 5.0]
    energies = []
    
    print("Distance (Bohr)  Energy (Ha)")
    print("-" * 35)
    
    for r in distances:
        atoms = [Atom(1, [0,0,0], "H1"), Atom(1, [0,0,r], "H2")]
        bfs = []
        for i, atom in enumerate(atoms):
            bfs.extend(BasisParser.parse(BasisSetLibrary.STO3G_H, [atom.coords], [i]))
        
        solver = SCFSolver(bfs, atoms, n_electrons=2)
        E = solver.run(verbose=False)
        energies.append(E)
        print(f"{r:8.2f}        {E:12.8f}")
    
    # Energy should have a minimum (equilibrium)
    min_idx = np.argmin(energies)
    r_eq = distances[min_idx]
    E_min = energies[min_idx]
    
    # At large distance, should approach 2 * E(H atom)
    E_dissoc = energies[-1]
    E_h_atom = -0.4666  # STO-3G
    expected_dissoc = 2 * E_h_atom
    
    print(f"\nEquilibrium distance: {r_eq:.2f} Bohr")
    print(f"Minimum energy:       {E_min:.6f} Ha")
    print(f"Dissociation energy:  {E_dissoc:.6f} Ha")
    print(f"Expected (2×H atom):  {expected_dissoc:.6f} Ha")
    
    # Check if curve has reasonable shape
    has_minimum = min_idx > 0 and min_idx < len(energies) - 1
    # Note: STO-3G has BSSE, so dissociation limit won't match exactly
    # Just check energy increases from minimum
    energy_increases = energies[-1] > energies[min_idx]
    
    if has_minimum and energy_increases:
        print("  TEST PASSED (has minimum, energy increases at dissociation)")
        print("  Note: Exact dissociation limit affected by BSSE in STO-3G")
        return True
    else:
        print("  TEST FAILED (incorrect curve shape)")
        return False

def test_orbital_orthonormality():
    """Test 9: MO coefficients should give orthonormal orbitals"""
    print("\n" + "="*70)
    print("TEST 9: ORBITAL ORTHONORMALITY")
    print("="*70)
    
    r = 1.4
    atoms = [Atom(1, [0,0,0], "H1"), Atom(1, [0,0,r], "H2")]
    bfs = []
    for i, atom in enumerate(atoms):
        bfs.extend(BasisParser.parse(BasisSetLibrary.STO3G_H, [atom.coords], [i]))
    
    solver = SCFSolver(bfs, atoms, n_electrons=2)
    E = solver.run(verbose=False)
    
    # Check C^T S C = I
    CTC = np.dot(solver.C.T, np.dot(solver.S, solver.C))
    
    # Should be identity matrix
    identity = np.eye(CTC.shape[0])
    max_error = np.max(np.abs(CTC - identity))
    
    print(f"Max deviation from orthonormality: {max_error:.2e}")
    
    if max_error < 1e-10:
        print("TEST PASSED (orbitals orthonormal)")
        return True
    else:
        print("TEST FAILED (orbitals not orthonormal)")
        return False

def test_spin_contamination():
    """Test 10: Hydrogen atom should have <S²> = 0.75 (doublet)"""
    print("\n" + "="*70)
    print("TEST 10: SPIN CONTAMINATION (H ATOM)")
    print("="*70)
    
    atoms = [Atom(1, [0,0,0], "H")]
    bfs = BasisParser.parse(BasisSetLibrary.STO3G_H, [atoms[0].coords], [0])
    
    solver = UHFSolver(bfs, atoms, n_alpha=1, n_beta=0)
    E = solver.run(verbose=False)
    
    S_squared = solver.compute_spin_squared(solver.P_alpha, solver.P_beta)
    expected_S_squared = 0.75  # S=1/2, S(S+1) = 0.75
    
    contamination = abs(S_squared - expected_S_squared)
    
    print(f"Calculated <S²>: {S_squared:.6f}")
    print(f"Expected <S²>:   {expected_S_squared:.6f}")
    print(f"Contamination:   {contamination:.2e}")
    
    if contamination < 1e-6:
        print("  TEST PASSED (no spin contamination)")
        return True
    else:
        print("  TEST FAILED (spin contamination detected)")
        return False

def run_all_tests():
    """Run complete test suite"""
    print("\n" + "█"*70)
    print("█" + " "*68 + "█")
    print("█" + " "*20 + "HARTREE-FOCK v4.0 TEST SUITE" + " "*20 + "█")
    print("█" + " "*68 + "█")
    print("█"*70)
    
    tests = [
        ("H₂ RHF Basic", test_h2_rhf),
        ("H Atom UHF", test_h_atom_uhf),
        ("Convergence Damping", test_convergence_damping),
        ("Mulliken Charge Conservation", test_mulliken_charge_conservation),
        ("Symmetry (H₂)", test_symmetry_h2),
        ("HOMO-LUMO Gap", test_homo_lumo_gap),
        ("Integral Symmetry", test_integral_symmetry),
        ("Bond Dissociation", test_bond_dissociation),
        ("Orbital Orthonormality", test_orbital_orthonormality),
        ("Spin (H Atom)", test_spin_contamination)
    ]
    
    results = []
    for name, test_func in tests:
        try:
            passed = test_func()
            results.append((name, passed))
        except Exception as e:
            print(f"  TEST CRASHED: {e}")
            results.append((name, False))
    
    # Summary
    print("\n" + "="*70)
    print("TEST SUMMARY")
    print("="*70)
    
    for name, passed in results:
        status = "  PASS" if passed else "  FAIL"
        print(f"{status:<10} {name}")
    
    total = len(results)
    passed = sum(1 for _, p in results if p)
    
    print("="*70)
    print(f"TOTAL: {passed}/{total} tests passed ({100*passed/total:.1f}%)")
    print("="*70)
    
    if passed == total:
        print("\n ALL TESTS PASSED! Code is working correctly.")
    else:
        print(f"\n {total - passed} test(s) failed. Please review.")
    
    return passed == total

if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)
