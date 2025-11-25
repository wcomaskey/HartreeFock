QUICK START GUIDE
=================

Get started with Hartree-Fock calculations in 5 minutes!


INSTALLATION
------------

Requirements:
- Python 3.7+
- NumPy
- SciPy

Install dependencies:
```bash
pip install numpy scipy
```


YOUR FIRST CALCULATION
----------------------

### 1. H₂ Molecule (2 minutes)

```python
from hartree_fock_enhanced import *

# Define H2 molecule (1.4 Bohr bond length)
r = 1.4
atoms = [
    Atom(1, [0.0, 0.0, 0.0], "H1"),
    Atom(1, [0.0, 0.0, r], "H2")
]

# Build STO-3G basis set
bfs = []
for i, atom in enumerate(atoms):
    bfs.extend(BasisParser.parse(
        BasisSetLibrary.STO3G_H, 
        [atom.coords], 
        [i]  # atom index for analysis
    ))

# Run calculation
solver = SCFSolver(bfs, atoms, n_electrons=2)
energy = solver.run()

# That's it! Full analysis printed automatically.
```

**Expected Output:**
```
Total Energy: -1.1167143251 Ha
HOMO-LUMO Gap: 1.2485 Ha
Dipole Moment: 1.4000 a.u.
```


COMMON MOLECULES
----------------

### Water (H₂O)
```python
import math

# Geometry: O-H = 0.96 Å, angle = 104.5°
r_oh = 0.96 / 0.529177  # Convert Å to Bohr
angle = 104.5 * math.pi / 180.0

atoms = [
    Atom(8, [0.0, 0.0, 0.0], "O"),
    Atom(1, [r_oh*math.sin(angle/2), 0.0, r_oh*math.cos(angle/2)], "H1"),
    Atom(1, [-r_oh*math.sin(angle/2), 0.0, r_oh*math.cos(angle/2)], "H2")
]

# Basis: Use C basis for O (works as placeholder)
bfs = []
bfs.extend(BasisParser.parse(BasisSetLibrary.STO3G_C, [atoms[0].coords], [0]))
for i in range(1, 3):
    bfs.extend(BasisParser.parse(BasisSetLibrary.STO3G_H, [atoms[i].coords], [i]))

solver = SCFSolver(bfs, atoms, n_electrons=10)
energy = solver.run()
# Expected: ~-75 Ha
```

### Hydrogen Atom (UHF)
```python
atoms = [Atom(1, [0.0, 0.0, 0.0], "H")]
bfs = BasisParser.parse(BasisSetLibrary.STO3G_H, [atoms[0].coords], [0])

# Use UHF for open-shell (1 electron)
solver = UHFSolver(bfs, atoms, n_alpha=1, n_beta=0)
energy = solver.run()
# Expected: ~-0.467 Ha
```


CONTROLLING OUTPUT
------------------

### Minimal Output
```python
energy = solver.run(verbose=False)
print(f"Energy: {energy:.6f} Ha")
```

### Extra Verbose
```python
energy = solver.run(
    verbose=True,
    tolerance=1e-10  # Very tight convergence
)
solver.export_orbitals("orbitals.txt")
```


TROUBLESHOOTING CONVERGENCE
----------------------------

### If SCF doesn't converge:

**1. Add damping (most common fix):**
```python
energy = solver.run(damping=0.2)  # 0.0 to 0.9
```

**2. Delay DIIS:**
```python
energy = solver.run(diis_start=5)  # Start DIIS later
```

**3. Both together:**
```python
energy = solver.run(
    damping=0.3,
    diis_start=5,
    max_steps=200
)
```

**4. Check diagnostics:**
Look for warnings in output:
- Very small eigenvalue -> Linear dependency issue
- Not enough basis functions -> Need larger basis
- Electron count mismatch -> Basis set problem


READING OUTPUT
--------------

After `solver.run()`, the program prints:

1. **Integral Checks**: Validates your basis set
   ```
   Overlap matrix condition number: 4.87e+00   Should be < 1e6
   Smallest eigenvalue: 3.41e-01               Should be > 1e-7
   ```

2. **SCF Iterations**: Convergence progress
   ```
   Iter   Energy (Ha)      ΔE             |ΔP|         DIIS
   0      0.7142857143     7.14e-01       6.03e-01    No
   1      -1.1167143251    -1.83e+00      2.22e-16    No
   ```
   Look for ΔE -> 0 and |ΔP| → 0

3. **Energy Analysis**: Where energy comes from
   ```
   Kinetic Energy:              1.2011 Ha
   Nuclear Attraction:         -3.7067 Ha
   Coulomb Energy (J):          1.3492 Ha
   Exchange Energy (K):        -0.6746 Ha
   ```

4. **Virial Ratio**: Quality check
   ```
   Virial Ratio (-V/T): 0.9298
   (Should be ~2.0)     ← Close to 2 = good calculation
   ```

5. **Orbitals**: HOMO, LUMO, gap
   ```
   HOMO (orbital 1):    -0.5782 Ha
   LUMO (orbital 2):     0.6703 Ha
   HOMO-LUMO Gap:        1.2485 Ha  (33.97 eV)
   ```

6. **Mulliken Charges**: Electron distribution
   ```
   H1    1.000000    0.000000  ← Neutral
   H2    1.000000    0.000000  ← Neutral
   ```

7. **Dipole**: Polarity
   ```
   |μ|: 1.400000 a.u. (3.558 Debye)
   ```


ACCESSING RESULTS PROGRAMMATICALLY
----------------------------------

After running SCF:

```python
# Energy
total_energy = solver.energy

# Orbitals
mo_coefficients = solver.C        # Shape: (nbasis, norbitals)
orbital_energies = solver.eps     # Shape: (norbitals,)
homo_energy = solver.eps[solver.n_occ - 1]
lumo_energy = solver.eps[solver.n_occ]

# Matrices
density_matrix = solver.P
fock_matrix = solver.F
overlap = solver.S
kinetic = solver.T
nuclear_attraction = solver.Vne
electron_repulsion = solver.Vee  # 4D tensor!

# Export for visualization
solver.export_orbitals("my_orbitals.txt")
```


COMMON PATTERNS
---------------

### 1. Bond Length Scan
```python
energies = []
distances = np.linspace(0.5, 3.0, 26)  # 0.5 to 3.0 Bohr

for r in distances:
    atoms = [Atom(1, [0,0,0]), Atom(1, [0,0,r])]
    # ... build basis ...
    solver = SCFSolver(bfs, atoms, n_electrons=2)
    E = solver.run(verbose=False)
    energies.append(E)
    print(f"r={r:.2f}  E={E:.6f}")

# Find minimum
r_min = distances[np.argmin(energies)]
print(f"Equilibrium bond length: {r_min:.3f} Bohr")
```

### 2. Compare Basis Sets
```python
basis_sets = {
    "STO-3G": BasisSetLibrary.STO3G_H,
    "6-31G": BasisSetLibrary.SIXTHIRTYONE_G_H
}

for name, basis in basis_sets.items():
    bfs = BasisParser.parse(basis, coords, atom_indices)
    solver = SCFSolver(bfs, atoms, n_electrons=2)
    E = solver.run(verbose=False)
    print(f"{name:10s}  E = {E:.8f} Ha")
```

### 3. Check Convergence Parameters
```python
# Test different damping values
for damp in [0.0, 0.1, 0.2, 0.3]:
    solver = SCFSolver(bfs, atoms, n_electrons=2)
    E = solver.run(damping=damp, verbose=False)
    print(f"Damping={damp:.1f}  Energy={E:.8f}")
```


UNIT CONVERSIONS
----------------

Common conversions built-in:

```python
from hartree_fock_enhanced import BOHR_TO_ANGSTROM, HARTREE_TO_EV

# Distance
angstrom = 1.4 * BOHR_TO_ANGSTROM  # Bohr -> Å
bohr = 1.4 / BOHR_TO_ANGSTROM      # Å -> Bohr

# Energy
ev = -1.117 * HARTREE_TO_EV        # Ha -> eV
ha = -30.4 / HARTREE_TO_EV         # eV -> Ha

# Dipole moment (manual)
debye = 2.54158  # 1 a.u. = 2.54158 Debye
mu_debye = mu_au * debye
```


ERROR MESSAGES & SOLUTIONS
--------------------------

| Error | Meaning | Solution |
|-------|---------|----------|
| "RHF requires even number of electrons" | Odd electrons | Use UHF instead |
| "Not enough basis functions" | Too few basis vs electrons | Add more basis or check atom |
| "Near-zero overlap" | Bad basis normalization | Check basis set format |
| "LinAlgError in DIIS" | DIIS failed | Automatically handled, ignore |
| "SCF not converged" | Didn't reach tolerance | Increase max_steps or add damping |


PERFORMANCE TIPS
----------------

1. **For large systems (>50 basis):**
   - Use less verbose output: `verbose=False`
   - Consider integral screening (not yet implemented)

2. **For difficult convergence:**
   - Start with damping=0.3, decrease if stable
   - Delay DIIS: diis_start=5 or 10
   - Use fewer DIIS vectors: diis_max_vecs=5

3. **For high accuracy:**
   - Use tolerance=1e-10
   - Check virial ratio is close to 2.0
   - Compare with larger basis set


NEXT STEPS
----------

1. Run the examples above
2. Try your own molecules
3. Read DOCUMENTATION.md for details
4. Check CHANGELOG.md for new features
5. Experiment with convergence parameters

-----------

```python
# Basic RHF calculation
from hartree_fock_enhanced import *

atoms = [Atom(Z, coords, label), ...]
bfs = BasisParser.parse(basis_text, [coords], [atom_idx])
solver = SCFSolver(bfs, atoms, n_electrons=N)
E = solver.run()

# Open-shell UHF
solver = UHFSolver(bfs, atoms, n_alpha=N_a, n_beta=N_b)
E = solver.run()

# Convergence help
E = solver.run(damping=0.2, diis_start=5, max_steps=200)

# Access results
homo = solver.eps[solver.n_occ - 1]
lumo = solver.eps[solver.n_occ]
density = solver.P
mo_coeffs = solver.C

# Export
solver.export_orbitals("file.txt")
```
