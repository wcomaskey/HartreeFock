# Jupyter Notebook Setup & Usage Guide

## Quick Start

### 1. Install Dependencies

```bash
# Basic requirements
pip install numpy scipy matplotlib jupyter

# For interactive 3D plots (optional)
pip install plotly

# For enhanced plotting (optional)
pip install seaborn pandas
```

### 2. Launch Notebook

```bash
# Navigate to the directory containing the notebook
cd /path/to/your/hartree_fock_directory

# Launch Jupyter
jupyter notebook HartreeFock_Interactive.ipynb
```

**Your browser will open automatically** showing the notebook.

### 3. Run the Notebook

**Option A: Run All Cells**
- Menu: `Kernel` → `Restart & Run All`
- This runs everything from start to finish

**Option B: Run Cell by Cell**
- Click on a cell
- Press `Shift + Enter` to run and move to next
- Or click the "Run" button in toolbar

**Option C: Run Selected Cells**
- Select cells you want to run
- Press `Shift + Enter` or click "Run"

## Notebook Structure

### Section 1: Historical Context
- **Read-only:** Pure markdown with history
- **No code to run**
- **Action needed:** Add your images to placeholders

### Section 2-5: Theory
- **Read-only:** Mathematical derivations
- **No code to run**
- **Interactive:** You can add LaTeX equations if desired

### Section 6-8: Interactive Examples
- **Code cells:** Run these!
- **Will perform calculations**
- **Generates plots**
- **Time:** 2-5 minutes total

### Section 9-10: Extensions & References
- **Read-only:** Discussion and bibliography
- **No code to run**

## Common Tasks

### Task 1: Run a Quick H₂ Calculation

Find this cell (Section 7):
```python
# Define H2 molecule
r = 1.4  # Bond length in Bohr
...
```

**Just run it!** You'll see:
- SCF iterations
- Final energy
- Mulliken populations
- Dipole moment
- All diagnostics

### Task 2: Generate Potential Energy Curve

Find this cell:
```python
# Scan bond lengths
distances = np.linspace(0.5, 4.0, 30)
...
```

**Run it** - takes ~30 seconds
- Calculates energy at 30 different bond lengths
- Finds equilibrium geometry
- Plots results

### Task 3: Visualize Molecular Orbitals

Find the orbital visualization cells and run them.
You'll get:
- HOMO and LUMO contour plots
- 2D cross-sections
- Bonding/antibonding character

### Task 4: Calculate Water Molecule

Find the H₂O example and run:
```python
# H2O geometry
r_oh = 0.96 / BOHR_TO_ANGSTROM
...
```

Takes ~10 seconds, gives full H₂O analysis.

## Modifying Examples

### Change Bond Length:
```python
# Original
r = 1.4

# Modified
r = 1.6  # Try different lengths
```

### Change Basis Set:
```python
# Current: STO-3G
BasisSetLibrary.STO3G_H

# Future: Add your own basis sets
```

### Adjust SCF Parameters:
```python
solver.run(
    tolerance=1e-8,    # Tighter convergence
    damping=0.2,       # Add damping if needed
    use_diis=True,     # DIIS acceleration
    verbose=True       # Show details
)
```

### Modify Plots:
```python
# Change figure size
plt.rcParams['figure.figsize'] = (16, 10)

# Change colors
plt.plot(x, y, 'r-')  # Red line

# Add grid
ax.grid(True, alpha=0.5)

# Change title
ax.set_title('My Custom Title', fontsize=18)
```

## Adding Images

### Method 1: Markdown Cell

```markdown
![Description](images/my_portrait.png)
*Caption: Douglas Hartree at Cambridge, 1935*
```

### Method 2: HTML (with sizing)

```html
<img src="images/hartree_portrait.jpg" width="400">
<p><em>Caption: Douglas Hartree (1897-1958)</em></p>
```

### Method 3: From URL

```markdown
![Description](https://example.com/image.png)
```

### Image Organization:

Create an `images/` folder:
```
Your Directory/
├── HartreeFock_Interactive.ipynb
├── images/
│   ├── portraits/
│   │   ├── hartree.jpg
│   │   ├── fock.jpg
│   │   ├── slater.jpg
│   │   └── ...
│   ├── equipment/
│   │   ├── mechanical_calculator.jpg
│   │   ├── ibm7090.jpg
│   │   └── ...
│   └── diagrams/
│       ├── scf_flowchart.png
│       ├── orbitals.png
│       └── ...
```

Then reference as:
```markdown
![Hartree](images/portraits/hartree.jpg)
```

## Keyboard Shortcuts

### Cell Operations:
- `Shift + Enter`: Run cell, move to next
- `Ctrl + Enter`: Run cell, stay on current
- `Alt + Enter`: Run cell, insert new below

### Cell Type:
- `Y`: Change to code cell
- `M`: Change to markdown cell
- `R`: Change to raw cell

### Cell Editing:
- `A`: Insert cell above
- `B`: Insert cell below
- `DD`: Delete cell (press D twice)
- `Z`: Undo delete
- `C`: Copy cell
- `V`: Paste cell

### Kernel Operations:
- `00`: Restart kernel (press 0 twice)
- `I I`: Interrupt kernel (press I twice)

### Mode:
- `Enter`: Enter edit mode (green border)
- `Esc`: Enter command mode (blue border)

## Troubleshooting

### Issue: "Module not found"

**Solution:**
```bash
# Make sure hartree_fock_enhanced.py is in same directory
# Or add to Python path
import sys
sys.path.insert(0, '/path/to/your/code')
```

### Issue: Plots don't show

**Solution:**
```python
# Add at top of notebook
%matplotlib inline

# Or for interactive plots
%matplotlib notebook
```

### Issue: Kernel dies during calculation

**Possible causes:**
- System out of memory
- Calculation too large
- Infinite loop

**Solutions:**
- Restart kernel: Kernel → Restart
- Use smaller basis set
- Check your code modifications

### Issue: Can't edit cells

**Solution:**
- You're in command mode (blue border)
- Press `Enter` to edit
- Or double-click the cell

### Issue: LaTeX not rendering

**Solution:**
```python
# Install required package
pip install jupyter_latex_envs

# Or use MathJax (usually works by default)
```

## Exporting

### To HTML:
```bash
jupyter nbconvert --to html HartreeFock_Interactive.ipynb
```

### To PDF:
```bash
# Requires LaTeX installation
jupyter nbconvert --to pdf HartreeFock_Interactive.ipynb
```

### To Python Script:
```bash
jupyter nbconvert --to script HartreeFock_Interactive.ipynb
```

### To Slides:
```bash
jupyter nbconvert --to slides HartreeFock_Interactive.ipynb --post serve
```

## Best Practices

### 1. Save Frequently
- `Ctrl + S` or click save icon
- Notebooks auto-save but better safe than sorry

### 2. Clear Outputs Before Committing
- Menu: `Cell` → `All Output` → `Clear`
- Reduces file size
- Removes old/incorrect results

### 3. Restart & Run All Periodically
- Ensures cells work in order
- Catches dependency issues
- Verifies reproducibility

### 4. Comment Your Modifications
```python
# Original value
# r = 1.4

# My modification: testing longer bond
r = 2.0  # Increased to study dissociation
```

### 5. Make Copies for Experiments
```bash
cp HartreeFock_Interactive.ipynb MyExperiment.ipynb
```

Then you can freely modify without losing original.

## Advanced Usage

### Running from Command Line:

```bash
# Execute notebook without opening browser
jupyter nbconvert --execute --to notebook \
    --inplace HartreeFock_Interactive.ipynb
```

### Using Virtual Environments:

```bash
# Create environment
python -m venv hf_env

# Activate
source hf_env/bin/activate  # Linux/Mac
# or
hf_env\Scripts\activate  # Windows

# Install packages
pip install numpy scipy matplotlib jupyter

# Run notebook
jupyter notebook
```

### Adding New Cells:

**Code cell for your own molecule:**
```python
# My molecule: Ammonia (NH3)
import numpy as np

# Define geometry
# ... your code ...
```

**Markdown cell for notes:**
```markdown
## My Notes

This is where I tested different parameters...
```

## Performance Tips

### For Large Calculations:

1. **Use verbose=False:**
```python
solver.run(verbose=False)  # Less output = faster
```

2. **Limit iterations:**
```python
solver.run(max_steps=20)  # Don't waste time if not converging
```

3. **Profile your code:**
```python
import time
start = time.time()
# ... your calculation ...
print(f"Time: {time.time() - start:.2f} seconds")
```

### For Many Plots:

```python
# Clear figures to save memory
plt.close('all')
```

## Getting Help

### In Notebook:
```python
# Get help on function
?solver.run

# Or
help(solver.run)

# View source
??solver.run
```

### Online Resources:
- Jupyter docs: https://jupyter.org/documentation
- Matplotlib gallery: https://matplotlib.org/stable/gallery/
- NumPy docs: https://numpy.org/doc/
- SciPy docs: https://docs.scipy.org/

### Common Questions:

**Q: How do I add my own molecule?**
A: See Example 3 (H₂O) as template. Define atoms and basis, then run SCF.

**Q: Can I use different basis sets?**
A: Yes! Add to BasisSetLibrary or parse from file. See DOCUMENTATION.md.

**Q: How do I compare with experimental data?**
A: Add your experimental values and plot alongside calculated values.

**Q: Can I export molecular orbitals?**
A: Yes! Use `solver.export_orbitals("filename.txt")`

## Summary

**To use this notebook:**
1. Install: `pip install numpy scipy matplotlib jupyter`
2. Launch: `jupyter notebook HartreeFock_Interactive.ipynb`
3. Run: `Kernel` → `Restart & Run All`
4. Explore: Modify parameters and re-run cells
5. Customize: Add your images and content
6. Share: Export to HTML or PDF

**The notebook is designed to be:**
- Self-contained (all code included)
- Interactive (modify and experiment)
- Educational (complete theory and history)
- Visual (plots and diagrams)
- Extensible (add your own content)

**Have fun exploring Hartree-Fock theory!** 

---

For more information, see:
- QUICKSTART.md - Quick tutorial
- DOCUMENTATION.md - Complete reference
- HartreeFock_History_Theory.md - Historical deep-dive
