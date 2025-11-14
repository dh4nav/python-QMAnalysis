# GitHub Copilot Instructions for QMAnalysis

## Project Overview
This is a Python application for quantum mechanics/computational chemistry analysis. It processes atom coordinate files (XYZ format) and CSV files with fixed values, performs measurements and calculations, and outputs results as CSV files and plots.

**Key Principle**: Users are NOT required to know how to code. They interact with the software through YAML configuration files.

## Project Structure

### Core Modules (`src/qmanalysis/`)
- `yamlreader.py` - Parse YAML input files with validation
- `xyzreader.py` - Read and parse XYZ coordinate files
- `globalconstantsreader.py` - Read constants from CSV files
- `containers.py` - Data structures for storing atoms, coordinates, and measurements
- `measure.py` - Calculate distances, angles, and dihedral angles
- `customcalculationrunner.py` - Execute user-defined calculations
- `gaussianoutreader.py` - Parse Gaussian output files
- `main.py` - Main entry point

### Configuration Files
- `inputfile.yaml` / `inputfile2.yaml` - User-facing configuration files
- `constants.csv` - Fixed values for calculations
- `globals.csv` - Global constants

## YAML Configuration Schema

Users define analysis workflows in YAML with these sections:

### 1. Files Section
```yaml
files:
  - path: "benzene*.xyz"  # Supports glob patterns
    type: xyz
    name: benz-            # Logical name prefix
    glob: True
```

### 2. Substitutions Section
Allows dynamic atom index mapping across multiple files:
```yaml
substitutions:
  - name: S1              # Substitution variable name
    entries:
      - file: benz-*      # File pattern (supports glob)
        atom_index: 5     # 1-based atom index
        glob: True
      - file: bxyz
        atom_index: 8
      - file: [Li, Ca, Be, Mg]  # List of file names
        atom_index: 3
```

**Note**: The `file` field can be:
- A single filename/pattern string
- A list of filenames (for applying the same atom_index to multiple files)
- When using a list, `glob: True` applies the pattern to each item in the list

### 3. Measurements Section
```yaml
measurements:
  distance:               # Distance between two atoms
    - name: O-H-bond
      a: 6                # Atom indices (1-based)
      b: 7
    - name: Subs-bond
      a: S1               # Can reference substitution variables
      b: 10
  angle:                  # Angle between three atoms
    - name: H-O-H
      a: 3
      b: 6
      c: 7
  dihedral:               # Dihedral angle between four atoms
    - name: torsion1
      a: 2
      b: 3
      c: 6
      d: 7
```

### 4. Output Section
```yaml
output:
  - file:
    - path: ./out.csv
      type: csv
  - graph:
    - type: scatter_plot
      x: Subs-bond        # Measurement name for x-axis
      y: H-O-N           # Measurement name for y-axis
      file: ./SB-HON.tiff
```

## Code Conventions

### Atom Indexing
- **Always use 1-based indexing** in YAML files and user-facing interfaces
- Convert to 0-based internally when working with Python lists/arrays
- Document clearly when functions use 0-based vs 1-based indexing

### File Handling
- Support glob patterns for batch processing
- Use `Path` objects from `pathlib` for cross-platform compatibility
- Handle missing files gracefully with informative error messages

### Measurements
- Distances in Angstroms (Å)
- Angles in degrees (°)
- Use numpy for vector calculations

### Data Flow
1. Parse YAML configuration
2. Load coordinate files (XYZ, Gaussian output)
3. Resolve substitutions for each file
4. Calculate measurements for each file
5. Run custom calculations if defined
6. Generate output (CSV, plots)

### Error Handling
- Validate YAML schema early
- Provide clear error messages for missing files, invalid atom indices
- Handle glob patterns that match zero files
- Validate measurement references (e.g., substitution variables exist)

## Testing Conventions
- Test data in `testdata/` directory
- Use pytest for unit tests
- Integration tests should use real YAML configs from `testdata/`
- Mock file I/O when appropriate

## Dependencies
- `numpy` - Numerical calculations
- `pandas` - Data manipulation and CSV output
- `matplotlib` - Plotting
- `pyyaml` - YAML parsing
- Standard library: `pathlib`, `glob`, `typing`

## Documentation
- Use NumPy-style docstrings
- Document YAML schema in `docs/yaml_format_reference.rst`
- Keep README.md updated with usage examples

## Special Notes
- The project evolved from Fortran code (`makexyz.f`)
- Legacy code in `oldQMAnalysis/` for reference
- Some test data from real research (Rudi80_Ursula datasets)
- Support for Gaussian quantum chemistry software outputs
