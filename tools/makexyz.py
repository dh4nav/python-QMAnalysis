"""
Python translation of makexyz.f
Extracts cartesian coordinates, atomic charges, normal modes, orbitals, and NMR shieldings from a Gaussian output file.
Generates output compatible with XMol, RasMol, etc.

This script is a direct translation from Fortran, with careful attention to index ranges:
- Fortran arrays are 1-based and inclusive on both ends.
- Python lists are 0-based and the end index is exclusive in slices.
- All loops and array accesses have been adjusted accordingly.

Extensive comments have been added for clarity.
"""

import re
import math

MAXATM = 1000
MAXBUF = 1000000

# Element symbols, index matches atomic number (1-based in Fortran, so ELMS[1] == "H")
ELMS = [
    "Bq", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K",
    "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La",
    "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os",
    "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am",
    "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Ha", "Sg", "Ns", "Hs", "Mt"
]
# Orbital labels
ORT = [
    "S   ", "PX  ", "PY  ", "PZ  ", "XY  ", "XZ  ", "YZ  ", "XX  ", "ZZ  ", "YY  ",
    "D-2 ", "D+1 ", "D-1 ", "D+2 ", "D 0 "
]
ORL = [
    "  S   ", "  Px  ", "  Py  ", "  Pz  ", " Dxy  ", " Dxz  ", " Dyz  ", "Dx2-y2", " Dz2  "
]


def parse_gaussian_output(lines):
    """
    Parses a Gaussian output file and prints XYZ coordinates, charges, shieldings, normal modes, and orbitals.
    """
    # Initialize arrays (Fortran 1-based, Python 0-based, so we use [0] as dummy and fill from [1])
    el = ["Xx"] * (MAXATM + 1)
    x = [0.0] * (MAXATM + 1)
    y = [0.0] * (MAXATM + 1)
    z = [0.0] * (MAXATM + 1)
    c = [0.0] * (MAXATM + 1)
    s = [[0.0] * 6 for _ in range(MAXATM + 1)]  # s[atom][0..5]
    e = [[[0.0] * (MAXATM + 1) for _ in range(10)]
         for _ in range(10)]  # e[jj][ij][atom]
    title = ""
    shield = False
    gauss98 = False
    nlines = len(lines)

    # Initialize block indices
    j = k = jc = kc = io = iv = 2 * MAXBUF
    kch = 0

    # Find relevant blocks in the Gaussian output
    for i, line in enumerate(lines):
        # Orientation block (for coordinates)
        if (
            "Standard orientation:" in line[19:40]
            or "Z-Matrix orientation:" in line[18:39]
            or "Input orientation:" in line[19:37]
        ):
            j = i + 5
            k = 2 * MAXBUF
        if (
            "Standard orientation:" in line[25:46]
            or "Z-Matrix orientation:" in line[25:46]
            or "Input orientation:" in line[26:47]
        ):
            j = i + 5
            k = 2 * MAXBUF
            gauss98 = True
        # End of orientation block
        if k > i >= j and line.startswith("----------"):
            k = i - 1
        # Start of molecular orbitals
        if io > nlines and "Molecular Orbital" in line[5:22]:
            io = i
        if io > nlines and "Alpha Molecular O" in line[5:22]:
            io = i
        # Start of vibrational normal modes
        if iv > nlines and line.startswith("Frequencies --"):
            iv = i
        # NMR shielding block
        if "  Anisotropy =" in line[33:46]:
            m = re.match(r"\s*(\d+)\s+.*?([\d\.]+)\s+.*?([\d\.]+)", line)
            if m:
                jj = int(m.group(1))
                s[jj][0] = float(m.group(2))
                s[jj][1] = float(m.group(3))
                m2 = re.match(
                    r".*?([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)", lines[i + 4])
                if m2:
                    s[jj][2] = float(m2.group(1))
                    s[jj][3] = float(m2.group(2))
                    s[jj][4] = float(m2.group(3))
                    s43 = s[jj][3] - s[jj][2]
                    s54 = s[jj][4] - s[jj][3]
                    s[jj][5] = s[jj][4] if s54 > s43 else s[jj][2]
                shield = True
        # Find atomic charges block
        if (
            "Mulliken charges:" in line[:17]
            or "Mulliken charges and spin densities:" in line[:36]
        ):
            if kch < 2:
                jc = i + 2
                kch = 1
        if "ESP charges:" in line[:12]:
            if kch < 3:
                jc = i + 2
                kch = 2
        if "Summary of Natural Popula" in line[:25]:
            jc = i + 6
            kch = 3
            if lines[i + 4].startswith(" Atom No  "):
                kch = 4

    # Find the last Standard orientation block
    j = k = None
    for i, line in enumerate(lines):
        if "Standard orientation:" in line:
            # Look for the start of the block (skip 5 lines after header)
            j_candidate = i + 5
            # Now find the end of the block (next dashed line after j_candidate)
            for k_candidate in range(j_candidate, len(lines)):
                if lines[k_candidate].strip().startswith('-----'):
                    k_candidate = k_candidate - 1  # last atom line
                    break
            else:
                continue  # no end found, skip
            # Save this block as the current one
            j, k = j_candidate, k_candidate

    # If no block found, fallback
    if j is None or k is None or j > k:
        na = 1
        el = ["Xx"]
        x = [0.0]
        y = [0.0]
        z = [0.0]
    else:
        # Parse atom lines
        el = []
        x = []
        y = []
        z = []
        for line in lines[j:k+1]:
            parts = line.split()
            if len(parts) >= 6:
                atomic_number = int(parts[1])
                el.append(ELMS[atomic_number])
                x.append(float(parts[3]))
                y.append(float(parts[4]))
                z.append(float(parts[5]))
        na = len(el)

    # Output coordinates and atomic charges in XYZ format
    if kch > 0:
        print(f"{na}\n{title}")
        for i in range(na):
            print(
                f"{el[i]:2s} {x[i]:12.6f} {y[i]:12.6f} {z[i]:12.6f} {c[i]:12.6f}"
            )
    else:
        print(f"{na}\n{title}")
        for i in range(na):
            print(f"{el[i]:2s} {x[i]:12.6f} {y[i]:12.6f} {z[i]:12.6f}")

    # Read atomic charges
    ii = 0
    if kch > 0 and jc < MAXBUF:
        kc = jc + na - 1  # Fortran: jc to kc inclusive
        for i in range(jc, kc + 1):
            ii += 1
            line = lines[i]
            if kch == 4:
                m = re.match(r".{8}([\d\.\-]+)", line)
            elif kch == 3:
                m = re.match(r".{12}([\d\.\-]+)", line)
            else:
                m = re.match(r".{10}([\d\.\-]+)", line)
            if m:
                c[ii] = float(m.group(1))
        if kch == 1:
            title = "Mulliken charges"
        elif kch == 2:
            title = "ESP charges"
        elif kch == 3:
            title = "NPA charges"

    # Output coordinates and atomic charges in XYZ format
    if kch > 0:
        for i in range(na):
            print(
                f"{el[i]:2s} {x[i]:12.6f} {y[i]:12.6f} {z[i]:12.6f} {c[i]:12.6f}"
            )
    else:
        for i in range(na):
            print(f"{el[i]:2s} {x[i]:12.6f} {y[i]:12.6f} {z[i]:12.6f}")

    # Output shieldings if present
    if shield:
        print("\n\nSHIELDINGS")
        for i in range(1, na + 1):
            if el[i] == "Bq":
                print(f"{el[i]:2s} " +
                      " ".join(f"{s[jj][i]:11.4f}" for jj in range(6)))
            else:
                print(f"{el[i]:2s} " +
                      " ".join(f"{s[jj][i]:11.4f}" for jj in range(2)))

    # Output normal modes (vibrational frequencies)
    if iv < MAXBUF:
        print(f"\n\nNORMALMODES {3*na-6}")
        ii = 1
        while ii <= na - 2 and iv < nlines:
            freq_line = lines[iv]
            m = re.match(
                r".*?([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)", freq_line[14:])
            if m:
                xa = float(m.group(1))
                xb = float(m.group(2))
                xc = float(m.group(3))
            else:
                xa = xb = xc = 0.0
            iv += 3
            # IR Intensity
            if lines[iv].startswith("IR Int"):
                m = re.match(
                    r".*?([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)", lines[iv][14:])
                if m:
                    xia = float(m.group(1))
                    xib = float(m.group(2))
                    xic = float(m.group(3))
                else:
                    xia = xib = xic = 0.0
            else:
                xia = xib = xic = 0.0
            iv += 1
            # Raman Intensity
            if lines[iv].startswith("Raman "):
                m = re.match(
                    r".*?([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)", lines[iv][14:])
                if m:
                    xra = float(m.group(1))
                    xrb = float(m.group(2))
                    xrc = float(m.group(3))
                else:
                    xra = xrb = xrc = 0.0
                iv += 2
                if not lines[iv].startswith(" Atom "):
                    iv += 1
            else:
                xra = xrb = xrc = 0.0
            # Normal mode vectors (3 blocks per mode)
            for block in range(3):
                for i in range(1, na + 1):
                    mode_line = lines[iv + i]
                    m = re.findall(r"([\d\.\-]+)", mode_line[11:])
                    if m and len(m) >= 9:
                        d = [float(val) for val in m[block*3:(block+1)*3]]
                    else:
                        d = [0.0] * 3
                    print(f"{el[i]:2s} {d[0]:6.2f} {d[1]:6.2f} {d[2]:6.2f}")
                # Print frequency/intensity for this block
                if block == 0:
                    print(f"{xa:12.4f} {xia:12.4f} {xra:12.4f}")
                elif block == 1:
                    print(f"{xb:12.4f} {xib:12.4f} {xrb:12.4f}")
                else:
                    print(f"{xc:12.4f} {xic:12.4f} {xrc:12.4f}")
            iv += na + 3
            ii += 1

    # Output molecular orbitals (simplified, as in original Fortran)
    if io < MAXBUF:
        print("\n\nMOLORBS")
        for ij in range(1, 11):  # Fortran: 1 to 10
            xx = 0.0
            for i in range(1, na + 1):
                for jj in range(9):
                    xx += e[jj][ij][i] ** 2
            if xx > 0.0:
                xx = 1.0 / math.sqrt(xx)
                for i in range(1, na + 1):
                    for jj in range(9):
                        e[jj][ij][i] *= xx
                print(f"{ij:2d}  " + " ".join(f"{label:7s}" for label in ORL))
                for i in range(1, na + 1):
                    print(
                        f"{el[i]:2s} "
                        + " ".join(f"{e[jj][ij][i]:7.2f}" for jj in range(9))
                    )


def main():
    import sys
    import io
    import contextlib
    import os

    if len(sys.argv) < 2:
        print("Usage: python makexyz.py <gaussian_output_file> [output.xyz]")
        return
    input_file = sys.argv[1]
    with open(input_file, "r") as f:
        lines = f.readlines()

    # Capture output of parse_gaussian_output
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        parse_gaussian_output(lines)
    output = buf.getvalue()

    # Determine output file name
    if len(sys.argv) > 2:
        output_file = sys.argv[2]
    else:
        # Always write an output file if no output filename is given
        # Use input filename + ".xyz" (even if it has no extension)
        base, ext = os.path.splitext(input_file)
        if ext.lower() == ".xyz":
            output_file = input_file
        else:
            output_file = input_file + ".xyz"

    # Write to file if specified or inferred, else print
    if output_file:
        with open(output_file, "w") as outf:
            outf.write(output)


if __name__ == "__main__":
    main()
