"""
Python translation of makexyz.f
Extracts cartesian coordinates, atomic charges, normal modes, orbitals, and NMR shieldings from a Gaussian output file.
Generates output compatible with XMol, RasMol, etc.
"""

import re
import math

MAXATM = 1000
MAXBUF = 1000000
ELMS = [
    "Bq",
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
    "Xe",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn",
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
    "Md",
    "No",
    "Lr",
    "Rf",
    "Ha",
    "Sg",
    "Ns",
    "Hs",
    "Mt",
]
ORT = [
    "S   ",
    "PX  ",
    "PY  ",
    "PZ  ",
    "XY  ",
    "XZ  ",
    "YZ  ",
    "XX  ",
    "ZZ  ",
    "YY  ",
    "D-2 ",
    "D+1 ",
    "D-1 ",
    "D+2 ",
    "D 0 ",
]
ORL = [
    "  S   ",
    "  Px  ",
    "  Py  ",
    "  Pz  ",
    " Dxy  ",
    " Dxz  ",
    " Dyz  ",
    "Dx2-y2",
    " Dz2  ",
]


def parse_gaussian_output(lines):
    # Initialize variables
    el = ["Xx"] * MAXATM
    x = [0.0] * MAXATM
    y = [0.0] * MAXATM
    z = [0.0] * MAXATM
    c = [0.0] * MAXATM
    s = [[0.0] * 6 for _ in range(MAXATM)]
    e = [[[0.0] * 10 for _ in range(10)] for _ in range(MAXATM)]
    title = ""
    shield = False
    gauss98 = False
    nlines = len(lines)
    # Find orientation block
    j = k = jc = kc = io = iv = 2 * MAXBUF
    kch = 0
    for i, line in enumerate(lines):
        if (
            "Standard orientation:" in line[19:39]
            or "Z-Matrix orientation:" in line[18:38]
            or "Input orientation:" in line[19:36]
        ):
            j = i + 5
            k = 2 * MAXBUF
        if (
            "Standard orientation:" in line[25:45]
            or "Z-Matrix orientation:" in line[25:45]
            or "Input orientation:" in line[26:46]
        ):
            j = i + 5
            k = 2 * MAXBUF
            gauss98 = True
        if k > i and i >= j and line.startswith("----------"):
            k = i - 1
        if io > nlines and "Molecular Orbital" in line[5:21]:
            io = i
        if io > nlines and "Alpha Molecular O" in line[5:21]:
            io = i
        if iv > nlines and line.startswith("Frequencies --"):
            iv = i
        if "  Anisotropy =" in line[33:46]:
            m = re.match(r"\s*(\d+)\s+.*?([\d\.]+)\s+.*?([\d\.]+)", line)
            if m:
                jj = int(m.group(1))
                s[jj][0] = float(m.group(2))
                s[jj][1] = float(m.group(3))
                m2 = re.match(r".*?([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)", lines[i + 4])
                if m2:
                    s[jj][2] = float(m2.group(1))
                    s[jj][3] = float(m2.group(2))
                    s[jj][4] = float(m2.group(3))
                    s43 = s[jj][3] - s[jj][2]
                    s54 = s[jj][4] - s[jj][3]
                    s[jj][5] = s[jj][4] if s54 > s43 else s[jj][2]
                shield = True
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
    # Read atoms
    ii = 0
    if j < MAXBUF:
        if gauss98:
            for i in range(j, k):
                ii += 1
                m = re.match(
                    r".{12}(\d+).{17}([\d\.\-]+).{12}([\d\.\-]+).{12}([\d\.\-]+)",
                    lines[i],
                )
                if m:
                    jj = int(m.group(1))
                    x[ii] = float(m.group(2))
                    y[ii] = float(m.group(3))
                    z[ii] = float(m.group(4))
                    el[ii] = ELMS[jj]
                    if jj < 0 or jj > 109:
                        ii -= 1
        else:
            for i in range(j, k):
                ii += 1
                m = re.match(
                    r".{12}(\d+).{6}([\d\.\-]+).{12}([\d\.\-]+).{12}([\d\.\-]+)",
                    lines[i],
                )
                if m:
                    jj = int(m.group(1))
                    x[ii] = float(m.group(2))
                    y[ii] = float(m.group(3))
                    z[ii] = float(m.group(4))
                    el[ii] = ELMS[jj]
                    if jj < 0 or jj > 109:
                        ii -= 1
    else:
        ii = 1
    na = ii
    # Read atomic charges
    ii = 0
    if kch > 0 and jc < MAXBUF:
        kc = jc + na - 1
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
    # Output coordinates and atomic charges
    print(f"{na}\n{title}")
    if kch > 0:
        for i in range(1, na + 1):
            print(f"{el[i]:2s} {x[i]:12.6f} {y[i]:12.6f} {z[i]:12.6f} {c[i]:12.6f}")
    else:
        for i in range(1, na + 1):
            print(f"{el[i]:2s} {x[i]:12.6f} {y[i]:12.6f} {z[i]:12.6f}")
    # Output shieldings
    if shield:
        print("\n\nSHIELDINGS")
        for i in range(1, na + 1):
            if el[i] == "Bq":
                print(f"{el[i]:2s} " + " ".join(f"{s[jj][i]:11.4f}" for jj in range(6)))
            else:
                print(f"{el[i]:2s} " + " ".join(f"{s[jj][i]:11.4f}" for jj in range(2)))
    # Output normal modes
    if iv < MAXBUF:
        print(f"\n\nNORMALMODES {3*na-6}")
        ii = 1
        while ii <= na - 2 and iv < nlines:
            freq_line = lines[iv]
            m = re.match(r".*?([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)", freq_line[14:])
            if m:
                xa = float(m.group(1))
                xb = float(m.group(2))
                xc = float(m.group(3))
            else:
                xa = xb = xc = 0.0
            iv += 3
            # IR Intensity
            if lines[iv].startswith("IR Int"):
                m = re.match(r".*?([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)", lines[iv][14:])
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
                m = re.match(r".*?([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)", lines[iv][14:])
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
            # Normal mode vectors
            for i in range(1, na + 1):
                mode_line = lines[iv + i]
                m = re.findall(r"([\d\.\-]+)", mode_line[11:])
                if m and len(m) >= 9:
                    d = [float(val) for val in m[:9]]
                else:
                    d = [0.0] * 9
                print(f"{el[i]:2s} {d[0]:6.2f} {d[1]:6.2f} {d[2]:6.2f}")
            print(f"{xa:12.4f} {xia:12.4f} {xra:12.4f}")
            for i in range(1, na + 1):
                mode_line = lines[iv + i]
                m = re.findall(r"([\d\.\-]+)", mode_line[11:])
                if m and len(m) >= 9:
                    d = [float(val) for val in m[3:6]]
                else:
                    d = [0.0] * 3
                print(f"{el[i]:2s} {d[0]:6.2f} {d[1]:6.2f} {d[2]:6.2f}")
            print(f"{xb:12.4f} {xib:12.4f} {xrb:12.4f}")
            for i in range(1, na + 1):
                mode_line = lines[iv + i]
                m = re.findall(r"([\d\.\-]+)", mode_line[11:])
                if m and len(m) >= 9:
                    d = [float(val) for val in m[6:9]]
                else:
                    d = [0.0] * 3
                print(f"{el[i]:2s} {d[0]:6.2f} {d[1]:6.2f} {d[2]:6.2f}")
            print(f"{xc:12.4f} {xic:12.4f} {xrc:12.4f}")
            iv += na + 3
            ii += 1
    # Output molecular orbitals
    if io < MAXBUF:
        print("\n\nMOLORBS")
        # This is a simplified version, full parsing may require more logic
        for ij in range(1, 11):
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

    if len(sys.argv) < 2:
        print("Usage: python makexyz.py <gaussian_output_file>")
        return
    with open(sys.argv[1], "r") as f:
        lines = f.readlines()
    parse_gaussian_output(lines)


if __name__ == "__main__":
    main()
