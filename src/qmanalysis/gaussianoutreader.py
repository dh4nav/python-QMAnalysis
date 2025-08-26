import pandas as pd
import re
from pathlib import Path


class GaussianOutFile:
    def __init__(self, atom_data, frame_data, file_path, file_name=None, timestep_name=None):
        self.atom_data = atom_data
        self.frame_data = frame_data
        self.file_path = file_path
        self.path = Path(self.file_path)
        self.file_name = file_name if file_name else self.path.name
        self.timestep_name = timestep_name
        self._read_gaussian_out()

    def _read_gaussian_out(self):
        orientation_start = None
        orientation_type = None
        atom_lines = []
        archive_lines = []
        raw_lines = []
        last_dash_idx = None
        with self.path.open('r') as f:
            for idx, line in enumerate(f):
                line = line.rstrip('\n')
                raw_lines.append(line)
                if re.search(r'(Standard|Input|Z-Matrix) orientation:', line):
                    orientation_start = idx
                    orientation_type = line.strip()
                if re.match(r'-{5,}', line):
                    last_dash_idx = idx
        # Atom block: only read if orientation_start found
        if orientation_start is None:
            raise ValueError(f"{self.file_path}: No orientation block found.")
        # Read atom table lines directly from file
        with self.path.open('r') as f:
            for idx, line in enumerate(f):
                if idx < orientation_start + 5:
                    continue
                if last_dash_idx is not None and idx >= last_dash_idx:
                    break
                if idx >= orientation_start + 5:
                    if line.startswith('-----'):
                        break
                    atom_lines.append(line.rstrip('\n'))
        atoms = []
        for line in atom_lines:
            tokens = line.split()
            if len(tokens) < 6:
                continue
            try:
                atomic_num = int(tokens[1])
            except (ValueError, IndexError):
                continue
            try:
                x, y, z = float(tokens[3]), float(tokens[4]), float(tokens[5])
            except (ValueError, IndexError):
                continue
            element = self._atomic_number_to_symbol(atomic_num)
            atom_index = len(atoms)
            alias = str(atom_index)
            charge = None
            atoms.append((element, x, y, z, alias, charge))
            self.atom_data.dataframe.loc[(self.file_name, self.file_path, self.timestep_name, atom_index)] = {
                "element": element,
                "x": x,
                "y": y,
                "z": z,
                "alias": alias,
                "charge": charge
            }
        # Archive block: read only the lines after last dashed line
        # Archive block starts with line beginning '1\'
        archive_block = ''
        with self.path.open('r') as f:
            lines = [line.rstrip('\n') for line in f]
        archive_start = None
        for i, line in enumerate(lines):
            if re.match(r'^\s*1\\1\\FAU-CCC', line):
                archive_start = i
                break
        if archive_start is not None:
            # Strip leading/trailing spaces from each line before joining
            archive_lines = [line.strip() for line in lines[archive_start:]]
            archive_block = '\\'.join(archive_lines)
            archive_block = archive_block.replace('\n', '').replace('\r', '')

        def is_plausible(val, key=None):
            # Check for plausible float values
            if val is pd.NA or val is None:
                return False
            try:
                fval = float(val)
                if not (-1e10 < fval < 1e10):
                    return False
                if str(val).lower() in ['nan', 'inf', '-inf', '']:
                    return False
                if re.match(r'.*e-$', str(val)):
                    return False
                return True
            except Exception:
                return False

        # Robust extraction of comment and charge/multiplicity
        file_comment = None
        charge = pd.NA
        multiplicity = pd.NA
        split_block = [s.strip()
                       for s in archive_block.split('\\') if s.strip()]
        # List of element symbols for detection
        element_symbols = {'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh',
                           'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U'}
        atom_block_start = None
        for i, field in enumerate(split_block):
            # Detect atom block start
            if ',' in field:
                first_token = field.split(',')[0]
                if first_token in element_symbols:
                    atom_block_start = i
                    break
        # Comment is the field before atom block
        if atom_block_start and atom_block_start > 0:
            file_comment = split_block[atom_block_start-1]
        # Charge/multiplicity is the field before the comment
        if atom_block_start and atom_block_start > 1:
            charge_mult_field = split_block[atom_block_start-2]
            parts = [p.strip() for p in charge_mult_field.split(',')]
            if len(parts) == 2:
                try:
                    charge = int(parts[0]) if parts[0].isdigit(
                    ) else float(parts[0])
                except Exception:
                    charge = pd.NA
                try:
                    multiplicity = int(parts[1])
                except Exception:
                    multiplicity = pd.NA

        def extract_archive_value(key, block, is_tuple=False):
            m = re.search(rf'{key}=([^\\]*)', block, re.IGNORECASE)
            if not m:
                return pd.NA if not is_tuple else (pd.NA, pd.NA, pd.NA)
            val = m.group(1).replace('\n', '').replace(' ', '').strip()
            if is_tuple:
                vals = [v.strip() for v in val.split(',')]
                tup = tuple(float(v) if is_plausible(v) else pd.NA for v in vals[:3]) if len(
                    vals) >= 3 else (pd.NA, pd.NA, pd.NA)
                return tup
            try:
                fval = float(val)
                return fval if is_plausible(fval, key) else pd.NA
            except Exception:
                return pd.NA
        # DEBUG: Print the archive block (commented out)
        # print("\n--- Archive Block ---\n", archive_block,
        #       "\n--- End Archive Block ---\n")

        # DEBUG: Print extracted values
        def debug_print(key, value):
            print(f"Extracted {key}: {value}")

        debug_print('Comment', file_comment)
        debug_print('Charge', charge)
        debug_print('Multiplicity', multiplicity)

        zeropoint = extract_archive_value('ZeroPoint', archive_block)
        debug_print('ZeroPoint', zeropoint)
        zpe = extract_archive_value('ZPE', archive_block)
        debug_print('ZPE', zpe)
        hf = extract_archive_value('HF', archive_block)
        debug_print('HF', hf)
        thermal = extract_archive_value('Thermal', archive_block)
        debug_print('Thermal', thermal)
        rmsd = extract_archive_value('RMSD', archive_block)
        debug_print('RMSD', rmsd)
        rmsf = extract_archive_value('RMSF', archive_block)
        debug_print('RMSF', rmsf)
        dipole = extract_archive_value('Dipole', archive_block, is_tuple=True)
        debug_print('Dipole', dipole)
        nimag = extract_archive_value('NIMag', archive_block)
        debug_print('NIMag', nimag)
        # Save to frame_data.dataframe
        row = {
            "raw_data": '\n'.join(raw_lines),
            "energy": hf,
            "zero-point energy": zeropoint,
            "file_comment": file_comment,
            "charge": charge,
            "multiplicity": multiplicity,
            "ZPE": zpe,
            "Thermal": thermal,
            "RMSD": rmsd,
            "RMSF": rmsf,
            "dipole": dipole,
            "nimag": nimag
        }
        self.frame_data.dataframe.loc[(
            self.file_name, self.file_path, self.timestep_name)] = row
        expected_cols = ["raw_data", "energy", "zero-point energy",
                         "file_comment", "charge", "multiplicity", "ZPE", "Thermal", "RMSD", "RMSF", "dipole", "nimag"]
        self.frame_data.dataframe = self.frame_data.dataframe.reindex(
            columns=expected_cols)

    @staticmethod
    def _atomic_number_to_symbol(num):
        # Complete periodic table up to uranium (Z=92)
        periodic_table = {
            1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne',
            11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca',
            21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn',
            31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr',
            41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn',
            51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce', 59: 'Pr', 60: 'Nd',
            61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb',
            71: 'Lu', 72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au', 80: 'Hg',
            81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th',
            91: 'Pa', 92: 'U'
        }
        return periodic_table.get(num, f'El{num}')
