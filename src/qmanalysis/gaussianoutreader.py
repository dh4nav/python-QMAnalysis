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
        # Read the file and store only the last orientation and archive block
        last_orientation_lines = []
        last_archive_lines = []
        in_orientation = False
        in_archive = False
        with self.path.open('r') as f:
            for line in f:
                line = line.rstrip('\n')
                # Detect orientation block start
                if re.search(r'(Standard|Input|Z-Matrix) orientation:', line):
                    in_orientation = True
                    last_orientation_lines = []
                if in_orientation:
                    last_orientation_lines.append(line)
                    if re.match(r'-{5,}', line):
                        in_orientation = False
                # Detect archive block start
                if re.match(r'^\s{1}1\\1\\', line):
                    in_archive = True
                    last_archive_lines = []
                if in_archive:
                    last_archive_lines.append(line)
                    # End archive block at empty line or line starting with '@'
                    if line.strip() == '' or line.lstrip().startswith('@'):
                        in_archive = False
        # Now process only the last orientation and archive block
        atoms = []
        for line in last_orientation_lines[5:]:
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
        archive_start = None
        archive_end = None
        for i, line in enumerate(last_archive_lines):
            if re.match(r'^\s{1}1\\1\\', line):
                archive_start = i
            # if archive_start is not None:
            #    print(f"End check line {i}: {repr(line)}")  # DEBUG
            if archive_start is not None and line.strip() == '':
                archive_end = i
                break
        if archive_start is not None and archive_end is not None:
            # Remove line breaks and one space at start/end of each line
            archive_lines = [line[1:-1] if line.startswith(' ') and line.endswith(
                ' ') else line.lstrip(' ').rstrip(' ') for line in last_archive_lines[archive_start:archive_end+1]]
            archive_block = ''.join(archive_lines)
            # Split at each '\', keeping empty fields
            split_block = archive_block.split('\\')
        else:
            split_block = []
        print("Archive split fields:")
        for idx, field in enumerate(split_block):
            print(f"[{idx}]: {field}")

        element_symbols = {'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh',
                           'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U'}

        def is_atom_line(field):
            tokens = field.split(',')
            if len(tokens) < 4:
                return False
            if tokens[0] not in element_symbols:
                return False
            try:
                float(tokens[1])
                float(tokens[2])
                float(tokens[3])
                return True
            except Exception:
                return False
        atom_block_start = None
        for i, field in enumerate(split_block):
            if is_atom_line(field):
                atom_block_start = i
                break

        # Process archive block fields
        archive_fields = []
        for line in last_archive_lines:
            archive_fields.extend(line.split('\\'))
        # Extract comment and charge/multiplicity
        file_comment = archive_fields[15] if len(archive_fields) > 15 else None
        if len(archive_fields) > 17:
            charge_multiplicity = archive_fields[17]
            print("CM: "+charge_multiplicity)
            charge, multiplicity = charge_multiplicity.split(
                ',') if charge_multiplicity else (pd.NA, pd.NA)

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
        nimag = extract_archive_value('NImag', archive_block)
        # Store NImag as integer if possible
        try:
            nimag = int(
                nimag) if nimag is not pd.NA and nimag is not None else pd.NA
        except Exception:
            nimag = pd.NA
        debug_print('NIMag', nimag)
        # Save to frame_data.dataframe
        row = {
            "raw_data": '\n'.join(raw_lines),
            "energy": pd.NA,
            "zero-point energy": zeropoint,
            "file_comment": file_comment,
            "charge": charge,
            "multiplicity": multiplicity,
            "HF": hf,
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
                         "file_comment", "charge", "multiplicity", "HF", "ZPE", "Thermal", "RMSD", "RMSF", "dipole", "nimag"]
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
