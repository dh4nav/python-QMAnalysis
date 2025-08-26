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
        with self.path.open('r') as f:
            lines = [line.rstrip('\n') for line in f]

        # Find the last orientation block (Standard/Input/Z-Matrix)
        orientation_start = None
        orientation_type = None
        for i, line in enumerate(lines):
            if re.search(r'(Standard|Input|Z-Matrix) orientation:', line):
                orientation_start = i
                orientation_type = line.strip()
        if orientation_start is None:
            raise ValueError(f"{self.file_path}: No orientation block found.")

        # Find the start and end of the atom table
        atom_table_start = orientation_start + 5
        atom_table_end = atom_table_start
        while atom_table_end < len(lines) and not lines[atom_table_end].startswith('-----'):
            atom_table_end += 1
        atom_lines = lines[atom_table_start:atom_table_end]

        atoms = []
        for line in atom_lines:
            tokens = line.split()
            # Skip header/separator lines and only parse lines with valid atom data
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

        # Archive block: find the last dashed line, then read until EOF
        archive_start = None
        for i in range(len(lines)-1, -1, -1):
            if re.match(r'-{5,}', lines[i]):
                archive_start = i + 1
                break
        archive_block = ''
        if archive_start is not None:
            archive_block = ''.join(lines[archive_start:]).replace('\n', '')
        # Extract values from archive block

        def extract_archive_value(key, block, is_tuple=False):
            m = re.search(rf'{key}=([^\\]+)', block)
            if not m:
                return pd.NA if not is_tuple else (pd.NA, pd.NA, pd.NA)
            val = m.group(1)
            if is_tuple:
                vals = val.split(',')
                return tuple(float(v) for v in vals[:3]) if len(vals) >= 3 else (pd.NA, pd.NA, pd.NA)
            try:
                return float(val)
            except Exception:
                return val
        zeropoint = extract_archive_value('ZeroPoint', archive_block)
        zpe = extract_archive_value('ZPE', archive_block)
        hf = extract_archive_value('HF', archive_block)
        thermal = extract_archive_value('Thermal', archive_block)
        rmsd = extract_archive_value('RMSD', archive_block)
        rmsf = extract_archive_value('RMSF', archive_block)
        dipole = extract_archive_value('Dipole', archive_block, is_tuple=True)
        nimag = extract_archive_value('NIMag', archive_block)

        # Save to frame_data.dataframe
        row = {
            "raw_data": '\n'.join(lines),
            "energy": hf,
            "zero-point energy": zeropoint,
            "file_comment": None,
            "ZPE": zpe,
            "Thermal": thermal,
            "RMSD": rmsd,
            "RMSF": rmsf,
            "dipole": dipole,
            "nimag": nimag
        }
        self.frame_data.dataframe.loc[(
            self.file_name, self.file_path, self.timestep_name)] = row
        # Ensure all columns are present
        expected_cols = ["raw_data", "energy", "zero-point energy",
                         "file_comment", "ZPE", "Thermal", "RMSD", "RMSF", "dipole", "nimag"]
        self.frame_data.dataframe = self.frame_data.dataframe.reindex(
            columns=expected_cols)

    @staticmethod
    def _atomic_number_to_symbol(num):
        # Minimal periodic table
        periodic_table = {
            1: 'H', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 19: 'K', 20: 'Ca',
            # Add more as needed
        }
        return periodic_table.get(num, f'El{num}')
