import pandas as pd
from pathlib import Path

class XYZFile:
    def __init__(self, atom_data, timestep_data, file_path):
        self.atom_data = atom_data
        self.timestep_data = timestep_data
        self.file_path = file_path

        self._read_xyz()

    def _read_xyz(self):
        path = Path(self.file_path)
        with path.open('r') as f:
            lines = [line.strip() for line in f if line.strip()]

        if len(lines) < 2:
            raise ValueError(f"{self.file_path}: File too short, missing atom count or comment line.")

        try:
            atom_count = int(lines[0])
        except ValueError:
            raise ValueError(f"{self.file_path}: First line must be an integer (atom count). Got: {lines[0]}")

        comment = lines[1]
        atom_lines = lines[2:]

        print(len(atom_lines))
        print(atom_count)


        if len(atom_lines) != atom_count:
            print(f"{self.file_path}: Expected {atom_count} atom lines, but got {len(atom_lines)}")
            raise IndexError(f"{self.file_path}: Expected {atom_count} atom lines, but got {len(atom_lines)}")

        # Add entry to timestep_data
        timestep_index = len(self.timestep_data.dataframe)
        self.timestep_data.dataframe.loc[timestep_index] = {
            "file_name": path.name,
            "file_path": str(path),
            "timestep_time": None,
            "timestep_index": timestep_index,
            "raw_data": "\n".join(lines),
            "energy": None,
            "zero-point energy": None,
            "file_comment": comment,
            "measurements": {}
        }

        file_index = timestep_index  # assuming 1 file = 1 timestep

        for atom_index, line in enumerate(atom_lines[:atom_count]):
            tokens = line.split()
            if len(tokens) < 4:
                raise ValueError(
                    f"{self.file_path}: Malformed atom line {atom_index + 3} (expected ≥4 tokens): {line}"
                )

            element, x, y, z = tokens[:4]
            try:
                x, y, z = float(x), float(y), float(z)
            except ValueError:
                raise ValueError(
                    f"{self.file_path}: Coordinates must be numeric in line {atom_index + 3}: {line}"
                )

            # Optional: alias and charge
            alias = tokens[4] if len(tokens) > 4 else None
            charge = float(tokens[5]) if len(tokens) > 5 else None

            self.atom_data.dataframe.loc[len(self.atom_data.dataframe)] = {
                "file_name": path.name,
                "file_path": str(path),
                "file_index": file_index,
                "atom_index": atom_index,
                "element": element,
                "x": x,
                "y": y,
                "z": z,
                "alias": alias,
                "charge": charge,
                "timestep_time": None,
                "timestep_index": timestep_index
            }
