import pandas as pd
from pathlib import Path


# class AtomData:
#     def __init__(self):
#         idx = pd.MultiIndex.from_tuples([], names=(
#             'file_name', 'file_path', 'timestep_name', 'atom_index'))
#         self.dataframe = pd.DataFrame(
#             index=ind, columns=["element", "alias", "charge", "x", "y", "z"])


# class FrameData:
#     def __init__(self):
#         idx = pd.MultiIndex.from_tuples([], names=(
#             'file_name', 'file_path', 'timestep_name'))
#         self.dataframe = pd.DataFrame(
#             columns=["raw_data", "energy", "zero-point energy", "file_comment"])

class XYZFile:
    def __init__(self, atom_data, frame_data, file_path, file_name=None, timestep_name=None):
        self.atom_data = atom_data
        self.timestep_data = frame_data
        self.file_path = file_path
        self.path = Path(self.file_path)
        if file_name:
            self.file_name = file_name
        else:
            self.file_name = self.path.name
        self.timestep_name = timestep_name

        self._read_xyz()

    def _read_xyz(self):
        with self.path.open('r') as f:
            lines = [line.strip() for line in f if line.strip()]

        if len(lines) < 2:
            raise ValueError(
                f"{self.file_path}: File too short, missing atom count or comment line.")

        try:
            atom_count = int(lines[0])
        except ValueError:
            raise ValueError(
                f"{self.file_path}: First line must be an integer (atom count). Got: {lines[0]}")

        comment = lines[1]
        atom_lines = lines[2:]

        print(len(atom_lines))
        print(atom_count)

        if len(atom_lines) < atom_count:
            print(
                f"{self.file_path}: Expected {atom_count} atom lines, but got {len(atom_lines)}")
            raise IndexError(
                f"{self.file_path}: Expected {atom_count} atom lines, but got {len(atom_lines)}")
        # Add entry to timestep_data

        self.timestep_data.dataframe.loc[(self.file_name, self.file_path, self.timestep_name)] = {
            "raw_data": "\n".join(lines),
            "energy": None,
            "zero-point energy": None,
            "file_comment": comment
        }

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
            alias = tokens[4] if len(tokens) > 4 else str(atom_index)
            charge = float(tokens[5]) if len(tokens) > 5 else None
            self.atom_data.dataframe.loc[(self.file_name, self.file_path, self.timestep_name, atom_index)] = {
                "element": element,
                "x": x,
                "y": y,
                "z": z,
                "alias": alias,
                "charge": charge
            }
