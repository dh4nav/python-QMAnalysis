from collections import OrderedDict
from pathlib import PurePath
import pandas as pd

class XYZFile():
    def __init__(self, atom_data, timestep_data, file_path, file_name = None):

        with open(file_path, "r", encoding="utf-8") as f:
            raw_data = f.read()

        if file_name is None:
            file_name = PurePath(file_path).name

        lines = raw_data.strip().splitlines()

        if len(lines) != int(lines[0])+2:
            raise IndexError("Number of Lines does not match number of atoms in line two")

        file_index = len(timestep_data.dataframe)

        timestep_data.dataframe.loc[file_index] = {
            "file_name": file_name, 
            "file_path": file_path, 
            "timestep_time": pd.NA, 
            "timestep_index": pd.NA, 
            "raw_data": raw_data, 
            "energy": pd.NA, 
            "zero-point energy": pd.NA,
            "file_comment": lines[1], 
            "measurements": {}
            }
        
        first_atom_index = len(atom_data.dataframe)

        for atom_index in range(int(lines[0])):
            elements = lines[atom_index+2].strip().split()
            if len(elements) != 4:
                raise ValueError(f"Line {atom_index+2} must have 4 fields: element x y z")
            atom_data.dataframe.loc[first_atom_index+atom_index] = {
            "file_name": file_name, 
            "file_path": file_path,
            "file_index": file_index,
            "atom_index": (atom_index),
            "element": elements[0],
            "x": float(elements[1]),
            "y": float(elements[2]),
            "z": float(elements[3]),
            "alias": str(atom_index),
            "charge": pd.NA,
            "timestep_time": pd.NA, 
            "timestep_index": pd.NA,     
            }




