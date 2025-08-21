import pandas as pd


class AtomData:
    def __init__(self):
        idx = pd.MultiIndex.from_tuples([], names=(
            'file_name', 'file_path', 'timestep_name', 'atom_index'))
        self.dataframe = pd.DataFrame(
            index=idx, columns=["element", "alias", "charge", "x", "y", "z"])


class FrameData:
    def __init__(self):
        idx = pd.MultiIndex.from_tuples([], names=(
            'file_name', 'file_path', 'timestep_name'))
        self.dataframe = pd.DataFrame(index=idx,
                                      columns=["raw_data", "energy", "zero-point energy", "file_comment"])


class MeasurementData:
    def __init__(self):
        idx = pd.MultiIndex.from_tuples([], names=(
            'file_name', 'file_path', 'timestep_name'))
        self.dataframe = pd.DataFrame(index=idx, columns=[])
