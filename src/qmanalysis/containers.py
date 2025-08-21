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


# class PerAtomData:
#     def __init__(self, csv):
#         self.dataframe = pd.read_csv(csv[0])
#         for f in csv[1:]:
#             temp_dataframe = pd.read_csv(csv)

#             # Confirm row count matches before concatenating
#             if len(self.dataframe) != len(temp_dataframe):
#                 raise ValueError(
#                     "Row count does not match between the existing DataFrame and the CSV file "+f)

#             # Concatenate along columns (axis=1)
#             self.dataframe = pd.concat(
#                 [self.dataframe, temp_dataframe], axis=1)


# class GlobalData:
#     def __init__(self, csv):
#         self.dataframe = pd.read_csv(csv[0])
#         for f in csv[1:]:
#             temp_dataframe = pd.read_csv(csv)

#             # Confirm row count matches before concatenating
#             if len(self.dataframe) != len(temp_dataframe):
#                 raise ValueError(
#                     "Row count does not match between the existing DataFrame and the CSV file "+f)

#             # Concatenate along columns (axis=1)
#             self.dataframe = pd.concat(
#                 [self.dataframe, temp_dataframe], axis=1)
