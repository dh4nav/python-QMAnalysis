import pandas as pd
from pathlib import Path


class GlobalConstantsFile:
    def __init__(self, file_path):
        self.df = None
        self.file_path = file_path
        self.base_path = path = Path(self.file_path)
        self.file_name = self.base_path.name

        self._read_csv()

    def _read_csv(self):
        read_df = pd.read_csv(self.file_path)
        self.df = read_df.set_index(
            "name").T.reset_index(drop=True)
        self.df.columns.name = None
        # pivot(
        #  columns="name", values="value")

        print(read_df)
        print(self.df)
