import pandas as pd
import numpy as np
import asteval as av

class DotNamespace:
    def __init__(self, data):
        self._data = data

    def __getattr__(self, name):
        if isinstance(self._data, pd.DataFrame):
            if name in self._data.columns:
                return self._data[name]
        elif isinstance(self._data, dict):
            if name in self._data:
                value = self._data[name]
                # Nest deeper if value is dict or dataframe
                if isinstance(value, (dict, pd.DataFrame)):
                    return DotNamespace(value)
                return value
        raise AttributeError(f"'{name}' not found in DotNamespace.")

class DataSwitch(av.Interpreter):
    def __init__(self, dataframe, atom_constants=None, global_constants=None):
        # Wrap each dataframe (or nested dicts of dfs) in DotNamespace
        d = {"data": dataframe}
        if atom_constants:
            d["atom_constants"] = atom_constants
        if global_constants:
            d["global_constants"] = global_constants#

        symtable = {k: DotNamespace(v) for k, v in d.items()}
        super().__init__(symtable=symtable)
