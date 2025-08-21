import pandas as pd
import numpy as np
import asteval as av
import re


class DfWrapper():
    def __init__(self, df):
        self._df = df

    def __getattr__(self, name):
        print("ga"+name)
        # if name in self._df:
        #     print("getx+ "+name)
        #     return self._df[name]
        # else:
        #     print("getx "+name)
        #     raise AttributeError("Missing attribute "+name)
        if hasattr(self._df, name):
            return getattr(self._df, name)
        raise AttributeError(f"{type(self).__name__} has no attribute {name}")

    def _wildcard_to_regex(self, pattern):
        # Escape all regex special chars first
        escaped = re.escape(pattern)
        # Replace escaped wildcard chars by regex equivalents
        regex = escaped.replace(r'\*', '.*').replace(r'\?', '.')
        # Add anchors to match the whole string
        regex = '^' + regex + '$'
        return regex

    def _filter_rows(self, df, name, value):
        print("fr")
        if isinstance(value, (int, float)):
            return self._get_single_scalar(df[df[name] == value])
        if isinstance(value, list):
            return self._get_single_scalar(df[df[name].isin(value)])
        if isinstance(value, str):
            return self._get_single_scalar(df[df[name].str.match(self._wildcard_to_regex(value))])

    def _get_single_scalar(self, df):
        print("df0")
        if df.shape == (1, 1) or df.shape == (1,):
            print("df1")
            print(df)
            return df[0, 0]
        else:
            return df

    def __call__(self, *args, **kwds):
        temp_df = self._df
        print("call"+str(args)+str(kwds))
        if "file_name" in kwds and kwds['file_name'] is not None:
            # print("requested file_name " + str(temp_df))
            temp_df = self._filter_rows(
                temp_df, "file_name", kwds['file_name'])
            # print("requested file_name " + str(temp_df))

        if "file_index" in kwds and kwds['file_index'] is not None:
            temp_df = self._filter_rows(
                temp_df, "file_index", kwds['file_index'])
            # print("requested file_index " + str(temp_df))

        if "ts_name" in kwds and kwds['ts_name'] is not None:
            temp_df = self._filter_rows(
                temp_df, "timestep_name", kwds['ts_name'])
            # print("requested ts_name " + str(temp_df))

        if "ts_index" in kwds and kwds['ts_index'] is not None:
            temp_df = self._filter_rows(
                temp_df, "timestep_index", kwds['ts_index'])
            # print("requested ts_index " + str(temp_df))

        return DfWrapper(temp_df)


class DataSwitch():
    def __init__(self, userdicts: dict = None):
        dfws = dict()
        for k in userdicts.keys():
            dfws[k] = DfWrapper(userdicts[k])
        self._interpreter = av.Interpreter(user_symbols=dfws)

    def __getattr__(self, name):
        pass

    def query(self, query):
        return self._interpreter(query)


# class DataMangler():
#     def __init__(self, userdicts: dict = None):
#         dfws = dict()
#         for k in userdicts.keys():
#             dfws[k] = userdicts[k]

#     def interpret(self, commands: list):
#         heap =
