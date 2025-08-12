import argparse
import strictyaml as sy
import qmanalysis.xyzreader as xr
import qmanalysis.yamlreader as yr
import qmanalysis.measure as mr
from qmanalysis.containers import AtomData, TimestepData, MeasurementData, PerAtomData, GlobalData
import matplotlib.pyplot as plt
from qmanalysis.dataswitch import DfWrapper, DataSwitch
from pathlib import Path
import fnmatch
import re
import asteval as av
from qmanalysis.globalconstantsreader import GlobalConstantsFile


def prepend_root_if_relative(file_path, root_path=None):
    """
    Prepend `root_path` to `file_path` if `file_path` is not absolute.
    Returns a Path object.
    """
    file_path = Path(file_path)
    if file_path.is_absolute() or root_path is None:
        return file_path
    else:
        return Path(root_path) / file_path


def prepend_root_if_relative_and_glob(file_path, root_path=None):
    """
    Prepend `root_path` to `file_path` if `file_path` is not absolute.
    Returns a Path object.
    """
    file_path = Path(file_path)
    if file_path.is_absolute() or root_path is None:
        raise NotImplementedError(
            "Globbing not supported for absolute paths yet")
    else:
        pattern = str(file_path)  # ensure it's a string
        return [p.resolve() for p in Path(root_path).glob(pattern)]


def main():

    testyaml = """
  name: Procedure Name
  comment: XYZ
  version: 111
  files:
    - path: "../benzene.xyz"
      type: xyz
      name: benz
    - path: "../benzene.xyz"
      type: xyz
      name: gurk  
    - path: "../benzene.xyz"
      type: xyz
      name: abcd  
  substitutions:
    - name: S1
      entries:
      - file: benz
        atom_index: 5
      - file: benz2
        atom_index: 4
      - file: benz3
        atom_index: 1
  measurements:
    distance:
      - name: O-H-bond
        a: 6
        b: 7
      - name: Subs-bond
        a: S1
        b: 10
    angle:
      - name: H-O-H
        a: 3
        b: 6
        c: 7
      - name: H-O-N
        a: 8
        b: 9
        c: 10
    dihedral:
      - name: torsion1
        a: 2
        b: 3
        c: 6
        d: 7
  output:
    - file:
      - path: ../out.csv
        type: csv
    - graph:
      - type: scatter_plot
        x: distance-Subs-bond
        y: angle-H-O-N 
        file: ../SB-HON.tiff
  """

    # Example usage
    try:
        parser = argparse.ArgumentParser()
        parser.add_argument("inputfile", help="Main command input file")

        # # Add -ac / --atom-constants argument (multiple allowed)
        # parser.add_argument(
        #     "-ac", "--atom_constants",
        #     action="append",
        #     type=str,
        #     help="Path to file defining atom-level constants. Can be specified multiple times."
        # )

        # # Add -gc / --global-constants argument (multiple allowed)
        # parser.add_argument(
        #     "-gc", "--global_constants",
        #     action="append",
        #     type=str,
        #     help="Path to file defining global constants. Can be specified multiple times."
        # )

        # Add -rp / --root-path argument (single allowed — default behavior)
        parser.add_argument(
            "-rp", "--root_path",
            type=str,
            help="Root path used for all files. If none is specified, the location of the input file is used as root path"
        )

        args = parser.parse_args()
        print(args)

        if args.root_path is None:
            # print((Path.cwd() / args.inputfile).resolve().parent )
            args.root_path = (Path.cwd() / args.inputfile).resolve().parent

        print(args)

        yamlparser = yr.YAMLFile()
        # yamlparser.load_string(testyaml)
        # prepend_root_if_relative(file_path=args.inputfile, root_path=args.root_path))
        yamlparser.load_file(args.inputfile)
        # validate_measurements(data)
    except (sy.YAMLError, ValueError) as e:
        print(f"Validation error: {e}")

    atom_data = AtomData()
    timestep_data = TimestepData()

    print(yamlparser.get_data())

    yamldata = yamlparser.get_data()

    for file in yamldata["files"]:
        if file["type"].lower() == "xyz":
            if "glob" in file and file["glob"]:
                for gfile in list(prepend_root_if_relative_and_glob(file_path=file["path"], root_path=args.root_path)):
                    xr.XYZFile(atom_data, timestep_data, file_path=gfile,
                               file_name=file["name"]+Path(gfile).stem)
            else:
                xr.XYZFile(atom_data, timestep_data, file_path=prepend_root_if_relative(
                    file_path=file["path"], root_path=args.root_path), file_name=file["name"])
        elif file["type"].lower() == "global_constants_csv":
            global_constants = GlobalConstantsFile(file_path=file["path"])
    # print(atom_data.dataframe)
    # print(timestep_data.dataframe)

    if "substitutions" in yamldata:
        for sub in yamldata["substitutions"]:
            for subfile in sub['entries']:
                # mask = (atom_data.dataframe["file_name"] == subfile['file']) & (atom_data.dataframe["atom_index"] == subfile['atom_index'])
                # atom_data.dataframe[mask, ["alias"]] == 777
                if "glob" in subfile:
                    # print(subfile['file'])
                    regex_pattern = fnmatch.translate(subfile['file'])
                    regex = re.compile(regex_pattern)
                    # print(regex)
                    atom_data.dataframe.loc[atom_data.dataframe["file_name"].str.match(regex) & (
                        atom_data.dataframe["atom_index"] == subfile['atom_index']), "alias"] = sub["name"]
                else:
                    atom_data.dataframe.loc[(atom_data.dataframe["file_name"] == subfile['file']) & (
                        atom_data.dataframe["atom_index"] == subfile['atom_index']), "alias"] = sub["name"]
    print("Finished reading")
    print(atom_data.dataframe)
    print(timestep_data.dataframe)

    timestep_names = timestep_data.dataframe.loc[:, "file_name"].to_numpy()

    # print(timestep_names)

    # from asteval import Interpreter
    # from types import MethodType

    # class Symboler:
    #   def __init__(self):
    #     self._store = {}

    #   def __getitem__(self, key):
    #     return self._store[key]

    #   def __setitem__(self, key, value):
    #     self._store[key] = value

    #   def __getattr__(self, name):
    #     # Called for undefined attributes
    #     return SymbolerAccessor(self, [name])

    #   def __call__(self, name, *args, **kwargs):
    #     print(f"Called {name} with args={args}, kwargs={kwargs}")
    #     return f"<result of {name}()>"

    # class SymbolerAccessor:
    #   def __init__(self, symboler, path):
    #     self._symboler = symboler
    #     self._path = path

    #   def __getattr__(self, name):
    #     return SymbolerAccessor(self._symboler, self._path + [name])

    #   def __call__(self, *args, **kwargs):
    #     full_name = ".".join(self._path)
    #     return self._symboler(full_name, *args, **kwargs)

    # class SymbolerInterpreter(Interpreter):
    #   def __init__(self, symboler, **kwargs):
    #     self.symboler = symboler
    #     super().__init__(usersyms={}, **kwargs)

    #   def run(self, expr, **kwargs):
    #     # sync all symboler keys to usersyms
    #     for k, v in self.symboler._store.items():
    #       self.symtable[k] = v
    #     return super().run(expr, **kwargs)
    # from collections import UserDict
    # from types import MethodType

    # * plain variables
    # * const is a dictionary, const.name returns constant name
    # * atoms, measurements, frames are all pd.DataFrame objects
    #   * attribute access in these (x as placeholder):
    #     * x.attribute returns one-column pd inside DfWrapper with respective column
    #     * x() returns filter result as pd in DfWrapper. Parameters: ts= for timestep, file= for file name. file name can contain wildcards or a list

    # class DfWrapper():
    #     def __init__(self, df):
    #         self._df = df

    #     def __getattr__(self, name):
    #         if name in self._df:
    #             print("getx+ "+name)
    #             return self._df[name]
    #         else:
    #             print("getx "+name)
    #             raise AttributeError("Missing attribute "+name)

    #     def _wildcard_to_regex(self, pattern):
    #         # Escape all regex special chars first
    #         escaped = re.escape(pattern)
    #         # Replace escaped wildcard chars by regex equivalents
    #         regex = escaped.replace(r'\*', '.*').replace(r'\?', '.')
    #         # Add anchors to match the whole string
    #         regex = '^' + regex + '$'
    #         return regex

    #     def _filter_rows(self, df, name, value):
    #         if isinstance(value, (int, float)):
    #             return df[df[name] == value]
    #         if isinstance(value, list):
    #             return df[df[name].isin(value)]
    #         if isinstance(value, str):
    #             rex = self._wildcard_to_regex(value)
    #             return df[df[name].str.match(rex)]

    #     def __call__(self, *args, **kwds):
    #         temp_df = self._df
    #         print(kwds)
    #         if "file_name" in kwds and kwds['file_name'] is not None:
    #             print("requested file_name " + str(temp_df))
    #             temp_df = self._filter_rows(
    #                 temp_df, "file_name", kwds['file_name'])
    #             print("requested file_name " + str(temp_df))

    #         if "file_index" in kwds and kwds['file_index'] is not None:
    #             temp_df = self._filter_rows(
    #                 temp_df, "file_index", kwds['file_index'])
    #             print("requested file_index " + str(temp_df))

    #         if "ts_name" in kwds and kwds['ts_name'] is not None:
    #             temp_df = self._filter_rows(
    #                 temp_df, "timestep_name", kwds['ts_name'])
    #             print("requested ts_name " + str(temp_df))

    #         if "ts_index" in kwds and kwds['ts_index'] is not None:
    #             temp_df = self._filter_rows(
    #                 temp_df, "timestep_index", kwds['ts_index'])
    #             print("requested ts_index " + str(temp_df))

    #         return DfWrapper(temp_df)

    # class Symboler(UserDict):
    #   def __init__(self, root_obj=None):
    #     self.root_obj = root_obj

    #   def __getitem__(self, key):
    #     if self.root_obj:
    #       try:
    #         return self.root_obj[key]
    #     else:
    #       try:
    #         if isinstance(self.data[key], pd.DataFrame):
    #           return Symboler(root_obj=self.data[key])
    #         else:
    #           return self.data[key]

    #   def __setitem__(self, key, value):
    #     if self.root_obj:
    #       try:
    #         self.root_obj[key] = value
    #     else:
    #       try:
    #         self.data[key]

    #   def __getattribute__(self, name):
    #     if self.root_obj:
    #       try:
    #         if hasattr(self.root_obj, name) and type(getattr(self.root_obj, name)) == types.MethodType:
    #           return
    #         return self.root_obj[name]
    #     else:
    #       try:
    #         if isinstance(self.data[name], pd.DataFrame):
    #           return Symboler(root_obj=self.data[name])
    #         else:
    #           return self.data[name]

    # # Example usage
    # symb = Symboler()
    # symb['x'] = 42
    # symb['y'] = 3.14
    # symb['atoms'] = atom_data.dataframe

    # test_df = DfWrapper(atom_data.dataframe)

    # Create an Interpreter with the custom symbol table
    aeval = DataSwitch(userdicts={
                       "const": global_constants.global_constants_data, "atoms": atom_data.dataframe})

    # Evaluate an expression using the custom symbol table
    result = aeval.query(
        """atoms(file_name='gurke').x / const.pi""")
    print(str(result))  # Output: 45.14

    measure = mr.Measure()

    md = MeasurementData(timestep_names)
    if args.atom_constants:
        ac = PerAtomData([prepend_root_if_relative(
            file_path=f, root_path=args.root_path) for f in args.atom_constants])
    else:
        ac = None

    if args.global_constants:
        gc = GlobalData([prepend_root_if_relative(
            file_path=f, root_path=args.root_path) for f in args.global_constants])
    else:
        gc = None

    data_switch = DataSwitch(
        dataframe=md, atom_constants=ac, global_constants=gc)
    print(timestep_names)

    temp_array = []

    if "measurements" in yamldata:
        if "distance" in yamldata["measurements"]:
            for distance in yamldata["measurements"]["distance"]:
                # print(distance['a'])
                temp_array = []
                for timestep_name in timestep_names:
                    # print(timestep_name)
                    # print(atom_data.dataframe[(atom_data.dataframe['file_name'] == timestep_name) & (atom_data.dataframe['alias'] == str(distance['a'])  )])
                    temp_array.append(measure.distance(
                        atom_data=atom_data,
                        atom_index1=atom_data.dataframe[(atom_data.dataframe['file_name'] == timestep_name) & (
                            atom_data.dataframe['alias'] == str(distance['a']))].index[0],
                        atom_index2=atom_data.dataframe[(atom_data.dataframe['file_name'] == timestep_name) & (atom_data.dataframe['alias'] == str(distance['b']))].index[0]))
                md.dataframe[distance["name"]] = temp_array
                # print(md.dataframe)

        if "angle" in yamldata["measurements"]:
            for angle in yamldata["measurements"]["angle"]:
                # print(angle)
                temp_array = []
                for timestep_name in timestep_names:
                    temp_array.append(measure.angle(
                        atom_data=atom_data,
                        atom_index1=atom_data.dataframe[(atom_data.dataframe['file_name'] == timestep_name) & (
                            atom_data.dataframe['alias'] == str(angle['a']))].index[0],
                        atom_index2=atom_data.dataframe[(atom_data.dataframe['file_name'] == timestep_name) & (
                            atom_data.dataframe['alias'] == str(angle['b']))].index[0],
                        atom_index3=atom_data.dataframe[(atom_data.dataframe['file_name'] == timestep_name) & (atom_data.dataframe['alias'] == str(angle['c']))].index[0]))
                md.dataframe[angle["name"]] = temp_array
        if "dihedral" in yamldata["measurements"]:
            for dihedral in yamldata["measurements"]["dihedral"]:
                # print(dihedral)
                temp_array = []
                for timestep_name in timestep_names:
                    temp_array.append(measure.dihedral(
                        atom_data=atom_data,
                        atom_index1=atom_data.dataframe[(atom_data.dataframe['file_name'] == timestep_name) & (
                            atom_data.dataframe['alias'] == str(dihedral['a']))].index[0],
                        atom_index2=atom_data.dataframe[(atom_data.dataframe['file_name'] == timestep_name) & (
                            atom_data.dataframe['alias'] == str(dihedral['b']))].index[0],
                        atom_index3=atom_data.dataframe[(atom_data.dataframe['file_name'] == timestep_name) & (
                            atom_data.dataframe['alias'] == str(dihedral['c']))].index[0],
                        atom_index4=atom_data.dataframe[(atom_data.dataframe['file_name'] == timestep_name) & (atom_data.dataframe['alias'] == str(dihedral['d']))].index[0]))
                md.dataframe[dihedral["name"]] = temp_array
    # print(md.dataframe)
    # print(atom_data.dataframe)
    # print(timestep_data.dataframe)
    # print(yamldata['output'])
    for one_output in yamldata['output']:
        if 'file' in one_output:
            # print(one_output)
            for file in one_output['file']:
                # print(file)
                if file['type'].lower() == "csv":
                    md.dataframe.to_csv(prepend_root_if_relative(
                        file_path=file['path'], root_path=args.root_path))
                else:
                    raise IndexError(f"{file['type']}: Unknown file type")
         # if "graphics" in yamldata['output']:
        # print("Z0")
        if 'graph' in one_output:
            # print("Z1")
            for graph in one_output['graph']:
                if graph['type'].lower() == "scatter_plot":
                    # print("Z")
                    # alt:
                    # fig, ax = plt.subplots(figsize=graph.get("figsize", (8, 6)))  # default figsize
                    # df.plot(x='x', y='y', ax=ax)

                    plot = md.dataframe.plot.scatter(
                        x=graph['x'], y=graph['y'])
                    # plt.show()
                    fig = plot.get_figure()
                    fig.savefig(prepend_root_if_relative(file_path=graph['file'], root_path=args.root_path), dpi=graph.get(
                        "dpi", 300), format=graph.get("file_format", "tiff"))

    # print(data_switch("data"))


if __name__ == "__main__":
    main()
