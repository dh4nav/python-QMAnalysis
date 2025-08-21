import argparse
import strictyaml as sy
import qmanalysis.xyzreader as xr
import qmanalysis.yamlreader as yr
import qmanalysis.measure as mr
from qmanalysis.containers import AtomData, FrameData, MeasurementData
import matplotlib.pyplot as plt
from qmanalysis.dataswitch import DfWrapper, DataSwitch
from pathlib import Path
import fnmatch
import re
import asteval as av
from qmanalysis.globalconstantsreader import GlobalConstantsFile
import pandas as pd
import numpy as np
import scipy

import qmanalysis.customcalculationrunner as ccr


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

        # Add -rp / --root-path argument (single allowed — default behavior)
        parser.add_argument(
            "-rp", "--root_path",
            type=str,
            help="Root path used for all files. If none is specified, the location of the input file is used as root path"
        )

        args = parser.parse_args()

        if args.root_path is None:
            # print((Path.cwd() / args.inputfile).resolve().parent )
            args.root_path = (Path.cwd() / args.inputfile).resolve().parent

        yamlparser = yr.YAMLFile()
        # yamlparser.load_string(testyaml)
        # prepend_root_if_relative(file_path=args.inputfile, root_path=args.root_path))
        yamlparser.load_file(args.inputfile)
        # validate_measurements(data)
    except (sy.YAMLError, ValueError) as e:
        print(f"Validation error: {e}")

    atom_data = AtomData()
    frame_data = FrameData()
    mesurement_data = MeasurementData()

    print(yamlparser.get_data())

    yamldata = yamlparser.get_data()

    for file in yamldata["files"]:
        if file["type"].lower() == "xyz":
            if "glob" in file and file["glob"]:
                for gfile in list(prepend_root_if_relative_and_glob(file_path=file["path"], root_path=args.root_path)):
                    xr.XYZFile(atom_data, frame_data, file_path=gfile,
                               file_name=file["name"]+Path(gfile).stem, timestep_name=file.get("timestep", None))
            else:
                xr.XYZFile(atom_data, frame_data, file_path=prepend_root_if_relative(
                    file_path=file["path"], root_path=args.root_path), file_name=file["name"], timestep_name=file.get("timestep", None))
        elif file["type"].lower() == "global_constants_csv":
            global_constants = GlobalConstantsFile(file_path=file["path"])
        elif file["type"].lower() == "per_file_constants_csv":
            per_file_constants = pd.read_csv(file_path=file["path"])
            frame_data = pd.merge(
                frame_data.df, per_file_constants, how='outer', on='file_name')
    # print(atom_data.dataframe)
    # print(timestep_data.dataframe)

    print("Atom data:")
    print(atom_data.dataframe)
    print("Timestep data:")
    print(frame_data.dataframe)

    if "substitutions" in yamldata:
        for sub in yamldata["substitutions"]:
            for subfile in sub['entries']:
                # Handle substitution entries with possible wildcards in file_name, file_path, timestep_name, and atom_index
                file_pattern = subfile.get('file', None)
                timestep_pattern = subfile.get('timestep', None)
                atom_index = subfile.get('atom_index', None)
                file_path_pattern = subfile.get('file_path', None)

                # Get MultiIndex levels
                idx = atom_data.dataframe.index

                # Build masks for each level
                masks = []
                # file_name
                if file_pattern is not None:
                    if any(char in file_pattern for char in ['*', '?', '[']):
                        masks.append([fnmatch.fnmatch(str(val), file_pattern)
                                     for val in idx.get_level_values('file_name')])
                    else:
                        masks.append(
                            [str(val) == file_pattern for val in idx.get_level_values('file_name')])
                # file_path
                if file_path_pattern is not None:
                    if any(char in file_path_pattern for char in ['*', '?', '[']):
                        masks.append([fnmatch.fnmatch(str(val), file_path_pattern)
                                     for val in idx.get_level_values('file_path')])
                    else:
                        masks.append(
                            [str(val) == file_path_pattern for val in idx.get_level_values('file_path')])
                # timestep_name
                if timestep_pattern is not None:
                    if any(char in timestep_pattern for char in ['*', '?', '[']):
                        masks.append([fnmatch.fnmatch(str(val), timestep_pattern)
                                     for val in idx.get_level_values('timestep_name')])
                    else:
                        masks.append(
                            [str(val) == timestep_pattern for val in idx.get_level_values('timestep_name')])
                # atom_index
                if atom_index is not None:
                    masks.append(
                        [val == atom_index for val in idx.get_level_values('atom_index')])

                # Combine all masks
                if masks:
                    final_mask = pd.Series([all(vals) for vals in zip(
                        *masks)], index=atom_data.dataframe.index)
                else:
                    final_mask = pd.Series(
                        [True] * len(atom_data.dataframe), index=atom_data.dataframe.index)

                atom_data.dataframe.loc[final_mask, "alias"] = sub["name"]

    print("Finished reading")
    print("Atom data:")
    print(atom_data.dataframe)
    print("Timestep data:")
    print(frame_data.dataframe)

    timestep_names = frame_data.dataframe.index.get_level_values(
        "timestep_name").unique().tolist()
    print(timestep_names)

    # Read CSV file with name/value pairs and add columns to frame_data.dataframe
    global_constants_csv_file = next(
        (f for f in yamldata.get("files", []) if f.get(
            "type", "").lower() == "global_constants_csv"), None
    )
    if global_constants_csv_file:
        csv_path = prepend_root_if_relative(
            file_path=global_constants_csv_file["path"], root_path=args.root_path
        )
        name_value_df = pd.read_csv(csv_path)
        for _, row in name_value_df.iterrows():
            col_name = row['name']
            col_value = row['value']
            frame_data.dataframe[col_name] = col_value

    # Read CSV files with frame constants and add columns to frame_data.dataframe
    frame_constants_csv_files = [
        f for f in yamldata.get("files", []) if f.get("type", "").lower() == "frame_constants_csv"
    ]
    for frame_constants_csv_file in frame_constants_csv_files:
        csv_path = prepend_root_if_relative(
            file_path=frame_constants_csv_file["path"], root_path=args.root_path
        )
        frame_constants_df = pd.read_csv(csv_path)
        # Identify label columns (file_name, file_path, timestep_name)
        label_cols = [col for col in frame_constants_df.columns if col in [
            "file_name", "file_path", "timestep_name"]]
        value_cols = [
            col for col in frame_constants_df.columns if col not in label_cols]
        # For each row in frame_constants_df, match to frame_data and set value
        for _, row in frame_constants_df.iterrows():
            mask = pd.Series([True] * len(frame_data.dataframe),
                             index=frame_data.dataframe.index)
            for label in label_cols:
                val = row[label]
                if pd.isna(val):
                    continue
                idx_vals = frame_data.dataframe.index.get_level_values(label)
                if any(char in str(val) for char in ['*', '?', '[']):
                    mask &= [fnmatch.fnmatch(str(idx_val), str(val))
                             for idx_val in idx_vals]
                else:
                    mask &= [str(idx_val) == str(val) for idx_val in idx_vals]
            for value_col in value_cols:
                frame_data.dataframe.loc[mask, value_col] = row[value_col]
        # Fill missing values with NaN (already default in pandas, but ensure columns exist)
        for value_col in value_cols:
            if value_col not in frame_data.dataframe.columns:
                frame_data.dataframe[value_col] = pd.NA

    measure = mr.Measure()

    def resolve_atom(label, timestep_name, m):
        file_val = m.get("file", None)
        timestep_val = timestep_name if timestep_name is not None else None

        def match_all(val):
            return val is None or (isinstance(val, float) and pd.isna(val)) or (isinstance(val, str) and val.strip() == "")

        try:
            idx = int(label)
            file_mask = np.ones(len(atom_data.dataframe), dtype=bool) if match_all(file_val) else (
                atom_data.dataframe.index.get_level_values("file_name") == file_val)
            timestep_mask = np.ones(len(atom_data.dataframe), dtype=bool) if match_all(timestep_val) else (
                atom_data.dataframe.index.get_level_values("timestep_name") == timestep_val)
            atom_index_mask = (
                atom_data.dataframe.index.get_level_values("atom_index") == idx)
            mask = file_mask & timestep_mask & atom_index_mask
        except ValueError:
            file_mask = np.ones(len(atom_data.dataframe), dtype=bool) if match_all(file_val) else (
                atom_data.dataframe.index.get_level_values("file_name") == file_val)
            timestep_mask = np.ones(len(atom_data.dataframe), dtype=bool) if match_all(timestep_val) else (
                atom_data.dataframe.index.get_level_values("timestep_name") == timestep_val)
            alias_mask = (atom_data.dataframe["alias"] == str(label))
            mask = file_mask & timestep_mask & alias_mask

        matches = atom_data.dataframe.index[mask]
        if len(matches) == 0:
            print(
                f"Could not resolve atom '{label}' for timestep '{timestep_name}' and file '{file_val}'")
            return None
        return matches[0]

    # Ensure frame_data.dataframe has a column for each measurement
    if "measurements" in yamldata:
        timestep_names = frame_data.dataframe.index.get_level_values(
            "timestep_name").unique().tolist()
        for mtype in ["distance", "angle", "dihedral"]:
            if mtype in yamldata["measurements"]:
                for m in yamldata["measurements"][mtype]:
                    results = []
                    for timestep_name in timestep_names:
                        try:
                            if mtype == "distance":
                                idx_a = resolve_atom(m["a"], timestep_name, m)
                                idx_b = resolve_atom(m["b"], timestep_name, m)
                                if idx_a is None or idx_b is None:
                                    val = None
                                else:
                                    val = measure.distance(
                                        atom_data, idx_a, idx_b)
                            elif mtype == "angle":
                                idx_a = resolve_atom(m["a"], timestep_name, m)
                                idx_b = resolve_atom(m["b"], timestep_name, m)
                                idx_c = resolve_atom(m["c"], timestep_name, m)
                                if None in (idx_a, idx_b, idx_c):
                                    val = None
                                else:
                                    val = measure.angle(
                                        atom_data, idx_a, idx_b, idx_c)
                            elif mtype == "dihedral":
                                idx_a = resolve_atom(m["a"], timestep_name, m)
                                idx_b = resolve_atom(m["b"], timestep_name, m)
                                idx_c = resolve_atom(m["c"], timestep_name, m)
                                idx_d = resolve_atom(m["d"], timestep_name, m)
                                if None in (idx_a, idx_b, idx_c, idx_d):
                                    val = None
                                else:
                                    val = measure.dihedral(
                                        atom_data, idx_a, idx_b, idx_c, idx_d)
                            results.append(val)
                        except Exception as e:
                            print(
                                f"Error in measurement '{m['name']}' for timestep '{timestep_name}': {e}")
                            results.append(None)
                    frame_data.dataframe[m["name"]] = results
    print("Finished measuring")
    print("Timestep data:")
    print(frame_data.dataframe)

    if "calc" in yamldata:
        runner = ccr.CustomCalculationRunner(frame_data)
        runner.run(yamldata["calc"])
        print("Custom calculations complete.")
        print(frame_data.dataframe)

    class FrameDataExporter:
        def __init__(self, frame_data):
            self.frame_data = frame_data

        def export_csv_tuples(self, file_path):
            df = self.frame_data.dataframe.copy()
            df = df.reset_index()
            df['tipe'] = list(
                df[self.frame_data.dataframe.index.names].apply(tuple, axis=1))
            cols = [
                'tipe'] + [col for col in df.columns if col not in self.frame_data.dataframe.index.names]
            df_export = df[cols]
            df_export.to_csv(file_path, index=False)

        def export_csv_multiindex(self, file_path):
            df = self.frame_data.dataframe.copy()
            df.to_csv(file_path)

    exporter = FrameDataExporter(frame_data)

    for one_output in yamldata.get('output', []):
        if 'file' in one_output:
            for file in one_output['file']:
                file_type = file.get('type', '').lower()
                file_path = prepend_root_if_relative(
                    file_path=file['path'], root_path=args.root_path)
                if file_type == "csv_tuples":
                    exporter.export_csv_tuples(file_path)
                elif file_type == "csv":
                    exporter.export_csv_multiindex(file_path)
                else:
                    raise IndexError(f"{file_type}: Unknown file type")

    # Visualization: Plot measurement n vs m with series grouped by file_name and timestep_name
    # Example YAML specification for such a plot:
    # output:
    #   - graph:
    #       - type: scatter_plot
    #         x: measurement_m
    #         y: measurement_n
    #         series_by: file_name      # group series by file_name
    #         parallel_by: timestep_name # parallel series by timestep_name
    #         file: ../series_plot.tiff
    #         file_format: tiff
    #         dpi: 300

    for one_output in yamldata.get('output', []):
        if 'graph' in one_output:
            for graph in one_output['graph']:
                if graph['type'].lower() == "scatter_plot":
                    x_col = graph['x']
                    y_col = graph['y']
                    series_by = graph.get('series_by', None)
                    parallel_by = graph.get('parallel_by', None)
                    df = frame_data.dataframe.reset_index()

                    if series_by and parallel_by:
                        # Group by parallel_by, then plot series_by within each group
                        unique_parallel = df[parallel_by].unique()
                        fig, ax = plt.subplots(
                            figsize=graph.get("figsize", (8, 6)))
                        for pval in unique_parallel:
                            subdf = df[df[parallel_by] == pval]
                            for sval in subdf[series_by].unique():
                                ssubdf = subdf[subdf[series_by] == sval]
                                ax.plot(ssubdf[x_col], ssubdf[y_col], marker='o',
                                        label=f"{series_by}: {sval}, {parallel_by}: {pval}")
                        ax.set_xlabel(x_col)
                        ax.set_ylabel(y_col)
                        ax.legend()
                        fig.savefig(prepend_root_if_relative(file_path=graph['file'], root_path=args.root_path), dpi=graph.get(
                            "dpi", 300), format=graph.get("file_format", "tiff"))
                    elif series_by:
                        fig, ax = plt.subplots(
                            figsize=graph.get("figsize", (8, 6)))
                        for sval in df[series_by].unique():
                            ssubdf = df[df[series_by] == sval]
                            ax.plot(ssubdf[x_col], ssubdf[y_col],
                                    marker='o', label=f"{series_by}: {sval}")
                        ax.set_xlabel(x_col)
                        ax.set_ylabel(y_col)
                        ax.legend()
                        fig.savefig(prepend_root_if_relative(file_path=graph['file'], root_path=args.root_path), dpi=graph.get(
                            "dpi", 300), format=graph.get("file_format", "tiff"))
                    else:
                        plot = df.plot.scatter(x=x_col, y=y_col)
                        fig = plot.get_figure()
                        fig.savefig(prepend_root_if_relative(file_path=graph['file'], root_path=args.root_path), dpi=graph.get(
                            "dpi", 300), format=graph.get("file_format", "tiff"))


if __name__ == "__main__":
    main()
