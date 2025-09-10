import argparse
from tokenize import group
import strictyaml as sy
import qmanalysis.xyzreader as xr
import qmanalysis.yamlreader as yr
import qmanalysis.measure as mr
from qmanalysis.containers import AtomData, FrameData, MeasurementData
import matplotlib.pyplot as plt
from pathlib import Path
import fnmatch
import re
import asteval as av
from qmanalysis.globalconstantsreader import GlobalConstantsFile
import pandas as pd
import numpy as np
import scipy

import qmanalysis.customcalculationrunner as ccr
from qmanalysis.gaussianoutreader import GaussianOutFile


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
        print(str(list(p.resolve() for p in Path(root_path).glob(pattern))))
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

        print(args.root_path)
        yamlparser = yr.YAMLFile()
        # yamlparser.load_string(testyaml)
        # prepend_root_if_relative(file_path=args.inputfile, root_path=args.root_path))
        yamlparser.load_file(prepend_root_if_relative(
            file_path=args.inputfile, root_path=args.root_path))
        # validate_measurements(data)
    except (sy.YAMLError, ValueError) as e:
        print(f"Validation error: {e}")

    atom_data = AtomData()
    frame_data = FrameData()
    mesurement_data = MeasurementData()

    print(yamlparser.get_data())

    yamldata = yamlparser.get_data()
    print(yamldata)
    for file in yamldata["files"]:
        ftype = file["type"].lower()
        # print("Reading file: "+file['name'])
        if ftype == "xyz":
            if "glob" in file and file["glob"]:
                globbed_files = list(prepend_root_if_relative_and_glob(
                    file_path=file["path"], root_path=args.root_path))
                if len(globbed_files) == 1:
                    xr.XYZFile(atom_data, frame_data, file_path=globbed_files[0],
                               file_name=file["name"], timestep_name=file.get("timestep", None))
                else:
                    for gfile in globbed_files:
                        xr.XYZFile(atom_data, frame_data, file_path=gfile,
                                   file_name=file["name"]+Path(gfile).stem, timestep_name=file.get("timestep", None))
            else:
                xr.XYZFile(atom_data, frame_data, file_path=prepend_root_if_relative(
                    file_path=file["path"], root_path=args.root_path), file_name=file["name"], timestep_name=file.get("timestep", None))
        elif ftype == "gaussian_out":
            if "glob" in file and file["glob"]:
                globbed_files = list(prepend_root_if_relative_and_glob(
                    file_path=file["path"], root_path=args.root_path))
                if len(globbed_files) == 1:
                    GaussianOutFile(atom_data, frame_data, file_path=globbed_files[0],
                                    file_name=file.get("name", None),
                                    timestep_name=file.get("timestep", None))
                else:
                    for gfile in globbed_files:
                        GaussianOutFile(atom_data, frame_data, file_path=gfile,
                                        file_name=file.get(
                                            "name", None) + Path(gfile).stem,
                                        timestep_name=file.get("timestep", None))
            else:
                GaussianOutFile(atom_data, frame_data, file_path=prepend_root_if_relative(
                    file_path=file["path"], root_path=args.root_path),
                    file_name=file.get("name", None),
                    timestep_name=file.get("timestep", None))
        elif ftype == "global_constants_csv":
            global_constants = GlobalConstantsFile(
                file_path=prepend_root_if_relative(file["path"], root_path=args.root_path))
        elif ftype == "per_file_constants_csv":
            per_file_constants = pd.read_csv(prepend_root_if_relative(
                file["path"], root_path=args.root_path))
            # Join per_file_constants to frame_data.dataframe on file_name, preserving multiindex
            frame_data.dataframe = frame_data.dataframe.join(
                per_file_constants.set_index('file_name'), on='file_name', rsuffix='_perfile')
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
        file_val = m
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
    # Meier
    # Ensure frame_data.dataframe has a column for each measurement

    # If no timestep_name is specified in the measurement, perform for all timestep_names
    if "measurements" in yamldata:
        all_timestep_names = frame_data.dataframe.index.get_level_values(
            "timestep_name").unique().tolist()
        for mtype in ["distance", "angle", "dihedral"]:
            if mtype in yamldata["measurements"]:
                for m in yamldata["measurements"][mtype]:
                    # If 'timestep' is specified in measurement, use only that, else use all
                    if "timestep" in m and m["timestep"] is not None:
                        timestep_names = [m["timestep"]]
                    else:
                        timestep_names = all_timestep_names
                    results = []
                    for idx, row in frame_data.dataframe.iterrows():
                        timestep_name = idx[frame_data.dataframe.index.names.index(
                            "timestep_name")]
                        if timestep_name not in timestep_names:
                            results.append(pd.NA)
                        else:
                            try:
                                if mtype == "distance":
                                    idx_a = resolve_atom(
                                        m["a"], None, idx[frame_data.dataframe.index.names.index(
                                            "file_name")])
                                    idx_b = resolve_atom(
                                        m["b"], None, idx[frame_data.dataframe.index.names.index(
                                            "file_name")])
                                    if idx_a is None or idx_b is None:
                                        val = pd.NA
                                    else:
                                        val = measure.distance(
                                            atom_data, idx_a, idx_b)
                                elif mtype == "angle":
                                    idx_a = resolve_atom(
                                        m["a"], None, idx[frame_data.dataframe.index.names.index(
                                            "file_name")])
                                    idx_b = resolve_atom(
                                        m["b"], None, idx[frame_data.dataframe.index.names.index(
                                            "file_name")])
                                    idx_c = resolve_atom(
                                        m["c"], None, idx[frame_data.dataframe.index.names.index(
                                            "file_name")])
                                    if None in (idx_a, idx_b, idx_c):
                                        val = pd.NA
                                    else:
                                        val = measure.angle(
                                            atom_data, idx_a, idx_b, idx_c)
                                elif mtype == "dihedral":
                                    idx_a = resolve_atom(
                                        m["a"], None, idx[frame_data.dataframe.index.names.index(
                                            "file_name")])
                                    idx_b = resolve_atom(
                                        m["b"], None, idx[frame_data.dataframe.index.names.index(
                                            "file_name")])
                                    idx_c = resolve_atom(
                                        m["c"], None, idx[frame_data.dataframe.index.names.index(
                                            "file_name")])
                                    idx_d = resolve_atom(
                                        m["d"], None, idx[frame_data.dataframe.index.names.index(
                                            "file_name")])
                                    if None in (idx_a, idx_b, idx_c, idx_d):
                                        val = pd.NA
                                    else:
                                        val = measure.dihedral(
                                            atom_data, idx_a, idx_b, idx_c, idx_d)
                                results.append(val)
                            except Exception as e:
                                print(
                                    f"Error in measurement '{m['name']}' for name '{idx[frame_data.dataframe.index.names.index(
                                        "file_name")]}': {e}")
                                results.append(None)
                    # If only one timestep, assign scalar, else assign list
                    if len(results) == 1:
                        frame_data.dataframe[m["name"]] = results[0]
                    else:
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

        def export_csv_tuples(self, file_path, include_raw_data=False):
            df = self.frame_data.dataframe.copy()
            df = df.reset_index()
            df['label'] = list(
                df[self.frame_data.dataframe.index.names].apply(tuple, axis=1))
            cols = [
                'label'] + [col for col in df.columns if col not in self.frame_data.dataframe.index.names]
            if not include_raw_data and 'raw_data' in cols:
                cols.remove('raw_data')
            df_export = df[cols]
            df_export.to_csv(file_path, index=False)

        def export_csv_multiindex(self, file_path, include_raw_data=False):
            df = self.frame_data.dataframe.copy()
            if not include_raw_data and 'raw_data' in df.columns:
                df = df.drop(columns=['raw_data'])
            df.to_csv(file_path)

        def export_xls_tuples(self, file_path, include_raw_data=False):
            df = self.frame_data.dataframe.copy()
            df = df.reset_index()
            df['label'] = list(
                df[self.frame_data.dataframe.index.names].apply(tuple, axis=1))
            cols = [
                'label'] + [col for col in df.columns if col not in self.frame_data.dataframe.index.names]
            if not include_raw_data and 'raw_data' in cols:
                cols.remove('raw_data')
            df_export = df[cols]
            df_export.to_excel(file_path, index=False)

        def export_xls_multiindex(self, file_path, include_raw_data=False):
            df = self.frame_data.dataframe.copy()
            if not include_raw_data and 'raw_data' in df.columns:
                df = df.drop(columns=['raw_data'])
            df.to_excel(file_path)

    exporter = FrameDataExporter(frame_data)

    for one_output in yamldata.get('output', []):
        if 'file' in one_output:
            for file in one_output['file']:
                file_type = file.get('type', '').lower()
                file_path = prepend_root_if_relative(
                    file_path=file['path'], root_path=args.root_path)
                # CSV export
                if file_type == 'csv':
                    if file.get('multiindex', False):
                        exporter.export_csv_multiindex(file_path)
                    else:
                        exporter.export_csv_tuples(file_path)
                # XLS export
                if file_type == 'xls' or file_type == 'xlsx':
                    if file.get('multiindex', False):
                        exporter.export_xls_multiindex(file_path)
                    else:
                        exporter.export_xls_tuples(file_path)
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

                    def resolve_columns(col_spec, df):
                        if isinstance(col_spec, list):
                            return [df[col] for col in col_spec]
                        elif isinstance(col_spec, str):
                            return [df[col_spec]]
                        elif isinstance(col_spec, int):
                            return [df.iloc[:, col_spec]]
                        else:
                            raise ValueError(
                                f"Unsupported column spec: {col_spec}")

                    x_spec = graph['x']
                    y_spec = graph['y']
                    x_label = graph.get('x_label', None)
                    y_label = graph.get('y_label', None)
                    series_by = graph.get('series_by', "file_name")
                    parallel_by = graph.get('parallel_by', None)
                    df = frame_data.dataframe.reset_index()

                    timestep_name = graph.get("timestep_name", None)
                    if timestep_name is not None:
                        df = df[df["timestep_name"] == timestep_name]

                    x_cols = resolve_columns(x_spec, df)
                    y_cols = resolve_columns(y_spec, df)

                    if len(x_cols) == 1 and len(y_cols) > 1:
                        x_cols = x_cols * len(y_cols)
                    if len(y_cols) == 1 and len(x_cols) > 1:
                        y_cols = y_cols * len(x_cols)

                    fig, ax = plt.subplots(
                        figsize=graph.get("figsize", (8, 6)))

                    # Get min/max for both axes
                    x_min = min([df[col.name].min() for col in x_cols])
                    x_max = max([df[col.name].max() for col in x_cols])
                    y_min = min([df[col.name].min() for col in y_cols])
                    y_max = max([df[col.name].max() for col in y_cols])

                    # Diagonal mode: square aspect ratio, identical axis scales, and diagonal line
                    if graph.get('diagonal', False):
                        ax.set_aspect('equal', adjustable='box')
                        x_min = min(x_min, y_min)
                        y_min = x_min
                        x_max = max(x_max, y_max)
                        y_max = x_max

                    # Extend axes by 5% on each side
                    x_axis_range = x_max - x_min
                    y_axis_range = y_max - y_min
                    x_pad = x_axis_range * 0.05
                    y_pad = y_axis_range * 0.05
                    ax.set_xlim(x_min - x_pad, x_max + x_pad)
                    ax.set_ylim(y_min - y_pad, y_max + y_pad)

                    if graph.get('diagonal', False):
                        # Add diagonal dashed line
                        ax.plot([x_min - x_pad, x_max + x_pad], [y_min - y_pad,
                                y_max + y_pad], linestyle='--', color='gray', linewidth=1)

                    label_offset_percentage = 0.03  # 3% of axis range

                    # Define marker symbols to cycle through
                    marker_symbols = ['x', '.', '+', '1', '2',
                                      '3', '4', '^', 'v', '<', '>', 'd', 's', 'p', '|', '_']
                    marker_fillstyles = {
                        '^': 'none', 'v': 'none', '<': 'none', '>': 'none',
                        'd': 'none', 's': 'none', 'p': 'none',
                        'x': 'full', '+': 'full', '|': 'full', '_': 'full',
                        '1': 'full', '2': 'full', '3': 'full', '4': 'full', '.': 'none'
                    }
                    unique_files = df["file_name"].unique()
                    file_marker_map = {fname: marker_symbols[i % len(
                        marker_symbols)] for i, fname in enumerate(unique_files)}

                    # Store all marker and label positions for optimization
                    # {marker_position: (x,y), marker_type: ".", label_position: (x,y), label_text: "label"}
                    marker_and_label_data = []

                    # Plot each file_name as a separate series with its marker
                    for fname in unique_files:
                        subdf = df[df["file_name"] == fname]
                        marker = file_marker_map[fname]
                        fillstyle = marker_fillstyles.get(marker, 'full')
                        label_bboxes = []  # Store bounding boxes of placed labels
                        group_threshold = 0.20  # threshold for grouping close markers
                        group_centers = []
                        label_offset_data = x_axis_range * label_offset_percentage
                        # 2.6% of axis range for stacking labels
                        label_stack_offset = x_axis_range * label_offset_percentage
                        # Collect all marker and label positions for optimization
                        all_marker_positions = [(row[xcol.name], row[ycol.name]) for i, (xcol, ycol) in enumerate(
                            zip(x_cols, y_cols)) for idx, row in subdf.iterrows()]
                        all_label_positions = []
                        label_texts = []
                        label_marker_indices = []
                        # First pass: calculate initial label positions
                        for i, (xcol, ycol) in enumerate(zip(x_cols, y_cols)):
                            for idx, row in subdf.iterrows():
                                x = row[xcol.name]
                                y = row[ycol.name]
                                label_text = str(row[series_by])
                                x_offset = x + label_offset_data
                                y_offset = y
                                all_label_positions.append(
                                    (x_offset, y_offset))
                                label_texts.append(label_text)
                                label_marker_indices.append((x, y))
                        used_positions = prune_close_positions(
                            all_marker_positions, group_threshold, x_axis_range, y_axis_range)

                        # Map used positions back to their original labels
                        for i, pos in enumerate(used_positions):
                            if pos:
                                marker_and_label_data.append({
                                    "marker_position": all_marker_positions[i],
                                    "marker_type": marker,
                                    "label_position": all_label_positions[i],
                                    "label_text": label_texts[i]
                                })
                            else:
                                marker_and_label_data.append({
                                    "marker_position": all_marker_positions[i],
                                    "marker_type": marker
                                })
                                # Calculate axis lengths and scales for optimal label placement
                                x_axis_length = x_max - x_min if x_max != x_min else 1.0
                                y_axis_length = y_max - y_min if y_max != y_min else 1.0
                                bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
                                width, height = bbox.width, bbox.height
                                aspect_ratio = width / height if height != 0 else 1.0
                                x_scale = aspect_ratio
                                y_scale = 1.0

                    # Second pass: optimize label positions
                    for i in range(len(marker_and_label_data)):
                        if "label_position" in marker_and_label_data[i]:
                            marker_pos = marker_and_label_data[i]["marker_position"]
                            label_pos = marker_and_label_data[i]["label_position"]
                            other_positions = [pos for j, pos in enumerate(
                                [d["marker_position"] for d in marker_and_label_data if "marker_position" in d] + [d["label_position"] for d in marker_and_label_data if "label_position" in d]) if j != i and pos != marker_pos and pos != label_pos]
                            new_label_pos = find_optimal_label_position(
                                marker_pos, label_pos, other_positions,
                                x_axis_length=x_axis_length,
                                y_axis_length=y_axis_length,
                                x_scale=x_scale,
                                y_scale=y_scale,
                                xlim=(x_min - x_pad, x_max + x_pad),
                                ylim=(y_min - y_pad, y_max + y_pad),
                                diagonal_line=graph.get('diagonal', False)
                            )
                            marker_and_label_data[i]["label_position"] = new_label_pos
                    # for i in range(len(marker_and_label_data)):
                    #     if "label_position" in marker_and_label_data[i]:
                    #         marker_pos = marker_and_label_data[i]["marker_position"]
                    #         label_pos = marker_and_label_data[i]["label_position"]
                    #         other_positions = [pos for j, pos in enumerate(
                    #             [d["marker_position"] for d in marker_and_label_data if "marker_position" in d] + [d["label_position"] for d in marker_and_label_data if "label_position" in d]) if j != i and pos != marker_pos and pos != label_pos]
                    #         new_label_pos = find_optimal_label_position(
                    #             marker_pos, label_pos, other_positions,
                    #             x_axis_length=x_axis_length,
                    #             y_axis_length=y_axis_length,
                    #             x_scale=x_scale,
                    #             y_scale=y_scale,
                    #             xlim=(x_min - x_pad, x_max + x_pad),
                    #             ylim=(y_min - y_pad, y_max + y_pad),
                    #             diagonal_line=graph.get('diagonal', False)
                    #         )
                    #         marker_and_label_data[i]["label_position"] = new_label_pos
                    # for i in range(len(marker_and_label_data)):
                    #     if "label_position" in marker_and_label_data[i]:
                    #         marker_pos = marker_and_label_data[i]["marker_position"]
                    #         label_pos = marker_and_label_data[i]["label_position"]
                    #         other_positions = [pos for j, pos in enumerate(
                    #             all_marker_positions + all_label_positions) if j != i and pos != marker_pos and pos != label_pos]
                    #         new_label_pos = find_optimal_label_position(
                    #             marker_pos, label_pos, other_positions)
                    #         marker_and_label_data[i]["label_position"] = new_label_pos

                    # optimized_label_positions = []
                    # for i, (marker_pos, label_pos) in enumerate(zip(label_marker_indices, all_label_positions)):
                    #     # Exclude current marker and label from other_positions
                    #     other_positions = [pos for j, pos in enumerate(
                    #         all_marker_positions + all_label_positions) if j != i and pos != marker_pos and pos != label_pos]
                    #     new_label_pos = find_optimal_label_position(
                    #         marker_pos, label_pos, other_positions)
                    #     optimized_label_positions.append(new_label_pos)
                    # Third pass: plot markers and optimized labels
                    label_bboxes = []
                    print(
                        f'Plotting markers and labels... ')
                    for i, marker_and_label in enumerate(marker_and_label_data):
                        print(f'{marker_and_label}')
                        (x, y) = marker_and_label["marker_position"]
                        if marker_fillstyles[marker_and_label["marker_type"]] == 'none':
                            ax.scatter(
                                x, y, marker=marker_and_label["marker_type"], facecolors='none', edgecolors='black', s=30, linewidths=0.5)
                        else:
                            ax.scatter(
                                x, y, marker=marker_and_label["marker_type"], color='black', s=30, linewidths=0.5)
                    # for i, marker_and_label in enumerate(marker_and_label_data):  #zip(label_texts, optimized_label_positions)):
                        if "label_position" in marker_and_label:
                            print(
                                f'Placing label: {marker_and_label["label_text"]} at {marker_and_label["label_position"]}')
                            (x_opt, y_opt) = marker_and_label["label_position"]
                            label_text = marker_and_label["label_text"]
                            print(
                                f'Placing label: {label_text} at {x_opt}, {y_opt}')
                            text_obj = ax.text(
                                x_opt, y_opt, label_text, fontsize=8, va='center', ha='center')
                            print(f'text object: {text_obj}')
                            plt.draw()
                            bbox = text_obj.get_window_extent(
                                renderer=fig.canvas.get_renderer())
                # Set axis labels
                if x_label:
                    ax.set_xlabel(', '.join(x_label) if isinstance(
                        x_label, list) else str(x_label))
                else:
                    ax.set_xlabel(', '.join([col.name for col in x_cols]))
                if y_label:
                    ax.set_ylabel(', '.join(y_label) if isinstance(
                        y_label, list) else str(y_label))
                else:
                    ax.set_ylabel(', '.join([col.name for col in y_cols]))

                if "title" in graph and graph["title"]:
                    ax.set_title(graph["title"])

                # fig.tight_layout()

                file_base = prepend_root_if_relative(
                    file_path=graph['file'], root_path=args.root_path)
                file_formats = graph.get("file_format", "tiff")
                if isinstance(file_formats, list):
                    formats_list = file_formats
                else:
                    formats_list = [file_formats]
                for fmt in formats_list:
                    ext = f".{fmt.lower()}"
                    file_out = str(file_base)
                    if not file_out.lower().endswith(ext):
                        file_out += ext
                    fig.savefig(
                        file_out,
                        dpi=graph.get("dpi", 300),
                        format=fmt
                    )


def prune_close_positions(positions, threshold, x_scale=1.0, y_scale=1.0):
    """
    Given a list of positions [(x1, y1), (x2, y2), ...] and a threshold,
    return a new list where positions closer than the threshold are grouped
    and replaced by their centroid.
    """

    positions = np.array(positions)
    used = np.ones(len(positions), dtype=bool)
    if (len(positions) == 1):
        return used

    for i in range(len(positions)):
        for j in range(len(positions)):
            if i != j:
                # print(
                #     f'comparing {i} and {j}: {positions[i]} vs {positions[j]} ({used[i]}, {used[j]}) => {positions[i] - positions[j]} => {((positions[i] - positions[j]) / np.array([x_scale, y_scale]))} => {np.linalg.norm((positions[i] - positions[j]) / np.array([x_scale, y_scale]))}')

                if used[i] and used[j] and np.linalg.norm((positions[i] - positions[j]) / np.array([x_scale, y_scale])) < threshold:
                    used[j] = False
    print(f"Pruning positions: {positions} -> {used}")
    return used


def find_optimal_label_position(marker_pos, label_pos, other_positions, x_axis_length=1.0, y_axis_length=1.0, x_scale=1.0, y_scale=1.0, xlim=None, ylim=None, diagonal_line=False):
    """
    Given marker_pos (x, y), label_pos (lx, ly), and a list of other_positions [(x, y), ...],
    calculate a circle around the marker with radius equal to the current label distance,
    and return the position on that circle farthest from any other element on the circle.
    The circle is scaled according to x_axis_length, y_axis_length, x_scale, and y_scale.
    Also avoids plot borders and diagonal line if present.
    """
    x, y = marker_pos
    lx, ly = label_pos
    # Compute scaled radius so that the circle is correct in image coordinates
    dx = (lx - x) * x_scale / x_axis_length
    dy = (ly - y) * y_scale / y_axis_length
    radius = np.hypot(dx, dy)
    if radius == 0:
        radius = 0.1
    # Calculate angles of all other elements relative to marker, with scaling
    angles = []
    for ox, oy in other_positions:
        ddx = (ox - x) * x_scale / x_axis_length
        ddy = (oy - y) * y_scale / y_axis_length
        if ddx == 0 and ddy == 0:
            continue
        angle = np.arctan2(ddy, ddx)
        angles.append(angle)
    # Add border avoidance: treat border as forbidden angles
    if xlim is not None and ylim is not None:
        border_angles = []
        # Sample points on the border and project to angles
        for bx in [xlim[0], xlim[1]]:
            for by in np.linspace(ylim[0], ylim[1], 20):
                ddx = (bx - x) * x_scale / x_axis_length
                ddy = (by - y) * y_scale / y_axis_length
                angle = np.arctan2(ddy, ddx)
                border_angles.append(angle)
        for by in [ylim[0], ylim[1]]:
            for bx in np.linspace(xlim[0], xlim[1], 20):
                ddx = (bx - x) * x_scale / x_axis_length
                ddy = (by - y) * y_scale / y_axis_length
                angle = np.arctan2(ddy, ddx)
                border_angles.append(angle)
        angles.extend(border_angles)
    # Add diagonal line avoidance if present
    if diagonal_line and xlim is not None and ylim is not None:
        diag_angles = []
        # Sample points along the diagonal line
        for t in np.linspace(0, 1, 40):
            bx = xlim[0] + t * (xlim[1] - xlim[0])
            by = ylim[0] + t * (ylim[1] - ylim[0])
            ddx = (bx - x) * x_scale / x_axis_length
            ddy = (by - y) * y_scale / y_axis_length
            angle = np.arctan2(ddy, ddx)
            diag_angles.append(angle)
        angles.extend(diag_angles)
    angles = np.sort(angles)
    if len(angles) == 0:
        label_angle = np.arctan2(dy, dx)
        new_lx = x + radius * np.cos(label_angle) * x_axis_length / x_scale
        new_ly = y + radius * np.sin(label_angle) * y_axis_length / y_scale
        return (new_lx, new_ly)
    gaps = []
    for i in range(len(angles)):
        a1 = angles[i]
        a2 = angles[(i + 1) % len(angles)]
        gap = (a2 - a1) % (2 * np.pi)
        gaps.append((gap, a1, a2))
    max_gap, a1, a2 = max(gaps, key=lambda t: t[0])
    mid_angle = (a1 + max_gap / 2) % (2 * np.pi)
    new_lx = x + radius * np.cos(mid_angle) * x_axis_length / x_scale
    new_ly = y + radius * np.sin(mid_angle) * y_axis_length / y_scale
    return (new_lx, new_ly)


if __name__ == "__main__":
    main()
