import argparse
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
            df['tipe'] = list(
                df[self.frame_data.dataframe.index.names].apply(tuple, axis=1))
            cols = [
                'tipe'] + [col for col in df.columns if col not in self.frame_data.dataframe.index.names]
            if not include_raw_data and 'raw_data' in cols:
                cols.remove('raw_data')
            df_export = df[cols]
            df_export.to_csv(file_path, index=False)

        def export_csv_multiindex(self, file_path, include_raw_data=False):
            df = self.frame_data.dataframe.copy()
            if not include_raw_data and 'raw_data' in df.columns:
                df = df.drop(columns=['raw_data'])
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

                    # Plot each file_name as a separate series with its marker
                    for fname in unique_files:
                        subdf = df[df["file_name"] == fname]
                        marker = file_marker_map[fname]
                        fillstyle = marker_fillstyles.get(marker, 'full')
                        label_bboxes = []  # Store bounding boxes of placed labels
                        group_threshold = 0.15  # threshold for grouping close markers
                        group_centers = []
                        label_offset_data = 0.027  # horizontal offset in data units
                        label_stack_offset_px = 16  # vertical offset in pixels for stacking labels
                        for i, (xcol, ycol) in enumerate(zip(x_cols, y_cols)):
                            for idx, row in subdf.iterrows():
                                x = row[xcol.name]
                                y = row[ycol.name]
                                if fillstyle == 'none':
                                    ax.scatter(
                                        x, y, marker=marker, facecolors='none', edgecolors='black', s=30, linewidths=0.5)
                                else:
                                    ax.scatter(
                                        x, y, marker=marker, color='black', s=30, linewidths=0.5)
                                label_text = str(row[series_by])
                                # Calculate label position: horizontal offset in data, vertical offset in pixels
                                display = ax.transData.transform((x, y))
                                display_offset = (
                                    display[0] + label_offset_data * ax.figure.dpi, display[1])
                                # Stack vertically in display coordinates
                                max_stack = 10
                                for stack_level in range(max_stack):
                                    y_disp_stacked = display_offset[1] + \
                                        stack_level * label_stack_offset_px
                                    label_disp_pos = (
                                        display_offset[0], y_disp_stacked)
                                    x_offset, y_offset = ax.transData.inverted().transform(label_disp_pos)
                                    text_obj = ax.text(
                                        x_offset, y_offset, label_text, fontsize=8, va='center', ha='left')
                                    plt.draw()  # Ensure renderer is updated
                                    bbox = text_obj.get_window_extent(
                                        renderer=fig.canvas.get_renderer())
                                    overlap = False
                                    for prev_bbox in label_bboxes:
                                        if bbox.overlaps(prev_bbox):
                                            overlap = True
                                            break
                                    if not overlap:
                                        label_bboxes.append(bbox)
                                        group_centers.append((x, y, marker))
                                        break
                                    else:
                                        text_obj.remove()  # Remove overlapping label before trying next position
                            # If all stack levels overlap, skip label
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

                fig.tight_layout()

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


if __name__ == "__main__":
    main()
