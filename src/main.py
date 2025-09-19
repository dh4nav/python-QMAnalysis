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
from scipy.optimize import minimize

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
    return Path(root_path) / file_path


def prune_close_positions(positions, threshold, x_scale=1.0, y_scale=1.0):
    # Given a list of positions [(x1, y1), (x2, y2), ...] and a threshold,
    # return a new list where positions closer than the threshold are grouped and replaced by their centroid.
    positions = np.array(positions)
    used = np.ones(len(positions), dtype=bool)
    if (len(positions) == 1):
        return used
    for i in range(len(positions)):
        for j in range(len(positions)):
            if i != j:
                if used[i] and used[j] and np.linalg.norm((positions[i] - positions[j]) / np.array([x_scale, y_scale])) < threshold:
                    used[j] = False
    return used


def circler(marker_positions, other_positions, radius, x_axis_start=0.0, y_axis_start=0.0, x_axis_end=1.0, y_axis_end=1.0,  diagonal_line=False):
    def downscaler(points, x_min, x_max, y_min, y_max):
        if len(points) == 0:
            return np.array(points)
        else:
            points = np.array(points)
            points[:, 0] = (points[:, 0] - x_min) / (x_max - x_min)
            points[:, 1] = (points[:, 1] - y_min) / (y_max - y_min)
            return points

    def upscaler(points, x_min, x_max, y_min, y_max):
        if len(points) == 0:
            return np.array(points)
        else:
            points = np.array(points)
            points[:, 0] = points[:, 0] * (x_max - x_min) + x_min
            points[:, 1] = points[:, 1] * (y_max - y_min) + y_min
            return points

    # === Circle setup ===
    print("Setting up circles...")
    centers = downscaler(
        marker_positions, x_axis_start, x_axis_end, y_axis_start, y_axis_end)
    others = downscaler(
        other_positions, x_axis_start, x_axis_end, y_axis_start, y_axis_end)
    radii = np.full(len(centers), radius)

    def g(point, others=others, diagonal_line=diagonal_line):
        x, y = point
        e = 0.0
        point_strength = 500.0
        for i, ((cx, cy), r) in enumerate(zip(centers, radii)):
            d = np.linalg.norm(point - np.array([cx, cy]))
            e += point_strength / (d**6 + 1e-9)  # repulsion from circle center
        wall_strength = 500000.0
        e += wall_strength / (x + 1e-9)
        e += wall_strength / (1.0 - x + 1e-9)
        e += wall_strength / (y + 1e-9)
        e += wall_strength / (1.0 - y + 1e-9)
        if diagonal_line:
            diag_strength = 500.0
            e += diag_strength / \
                ((np.linalg.norm(
                    np.array([x, y]) - np.array([(x+y)/2.0, (x+y)/2.0]))**4) + 1e-9)
        return e

    def get_points(thetas, centers=centers, radii=radii):
        return np.array([
            [cx + r * np.cos(theta), cy + r * np.sin(theta)]
            for (cx, cy), r, theta in zip(centers, radii, thetas)
        ])

    def energy(thetas, centers=centers, radii=radii):
        points = get_points(thetas, centers, radii)
        n = len(points)
        e = 0.0
        for i in range(n):
            for j in range(i + 1, n):
                d = np.linalg.norm(points[i] - points[j])
                e += 1.0 / (d**2 + 1e-9)
        for p in points:
            e += g(p)
        return e

    def find_best_config(n_starts=20):
        best_result = None
        best_energy = np.inf
        for _ in range(n_starts):
            theta0 = np.random.rand(len(centers)) * 2 * np.pi
            res = minimize(energy, theta0, method="BFGS")
            if res.fun < best_energy:
                best_energy = res.fun
                best_result = res
        return best_result

    res = find_best_config(n_starts=30)
    print("Best energy:", res.fun)
    print("Optimal angles:", res.x)
    return upscaler(get_points(res.x), x_axis_start, x_axis_end, y_axis_start, y_axis_end)


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("inputfile", help="Main command input file")

    # Add -rp / --root-path argument (single allowed — default behavior)
    parser.add_argument(
        "-rp", "--root_path",
        type=str,
        help="Root path used for all files. If none is specified, the location of the input file is used as root path"
    )

    args = parser.parse_args()

    # --- BEGIN MOVED CODE ---
    # The following code was previously at the top-level and is now inside main()
    yamlparser = yr.YAMLFile()
    yamlparser.load_file(prepend_root_if_relative(
        file_path=args.inputfile, root_path=args.root_path))
    atom_data = AtomData()
    frame_data = FrameData()
    mesurement_data = MeasurementData()
    yamldata = yamlparser.get_data()

    # Read files
    for file in yamldata["files"]:
        ftype = file["type"].lower()
        if ftype == "xyz":
            if "glob" in file and file["glob"]:
                from glob import glob
                import os
                globbed_files = list(
                    glob(str(prepend_root_if_relative(file["path"], args.root_path))))
                if len(globbed_files) == 1:
                    xr.XYZFile(atom_data, frame_data, file_path=globbed_files[0],
                               file_name=file["name"], timestep_name=file.get("timestep", None))
                else:
                    for gfile in globbed_files:
                        xr.XYZFile(atom_data, frame_data, file_path=gfile,
                                   file_name=file["name"]+os.path.splitext(os.path.basename(gfile))[0], timestep_name=file.get("timestep", None))
            else:
                xr.XYZFile(atom_data, frame_data, file_path=prepend_root_if_relative(
                    file_path=file["path"], root_path=args.root_path), file_name=file["name"], timestep_name=file.get("timestep", None))
        elif ftype == "gaussian_out":
            if "glob" in file and file["glob"]:
                from glob import glob
                import os
                globbed_files = list(
                    glob(str(prepend_root_if_relative(file["path"], args.root_path))))
                if len(globbed_files) == 1:
                    GaussianOutFile(atom_data, frame_data, file_path=globbed_files[0],
                                    file_name=file.get("name", None), timestep_name=file.get("timestep", None))
                else:
                    for gfile in globbed_files:
                        GaussianOutFile(atom_data, frame_data, file_path=gfile,
                                        file_name=(
                                            file.get("name", "") or "")+os.path.splitext(os.path.basename(gfile))[0],
                                        timestep_name=file.get("timestep", None))
            else:
                GaussianOutFile(atom_data, frame_data, file_path=prepend_root_if_relative(
                    file_path=file["path"], root_path=args.root_path), file_name=file.get("name", None), timestep_name=file.get("timestep", None))
        elif ftype == "global_constants_csv":
            global_constants = GlobalConstantsFile(
                file_path=prepend_root_if_relative(file["path"], root_path=args.root_path))
        elif ftype == "per_file_constants_csv":
            per_file_constants = pd.read_csv(prepend_root_if_relative(
                file["path"], root_path=args.root_path))
            frame_data.dataframe = frame_data.dataframe.join(
                per_file_constants.set_index('file_name'), on='file_name', rsuffix='_perfile')

    # Substitutions
    if "substitutions" in yamldata:
        for sub in yamldata["substitutions"]:
            for subfile in sub['entries']:
                file_pattern = subfile.get('file', None)
                timestep_pattern = subfile.get('timestep', None)
                atom_index = subfile.get('atom_index', None)
                file_path_pattern = subfile.get('file_path', None)
                idx = atom_data.dataframe.index
                masks = []
                if file_pattern is not None:
                    if any(char in file_pattern for char in ['*', '?', '[']):
                        masks.append([fnmatch.fnmatch(str(val), file_pattern)
                                     for val in idx.get_level_values('file_name')])
                    else:
                        masks.append(
                            [str(val) == file_pattern for val in idx.get_level_values('file_name')])
                if file_path_pattern is not None:
                    if any(char in file_path_pattern for char in ['*', '?', '[']):
                        masks.append([fnmatch.fnmatch(str(val), file_path_pattern)
                                     for val in idx.get_level_values('file_path')])
                    else:
                        masks.append(
                            [str(val) == file_path_pattern for val in idx.get_level_values('file_path')])
                if timestep_pattern is not None:
                    if any(char in timestep_pattern for char in ['*', '?', '[']):
                        masks.append([fnmatch.fnmatch(str(val), timestep_pattern)
                                     for val in idx.get_level_values('timestep_name')])
                    else:
                        masks.append(
                            [str(val) == timestep_pattern for val in idx.get_level_values('timestep_name')])
                if atom_index is not None:
                    masks.append(
                        [val == atom_index for val in idx.get_level_values('atom_index')])
                if masks:
                    final_mask = pd.Series([all(vals) for vals in zip(
                        *masks)], index=atom_data.dataframe.index)
                else:
                    final_mask = pd.Series(
                        [True] * len(atom_data.dataframe), index=atom_data.dataframe.index)
                atom_data.dataframe.loc[final_mask, "alias"] = sub["name"]

    # Global constants
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

    # Frame constants
    frame_constants_csv_files = [
        f for f in yamldata.get("files", []) if f.get("type", "").lower() == "frame_constants_csv"
    ]
    for frame_constants_csv_file in frame_constants_csv_files:
        csv_path = prepend_root_if_relative(
            file_path=frame_constants_csv_file["path"], root_path=args.root_path
        )
        frame_constants_df = pd.read_csv(csv_path)
        label_cols = [col for col in frame_constants_df.columns if col in [
            "file_name", "file_path", "timestep_name"]]
        value_cols = [
            col for col in frame_constants_df.columns if col not in label_cols]
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
        for value_col in value_cols:
            if value_col not in frame_data.dataframe.columns:
                frame_data.dataframe[value_col] = pd.NA

    # --- Measurements ---
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

    if "measurements" in yamldata:
        all_timestep_names = frame_data.dataframe.index.get_level_values(
            "timestep_name").unique().tolist()
        for mtype in ["distance", "angle", "dihedral"]:
            if mtype in yamldata["measurements"]:
                for m in yamldata["measurements"][mtype]:
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
                                        m["a"], None, idx[frame_data.dataframe.index.names.index("file_name")])
                                    idx_b = resolve_atom(
                                        m["b"], None, idx[frame_data.dataframe.index.names.index("file_name")])
                                    if idx_a is None or idx_b is None:
                                        val = pd.NA
                                    else:
                                        val = measure.distance(
                                            atom_data, idx_a, idx_b)
                                elif mtype == "angle":
                                    idx_a = resolve_atom(
                                        m["a"], None, idx[frame_data.dataframe.index.names.index("file_name")])
                                    idx_b = resolve_atom(
                                        m["b"], None, idx[frame_data.dataframe.index.names.index("file_name")])
                                    idx_c = resolve_atom(
                                        m["c"], None, idx[frame_data.dataframe.index.names.index("file_name")])
                                    if None in (idx_a, idx_b, idx_c):
                                        val = pd.NA
                                    else:
                                        val = measure.angle(
                                            atom_data, idx_a, idx_b, idx_c)
                                elif mtype == "dihedral":
                                    idx_a = resolve_atom(
                                        m["a"], None, idx[frame_data.dataframe.index.names.index("file_name")])
                                    idx_b = resolve_atom(
                                        m["b"], None, idx[frame_data.dataframe.index.names.index("file_name")])
                                    idx_c = resolve_atom(
                                        m["c"], None, idx[frame_data.dataframe.index.names.index("file_name")])
                                    idx_d = resolve_atom(
                                        m["d"], None, idx[frame_data.dataframe.index.names.index("file_name")])
                                    if None in (idx_a, idx_b, idx_c, idx_d):
                                        val = pd.NA
                                    else:
                                        val = measure.dihedral(
                                            atom_data, idx_a, idx_b, idx_c, idx_d)
                                results.append(val)
                            except Exception as e:
                                print(
                                    f"Error in measurement '{m['name']}' for name '{idx[frame_data.dataframe.index.names.index('file_name')]}': {e}")
                                results.append(None)
                    if len(results) == 1:
                        frame_data.dataframe[m["name"]] = results[0]
                    else:
                        frame_data.dataframe[m["name"]] = results

    # --- Calculations ---
    if "calc" in yamldata:
        runner = ccr.CustomCalculationRunner(frame_data)
        runner.run(yamldata["calc"])
        print("Custom calculations complete.")
        print(frame_data.dataframe)

    # --- Exporting ---
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
                if file_type == 'csv':
                    if file.get('multiindex', False):
                        exporter.export_csv_multiindex(file_path)
                    else:
                        exporter.export_csv_tuples(file_path)
                if file_type == 'xls' or file_type == 'xlsx':
                    if file.get('multiindex', False):
                        exporter.export_xls_multiindex(file_path)
                    else:
                        exporter.export_xls_tuples(file_path)

    # --- Plotting ---
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
                    x_min = min([df[col.name].min() for col in x_cols])
                    x_max = max([df[col.name].max() for col in x_cols])
                    y_min = min([df[col.name].min() for col in y_cols])
                    y_max = max([df[col.name].max() for col in y_cols])
                    if graph.get('diagonal', False):
                        ax.set_aspect('equal', adjustable='box')
                        x_min = min(x_min, y_min)
                        y_min = x_min
                        x_max = max(x_max, y_max)
                        y_max = x_max
                        x_axis_range = x_max - x_min
                        y_axis_range = y_max - y_min
                        x_pad = x_axis_range * 0.07
                        y_pad = y_axis_range * 0.07
                        ax.set_xlim(x_min - x_pad, x_max + x_pad)
                        ax.set_ylim(y_min - y_pad, y_max + y_pad)
                        for label in ax.get_xticklabels():
                            label.set_horizontalalignment('center')
                        for label in ax.get_yticklabels():
                            label.set_verticalalignment('center')
                    else:
                        x_axis_range = x_max - x_min
                        y_axis_range = y_max - y_min
                        x_pad = x_axis_range * 0.1
                        y_pad = y_axis_range * 0.1
                        ax.set_xlim(x_min - x_pad, x_max + x_pad)
                        ax.set_ylim(y_min - y_pad, y_max + y_pad)
                    ax.set_xticklabels(
                        [f"{tick:g}" for tick in ax.get_xticks()], fontdict={'family': graph.get('xticksfont', graph.get('font', {'family': 'serif'})).get('family', 'serif'), 'size': graph.get('xticksfont', graph.get('font', {'size': 9})).get('size', 9)})
                    ax.set_yticklabels(
                        [f"{tick:g}" for tick in ax.get_yticks()], fontdict={'family': graph.get('yticksfont', graph.get('font', {'family': 'serif'})).get('family', 'serif'), 'size': graph.get('yticksfont', graph.get('font', {'size': 9})).get('size', 9)})
                    label_offset_percentage = 0.045
                    if graph.get('diagonal', False):
                        label_offset_percentage = 0.045
                        ax.plot([x_min - x_pad, x_max + x_pad], [y_min - y_pad,
                                y_max + y_pad], linestyle='--', color='gray', linewidth=1)
                    marker_symbols = ['x', '.', '+', '1', '2', '3',
                                      '4', '^', 'v', '<', '>', 'd', 's', 'p', '|', '_', '*', 'o']
                    marker_fillstyles = {'^': 'none', 'v': 'none', '<': 'none', '>': 'none', 'd': 'none', 's': 'none', 'p': 'none',
                                         'x': 'full', '+': 'full', '|': 'full', '_': 'full', '1': 'full', '2': 'full', '3': 'full', '4': 'full', '.': 'none', '*': 'none', 'o': 'none'}
                    unique_files = df["file_name"].unique()
                    file_marker_map = {fname: marker_symbols[i % len(
                        marker_symbols)] for i, fname in enumerate(unique_files)}
                    marker_and_label_data = []
                    for fname in unique_files:
                        subdf = df[df["file_name"] == fname]
                        marker_map = graph.get('marker_map', {})
                        column_marker_map = graph.get('column_marker_map', {})
                        name_column_marker_map = graph.get(
                            'name_column_marker_map', [])
                        marker_base = marker_map.get(fname, {}).get(
                            'marker', file_marker_map[fname])
                        label_override = marker_map.get(
                            fname, {}).get('label', None)
                        fillstyle = marker_fillstyles.get(marker_base, 'full')
                        group_threshold = 0.05
                        label_offset_data = x_axis_range * label_offset_percentage
                        label_stack_offset = x_axis_range * label_offset_percentage
                        all_marker_positions = [(row[xcol.name], row[ycol.name]) for i, (xcol, ycol) in enumerate(
                            zip(x_cols, y_cols)) for idx, row in subdf.iterrows()]
                        all_label_positions = []
                        label_texts = []
                        marker_list = []
                        label_marker_indices = []
                        for i, (xcol, ycol) in enumerate(zip(x_cols, y_cols)):
                            for idx, row in subdf.iterrows():
                                x = row[xcol.name]
                                y = row[ycol.name]
                                x_offset = x + label_offset_data
                                y_offset = y
                                label_text = label_override if label_override is not None else str(
                                    row[series_by])
                                marker = marker_base
                                col_keys = [xcol.name, ycol.name]
                                marker_and_label_written = False
                                for col_key in col_keys:
                                    marker_and_label_written = False
                                    for ncm in name_column_marker_map:
                                        if (fname in ncm.get('name', [])) and (col_key in ncm.get('columns', [])):
                                            if 'labeladd' in ncm['substitution']:
                                                label_texts.append(
                                                    label_text + ncm['substitution']['labeladd'])
                                            else:
                                                label_texts.append(
                                                    ncm['substitution'].get('label', label_text))
                                            marker_list.append(
                                                ncm['substitution'].get('marker', marker))
                                            marker_and_label_written = True
                                            break
                                    if not marker_and_label_written and col_key in column_marker_map:
                                        if col_key in column_marker_map:
                                            label_texts.append(
                                                column_marker_map[col_key].get('label', label_text))
                                            marker_list.append(
                                                column_marker_map[col_key].get('marker', marker))
                                            marker_and_label_written = True
                                            break
                                if not marker_and_label_written:
                                    label_texts.append(label_text)
                                    marker_list.append(marker)
                                all_label_positions.append(
                                    (x_offset, y_offset))
                                label_marker_indices.append((x, y))

                        used_positions = prune_close_positions(
                            all_marker_positions, group_threshold, x_axis_range, y_axis_range)
                        for i, pos in enumerate(used_positions):
                            if pos:
                                marker_and_label_data.append({
                                    "marker_position": all_marker_positions[i],
                                    "marker_type": marker_list[i],
                                    "label_position": all_label_positions[i],
                                    "label_text": label_texts[i]
                                })
                            else:
                                marker_and_label_data.append({
                                    "marker_position": all_marker_positions[i],
                                    "marker_type": marker_list[i]
                                })
                                x_axis_length = x_max - x_min if x_max != x_min else 1.0
                                y_axis_length = y_max - y_min if y_max != y_min else 1.0
                                bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
                                width, height = bbox.width, bbox.height
                                aspect_ratio = width / height if height != 0 else 1.0
                                x_scale = aspect_ratio
                                y_scale = 1.0
                    label_indices = [i for i, d in enumerate(
                        marker_and_label_data) if "label_position" in d]
                    if label_indices:
                        marker_positions = [
                            marker_and_label_data[i]["marker_position"] for i in label_indices]
                        initial_label_positions = [
                            marker_and_label_data[i]["label_position"] for i in label_indices]
                        other_positions = [d["marker_position"] for i, d in enumerate(
                            marker_and_label_data) if i not in label_indices]
                        opt_label_positions = circler(marker_positions, other_positions, label_offset_percentage, x_axis_start=x_min,
                                                      y_axis_start=y_min, x_axis_end=x_max, y_axis_end=y_max,  diagonal_line=graph.get('diagonal', False))
                        for idx, opt_pos in zip(label_indices, opt_label_positions):
                            marker_and_label_data[idx]["label_position"] = opt_pos
                    label_bboxes = []
                    print(f'Plotting markers and labels... ')
                    for i, marker_and_label in enumerate(marker_and_label_data):
                        print(f'M1: {marker_and_label}')
                        (x, y) = marker_and_label["marker_position"]
                        if marker_fillstyles[marker_and_label["marker_type"]] == 'none':
                            ax.scatter(
                                x, y, marker=marker_and_label["marker_type"], facecolors='none', edgecolors='black', s=30, linewidths=0.5)
                        else:
                            ax.scatter(
                                x, y, marker=marker_and_label["marker_type"], color='black', s=30, linewidths=0.5)
                        if "label_position" in marker_and_label:
                            print(
                                f'Placing label: {marker_and_label["label_text"]} at {marker_and_label["label_position"]}')
                            (x_opt, y_opt) = marker_and_label["label_position"]
                            label_text = marker_and_label["label_text"]
                            print(
                                f'Placing label: {label_text} at {x_opt}, {y_opt}')
                            text_obj = ax.text(
                                x_opt, y_opt, label_text, fontdict={'family': graph.get('labelfont', graph.get('font', {'family': 'serif'})).get('family', 'serif'), 'size': graph.get('labelfont', graph.get('font', {'size': 8})).get('size', 8)}, va='center', ha='center')
                            print(f'text object: {text_obj}')
                            plt.draw()
                            bbox = text_obj.get_window_extent(
                                renderer=fig.canvas.get_renderer())
                    if x_label:
                        ax.set_xlabel(', '.join(x_label) if isinstance(
                            x_label, list) else str(x_label), fontdict={'family': graph.get('xlabelfont', graph.get('font', {'family': 'serif'})).get('family', 'serif'), 'size': graph.get('xlabelfont', graph.get('font', {'size': 11})).get('size', 11)})
                    else:
                        ax.set_xlabel(
                            ', '.join([col.name for col in x_cols]), fontdict={'family': graph.get('xlabelfont', graph.get('font', {'family': 'serif'})).get('family', 'serif'), 'size': graph.get('xlabelfont', graph.get('font', {'size': 11})).get('size', 11)})
                    if y_label:
                        ax.set_ylabel(', '.join(y_label) if isinstance(
                            y_label, list) else str(y_label), fontdict={'family': graph.get('ylabelfont', graph.get('font', {'family': 'serif'})).get('family', 'serif'), 'size': graph.get('ylabelfont', graph.get('font', {'size': 11})).get('size', 11)})
                    else:
                        ax.set_ylabel(
                            ', '.join([col.name for col in y_cols]), fontdict={'family': graph.get('ylabelfont', graph.get('font', {'family': 'serif'})).get('family', 'serif'), 'size': graph.get('ylabelfont', graph.get('font', {'size': 11})).get('size', 11)})
                    if "title" in graph and graph["title"]:
                        ax.set_title(graph["title"],
                                     fontdict={'family': graph.get('titlefont', graph.get('font', {'family': 'serif'})).get('family', 'serif'), 'size': graph.get('titlefont', graph.get('font', {'size': 12})).get('size', 12)})
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
                        fig.savefig(file_out, dpi=graph.get(
                            "dpi", 300), format=fmt)

        # Place beep at the very end, after all processing and exporting
    if yamldata.get('ping', False):
        try:
            import beepy
            beepy.beep()
        except ImportError:
            print('beepy not installed, cannot beep')


def prune_close_positions(positions, threshold, x_scale=1.0, y_scale=1.0):
    # Given a list of positions [(x1, y1), (x2, y2), ...] and a threshold,
    # return a new list where positions closer than the threshold are grouped and replaced by their centroid.

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
    # print(f"Pruning positions: {positions} -> {used}")
    return used


if __name__ == "__main__":
    main()
