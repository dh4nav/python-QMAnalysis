import pytest
import pandas as pd
from qmanalysis.containers import AtomData, FrameData, MeasurementData


def test_atomdata_init_empty():
    atom = AtomData()
    # Check index names
    assert atom.dataframe.index.names == [
        'file_name', 'file_path', 'timestep_name', 'atom_index']
    # Check columns
    expected_cols = ["element", "alias", "charge", "x", "y", "z"]
    assert list(atom.dataframe.columns) == expected_cols
    # Should be empty
    assert atom.dataframe.empty


def test_atomdata_add_row():
    atom = AtomData()
    idx = ('file1', '/path/to/file1', 'init', 1)
    atom.dataframe.loc[idx] = ["C", "CA", 0, 0.0, 1.0, 2.0]
    assert atom.dataframe.loc[idx, "element"] == "C"
    assert atom.dataframe.loc[idx, "alias"] == "CA"
    assert atom.dataframe.loc[idx, "x"] == 0.0


def test_framedata_init_empty():
    frame = FrameData()
    assert frame.dataframe.index.names == [
        'file_name', 'file_path', 'timestep_name']
    expected_cols = ["raw_data", "energy", "zero-point energy", "file_comment"]
    assert list(frame.dataframe.columns) == expected_cols
    assert frame.dataframe.empty


def test_framedata_add_row():
    frame = FrameData()
    idx = ('file1', '/path/to/file1', 'init')
    frame.dataframe.loc[idx] = ["raw", 1.23, 0.45, "comment"]
    assert frame.dataframe.loc[idx, "energy"] == 1.23
    assert frame.dataframe.loc[idx, "file_comment"] == "comment"


def test_measurementdata_init_empty():
    meas = MeasurementData()
    assert meas.dataframe.index.names == [
        'file_name', 'file_path', 'timestep_name']
    assert meas.dataframe.empty
    # Should have no columns
    assert len(meas.dataframe.columns) == 0


def test_measurementdata_add_row_and_column():
    meas = MeasurementData()
    idx = ('file1', '/path/to/file1', 'init')
    meas.dataframe.loc[idx, "distance"] = 1.5
    meas.dataframe.loc[idx, "angle"] = 120.0
    assert meas.dataframe.loc[idx, "distance"] == 1.5
    assert meas.dataframe.loc[idx, "angle"] == 120.0


def test_atomdata_multiple_rows():
    atom = AtomData()
    idx1 = ('file1', '/path/to/file1', 'init', 1)
    idx2 = ('file1', '/path/to/file1', 'init', 2)
    atom.dataframe.loc[idx1] = ["C", "CA", 0, 0.0, 1.0, 2.0]
    atom.dataframe.loc[idx2] = ["H", "HA", 0, 1.0, 2.0, 3.0]
    assert len(atom.dataframe) == 2
    assert set(atom.dataframe["element"]) == {"C", "H"}


def test_framedata_multiple_rows():
    frame = FrameData()
    idx1 = ('file1', '/path/to/file1', 'init')
    idx2 = ('file2', '/path/to/file2', 'ts')
    frame.dataframe.loc[idx1] = ["raw", 1.23, 0.45, "comment"]
    frame.dataframe.loc[idx2] = ["raw2", 2.34, 0.67, "comment2"]
    assert len(frame.dataframe) == 2
    assert set(frame.dataframe["file_comment"]) == {"comment", "comment2"}


def test_measurementdata_multiple_rows_and_columns():
    meas = MeasurementData()
    idx1 = ('file1', '/path/to/file1', 'init')
    idx2 = ('file2', '/path/to/file2', 'ts')
    meas.dataframe.loc[idx1, "distance"] = 1.5
    meas.dataframe.loc[idx2, "distance"] = 2.5
    meas.dataframe.loc[idx1, "angle"] = 120.0
    meas.dataframe.loc[idx2, "angle"] = 130.0
    assert len(meas.dataframe) == 2
    assert set(meas.dataframe["distance"]) == {1.5, 2.5}
    assert set(meas.dataframe["angle"]) == {120.0, 130.0}
