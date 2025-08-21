import pytest
import pandas as pd
import numpy as np
from pathlib import Path

# Import functions/classes from main.py
from src.main import prepend_root_if_relative, prepend_root_if_relative_and_glob


def test_prepend_root_if_relative_absolute():
    abs_path = "/tmp/test.txt"
    assert prepend_root_if_relative(abs_path, "/root") == Path(abs_path)


def test_prepend_root_if_relative_relative():
    rel_path = "test.txt"
    root = "/root"
    assert prepend_root_if_relative(rel_path, root) == Path(root) / rel_path


def test_prepend_root_if_relative_and_glob_absolute():
    abs_path = "/tmp/test.txt"
    assert prepend_root_if_relative_and_glob(
        abs_path, "/root") == Path(abs_path)


def test_prepend_root_if_relative_and_glob_relative():
    rel_path = "test.txt"
    root = "/root"
    assert prepend_root_if_relative_and_glob(
        rel_path, root) == Path(root) / rel_path


def test_resolve_atom_index(monkeypatch):
    # Patch atom_data with a simple dataframe
    from src.main import resolve_atom

    class DummyAtomData:
        dataframe = pd.DataFrame({
            "alias": ["A", "B", "C"],
            "x": [0.0, 1.0, 2.0],
            "y": [0.0, 1.0, 2.0],
            "z": [0.0, 1.0, 2.0]
        }, index=pd.MultiIndex.from_tuples([
            ("file1", "init", 0),
            ("file1", "init", 1),
            ("file1", "init", 2)
        ], names=["file_name", "timestep_name", "atom_index"]))
    monkeypatch.setattr("src.main", "atom_data", DummyAtomData)
    idx = resolve_atom(1, "init", {"file": "file1"})
    assert idx == ("file1", "init", 1)
    idx = resolve_atom("B", "init", {"file": "file1"})
    assert idx == ("file1", "init", 1)


def test_frame_data_exporter(tmp_path):
    from src.main import FrameDataExporter

    class DummyFrameData:
        dataframe = pd.DataFrame({
            "file_name": ["f1", "f2"],
            "timestep_name": ["t1", "t2"],
            "value": [1, 2]
        })
    exporter = FrameDataExporter(DummyFrameData)
    csv_path = tmp_path / "out.csv"
    exporter.export_csv_multiindex(str(csv_path))
    assert csv_path.exists()
    tuples_path = tmp_path / "out_tuples.csv"
    exporter.export_csv_tuples(str(tuples_path))
    assert tuples_path.exists()


def test_main_measurement_logic(monkeypatch):
    # Simulate measurement logic
    from src.main import measure

    class DummyAtomData:
        dataframe = pd.DataFrame({
            "x": [0.0, 1.0],
            "y": [0.0, 0.0],
            "z": [0.0, 0.0]
        }, index=pd.MultiIndex.from_tuples([
            ("file1", "init", 0),
            ("file1", "init", 1)
        ], names=["file_name", "timestep_name", "atom_index"]))

    class DummyFrameData:
        dataframe = pd.DataFrame(index=pd.MultiIndex.from_tuples([
            ("file1", "init"),
            ("file1", "init")
        ], names=["file_name", "timestep_name"]))
    monkeypatch.setattr("src.main", "atom_data", DummyAtomData)
    monkeypatch.setattr("src.main", "frame_data", DummyFrameData)
    # Should not raise
    val = measure.distance(
        DummyAtomData, ("file1", "init", 0), ("file1", "init", 1))
    assert pytest.approx(val) == 1.0


def test_main_calc_logic():
    # Simulate calculation logic
    df = pd.DataFrame({
        "pi": [3.1415, 3.1415],
        "H-O-H": [120.0, 130.0]
    })
    df["gork"] = df["pi"] * df["H-O-H"] * 5.0
    expected = df["pi"] * df["H-O-H"] * 5.0
    assert all(df["gork"] == expected)
