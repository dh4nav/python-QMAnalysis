import pytest
import pandas as pd
import tempfile
from pathlib import Path
from qmanalysis.xyzreader import XYZFile
from qmanalysis.containers import AtomData, FrameData


@pytest.fixture
def xyz_file_content():
    return """3
benzene molecule
C 0.0 0.0 0.0 CA 0.0
H 1.0 0.0 0.0 HA 0.1
O 0.0 1.0 0.0 OA -0.2
"""


def test_xyzfile_reads_atoms_and_frame(tmp_path, xyz_file_content):
    file_path = tmp_path / "benzene.xyz"
    file_path.write_text(xyz_file_content)
    atom_data = AtomData()
    frame_data = FrameData()
    xyz = XYZFile(atom_data, frame_data, str(file_path),
                  file_name="benz", timestep_name="init")
    # Check frame_data
    idx = ("benz", str(file_path), "init")
    assert idx in frame_data.dataframe.index
    assert frame_data.dataframe.loc[idx, "file_comment"] == "benzene molecule"
    # Check atom_data
    for i, (element, alias, charge) in enumerate([("C", "CA", 0.0), ("H", "HA", 0.1), ("O", "OA", -0.2)]):
        atom_idx = ("benz", str(file_path), "init", i)
        assert atom_idx in atom_data.dataframe.index
        assert atom_data.dataframe.loc[atom_idx, "element"] == element
        assert atom_data.dataframe.loc[atom_idx, "alias"] == alias
        assert atom_data.dataframe.loc[atom_idx, "charge"] == charge


def test_xyzfile_missing_atom_lines(tmp_path):
    content = """2
comment
C 0.0 0.0 0.0
"""
    file_path = tmp_path / "bad.xyz"
    file_path.write_text(content)
    atom_data = AtomData()
    frame_data = FrameData()
    with pytest.raises(IndexError, match="Expected 2 atom lines, but got 1"):
        XYZFile(atom_data, frame_data, str(file_path),
                file_name="bad", timestep_name="init")


def test_xyzfile_non_integer_atom_count(tmp_path):
    content = """notanint
comment
C 0.0 0.0 0.0
H 1.0 0.0 0.0
"""
    file_path = tmp_path / "bad.xyz"
    file_path.write_text(content)
    atom_data = AtomData()
    frame_data = FrameData()
    with pytest.raises(ValueError, match="First line must be an integer"):
        XYZFile(atom_data, frame_data, str(file_path),
                file_name="bad", timestep_name="init")


def test_xyzfile_malformed_atom_line(tmp_path):
    content = """2
comment
C 0.0 0.0
H 1.0 0.0 0.0
"""
    file_path = tmp_path / "bad.xyz"
    file_path.write_text(content)
    atom_data = AtomData()
    frame_data = FrameData()
    with pytest.raises(ValueError, match="Malformed atom line"):
        XYZFile(atom_data, frame_data, str(file_path),
                file_name="bad", timestep_name="init")


def test_xyzfile_non_numeric_coordinates(tmp_path):
    content = """1
comment
C x y z
"""
    file_path = tmp_path / "bad.xyz"
    file_path.write_text(content)
    atom_data = AtomData()
    frame_data = FrameData()
    with pytest.raises(ValueError, match="Coordinates must be numeric"):
        XYZFile(atom_data, frame_data, str(file_path),
                file_name="bad", timestep_name="init")


def test_xyzfile_default_alias_and_charge(tmp_path):
    content = """1
comment
C 0.0 0.0 0.0
"""
    file_path = tmp_path / "one.xyz"
    file_path.write_text(content)
    atom_data = AtomData()
    frame_data = FrameData()
    XYZFile(atom_data, frame_data, str(file_path),
            file_name="one", timestep_name="init")
    idx = ("one", str(file_path), "init", 0)
    assert atom_data.dataframe.loc[idx, "alias"] == "0"
    assert pd.isna(atom_data.dataframe.loc[idx, "charge"])
