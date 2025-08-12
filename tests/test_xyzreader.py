import pytest
import pandas as pd
from pathlib import Path
from qmanalysis.xyzreader import XYZFile  # Adjust this import as needed

# --- Fixtures ---

@pytest.fixture
def dummy_timestep_data():
    class DummyData:
        def __init__(self):
            self.dataframe = pd.DataFrame(columns=[
                "file_name", "file_path", "timestep_time", "timestep_index",
                "raw_data", "energy", "zero-point energy", "file_comment", "measurements"
            ])
    return DummyData()

@pytest.fixture
def dummy_atom_data():
    class DummyAtomData:
        def __init__(self):
            self.dataframe = pd.DataFrame(columns=[
                "file_name", "file_path", "file_index", "atom_index",
                "element", "x", "y", "z", "alias", "charge",
                "timestep_time", "timestep_index"
            ])
    return DummyAtomData()

@pytest.fixture
def xyz_file_content():
    return """3
Comment line
H 0.0 0.0 0.0
C 1.0 0.0 0.0
O 0.0 1.0 0.0
"""

# --- Tests ---

def test_xyzfile_parsing(tmp_path, xyz_file_content, dummy_atom_data, dummy_timestep_data):
    file_path = tmp_path / "test.xyz"
    file_path.write_text(xyz_file_content)

    XYZFile(dummy_atom_data, dummy_timestep_data, str(file_path))

    # Check timestep entry
    assert len(dummy_timestep_data.dataframe) == 1
    timestep_row = dummy_timestep_data.dataframe.iloc[0]
    assert timestep_row["file_path"] == str(file_path)
    assert timestep_row["file_comment"] == "Comment line"
    assert timestep_row["measurements"] == {}

    # Check atom entries
    df = dummy_atom_data.dataframe
    assert len(df) == 3
    assert list(df["element"]) == ["H", "C", "O"]
    assert df.iloc[1]["x"] == 1.0
    assert df.iloc[0]["alias"] == str(df.iloc[0]["atom_index"])
    assert pd.isna(df.iloc[0]["charge"])

def test_optional_alias_and_charge(tmp_path, dummy_atom_data, dummy_timestep_data):
    content = """2
Water molecule
H 0.0 0.0 0.0 H1 0.25
O 1.0 0.0 0.0 O1 -0.5
"""
    file_path = tmp_path / "alias_charge.xyz"
    file_path.write_text(content)

    XYZFile(dummy_atom_data, dummy_timestep_data, str(file_path))

    df = dummy_atom_data.dataframe
    assert df.loc[0, "alias"] == "H1"
    assert df.loc[0, "charge"] == 0.25
    assert df.loc[1, "alias"] == "O1"
    assert df.loc[1, "charge"] == -0.5

def test_wrong_atom_count_too_low(tmp_path, dummy_atom_data, dummy_timestep_data):
    content = """4
Comment line
H 0.0 0.0 0.0
C 1.0 0.0 0.0
O 0.0 1.0 0.0
"""
    file_path = tmp_path / "wrong_count.xyz"
    file_path.write_text(content)

    with pytest.raises(IndexError, match="Expected 4 atom lines, but got 3"):
        XYZFile(dummy_atom_data, dummy_timestep_data, str(file_path))

def test_wrong_atom_count_too_high(tmp_path, dummy_atom_data, dummy_timestep_data):
    content = """2
Comment line
H 0.0 0.0 0.0
C 1.0 0.0 0.0
O 0.0 1.0 0.0
"""
    file_path = tmp_path / "wrong_count.xyz"
    file_path.write_text(content)

    with pytest.raises(IndexError, match="Expected 2 atom lines, but got 3"):
        XYZFile(dummy_atom_data, dummy_timestep_data, str(file_path))

def test_malformed_atom_line(tmp_path, dummy_atom_data, dummy_timestep_data):
    content = """2
Bad line
H 0.0 0.0
C 1.0 0.0 0.0
"""
    file_path = tmp_path / "bad_line.xyz"
    file_path.write_text(content)

    with pytest.raises(ValueError, match="Malformed atom line 3"):
        XYZFile(dummy_atom_data, dummy_timestep_data, str(file_path))

def test_non_numeric_coordinates(tmp_path, dummy_atom_data, dummy_timestep_data):
    content = """1
Non-numeric
H A B C
"""
    file_path = tmp_path / "non_numeric.xyz"
    file_path.write_text(content)

    with pytest.raises(ValueError, match="Coordinates must be numeric in line 3"):
        XYZFile(dummy_atom_data, dummy_timestep_data, str(file_path))

def test_missing_header_line(tmp_path, dummy_atom_data, dummy_timestep_data):
    content = "2\n"
    file_path = tmp_path / "missing_comment.xyz"
    file_path.write_text(content)

    with pytest.raises(ValueError, match="File too short"):
        XYZFile(dummy_atom_data, dummy_timestep_data, str(file_path))

def test_invalid_atom_count_header(tmp_path, dummy_atom_data, dummy_timestep_data):
    content = "Two\nComment\nH 0.0 0.0 0.0"
    file_path = tmp_path / "bad_header.xyz"
    file_path.write_text(content)

    with pytest.raises(ValueError, match="First line must be an integer"):
        XYZFile(dummy_atom_data, dummy_timestep_data, str(file_path))
