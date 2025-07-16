import pytest
import pandas as pd
from pathlib import PurePath
from tempfile import NamedTemporaryFile
from qmanalysis.xyzreader import XYZFile  # replace with actual module name

class DummyData:
    def __init__(self):
        self.dataframe = pd.DataFrame(columns=[
            "file_name", "file_path", "timestep_time", "timestep_index",
            "raw_data", "energy", "zero-point energy", "file_comment", "measurements"
        ])

class DummyAtomData:
    def __init__(self):
        self.dataframe = pd.DataFrame(columns=[
            "file_name", "file_path", "file_index", "atom_index",
            "element", "x", "y", "z", "alias", "charge",
            "timestep_time", "timestep_index"
        ])

@pytest.fixture
def xyz_file_content():
    return """3
Comment line
H 0.0 0.0 0.0
C 1.0 0.0 0.0
O 0.0 1.0 0.0
"""


def test_xyzfile_parsing(xyz_file_content):
    with NamedTemporaryFile(delete=False, mode='w', suffix='.xyz') as f:
        f.write(xyz_file_content)
        file_path = f.name

    atom_data = DummyAtomData()
    timestep_data = DummyData()

    XYZFile(atom_data, timestep_data, file_path)

    # Check timestep entry
    assert len(timestep_data.dataframe) == 1
    timestep_row = timestep_data.dataframe.iloc[0]
    assert timestep_row["file_path"] == file_path
    assert timestep_row["file_comment"] == "Comment line"
    assert timestep_row["measurements"] == {}

    # Check atom entries
    assert len(atom_data.dataframe) == 3
    assert all(atom_data.dataframe["element"] == ["H", "C", "O"])
    assert atom_data.dataframe.iloc[1]["x"] == 1.0


def test_wrong_atom_count():
    content = """2
Comment line
H 0.0 0.0 0.0
C 1.0 0.0 0.0
O 0.0 1.0 0.0
"""  # says 2 atoms, gives 3
    with NamedTemporaryFile(delete=False, mode='w', suffix='.xyz') as f:
        f.write(content)
        file_path = f.name

    atom_data = DummyAtomData()
    timestep_data = DummyData()

    with pytest.raises(IndexError):
        XYZFile(atom_data, timestep_data, file_path)


def test_malformed_atom_line():
    content = """2
Comment line
H 0.0 0.0
C 1.0 0.0 0.0
"""  # first atom line has only 3 fields
    with NamedTemporaryFile(delete=False, mode='w', suffix='.xyz') as f:
        f.write(content)
        file_path = f.name

    atom_data = DummyAtomData()
    timestep_data = DummyData()

    with pytest.raises(ValueError):
        XYZFile(atom_data, timestep_data, file_path)
