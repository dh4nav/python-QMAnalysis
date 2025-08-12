import pytest
import pandas as pd
import numpy as np
from qmanalysis.measure import Measure
from qmanalysis.containers import AtomData, TimestepData, MeasurementData


@pytest.fixture
def mock_atom_data():
    atom = AtomData()
    atom.dataframe = pd.DataFrame([
        {"x": 0.0, "y": 0.0, "z": 0.0},
        {"x": 1.0, "y": 0.0, "z": 0.0},
        {"x": 1.0, "y": 1.0, "z": 0.0},
        {"x": 1.0, "y": 1.0, "z": 1.0}
    ])
    return atom

def test_atomdata_columns():
    atom = AtomData()
    expected = ["file_name", "file_index", "timestep_time", "timestep_index",
                "atom_index", "element", "alias", "charge", "x", "y", "z"]
    assert list(atom.dataframe.columns) == expected

def test_timestepdata_columns():
    timestep = TimestepData()
    expected = ["file_name", "file_path", "timestep_time", "timestep_index",
                "raw_data", "energy", "zero-point energy", "file_comment", "measurements"]
    assert list(timestep.dataframe.columns) == expected

@pytest.fixture
def mock_yaml_data():
    return {
        "measurements": {
            "energy": ["E_total", "E_kinetic"],
            "geometry": ["bond_length", "angle"]
        }
    }

def test_measurementdata(mock_yaml_data):
    names = ["timestep1", "timestep2"]
    # Simulate the logic you commented out in MeasurementData
    expected_index = []
    for mtype in mock_yaml_data["measurements"]:
        for mname in mock_yaml_data["measurements"][mtype]:
            expected_index.append(f"{mtype}-{mname}")

    # Use just names in constructor (like in your class)
    measurement = MeasurementData(names)
    assert list(measurement.dataframe["file_name"]) == names
    assert len(measurement.dataframe) == len(names)

def test_distance(mock_atom_data):
    m = Measure()
    dist = m.distance(mock_atom_data, 0, 1)
    assert pytest.approx(dist) == 1.0

def test_angle(mock_atom_data):
    m = Measure()
    angle = m.angle(mock_atom_data, 0, 1, 2)
    assert pytest.approx(angle, abs=0.1) == 90.0

def test_plane_normal(mock_atom_data):
    m = Measure()
    normal = m.plane_normal(mock_atom_data, 0, 1, 2)
    expected = np.array([0.0, 0.0, 1.0])
    np.testing.assert_allclose(normal, expected, atol=1e-6)

def test_dihedral(mock_atom_data):
    m = Measure()
    angle = m.dihedral(mock_atom_data, 0, 1, 2, 3)
    assert pytest.approx(angle, abs=0.1) == 90.0

def test_check_equal_raises():
    m = Measure()
    df = pd.DataFrame([
        {"x": 1.0, "y": 1.0, "z": 1.0},
        {"x": 1.0, "y": 1.0, "z": 1.0}
    ])
    with pytest.raises(ValueError, match="Atom coordinates for atom 0 and atom 1 can not be equal"):
        m._checkEqual(df)
