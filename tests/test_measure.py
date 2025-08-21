import pytest
import numpy as np
import pandas as pd
from qmanalysis.measure import Measure


@pytest.fixture
def atom_data():
    class DummyAtomData:
        pass
    atom = DummyAtomData()
    # Create a DataFrame with 4 atoms in 3D space
    atom.dataframe = pd.DataFrame([
        {"x": 0.0, "y": 0.0, "z": 0.0},
        {"x": 1.0, "y": 0.0, "z": 0.0},
        {"x": 1.0, "y": 1.0, "z": 0.0},
        {"x": 1.0, "y": 1.0, "z": 1.0}
    ])
    return atom


def test_distance(atom_data):
    m = Measure()
    dist = m.distance(atom_data, 0, 1)
    assert pytest.approx(dist) == 1.0


def test_angle(atom_data):
    m = Measure()
    angle = m.angle(atom_data, 0, 1, 2)
    assert pytest.approx(angle, abs=0.1) == 90.0


def test_plane_normal(atom_data):
    m = Measure()
    normal = m.plane_normal(atom_data, 0, 1, 2)
    expected = np.array([0.0, 0.0, 1.0])
    np.testing.assert_allclose(normal, expected, atol=1e-6)


def test_dihedral(atom_data):
    m = Measure()
    angle = m.dihedral(atom_data, 0, 1, 2, 3)
    assert pytest.approx(angle, abs=0.1) == 90.0


def test_check_equal_raises():
    m = Measure()
    df = pd.DataFrame([
        {"x": 1.0, "y": 1.0, "z": 1.0},
        {"x": 1.0, "y": 1.0, "z": 1.0}
    ])
    with pytest.raises(ValueError, match="Atom coordinates for atom 0 and atom 1 can not be equal"):
        m._checkEqual(df)


def test_get_vector_norm_false(atom_data):
    m = Measure()
    df = atom_data.dataframe.loc[[0, 1], ["x", "y", "z"]]
    vec = m._getVector(df, norm=False)
    assert np.allclose(vec, np.array([1.0, 0.0, 0.0]))


def test_get_vector_norm_true(atom_data):
    m = Measure()
    df = atom_data.dataframe.loc[[0, 1], ["x", "y", "z"]]
    vec = m._getVector(df, norm=True)
    assert np.allclose(vec, np.array([1.0, 0.0, 0.0]))


def test_angle_colinear():
    m = Measure()

    class DummyAtomData:
        pass
    atom = DummyAtomData()
    atom.dataframe = pd.DataFrame([
        {"x": 0.0, "y": 0.0, "z": 0.0},
        {"x": 1.0, "y": 0.0, "z": 0.0},
        {"x": 2.0, "y": 0.0, "z": 0.0}
    ])
    angle = m.angle(atom, 0, 1, 2)
    assert pytest.approx(angle, abs=0.1) == 180.0


def test_plane_normal_nonplanar():
    m = Measure()

    class DummyAtomData:
        pass
    atom = DummyAtomData()
    atom.dataframe = pd.DataFrame([
        {"x": 0.0, "y": 0.0, "z": 0.0},
        {"x": 1.0, "y": 0.0, "z": 1.0},
        {"x": 0.0, "y": 1.0, "z": 1.0}
    ])
    normal = m.plane_normal(atom, 0, 1, 2)
    assert normal.shape == (3,)
    assert not np.allclose(normal, np.zeros(3))


def test_dihedral_linear():
    m = Measure()

    class DummyAtomData:
        pass
    atom = DummyAtomData()
    atom.dataframe = pd.DataFrame([
        {"x": 0.0, "y": 0.0, "z": 0.0},
        {"x": 1.0, "y": 0.0, "z": 0.0},
        {"x": 2.0, "y": 0.0, "z": 0.0},
        {"x": 3.0, "y": 0.0, "z": 0.0}
    ])
    angle = m.dihedral(atom, 0, 1, 2, 3)
    assert pytest.approx(angle, abs=0.1) == 0.0
