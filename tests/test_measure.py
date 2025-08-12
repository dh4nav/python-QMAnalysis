import unittest
import numpy as np
import pandas as pd
from qmanalysis.measure import Measure  # replace with actual module name

class AtomData:
    def __init__(self, dataframe):
        self.dataframe = dataframe

class TestMeasure(unittest.TestCase):

    def setUp(self):
        self.measure = Measure()
        self.data = pd.DataFrame([
            {"x": 0, "y": 0, "z": 0},   # Atom 0
            {"x": 1, "y": 0, "z": 0},   # Atom 1
            {"x": 1, "y": 1, "z": 0},   # Atom 2
            {"x": 1, "y": 1, "z": 1}    # Atom 3
        ])
        self.atom_data = AtomData(self.data)

    def test_distance(self):
        dist = self.measure.distance(self.atom_data, 0, 1)
        self.assertAlmostEqual(dist, 1.0)

    def test_angle_90_deg(self):
        angle = self.measure.angle(self.atom_data, 0, 1, 2)
        self.assertAlmostEqual(angle, 90.0)

    def test_plane_normal(self):
        normal = self.measure.plane_normal(self.atom_data, 0, 1, 2)
        expected = np.array([0, 0, 1])
        np.testing.assert_array_almost_equal(normal, expected)

    def test_dihedral_90_deg(self):
        dih = self.measure.dihedral(self.atom_data, 0, 1, 2, 3)
        self.assertAlmostEqual(dih, 90.0)

    def test_checkEqual_raises(self):
        duplicate_data = pd.DataFrame([
            {"x": 0, "y": 0, "z": 0},
            {"x": 0, "y": 0, "z": 0}
        ])
        with self.assertRaises(ValueError):
            self.measure._checkEqual(duplicate_data)

    def test_getVector_normalized(self):
        vector = self.measure._getVector(self.data.loc[[0, 1], ["x", "y", "z"]], norm=True)
        expected = np.array([1.0, 0.0, 0.0])
        np.testing.assert_array_almost_equal(vector, expected)

    def test_getVector_non_normalized(self):
        vector = self.measure._getVector(self.data.loc[[0, 1], ["x", "y", "z"]], norm=False)
        expected = np.array([1.0, 0.0, 0.0])
        np.testing.assert_array_almost_equal(vector, expected)

if __name__ == '__main__':
    unittest.main()


# import unittest
# from collections import OrderedDict
# import tempfile
# import os
# import numpy as np

# # Assuming the XYZFile class is defined in xyzreader.py
# from qmanalysis.measure import Measure


# class TestMeasure(unittest.TestCase):

#     def test_distance(self):
#         measure = Measure()
#         self.assertAlmostEqual(measure.distance([1.0,1.0,1.0],[0.0,1.0,1.0]), 1.0)
#         self.assertAlmostEqual(measure.distance([1.0,1.0,1.0],[1.0,1.0,1.0]), 0.0)
#         self.assertAlmostEqual(measure.distance([0.0,0.0,0.0],[0.0,0.0,0.0]), 0.0)

#     def test_angle(self):
#         measure = Measure()
#         self.assertAlmostEqual(measure.angle([1.0,0.0,0.0],[0.0,0.0,0.0],[1.0,0.0,1.0]), 45.0)
#         self.assertAlmostEqual(measure.angle([0.0,0.0,1.0],[0.0,0.0,0.0],[0.0,0.0,2.0]), 0.0)
#         self.assertAlmostEqual(measure.angle([1.0,0.0,0.0],[0.0,0.0,0.0],[-1.0,0.0,0.0]), 180.0)
#         with self.assertRaises(ValueError):
#             measure.angle([0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0])
#         with self.assertRaises(ValueError):
#             measure.angle([2.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0])
#         with self.assertRaises(ValueError):
#             measure.angle([0.0,0.0,0.0],[0.0,2.0,0.0],[0.0,0.0,0.0])
#         with self.assertRaises(ValueError):
#             measure.angle([0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,2.0])

#     def test_plane_normal(self):
#         measure = Measure()
#         self.assertTrue(np.allclose(measure.plane_normal([1.0,0.0,0.0],[0.0,0.0,0.0],[0.0,1.0,0.0]), np.array([0.0,0.0,1.0])))
#         self.assertTrue(np.allclose(measure.plane_normal([2.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.5,0.0]), np.array([0.0,0.0,1.0])))
#         with self.assertRaises(ValueError):
#             measure.plane_normal([0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0])
#         with self.assertRaises(ValueError):
#             measure.plane_normal([2.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0])
#         with self.assertRaises(ValueError):
#             measure.plane_normal([0.0,0.0,0.0],[0.0,2.0,0.0],[0.0,0.0,0.0])
#         with self.assertRaises(ValueError):
#             measure.plane_normal([0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,2.0])

#     def test_dihedral(self):
#         measure = Measure()
#         self.assertAlmostEqual(measure.dihedral([1.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,1.0],[1.0,0.0,1.0]), 0.0)
#         self.assertAlmostEqual(measure.dihedral([2.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,2.0],[1.0,1.0,2.0]), 45.0)
#         self.assertAlmostEqual(measure.dihedral([1.0,0.0,0.0],[0.0,0.0,0.0],[0.0,1.0,0.0],[-1.0,1.0,1.0]), 135.0)
#         with self.assertRaises(ValueError):
#             measure.dihedral([0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0])
#         with self.assertRaises(ValueError):
#             measure.dihedral([2.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0])
#         with self.assertRaises(ValueError):
#             measure.dihedral([0.0,0.0,0.0],[0.0,2.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0])
#         with self.assertRaises(ValueError):
#             measure.dihedral([0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,2.0],[0.0,0.0,0.0])
#         with self.assertRaises(ValueError):
#             measure.dihedral([0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,1.0,0.0])
#         with self.assertRaises(ValueError):
#             measure.dihedral([0.0,2.0,0.0],[0.0,2.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0])
#         with self.assertRaises(ValueError):
#             measure.dihedral([0.0,2.0,0.0],[0.0,0.0,0.0],[0.0,2.0,0.0],[0.0,0.0,0.0])
#         with self.assertRaises(ValueError):
#             measure.dihedral([0.0,2.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,2.0,0.0]) 
#         with self.assertRaises(ValueError):
#             measure.dihedral([0.0,0.0,0.0],[0.0,0.0,1.0],[0.0,0.0,1.0],[0.0,0.0,0.0]) 
#         with self.assertRaises(ValueError):
#             measure.dihedral([0.0,0.0,0.0],[0.0,0.0,1.0],[0.0,0.0,0.0],[0.0,0.0,1.0])
#         with self.assertRaises(ValueError):
#             measure.dihedral([0.0,0.0,0.0],[0.0,0.0,0.0],[1.0,0.0,0.0],[1.0,0.0,0.0])

# if __name__ == '__main__':
#     unittest.main()
