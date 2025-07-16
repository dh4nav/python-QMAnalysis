import pytest
from qmanalysis.yamlreader import YAMLFile  # Passe den Import an dein Projekt an

valid_yaml = """
name: test_project
comment: This is a test
version: 1

files:
  - path: "data/file1.xyz"
    type: "xyz"

graphics:
  - plot:
      name: "main_plot"
      identifier: "S1"
"""

yaml_with_measurements = """
name: test
comment: test comment
version: 1

files:
  - path: "example.xyz"
    type: "xyz"

graphics:
  - plot:
      name: "plot1"
      identifier: 1

measurements:
  distance:
    - name: "AB"
      a: "S1"
      b: 2
  angle:
    - name: "ABC"
      a: 1
      b: 2
      c: "S3"
  dihedral:
    - name: "ABCD"
      a: "X"
      b: 2
      c: "Y"
      d: 4
"""

invalid_yaml_missing_required = """
comment: This is missing name and version

files:
  - path: "data/file1.xyz"
    type: "xyz"

graphics:
  - plot:
      name: "main_plot"
      identifier: "S1"
"""

def test_load_string_valid_yaml():
    yml = YAMLFile().load_string(valid_yaml)
    data = yml.get_data()
    assert data["name"] == "test_project"
    assert data["files"][0]["type"] == "xyz"
    assert data["graphics"][0]["plot"]["identifier"] == "S1"

def test_load_string_with_measurements():
    yml = YAMLFile().load_string(yaml_with_measurements)
    data = yml.get_data()
    assert "measurements" in data
    assert "distance" in data["measurements"]
    assert data["measurements"]["distance"][0]["name"] == "AB"
    assert data["measurements"]["angle"][0]["c"] == "S3"

def test_str_method():
    yml = YAMLFile().load_string(valid_yaml)
    assert str(yml) == valid_yaml

def test_invalid_yaml_raises():
    with pytest.raises(Exception):
        YAMLFile().load_string(invalid_yaml_missing_required)

def test_data_return_type():
    yml = YAMLFile().load_string(valid_yaml)
    data = yml.get_data()
    assert isinstance(data, dict)


# import unittest
# from collections import OrderedDict
# import tempfile
# import os
# import strictyaml as sy

# # Assuming the YAMLFile class is defined in yamlreader.py
# from qmanalysis.yamlreader import YAMLFile


# class TestYAMLFile(unittest.TestCase):

#     def setUp(self):
#         # Create valid yaml content
#         self.valid_yaml = testyaml="""name: Procedure Name
# comment: XYZ
# version: 111
# substitutions:
#   - S1:
#     - 
#       file: abc/cde.xyz
#       atom: 25 
#     - 
#       file: abab/fff.xyz
#       atom: 19
# sequences:
#   - 1:
#     - abcd
#     - xyza.bla
# measurements:
#   distance:
#     - name: O-H bond
#       a: 25
#       b: S1
#   angle:
#     - name: H-O-H
#       a: S2
#       b: 3
#       c: 4
#     - name: H-O-N
#       a: S2
#       b: S5
#       c: 4
#   dihedral:
#     - name: torsion1
#       a: 1
#       b: 2
#       c: 3
#       d: S3
# """

#         self.bad_yaml = testyaml="""name: Procedure Name
# comment: XYZ
# version: 111
# substitutions:
#   - S1:
#     - 
#       file: abc/cde.xyz
#       atom: 25 
#     - 
#       file: abab/fff.xyz
#       atom: 19
# sequences:
#   - 1:
#     - abcd
#     - xyza.bla
# measurements:
#   -distance:
#     - name: O-H bond
#       a: 25
#       b: S1
#   angle:
#     - name: H-O-H
#       a: S2
#       b: 3
#       c: 4
#     - name: H-O-N
#       a: S2
#       b: S5
#       c: 4
#   dihedral:
#     - name: torsion1
# a: 1
#       b: 2
#       c: 3
#       d: S3
# """
#         self.tempfile = tempfile.NamedTemporaryFile(delete=False, mode='w+', suffix='.yaml')
#         self.tempfile.write(self.valid_yaml)
#         self.tempfile.close()

#     def tearDown(self):
#         os.remove(self.tempfile.name)

#     def test_str_file(self):
#         yaml = YAMLFile()
#         yaml.load_file(self.tempfile.name)
#         self.assertEqual(str(yaml), self.valid_yaml)

#     def test_str_string(self):
#         yaml = YAMLFile()
#         yaml.load_string(self.valid_yaml)
#         self.assertEqual(str(yaml), self.valid_yaml)

#     def test_parsed_file(self):
#         yaml = YAMLFile()
#         yaml.load_file(self.tempfile.name)
#         expected = dict({'name': 'Procedure Name', 'comment': 'XYZ', 'version': 111, 'substitutions': [{'S1': [{'file': 'abc/cde.xyz', 'atom': 25}, {'file': 'abab/fff.xyz', 'atom': 19}]}], 'sequences': [{'1': ['abcd', 'xyza.bla']}], 'measurements': {'distance': [{'name': 'O-H bond', 'a': '25', 'b': 'S1'}], 'angle': [{'name': 'H-O-H', 'a': 'S2', 'b': '3', 'c': '4'}, {'name': 'H-O-N', 'a': 'S2', 'b': 'S5', 'c': '4'}], 'dihedral': [{'name': 'torsion1', 'a': '1', 'b': '2', 'c': '3', 'd': 'S3'}]}})
        
#         self.assertEqual(yaml.get_data(), expected)

#     def test_parsed_string(self):
#         yaml = YAMLFile()
#         yaml.load_string(self.valid_yaml)
#         expected = dict({'name': 'Procedure Name', 'comment': 'XYZ', 'version': 111, 'substitutions': [{'S1': [{'file': 'abc/cde.xyz', 'atom': 25}, {'file': 'abab/fff.xyz', 'atom': 19}]}], 'sequences': [{'1': ['abcd', 'xyza.bla']}], 'measurements': {'distance': [{'name': 'O-H bond', 'a': '25', 'b': 'S1'}], 'angle': [{'name': 'H-O-H', 'a': 'S2', 'b': '3', 'c': '4'}, {'name': 'H-O-N', 'a': 'S2', 'b': 'S5', 'c': '4'}], 'dihedral': [{'name': 'torsion1', 'a': '1', 'b': '2', 'c': '3', 'd': 'S3'}]}})
        
#         self.assertEqual(yaml.get_data(), expected)


#     def test_invalid_yaml_file(self):
#         with tempfile.NamedTemporaryFile(delete=False, mode='w+') as f:
#             f.write(self.bad_yaml)
#             f.close()
#             with self.assertRaises(sy.ruamel.scanner.ScannerError):
#                 yaml = YAMLFile()
#                 yaml.load_file(f.name)
#             os.remove(f.name)

#     def test_invalid_yaml_string(self):
#         with self.assertRaises(sy.ruamel.scanner.ScannerError):
#             yaml = YAMLFile()
#             yaml.load_string(self.bad_yaml)


# if __name__ == '__main__':
#     unittest.main()
