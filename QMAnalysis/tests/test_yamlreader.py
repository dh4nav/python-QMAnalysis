import pytest
from qmanalysis.yamlreader import YAMLFile  # adjust to your module path
import strictyaml as sy

# --- Valid test YAML strings ---

valid_yaml = """
name: test_project
comment: This is a test
version: 1

files:
  - path: "data/file1.xyz"
    type: "xyz"

output:
  - graphics:
      - plot:
          name: "main_plot"
          identifier: "S1"
    file:
      - path: "output/results.dat"
        type: "text"
"""

yaml_with_measurements = """
name: test
comment: test comment
version: 1

files:
  - path: "example.xyz"
    type: "xyz"

output:
  - graphics:
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

output:
  - graphics:
      - plot:
          name: "main_plot"
          identifier: "S1"
"""

# --- Tests ---

def test_load_string_valid_yaml():
    yml = YAMLFile().load_string(valid_yaml)
    data = yml.get_data()
    assert data["name"] == "test_project"
    assert data["files"][0]["type"] == "xyz"
    assert data["output"][0]["graphics"][0]["plot"]["identifier"] == "S1"

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
    with pytest.raises(sy.YAMLValidationError):
        YAMLFile().load_string(invalid_yaml_missing_required)

def test_data_return_type():
    yml = YAMLFile().load_string(valid_yaml)
    data = yml.get_data()
    assert isinstance(data, dict)

def test_loading_from_file(tmp_path):
    file = tmp_path / "test.yaml"
    file.write_text(valid_yaml)
    yml = YAMLFile().load_file(file)
    data = yml.get_data()
    assert data["version"] == 1
    assert str(yml) == valid_yaml

def test_malformed_yaml_raises():
    bad_yaml = "name: test\nfiles: [\n - path: \"bad\""  # missing closing brackets
    with pytest.raises(sy.YAMLError):
        YAMLFile().load_string(bad_yaml)

def test_output_file_section():
    yml = YAMLFile().load_string(valid_yaml)
    output = yml.get_data()["output"][0]
    assert "file" in output
    assert output["file"][0]["path"] == "output/results.dat"
    assert output["file"][0]["type"] == "text"
