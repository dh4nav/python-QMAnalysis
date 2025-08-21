import pytest
from qmanalysis.yamlreader import YAMLFile


def test_load_string_and_get_data():
    yaml_str = """
name: Test Procedure
comment: Just a test
version: 1
files:
  - path: "constants.csv"
    type: global_constants_csv
  - path: "mol.xyz"
    type: xyz
    name: testmol
    timestep: "init"
measurements:
  distance:
    - name: bond1
      a: 1
      b: 2
      description: "Test bond"
  angle:
    - name: angle1
      a: 1
      b: 2
      c: 3
substitutions:
  - name: S1
    entries:
      - file: testmol
        atom_index: 1
output:
  - file:
      - path: ./out.csv
        type: csv
  - graph:
      - type: scatter_plot
        x: bond1
        y: angle1
        file: ./plot.tiff
        file_format: tiff
        dpi: 300
"""
    yamlfile = YAMLFile()
    yamlfile.load_string(yaml_str)
    data = yamlfile.get_data()
    assert data["name"] == "Test Procedure"
    assert data["version"] == 1
    assert isinstance(data["files"], list)
    assert "measurements" in data
    assert "distance" in data["measurements"]
    assert data["measurements"]["distance"][0]["name"] == "bond1"
    assert "output" in data
    assert isinstance(data["output"], list)


def test_load_file(tmp_path):
    yaml_content = """
name: FileTest
comment: File test
version: 2
files:
  - path: "constants.csv"
    type: global_constants_csv
"""
    yaml_path = tmp_path / "test.yaml"
    yaml_path.write_text(yaml_content)
    yamlfile = YAMLFile()
    yamlfile.load_file(str(yaml_path))
    data = yamlfile.get_data()
    assert data["name"] == "FileTest"
    assert data["version"] == 2
    assert isinstance(data["files"], list)


def test_schema_accepts_description():
    yaml_str = """
name: DescTest
comment: Description test
version: 1
files:
  - path: "constants.csv"
    type: global_constants_csv
measurements:
  distance:
    - name: bond2
      a: 1
      b: 2
      description: "A bond description"
"""
    yamlfile = YAMLFile()
    yamlfile.load_string(yaml_str)
    data = yamlfile.get_data()
    assert data["measurements"]["distance"][0]["description"] == "A bond description"


def test_str_method():
    yaml_str = """
name: StrTest
comment: Str test
version: 1
files:
  - path: "constants.csv"
    type: global_constants_csv
"""
    yamlfile = YAMLFile()
    yamlfile.load_string(yaml_str)
    s = str(yamlfile)
    assert "StrTest" in s or "strtest" in s.lower()


def test_schema_missing_required_fields():
    yaml_str = """
comment: Missing name
version: 1
files:
  - path: "constants.csv"
    type: global_constants_csv
"""
    yamlfile = YAMLFile()
    with pytest.raises(Exception):
        yamlfile.load_string(yaml_str)
