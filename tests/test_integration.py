import os
import sys
import tempfile
import shutil
import subprocess
import pytest
from pathlib import Path

TEST_YAML_MINIMAL = """
name: test
comment: Minimal test
version: 1
files:
  - path: test.xyz
    type: xyz
output: []
"""

TEST_XYZ = """
3
comment
H 0.0 0.0 0.0
H 0.0 0.0 1.0
O 0.0 1.0 0.0
"""


def write_file(path, content):
    with open(path, "w") as f:
        f.write(content)


def run_cli(yaml_path, root=None):
    cmd = [sys.executable, "-m", "src.main", yaml_path]
    if root:
        cmd.extend(["-rp", root])
    return subprocess.run(cmd, capture_output=True, text=True)


def test_minimal_xyz(tmp_path):
    yaml_path = tmp_path / "test.yaml"
    xyz_path = tmp_path / "test.xyz"
    write_file(yaml_path, TEST_YAML_MINIMAL)
    write_file(xyz_path, TEST_XYZ)
    result = run_cli(str(yaml_path), root=str(tmp_path))
    assert result.returncode == 0
    assert "error" not in result.stderr.lower()


def test_invalid_yaml(tmp_path):
    yaml_path = tmp_path / "bad.yaml"
    yaml_path.write_text("name: [unclosed_list\n")
    result = run_cli(str(yaml_path), root=str(tmp_path))
    assert result.returncode != 0 or "error" in result.stderr.lower()


def test_ping_option(monkeypatch, tmp_path):
    yaml_path = tmp_path / "test.yaml"
    xyz_path = tmp_path / "test.xyz"
    yaml = TEST_YAML_MINIMAL + "ping: true\n"
    write_file(yaml_path, yaml)
    write_file(xyz_path, TEST_XYZ)
    called = {}
    monkeypatch.setitem(sys.modules, "beepy", type(
        "B", (), {"beep": lambda: called.setdefault("beep", True)})())
    result = run_cli(str(yaml_path), root=str(tmp_path))
    assert result.returncode == 0
    assert called.get("beep")


def test_missing_file(tmp_path):
    yaml_path = tmp_path / "test.yaml"
    yaml = TEST_YAML_MINIMAL.replace("test.xyz", "missing.xyz")
    write_file(yaml_path, yaml)
    result = run_cli(str(yaml_path), root=str(tmp_path))
    assert result.returncode != 0 or "missing" in result.stderr.lower()


def test_invalid_file_type(tmp_path):
    yaml_path = tmp_path / "test.yaml"
    yaml = TEST_YAML_MINIMAL.replace("type: xyz", "type: unknown")
    xyz_path = tmp_path / "test.xyz"
    write_file(yaml_path, yaml)
    write_file(xyz_path, TEST_XYZ)
    result = run_cli(str(yaml_path), root=str(tmp_path))
    assert result.returncode != 0 or "unknown" in result.stderr.lower()

# More integration tests can be added for plotting, globbing, etc.
