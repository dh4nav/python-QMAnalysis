import pytest
import pandas as pd
from qmanalysis.customcalculationrunner import CustomCalculationRunner


@pytest.fixture
def frame_data():
    class DummyFrameData:
        pass
    df = pd.DataFrame({
        "pi": [3.1415, 3.1415],
        "H-O-H": [120.0, 130.0],
        "böp": [5.0, 6.0]
    })
    DummyFrameData.dataframe = df
    return DummyFrameData()


def test_run_simple_expression(frame_data):
    runner = CustomCalculationRunner(frame_data)
    calc = [{"name": "gork", "expr": "pi * `H-O-H` * 5.0"}]
    runner.run(calc)
    expected = frame_data.dataframe["pi"] * frame_data.dataframe["H-O-H"] * 5.0
    assert all(frame_data.dataframe["gork"] == expected)


def test_run_numpy_expression(frame_data):
    runner = CustomCalculationRunner(frame_data)
    calc = [{"name": "sqrt_gork", "expr": "np.sqrt(pi * `H-O-H`)"}]
    runner.run(calc)
    import numpy as np
    expected = np.sqrt(
        frame_data.dataframe["pi"] * frame_data.dataframe["H-O-H"])
    assert all(frame_data.dataframe["sqrt_gork"] == expected)


def test_run_chained_expression(frame_data):
    runner = CustomCalculationRunner(frame_data)
    calc = [
        {"name": "böps", "expr": "pi * `H-O-H`"},
        {"name": "böpn", "expr": "böp * böp"}
    ]
    runner.run(calc)
    expected_böps = frame_data.dataframe["pi"] * frame_data.dataframe["H-O-H"]
    expected_böpn = frame_data.dataframe["böp"] * frame_data.dataframe["böp"]
    assert all(frame_data.dataframe["böps"] == expected_böps)
    assert all(frame_data.dataframe["böpn"] == expected_böpn)


def test_run_invalid_expression(frame_data):
    runner = CustomCalculationRunner(frame_data)
    calc = [{"name": "fail", "expr": "not_a_column * 2"}]
    with pytest.raises(Exception):
        runner.run(calc)
