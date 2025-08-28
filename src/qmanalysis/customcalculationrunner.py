# --- Custom Calculation Section ---

import numpy as np
import pandas as pd
import scipy    # Import scipy if needed, otherwise it will be optional


class CustomCalculationRunner:
    """
    Runs user-defined calculations on frame_data using numpy, scipy, and user-supplied values.
    Stores results in new columns as specified in the YAML input.
    """

    def __init__(self, frame_data, extra_globals=None):
        self.frame_data = frame_data
        # Allowed modules/functions for eval
        self.safe_globals = {
            "np": np,
            "numpy": np,
            "pd": pd,
            "max": max,
            "min": min,
            "abs": abs,
            "sum": sum,
            "len": len,
            "round": round,
            "float": float,
            "int": int,
            "str": str
        }

        # Optionally add scipy if needed
        try:
            self.safe_globals["scipy"] = scipy
        except ImportError:
            pass

        def pivot(file_name, col_name):
            # Find the row with the given file_name in the MultiIndex
            # Assumes file_name is in the first level of the MultiIndex
            try:
                # Get all rows with the given file_name
                df = self.frame_data.dataframe.xs(file_name, level="file_name")
                # If multiple rows, you may need to further filter by other index levels
                # For now, just take the first match
                value = df.iloc[0][col_name]
                return value
            except Exception as e:
                print(f"pivot error: {e}")
                return None

        self.safe_globals["pivot"] = pivot

        if extra_globals:
            self.safe_globals.update(extra_globals)

    def run(self, calculations):
        """
        calculations: list of dicts, each with keys:
            - name: output column name
            - expr: string, Python expression using allowed functions and frame_data columns
        """
        for calc in calculations:
            name = calc["name"]
            expr = calc["expr"]
            # Prepare local variables: all columns as Series
            local_vars = {
                col: self.frame_data.dataframe[col] for col in self.frame_data.dataframe.columns}
            # Allow user-supplied values
            if "values" in calc:
                local_vars.update(calc["values"])
            # Evaluate the expression for each row
            try:
                import re

                def safe_eval(row):
                    # Replace column names in expr with their values from the row, using bracket notation
                    expr_safe = expr
                    for col in self.frame_data.dataframe.columns:
                        # Only replace if col is in expr as a variable (not inside quotes)
                        # This regex matches col name not inside single or double quotes
                        def replacer(match):
                            # Check if match is inside quotes
                            before = expr_safe[:match.start()]
                            num_single = before.count("'")
                            num_double = before.count('"')
                            # If odd number of quotes before, we're inside quotes
                            if num_single % 2 == 1 or num_double % 2 == 1:
                                return match.group(0)
                            return f'row["{col}"]'
                        expr_safe = re.sub(
                            rf'\b{re.escape(col)}\b', replacer, expr_safe)
                    try:
                        return eval(expr_safe, self.safe_globals, {"row": row, **local_vars})
                    except Exception as e:
                        print(
                            f"Error in custom calculation '{name}' for row: {e}")
                        return np.nan
                self.frame_data.dataframe[name] = self.frame_data.dataframe.apply(
                    safe_eval, axis=1)
            except Exception as e:
                print(f"Error in custom calculation '{name}': {e}")
                self.frame_data.dataframe[name] = np.nan

# --- Run custom calculations if specified in YAML ---
# Example YAML for two calculations:
# custom_calculations:
#   - name: calc1
#     expr: "distance_O_H_bond * 2"
#   - name: calc2
#     expr: "calc1 + angle_H_O_H"
