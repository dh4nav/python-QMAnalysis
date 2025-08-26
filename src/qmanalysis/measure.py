import numpy as np
import pandas as pd
from itertools import combinations


class Measure:

    def _checkEqual(self, dataframe):
        for comb in combinations(range(len(dataframe.index)), 2):
            if dataframe.iloc[comb[0]].equals(dataframe.iloc[comb[1]]):
                raise ValueError("Atom coordinates for atom " +
                                 str(comb[0]) + " and atom " + str(comb[1]) + " can not be equal")

    def _getVector(self, dataframe, norm=False):
        vector1 = dataframe.iloc[1] - dataframe.iloc[0]
        if norm:
            vector1 = vector1/np.linalg.norm(vector1)
        return vector1

    def distance(self, atom_data, atom_index1, atom_index2):
        print(np.linalg.norm((atom_data.dataframe.loc[atom_index1, [
              "x", "y", "z"]] - atom_data.dataframe.loc[atom_index2, ["x", "y", "z"]]).to_numpy()))
        return np.linalg.norm((atom_data.dataframe.loc[atom_index1, ["x", "y", "z"]] - atom_data.dataframe.loc[atom_index2, ["x", "y", "z"]]).to_numpy())

    def angle(self, atom_data, atom_index1, atom_index2, atom_index3):
        self._checkEqual(atom_data.dataframe.loc[[
                         atom_index1, atom_index2, atom_index3], ["x", "y", "z"]])

        vector1 = self._getVector(
            atom_data.dataframe.loc[[atom_index2, atom_index1], ["x", "y", "z"]], norm=True)
        vector2 = self._getVector(
            atom_data.dataframe.loc[[atom_index2, atom_index3], ["x", "y", "z"]], norm=True)

        dotprod = np.dot(vector1, vector2)
        angle_rad = np.arccos(dotprod)
        angle_deg = np.degrees(angle_rad)
        print(angle_deg)
        return angle_deg

    def plane_normal(self, atom_data, atom_index1, atom_index2, atom_index3):

        self._checkEqual(atom_data.dataframe.loc[[
                         atom_index1, atom_index2, atom_index3], ["x", "y", "z"]])

        vector1 = self._getVector(
            atom_data.dataframe.loc[[atom_index2, atom_index1], ["x", "y", "z"]], norm=True)
        vector2 = self._getVector(
            atom_data.dataframe.loc[[atom_index2, atom_index3], ["x", "y", "z"]], norm=True)

        normal = np.cross(vector2, vector1)

        return (normal/np.linalg.norm(normal))

    def dihedral(self, atom_data, atom_index1, atom_index2, atom_index3, atom_index4):

        self._checkEqual(atom_data.dataframe.loc[[
                         atom_index1, atom_index2, atom_index3, atom_index4], ["x", "y", "z"]])

        plane_normal1 = self.plane_normal(
            atom_data, atom_index1, atom_index2, atom_index3)
        plane_normal2 = self.plane_normal(
            atom_data, atom_index2, atom_index3, atom_index4)

        dotprod = np.dot(plane_normal1, plane_normal2)
        angle_rad = np.arccos(dotprod)
        angle_deg = np.degrees(angle_rad)

        return angle_deg
