import pandas as pd

class AtomData:
    def __init__(self):
        self.dataframe =  pd.DataFrame(columns = ["file_name", "file_index", "timestep_time", "timestep_index", "atom_index", "element", "alias", "charge", "x", "y", "z" ])   

class TimestepData:
    def __init__(self):
        self.dataframe =  pd.DataFrame(columns = ["file_name", "file_path", "timestep_time", "timestep_index", "raw_data", "energy", "zero-point energy", "file_comment", "measurements" ])     

# class MeasurementData:
#     __init__(self):
#         self.dataframe = pd.dataframe(columns = ["file_name"], )


# class Atom:
    
#     def __init__(self, element=None, coordinates=None):
#         self.element = element
#         self.coordinates = list(coordinates)

#     def __str__(self):
#         return OrderedDict({'element': self.element, 'cordinates': [self.coordinates[0], self.coordinates[1], self.coordinates[2]]})
    
#     def __eq__(self, other): 
#         if not isinstance(other, Atom):
#             # don't attempt to compare against unrelated types
#             return NotImplemented

#         return self.element == other.element and self.coordinates == other.coordinates



# class SimulationStep:
    
#     def __init__(self):
#         self.name = ""
#         self.raw = ""
#         self.atoms = pd.Series()
#         self.metadata = dict()

#     def __len__(self):
#         return len(self.atoms)

#     def __str__(self):
#         return self.raw

#     def get_atoms(self):
#         return self.atoms
    
#     def get_metadata(self):
#         return self.metadata
    
#     def __getitem__(self, index):
#         return self.atoms[index]
    
# class SimulationTrajectory:
    
#     def __init__(self, steps=None, metadata=None):
#         self.steps = pd.Series()
#         self.metadata = pd.Series()

#     def __len__(self):
#         return len(self.steps)

#     def __str__(self):
#         return self.steps.__str__()
    
#     def get_metadata(self):
#         return self.metadata
    
#     def __getitem__(self, index):
#         return self.atoms[index]