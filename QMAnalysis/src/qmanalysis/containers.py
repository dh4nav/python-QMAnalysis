import pandas as pd

class AtomData:
    def __init__(self):
        self.dataframe =  pd.DataFrame(columns = ["file_name", "file_index", "timestep_time", "timestep_index", "atom_index", "element", "alias", "charge", "x", "y", "z" ])   

class TimestepData:
    def __init__(self):
        self.dataframe =  pd.DataFrame(columns = ["file_name", "file_path", "timestep_time", "timestep_index", "raw_data", "energy", "zero-point energy", "file_comment", "measurements" ])     

class MeasurementData:
    def __init__(self, timestep_names):
        # indexlist = []
        # for measurement_type in yamldata["measurements"]:
        #     for measurement_name in yamldata["measurements"][measurement_type]:
        #         indexlist.append(measurement_type + "-" + measurement_name)

        # print(indexlist)

        self.dataframe = pd.DataFrame({"file_name": timestep_names}, index=list(range(len(timestep_names)))) #indexlist, )
    
        print(self.dataframe)
