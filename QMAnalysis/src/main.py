import argparse
import strictyaml as sy
import qmanalysis.xyzreader as xr
import qmanalysis.yamlreader as yr
import qmanalysis.measure as mr
from qmanalysis.containers import AtomData, TimestepData, MeasurementData

def main():

  testyaml="""
  name: Procedure Name
  comment: XYZ
  version: 111
  files:
    - path: "../benzene.xyz"
      type: xyz
      name: benz
    - path: "../benzene.xyz"
      type: xyz
      name: benz2  
    - path: "../benzene.xyz"
      type: xyz
      name: benz3  
  substitutions:
    - name: S1
      entries:
      - file: benz
        atom_index: 5
      - file: benz2
        atom_index: 4
      - file: benz3
        atom_index: 1
  measurements:
    distance:
      - name: O-H-bond
        a: 6
        b: 7
      - name: Subs-bond
        a: S1
        b: 10
    angle:
      - name: H-O-H
        a: 3
        b: 6
        c: 7
      - name: H-O-N
        a: 8
        b: 9
        c: 10
    dihedral:
      - name: torsion1
        a: 2
        b: 3
        c: 6
        d: 7
  """


  # Example usage
  try:
      parser=argparse.ArgumentParser()
      parser.add_argument("inputfile", help="Main command input file", )
      args = parser.parse_args()
      yamlparser = yr.YAMLFile()
      yamlparser.load_string(testyaml)
      #validate_measurements(data)
  except (sy.YAMLError, ValueError) as e:
      print(f"Validation error: {e}")


  atom_data = AtomData()
  timestep_data = TimestepData()

  print(yamlparser.get_data())

  yamldata = yamlparser.get_data()

  for file in yamldata["files"]:
    if file["type"].lower() == "xyz":
      xr.XYZFile(atom_data, timestep_data, file_path = file["path"], file_name = file["name"])
  print(atom_data.dataframe)
  print(timestep_data.dataframe)

  if "substitutions" in yamldata:
    for sub in yamldata["substitutions"]:
       for subfile in sub['entries']:
          #mask = (atom_data.dataframe["file_name"] == subfile['file']) & (atom_data.dataframe["atom_index"] == subfile['atom_index'])
          #atom_data.dataframe[mask, ["alias"]] == 777
          atom_data.dataframe.loc[(atom_data.dataframe["file_name"] == subfile['file']) & (atom_data.dataframe["atom_index"] == subfile['atom_index']), "alias"] = sub["name"]
  
  print(atom_data.dataframe)

  timestep_names = timestep_data.dataframe.loc[:,"file_name"].to_numpy()

  print(timestep_names)

  measure = mr.Measure()

  md = MeasurementData(timestep_names)

  print("X")

  temp_array = []

  if "measurements" in yamldata:
    if "distance" in yamldata["measurements"]:
       for distance in yamldata["measurements"]["distance"]:
          print(distance['a'])
          temp_array = []
          for timestep_name in timestep_names:
            print(timestep_name)
            print(atom_data.dataframe[(atom_data.dataframe['file_name'] == timestep_name) & (atom_data.dataframe['alias'] == str(distance['a'])  )])
            temp_array.append(measure.distance(
                atom_data=atom_data, 
                atom_index1=atom_data.dataframe[(atom_data.dataframe['file_name'] == timestep_name) & (atom_data.dataframe['alias'] == str(distance['a'])  )].index[0],
                atom_index2=atom_data.dataframe[(atom_data.dataframe['file_name'] == timestep_name) & (atom_data.dataframe['alias'] == str(distance['b'])  )].index[0]))
          md.dataframe["distance-"+distance["name"]] = temp_array
          print(md.dataframe)
          
    if "angle" in yamldata["measurements"]:
       for angle in yamldata["measurements"]["angle"]:
          print(angle)
          temp_array = []
          for timestep_name in timestep_names:
            temp_array.append(measure.angle(
                atom_data=atom_data, 
                atom_index1=atom_data.dataframe[(atom_data.dataframe['file_name'] == timestep_name) & (atom_data.dataframe['alias'] == str(angle['a']))].index[0],
                atom_index2=atom_data.dataframe[(atom_data.dataframe['file_name'] == timestep_name) & (atom_data.dataframe['alias'] == str(angle['b']))].index[0],
                atom_index3=atom_data.dataframe[(atom_data.dataframe['file_name'] == timestep_name) & (atom_data.dataframe['alias'] == str(angle['c']))].index[0]))
          md.dataframe["angle-"+angle["name"]] = temp_array
    if "dihedral" in yamldata["measurements"]:
       for dihedral in yamldata["measurements"]["dihedral"]:
          print(dihedral)
          temp_array = []
          for timestep_name in timestep_names:
            temp_array.append(measure.dihedral(
                atom_data=atom_data, 
                atom_index1=atom_data.dataframe[(atom_data.dataframe['file_name'] == timestep_name) & (atom_data.dataframe['alias'] == str(dihedral['a']))].index[0],
                atom_index2=atom_data.dataframe[(atom_data.dataframe['file_name'] == timestep_name) & (atom_data.dataframe['alias'] == str(dihedral['b']))].index[0],
                atom_index3=atom_data.dataframe[(atom_data.dataframe['file_name'] == timestep_name) & (atom_data.dataframe['alias'] == str(dihedral['c']))].index[0],
                atom_index4=atom_data.dataframe[(atom_data.dataframe['file_name'] == timestep_name) & (atom_data.dataframe['alias'] == str(dihedral['d']))].index[0]))
          md.dataframe["dihedral-"+dihedral["name"]] = temp_array
  print(md.dataframe)
  print(atom_data.dataframe)
  print(timestep_data.dataframe)
  plot = md.dataframe.plot.scatter(x="distance-O-H-bond", y="distance-Subs-bond")

