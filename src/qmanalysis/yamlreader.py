import strictyaml as sy
from collections import OrderedDict


class YAMLFile:
    def load_file(self, filepath: str):
        self.filepath = filepath
        with open(filepath, "r", encoding="utf-8") as f:
            self.raw = f.read()

        self.parsed = sy.load(self.raw, self._schema())
        self.data = self.parsed.data
        # self._validate_measurements()
        return self

    def load_string(self, yamlstring):
        self.raw = yamlstring
        self.parsed = sy.load(self.raw, self._schema())
        self.data = self.parsed.data
        # self._validate_measurements()
        return self

    def __str__(self):
        return self.raw

    def get_data(self):
        return self.data

    def _schema(self):
        from strictyaml import Map, Seq, Str, Int, MapPattern, Optional, Regex, Bool

        identifier = Str() | Int()
        valid_substitution = Regex(r"S\d+")

        return Map({
            "name": Str(),
            "comment": Str(),
            "version": Int(),

            "files": Seq(
                Map({
                    "path": Str(),
                    Optional("name"): Str(),
                    "type": Str(),
                    # Support for global_constants_csv and frame_constants_csv
                    Optional("options"): Str() | Int() | MapPattern(Str(), Str()),
                    Optional("glob"): Bool(),
                    Optional("timestep"): Str()
                })
            ),

            Optional("substitutions"): Seq(
                Map({
                    "name": Str(),
                    "entries": Seq(
                        Map({
                            "file": Str(),
                            "atom_index": Int(),
                            Optional("glob"): Bool()
                        })
                    )
                })
            ),

            Optional("sequences"): Seq(
                Map({
                    "name": Str(),
                    "entries": Seq(identifier)
                })
            ),

            Optional("measurements"): Map({
                Optional("distance"): Seq(
                    Map({
                        "name": Str(),
                        "a": identifier,
                        "b": identifier,
                        Optional("timestep"): Str
                    })
                ),
                Optional("angle"): Seq(
                    Map({
                        "name": Str(),
                        "a": identifier,
                        "b": identifier,
                        "c": identifier,
                        Optional("timestep"): Str
                    })
                ),
                Optional("dihedral"): Seq(
                    Map({
                        "name": Str(),
                        "a": identifier,
                        "b": identifier,
                        "c": identifier,
                        "d": identifier,
                        Optional("timestep"): Str
                    })
                )
            }),

            # Add calc logic
            Optional("calc"): Seq(
                Map({
                    "name": Str(),
                    "expr": Str()
                })
            ),

            "output": Seq(
                Map({
                    Optional("graph"): Seq(
                        Map({
                            "type": Str(),
                            "x": identifier,
                            "y": identifier,
                            Optional("series_by"): Str(),
                            Optional("parallel_by"): Str(),
                            "file": Str(),
                            Optional("file_format"): Str(),
                            Optional("dpi"): Int(),
                            Optional("figsize"): Seq(Int())
                        })
                    ),
                    Optional("file"): Seq(
                        Map({
                            "path": Str(),
                            "type": Str()
                        })
                    )
                })
            )

        })

    # def _validate_measurements(self):
    #     data = self.data
    #     defined_subs = {}
    #     defined_atoms = set()

    #     for item in data.get("substitutions", []):
    #         for key, entries in item.items():
    #             defined_subs[key] = {entry["atom"] for entry in entries}
    #             defined_atoms.update(defined_subs[key])

    #     def is_valid_identifier(x):
    #         return isinstance(x, int) and x in defined_atoms or isinstance(x, str) and x in defined_subs

    #     def assert_unique_ids(measure, keys):
    #         values = [measure[k] for k in keys]
    #         if len(set(values)) != len(values):
    #             raise ValueError(
    #                 f"[{measure['name']}] has duplicate identifiers: {values}"
    #             )

    #     if "measurements" in data:
    #         for kind, keys in {
    #             "distance": ["a", "b"],
    #             "angle": ["a", "b", "c"],
    #             "dihedral": ["a", "b", "c", "d"]
    #         }.items():
    #             for m in data["measurements"].get(kind, []):
    #                 for k in keys:
    #                     if not is_valid_identifier(m[k]):
    #                         raise ValueError(
    #                             f"[{m['name']}] {kind}: '{m[k]}' is not a valid atom or substitution reference."
    #                         )
    #                 assert_unique_ids(m, keys)
