import yaml
from collections import OrderedDict


class YAMLFile:
    def _validate_data(self, data):
        """Basic validation of required fields"""
        required_fields = ['name', 'comment', 'version', 'files']
        for field in required_fields:
            if field not in data:
                raise ValueError(f"Missing required field: {field}")

        # Validate files structure
        if not isinstance(data['files'], list):
            raise ValueError("'files' must be a list")

        for file_entry in data['files']:
            if 'path' not in file_entry:
                raise ValueError("Each file entry must have a 'path' field")
            if 'type' not in file_entry:
                raise ValueError("Each file entry must have a 'type' field")

    def load_file(self, filepath: str):
        self.filepath = filepath
        with open(filepath, "r", encoding="utf-8") as f:
            self.raw = f.read()

        self.data = yaml.safe_load(self.raw)
        self._validate_data(self.data)
        return self

    def load_string(self, yamlstring):
        self.raw = yamlstring
        self.data = yaml.safe_load(self.raw)
        self._validate_data(self.data)
        return self

    def __str__(self):
        return self.raw

    def get_data(self):
        return self.data
