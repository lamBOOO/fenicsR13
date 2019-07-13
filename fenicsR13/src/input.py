"Class to handle input file"

import yaml
import schema

class Input:
    "Class to store input data as dict"

    def __init__(self, yaml_file):
        with open(yaml_file, 'r') as stream:
            self.dict = yaml.safe_load(stream)

        conf_schema = schema.Schema({
            'meshes': [str],
            'nsd': schema.And(int, lambda n: n == 2, error='nsd must be 2')
        })

        # Validate
        conf_schema.validate(self.dict)
