"Class to handle input file"

import yaml
from cerberus import Validator

from schema import Schema, And, Use, Or



class Input:
    "Class to store input data as dict"

    def __init__(self, yaml_file):
        with open(yaml_file, "r") as stream:
            self.dict = yaml.safe_load(stream)

        v = Validator()
        input_schema = {
            "meshes": {
                "type": "list",
                "required": True,
                "schema": {"type": "string"}
            },
            "nsd": {
                "type": "integer",
                "required": True,
                "allowed": [2]
            },
            "stabilization": {
                "type": "dict",
                "required": True,
                "keysrules": {"type": "string", "regex": "cip"},
                "valueschema": {
                    "type": "dict",
                    "schema": {
                        "enable": {
                            "type": "boolean",
                            "required": True
                        },
                        "delta_1": {
                            "type": "float",
                            "required": True
                        },
                        "delta_2": {
                            "type": "float",
                            "required": True
                        },
                        "delta_3": {
                            "type": "float",
                            "required": True
                        },
                    }
                }
            }
        }

        if not v.validate(self.dict, input_schema):
            print(v.errors)
            raise Exception("Parsing error")

# test = Input("input.yml")