"Class to handle input file"

from json import dumps
import yaml
from cerberus import Validator

class Input:
    "Class to store input data as dict"

    def __init__(self, yaml_file):
        with open(yaml_file, "r") as stream:
            self.dict = yaml.safe_load(stream)

        val = Validator()
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
            "mode": {
                "type": "string",
                "required": True,
                "allowed": ["heat", "stress", "coupled"]
            },
            "use_coeffs": {
                "type": "boolean",
                "required": True,
            },
            "tau": {
                "type": "float",
                "required": True,
                "min": 0.000000001
            },
            "xi_tilde": {
                "type": "float",
                "required": True,
                "min": 0.000000001
            },
            "theta_w_inner": {
                "type": "float",
                "required": True,
                "min": 0.000000001
            },
            "theta_w_outer": {
                "type": "float",
                "required": True,
                "min": 0.000000001
            },
            "v_t_inner": {
                "type": "float",
                "required": True,
            },
            "v_t_outer": {
                "type": "float",
                "required": True,
            },
            "heat_source": {
                "type": "string",
                "required": True,
            },
            "mass_source": {
                "type": "string",
                "required": True,
            },
            "convergence_study": {
                "type": "dict",
                "required": True,
                "schema": {
                    "enable": {
                        "type": "boolean",
                        "required": True
                    },
                    "exact_solution": {
                        "type": "string",
                        "required": True
                    },
                    "plot": {
                        "type": "boolean",
                        "required": True
                    },
                }
            },
            "output_folder": {
                "type": "string",
                "required": True,
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

        if not val.validate(self.dict, input_schema):
            print(val.errors)
            raise Exception("Parsing error")

        print("Input:\n" + dumps(self.dict, indent=4))

# test = Input("input.yml")
