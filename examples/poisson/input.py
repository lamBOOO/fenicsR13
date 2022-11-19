"""
Module for input related Classes.

Contains the Input class.
"""

from json import dumps
import yaml
from cerberus import Validator


class Input:
    """
    Class to handle the input file in YAML_ format.

    .. _YAML: https://en.wikipedia.org/wiki/YAML

    Content:

        # General
        # =======
        # - output_folder: Used as output folder

        # Meshes
        # ======
        # - meshes: List of input meshes in h5 format to run simulations on
        
        # Formulation Parameters
        # ======================
        # - nsd: Number of spatial dimensions == 2
        # - heat_source: Heat source function
        # - mass_source: Mass source function

        # Region Parameters
        # =================
        # - regs: Dictionary of all mesh regions
        #   - reg_id: Must contain the following parameters:
        #     - kn: Knudsen number

        # Boundary Conditions
        # ===================
        # - bcs: Dictionary of all boundary IDs from mesh
        #   - bc_id: must contain the following parameters
        #     - chi_tilde: Refaction coefficient in Maxwell accomodation model
        #     - theta_w: Value for temperature at wall
        #     - p_w: Value for pressure at wall

    """

    def __init__(self, yaml_file):
        """Construct the Input class."""
        try:
            with open(yaml_file, "r") as stream:
                self.dict = yaml.safe_load(stream)
        except FileNotFoundError:
            print(f"[Error] The input file {yaml_file} is not found.\n")
            raise

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
            "regs": {
                "type": "dict",
                "required": True,
                "keysrules": {"type": "integer"},
                "valuesrules": {
                    "type": "dict",
                    "schema": {
                        "kn": {
                            "anyof": [
                                {"type": "string"},
                                {"type": "float", "min": 1E-10}
                            ],
                            "required": True,
                        },
                    }
                }
            },
            "bcs": {
                "type": "dict",
                "required": True,
                "keysrules": {"type": "integer"},
                "valuesrules": {
                    "type": "dict",
                    "schema": {
                        "chi_tilde": {
                            "anyof": [
                                {"type": "string"},
                                {"type": "float", "min": 1E-10}
                            ],
                            "required": True,
                        },
                        "theta_w": {
                            "anyof": [{"type": "string"}, {"type": "float"}],
                            "required": True,
                        },
                        "p_w": {
                            "anyof": [{"type": "string"}, {"type": "float"}],
                            "required": True
                        },
                    }
                }
            },
            "heat_source": {
                "anyof": [{"type": "string"}, {"type": "float"}],
                "required": True,
            },
            "mass_source": {
                "anyof": [{"type": "string"}, {"type": "float"}],
                "required": True,
            },
            "output_folder": {
                "type": "string",
                "required": True,
            },
         }

        if not val.validate(self.dict, input_schema):
            print("! Parsing Error: \n" + str(val.errors) + "\n")
            raise Exception("Parsing error")

        print("Input:\n" + dumps(self.dict, indent=1))
