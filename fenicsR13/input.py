"""
Module for input related Classes.

Contains the Input class.
"""

from functools import reduce
import operator
from json import dumps
import yaml
from cerberus import Validator


class Input:
    """
    Class to handle the input file in YAML_ format.

    .. _YAML: https://en.wikipedia.org/wiki/YAML

    An example input file could look like:

    .. code-block:: yaml

        # General
        # =======
        # - output_folder: Used as output folder
        output_folder: r13_1_coeffs_sources_rot_noinflow_p2p2p2p2p2_stab

        # Meshes
        # ======
        # - meshes: List of input meshes in h5 format to run simulations on
        meshes:
          - ../mesh/ring0.h5
          - ../mesh/ring1.h5
          - ../mesh/ring2.h5
          - ../mesh/ring3.h5
          - ../mesh/ring4.h5
          # - ../mesh/ring5.h5
          # - ../mesh/ring6.h5
          # - ../mesh/ring7.h5
          # - ../mesh/ring8.h5

        # Numerical Parameters
        # ====================
        # - stabilization: Must contain cip
        #   - cip: Collection of Continous Interior Penalty (CIP) parameters
        #     - enable: Enable CIP stabilization
        #     - delta_theta: Stabilization of grad(T)*grad(T_test) over edge
        #     - delta_u: Stabilization of grad(u)*grad(u_test) over edge
        #     - delta_p: Stabilization of grad(p)*grad(p_test) over edge
        stabilization:
          cip:
            enable: True
            delta_theta: 1.0
            delta_u: 1.0
            delta_p: 0.01

        # Formulation Parameters
        # ======================
        # - nsd: Number of spatial dimensions == 2
        # - kn: Knudsen numberkn
        # - heat_source: Heat source function
        # - mass_source: Mass source function
        # - body_force: Body force
        nsd: 2
        heat_source: 0
        mass_source: 1.0 * (1.0 - (5.0*pow(R,2))/(18.0*pow(0.1,2))) * cos(phi)
        body_force: [0,0]

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
        #     - u_t_w: Value for tangential velocity at wall
        #     - u_n_w: Value for normal velocity at wall
        #     - p_w: Value for pressure at wall
        #     - epsilon_w: Inflow-model parameter <=> Weight of pressure
        bcs:
          3000:
            chi_tilde: 1.0
            theta_w: 1.0
            u_t_w: -10
            u_n_w: 0
            p_w: 0
            epsilon_w: 0
          3100:
            chi_tilde: 1.0
            theta_w: 0.5
            u_t_w: 0
            u_n_w: 0
            p_w: 0
            epsilon_w: 0

    Further input examples can be found in the ``tests`` or ``examples``
    folders.

    Examples
    --------

    >>> corrupt_input = Input("etc/doctest_data/corrupt_input.yml")
    Traceback (most recent call last):
    ...
    Exception: Parsing error

    >>> working_input = Input("tests/heat/inputs/heat_01_coeffs_p1p1_stab.yml")
    Input:...

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
                        "u_t_w": {
                            "anyof": [{"type": "string"}, {"type": "float"}],
                            "required": True
                        },
                        "u_n_w": {
                            "anyof": [{"type": "string"}, {"type": "float"}],
                            "required": True
                        },
                        "p_w": {
                            "anyof": [{"type": "string"}, {"type": "float"}],
                            "required": True
                        },
                        "epsilon_w": {
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
            "body_force": {
                "type": "list",
                "required": True,
                "schema": {"anyof": [{"type": "string"}, {"type": "float"}]}
            },
            "output_folder": {
                "type": "string",
                "required": True,
            },
            "stabilization": {
                "type": "dict",
                "required": True,
                "schema": {
                    "cip": {
                        "type": "dict",
                        "required": True,
                        "schema": {
                            "enable": {
                                "type": "boolean",
                                "required": True
                            },
                            "delta_theta": {
                                "type": "float",
                                "required": True
                            },
                            "delta_u": {
                                "type": "float",
                                "required": True
                            },
                            "delta_p": {
                                "type": "float",
                                "required": True
                            },
                        }
                    },
                }
            },
         }

        if not val.validate(self.dict, input_schema):
            print("! Parsing Error: \n" + str(val.errors) + "\n")
            raise Exception("Parsing error")

        print("Input:\n" + dumps(self.dict, indent=1))

    def get_from_input(self, map_list):
        """
        Get value of the input dict by using a list of stacked keys.

        Source: https://stackoverflow.com/questions/14692690/..
                ..access-nested-dictionary-items-via-a-list-of-keys/37704379
        """
        try:
            return reduce(operator.getitem, map_list, self.dict)
        except Exception:
            raise Exception("Dict has no entry with the key:", map_list)

    def set_in_input(self, map_list, value):
        """Change value of the input dict by using a list of stacked keys."""
        self.get_from_input(map_list[:-1])[map_list[-1]] = value
