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
        # - case_name: Used as output folder
        case_name: r13_1_coeffs_sources_rot_noinflow_p2p2p2p2p2_stab

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
        # - elements: Must contain the fields: theta, s, p, u, sigma
        #   - fields: List of FEM parameters (shape, degree)
        #     - shape: Element shape, e.g. Lagrange
        #     - degree: Element degree, e.g. 2
        # - stabilization: Must contain cip
        #   - cip: Collection of Continous Interior Penalty (CIP) parameters
        #     - enable: Enable CIP stabilization
        #     - delta_1: Stabilization of grad(T)*grad(T_test) over edge
        #     - delta_2: Stabilization of grad(u)*grad(u_test) over edge
        #     - delta_3: Stabilization of grad(p)*grad(p_test) over edge
        elements:
          theta:
            shape: Lagrange
            degree: 2
          s:
            shape: Lagrange
            degree: 2
          p:
            shape: Lagrange
            degree: 2
          u:
            shape: Lagrange
            degree: 2
          sigma:
            shape: Lagrange
            degree: 2
        stabilization:
          cip:
            enable: True
            delta_1: 1.0
            delta_2: 1.0
            delta_3: 0.01

        # Formulation Parameters
        # ======================
        # - nsd: Number of spatial dimensions == 2
        # - mode: Formulation mode, one of heat, stress, r13
        # - use_coeffs: Use real R13 coefficients, False only valid in heat
        # - kn: Knudsen numberkn
        # - xi_tilde: Refaction coefficient in Maxwell accomodation model
        # - heat_source: Heat source function for mode==heat||r13
        # - mass_source: Mass source function for mode==stress||r13
        nsd: 2
        mode: r13
        use_coeffs: True
        kn: 1.0
        xi_tilde: 1.0
        heat_source: 0
        mass_source: 1.0 * (1.0 - (5.0*pow(R,2))/(18.0*pow(kn,2))) * cos(phi)

        # Boundary Conditions
        # ===================
        # - bcs: Dictionary of all boundary IDs from mesh
        #   - bc_id: must contain theta_w, u_t_w, u_n_w, p_w, epsilon_w
        #     - theta_w: Value for temperature at wall
        #     - u_t_w: Value for tangential velocity at wall
        #     - u_n_w: Value for normal velocity at wall
        #     - p_w: Value for pressure at wall
        #     - epsilon_w: Inflow-model parameter <=> Weight of pressure
        bcs:
          3000:
            theta_w: 1.0
            u_t_w: -10
            u_n_w: 0
            p_w: 0
            epsilon_w: 0
          3100:
            theta_w: 0.5
            u_t_w: 0
            u_n_w: 0
            p_w: 0
            epsilon_w: 0

        # Convergence Study
        # =================
        # - enable: Enable convergence study on given meshes
        # - exact_solution: Path to exact solution in cpp-format to compare
        # - plot: Show errors in matplotlib window. PDF output is default
        # - write_systemmatrix: Writes out systemmatrix (LHS)
        # - rescale_pressure: Rescale numerical pressure result for zero mean
        # - relative_errors: Use relative errors. If esol=0, use absolute.
        convergence_study:
          enable: True
          exact_solution: esols/1_coeffs_sources_rot_noinflow.cpp
          plot: False
          write_systemmatrix: False
          rescale_pressure: True
          relative_error: True

        # Postprocessing
        # ==============
        # - write_pdfs: Write all solution fields as PDF plot
        # - massflow: List of BC IDs to compute massflow J=int_bc dot(u,n) ds
        postprocessing:
        write_pdfs: True

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
                "allowed": ["heat", "stress", "r13"]
            },
            "use_coeffs": {
                "type": "boolean",
                "required": True,
            },
            "kn": {
                "type": "float",
                "required": True,
                "min": 0.000000001
            },
            "xi_tilde": {
                "type": "float",
                "required": True,
                "min": 0.000000001
            },
            "bcs": {
                "type": "dict",
                "required": True,
                "keysrules": {"type": "integer"},
                "valuesrules": {
                    "type": "dict",
                    "schema": {
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
            "postprocessing": {
                "type": "dict",
                "required": True,
                "schema": {
                    "write_pdfs": {
                        "type": "boolean",
                        "required": True
                    },
                    "massflow": {
                        "type": "list",
                        "required": True,
                        "schema": {"type": "integer"}
                    },
                }
            },
            "parameter_study": {
                "type": "dict",
                "required": True,
                "schema": {
                    "enable": {
                        "type": "boolean",
                        "required": True
                    },
                    "parameter_key": {
                        "type": "list",
                        "required": True,
                    },
                    "parameter_values": {
                        "type": "list",
                        "required": True,
                    },
                }
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
                    "write_systemmatrix": {
                        "type": "boolean",
                        "required": True
                    },
                    "rescale_pressure": {
                        "type": "boolean",
                        "required": True
                    },
                    "relative_error": {
                        "type": "boolean",
                        "required": True
                    },
                }
            },
            "case_name": {
                "type": "string",
                "required": True,
            },
            "stabilization": {
                "type": "dict",
                "required": True,
                "keysrules": {"type": "string", "regex": "cip"},
                "valuesrules": {
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
            },
            "elements": {
                "type": "dict",
                "required": True,
                "schema": {
                    "theta": {
                        "type": "dict",
                        "required": True,
                        "schema": {
                            "shape": {
                                "type": "string",
                                "required": True
                            },
                            "degree": {
                                "type": "integer",
                                "required": True,
                                "min": 0
                            }
                        }
                    },
                    "s": {
                        "type": "dict",
                        "required": True,
                        "schema": {
                            "shape": {
                                "type": "string",
                                "required": True
                            },
                            "degree": {
                                "type": "integer",
                                "required": True,
                                "min": 0
                            }
                        }
                    },
                    "p": {
                        "type": "dict",
                        "required": True,
                        "schema": {
                            "shape": {
                                "type": "string",
                                "required": True
                            },
                            "degree": {
                                "type": "integer",
                                "required": True,
                                "min": 0
                            }
                        }
                    },
                    "u": {
                        "type": "dict",
                        "required": True,
                        "schema": {
                            "shape": {
                                "type": "string",
                                "required": True
                            },
                            "degree": {
                                "type": "integer",
                                "required": True,
                                "min": 0
                            }
                        }
                    },
                    "sigma": {
                        "type": "dict",
                        "required": True,
                        "schema": {
                            "shape": {
                                "type": "string",
                                "required": True
                            },
                            "degree": {
                                "type": "integer",
                                "required": True,
                                "min": 0
                            }
                        }
                    },
                }
            }
        }

        if not val.validate(self.dict, input_schema):
            print(val.errors)
            raise Exception("Parsing error")

        print("Input:\n" + dumps(self.dict, indent=None))

    def get_from_input(self, map_list):
        """
        Get value of the input dict by using a list of stacked keys.

        Source: https://stackoverflow.com/questions/14692690/..
                ..access-nested-dictionary-items-via-a-list-of-keys/37704379
        """
        try:
            return reduce(operator.getitem, map_list, self.dict)
        except:
            raise Exception("Dict has no entry with the key:", map_list)

    def set_in_input(self, map_list, value):
        """Change value of the input dict by using a list of stacked keys."""
        self.get_from_input(map_list[:-1])[map_list[-1]] = value
