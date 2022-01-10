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
        # - elements: Must contain the fields: theta, s, p, u, sigma
        #   - fields: List of FEM parameters (shape, degree)
        #     - shape: Element shape, e.g. Lagrange
        #     - degree: Element degree, e.g. 2
        # - stabilization: Must contain cip and gls
        #   - cip: Collection of Continous Interior Penalty (CIP) parameters
        #     - enable: Enable CIP stabilization
        #     - delta_theta: Stabilization of grad(T)*grad(T_test) over edge
        #     - delta_u: Stabilization of grad(u)*grad(u_test) over edge
        #     - delta_p: Stabilization of grad(p)*grad(p_test) over edge
        #   - gls: Collection of Garlerkin Least Squares (GLS) parameters
        #     - enable: Enable GLS stabilization
        #     - tau_energy: Stabilization with energy eq. residual
        #     - tau_heatflux: Stabilization with heatflux eq. residual
        #     - tau_mass: Stabilization with mass eq. residual
        #     - tau_momentum: Stabilization with momentum eq. residual
        #     - tau_stress: Stabilization with stress eq. residual
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
            delta_theta: 1.0
            delta_u: 1.0
            delta_p: 0.01
          gls:
            enable: False
            tau_energy: 0.001
            tau_heatflux: 0.001
            tau_mass: 0.01
            tau_momentum: 0.01
            tau_stress: 0.01

        # Formulation Parameters
        # ======================
        # - nsd: Number of spatial dimensions == 2
        # - mode: Formulation mode, one of heat, stress, r13
        # - kn: Knudsen numberkn
        # - heat_source: Heat source function for mode==heat||r13
        # - mass_source: Mass source function for mode==stress||r13
        # - body_force: Body force for mode==stress||r13
        nsd: 2
        mode: r13
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
        # - write_vecs: Write all solution fields as vectors
        # - massflow: List of BC IDs to compute massflow J=int_bc dot(u,n) ds
        # - line_integrals: List of line integral dicts:
        #   - name: Name for output
        #   - expr: Expression to evaluate
        #   - start: Start point
        #   - end: End point
        #   - res: Sampling resolution of line
        postprocessing:
        write_pdfs: True
        write_vecs: True
        massflow: []
        line_integrals:
          - name: "avg(abs(uy(x,y=0.5)))"
            expr: abs(uy)/8
            start: [0, 0.5]
            end: [8, 0.5]
            res: 10000
          - name: "avg(abs(uy(x,y=4.5)))"
            expr: abs(uy)/8
            start: [0, 4.5]
            end: [8, 4.5]
            res: 10000

        # Parameter Study
        # ==============
        # - enable: Repeat simulation with different p. values (study)
        # - parameter_key: Key as list, e.g. ["elemenets", "p", "degree"]
        # - parameter_values: List of value for parameter, e.g. [0.01,0.1,1,10]
        parameter_study:
        enable: True
        parameter_key: ["kn"]
        parameter_values: [1,2,3]

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
                "allowed": [2, 3]
            },
            "mode": {
                "type": "string",
                "required": True,
                "allowed": ["heat", "stress", "r13"]
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
            "polar_coord_syst": {
                "type": "boolean",
                "required": True,

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
                        "u_x_w": {
                            "anyof": [{"type": "string"}, {"type": "float"}],
                            "required": True
                        },
                        "u_y_w": {
                            "anyof": [{"type": "string"}, {"type": "float"}],
                            "required": True
                        },
                        "u_z_w": {
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
            "solver": {
                "type": "dict",
                "required": True,
                "schema": {
                    "solver_name": {
                        "type": "string",
                        "required": True,
                    },
                    "preconditioner": {
                        "type": "string",
                        "required": True,
                    }

                }
            },
            "postprocessing": {
                "type": "dict",
                "required": True,
                "schema": {
                    "write_pdfs": {
                        "type": "boolean",
                        "required": True
                    },
                    "write_vecs": {
                        "type": "boolean",
                        "required": True
                    },
                    "massflow": {
                        "type": "list",
                        "required": True,
                        "schema": {"type": "integer"}
                    },
                    "line_integrals": {
                        "type": "list",
                        "required": True,
                        "schema": {
                            "type": "dict",
                            "schema": {
                                "name": {
                                    "type": "string",
                                    "required": True,
                                },
                                "expr": {
                                    "anyof": [
                                        {"type": "string"}, {"type": "float"}
                                    ],
                                    "required": True,
                                },
                                "start": {
                                    "type": "list",
                                    "required": True,
                                    "schema": {"anyof": [
                                        {"type": "string"}, {"type": "float"}
                                    ]}
                                },
                                "end": {
                                    "type": "list",
                                    "required": True,
                                    "schema": {"anyof": [
                                        {"type": "string"}, {"type": "float"}
                                    ]}
                                },
                                "res": {
                                    "type": "integer",
                                    "required": True,
                                    "min": 1
                                },
                            },
                        },
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
                    "gls": {
                        "type": "dict",
                        "required": True,
                        "schema": {
                            "enable": {
                                "type": "boolean",
                                "required": True
                            },
                            "tau_energy": {
                                "type": "float",
                                "required": True
                            },
                            "tau_heatflux": {
                                "type": "float",
                                "required": True
                            },
                            "tau_mass": {
                                "type": "float",
                                "required": True
                            },
                            "tau_momentum": {
                                "type": "float",
                                "required": True
                            },
                            "tau_stress": {
                                "type": "float",
                                "required": True
                            },
                        }
                    },
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
            print("! Parsing Error: \n" + str(val.errors) + "\n")
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
        except Exception:
            raise Exception("Dict has no entry with the key:", map_list)

    def set_in_input(self, map_list, value):
        """Change value of the input dict by using a list of stacked keys."""
        self.get_from_input(map_list[:-1])[map_list[-1]] = value
