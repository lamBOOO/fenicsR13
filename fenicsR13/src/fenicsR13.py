#!/usr/bin/env python3

# pylint: disable=invalid-name
# pylint: disable=unsubscriptable-object

"""
Program to solve the decoupled (removed coupling term) heat system of the
linearized R13 equations
"""

import sys
import dolfin as df

import meshes
from input import Input
from solver import Solver
from postprocessor import Postprocessor

def main():
    "Main Program"

    # Dolfin settings
    df.set_log_level(100) # 1: all logs
    df.parameters["ghost_mode"] = "shared_vertex"

    inputfile = sys.argv[1] if len(sys.argv) == 2 else "input.yml"

    params = Input(inputfile).dict
    mesh_names = params["meshes"]

    convergence_study = params["convergence_study"]["enable"]
    plot = params["convergence_study"]["plot"]

    data = []

    for p, mesh_name in enumerate(mesh_names):

        print("Mesh: " + mesh_name)

        mesh_name = mesh_names[p]

        current_mesh = meshes.H5Mesh(mesh_name)
        solver = Solver(params, current_mesh, p)

        solver.setup_function_spaces()
        solver.assemble()
        solver.solve()
        solver.write_solutions()
        solver.write_parameters()

        if convergence_study:

            solver.load_exact_solution()
            solver.calc_errors()

            errors = solver.errors

            if params["mode"] == "heat":
                # FIXME: Resolve this if statement proberly
                data.append({
                    "h": current_mesh.mesh.hmax(),
                    "theta": {
                        "L_2": errors["f"]["l2"]["theta"],
                        "l_inf": errors["v"]["linf"]["theta"],
                    },
                    "sx": {
                        "L_2": errors["f"]["l2"]["s"][0],
                        "l_inf": errors["v"]["linf"]["s"][0],
                    },
                    "sy": {
                        "L_2": errors["f"]["l2"]["s"][1],
                        "l_inf": errors["v"]["linf"]["s"][1],
                    },
                })
            elif params["mode"] == "stress":
                data.append({
                    "h": current_mesh.mesh.hmax(),
                    "p": {
                        "L_2": errors["f"]["l2"]["p"],
                        "l_inf": errors["v"]["linf"]["p"],
                    },
                    "ux": {
                        "L_2": errors["f"]["l2"]["u"][0],
                        "l_inf": errors["v"]["linf"]["u"][0],
                    },
                    "uy": {
                        "L_2": errors["f"]["l2"]["u"][1],
                        "l_inf": errors["v"]["linf"]["u"][1],
                    },
                    "sigmaxx": {
                        "L_2": errors["f"]["l2"]["sigma"][0],
                        "l_inf": errors["v"]["linf"]["sigma"][0],
                    },
                    "sigmaxy": {
                        "L_2": errors["f"]["l2"]["sigma"][1],
                        "l_inf": errors["v"]["linf"]["sigma"][1],
                    },
                    "sigmayy": {
                        "L_2": errors["f"]["l2"]["sigma"][2],
                        "l_inf": errors["v"]["linf"]["sigma"][2],
                    }
                })

            if p == len(mesh_names)-1: # after last mesh
                postp = Postprocessor(data)
                postp.write_errors()
                if plot:
                    postp.plot_errors()

if __name__ == '__main__':
    main()
