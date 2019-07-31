#!/usr/bin/env python3

# pylint: disable=invalid-name
# pylint: disable=unsubscriptable-object

"""
Program to solve the decoupled (removed coupling term) heat system of the
linearized R13 equations
"""

import sys
import gc
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

            if params["convergence_study"]["write_systemmatrix"]:
                solver.write_systemmatrix()

            solver.load_exact_solution()
            solver.calc_errors()

            errors = solver.errors

            data.append({
                "h": current_mesh.mesh.hmax(),
                **errors
            })

            if p == len(mesh_names)-1: # after last mesh
                postp = Postprocessor(data, params["case_name"])
                postp.write_errors()
                if plot:
                    postp.plot_errors()

        solver = None
        gc.collect()

if __name__ == '__main__':
    main()
