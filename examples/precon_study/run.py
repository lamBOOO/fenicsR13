import os
import csv
import re
import json
import subprocess
import itertools


def extract_log_info(logfile_path):

    solved_times = []
    iteration_values = []
    last_iteration = None

    with open(logfile_path, 'r') as logfile:
        for line in logfile:
            # Suche in jeder Zeile nach einem Iterationswert
            iteration_match = re.search(
                r'(\d+)\s+KSP (?:preconditioned|unpreconditioned) resid norm',
                line,
                re.IGNORECASE
            )
            if iteration_match:
                last_iteration = int(iteration_match.group(1))

            # Suche nach einem "Finished solve:"-Eintrag
            solve_match = re.search(
                r'Finished solve:\s*(\d+\.?\d*)',
                line,
                re.IGNORECASE
            )
            if solve_match:
                solved_time = float(solve_match.group(1))
                solved_times.append(solved_time)
                iteration_values.append(last_iteration)
                last_iteration = None

    return solved_times, iteration_values


def collect_log_data(base_dir):

    log_data = []

    for root, dirs, files in os.walk(base_dir):
        for file in files:
            if file.endswith('.log'):
                logfile_path = os.path.join(root, file)
                folder_name = os.path.basename(root)

                solved_times, iterations = extract_log_info(logfile_path)

                if solved_times and iterations:
                    for solved_time, iteration in zip(solved_times, iterations):
                        log_data.append((folder_name, solved_time, iteration))
                    print(f"Valid log found: {logfile_path}")
                else:
                    print(f"Invalid log skipped: {logfile_path}")

    return log_data


def write_to_csv(output_path, log_data):

    # Gruppiere die Logdaten nach Folder Name
    grouped_data = {}
    for folder_name, solved_time, iteration in log_data:
        if folder_name not in grouped_data:
            grouped_data[folder_name] = {"solved_times": [], "iterations": []}
        grouped_data[folder_name]["solved_times"].append(solved_time)
        grouped_data[folder_name]["iterations"].append(iteration)

    # Schreibe die gruppierten Daten in die CSV-Datei
    with open(output_path, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        # Schreibe den Header
        csvwriter.writerow(['Folder Name', 'Solved Times', 'Iterations'])

        for folder, data in grouped_data.items():
            solved_times_str = ", ".join(map(str, data["solved_times"]))
            iterations_str = ", ".join(map(str, data["iterations"]))
            csvwriter.writerow([folder, solved_times_str, iterations_str])


def update_and_run_simulations():

    # Pfade zur Original-YAML-Datei und zum Verzeichnis für Ausgaben
    original_yaml = (
        "../../tests/square_manufactured_solution/inputs/"
        "square_manufactured_solution_2d_r13.yml"
    )
    # Verzeichnis für Logs und Plots
    output_dir = "../../examples/precon_study/Simulation_Output_v1"

    # Optionen für das Kartesische Produkt
    ksp_types = ["gmres", "fgmres"]  # , "tfqmr"]
    pc_types = ["fieldsplit"]
    fieldsplit_types = ["schur"]  # , "additive",
    # "multiplicative", "symmetric_multiplicative"]
    schur_fact_types = ["lower"]  # ,"full", "diag", "upper"]
    fieldsplit_0_types = ["hypre"]  # ,"jacobi", "bjacobi", "sor","lu", "qr",
    # "shell", "eisenstat", "ilu", "icc", "asm",
    # "gasm","ksp","composite","cholesky", "gamg", "bddc"]
    fieldsplit_1_types = ["jacobi"]  # ,"none"

    # Konstante Optionen, die für alle Simulationen beibehalten werden
    constant_options = {
        "petsc_options.fieldsplit_0_ksp_type": "preonly",
        "petsc_options.fieldsplit_1_ksp_type": "preonly",
        "petsc_options.fieldsplit_0_pc_gamg_type": "classical"
    }

    # Erzeuge alle Kombinationen der Optionen
    combinations = []
    for ksp in ksp_types:
        for pc in pc_types:
            if pc == "fieldsplit":
                for f in fieldsplit_types:
                    for f0 in fieldsplit_0_types:
                        for f1 in fieldsplit_1_types:
                            if f == "schur":
                                for s in schur_fact_types:
                                    combinations.append((ksp, pc, f, s, f0, f1))
                            else:
                                combinations.append((ksp, pc, f, None, f0, f1))
            else:
                # Für andere pc-Typen wird das volle Produkt verwendet.
                for f, s in itertools.product(
                        fieldsplit_types, schur_fact_types):
                    combinations.append((ksp, pc, f, s, f0, f1))

    for idx, (
        ksp_type,
        pc_type,
        fieldsplit_type,
        schur_fact_type, fieldsplit_0_type, fieldsplit_1_type
    ) in enumerate(combinations):
        # Erstellen des spezifischen Output-Ordners im precon_study-Verzeichnis
        if schur_fact_type is None:
            output_folder_name = (
                f"{ksp_type}_{fieldsplit_type}_"
                f"{fieldsplit_0_type}_{fieldsplit_1_type}"
            )
        else:
            output_folder_name = (
                f"{ksp_type}_{fieldsplit_type}_{schur_fact_type}_"
                f"{fieldsplit_0_type}_{fieldsplit_1_type}")
        output_folder = os.path.join(output_dir, output_folder_name)
        os.makedirs(output_folder, exist_ok=True)

        if pc_type == "fieldsplit":
            if fieldsplit_type == "schur":
                updates = {
                    "output_folder": output_folder,
                    "petsc_options.ksp_type": ksp_type,
                    "petsc_options.pc_type": pc_type,
                    "petsc_options.pc_fieldsplit_type": fieldsplit_type,
                    "petsc_options.pc_fieldsplit_schur_fact_type": (
                        schur_fact_type
                    ),
                    "petsc_options.fieldsplit_0_pc_type": fieldsplit_0_type,
                    "petsc_options.fieldsplit_1_pc_type": fieldsplit_1_type
                }
                updates.update(constant_options)
            else:
                updates = {
                    "output_folder": output_folder,
                    "petsc_options.ksp_type": ksp_type,
                    "petsc_options.pc_type": pc_type,
                    "petsc_options.pc_fieldsplit_type": fieldsplit_type,
                    "petsc_options.fieldsplit_0_pc_type": fieldsplit_0_type,
                    "petsc_options.fieldsplit_1_pc_type": fieldsplit_1_type
                }
                updates.update(constant_options)

        print(
            f"Updating YAML with the following options: "
            f"{json.dumps(updates, indent=2)}"
        )

        output_file_name = f"config_update_{idx + 1}.yml"
        yaml_output_file_path = os.path.join(
            '../../tests/square_manufactured_solution/inputs/',
            output_file_name)

        subprocess.run(["python3", "../../examples/precon_study/update_yml.py",
                        original_yaml, yaml_output_file_path,
                        json.dumps(updates)], check=True)

        log_file_name = f"update_{idx + 1}.log"
        log_file_path = os.path.join(output_folder, log_file_name)

        subprocess.run(
            ["fenicsR13",
             f"../../tests/square_manufactured_solution/inputs/"
             f"{output_file_name}"],
            stdout=open(log_file_path, "w"),
            stderr=subprocess.STDOUT
        )

        # Lösche alle Dateien außer dem Log im Ordner
        for file in os.listdir(output_folder):
            file_path = os.path.join(output_folder, file)
            if file != f"update_{idx + 1}.log" and os.path.isfile(file_path):
                os.remove(file_path)

        plot_result = subprocess.run([
            "python3",
            "../../examples/precon_study/convergence_plots.py",
            log_file_path
        ])
        if plot_result.returncode != 0:
            print(f"Error creating plot for log file: {log_file_path}")


def main():
    base_dir = "../../examples/precon_study/"
    output_csv = os.path.join(base_dir, "log_data_gmres_fgmres.csv")

    update_and_run_simulations()
    log_data = collect_log_data(base_dir)
    write_to_csv(output_csv, log_data)
    print(f"Daten wurden erfolgreich in {output_csv} gespeichert.")


if __name__ == '__main__':
    main()
