# -*- coding: utf-8 -*-
import re
import matplotlib.pyplot as plt
import os
import sys

# 1. Überprüfen, ob ein Argument übergeben wurde
if len(sys.argv) < 2:
    print("Bitte geben Sie den Pfad zur Logdatei an.")
    sys.exit(1)

logfile = os.path.abspath(sys.argv[1])

# Überprüfen, ob die Logdatei existiert
if not os.path.isfile(logfile):
    print(f"Logdatei nicht gefunden: {logfile}")
    exit(1)

# Verzeichnis für gespeicherte Plots
output_dir = os.path.dirname(logfile)
os.makedirs(output_dir, exist_ok=True)

# 2. Daten für jede Mesh-Größe speichern
data = {}

with open(logfile, 'r') as file:
    current_mesh_size = None
    for line in file:
        # 3. Muster anpassen: Mesh-Größe erfassen
        mesh_match = re.search(r'Mesh:\s*([^\s]+)', line)
        if mesh_match:
            current_mesh_size = mesh_match.group(1)
            if current_mesh_size not in data:
                data[current_mesh_size] = ([], [])

        if current_mesh_size:  # Nur wenn eine Mesh-Größe gesetzt ist
            residual_match = re.search(
                r"^\s*(\d+)\s+KSP.*\|\|r\(i\)\|\|/"
                r".*\|\|b\|\|\s+([\d\.eE+-]+)",
                line
            )
            if residual_match:
                # Iterationsnummer und relatives Residuum speichern
                iterations, relative_residues = data[current_mesh_size]
                iterations.append(int(residual_match.group(1)))
                relative_residues.append(float(residual_match.group(2)))

# 4. Plot erstellen
plt.figure(figsize=(10, 8))

for mesh_size, (iterations, relative_residues) in data.items():
    # Extrahiere nur den Dateinamen aus dem Pfad
    mesh_name = os.path.basename(mesh_size)

    plt.plot(
        iterations,
        relative_residues,
        marker='o',
        linestyle='-',
        label=f'Mesh: {mesh_name}'
    )

plt.yscale('log')  # Logarithmische Skala für Residuen
plt.xlabel(r'Iteration')
plt.ylabel(r'Rel. Res. Norm')
plt.title(r'Convergence Plot for Different Mesh Sizes')
plt.grid(True)
plt.legend()  # Legende anzeigen

# Plot speichern
output_dir_name = os.path.basename(output_dir.rstrip(os.sep))

output_file = os.path.join(
    output_dir,
    f"convergence_plot_{output_dir_name}.png"
)
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"Plot gespeichert in: {output_file}")
