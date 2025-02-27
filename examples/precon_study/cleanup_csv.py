import pandas as pd
import numpy as np

# CSV-Datei laden (passe den Dateinamen bzw. Pfad an)
df = pd.read_csv("log_data_gmres_fgmres.csv", delimiter=",")

# ------------------------------
# 1. Spalte "Solved Times" in separate Spalten aufteilen
# ------------------------------
def process_solved_times(x):
    """
    - Entfernt eventuelle Anführungszeichen.
    - Teilt den String anhand von ", " in einzelne Werte.
    - Wandelt jeden Wert in einen Float um, rundet ihn auf 3 Nachkommastellen und formatiert ihn als String.
    - Falls weniger als 4 Werte vorhanden sind, wird mit "0.000" aufgefüllt.
    - Liefert eine Liste mit genau 4 Einträgen zurück.
    """
    # Entferne Anführungszeichen und splitte den String
    times = x.strip('"').split(", ")
    # Konvertiere jeden Eintrag, runde ihn und formatiere ihn als String
    times = [f"{np.round(float(t), 3):.3f}" for t in times if t != ""]
    # Wenn weniger als 4 Werte vorhanden, mit "-" auffüllen
    if len(times) < 4:
        times.extend(["-"] * (4 - len(times)))
    return times


df["Folder Name"] = df["Folder Name"].apply(lambda x: x.replace("_", r"\_"))

# Wende die Funktion an, sodass in der Spalte "Solved Times" Listen gespeichert werden
df["Solved Times"] = df["Solved Times"].apply(process_solved_times)

# Erzeuge für jeden Solved Time-Wert eine eigene Spalte (Solved Time 1, Solved Time 2, etc.)
for i in range(4):
    df[f"Solved Time {i+1}"] = df["Solved Times"].apply(lambda x: x[i])

# Entferne die ursprüngliche Spalte "Solved Times"
df.drop(columns=["Solved Times"], inplace=True)

# ------------------------------
# 2. Spalte "Iterations" in separate Spalten aufteilen
# ------------------------------
# Hier wird davon ausgegangen, dass in der Spalte "Iterations" Werte als String stehen.
iterations_split = df["Iterations"].apply(lambda x: x.strip('"').split(", "))
# Ermittle, wie viele Iterationswerte in einer Zeile maximal vorkommen
max_iter = iterations_split.apply(len).max()

# Für jede Iteration eine eigene Spalte erzeugen
for i in range(max_iter):
    df[f"Iteration {i+1}"] = iterations_split.apply(lambda x: x[i] if i < len(x) else "-")

# Entferne die ursprüngliche Spalte "Iterations"
df.drop(columns=["Iterations"], inplace=True)

# ------------------------------
# 3. Bereinigte CSV speichern
# ------------------------------
df_schur = df[df["Folder Name"].str.contains("schur", case=False, regex=False)]

# Alle übrigen Zeilen in ein separates DataFrame
df_additive = df[df["Folder Name"].str.contains("additive", case=False, regex=False)]
df_multiplicative = df[df["Folder Name"].str.contains("multiplicative", case=False, regex=False)]
df_gmres = df[df["Folder Name"].str.contains("gmres", case=False, regex=False)]
# Speichern der beiden neuen CSV-Dateien
df_schur.to_csv("log_schur.csv", index=False)
df_additive.to_csv("log_additive.csv", index=False)
df_multiplicative.to_csv("log_mult.csv", index=False)
df_gmres.to_csv("log_gmres.csv", index=False)

print("Bereinigte CSV-Dateien wurden gespeichert!")
