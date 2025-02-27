import yaml
import sys
import os
import json

def load_yaml(file_path):
    """Lädt die YAML-Datei."""
    with open(file_path, 'r') as file:
        return yaml.safe_load(file)

def save_yaml(data, file_path):
    """Speichert die YAML-Daten in einer Datei."""
    with open(file_path, 'w') as file:
        yaml.safe_dump(data, file, default_flow_style=False, sort_keys=False)  # Sort keys optional
    print(f"Die Datei wurde gespeichert unter: {file_path}")

def update_yaml(data, updates):
    """Aktualisiert die YAML-Daten basierend auf den Updates."""
    for key, value in updates.items():
        # Unterpfade wie petsc_options.ksp_type handhaben
        keys = key.split('.')
        ref = data
        for k in keys[:-1]:
            ref = ref.setdefault(k, {})  # Gehe tiefer in die Struktur oder erstelle Schlüssel
        ref[keys[-1]] = value  # Setze den Wert
    return data

def main():
    # Überprüfe die Anzahl der Argumente
    if len(sys.argv) < 4:
        print("Nicht genügend Argumente übergeben. Benutzung:")
        print("python3 update_yml_v2.py <original_yaml> <output_file> <updates_json>")
        sys.exit(1)

    # Argumente vom Bash-Skript
    original_file = sys.argv[1]
    output_file = sys.argv[2]
    updates_json = sys.argv[3]

    # YAML-Datei laden
    yaml_data = load_yaml(original_file)

    # Updates anwenden
    updates = json.loads(updates_json)  # Konvertiere den JSON-String in ein Dictionary
    updated_yaml = update_yaml(yaml_data, updates)

    # Neue Datei speichern
    save_yaml(updated_yaml, output_file)

if __name__ == "__main__":
    main()