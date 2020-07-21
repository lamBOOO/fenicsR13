import yaml
import os
import subprocess


# setup
input_file = "input.yml"
output_folder = "inputs/"

# setup parameters
kn_range = [(i + 1) / 10 for i in range(10)]

# read template file
try:
    with open(input_file, "r") as stream:
        data = yaml.safe_load(stream)
except FileNotFoundError:
    print(f"[Error] The input file {input_file} is not found.\n")
    raise

# create output folder if not exists
os.makedirs(os.path.dirname(output_folder), exist_ok=True)

# write all input_files
for kn1 in kn_range:
    data["regs"][4000]["kn"] = kn1
    data["output_folder"] = "{}_".format(kn1)

    output_path = output_folder + "{}.yml".format(kn1)
    with open(output_path, "w") as stream:
        yaml.dump(data, stream)

# run study
for kn1 in kn_range:
    subprocess.check_call([
        "fenicsR13", output_folder + "{}.yml".format(kn1)
    ], cwd=".")
