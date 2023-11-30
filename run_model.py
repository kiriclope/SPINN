import subprocess
import shutil
import yaml


def run_cpp(session, bin_path="/home/leon/models/lif_cpp/bin/LifNet", conf_path="/home/leon/models/lif_cpp/conf/config_EI.yml"):
    cmd = ["screen", "-dmS", session, bin_path, conf_path]
    # cmd = ["screen", "-dmS", session, bin_path, conf_path, "&&", "notify-send", "Job Finished", session, "has finished"]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    # Use these lines to print stdout and stderr in real-time
    for line in proc.stdout:
        print(line.decode().strip())

    for line in proc.stderr:
        print("Error: " + line.decode().strip())

    # Wait for the process to terminate.
    proc.communicate()


def update_conf(conf_name, key_to_change, new_value, axis=None):
    with open(conf_name + ".yml", "r") as stream:
        try:
            conf = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    # Now the file is implemented as a Python dictionary
    if axis is None:
        conf[key_to_change] = new_value
    else:
        conf[key_to_change][axis] = new_value

    with open(conf_name + ".yml", "w") as file:
        documents = yaml.dump(conf, file)


def restore_conf(src, dest):
    try:
        shutil.copy(src, dest)
        print("File moved successfully!")
    except FileNotFoundError:
        print("Error: Source file not found!")
    except PermissionError:
        print("Error: You do not have permission to move the file!")
    except shutil.Error as e:
        print(f"Error: Failed to move the file. Reason: {e}")

