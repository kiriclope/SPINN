import subprocess


def run_cpp(bin_path="../bin/LifNet", conf_path="../conf/config_EI.yml"):
    cmd = [bin_path, conf_path]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Use these lines to print stdout and stderr in real-time
    for line in proc.stdout:
        print(line.decode().strip())

    for line in proc.stderr:
        print("Error: " + line.decode().strip())

    # Wait for the process to terminate.
    proc.communicate()
