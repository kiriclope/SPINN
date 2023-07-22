import subprocess


def run_cpp():
    cmd = ["./bin/LifNet", "./conf/config_EI.yml"]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Use these lines to print stdout and stderr in real-time
    for line in proc.stdout:
        print(line.decode().strip())

    for line in proc.stderr:
        print("Error: " + line.decode().strip())

    # Wait for the process to terminate.
    proc.communicate()


run_cpp()
