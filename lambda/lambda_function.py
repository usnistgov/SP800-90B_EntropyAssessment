import subprocess
import os
import base64

def store_file(file_b64):
    """
    Takes the base64 string form of a file, writes the file to a 
    temporary file path, and then returns the file path.

    Args:
        file_b64: the base64 respresentation of the file as a string

    Returns:
        string: The location where the file was written.
    """

    # The temporary file path used
    tmp_file_path = "/tmp/input_data.bin"

    # Write the file out
    with open(tmp_file_path, "wb") as f:
        f.write(base64.b64decode(file_b64))

    # Return the file path.
    return tmp_file_path

def _run_command(command):
    """
    Takes a command and runs it, returning the response JSON

    Args:
        command: array of strings respresenting the command to be run.

    Returns:
        object: An object representing the response JSON
    """

    # Load the current state of the environment variables, then add the lib and bin directories from the 
    # 90B tests to the env LD_LIBRARY_PATH and PATH variables, respectively
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = "/var/task/lib:" + env.get("LD_LIBRARY_PATH", "")
    env["PATH"] = "/var/task/bin/:" + env.get("PATH", "")
    print(f"DEBUG: Starting command invocation: {command}")
    # Invoke the command and return the response body
    try:
        result = subprocess.run(
            command, 
            env=env,
            capture_output=True,
            text=True,
            check=True,
        )
        print("DEBUG: Command completed successfully!")
        print(result.stdout)
        return {
            "statusCode": 200,
            "cmd": command,
            "stdout": result.stdout,
            "stderr": result.stderr,
        }
    except subprocess.CalledProcessError as e:
        return {
            "statusCode": 500,
            "cmd": command,
            "returncode": e.returncode,
            "stdout": e.stdout,
            "stderr": e.stderr,
        }


def raw_noise(event, context):
    #
    # This function expects JSON of the format:
    # {
    #    "iid": true,
    #    "bits_per_sample": 8,
    #    "file_b64": "<base64 encoded file content>",
    # }
    #
    command = []
    
    if(event.get("iid") == True):
        command.append("ea_iid")
    else:
        command.append("ea_non_iid")

    #command += "-o -a -i -q"
    command.extend(["-o", "-a", "-i", "-q"])

    file_path = store_file(event.get("file_b64"))
    command.append(file_path)

    command.append(str(event.get("bits_per_sample")))

    return _run_command(command)

def restart_noise(event, context):
    #
    # This function expects JSON of the format:
    # {
    #    "iid": true,
    #    "h_i": 8,
    #    "bits_per_sample": 8,
    #    "file_b64": "<base64 encoded file content>",
    # }
    #
    command = ["ea_restart"]
    
    if(event.get("iid") == True):
        command.append("-i")
    else:
        command.append("-n")
    
    file_path = store_file(event.get("file_b64"))
    command.append(file_path)
    command.append(str(event.get("bits_per_sample")))
    command.append(str(event.get("h_i")))
    command.extend(["-q", "-o"])

    return _run_command(command)

def conditioning_component(event, context):
    #
    # This function expects JSON of the format:
    # {
    #    "iid": true,
    #    "n_in": 8,
    #    "n_out": 8,
    #    "n_w": 8,
    #    "h_in": 8,
    #    "bits_per_sample": 8,
    #    "file_b64": "<base64 encoded file content>",
    # }
    #

    # In keeping with the existing run-conditioning.sh, we assume non-vetted for all ESV CCs
    # Since we don't accept them otherwise
    command = ["ea_conditioning", "-n"]

    # Passing the usual variables
    command.append(str(event.get("n_in")))
    command.append(str(event.get("n_out")))
    command.append(str(event.get("n_w")))
    command.append(str(event.get("h_in")))

    # Provide the input file, using -i and the file path
    command.append("-i")
    file_path = store_file(event.get("file_b64"))
    command.append(file_path)

    # Disable normal output (-q), and ouput to JSON in stdout (-o, provided without a filename)
    command.extend(["-q", "-o"])
    
    return _run_command(command)