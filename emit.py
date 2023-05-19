# defines logic for process to process communication in nextflow.
# Supports two modes:
# 1. by stdout: print outputs to stdout
# 2. by file: print outputs to file

import json


def emit_stdout(output):
    print(json.dumps(output))


def emit_file(output, output_file):
    with open(output_file, "w") as f:
        json.dump(output, f)


def emit(output, output_file=None):
    if output_file:
        emit_file(output, output_file)
    else:
        emit_stdout(output)
