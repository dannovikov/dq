import sys
from emit import emit


def clean_nucleotides(id, seq):
    iupac_code = set(["A", "C", "G", "T", "U", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N", "-", "*"])
    ambig_code = set(["R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N"])

    flags = {"unknownChars": 0, "ambig5percent": 0, "nString": 0}

    consecutive_ns = 0
    ambig_count = 0
    clean_seq = ""

    for nucl in seq:
        if nucl not in iupac_code:
            flags["unknownChars"] = 1
        else:
            clean_seq += nucl

        if nucl in ambig_code:
            ambig_count += 1

        if nucl == "N":
            consecutive_ns += 1
            if consecutive_ns >= 5:
                flags["nString"] = 1
        else:
            consecutive_ns = 0

    if ambig_count / len(seq) >= 0.05:
        flags["ambig5percent"] = 1

    emit({"id": id, "clean_seq": clean_seq, "flags": flags})


if __name__ == "__main__":
    clean_nucleotides(*(open(sys.argv[1]).read().splitlines()))
