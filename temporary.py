import sys

outfile = sys.argv[1]

with open(outfile, "w") as f:
    f.write(sys.stdin.read())
