import subprocess
import tempfile

INPUT = "/cluster/home/mollyzhang/release1.0/merged.csv"
OUTPUT = "/cluster/home/mollyzhang/release1.0/data/VEP/vep_input.txt"


def main():
    out = subprocess.check_output(["cut", "-d,", "-f3", INPUT])
    genome_coors = out.split("\n")[1:-1]
    f = open(OUTPUT, "w")
    for genome_coor in genome_coors:
        items = genome_coor[3:].replace("-", ".").split(":")
        items = items[0:2] + ['.'] + items[2].split(">")
        f.write(" ".join(items) + "\n")
        

if __name__ == "__main__":
    main()

