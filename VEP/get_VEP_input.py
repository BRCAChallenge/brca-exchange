import subprocess
import tempfile

INPUT = "/cluster/home/mollyzhang/release1.0/merged.csv"
OUTPUT = "/cluster/home/mollyzhang/release1.0/data/VEP/vep_input.txt"


def main():
    f = tempfile.tempfile()
    subprocess.call(["cut", "-d,", "-f3", INPUT], stdout=f)
    for line in f:
        print line.strip()

if __name__ == "__main__":
    main()

