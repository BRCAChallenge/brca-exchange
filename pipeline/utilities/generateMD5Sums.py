#!/usr/bin/env python
import os
import argparse
import hashlib


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputDir", help="Input directory for generating md5sums")
    parser.add_argument("-o", "--outputFile", type=argparse.FileType('w'),
                        help="Output file for md5sums")

    args = parser.parse_args()

    # Recurses through a directory and it's subdirectories and generates md5 hashes for each file
    # All hashes are appended to an output file as specified with the -o flag.
    for subdir, dirs, files in os.walk(args.inputDir):
        for file in files:
            filename = os.path.join(subdir, file)
            md5hash = hashlib.md5(open(filename, 'rb').read()).hexdigest()
            args.outputFile.write(file + ": " + md5hash + '\n')

if __name__ == "__main__":
    main()
