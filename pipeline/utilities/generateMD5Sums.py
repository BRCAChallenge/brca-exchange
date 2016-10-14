#!/usr/bin/env python
import os
import argparse
import hashlib


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputDir", help="Input directory for generating md5sums")
    parser.add_argument("-o", "--outputFile", help="Output file for md5sums")

    args = parser.parse_args()

    output_file_name = args.outputFile
    # Recurses through a directory and it's subdirectories and generates md5 hashes for each file
    # All hashes are appended to an output file as specified with the -o flag.
    with open(output_file_name, 'w') as f_out:
        for subdir, dirs, files in os.walk(args.inputDir):
            for file in files:
                # Don't hash the output file
                if file == output_file_name.split('/')[-1]:
                    continue
                filename = os.path.join(subdir, file)
                md5hash = hashlib.md5(open(filename, 'rb').read()).hexdigest()
                f_out.write(file + ": " + md5hash + '\n')

if __name__ == "__main__":
    main()
