#!/usr/bin/env python
import argparse
import hashlib
import logging
import os
import sys
from pathlib import Path
from typing import Set


def read_file_list(path: Path) -> Set[Path]:
    with open(path, "r") as f:
        return {Path(line.strip()) for line in f.readlines()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputDir", help="Input directory for generating md5sums")
    parser.add_argument("-o", "--outputFile", help="Output file for md5sums")
    parser.add_argument("-f", "--keepListFilePath", help="only consider files whose paths are given in this file."
                                                         "Paths are expected to be relative to --inputDir")
    parser.add_argument("-d", "--discardListFilePath",
                        help="paths in this file are ignored, i.e. not included in final tarball. Paths are expected to be relative to --inputDir")

    args = parser.parse_args()

    logging.getLogger().setLevel(logging.INFO)

    keep_list = read_file_list(args.keepListFilePath)
    discard_list = read_file_list(args.discardListFilePath)

    lists_intersec = keep_list.intersection(discard_list)
    if lists_intersec:
        sys.exit(f"Keep list and discard list are not disjoint! Got {', '.join([str(p) for p in lists_intersec])}")

    input_dir = Path(args.inputDir)
    output_file_name = args.outputFile

    # Recurses through a directory and it's subdirectories and generates md5 hashes for each file in keeplist
    # All hashes are appended to an output file as specified with the -o flag.
    with open(output_file_name, 'w') as f_out:
        for subdir, _, files in os.walk(input_dir):
            subdir_path = Path(subdir).relative_to(input_dir)

            for file in files:
                # Don't hash the output file
                if file == output_file_name.split('/')[-1]:
                    continue

                filepath = subdir_path / file
                if filepath not in keep_list:
                    logging.info(f"Won't include file {filepath} in tarball")

                    if filepath not in discard_list:
                        # If a file is neither in the keep list nor the discard,list, fail, as it is an unexpected file.
                        # This way we are making sure we don't forget to include a newly created file into the tarball
                        sys.exit(f"Found found file {filepath} neither in keep list nor discard list. ")
                else:
                    md5hash = hashlib.md5(open(Path(subdir) / file, 'rb').read()).hexdigest()
                    f_out.write(f"{md5hash}  {filepath}\n")  # 2 whitespaces in order to be compatible with GNU md5sum


if __name__ == "__main__":
    main()
