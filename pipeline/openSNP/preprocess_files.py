"""find out file statistics of opensnp data dump and seperate files into
folders based on their content"""


import glob
import os
import shutil

RAW_DATA = "data/datadump_opensnp_9_8_2016/"


def main():
    """nothing but the main function"""
    file_dict = file_stats()
    organize_files(file_dict)

def organize_files(file_dict):
    """move files into seperate folders, each folder is named after the origin
       of data, like 23andme, ancetry, fitbit, etc"""

    path = "data/organized_opensnp_filedump/"
    try:
        os.mkdir(path)
    except OSError:
        pass
    for folder in file_dict.keys():
        print "moving", folder, "files"
        if "readme" not in folder:
            try:
                os.mkdir(path + folder)
            except OSError:
                pass
            for filepath in file_dict[folder]:
                filename = filepath.split("/")[-1]
                shutil.copyfile(filepath, path+folder+"/"+filename)
        else:
            for filepath in file_dict[folder]:
                filename = filepath.split("/")[-1]
                shutil.copyfile(filepath, path + filename)
    print "files moved and organized at data/organized_opensnp_filedump/"

def file_stats():
    """find out how many files are of each type"""

    file_dict = {"picture": [],
                 "fitbit": [],
                 "23andme": [],
                 "ancestry": [],
                 "ftdna" :[],
                 "IYG": [],
                 "decodeme": [],
                 "phenotype_readme": []}
    subfix = {}
    for index, filename in enumerate(glob.glob(RAW_DATA + "*")):
        # find out file subfix distribution
        filetype = filename.split(".")[-1]
        if filetype in subfix.keys():
            subfix[filetype] += 1
        else:
            subfix[filetype] = 1

        # seperate files based on content
        found = False
        for key in file_dict.keys():
            if key in filename:
                found = True
                file_dict[key].append(filename)
        if not found:
            file_dict["phenotype_readme"].append(filename)


    print subfix
    sum_file = 0
    for key, value in file_dict.iteritems():
        print key, ":", len(value)
        sum_file += len(value)
    print "total number of files", index+1
    print "sum of files organized:", sum_file

    return file_dict


if __name__ == "__main__":
    main()

