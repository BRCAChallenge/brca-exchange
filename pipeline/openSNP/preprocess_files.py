import glob
import os


RAW_DATA = "data/datadump_opensnp_9_8_2016/"


def main():
    file_dict = file_stats()
    organize_files(file_dict)

def organize_files(file_dict):
    path = "data/organized_opensnp_filedump/"
    try: os.mkdir(path)
    except OSError: pass
    for folder in file_dict.keys():
        if "readme" not in folder:
            try: os.mkdir(path + folder)
            except OSError: pass
        else:
            pass


def file_stats():
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
    sum = 0
    for key, value in file_dict.iteritems():
        print key, ":", len(value)
        sum += len(value)
    print "sum:", sum

    return file_dict


if __name__ == "__main__":
    main()

