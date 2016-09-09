import glob

RAW_DATA = "data/datadump_opensnp_9_8_2016/"


def main():
    file_stats()



def file_stats():
    file_dict = {"picture": [],
                 "fitbit": [],
                 "23andme": [],
                 "ancestry": [],
                 "ftdna" :[],
                 "IYG": [],
                 "decodeme": []}
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
            print filename


    print subfix
    sum = 0
    for key, value in file_dict.iteritems():
        print key
        print len(value)
        sum += len(value)
    print sum




if __name__ == "__main__":
    main()

