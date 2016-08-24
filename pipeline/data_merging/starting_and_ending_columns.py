import csv
with open('../../website/django/data/resources/aggregated.tsv') as f:
    with open('../../website/django/data/resources/written.tsv', 'w') as w:
        cntr = 0
        for line in csv.reader(f, dialect="excel-tab"):
            if cntr != 0 :
                Chrm36, StartHg36, Alt36 = line[61].split(':')
                x36, y36 = Alt36.split(">")
                EndHg36 = int(StartHg36) + len(y36)

                Chrm37, StartHg37, Alt37 = line[60].split(':')
                x37, y37 = Alt37.split(">")
                EndHg37 = int(StartHg37) + len(y37)

                Chrm38, StartHg38, Alt38 = line[59].split(':')
                x38, y38 = Alt38.split(">")
                EndHg38 = int(StartHg38) + len(y38)
                for i in line:
                    w.write(i+"\t")
                w.write(StartHg36+"\t"+str(EndHg36)+"\t"+StartHg37+"\t"+str(EndHg37)+"\t"+StartHg38+"\t"+str(EndHg38)+"\n")
            else:
                for i in line:
                    w.write(i+"\t")
                w.write("Hg36_Start\tHg36_End\tHg37_Start\tHg37_End\tHg38_Start\tHg38_End\n")
                cntr += 1
