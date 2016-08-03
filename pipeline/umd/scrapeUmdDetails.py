#!/usr/bin/env python

# convert UMD details pages to tab-sep file
# first argument is the index file
# other arguments is a list of files e.g. *.html
import glob, sys
import lxml.html as et
from collections import OrderedDict
from os.path import basename, splitext

def cleanField(str):
    " remove tab, \n \m, ---- from field "
    return str.replace("\n", "").replace("\m","").lstrip("-")

def parseIndexPage(fname):
    " return a list of dicts with header -> value "
    data = open(fname).read()
    root = et.document_fromstring(data)
    i = 0
    allData = {}
    for tEl in root.findall(".//table"):
        if i == 0:
            i += 1
            continue

        headers = []
        for tdEl in tEl.findall(".//th"):
            headers.append(tdEl.text_content())

        #print "headers", headers

        for trEl in tEl.findall(".//tr"):
            mutRow = {}
            col = 0
            for tdEl in trEl.findall(".//td"):
                #print headers, col, tdEl.text_content()
                mutRow[headers[col]] = tdEl.text_content().lstrip("-").rstrip("\n")
                col += 1
            if len(mutRow)!=0:
                allData[mutRow["UMD_id"]] = mutRow
    return allData

indexFname = sys.argv[1]
otherFnames = sys.argv[2:]

indexData = parseIndexPage(indexFname)

# make sure all index data has the same fields
indexHeaders = set(indexData.values()[0].keys())
for umdId, mutData in indexData.iteritems():
    assert(len(set(mutData) - indexHeaders) == 0)

indexHeaders = sorted(list(indexHeaders))

allKeys = OrderedDict()
allData = []
for fname in otherFnames:
    rowId = splitext(basename(fname))[0]
    sys.stderr.write("Parsing %s\n" % fname)
    data = open(fname).read()
    #data = data.replace("<table class='rien' width='80'>","")
    root = et.document_fromstring(data)
    #root = tree.getroot()
    data = {}
    # find the tables, get the first row as the keys
    # get the 2nd row as values and put everything into a dict as key->value
    # [@id="mydiv and @class!="exclass"]
    for tEl in root.findall(".//table"):
        keys = ["id"]
        vals = [rowId]
        #print "element", tEl
        #print "new table", tEl.text_content().encode("utf8")
        if tEl.get("class")=="rien":
            continue
        rowIdx = 0
        for rowEl in tEl.findall("tr"):
            skipRow = False
            ###print "row", rowIdx, rowEl.text_content(), skipRow
            for tdEl in rowEl.findall("td"):
                tdText = tdEl.text_content()
                if tdText!=None and "Impact on splicing" in tdText:
                    skipRow = True
                    #print "SKIP", rowIdx
                    break

                if tdText==None:
                    tdText=""
                else:
                    tdText = tdText.strip().encode("utf8")

                if rowIdx==0:
                    assert(tdText!="58.5")
                    # put this into a special x-link field?
                    if tdText.startswith("Other variation(s)"):
                        continue
                    keys.append(tdText)
                elif rowIdx==1:
                    vals.append(tdText)
                else:
                    assert(False)
            if skipRow:
                continue
            else:
                rowIdx+=1

        if len(keys)==len(vals):
            data.update(zip(keys, vals))
            for key in keys:
                if key not in allKeys:
                    allKeys[key]=True
    allData.append(data)
        #for key, val in data.items():
           #print key, val

# print header
headers = list(allKeys)
headers.extend(indexHeaders)
print "#"+"\t".join(headers)

for data in allData:
    if len(data)==0:
        continue
    row = []
    umdId = data["id"]
    for key in allKeys:
        row.append(cleanField(data.get(key, "")))
    # append the index data
    for indexKey in indexHeaders:
        row.append(cleanField(indexData[umdId][indexKey]))

    #row = [x.encode("utf8") for x in row]
    print "\t".join(row)
