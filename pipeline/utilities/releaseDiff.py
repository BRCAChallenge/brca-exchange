#!/usr/bin/env python

import argparse
import csv
import re


class transformer(object):
    """
    Make the expected changes to update data from one version to another
    """
    _renamedColumns = {}

    _makeExpectedChanges = {}

    def __init__(self, oldColumns, newColumns):
        (self._oldColumnsRemoved, self._newColumnsAdded, 
         self._newColumnNameToOld) = self._mapColumnNames(oldColumns, 
                                                          newColumns)


    def _mapColumnNames(self, oldColumns, newColumns):
        """
        Given a set of the old column names and a set of the new column 
        names, plus a hash mapping old to new names, return three things:
        - a list of new columns not represented in the old
        - a list of old columns not represented in the new
        - a dictionary describing columns represented in both old and
          new, and mapping the new names to the corresponding old names.
          """
        newToOldNameMapping = {}
        oldColumnsRemoved = []
        newColumnsAdded = []
        for ocol in oldColumns:
            if ocol in newColumns:
                newToOldNameMapping[ocol] = ocol
            elif self._renamedColumns.has_key(ocol):
                newToOldNameMapping[self._renamedColumns[ocol]] = ocol
            else:
                oldColumnsRemoved.append(ocol)
        for ncol in newColumns:
            if not ncol in oldColumns:
                if not ncol in self._renamedColumns.values():
                    newColumnsAdded.append(ncol)
        return (oldColumnsRemoved, newColumnsAdded, newToOldNameMapping)

    def _consistentDelimitedLists(self, oldValues, newValues):
        """Determine if the old and new values are comma-separated
        lists in which the same elements have been assembled in
        a differnt order
        """
        listsAreConsisten = False
        oldTokens = oldValues.split(",")
        newTokens = newValues.split(",")
        numberSharedTokens = 0
        for token in oldTokens:
            if token in newTokens:
                numberSharedTokens += 1
        if numberSharedTokens == len(newValues) and \
                numberSharedTokens == len(oldValues):
            listsAreConsistent = True
        return listsAreConsistent

    def compareField(self, oldRow, newRow, field):
        """
        Compare the old and new versions of a specific field.  Report if
        the field is added, has cosmetic changes, has major changes, or
        is unchanged.
        """
        if field in self._newColumnsAdded:
            return "added"
        else:
            newValue = newRow[field]
            oldValue = oldRow[self._newColumnNameToOld[field]]
            if oldValue == newValue:
                return "unchanged"
            elif self._consistentDelimitedLists(oldValue, newValue):
                return "unchanged"
            elif self._makeExpectedChanges.has_key(field):
                updatedOldValue = self._makeExpectedChanges[field](oldValue)
                if updatedOldValue == newValue:
                    return "minor change: %s %s" % (oldValue, newValue)
                else:
                    return "major change: %s %s" % (oldValue, newValue)
            else:
                return "major change: %s %s" % (oldValue, newValue)
        

    def compareRow(self, oldRow, newRow):
        """
        Compare the contents of an old row to a new row.  Indicate any minor
        (cosmetic) changes, major changes, or new values
        """
        for field in newRow.keys():
            result = self.compareField(oldRow, newRow, field)
            if re.search("major change", result):
                print field, result
            
        


class v1ToV2(transformer):
    _renamedColumns = {"SIFT_VEP" : "Sift_Prediction",
                       "PolyPhen_VEP" : "Polyphen_Prediction", 
                       "BIC_Identifier" : "BIC_Nomenclature", 
                       "Pathogenicity_default" : "Pathogenicity_expert",
                       "Pathogenicity_research" : "Pathogenicity_all" }


    _makeExpectedChanges = {
        "Synonyms" : (lambda xx: re.sub("^,", "", xx))
        "HGVS_Protein" : (lambda xx: re.sub("$", ")", re.sub("p.", "p.(", xx)))
        }



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--v2", default="built.tsv", 
                        help="File with version 2 of the data")
    parser.add_argument("--v1", default="aggregated.tsv",
                        help="File with Version 1 of the data")
    parser.add_argument("--removed", default="removed.tsv",
                        help="Variants from the old file that were removed in the new file")
    parser.add_argument("--added", default="added.tsv",
                        help="Variants that were added in the new file")

    args = parser.parse_args()
    v1In = csv.DictReader(open(args.v1, "r"), delimiter="\t")
    v2In = csv.DictReader(open(args.v2, "r"), delimiter="\t")
    removed = csv.DictWriter(open(args.removed, "w"), delimiter="\t",
                             fieldnames=v1In.fieldnames)
    added = csv.DictWriter(open(args.added, "w"), delimiter="\t",
                           fieldnames=v2In.fieldnames)
    v1v2 = v1ToV2(v1In.fieldnames, v2In.fieldnames)
    #
    # Save the old variants in a dictionary for which the hg38 genomic
    # HGVS string is the key, and for which the value is the full row.
    oldData = {}
    for oldRow in v1In:
        oldData[oldRow["Genomic_Coordinate_hg38"]] = oldRow
    newData = {}
    for newRow in v2In:
        newData[newRow["Genomic_Coordinate_hg38"]] = newRow
    for oldVariant in oldData.keys():
        if not newData.has_key(oldVariant):
            removed.writerow(oldData[oldVariant])
    for newVariant in newData.keys():
        if not oldData.has_key(newVariant):
            added.writerow(newData[newVariant])
        else:
            v1v2.compareRow(oldData[newVariant], newData[newVariant])
    

if __name__ == "__main__":
    #print "hello world"
    main()


