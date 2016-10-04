#!/usr/bin/env python

import argparse
import csv
import re
import pdb

added_data = None
diff = None
total_variants_with_changes = 0
total_variants_with_additions = 0

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
        listsAreConsistent = False
        if oldValues is None or newValues is None:
            return False
        elif re.search(",", oldValues) and re.search(",", newValues):
            oldTokens = oldValues.split(",")
            newTokens = newValues.split(",")
            numberSharedTokens = 0
            for token in oldTokens:
                if token in newTokens:
                    numberSharedTokens += 1
            if numberSharedTokens == len(newTokens) and \
                    numberSharedTokens == len(oldTokens):
                listsAreConsistent = True
        return listsAreConsistent

    def _normalize(self, value):
        """Make all values similar for improved comparison"""

        # Replace all blank values with dashes for easier comparison
        if value == "" or value is None:
            return "-"
        # Some values start with ", " which throws off the comparison -- overwrite it.
        elif value[:2] == ", ":
            return value[2:]
        # Some values end with "," which throws off the comparison -- overwrite it.
        elif value[len(value)-1] == ",":
            return value[:len(value)-1]
        else:
            return value

    def compareField(self, oldRow, newRow, field):
        """
        Compare the old and new versions of a specific field.  Report if
        the field is added, has cosmetic changes, has major changes, or
        is unchanged.
        """
        global added_data
        if field in self._newColumnsAdded:
            return "added"
        else:
            newValue = self._normalize(newRow[field])
            oldValue = self._normalize(oldRow[self._newColumnNameToOld[field]])
            if oldValue == newValue:
                return "unchanged"
            elif oldValue == "-" or oldValue in newValue:
                variant = newRow["pyhgvs_Genomic_Coordinate_38"]
                return "added data: %s | %s" % (oldValue, newValue)
            elif self._consistentDelimitedLists(oldValue, newValue):
                return "unchanged"
            elif self._makeExpectedChanges.has_key(field):
                updatedOldValue = self._normalize(self._makeExpectedChanges[field](oldValue))
                if updatedOldValue == newValue:
                    return "minor change: %s | %s" % (oldValue, newValue)
                else:
                    return "major change: %s | %s" % (oldValue, newValue)
            else:
                return "major change: %s | %s" % (oldValue, newValue)


    def compareRow(self, oldRow, newRow):
        """
        Compare the contents of an old row to a new row.  Indicate any minor
        (cosmetic) changes, major changes, or new values
        """

        global added_data
        global diff
        global total_variants_with_additions
        global total_variants_with_changes

        # Uncomment if using old data schema (e.g. pre pyhgvs_Genomic_Coordinate_38)
        columns_to_ignore = [
                            # "Genomic_Coordinate_hg36",
                            # "Genomic_Coordinate_hg37",
                            # "Genomic_Coordinate_hg38",
                            # "Hg37_Start",
                            # "Hg37_End",
                            # "Hg36_Start",
                            # "Hg36_End",
                            # "HGVS_cDNA",
                            # "HGVS_Protein"
                            ]

        # Header to group all logs the same variant
        variant_intro = "\n\n %s \n Old Source: %s \n New Source: %s \n\n" % (newRow["pyhgvs_Genomic_Coordinate_38"],
                                                                              oldRow["Source"], newRow["Source"])

        changeset = ""
        added_data_str = ""

        for field in newRow.keys():
            if field not in columns_to_ignore:
                result = self.compareField(oldRow, newRow, field)
                if re.search("major change", result):
                    result = re.sub("major change: ", "", result)
                    changeset += "%s: %s \n" % (field, result)
                if re.search("added data", result):
                    result = re.sub("added data: ", "", result)
                    added_data_str += "%s: %s \n" % (field, result)

        # If there are any changes, log them in the diff
        if len(changeset) > 0:
            diff.write(variant_intro)
            diff.write(changeset)
            total_variants_with_changes += 1

        # If there are any additions, log them in added_data
        if len(added_data_str) > 0:
            added_data.write(variant_intro)
            added_data.write(added_data_str)
            total_variants_with_additions += 1





class v1ToV2(transformer):
    # Here are columns that were renamed between the April 2016 release and the September 2016
    # release. In this dictionary, the key is the old name, and the value is the new name.
    _renamedColumns = {"SIFT_VEP": "Sift_Prediction",
                       "PolyPhen_VEP": "Polyphen_Prediction",
                       "BIC_Identifier": "BIC_Nomenclature",
                       "Pathogenicity_default": "Pathogenicity_expert",
                       "Pathogenicity_research": "Pathogenicity_all"}

    #
    # This dictionary documents and implements some expected formatting changes between the
    # April 2016 release and the September 2016 release.  For each named field, there is a
    # lambda function that if applied to the old value, would generate the equivalent new value.
    #
    _makeExpectedChanges = {
        "Synonyms": (lambda xx: re.sub("^,", "", xx)),
        "HGVS_Protein": (lambda xx: re.sub(".p.", ":p.",
                                           re.sub("$", ")",
                                                  re.sub("p.", "p.(",
                                                         re.sub("NM_000059", "NP_000050.2",
                                                                xx))))),
        "Reference_Sequence": (lambda xx: re.sub("NM_000059", "NM_000059.3",
                                                 re.sub("NM_007294", "NM_007294.3", xx))),
        "Allele_Frequency": (lambda xx: re.sub("\(ExAC", "(ExAC)", xx)),
        "Polyphen_Prediction": (lambda xx: re.sub("\(*$", "", xx)),
        "Sift_Prediction": (lambda xx: re.sub("\(*$", "", xx)),
        "Clinical_significance_citations_ENIGMA": (lambda xx: re.sub("", "-", xx)),
        "Date_last_evaluated_ENIGMA": (lambda xx: re.sub("/15$", "/2015", xx)),
        "Submitter_ClinVar": (lambda xx: re.sub("Invitae_", "Invitae", xx))
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
    parser.add_argument("--added_data", default="added_data.tsv",
                        help="Variants with data added in version 2")
    parser.add_argument("--diff", default="diff.txt",
                        help="Variant diff output file")

    args = parser.parse_args()
    v1In = csv.DictReader(open(args.v1, "r"), delimiter="\t")
    v2In = csv.DictReader(open(args.v2, "r"), delimiter="\t")
    removed = csv.DictWriter(open(args.removed, "w"), delimiter="\t",
                             fieldnames=v1In.fieldnames)
    added = csv.DictWriter(open(args.added, "w"), delimiter="\t",
                           fieldnames=v2In.fieldnames)
    global added_data
    added_data = open(args.added_data, "w")
    global diff
    diff = open(args.diff, "w")
    v1v2 = v1ToV2(v1In.fieldnames, v2In.fieldnames)
    #
    # Save the old variants in a dictionary for which the hg38 genomic
    # HGVS string is the key, and for which the value is the full row.
    oldData = {}
    for oldRow in v1In:
        oldData[oldRow["pyhgvs_Genomic_Coordinate_38"]] = oldRow
    newData = {}
    for newRow in v2In:
        newData[newRow["pyhgvs_Genomic_Coordinate_38"]] = newRow
    for oldVariant in oldData.keys():
        if not newData.has_key(oldVariant):
            removed.writerow(oldData[oldVariant])
    for newVariant in newData.keys():
        if not oldData.has_key(newVariant):
            added.writerow(newData[newVariant])
        else:
            v1v2.compareRow(oldData[newVariant], newData[newVariant])

    print "Number of variants with additions: " + str(total_variants_with_additions)
    print "Number of variants with changes: " + str(total_variants_with_changes)

if __name__ == "__main__":
    main()
