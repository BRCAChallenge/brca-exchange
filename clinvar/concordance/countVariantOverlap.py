#!/usr/bin/env python
"""
countVariantOverlap: count how often the same variant is reported by 
different submittors, and how concordant those submitters are or are not.
"""
import argparse
import csv
import re

class submitterPair:

    def __init__(self):
        self.totalVariants = 0
        self.conflictingVariants = 0

class variantData:

    def __init__(self, id, submitter, clinicalSignificance, row):
        self.id = id
        self.submitters = list()
        self.pathogenicSubmissions = list()
        self.nonPathogenicSubmissions = list()
        self.rows = list()
        self.newSubmission(submitter, clinicalSignificance, row)

    def isPathogenic(self, clinicalSignificance):
        if re.search(r'(?i)Pathogenic', clinicalSignificance):
            return True
        else:
            return False

    def newSubmission(self, submitter, clinicalSignificance, row):
        self.submitters.append(submitter)
        self.rows.append(row)
        if self.isPathogenic(clinicalSignificance):
            self.pathogenicSubmissions.append(submitter)
        else:
            self.nonPathogenicSubmissions.append(submitter)

    def hasDiscordantSubmissions(self):
        if len(self.pathogenicSubmissions) != len(self.submitters) and \
           len(self.nonPathogenicSubmissions) != len(self.submitters):
            return True
        else:
            return False

    def report(self, verbose=False):
        if self.hasDiscordantSubmissions() or verbose:
            print self.id, "DISCORDANT", 
            print "pathogenic submissions:", self.pathogenicSubmissions,
            print "nonpathogenic submission", self.nonPathogenicSubmissions
    
def headerLineToColumnIndex(contents):
    """Given a line of header information, return a dictionary indicating
    which column is in which position"""
    columnIndex = dict()
    for ii in range(0, len(contents)):
        columnIndex[contents[ii]] = ii
    return columnIndex

def addNewSubmitter(submitterCounts, submitter):
    """Add a new entry in this dict of dicts for this new submitter"""
    if len(submitterCounts) == 0:
        # Here is the case where the first submitter is being entered
        submitterCounts[submitter] = dict()
        submitterCounts[submitter][submitter] = list()
    else:
        submitterCounts[submitter] = dict()
        submitterCounts[submitter][submitter] = list()
        for existingSubmitter in submitterCounts.keys():
            submitterCounts[existingSubmitter][submitter] = list()
            submitterCounts[submitter][existingSubmitter] = list()
    return submitterCounts

def updateCounts(currentCounts, submitterA, submitterB, row):
    """Update the number of counts for variants assessed (possibly differently)
    by two different submitters"""
    currentCounts[submitterA][submitterB].append(row)
    currentCounts[submitterB][submitterA].append(row)
    return currentCounts

def conflictingAssessment(assessmentA, assessmentB):
    """Return True if one assessment is pathogenic or likely pathogenic, and 
    the other is either uncertain or benign or likely benign"""

    if assessmentA == "Benign" or assessmentA == "Likely benign" or assessmentA == "Uncertain significance":
        if assessmentB == "Pathogenic" or assessmentB == "Likely Pathogenic":
            return True
    if assessmentA == "Pathogenic" or assessmentA == "Likely pathogenic":
        if assessmentB == "Benign" or assessmentB == "Likely benign" or assessmentB == "Uncertain significance":
            return True
    return False

def printTotalCounts(variantsBySubmitter, submitters, sharedVariants, columnIndex):
    """Output the total number of variants per submitter"""
    print "\nSubmitter\tSubmissions\tShared Submissions"
    for submitterA in submitters:
        sharedSubmissionsThisSubmitter = dict()
        for submitterB in submitters:
            if submitterB != submitterA:
                for submission in sharedVariants[submitterA][submitterB]:
                    variantName = submission[columnIndex["HGVS"]]
                    sharedSubmissionsThisSubmitter[variantName] = 1
        totalSharedSubmissions = len(sharedSubmissionsThisSubmitter.keys())
        print "%s\t%d\t%d" % (submitterA, variantsBySubmitter[submitterA],
                              totalSharedSubmissions)


def printJointCounts(sharedVariants, conflictingVariants, submitters):
    """Output a matrix of joint counts by submitter A and submitter B"""
    print "\nShared Variants by Submitters: Total (Discordant)"
    for submitter in submitters:
        print "\t%s" % (submitter.split()[0]),
    print
    for submitterA in submitters:
        print submitterA.split()[0],
        for submitterB in submitters:
            print "\t%d (%d)" % (len(sharedVariants[submitterA][submitterB]),
                                 len(conflictingVariants[submitterA][submitterB])),
        print

def printDiscordantVariants(variants):
    """Output a summary of all conflicting variants, with the labs that 
    submitted them and their assessment"""
    print "\nAll discordant variants\n";
    for variant in sorted(variants.keys()):
        variants[variant].report()




def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("inputFile")
    parser.add_argument("-d", "--debug", default=False)
    args = parser.parse_args()

    variantsBySubmitter = dict()
    variants = dict()
    sharedVariants = dict()
    confictingVariantsBySubmitter = dict()

    reader = csv.reader(open(args.inputFile), delimiter='\t')
    columnIndex = headerLineToColumnIndex(reader.next())
    for row in reader:
        variant = row[columnIndex["HGVS"]]
        submitter = row[columnIndex["Submitter"]]
        submitter = re.sub("(\s)*$", "", submitter)  # remove any trailing whitespace
        clinicalSignificance = row[columnIndex["ClinicalSignificance"]]
        clinicalSignificance = re.sub("(\s)*$", "", clinicalSignificance)  
        #
        # If this is a new submitter, give this submitter a slot in the
        # count data.
        if not variantsBySubmitter.has_key(submitter):
            sharedVariants = addNewSubmitter(sharedVariants, submitter)
            confictingVariantsBySubmitter = addNewSubmitter(confictingVariantsBySubmitter,
                                                            submitter)
            variantsBySubmitter[submitter] = 0
        variantsBySubmitter[submitter] += 1
        #
        # If this variant hasn't been seen before, record seeing it for 
        # the first time.  Else, it's been submitted by another lab,
        # and then update the counts for sharedVariants and for
        # confictingVariantsBySubmitter if appropriate.
        if not variants.has_key(variant):
            variants[variant] = variantData(variant, submitter, clinicalSignificance, row)
        else:
            for previousRow in variants[variant].rows:
                prevSubmitter = previousRow[columnIndex["Submitter"]]
                prevSubmitter = re.sub("(\s)*$", "", prevSubmitter) 
                prevSignificance = previousRow[columnIndex["ClinicalSignificance"]]
                prevSignificance = re.sub("(\s)*$", "", prevSignificance)  
                sharedVariants = updateCounts(sharedVariants, submitter, 
                                              prevSubmitter, row)
                if conflictingAssessment(clinicalSignificance, 
                                         prevSignificance):
                    confictingVariantsBySubmitter = updateCounts(confictingVariantsBySubmitter,
                                                                 submitter, prevSubmitter, row)
                    if args.debug:
                        print("%s: conflicting assessments of %s (%s) and %s (%s)" 
                              % (variant, clinicalSignificance, submitter, prevSignificance,
                                 prevSubmitter))
            variants[variant].newSubmission(submitter, clinicalSignificance, row)

    submitters = sorted(variantsBySubmitter.keys())
    printDiscordantVariants(variants)
    printTotalCounts(variantsBySubmitter, submitters, sharedVariants, 
                     columnIndex)
    printJointCounts(sharedVariants, confictingVariantsBySubmitter, submitters)


if __name__ == "__main__":
    # execute only if run as a script
    main()
