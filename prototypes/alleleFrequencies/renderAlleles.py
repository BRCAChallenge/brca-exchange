import vcfAlleleLibrary
import matplotlib.pyplot as plt
import numpy as np

"""Given a list of ExAC alleles, a bar chart of allele frequencies is plotted, with 
subpopulations on the x axis and frequencies on the y axis"""


def plotAlleles(alleleList):
    """Given an alleleList, will produce 2 bar charts rendering allele freq
    data"""
    figureNum = 0 #will allow "figure 1" and "figure 2" to be made
    for allele in alleleList: #allele is an object with attributes
        figureNum += 1
        if len(allele.pops) == 5:
            color = 'cyan'
            title = '1000 Genomes Frequencies'
            plotAllele(allele, title, color)
            truncPlot(allele, title, color)
        if len(allele.pops) == 7:
            color = 'salmon'
            title = 'ExAC Frequencies'
            plotAllele(allele, title, color)
            truncPlot(allele, title, color)
        else:
            print 'Vcf.gz not recognized'


def plotAllele(allele, title, color):
    """given an allele, will create a bar chart for its allele frequencies"""
    labels, AFs = makeFrequencyVectors(allele)
    sigLine = np.full(len(labels), 0.05)
    xAlign = np.arange(len(labels))
    fig, ax = plt.subplots()
    ax.bar(xAlign, AFs, color= color)
    ax.plot(xAlign, sigLine, color = 'black')
    ax.set_ylabel('Allele Frequencies')
    ax.set_title(title)
    ax.set_xticks(xAlign) #this is the line that made a difference for label alignment
    ax.set_xticklabels(labels)
    #ax.set_ylim((0,1))
    plt.show()

def truncPlot(allele, title, color):
    labels, AFs = makeFrequencyVectors(allele)
    sigLine = np.full(len(labels), 0.05)
    xAlign = np.arange(len(labels))
    fig, (ax, ax2) = plt.subplots(2, 1, sharex=True)
    ax.plot(xAlign, sigLine, color = 'black')
    #ax2.plot(xAlign, sigLine, color = 'black')
    ax.set_title(title)
    ax.bar(xAlign, AFs, color= color)
    ax2.bar(xAlign, AFs, color = color)
    ax.set_ylim(.02, 1.)
    ax2.set_ylim(0, .001) #figure out how to let this autofit
    ax.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax.xaxis.tick_top()
    ax.tick_params(labeltop='off')
    ax2.xaxis.tick_bottom()

    d = .015
    kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    ax.plot((-d, +d), (-d, +d), **kwargs)  # top-left diagonal
    ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal
    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
    ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

    ax2.set_ylabel('Allele Frequencies')
    ax2.set_xlabel('populations')
    ax.set_xticks(xAlign)  # this is the line that made a difference for label alignment
    ax.set_xticklabels(labels)
    plt.show()

def findYLimit(allele):
    for subpop in allele.pops:
        sps = vars(subpop)

def plotTgAllele(allele, figureNum, color):
    labels, AFs = makeFrequencyVectors(allele)
    sigLine = [0.05, 0.05, 0.05, 0.05, 0.05]
    xAlign = np.arange(len(labels))
    fig, ax = plt.subplots()
    ax.bar(xAlign, AFs, color = 'cyan')
    ax.plot(xAlign, sigLine)
    ax.set_ylabel('Allele Frequency')
    ax.set_title('1000Genomes Frequencies')
    ax.set_xticks(xAlign)
    ax.set_xticklabels(labels)
    plt.show()


def makeFrequencyVectors(allele):
    """given an allele, will return a labels and an allele frequencies vector"""
    labels = []
    AFs = []
    subpopTups = []
    for sp in sorted(allele.pops, key=lambda p: (p.namepop,)): #FIX THIS, YOU DONT NEED subpopTups #keys, given one element in the list, returns value youa re supposed to sort on
        labels.append(sp.namepop)
        AFs.append(sp.af)
    return labels, AFs

# def calcAF(subpop):
#     """Given a subpop, will generate a tuple of (Name, AF)"""
#     ac = float(subpop.ac)
#     an = float(subpop.an)
#     af = (ac/an)
#     spTup = (subpop.namepop, af)
#     print spTup
#     return spTup


