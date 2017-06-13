
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

"""Given a list of ExAC alleles, a bar chart of allele frequencies is plotted, with 
subpopulations on the x axis and frequencies on the y axis"""


def plotAlleles(alleleList, filename):
    """Given an alleleList, will produce 2 bar charts rendering allele freq
    data"""
    figureNum = 0 #will allow "figure 1" and "figure 2" to be made
    for allele in alleleList: #allele is an object with attributes
        if len(allele.pops) == 5:
            color = 'cyan'
            title = '1000 Genomes Frequencies'
            chromPos = 'chrom ' + str(allele.chrom) + ', ' + 'position ' + str(allele.pos)
            plotAllele(allele, title, color, filename)
            truncPlot(allele, title, color, filename)
            for sp in allele.pops:
                print sp.namepop, sp.af
            figureNum += 1
            print 'plotdone'
        elif len(allele.pops) == 7:
            color = 'salmon'
            title = 'ExAC Frequencies'
            chromPos = 'chrom ' + str(allele.chrom)+', '+'position '+str(allele.pos)
            plotAllele(allele, title, color, filename)
            truncPlot(allele, title, color, filename)
            for sp in allele.pops:
                print sp.namepop, sp.af
            figureNum += 1
            print 'plotdone'
        else:
            print "vcf.gz file not compatible"


def plotAllele(allele, title, color, filename):
    """given an allele, will create a bar chart for its allele frequencies"""
    #pp = PdfPages(filename)
    labels, AFs = makeFrequencyVectors(allele)
    sigLine = np.full(len(labels), 0.01) #creates a line across all points of the histogram
    xAlign = np.arange(len(labels))
    fig, ax = plt.subplots()
    ax.bar(xAlign, AFs, color= color)
    ax.plot(xAlign, sigLine, color = 'black')
    ax.set_ylabel('Allele Frequencies')
    ax.set_title(title)
    ax.set_xticks(xAlign) #this is the line that made a difference for label alignment
    ax.set_xticklabels(labels)
    #ax.set_ylim((0,1))
    #pp.savefig()
    #pp.close()
    plt.show()

def truncPlot(allele, title, color, filename):
    #pp = PdfPages(filename) #do some stripping?
    labels, AFs = makeFrequencyVectors(allele)
    sigLine = np.full(len(labels), 0.01)
    xAlign = np.arange(len(labels))
    fig, (ax, ax2) = plt.subplots(2, 1, sharex=True)
    ax.plot(xAlign, sigLine, color = 'black')
    #ax2.plot(xAlign, sigLine, color = 'black')
    ax.set_title(title)
    ax.bar(xAlign, AFs, color= color)
    ax2.bar(xAlign, AFs, color = color)
    yLim = findYLim(allele) #THERE IS A CONCERN I HAVE WITH THE AUTOFIT
    ax.set_ylim(.005, 1.)
    ax2.set_ylim(0, yLim) #figure out how to get this to autofit the data through another function
    ax.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax.xaxis.tick_top()
    ax.tick_params(labeltop='off')
    ax2.xaxis.tick_bottom()

    dd = .015
    kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    ax.plot((-dd, +dd), (-dd, +dd), **kwargs)  # top-left diagonal
    ax.plot((1 - dd, 1 + dd), (-dd, +dd), **kwargs)  # top-right diagonal
    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-dd, +dd), (1 - dd, 1 + dd), **kwargs)  # bottom-left diagonal
    ax2.plot((1 - dd, 1 + dd), (1 - dd, 1 + dd), **kwargs)  # bottom-right diagonal

    ax2.set_ylabel('Allele Frequencies')
    ax2.set_xlabel('populations')
    ax.set_xticks(xAlign)  # this is the line that made a difference for label alignment
    ax.set_xticklabels(labels)

    #pp.savefig()
    #pp.close()
    plt.show()

def findYLim(allele):
    afs = []
    for subpop in allele.pops:
        afs.append(subpop.af)
    if max(afs) < 0.0004:
        yLim = max(afs) + max(afs)/2
    else:
        yLim = 0.003
    return yLim

def makePdfFilename(filename):
    newFn = filename.strip('.vcf')
    newerFn = filename.strip


def makeFrequencyVectors(allele):
    """given an allele, will return a labels and an allele frequencies vector"""
    labels = []
    AFs = []
    for sp in sorted(allele.pops, key=lambda p: (p.namepop,)): #FIX THIS, YOU DONT NEED subpopTups #keys, given one element in the list, returns value youa re supposed to sort on
        labels.append(sp.namepop)
        AFs.append(sp.af)
    return labels, AFs




