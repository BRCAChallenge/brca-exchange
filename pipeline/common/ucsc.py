import requests
import json


GENE_SYMBOL_COLUMN=18
NCBI_ID_COLUMN=21
MANE_STAT_COLUMN=24
MANE_SELECT="MANE Select"

API_URL="https://api.genome.ucsc.edu/getData/track"

def track_data(track, genome, maxItemsOutput=1000000):
    """
    Use the genome browser API to retrieve the track data for the track and genome of interest.
    There are additional parameters supported to filter by genomic coordinates, not used here.
    See https://genome.ucsc.edu/goldenpath/help/api.html
    """
    # Set up the query
    params = {
        "genome": genome,
        "track": track,
        "maxItemsOutput": maxItemsOutput,
        "jsonOutputArrays": 1
    }
    response = requests.get(API_URL, params=params)
    if response.status_code == 200:
        #
        # Search for the row with the maneStat field of 'MANE Select'
        # and the desired gene symbol.  
        data = response.json()
        return(data)
    return(None)


def symbols_to_mane_transcripts(gene_symbols):
    """
    Given a gene symbol, query the genome browser's API to
    retrieve the NCBI (RefSeq) accession its the MANE SELECT transcript
    """
    mane_select_transcripts = dict()
    for symbol in gene_symbols:
        mane_select_transcripts[symbol] = None
    data = track_data("mane", "hg38")
    if data is not None:
        for row in data["mane"]:
            if row[MANE_STAT_COLUMN] == MANE_SELECT:
                for symbol in gene_symbols:
                    if row[GENE_SYMBOL_COLUMN].upper() == symbol.upper():
                        mane_select_transcripts[symbol] = row[NCBI_ID_COLUMN]
    return(mane_select_transcripts)

    
