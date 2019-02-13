import pandas as pd
from collections import namedtuple
# TODO: move elsewhere?

#header = ['symbol', 'chr',  'start_hg37', 'end_hg37', 'start_hg38', 'end_hg38']

#GeneConfig = namedtuple("GeneConfig", header)

# TODO: use dataframe or some intermediary structure?

def load_config(path):
    print("loading " + path)
    df = pd.read_csv(path, sep=',', header=0)
    # TODO validate against GeneConfig?
    #return {x.gene_name: x for x in df.itertuples(False, 'GeneConfig')}
    return df.set_index('symbol', drop=False)



# tweak path
#d = load_config('/Users/marc/git/brca-exchange/pipeline/luigi/gene_config_brca_only.txt')

#print(d['BRCA2'])