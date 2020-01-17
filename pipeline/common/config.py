import pandas as pd

HGVS_CDNA_DEFAULT_AC = 'hgvs_cdna_default_ac'
SYMBOL_COL = 'symbol'
SYNONYM_AC_COL = 'synonyms_ac_col'

def load_config(path):
    df = pd.read_csv(path, sep=',', header=0)
    return df.set_index(SYMBOL_COL, drop=False)
