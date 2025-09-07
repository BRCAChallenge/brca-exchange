# Here are the canonical BRCA transcripts in ENSEMBL nomenclature
BRCA1_CANONICAL = "ENST00000357654"
BRCA2_CANONICAL = "ENST00000380152"

# clinically important domain boundaries
brca1CIDomains = {"enigma": {"ring": {"domStart": 43124096,
                                      "domEnd": 43104260},
                             "brct": {"domStart": 43070966,
                                      "domEnd": 43045705}},
                  "priors": {"initiation": {"domStart": 43124096,
                                            "domEnd": 43124094},
                             "ring": {"domStart": 43124084,
                                      "domEnd": 43104875},
                             "brct": {"domStart": 43070966,
                                      "domEnd": 43045705}}}

brca2CIDomains = {"enigma": {"dnb": {"domStart": 32356433,
                                     "domEnd": 32396954}},
                  "priors": {"initiation": {"domStart": 32316461,
                                            "domEnd": 32316463},
                             "palb2": {"domStart": 32316491,
                                       "domEnd": 32319108},
                             "dnb": {"domStart": 32356433,
                                     "domEnd": 32396954},
                             "tr2/rad5": {"domStart": 32398318,
                                          "domEnd": 32398428}}}

# BRCA1/BRCA2 grey zone boundaries
greyZones = {"BRCA2": {"greyZoneStart": 32398438,
                       "greyZoneEnd": 32398488}}

# standard window sizes for splice donor and acceptor
# 3/5/18, defined by Sean Tavtigian and Michael Parsons
# also found in MaxEntScan score splice site definitions (Yeo and Burge 2004)
STD_DONOR_SIZE = 9
STD_DONOR_INTRONIC_LENGTH = 6
STD_DONOR_EXONIC_LENGTH = 3
STD_ACC_SIZE = 23
STD_ACC_INTRONIC_LENGTH = 20
STD_ACC_EXONIC_LENGTH = 3

# standard exonic portion size and de novo acceptor length
# 3/5/18, defined by Sean Tavtigian and Michael Parsons
# stdExonicPortion also found in MaxEntScan score splice site definitions (Yeo and Burge 2004)
STD_EXONIC_PORTION = 3
STD_DE_NOVO_LENGTH = 10

# standard de novo offset
# subject to change if values above change
# default value as of 3/5/18 is 7
STD_DE_NOVO_OFFSET = STD_DE_NOVO_LENGTH - STD_EXONIC_PORTION

# Canonical BRCA transcripts in RefSeq nomenclature
BRCA1_RefSeq = "NM_007294.3"
BRCA2_RefSeq = "NM_000059.3"

# Probability constants from SVT and MP valid as of 5/24/18
LOW_PROBABILITY = 0.02
LOW_SPLICING_PROBABILITY = 0.04
MODERATE_DENOVO_PROBABILITY = 0.3
MODERATE_SPLICING_PROBABILITY = 0.34
CAPPED_PROBABILITY = 0.5
HIGH_PROBABILITY = 0.64
HIGH_SPLICING_PROBABILITY = 0.97
PATHOGENIC_PROBABILITY = 0.99

# Splicing cutoff constants from SVT and MP valid as of 5/24/18
LOW_MES_CUTOFF = 6.2
HIGH_MES_CUTOFF = 8.5
MIN_REF_ALT_ZDIFF = 0.5
REF_DONOR_REFZ_CUTOFF = -1.5
REF_DONOR_HIGH_CUTOFF = 0.0
REF_DONOR_LOW_CUTOFF = -2.0
REF_ACC_REFZ_CUTOFF = -1.0
REF_ACC_HIGH_CUTOFF = 0.5
REF_ACC_LOW_CUTOFF = -1.5
DE_NOVO_DONOR_HIGH_CUTOFF = 0.0
DE_NOVO_DONOR_LOW_CUTOFF = -2.0

# summary statistics for z-score computation
BRCA_ZSCORES = {
    "donors": { "std": 2.3289956850167082, "mean": 7.9380909090909073 },
    "acceptors": { "std": 2.4336623152078452, "mean": 7.984909090909091 }
}
