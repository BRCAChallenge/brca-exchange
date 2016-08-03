import pandas as pd
import string_comp
import add_uniq_id as aui


ONEKG = "../allele_frequency/1000g_af.txt"

df = pd.read_csv(ONEKG, sep="\t", header=None, names=["genome", "af"])

print df.head()
print df.columns

new_columns = ["af", "genome"]

df = df[new_columns]
print df.head()

def three_difference():
    f1 = open("../result_comparison/invitae_discord", "r").read()

    f2 = open("../result_comparison/molly_discord", "r").read()

    set1 = set(f1.split("\n"))
    set2 = set(f2.split("\n"))

    print "NM_007294.3:c.5089T>C" in set1
    print "NM_007294.3:c.5089T>C" in set2
    print set1-set2