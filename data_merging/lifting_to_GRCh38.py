"""this script do the conversion of enigma variant from GRCh37 to 38"""
import pandas as pd
from pprint import pprint


ENIGMA_37 = "/hive/groups/cgl/brca/enigma_variants_9-29-2015.tsv"
ENIGMA_38 = "/hive/groups/cgl/brca/enigma_variants_GRCh38_2-27-2016.tsv"
CLINVAR_38 = "/hive/groups/cgl/brca/release1.0/data/ClinVar/ClinVarBrca.txt"



def main():
    clinvar_df = pd.DataFrame.from_csv(CLINVAR_38, sep="\t", index_col=None)
    enigma_in_clinvar = clinvar_df[clinvar_df['Submitter'].str.contains("ENIGMA")]
    lift_dict = enigma_in_clinvar[["SCV","Genomic_Coordinate"]].set_index("SCV").to_dict()['Genomic_Coordinate']
    enigma_df = pd.DataFrame.from_csv(ENIGMA_37, sep="\t", index_col=None)
    enigma_df["Genomic_Coordinate_38"] = ""
    for index, row in enigma_df.iterrows():
        if row["ClinVarAccession"] not in lift_dict.keys():
            raise Exception("bad time") 
        row["Genomic_Coordinate_38"] = lift_dict[row["ClinVarAccession"]]
    enigma_38 = enigma_df.drop("Genomic_Coordinate", axis=1)
    enigma_38.rename(columns={'Genomic_Coordinate_38': 'Genomic_Coordinate'}, inplace=True)
    cols = list(enigma_38.columns.values)
    cols = [cols[0]] + [cols[-1]] + cols[1:-1]
    enigma_38 = enigma_38[cols]
    enigma_38.to_csv(ENIGMA_38, sep="\t", index=False)

if __name__ == "__main__":
    main()
