

vcf_mandatory_cols = ['contigName', 'start', 'end', 'names', 'referenceAllele', 'alternateAlleles', 'qual', 'filters']

pop_names = ['AFR', 'AMR', 'EAS', 'NFE', 'SAS']
col_prefixes = ['INFO_faf95', 'INFO_AC', 'INFO_AF', 'INFO_AN']
additional_cols = ['INFO_popmax', 'INFO_AC_popmax', 'INFO_AF_popmax', 'INFO_AN_popmax', 'INFO_vep'] + [ f'{prefix}_{p.lower()}' for p in pop_names for prefix in col_prefixes]

faf95_col_names = [f'faf95_{p.lower()}' for p in pop_names]
