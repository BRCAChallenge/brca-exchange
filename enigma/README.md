this folder documents the process of cleaning and processing the excel file containing the variants classified by ENIGMA

ENIGMA excel files at excel_raw_files:

raw_files/ENIGMA_for_BRCAsite_9.21.2015.xlsx
raw_files/ENIGMA_SubmissionClinVar_2016-05-31_ncbi.xlsx

1. open excel file with google spreadsheet, go to the "Variant" tab
2. delete the following rows:
    row1: ("This worksheet is REQUIRED for submission. All assertions......")
    row2: (Gene, Variant definition...)
    row4: (Optional but highly recommended....)
    row5: (Please start your submission in ....)
3. delete the columns with the following title:
    ##Local ID
    Linking ID
    Mode of Inheritance
    Citations or URLs for  clinical significance without database identifiers
    Explanation if clinical significance is other or drug response
    Drug response condition
    Functional consequence 
    Comment on functional consequence
    click '+' above this column to add more clinical significance data
    Affected status
4.  file -> download as -> Tab-separated values (tsv, .current sheet)
    file renamed and saved as raw_files/ENIGMA_variants_batch1_9_21_2015.tsv
