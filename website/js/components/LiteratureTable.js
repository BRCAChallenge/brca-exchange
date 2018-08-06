'use strict';

const React = require('react'),
    { Table, OverlayTrigger, Popover } = require('react-bootstrap');
    //util = require('./util'),
    //_ = require('lodash');class

const fakeLit = [
    {
        "pmid": "29684080",
        "title": "Unexpected cancer-predisposition gene variants in Cowden syndrome and Bannayan-Riley-Ruvalcaba syndrome patients without underlying germline PTEN mutations.",
        "url": "http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007352",
        "journal": "PLoS genetics",
        "year": "2018",
        "authors": "Yehia, Lamis; Ni, Ying; Sesock, Kaitlin; Niazi, Farshad; Fletcher, Benjamin; Chen, Hannah Jin Lian; LaFramboise, Thomas; Eng, Charis",
        "keywords": "Adolescent/Adult/Aged/Child/Child, Preschool/DNA Glycosylases/DNA Mutational Analysis/Female/Genes, BRCA1/Genes, BRCA2/Genetic Predisposition to Disease/Genetic Variation/Germ-Line Mutation/Hamartoma Syndrome, Multiple/Humans/Infant/Male/Middle Aged/Neoplasms, Multiple Primary/Oncogenes/PTEN Phosphohydrolase/Prospective Studies/Proto-Oncogene Proteins c-ret/Proto-Oncogene Proteins p21(ras)/Tumor Suppressor Proteins/Whole Exome Sequencing/Xeroderma Pigmentosum Group D Protein/Young Adult",
        "mentions": 1,
        "abstract": "Patients with heritable cancer syndromes characterized by germline PTEN mutations (termed PTEN hamartoma tumor syndrome, PHTS) benefit from PTEN-enabled cancer risk assessment and clinical management. PTEN-wildtype patients (~50%) remain at increased risk of developing certain cancers. Existence of germline mutations in other known cancer susceptibility genes has not been explored in these patients, with implications for different medical management. We conducted a 4-year multicenter prospective study of incident patients with features of Cowden/Cowden-like (CS/CS-like) and Bannayan-Riley-Ruvalcaba syndromes (BRRS) without PTEN mutations. Exome sequencing and targeted analysis were performed including 59 clinically actionable genes from the American College of Medical Genetics and Genomics (ACMG) and 24 additional genes associated with inherited cancer syndromes. Pathogenic or likely pathogenic cancer susceptibility gene alterations were found in 7 of the 87 (8%) CS/CS-like and BRRS patients and included MUTYH, RET, TSC2, BRCA1, BRCA2, ERCC2 and HRAS. We found classic phenotypes associated with the identified genes in 5 of the 7 (71.4%) patients. Variant positive patients were enriched for the presence of second malignant neoplasms compared to patients without identified variants (OR = 6.101, 95% CI 1.143-35.98, p = 0.035). Germline variant spectrum and frequencies were compared to The Cancer Genome Atlas (TCGA), including 6 apparently sporadic cancers associated with PHTS. With comparable overall prevalence of germline variants, the spectrum of mutated genes was different in our patients compared to TCGA. Intriguingly, we also found notable enrichment of variants of uncertain significance (VUS) in our patients (OR = 2.3, 95% CI 1.5-3.5, p = 0.0002). Our data suggest that only a small subset of PTEN-wildtype CS/CS-like and BRRS patients could be accounted for by germline variants in some of the known cancer-related genes. Thus, the existence of alterations in other and more likely non-classic cancer-associated genes is plausible, reflecting the complexity of these heterogeneous hereditary cancer syndromes."
    },
    {
        "pmid": "29566657",
        "title": "Germline breast cancer susceptibility gene mutations and breast cancer outcomes.",
        "url": "https://bmccancer.biomedcentral.com/articles/10.1186/s12885-018-4229-5",
        "journal": "BMC cancer",
        "year": "2018",
        "authors": "Wang, Yong Alison; Jian, Jhih-Wei; Hung, Chen-Fang; Peng, Hung-Pin; Yang, Chi-Fan; Cheng, Hung-Chun Skye; Yang, An-Suei",
        "keywords": "",
        "mentions": 1,
        "abstract": "<p>BACKGROUND</p> It is unclear whether germline breast cancer susceptibility gene mutations affect breast cancer related outcomes. We wanted to evaluate mutation patterns in 20 breast cancer susceptibility genes and correlate the mutations with clinical characteristics to determine the effects of these germline mutations on breast cancer prognosis.<p>METHODS</p> The study cohort included 480 ethnic Chinese individuals in Taiwan with at least one of the six clinical risk factors for hereditary breast cancer: family history of breast or ovarian cancer, young age of onset for breast cancer, bilateral breast cancer, triple negative breast cancer, both breast and ovarian cancer, and male breast cancer. PCR-enriched amplicon-sequencing on a next generation sequencing platform was used to determine the germline DNA sequences of all exons and exon-flanking regions of the 20 genes. Protein-truncating variants were identified as pathogenic.<p>RESULTS</p> We detected a 13.5% carrier rate of pathogenic germline mutations, with BRCA2 being the most prevalent and the non-BRCA genes accounting for 38.5% of the mutation carriers. BRCA mutation carriers were more likely to be diagnosed of breast cancer with lymph node involvement (66.7% vs 42.6%; P = 0.011), and had significantly worse breast cancer specific outcomes. The 5-year disease-free survival was 73.3% for BRCA mutation carriers and 91.1% for non-carriers (hazard ratio for recurrence or death 2.42, 95% CI 1.29-4.53; P = 0.013). After adjusting for clinical prognostic factors, BRCA mutation remained an independent poor prognostic factor for cancer recurrence or death (adjusted hazard ratio 3.04, 95% CI 1.40-6.58; P = 0.005). Non-BRCA gene mutation carriers did not exhibit any significant difference in cancer characteristics or outcomes compared to those without detected mutations. Among the risk factors for hereditary breast cancer, the odds of detecting a germline mutation increased significantly with having bilateral breast cancer (adjusted odds ratio 3.27, 95% CI 1.64-6.51; P = 0.0008) or having more than one risk factor (odds ratio 2.07, 95% CI 1.22-3.51; P = 0.007).<p>CONCLUSIONS</p> Without prior knowledge of the mutation status, BRCA mutation carriers had more advanced breast cancer on initial diagnosis and worse cancer-related outcomes. Optimal approach to breast cancer treatment for BRCA mutation carriers warrants further investigation."
    },
];

function limitFieldLength(string, maxLength) {
    if (string.length <= maxLength) {
        return string;
    }

    let popper = (<Popover>{string}</Popover>);
    return (
        <span>
            {string.substring(0, maxLength)}
            <OverlayTrigger placement='bottom' overlay={popper}><span>...</span></OverlayTrigger>
        </span>
    );
}

function formatKeywords(keywords) {
    return keywords.split(",").map(elem => (<div>{elem.replace(/\//g, " / ")}</div>));
}

function LiteratureRow(props) {
}

class LiteratureTable extends React.Component {
    constructor(props) {
        super(props);
    }

    render() {
        let litRows = fakeLit.map(({title, authors, journal, year, keywords, pmid}) => (
            <tr key={`lit_${pmid}`}>
                <td><div>{limitFieldLength(title, 100)}</div></td>
                <td><div>{authors}</div></td>
                <td><div>{journal}</div></td>
                <td><div>{year}</div></td>
                <td><div>{formatKeywords(keywords)}</div></td>
                <td><div><a href={`https://www.ncbi.nlm.nih.gov/pubmed/${pmid}`} target='_blank'>{pmid}</a></div></td>
            </tr>
        ));
        return (
            <div>
                <h3>BRCA Exchange Literature Search Results</h3>
                <Table className='variant-literature nopointer' responsive bordered>
                    <thead>
                        <tr className='active'>
                            <th>Title</th>
                            <th>Authors</th>
                            <th>Journal</th>
                            <th>Date</th>
                            <th>Keywords</th>
                            <th>PMID</th>
                        </tr>
                    </thead>
                    <tbody>
                        {litRows}
                    </tbody>
                </Table>
            </div>
        );
    }
}

export default LiteratureTable;
