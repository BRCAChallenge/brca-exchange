*In Silico* Probabilities of Pathogenicity can be used to predict whether a variant will have a pathogenic effect. Calculations do not require clinical data for every variant; thus it can be applied for almost every variant based on type or location, regardless of how rare the variant is or what variant-level clinical data is available.  Though each calculation does not require clinical data, the assigned probabilities are based on bioinformatic categories that were thoroughly developed and calibrated using clinical data sets (Tavtigian 2008). **Please consider these probabilities with caution, accounting for the amount of clinical data that is or is not available for the variant in question.** Variants with clinical data will have ENIGMA, ClinVar, or Multifactorial Likelihood Analysis data on the Variant Details Page, all of which provide clinical information more directly.

A variant’s risk can be associated with impact on protein translation or interference with mRNA splicing. The [*In Silico*](https://en.wikipedia.org/wiki/In_silico) Probability of Pathogenicity accounts for both mechanisms, and uses the type and the location of the variant to make predictions. The higher of two estimations, the protein-level estimation or the splicing-level estimation, is ultimately assigned as the *In Silico* Probability of Pathogenicity. By choosing the higher estimation, our model ensures that pathogenicity is predicted conservatively, assigning risk according to highest possible variant impact.


<figure style="width: 80%; margin: 0 auto 1em auto;">
    <img src="decision.png" />
    <figcaption>In this example, the Splicing-level estimation is 0.3, which is greater than the Protein-level estimation 0.02. Thus, the *In Silico* Probability of Pathogenicity is assigned as 0.3 due to predicted splicing-level impact. In other words, splicing impact introduces the most risk for this variant.</figcaption>
</figure>

In the tile, an *In Silico* Probability of Pathogenicity is displayed on a scale of 0.00 to 1.00. An arrow and a line indicates the probability, with red indicating a higher probability of impact and blue a lower. It also displays variant location, variant type, and a link to the [HCI Breast Cancer Genes Prior probabilities website](http://priors.hci.utah.edu/PRIORS/), which provides more information about these probabilities.

The protein-level and splicing-level estimation nested tiles can be expanded by clicking on the header. For the protein-level and splicing-level nested tiles, probabilities are assigned in order to estimate a variant’s predicted impact on both protein translation and mRNA splicing respectively.

The Protein-level Estimation checks whether the variant is in a Clinically Important Functional Domain (CI Domain) and the predicted severity of the variant using Align GV-GD scores. Align GV-GD scores assess the difference in physicochemical properties for the amino acid change and the level of conservation across species.

A splicing-level estimation is also calculated. A variant can impact splicing in two ways: either it alters an existing (wild-type) splice site, or it creates a new (de novo) splice site that is used instead of the wild-type site. A splice site can refer to either a splice donor or splice acceptor site. The splicing-level estimation calculates probabilities for splice likelihood on wild-type donors, wild-type acceptors, and de novo donors. Wild-type sites may be altered, and de novo donors may impact splicing by creating new splice sites. Probabilities for each type of splice site are displayed across the three tabs in the nested tile. Because these probabilities are often mutually exclusive, it is common for at least one value to be Not Applicable (N/A).

In general, green, yellow, and red indicate low, moderate, and high probability of pathogenicity respectively. However, there are some specific color schemes used in each splicing-level tab. For variants in wild-type donor and acceptor regions, a blue row labeled “Improved” indicates that the wild-type splice site score was improved by the variant. A gray row labeled “Outside consensus region” indicates that the variant is not in a splice site. Both of these colors have low probabilities associated with them. In the De Novo Donor tab, however, increased splice likelihood actually heightens risk of pathogenicity; increased likelihood of de novo splicing increases the likelihood for altered splicing. Thus, in the case of De Novo Donors, a red row labeled “Increased” is used to indicate that the variant increases probability of pathogenicity.

Please be reminded that although these estimations are well calibrated using clinical data sets, the estimations do not directly utilize clinical evidence collected from patients with the specific variant in question.

For more information about these probabilities, please visit our resource detailing [Scoring Systems for *In Silico* Probabilities of Pathogenicity](/about/insilicoScoring).

##### References
1. Vallée, MP. 2016. [PMID 26913838](https://www.ncbi.nlm.nih.gov/pubmed/26913838).

