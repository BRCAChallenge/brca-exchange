#### Computational Predictions

Computational predictions are algorithmic estimates of the impact of the variant on BRCA protein function or abundance, or mRNA splicing. These predictors can inform variant interpretation, but should not be used in place of expert interpretation.

The Computational Prediction tile contains the results of computational methods that were selected by the ENIGMA Consortium as informative, and serve as inputs to the classificaiton rules of the ClinGen ENIGMA/BRCA Variant Curation Expert Panel.  Along with the results of these machine learning methods, the tile contains information on the type of variant and whether or not it overlaps a clinically-important protein domain.  These are contextual attributes of the variant which inform curators on which computational prediction methods are relevant

For each computational predictor, we present a set of general fields that are common to all variants (**Method**, **Description**, **Publication**) as well as fields that are specific to each prediction method.

----
#### Variant Type ((varType))

Indicates the type of impact that the variant has on the reference sequence: a substition, insertion, deletion, or a change that's more complex.

----
#### Variant Location ((varLoc))

Indicates which clinically-important functional domain the variant overlaps, if any.  If the variant does not overlap a clinically-important domain, then estimates of the impact of the variant on the protein function are not relevant.


----
#### General Fields

These fields are common to each of computational prediction methods, and are constant across all variants.

* #### Method ((Method))
    The name of the computational prediction method
* #### Description ((Description))
    A brief description of the method
* #### Publication ((Publication))
    PubMed ID and hyperlink to a publication in PubMed that describes the method

----
#### SpliceAI

SpliceAI is a deep neural network method that predicts the effect of genetic variants on pre-mRNA splicing.  There are four ways in which a variant could impact splicing:

1. By weakening or removing an existing splice site at the start of an intron (*Donor Loss*)

2. By weakening or removing an existing splice site at the start of an exon (*Acceptor Loss*)

3. By adding a new splice site that effectively starts a new intron prematurely as compared to the reference sequence  (*Donor Gain*)

4. By adding a new splice site that effectively terminates an intron prematurely as compared to the reference sequence (*Acceptor Gain*)

For any given variant, SpliceAI first estimates the probability that the *wild type allele* (the reference sequence without the variant) at the location of the variant would be part of a splice donor (intron start), a splice acceptor (intron end) or neither.  Then it considers the variant, and re-estimates the probability that the *alternative allele* (the reference sequence as modified by the variant) would be a splice donor, a splice acceptor, or neither.  In comparing these probabilities, it estimates how the variant would impact the probability of a loss of an existing splice site (*donor* or *acceptor loss*) or the creation of a new splice site (*donor* or *acceptor gain*).  It also indicates the location of the splice site loss or gain relative to the variant, where negative positions indicate a splice site upstream from the variant, positive positions indicate a splice site downstream from the variant, and a position of 0 indicates a splice site at the variant itself.

The largest of SpliceAI's four delta scores (representing the gain/loss of a donor/acceptor splice site) is widely used as an overall estimate of the probability that the variant will impact pre-mRNA splicing.

* #### Delta Score Donor Loss ((DS_DL_spliceAI))
    The Delta Score Donor Loss indicates the likelihood that the variant will impact splicing by weakening or removing an existing donor splice site, an exon/intron border.

* #### Delta Score Acceptor Loss ((DS_AL_spliceAI))
    The Delta Score Acceptor Loss indicates the likelihood that the variant will impact splicing by weakening or removing an existing acceptor splice site, an intron/exon border.

* #### Delta Score Donor Gain ((DS_DG_spliceAI))
    The Delta Score Donor Gain indicates the likelihood that the variant will impact splicing by introducing a new donor splice site, an exon/intron border.

* #### Delta Score Acceptor Gain ((DS_AG_spliceAI))
    The Delta Score Acceptor Gain indicates the likelihood that the variant will impact splicing by introducing a new acceptor splice site, an intron/exon border.

* #### Delta Position Donor Loss ((DP_DL_spliceAI))
    For variants predicted to weaken or remove an existing donor splice site, the delta position donor loss indicates the number of nucleotides in the pre-mRNA transcript between the variant and the donor splice site in question.  Negative values indicate transcript positions downstream of the variant, and positive values indicate transcript positions upstream of the variant.

* #### Delta Position Acceptor Loss ((DP_AL_spliceAI))
    For variants predicted to weaken or remove an existing acceptor splice site, the Delta Position Acceptor Loss indicates the number of nucleotides in the pre-mRNA transcript between the variant and the acceptor splice site in question.  Negative values indicate transcript positions downstream of the variant, and positive values indicate transcript positions upstream of the variant.


* #### Delta Position Donor Gain ((DP_DG_spliceAI))
    For variants predicted to create a new donor splice site, the Delta Position Donor Gain indicates the number of nucleotides in the pre-mRNA transcript between the variant and this new donor splice site.  Negative values indicate transcript positions downstream of the variant, and positive values indicate transcript positions upstream of the variant.


* #### Delta Position Acceptor Gain ((DP_AG_spliceAI))
    For variants predicted to create a new acceptor splice site, the Delta Position Acceptor Gain indicates the number of nucleotides in the pre-mRNA transcript between the variant and this new acceptor splice site.  Negative values indicate transcript positions downstream of the variant, and positive values indicate transcript positions upstream of the variant.

* #### Result ((result_spliceai))
    The largest of the four delta scores (Delta Score Donor Loss, Delta Score Acceptor Loss, Delta Score Donor Gain, Delta Score Acceptor Gain) represents an overall estimate of the overall likelihood that the variant will impact pre-mRNA splicing.


----
#### BayesDel

BayesDel is a machine learning method that predicts the impact of the variant on protein function or abundance.  It produces a numeric score, where a higher score represents a greater likelihood that the variant will have a significant functional impact.

In *BRCA1* and *BRCA2*, variants that fall outside clinically-important functional domains are rarely pathogenic, even if they are predicted to alter protein function.  This is because these protein regions are unstructured, and any changes to the protein sequence in these regions has little effect on protein function.  For this reason, BayesDel predictions are not relevant outside clinically-important domains.

* #### Result ((BayesDel_nsfp33a_noAF))
    The result of BayesDel is a numeric score, where negative scores suggest no functional impact, positive scores suggest possible functional impact, and larger positive scores predict functional impact wiht a greater likelihood.
