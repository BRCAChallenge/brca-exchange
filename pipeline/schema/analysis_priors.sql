-- In silico prior probability fields shown in the UI, one row per variant.

CREATE TABLE IF NOT EXISTS pipeline.analysis_priors (
    "VRS_Digest"                    text PRIMARY KEY
        REFERENCES pipeline.variant("VRS_Digest") ON DELETE CASCADE,

    -- location / summary
    "varLoc"                        text,

    -- protein-level
    "proteinPrior"                  text,

    -- donor splice site
    "refDonorPrior"                 text,
    "refRefDonorMES"                text,
    "refRefDonorZ"                  text,
    "altRefDonorMES"                text,
    "altRefDonorZ"                  text,
    "refRefDonorSeq"                text,
    "altRefDonorSeq"                text,

    -- de novo donor
    "deNovoDonorPrior"              text,
    "refDeNovoDonorMES"             text,
    "refDeNovoDonorZ"               text,
    "altDeNovoDonorMES"             text,
    "altDeNovoDonorZ"               text,
    "refDeNovoDonorSeq"             text,
    "altDeNovoDonorSeq"             text,
    "deNovoDonorGenomicSplicePos"   text,
    "deNovoDonorTranscriptSplicePos" text,
    "closestDonorGenomicSplicePos"  text,
    "closestDonorTranscriptSplicePos" text,
    "closestDonorRefMES"            text,
    "closestDonorRefZ"              text,
    "closestDonorRefSeq"            text,
    "closestDonorAltMES"            text,
    "closestDonorAltZ"              text,
    "closestDonorAltSeq"            text,

    -- acceptor splice site
    "refAccPrior"                   text,
    "refRefAccMES"                  text,
    "refRefAccZ"                    text,
    "altRefAccMES"                  text,
    "altRefAccZ"                    text,
    "refRefAccSeq"                  text,
    "altRefAccSeq"                  text
);
