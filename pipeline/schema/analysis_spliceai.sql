-- SpliceAI scores, one row per variant.

CREATE TABLE IF NOT EXISTS pipeline.analysis_spliceai (
    "VRS_Digest" text PRIMARY KEY
        REFERENCES pipeline.variant("VRS_Digest") ON DELETE CASCADE,
    "DS_AG"  text,
    "DS_AL"  text,
    "DS_DG"  text,
    "DS_DL"  text,
    "DP_AG"  text,
    "DP_AL"  text,
    "DP_DG"  text,
    "DP_DL"  text,
    result   text
);
