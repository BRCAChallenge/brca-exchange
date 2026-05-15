-- BayesDel scores, one row per variant.

CREATE TABLE IF NOT EXISTS pipeline.analysis_bayesdel (
    "VRS_Digest"             text PRIMARY KEY
        REFERENCES pipeline.variant("VRS_Digest") ON DELETE CASCADE,
    "BayesDel_nsfp33a_noAF" text
);
