-- VEP annotation results, one row per variant.

CREATE TABLE IF NOT EXISTS pipeline.analysis_vep (
    "VRS_Digest"  text PRIMARY KEY
        REFERENCES pipeline.variant("VRS_Digest") ON DELETE CASCADE,
    variant_class text
);
