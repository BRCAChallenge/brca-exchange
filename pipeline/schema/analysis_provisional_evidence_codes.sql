-- Provisional evidence codes, one row per variant.

CREATE TABLE IF NOT EXISTS pipeline.analysis_provisional_evidence_codes (
    "VRS_Digest"        text PRIMARY KEY
        REFERENCES pipeline.variant("VRS_Digest") ON DELETE CASCADE,
    popfreq_code        text,
    popfreq_description text
);
