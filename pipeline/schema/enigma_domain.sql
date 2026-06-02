-- ENIGMA Consortium functional domains of potential clinical importance.
-- Coordinates are GRCh38; start < end for both genes.
-- Source: website/js/SplicingData.js (brca1CIDomains / brca2CIDomains)

CREATE TABLE IF NOT EXISTS pipeline.enigma_domain (
    id       serial  PRIMARY KEY,
    gene     text    NOT NULL,
    name     text    NOT NULL,
    chrom    text    NOT NULL,
    assembly text    NOT NULL DEFAULT 'GRCh38',
    start    integer NOT NULL,
    "end"    integer NOT NULL
);

INSERT INTO pipeline.enigma_domain (gene, name, chrom, assembly, start, "end") VALUES
    ('BRCA1', 'ring',          '17', 'GRCh38', 43104260, 43124093),
    ('BRCA1', 'coiled-coil',   '17', 'GRCh38', 43082489, 43090958),
    ('BRCA1', 'brct',          '17', 'GRCh38', 43045699, 43070966),
    ('BRCA2', 'palb2 binding', '13', 'GRCh38', 32316488, 32319129),
    ('BRCA2', 'dna binding',   '13', 'GRCh38', 32356433, 32396954);
