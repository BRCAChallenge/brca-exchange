# Project: brca-exchange

## Stack
- Django 5.2 / Python 3.13 / PostgreSQL 18 (Ubuntu 24.04, PGDG apt repo)
- DRF + drf-spectacular for OpenAPI schema generation
- Node.js 22 for any frontend tooling
- Deployed and developed via VSCode Remote-SSH; DBeaver connects through an SSH tunnel for ad hoc schema exploration

## Domain context (assume this knowledge, don't re-explain it to me)
- Core entity is the `variant` table: genomic coordinates, pathogenicity calls, cross-references to ENIGMA, ClinVar, LOVD, exLOVD, gnomAD
- Supporting tables: reports, in silico predictions, protein structures, data releases
- Working genes of interest: BRCA1/BRCA2
- External APIs/resources in regular use: ClinGen Allele Registry (HGVS/VRS lookups), GA4GH VRS (with SeqRepo/UTA), gnomAD v4 (LCR intervals, exome capture regions, coverage fields like `total_dp`, indel size boundaries)
- ENIGMA VCEP evidence rules (BA1/BS1/PM2, gnomAD non-cancer FAF) are relevant when touching pathogenicity/classification logic

## Conventions
- Prefer Django ORM migrations over raw SQL; flag if a change needs a manual migration
- Match existing DRF serializer/viewset patterns rather than introducing new ones
- Any script hitting external APIs (ClinGen, gnomAD) should be idempotent/rate-limit aware — treat these as slow, rate-limited, sometimes flaky
- Luigi is the automation framework for PostgreSQL/ETL-style tasks — extend existing tasks rather than writing standalone scripts when the work is pipeline-shaped

## Commands
- Test: `<fill in your test command>`
- Lint/format: `<fill in>`
- Local server: `<fill in>`
- Migrations: `python manage.py makemigrations && python manage.py migrate`

## Things to always ask before doing
- Any destructive DB operation on `brca-prod` (this is treated as a production server)
- Schema changes touching the `variant` table or its cross-reference tables
- Anything that would re-download large reference datasets (gnomAD, ClinVar) — confirm scope first, these are large

## Style
- Keep explanations terse; I'm the domain expert on the genetics/bioinformatics side, you're the domain expert on Django/infra side — no need to explain BRCA biology or basic Django patterns to me
