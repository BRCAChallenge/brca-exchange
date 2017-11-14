'use strict';

module.exports = {
    sources: [
        {
            name: 'ClinVar',
            submissions: [
                {
                    submitter: 'GeneDx',
                    cols: [
                        { prop: 'clinical-significance', title: 'Clinical Significance', value: 'Benign' },
                        { prop: 'analysis-method', title: 'Analysis Method', value: 'clinical testing' },
                        { prop: 'date-last-updated', title: 'Date Last Updated', value: 'Feb 21, 2014' },
                        { prop: 'scv-accession', title: 'SCV Accession', value: 'SCV000167236.10' },
                        {
                            prop: 'supporting-observations',
                            title: 'Supporting Observations',
                            value: 'This variant is considered likely benign or benign based on one or more of the following criteria: it is a conservative change, it occurs at a poorly conserved position in the protein, it is predicted to be benign by multiple in silico algorithms, and/or has population frequency not consistent with disease.'
                        }
                    ]
                },
                {
                    submitter: 'Invitae',
                    cols: [
                        { prop: 'clinical-significance', title: 'Clinical Significance', value: 'Benign' },
                        { prop: 'analysis-method', title: 'Analysis Method', value: 'clinical testing' }
                    ]
                }
            ]
        },
        {
            name: 'LOVD',
            submissions: [
                {
                    submitter: 'Hans Gille (Amsterdam, NL)',
                    cols: [
                        { prop: 'variant-effect', title: 'Variant Effect', value: '?/?' },
                        { prop: 'genetic-origin', title: 'Genetic Origin', value: 'Germline' },
                        { prop: 'individuals', title: 'Individuals', value: '12' },
                        { prop: 'submission-id', title: 'Submission ID', value: '00011551' },
                        { prop: 'variant-haplotype', title: 'Variant Haplotype', value: '-' }
                    ]
                }
            ]
        }
    ]
};
