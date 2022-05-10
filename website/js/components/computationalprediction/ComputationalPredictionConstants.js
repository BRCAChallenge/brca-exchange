/*global module: false, require: false, URL: false, Blob: false */
'use strict';

const ComputationalPredictionConstants = [{
        Method: "BayesDel",
        Description: "BayesDel is a meta-predictor of variant pathogenicity",
        Publication: "27995669",
    },
    {
        Method: "SpliceAI",
        Description: "SpliceAI estimates the likelihood that the variant impacts splicing in one of four possible ways.  The strongest of the four scores is generally used to indicate the probability of impact on splicing.",
        Publication: "30661751",
    }
];

module.exports = ComputationalPredictionConstants;
