'use strict';

/*
THIS FILE CONTAINS DATA DETAILING BRCA1/BRCA2 EXONS,
SPLICE SITES, AND CLINICALLY IMPORTANT REGIONS
*/

/*
NOTE: BRCA1 is transcribed on the minus strand, meaning it is read from right
to left instead of the usual left to right. Genomic coordinates are always
based on the plus strand so the start coordinate for BRCA1 occurs to the right
of the end coordinate in genomic position (so it has a greater value)
*/

/*
Exon boundary data
*/

const brca1Exons = {
    'exon1': {
        'end': 43125270,
        'start': 43125483
    },
    'exon2': {
        'end': 43124016,
        'start': 43124115
    },
    'exon3': {
        'end': 43115725,
        'start': 43115779
    },
    'exon5': {
        'end': 43106455,
        'start': 43106533
    },
    'exon6': {
        'end': 43104867,
        'start': 43104956
    },
    'exon7': {
        'end': 43104121,
        'start': 43104261
    },
    'exon8': {
        'end': 43099774,
        'start': 43099880
    },
    'exon9': {
        'end': 43097243,
        'start': 43097289
    },
    'exon10': {
        'end': 43095845,
        'start': 43095922
    },
    'exon11': {
        'end': 43091434,
        'start': 43094860
    },
    'exon12': {
        'end': 43090943,
        'start': 43091032
    },
    'exon13': {
        'end': 43082403,
        'start': 43082575
    },
    'exon14': {
        'end': 43076487,
        'start': 43076614
    },
    'exon15': {
        'end': 43074330,
        'start': 43074521
    },
    'exon16': {
        'end': 43070927,
        'start': 43071238
    },
    'exon17': {
        'end': 43067607,
        'start': 43067695
    },
    'exon18': {
        'end': 43063873,
        'start': 43063951
    },
    'exon19': {
        'end': 43063332,
        'start': 43063373
    },
    'exon20': {
        'end': 43057051,
        'start': 43057135
    },
    'exon21': {
        'end': 43051062,
        'start': 43051117
    },
    'exon22': {
        'end': 43049120,
        'start': 43049194
    },
    'exon23': {
        'end': 43047642,
        'start': 43047703
    },
    'exon24': {
        'end': 43044294,
        'start': 43045802
    }
};

const brca2Exons = {
    'exon1': {
        'end': 32315667,
        'start': 32315479
    },
    'exon2': {
        'end': 32316527,
        'start': 32316421
    },
    'exon3': {
        'end': 32319325,
        'start': 32319076
    },
    'exon4': {
        'end': 32325184,
        'start': 32325075
    },
    'exon5': {
        'end': 32326150,
        'start': 32326100
    },
    'exon6': {
        'end': 32326282,
        'start': 32326241
    },
    'exon7': {
        'end': 32326613,
        'start': 32326498
    },
    'exon8': {
        'end': 32329492,
        'start': 32329442
    },
    'exon9': {
        'end': 32331030,
        'start': 32330918
    },
    'exon10': {
        'end': 32333387,
        'start': 32332271
    },
    'exon11': {
        'end': 32341196,
        'start': 32336264
    },
    'exon12': {
        'end': 32344653,
        'start': 32344557
    },
    'exon13': {
        'end': 32346896,
        'start': 32346826
    },
    'exon14': {
        'end': 32355288,
        'start': 32354860
    },
    'exon15': {
        'end': 32356609,
        'start': 32356427
    },
    'exon16': {
        'end': 32357929,
        'start': 32357741
    },
    'exon17': {
        'end': 32362693,
        'start': 32362522
    },
    'exon18': {
        'end': 32363533,
        'start': 32363178
    },
    'exon19': {
        'end': 32370557,
        'start': 32370401
    },
    'exon20': {
        'end': 32371100,
        'start': 32370955
    },
    'exon21': {
        'end': 32376791,
        'start': 32376669
    },
    'exon22': {
        'end': 32379515,
        'start': 32379316
    },
    'exon23': {
        'end': 32379913,
        'start': 32379749
    },
    'exon24': {
        'end': 32380145,
        'start': 32380006
    },
    'exon25': {
        'end': 32394933,
        'start': 32394688
    },
    'exon26': {
        'end': 32397044,
        'start': 32396897
    },
    'exon27': {
        'end': 32399672,
        'start': 32398161
    }
};

/*
Reference splice donor boundaries
*/
// See note at top of file for explanation of start > end.
const brca1RefSpliceDonorBounds = {
    'exon1': {
        'start': 43125273,
        'end': 43125265
    },
    'exon2': {
        'start': 43124019,
        'end': 43124011
    },
    'exon3': {
        'start': 43115728,
        'end': 43115720
    },
    'exon5': {
        'start': 43106458,
        'end': 43106450
    },
    'exon6': {
        'start': 43104870,
        'end': 43104862
    },
    'exon7': {
        'start': 43104124,
        'end': 43104116
    },
    'exon8': {
        'start': 43099777,
        'end': 43099769
    },
    'exon9': {
        'start': 43097246,
        'end': 43097238
    },
    'exon10': {
        'start': 43095848,
        'end': 43095840
    },
    'exon11': {
        'start': 43091437,
        'end': 43091429
    },
    'exon12': {
        'start': 43090946,
        'end': 43090938
    },
    'exon13': {
        'start': 43082406,
        'end': 43082398
    },
    'exon14': {
        'start': 43076490,
        'end': 43076482
    },
    'exon15': {
        'start': 43074333,
        'end': 43074325
    },
    'exon16': {
        'start': 43070930,
        'end': 43070922
    },
    'exon17': {
        'start': 43067610,
        'end': 43067602
    },
    'exon18': {
        'start': 43063876,
        'end': 43063868
    },
    'exon19': {
        'start': 43063335,
        'end': 43063327
    },
    'exon20': {
        'start': 43057054,
        'end': 43057046
    },
    'exon21': {
        'start': 43051065,
        'end': 43051057
    },
    'exon22': {
        'start': 43049123,
        'end': 43049115
    },
    'exon23': {
        'start': 43047645,
        'end': 43047637
    }
};

const brca2RefSpliceDonorBounds = {
    'exon1': {
        'start': 32315665,
        'end': 32315673
    },
    'exon2': {
        'start': 32316525,
        'end': 32316533
    },
    'exon3': {
        'start': 32319323,
        'end': 32319331
    },
    'exon4': {
        'start': 32325182,
        'end': 32325190
    },
    'exon5': {
        'start': 32326148,
        'end': 32326156
    },
    'exon6': {
        'start': 32326280,
        'end': 32326288
    },
    'exon7': {
        'start': 32326611,
        'end': 32326619
    },
    'exon8': {
        'start': 32329490,
        'end': 32329498
    },
    'exon9': {
        'start': 32331028,
        'end': 32331036
    },
    'exon10': {
        'start': 32333385,
        'end': 32333393
    },
    'exon11': {
        'start': 32341194,
        'end': 32341202
    },
    'exon12': {
        'start': 32344651,
        'end': 32344659
    },
    'exon13': {
        'start': 32346894,
        'end': 32346902
    },
    'exon14': {
        'start': 32355286,
        'end': 32355294
    },
    'exon15': {
        'start': 32356607,
        'end': 32356615
    },
    'exon16': {
        'start': 32357927,
        'end': 32357935
    },
    'exon17': {
        'start': 32362691,
        'end': 32362699
    },
    'exon18': {
        'start': 32363531,
        'end': 32363539
    },
    'exon19': {
        'start': 32370555,
        'end': 32370563
    },
    'exon20': {
        'start': 32371098,
        'end': 32371106
    },
    'exon21': {
        'start': 32376789,
        'end': 32376797
    },
    'exon22': {
        'start': 32379513,
        'end': 32379521
    },
    'exon23': {
        'start': 32379911,
        'end': 32379919
    },
    'exon24': {
        'start': 32380143,
        'end': 32380151
    },
    'exon25': {
        'start': 32394931,
        'end': 32394939
    },
    'exon26': {
        'start': 32397042,
        'end': 32397050
    }
};
/*
Reference splice acceptor boundaries
*/
// See note at top of file for explanation of start > end.
const brca1RefSpliceAcceptorBounds = {
    'exon2': {
        'start': 43124135,
        'end': 43124113
    },
    'exon3': {
        'start': 43115799,
        'end': 43115777
    },
    'exon5': {
        'start': 43106553,
        'end': 43106531
    },
    'exon6': {
        'start': 43104976,
        'end': 43104954
    },
    'exon7': {
        'start': 43104281,
        'end': 43104259
    },
    'exon8': {
        'start': 43099900,
        'end': 43099878
    },
    'exon9': {
        'start': 43097309,
        'end': 43097287
    },
    'exon10': {
        'start': 43095942,
        'end': 43095920
    },
    'exon11': {
        'start': 43094880,
        'end': 43094858
    },
    'exon12': {
        'start': 43091052,
        'end': 43091030
    },
    'exon13': {
        'start': 43082595,
        'end': 43082573
    },
    'exon14': {
        'start': 43076634,
        'end': 43076612
    },
    'exon15': {
        'start': 43074541,
        'end': 43074519
    },
    'exon16': {
        'start': 43071258,
        'end': 43071236
    },
    'exon17': {
        'start': 43067715,
        'end': 43067693
    },
    'exon18': {
        'start': 43063971,
        'end': 43063949
    },
    'exon19': {
        'start': 43063393,
        'end': 43063371
    },
    'exon20': {
        'start': 43057155,
        'end': 43057133
    },
    'exon21': {
        'start': 43051137,
        'end': 43051115
    },
    'exon22': {
        'start': 43049214,
        'end': 43049192
    },
    'exon23': {
        'start': 43047723,
        'end': 43047701
    },
    'exon24': {
        'start': 43045822,
        'end': 43045800
    }
};

const brca2RefSpliceAcceptorBounds = {
    'exon2': {
        'start': 32316402,
        'end': 32316424
    },
    'exon3': {
        'start': 32319057,
        'end': 32319079
    },
    'exon4': {
        'start': 32325056,
        'end': 32325078
    },
    'exon5': {
        'start': 32326081,
        'end': 32326103
    },
    'exon6': {
        'start': 32326222,
        'end': 32326244
    },
    'exon7': {
        'start': 32326479,
        'end': 32326501
    },
    'exon8': {
        'start': 32329423,
        'end': 32329445
    },
    'exon9': {
        'start': 32330899,
        'end': 32330921
    },
    'exon10': {
        'start': 32332252,
        'end': 32332274
    },
    'exon11': {
        'start': 32336245,
        'end': 32336267
    },
    'exon12': {
        'start': 32344538,
        'end': 32344560
    },
    'exon13': {
        'start': 32346807,
        'end': 32346829
    },
    'exon14': {
        'start': 32354841,
        'end': 32354863
    },
    'exon15': {
        'start': 32356408,
        'end': 32356430
    },
    'exon16': {
        'start': 32357722,
        'end': 32357744
    },
    'exon17': {
        'start': 32362503,
        'end': 32362525
    },
    'exon18': {
        'start': 32363159,
        'end': 32363181
    },
    'exon19': {
        'start': 32370382,
        'end': 32370404
    },
    'exon20': {
        'start': 32370936,
        'end': 32370958
    },
    'exon21': {
        'start': 32376650,
        'end': 32376672
    },
    'exon22': {
        'start': 32379297,
        'end': 32379319
    },
    'exon23': {
        'start': 32379730,
        'end': 32379752
    },
    'exon24': {
        'start': 32379987,
        'end': 32380009
    },
    'exon25': {
        'start': 32394669,
        'end': 32394691
    },
    'exon26': {
        'start': 32396878,
        'end': 32396900
    },
    'exon27': {
        'start': 32398142,
        'end': 32398164
    }
};

/*
Clinically important regions
*/
// See note at top of file for explanation of start > end.
const brca1CIDomains = {
    "ENIGMA Consortium": { // 'enigma'
        "code": 'enigma',
        "label": "Clinically Important Functional Domains (ENIGMA Consortium)",
        "domains": [
            {
                "name": "ring",
                "start": 43124093,
                "end": 43104260
            },
            {   "name": "coiled-coil",
                "start": 43090958,
                "end": 43082489
            },
            {
                "name": "brct",
                "start": 43070966,
                "end": 43045699
            }
        ]
    },
    "Huntsman Cancer Institute": { // 'priors'
        "code": 'huntsman',
        "label": "Huntsman Cancer Institute Functional Domains",
        "domains": [
            {
                "name": "initiation",
                "start": 43124096,
                "end": 43124094
            },
            /*
            {
                "name": "ring",
                "start": 43124084,
                "end": 43104875
            },
            {
                "name": "brct",
                "start": 43070966,
                "end": 43045705
            }
            */
        ]
    }
};

const brca2CIDomains = {
    "ENIGMA Consortium": { // 'enigma'
        "code": 'enigma',
        "label": "Clinically Important Functional Domains (ENIGMA Consortium)",
        "domains": [
            {
                "name": "palb2 binding",
                "start": 32316488,
                "end": 32319129
            },
            {
                "name": "dna binding",
                "start": 32356433,
                "end": 32396954
            }
        ]
    },
    "Huntsman Cancer Institute": { // 'priors'
        "code": 'huntsman',
        "label": "Huntsman Cancer Institute Functional Domains",
        "domains": [
            {
                "name": "initiation",
                "start": 32316461,
                "end": 32316463
            },
            {
                "name": "palb2",
                "start": 32316491,
                "end": 32319108
            },
            /*
            {
                "name": "dnb",
                "start": 32356433,
                "end": 32396954
            },
            */
            {
                "name": "tr2/rad5",
                "start": 32398318,
                "end": 32398428
            }
        ]
    }
};

// contains all the metadata for each gene, so we can load it easily at once
// sense field values: if '-' (aka antisense), start > end; if '+', start < end
const geneMeta = {
    'BRCA1': {
        strand: '-',
        exons: brca1Exons,
        spliceDonors: brca1RefSpliceDonorBounds,
        spliceAcceptors: brca1RefSpliceAcceptorBounds,
        CIDomains: brca1CIDomains
    },
    'BRCA2': {
        strand: '+',
        exons: brca2Exons,
        spliceDonors: brca2RefSpliceDonorBounds,
        spliceAcceptors: brca2RefSpliceAcceptorBounds,
        CIDomains: brca2CIDomains
    }
};

module.exports = {
    brca1Exons,
    brca2Exons,
    brca1RefSpliceDonorBounds,
    brca2RefSpliceDonorBounds,
    brca1RefSpliceAcceptorBounds,
    brca2RefSpliceAcceptorBounds,
    brca1CIDomains,
    brca2CIDomains,
    geneMeta
};
