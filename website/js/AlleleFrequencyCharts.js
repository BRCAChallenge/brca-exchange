'use strict';

var React = require('react'),
    BarChart = require('./BarChart');

// copied from index.js, move to utility file
function isEmptyField(value) {
    if (Array.isArray(value)) {
        value = value[0];
    }

    if (value === null || (typeof value === 'undefined')) {
        return true;
    }

    var v = value.trim();
    return v === '' || v === '-' || v === 'None';
}

var alleleFrequencyCharts = function (variant, prop) {
    let frequencyProps, title, pointFormat;
    if (prop === 'Allele_Frequency_Charts_1000_Genomes') {
        if (!variant['Variant_in_1000_Genomes']) { // eslint-disable-line dot-notation
            return false;
        }
        frequencyProps = [
            {label: 'AFR', prop: 'AFR_Allele_frequency_1000_Genomes'},
            {label: 'AMR', prop: 'AMR_Allele_frequency_1000_Genomes'},
            {label: 'EAS', prop: 'EAS_Allele_frequency_1000_Genomes'},
            {label: 'EUR', prop: 'EUR_Allele_frequency_1000_Genomes'},
            {label: 'SAS', prop: 'SAS_Allele_frequency_1000_Genomes'},
        ];
        title = '1000 Genomes';
        pointFormat =  "{point.y}";
    }
    else if (prop === 'Allele_Frequency_Charts_ExAC') {
        if (!variant['Variant_in_ExAC']) {
            return false;
        }
        frequencyProps = [
            {label: 'AFR', prop: 'Allele_frequency_AFR_ExAC'},
            {label: 'AMR', prop: 'Allele_frequency_AMR_ExAC'},
            {label: 'EAS', prop: 'Allele_frequency_EAS_ExAC'},
            {label: 'FIN', prop: 'Allele_frequency_FIN_ExAC'},
            {label: 'NFE', prop: 'Allele_frequency_NFE_ExAC'},
            {label: 'OTH', prop: 'Allele_frequency_OTH_ExAC'},
            {label: 'SAS', prop: 'Allele_frequency_SAS_ExAC'},
        ];
        title = 'ExAC';
        pointFormat =  "{point.y}<br /><em>({point.count} of {point.number})</em>";
    } else {
        return false;
    }
    /*
    let frequencyProps = [
        //{label: 'Maximum Allele Frequency (1000 Genomes and ESP)', prop: 'Max_Allele_Frequency'},
        //{label: 'Allele Frequency (1000 Genomes)', prop: 'Allele_frequency_1000_Genomes'},
        // ESP
        {label: 'EA', prop: 'EA_Allele_Frequency_ESP'},
        {label: 'AA', prop: 'AA_Allele_Frequency_ESP'},
        {label: 'Allele Frequency (ESP)', prop: 'Allele_Frequency_ESP'},
        // ExAC
        {label: 'Allele frequency (ExAC minus TCGA)', prop: 'Allele_frequency_ExAC'},
    ];
    */

    let categories = [];
    let data = [];

    for (let frequencyType of frequencyProps) {
        if (!isEmptyField(variant[frequencyType.prop])) {
            categories.push(frequencyType.label);
            data.push({
                y: parseFloat(variant[frequencyType.prop]),
                number: variant[frequencyType.prop.replace("frequency", "number")],
                count: variant[frequencyType.prop.replace("frequency", "count")]
            });
        }
    }

    let chart2Max = Math.max(parseFloat(variant['Max_Allele_Frequency']), 0.01);
    let fullscaleChartOptions = {
        title: { text: title},
        legend: { enabled: false },
        plotOptions: { column: { tooltip: { pointFormat: pointFormat } } },
        chart: { spacing: [ 8, 4, 14, 8 ], height: 300 },
        yAxis: {
            max: 1.0,
            tickInterval: null,
            labels: {
                rotation: 45,
                style: { fontSize: "8px" },
                x: -2,
                y: 2
            },
            plotLines: [{
                color: 'red',
                dashStyle: 'dash',
                value: chart2Max,
                width: 1,
                zIndex: 5
            }]
        },
        xAxis: { categories: categories },
        series: [{ data: data }]
    };

    let scaledChartOptions = {
        title: { text: `${title} (scaled)`},
        legend: { enabled: false },
        plotOptions: { column: { tooltip: { pointFormat: pointFormat } } },
        chart: { spacing: [ 8, 10, 14, 0 ], height: 300 },
        yAxis: {
            max: chart2Max,
            tickInterval: null,
            labels: {
                rotation: 45,
                style: { fontSize: "8px" },
                x: -2,
                y: 2
            }
        },
        xAxis: { categories: categories },
        series: [{ data: data }]
    };
    return [
        <div className='alleleFrequencyChart'>
            <BarChart container={`${prop}_alleleFreq1`} options={fullscaleChartOptions} />
        </div>,
        <div className='alleleFrequencyChart'>
            <BarChart container={`${prop}_alleleFreq2`} options={scaledChartOptions} />
        </div>,
    ];
};

module.exports = alleleFrequencyCharts;
