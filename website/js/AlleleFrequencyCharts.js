'use strict';

var React = require('react'),
    BarChart = require('./BarChart'),
    util = require('./util'),
    _ = require('lodash');

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

    if (_.every(_.values(_.pick(variant, _.map(frequencyProps, e => e.prop))), util.isEmptyField)) {
        return false;
    }

    let categories = [];
    let data = [];

    for (let frequencyType of frequencyProps) {
        if (!util.isEmptyField(variant[frequencyType.prop])) {
            categories.push(frequencyType.label);
            data.push({
                y: parseFloat(variant[frequencyType.prop]),
                number: variant[frequencyType.prop.replace("frequency", "number")],
                count: variant[frequencyType.prop.replace("frequency", "count")]
            });
        }
    }

    let chart2Max = Math.max.apply(null, _.filter(_.map(_.values(_.pick(variant, _.map(frequencyProps, e => e.prop))), parseFloat), e => isFinite(e)));
    let lineColor = '#D00';
    let fullscaleChartOptions = {
        title: { text: title},
        legend: { enabled: false },
        plotOptions: {
            column: {
                tooltip: { pointFormat: pointFormat },
                dataLabels: { enabled: false }
            }
        },
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
                color: lineColor,
                dashStyle: 'Dash',
                value: 0.01,
                width: 1,
                zIndex: 5,
                label: {
                    text: '1%',
                    align: 'left',
                    style: { color: lineColor }
                }
            }]
        },
        xAxis: { categories: categories },
        series: [{ data: data }]
    };

    let scaledChartOptions = {
        title: { text: `${title} (scaled)`},
        legend: { enabled: false },
        plotOptions: {
            column: {
                tooltip: { pointFormat: pointFormat },
                dataLabels: { enabled: false }
            }
        },
        chart: { spacing: [ 8, 10, 14, 0 ], height: 300 },
        yAxis: {
            max: chart2Max,
            tickInterval: null,
            labels: {
                rotation: 45,
                style: { fontSize: "8px" },
                x: -2,
                y: 2
            },
            plotLines: [{
                color: lineColor,
                dashStyle: 'Dash',
                value: 0.01,
                width: 1,
                zIndex: 5,
                label: {
                    text: '1%',
                    align: 'left',
                    style: { color: lineColor }
                }
            }]
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
