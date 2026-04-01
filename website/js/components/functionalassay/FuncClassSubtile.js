'use strict';
import React from 'react';
import * as d3 from 'd3';

import varScoresArray from './mockdata/funcscores_array';
import {impacts} from "./FunctionalAssayTile";

export default class FuncClassSubtile extends React.Component {
    constructor(props) {
        super(props);
        this.createBarChart = this.createBarChart.bind(this);
    }

    componentDidMount() {
        this.createBarChart();
    }

    componentDidUpdate() {
        this.createBarChart();
    }

    createBarChart() {
        const {score} = this.props;
        // FIXME: 'varScoresArray' was generated from a copy of the findlay functional scores.
        //  ideally this should be generated/exported by the pipeline to the front-end, so it stays in sync with the source.
        const values = varScoresArray;

        const margin = { top: 0, bottom: 80, left: 65, right: 20 };
        const width = 500 - margin.left - margin.right;
        const height = 150 - margin.top - margin.bottom;

        const max = d3.max(values);
        const min = d3.min(values);
        const x = d3.scaleLinear()
            .domain([min, max])
            .range([0, width])
            .clamp(true);

        // Generate a histogram using twenty uniformly-spaced bins.
        const data = d3.histogram()
            .domain(x.domain())
	    .thresholds(x.ticks(40))
            (values);

        const yMax = d3.max(data, d => d.length);
        // const yMin = d3.min(data, d => d.length);

        const y = d3.scaleLinear()
            .domain([0, yMax])
            .range([height, 0]);

        const xAxis = d3.axisBottom(x)
	    .ticks(20);

        const yAxis = d3.axisLeft(y)
            .ticks(4)
            .tickSize(0);

        // gridlines along y (separate axis instance with extended tickSize)
        const yAxisGrid = d3.axisLeft(y)
            .ticks(4)
	    .tickSize(-width)
            .tickFormat(() => "");

        // const svgElem = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
        const svg = d3.select("#func-assay-obj")
            .attr("class", "func-assay")
            .append("g")
            .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

        // add horizontal grid lines
        svg.append("g")
            .attr("class", "grid")
            .call(yAxisGrid);

        // draw histogram bars
        const bar = svg.selectAll(".bar")
            .data(data)
            .enter().append("g")
            .attr("class", "bar")
            .attr("transform", d => "translate(" + x(d.x0) + "," + y(d.length) + ")");
        bar.append("rect")
            .attr("x", 1)
            .attr("width", d => (x(d.x1) - x(d.x0)) - 1)
            .attr("height", d => height - y(d.length))
            .attr("fill", "#7eb6ea");

        // draw the x-axis...
        svg.append("g")
            .attr("class", "x axis")
            .attr("transform", "translate(0," + height + ")")
            .call(xAxis)
            .append("text") // add x-axis "Function Score" label
            .attr("class", "axis-label")
            .attr("y", 35)
            .attr("x", 170)
            .style("text-anchor", "middle")
            .text("Function Score");

        // ...and the y-axis
        svg.append("g")
            .attr("class", "y axis")
            .call(yAxis)
            .append("text") // add y-axis "SNVs" label
            .attr("class", "axis-label")
            .attr("transform", "rotate(-90)")
            .attr("y", -30)
            .attr("x", -70)
            .style("text-anchor", "middle")
            .text("SNVs");

        // draw classification regions below chart
        const regions = svg.selectAll(".region")
            .data(impacts)
            .enter().append("g")
            .attr("transform", d => `translate(${x(d.range[0])},${height + 45})`);

        regions.append("rect")
            .attr("width", d => x(Math.min(d.range[1], max)) - x(Math.max(d.range[0], min)))
            .attr("height", 10)
            .attr("fill", d => d.color);

        regions.append("text")
            .attr("width", d => x(Math.min(d.range[1], max)) - x(Math.max(d.range[0], min)))
            .attr("x", (d, i) => (i === 0) ? 0 : x(Math.min(d.range[1], max)) - x(Math.max(d.range[0], min)))
            .attr("dy", 30)
            .attr("style", "font-size: 12px")
            .attr("text-anchor", (d, i) => (i === 0) ? "start" : "end")
            .text(d => d.label);

        // draw caret on classification chart (only if score is valid)
	const scoreX = x(score);
	if (scoreX !== undefined && scoreX !== null && !Number.isNaN(scoreX)) {
        svg.append("g")
            .attr("transform", `translate(${x(score)}, ${height + 50})`)
            .append("polygon")
            .attr("stroke", "none")
            .attr("fill", "black")
            .attr("points", "0,0 5,10 -5,10");

        svg.append("g")
            .attr("transform", `translate(${x(score)}, ${height + 50})`)
            .append("line")
            .attr("stroke", "black")
            .attr("stroke-dasharray", 9)
            .attr("stroke-width", 1)
            .attr("fill", "black")
            .attr("y", 0)
            .attr("y2", -yMax);
	}
    }

    render() {
        return (<svg id="func-assay-obj" width="100%" viewBox="0 0 500 150" />);
    }
}
