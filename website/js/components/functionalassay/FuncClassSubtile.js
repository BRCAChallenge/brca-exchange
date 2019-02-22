'use strict';
import React from 'react';
import d3 from 'd3';

import varScores from './mockdata/variants_to_funcscores';
import {impacts} from "./FunctionalAssayTile";

export default class FuncClassSubtile extends React.Component {
    constructor(props) {
        super(props);
        this.data = Object.values(varScores);
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
        var values = this.data;

        const margin = { top: 0, bottom: 80, left: 65, right: 20 };
        const width = 500 - margin.left - margin.right;
        const height = 150 - margin.top - margin.bottom;

        var max = d3.max(values);
        var min = d3.min(values);
        var x = d3.scale.linear()
            .domain([min, max])
            .range([0, width])
            .clamp(true);

        // Generate a histogram using twenty uniformly-spaced bins.
        var data = d3.layout.histogram()
            .bins(x.ticks(40))
            (values);

        var yMax = d3.max(data, d => d.length);
        // var yMin = d3.min(data, d => d.length);

        var y = d3.scale.linear()
            .domain([0, yMax])
            .range([height, 0]);

        var xAxis = d3.svg.axis()
            .scale(x)
            .ticks(20)
            .orient("bottom");

        var yAxis = d3.svg.axis()
            .scale(y)
            .ticks(4)
            .tickSize(0)
            .orient("left");

        // FIXME: instead of duplicating the axis object, figure out how to clone yAxis
        var yAxis2 = d3.svg.axis()
            .scale(y)
            .ticks(4)
            .orient("left");

        // var svgElem = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
        var svg = d3.select("#func-assay-obj")
            .attr("class", "func-assay")
            .append("g")
            .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

        // add horizontal grid lines
        svg.append("g")
            .attr("class", "grid")
            .call(yAxis2
                .tickSize(-width)
                .tickFormat("")
            );

        // draw histogram bars
        var bar = svg.selectAll(".bar")
            .data(data)
            .enter().append("g")
            .attr("class", "bar")
            .attr("transform", d => "translate(" + x(d.x) + "," + y(d.y) + ")");
        bar.append("rect")
            .attr("x", 1)
            .attr("width", (x(data[0].dx) - x(0)) - 1)
            .attr("height", d => height - y(d.y))
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
            .style("text-anchor", "center")
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
            .style("text-anchor", "center")
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

        // draw caret on classification chart
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
            .attr("stroke-width", 3)
            .attr("fill", "black")
            .attr("y", 0)
            .attr("y2", -yMax)
            .attr("points", "0,0 5,10 -5,10");

    }

    render() {
        return (<svg id="func-assay-obj" width="100%" height="auto" viewBox="0 0 500 150" />);
    }
}
