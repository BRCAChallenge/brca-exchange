'use strict';

var React = require('react'),
    _ = require('lodash'),
    brca1Exons = JSON.parse(require('raw!../content/brca1LollipopDomain.json')),
    brca2Exons = JSON.parse(require('raw!../content/brca2LollipopDomain.json'));

var intronWidth = 40,
    maxUnscaledExonWidth = 150, // base count
    exonFill = "#92c3e5",
    exonStroke = "rgba(0, 0, 0, 0.5)",
    intronFill = "#b5c4d2",
    intronStroke = "rgba(0, 0, 0, 0.5)",
    highlightFill = "#ffffbb",
    highlightStroke = "rgba(0, 0, 0, 0.5)",
    // made up values, how do I get UTR sizes for BRCA1 and BRCA2?
    leaderSize = 75,
    tailSize = 50,
    zoomMargin = 200;

function exonSizeTx (bases) {
    return bases > 150 ? 150 + Math.sqrt(bases - 150) : bases;
}

class Exon extends React.Component {
    constructor(props) {
        super(props);
    }
    componendDidMount() {
    }
    render() {
        let { txStart,
              txEnd,
              xScale,
              width,
              height,
              x,
              variant,
              highlight
        } = this.props;

        let variantStart, variantEnd, variantX, variantWidth;

        if (variant) {
            variantStart = variant["Hg38_Start"];
            variantEnd = variant["Hg38_End"];
            variantX = x + width * (variantStart - txStart) / (txEnd - txStart)
            variantWidth = width * (variantEnd - variantStart) / (txEnd - txStart);
            console.log(variantWidth);
            variantWidth = Math.max(variantWidth, 2);
        }

        return (
            <g>
                <rect x={x} width={width} height={height} rx={4} ry={4} fill={highlight ? highlightFill : exonFill} stroke={highlight ? highlightStroke : exonStroke} />
                { variant !== undefined &&
                    <rect x={variantX} width={variantWidth} height={height} fill="red" /> }
            </g>
        );
    }
}

class Intron extends React.Component {
    render() {
        let { txStart,
              txEnd,
              x,
              height,
              width,
              variant,
              highlight
        } = this.props;

        let variantStart, variantEnd, variantX, variantWidth;

        if (variant) {
            variantStart = variant["Hg38_Start"];
            variantEnd = variant["Hg38_End"];
            variantX = x + width * (variantStart - txStart) / (txEnd - txStart)
            variantWidth = width * (variantEnd - variantStart) / (txEnd - txStart);
            console.log(variantWidth);
            variantWidth = Math.max(variantWidth, 2);
        }

        //let width = variant !== undefined ? 2 * intronWidth : intronWidth;
        return (
            <g>
                <rect x={x} y={10} width={width} height={height} fill={highlight ? highlightFill : intronFill} stroke={highlight ? highlightStroke : intronStroke} />
                { variant !== undefined &&
                    <rect x={variantX} y={10} width={variantWidth} height={height} fill="red" /> }
            </g>
        );
    }
}

class Transcript extends React.Component {
    render() {
        let { variant,
              exons,
              preceding,
              following,
              width
        } = this.props;

        let scale,
            blocks = [],
            x,
            zoomLineStartLeft,
            zoomLineStartRight;

        // precalculate scale
        {
            let totalWidth = 0;
            for (let i = 0; i < exons.length; i++) {
                totalWidth += exonSizeTx(exons[i][1] - exons[i][0]);
                if (i < exons.length - 1) {
                    totalWidth += intronWidth;
                }
            }
            totalWidth += leaderSize + tailSize;
            // variant is intronic
            if (following - preceding === 1) {
                totalWidth += intronWidth;
            }
            scale = (this.props.width - 2) / totalWidth;
        }

        x = leaderSize * scale;

        for (let i = 0; i < exons.length; i++) {
            let exon = exons[i];
            let exonWidth = scale * exonSizeTx(exon[1] - exon[0]);

            // draw exon
            {
                let highlight, _variant;
                if (i >= preceding && i <= following) { // exon is adjacent to, or contains, variant
                    highlight = true;
                    if (variant["Hg38_Start"] >= exon[0] && variant["Hg38_Start"] <= exon[1]) { // exon contains variant
                        _variant = variant;
                    }
                }
                blocks.push(<Exon x={x} width={exonWidth} height={40} txStart={exon[0]} txEnd={exon[1]} highlight={highlight} variant={_variant} />);
            }

            if (i == preceding) {
                zoomLineStartLeft = x;
            } else if (i == following) {
                zoomLineStartRight = x + exonWidth;
            }

            x += exonWidth;

            // draw intron
            if (i < exons.length - 1) {
                let highlight, _variant;
                if (i >= preceding && i < following) {  // intron is adjacent to, or contains, variant
                    highlight = true;
                    if (variant["Hg38_Start"] > exon[1] && variant["Hg38_Start"] < exons[i+1][0]) { // intron contains variant
                        _variant = variant;
                    }
                }
                blocks.push(<Intron txStart={exon[1]} txEnd={exons[i+1][0]} x={x} width={scale * intronWidth} height={20} highlight={highlight} variant={_variant} />);
                x += scale * intronWidth;
            }
        }
        return (
            <g>
                <g transform="translate(0, 10)">
                    <rect x={1} y={8} width={scale * leaderSize} height={24} fill={intronFill} stroke={intronStroke} />
                    { blocks }
                    <rect x={x} y={8} width={scale * tailSize} height={24} fill={intronFill} stroke={intronStroke} />
                </g>
                <line x1={zoomLineStartLeft} y1={50} x2={zoomMargin} y2={95} stroke="#ccc" strokeWidth="3" strokeDasharray="8,3" fill="none" />
                <line x1={zoomLineStartRight} y1={50} x2={width - zoomMargin - 2.5} y2={95} stroke="#ccc" strokeWidth="3" strokeDasharray="8,3" fill="none" />
            </g>
        );
    }
}

class Zoom extends React.Component {
    render() {
        let { variant,
              exons,
              preceding,
              following,
              width
        } = this.props; 

        let blocks = [],
            scale,
            x = zoomMargin;


        // precalculate scale
        {  
            let totalWidth = 0;
            for (let i = preceding; i <= following; i++) {
                totalWidth += exonSizeTx(exons[i][1] - exons[i][0]);
                if (i < exons.length - 1) {
                    totalWidth += intronWidth;
                }
            }
            // variant is intronic
            if (following - preceding === 1) {
                totalWidth += intronWidth;
            }
            scale = (18 + width - 2 * zoomMargin) / totalWidth;
        }

        for (let i = preceding; i <= following; i++) {
            let exon = exons[i],
                exonWidth = scale * exonSizeTx(exon[1] - exon[0]);

            {
                let _variant;
                if (variant["Hg38_Start"] > exon[0] && variant["Hg38_Start"] < exon[1]) { // exon contains variant
                    _variant = variant;
                }
                blocks.push(<Exon txStart={exon[0]} txEnd={exon[1]} x={x} width={exonWidth} height={80} highlight={true} variant={_variant} />);
            }

            x += exonWidth;
            if (i != following) {
                let _variant;
                if (variant["Hg38_Start"] > exon[1] && variant["Hg38_Start"] < exons[i+1][0]) { // intron contains variant
                    _variant = variant;
                }
                blocks.push(<Intron txStart={exon[1]} txEnd={exons[i+1][0]} x={x} width={intronWidth} height={60} highlight={true} variant={_variant} />);
                x += intronWidth;
            }
        }
        return (
            <g transform="translate(0, 100)">
                { blocks }
            </g>
        );
    }
}

class Splicing extends React.Component {
    /* props:
       variant, width, height
    */
    componentDidMount() {
    }
    componentWillUnmount() {
    }
    render() {
        let { variant,
              width,
              height
        } = this.props;
        let variantStart = variant["Hg38_Start"],
            exons = _.map(variant["Gene_Symbol"] === "BRCA1" ? brca1Exons : brca2Exons, e => e.coord.split('-')),
            precedingExonIndex, followingExonIndex;

        for (let i = 0; i < exons.length; i++) {
            if (exons[i][0] > variantStart) {
                followingExonIndex = i;
                if (exons[i - 1][1] < variantStart) {
                    precedingExonIndex = i - 1;
                } else {
                    precedingExonIndex = i - 2;
                }
                break;
            }
        }
        return (
            <svg width={width} height={height}>
                <Transcript variant={variant} exons={exons} preceding={precedingExonIndex} following={followingExonIndex} width={width} />
                <Zoom variant={variant} exons={exons} preceding={precedingExonIndex} following={followingExonIndex} width={width} />
            </svg>
        );
    }
/*
        let variant = this.props.variant,
            exons = _.map(variant["Gene_Symbol"] === "BRCA1" ? brca1Exons : brca2Exons, e => e.coord.split('-')),
            exonBlocks = [],
            exonLabels = [],
            layout = [],
            lastExonEnd = 0,
            unscaledSize = leaderSize,
            variantLocation = variant["Hg38_Start"],
            precedingExonIndex, followingExonIndex;

        // determine exon/intron that contains variant
        for (let i = 0; i < exons.length; i++) {
            // ??? do we need to anticipate variants in the leader/tail?
            if (exons[i][0] > variantLocation) {
                followingExonIndex = i;
                if (exons[i - 1][1] < variantLocation) {
                    // variant is intronic
                    precedingExonIndex = i - 1;
                } else {
                    // variant is exonic
                    precedingExonIndex = i - 2;
                }
                break;
            }
        }

        // compute total width first, to allow for clean scaling
        for (let i = 0; i < exons.length; i++) {
            let exonWidth = exons[i][1] - exons[i][0];

            if (exonWidth > maxUnscaledExonWidth) {
                exonWidth = maxUnscaledExonWidth + Math.sqrt(exonWidth);
            }
            unscaledSize += intronWidth + exonWidth;

            // intronic variant
            if (precedingExonIndex === i && followingExonIndex === i + 1) {
                unscaledSize += intronWidth;
            }
        }
        unscaledSize += tailSize;

        let scaleFactor = 799 / unscaledSize;
        lastExonEnd = scaleFactor * leaderSize;
        let variantBlock;
        for (let i = 0; i < exons.length; i++) {
            let exonWidth = exons[i][1] - exons[i][0];

            if (exonWidth > maxUnscaledExonWidth) {
                exonWidth = maxUnscaledExonWidth + Math.sqrt(exonWidth);
            }

            exonLabels.push(<text x={lastExonEnd + (scaleFactor * exonWidth)/2} textAnchor='middle'>{i+1}</text>);

            if (variantLocation > exons[i][0] && variantLocation < exons[i][1]) {
                let start = scaleFactor * exonWidth * (variantLocation - exons[i][0])/(exons[i][1] - exons[i][0]);
                let width = scaleFactor * exonWidth * (variant["Hg38_End"] - variantLocation)/(exons[i][1] - exons[i][0]);
                variantBlock = <rect x={lastExonEnd + start} width={Math.max(width, 1)} height={40} fill="red" />;
            } else if (precedingExonIndex === i && followingExonIndex === i + 1) {
                let start = scaleFactor * intronWidth * 2 * (variantLocation - exons[i][1])/(exons[i + 1][0] - exons[i][1]);
                let width = scaleFactor * intronWidth * 2 * (variant["Hg38_End"] - variantLocation)/(exons[i + 1][0] - exons[i][1]);
                variantBlock = <rect x={lastExonEnd + scaleFactor * exonWidth + start} width={Math.max(width, 1)} y={15} height={10} fill="red" />;
            }

            lastExonEnd += scaleFactor * (intronWidth + exonWidth);
            // intronic variant
            if (precedingExonIndex === i && followingExonIndex === i + 1) {
                lastExonEnd += scaleFactor * intronWidth;
            }
        }
        return (
            <svg width={800} height={300}>
                <rect x={exonBlocks[0].props.x + 1} width={exonBlocks[exonBlocks.length - 1].props.x + 1} fill={intronFillColor} stroke={intronStrokeColor} y={35} height={10} />
                <g fill={exonFillColor} stroke={exonStrokeColor} strokeWidth={0.5} transform='translate(0,20)'>
                    <rect x={1} width={scaleFactor * leaderSize} y={10} height={20}  fill={intronFillColor}/>
                    { exonBlocks }
                    <rect x={exonBlocks[exonBlocks.length - 1].props.x + exonBlocks[exonBlocks.length - 1].props.width} width={scaleFactor * tailSize} y={10} height={20} fill={intronFillColor} />
                </g>
                <g transform='translate(0,20)'>
                    { variantBlock }
                </g>
                <g transform='translate(0,74)'>
                    { exonLabels }
                </g>
            </svg>
        );
    }
    */
}

module.exports = Splicing;
