'use strict';

var React = require('react');
var ReactDOM = require('react-dom');
var _ = require('lodash');

import {CollapsableMixin} from 'react-bootstrap';
import classNames from "classnames";
import * as d3s from 'd3-scale';
import update from 'immutability-helper';

import {geneMeta} from './SplicingData.js';

// for some bizarre reason svg clipping masks don't pick this up from the CSS
const exonBorderRadius = 4;

const intronWidth = 40;
const leaderSize = 75, tailSize = 50; // made up values, how do I get UTR sizes for BRCA1 and BRCA2?
const zoomMargin = 200;
const intronMag = 2; // factor by which the intron for a fully-intronic variant is scaled

const segmentDims = {
    transcript: {
        intron: { height: 13, yOffset: 8 },
        exon: { height: 30, }
    },
    zoom: {
        intron: { height: 25, yOffset: 7 },
        exon: { height: 40 }
    },
};


// --------------------------------------------------------------------------------------------------------------
// --- supporting methods
// --------------------------------------------------------------------------------------------------------------

function exonSizeTx(bases) {
    return bases > 150 ? 150 + Math.sqrt(bases - 150) : bases;
}

function intronSizeTx() {
    // like exonSizeTx() this method receives bases as an argument,
    // in case we ever want to scale by the real size of the intron
    return intronWidth;
}

function variantInfo(variant) {
    // let [before, after] = _.map(variant.Genomic_Coordinate_hg38.split(':').pop().split('>'), e => e.length);
    let before = variant.Ref.length, after = variant.Alt.length;

    if (before === 1 && after === 1) {
        return { changed: 1, inserted: 0, deleted: 0 };
    }
    before--; after--; // since CG>G is a removal of one base, not a removal of two and insertion of one
    return {
        changed: Math.min(before, after),
        inserted: after > before ? after - before : 0,
        deleted: before > after ? before - after : 0
    };
}

function sortCoords([a, b]) {
    return [a < b ? a : b, b > a ? b : a];
}

// given the intervals (a: [a1, a2], b: [b1, b2]), return true if the two overlap (inclusive)
function overlaps(a, b) {
    // ensure the pairs are internally sorted
    const as = (a[0] <= a[1]) ? a : [a[1], a[0]];
    const bs = (b[0] <= b[1]) ? b : [b[1], b[0]];
    // sort each pair to reduce the number of cases we need to check
    const [first, second] = (as[0] < bs[0]) ? [as, bs] : [bs, as];

    // they can only overlap if first extends into (including surpassing) second
    return first[1] >= second[0];
}

// given an array, returns pairs of successive elements; e.g. [1,2,3] produces [[1,2],[2,3]]
// (used to create introns between pairs of exons later in the code)
function pairwise(seq) {
    return _.zip(_.take(seq, seq.length - 1), _.tail(seq));
}


// --------------------------------------------------------------------------------------------------------------
// --- components
// --------------------------------------------------------------------------------------------------------------

/**
 * A colored overlay on a specific exon or intron (aka a segment). Entities that span segment boundaries (e.g. variants,
 * acceptor/donor sites) are split into one region per overlapped segment.
 */
class Region extends React.Component {
    render() {
        let { region, width, height, txStart, txEnd, scale, mask, selected, nudgeable, className } = this.props;

        // txStart, txEnd are the parent exon/intron's span in bases
        // height is the pixel height of the parent element
        // scale is a d3 scaling object that converts from genomic coords to pixel positions in the parent

        // no point drawing an invisible zero-width event or one that doesn't overlap this container
        if (Math.abs(region.end - region.start) <= 0 || !overlaps([region.start, region.end], [txStart, txEnd])) {
            return null;
        }

        // project genomic coords into pixel positions
        const bpStartPx = scale(region.start);
        const bpEndPx = scale(region.end);

        // convert coordinates (minning b/c they may be swapped) into a pos and a width
        let bpMinPx = Math.min(bpStartPx, bpEndPx);
        let widthPx = Math.abs(bpStartPx - bpEndPx);

        // constrain to a minimum width, shifting away from the margin if necessary
        if (nudgeable && widthPx < 2) {
            widthPx = 2;

            if (bpMinPx + widthPx > width) {
                bpMinPx -= widthPx;
            }
        }

        return (
            <g>
                <rect x={bpMinPx} width={widthPx} height={height}
                    className={className}
                    clipPath={mask && `url(#${mask})`}
                />

                {/*
                draw the outline separately, ignoring the mask, so we can see it around variants
                */}
                { selected &&
                    <rect x={bpMinPx} width={widthPx} height={height}
                        rx={exonBorderRadius}
                        fill="transparent"
                        className={`selected-ci-path`}
                    />
                }
            </g>
        );
    }
}

/**
 * A set of Regions, which can potentially render multiple regions for a variant onto the given segment.
 *
 * It reads the variant info and possibly creates the following regions:
 * - insertion:
 *      if the alt sequence is longer than the ref, it's an insertion and renders an inserted area
 * - deletion:
 *      if the ref sequence is longer than alt, it's a deletion and renders a deleted area
 * - change:
 *      if alt and ref both have lengths > 1, there's an overlapping 'changed' region of width
 *      min(length(ref),length(alt)), rendered before the insertion/deletion region
 * - single-base change:
 *      if alt and ref both have lengths == 1, then only a changed region of width 1 is rendered
 */
class Variant extends React.Component {
    render() {
        let { variant, width, height, txStart, txEnd, scale, mask } = this.props;

        const variantStart = variant.Hg38_Start;
        const delta = variantInfo(variant);

        // the colors below are defined in the top-level svg's <defs> section
        const events = {
            changed: {
                widthBP: delta.changed,
                span: {start: variantStart, end: variantStart + delta.changed}
            },
            deleted: {
                widthBP: delta.deleted,
                span: {start: variantStart + delta.changed, end: variantStart + delta.changed + delta.deleted}
            },
            inserted: {
                widthBP: delta.inserted,
                // FIXME: should insertions be drawn as points, not intervals, since there's no corresponding region in the source to annotate?
                span: {start: variantStart + delta.changed, end: variantStart + delta.changed + delta.inserted}
            },
        };

        return (
            <g>
            {
                _.toPairs(events)
                    .filter((keyAndEvent) => keyAndEvent[1].widthBP > 0)
                    .map(([key, event]) =>
                        <Region key={`event_${key}`} region={event.span}
                            className={`variant-region ${key}`}
                            x={0} width={width} height={height} txStart={txStart} txEnd={txEnd} scale={scale}
                            mask={mask} nudgeable={true}
                        />
                    )
            }
            </g>
        );
    }
}

/**
 * Generates regions as overlays on a given segment.
 *
 * Current region types rendered:
 * - donor, acceptor sites
 * - variant inserted, deleted, changed regions
 */
class SegmentRegions extends React.Component {
    render() {
        const {
            mask, variant, donors, acceptors, CIDomains, width, height, txStart, txEnd, zoomed, isFlipped, isIntron
        } = this.props;
        const n = this.props.n || 'intron'; // section indicator, used for creating region keys
        const flatDomains = _.flatMap(CIDomains, (v, k) => v.domains.map((x) => ({
            org: k, name: x.name, code: v.code,
            span: {start: x.start, end: x.end}
        })));

        // regions drawn within this segment will use this scale to convert
        // from BP positions to (segment-relative) pixel positions
        const scale = d3s.scaleLinear()
            .clamp(true);

        // we'll append to these going forward
        const [txS, txE] = sortCoords([txStart, txEnd]);
        let domain = [txS];
        let range = [0];

        // exons and introns alternate in starting with a donor or acceptor,
        // and which site comes first depends on the strandedness (i.e. 'isFlipped')
        const sites = ((!isFlipped && isIntron) || (isFlipped && !isIntron)) ? [donors, acceptors] : [acceptors, donors];

        // after having ordered the sites, we can build the scale's intervals in a (somewhat) order-independent way
        sites.forEach((site, idx) => {
            // our start and end exons are missing a donor/acceptor, thus the following check
            if (site.length > 0) {
                const isFirst = idx === 0;
                const [siteStart, siteEnd] = sortCoords([site[0].start, site[0].end]);

                if (isFirst) {
                    domain.push(siteEnd);
                    range.push(isIntron ? width * (!zoomed ? 0.3 : 0.1) : 5);
                }
                else {
                    domain.push(siteStart);
                    range.push(isIntron ? width * (!zoomed ? 0.7 : 0.9) : width - 5);
                }
            }
        });

        // cap them off
        domain.push(txE);
        range.push(width);

        // apply the intervals we've constructed to the scale
        scale.domain(domain).range(range);

        // note that the order of elements defines their rendering order;
        // later elements will be rendered on top of earlier ones

        return (
            <g>
            {
                flatDomains
                    .filter(({span}) => overlaps([span.start, span.end], [txStart, txEnd]))
                    .map(({org, name, code, span}, idx) =>
                        <Region key={`cidomain_${org}_${name}_${idx}`}
                            className={`region cidomain domain-${code}`}
                            region={span} label={zoomed && `${org}: ${name}`}
                            x={0} width={width} height={height}
                            txStart={txStart} txEnd={txEnd} scale={scale}
                            mask={mask} // fill={CIDomainFills[org]}
                            selected={this.props.selectedDomain === `${org}_${name}`}
                        />
                    )
            }

            {
                donors.map((donorSpan, idx) => (
                    <Region key={`donor_${n}_${idx}`} region={donorSpan}
                        className="region donor"
                        x={0} width={width} height={height} txStart={txStart} txEnd={txEnd} scale={scale}
                        mask={mask}
                    />
                ))
            }

            {
                acceptors.map((acceptorSpan, idx) => (
                    <Region key={`acceptor_${n}_${idx}`} region={acceptorSpan}
                        className="region acceptor"
                        x={0} width={width} height={height} txStart={txStart} txEnd={txEnd} scale={scale}
                        mask={mask}
                    />
                ))
            }

            {
                /* draw variants over everything else */
                variant && (
                    <Variant variant={variant}
                        x={0} width={width} height={height} txStart={txStart} txEnd={txEnd}
                        scale={scale}
                        mask={mask}
                    />
                )
            }
            </g>
        );
    }
}

class Exon extends React.Component {
    render() {
        const {
            n, txStart, txEnd, width, height, x, variant, zoomed, highlight, isFlipped,
            donors, acceptors, CIDomains
        } = this.props;

        // the clip mask allows us to draw variants within the rounded-rectangle exon
        // we need to assign different mask IDs for each exon, in either zoomed-in or full-transcript mode
        const clipMaskID = `exon-${n}-clip-${zoomed ? 'zoomed' : 'full'}`;

        return (
            <g>
                <text x={x + width / 2} y={height + 14} textAnchor="middle">{n}</text>

                <g transform={isFlipped ? `translate(${x + width}) scale(-1,1)` : `translate(${x})`}>
                    <defs>
                        <clipPath id={clipMaskID}>
                            <rect
                                x={0} width={width} height={height} rx={exonBorderRadius} ry={exonBorderRadius}
                                className={`segment exon`}
                            />
                        </clipPath>
                    </defs>

                    <rect x={0} className={`segment exon ${highlight && 'highlighted'}`} width={width} height={height} rx={exonBorderRadius} ry={exonBorderRadius} />

                    <SegmentRegions n={n} variant={variant} donors={donors} acceptors={acceptors} CIDomains={CIDomains}
                        selectedDomain={this.props.selectedDomain}
                        txStart={txStart} txEnd={txEnd} zoomed={zoomed} isFlipped={isFlipped} isIntron={false}
                        width={width} height={height} mask={clipMaskID} />

                    {/* extra unfilled rect overlay, used to re-draw the outline on top of the regions */}
                    <rect x={0} className="segment exon outline" width={width} height={height} rx={exonBorderRadius} ry={exonBorderRadius} />
                </g>
            </g>
        );
    }
}

class Intron extends React.Component {
    render() {
        const {
            txStart, txEnd, x, height, width, variant, zoomed, highlight, isFlipped,
            donors, acceptors, CIDomains
        } = this.props;

        return (
            <g transform={`translate(0, ${zoomed ? segmentDims.zoom.intron.yOffset : segmentDims.transcript.intron.yOffset})`}>
                <g transform={isFlipped ? `translate(${x + width}) scale(-1,1)` : `translate(${x})`}>
                    <rect x={0} className={`segment intron ${highlight && 'highlighted'}`} width={width} height={height} />

                    <SegmentRegions variant={variant} donors={donors} acceptors={acceptors} CIDomains={CIDomains}
                        selectedDomain={this.props.selectedDomain}
                        txStart={txStart} txEnd={txEnd} zoomed={zoomed} isFlipped={isFlipped} isIntron={true}
                        width={width} height={height} />

                    {/* extra unfilled rect overlay, used to re-draw the outline on top of the regions */}
                    <rect x={0} className="segment intron outline" width={width} height={height} />
                </g>
            </g>
        );
    }
}

class Transcript extends React.Component {
    render() {
        const {variant, donors, acceptors,  CIDomains, segments, width, isFlipped} = this.props;


        // ------------------------------------------------
        // --- precalculate scale
        // --- TODO: scales are slightly different between exonic / intronic variants. fix.
        // ------------------------------------------------

        let totalWidth = leaderSize + tailSize + _.sum(segments.map((segment) => {
            return (segment.type === 'exon')
                ? exonSizeTx(segment.span.end - segment.span.start)
                : intronSizeTx(segment.span.end - segment.span.start);
        }));

        // set per-element scale from totalWidth
        const scale = (this.props.width - 2) / totalWidth;


        // ------------------------------------------------
        // --- create visible elements for each segment
        // ------------------------------------------------

        // we precalculate the element offsets so that we don't have to keep state in our segment mapping below
        // we also need to remember the position + width of the last element to place the 'tail' intron
        const agg = segments.reduce((agg, segment) => {
            // figure out the width of the current element
            const curWidth = (segment.type === 'exon')
                ? scale * exonSizeTx(segment.span.end - segment.span.start)
                : scale * intronSizeTx(segment.span.end - segment.span.start);

            // add the offset of this segment and move our cursor to the next position
            agg.offsets.push({x: agg.curX, width: curWidth, highlighted: segment.highlighted}); agg.curX += curWidth;

            return agg;
        }, {offsets: [], curX: leaderSize * scale});
        const tailPos = agg.curX;
        const offsets = agg.offsets;

        const blocks = segments.map((segment, i) => {
            const passedVariant = overlaps([variant.Hg38_Start, variant.Hg38_End], [segment.span.start, segment.span.end]) ? variant : null;

            const relevantDonors = Object.values(donors)
                    .filter(site => overlaps([site.start, site.end], [segment.span.start, segment.span.end]));
            const relevantAcceptors = Object.values(acceptors)
                    .filter(site => overlaps([site.start, site.end], [segment.span.start, segment.span.end]));

            if (segment.type === 'exon') {
                return (
                    <Exon key={`exon_${segment.id}`} n={segment.id}
                        variant={passedVariant} donors={relevantDonors} acceptors={relevantAcceptors} CIDomains={CIDomains}
                        selectedDomain={this.props.selectedDomain}
                        txStart={segment.span.start} txEnd={segment.span.end}
                        x={offsets[i].x} width={offsets[i].width} height={segmentDims.transcript.exon.height}
                        highlight={segment.highlighted} zoomed={false} isFlipped={isFlipped}
                    />
                );
            }
            else if (segment.type === 'intron') {
                return (
                    <Intron key={`intron_${segment.id}`} n={segment.id}
                        variant={passedVariant} donors={relevantDonors} acceptors={relevantAcceptors} CIDomains={CIDomains}
                        selectedDomain={this.props.selectedDomain}
                        txStart={segment.span.start} txEnd={segment.span.end}
                        x={offsets[i].x} width={offsets[i].width} height={segmentDims.transcript.intron.height}
                        highlight={segment.highlighted} zoomed={false} isFlipped={isFlipped}
                    />
                );
            }
        });

        // zoom lines should always be determined by the first and last highlighted element
        const highlightedOffsets = offsets.filter(x => x.highlighted);
        const firstHighlight = _.first(highlightedOffsets), lastHighlight = _.last(highlightedOffsets);

        const zoomLineStartLeft = firstHighlight.x;
        const zoomLineStartRight = lastHighlight.x + lastHighlight.width;

        return (
            <g>
                <line x1={zoomLineStartLeft} y1={50} x2={zoomMargin} y2={65} className="zoomline" />
                <line x1={zoomLineStartRight} y1={50} x2={width - zoomMargin - 2.5} y2={65} className="zoomline" />

                <rect x={0} y={53} width={820} height={14} fill="rgba(255,255,255,0.6)" />

                <g transform="translate(0, 10)">
                    <rect x={1} y={segmentDims.transcript.intron.yOffset}
                        width={scale * leaderSize} height={segmentDims.transcript.intron.height}
                        className="segment intron cap" />
                    { blocks }
                    <rect x={tailPos} y={segmentDims.transcript.intron.yOffset}
                        width={scale * tailSize} height={segmentDims.transcript.intron.height}
                        className="segment intron cap" />
                </g>
            </g>
        );
    }
}

class Zoom extends React.Component {
    render() {
        const { variant, donors, acceptors, CIDomains, segments, width, isFullyIntronic, isFlipped } = this.props;


        // ------------------------------------------------
        // --- precalculate scale
        // ------------------------------------------------

        let totalWidth = _.sum(segments.map((segment) => {
            return (segment.type === 'exon')
                ? exonSizeTx(segment.span.end - segment.span.start)
                : intronSizeTx(segment.span.end - segment.span.start) * (isFullyIntronic ? intronMag : 1);
        }));

        // set per-element scale from totalWidth
        const scale = (width - 2 * zoomMargin) / totalWidth;


        // ------------------------------------------------
        // --- create visible elements for each segment
        // ------------------------------------------------

        // we precalculate the element offsets so that we don't have to keep state in our segment mapping below
        const offsets = segments.reduce((agg, segment) => {
            // figure out the width of the current element
            const curWidth = (segment.type === 'exon')
                ? scale * exonSizeTx(segment.span.end - segment.span.start)
                : scale * intronSizeTx(segment.span.end - segment.span.start) * (isFullyIntronic ? intronMag : 1 );

            // add the offset of this segment and move our cursor to the next position
            agg.offsets.push({x: agg.curX, width: curWidth}); agg.curX += curWidth;

            return agg;
        }, {offsets: [], curX: zoomMargin}).offsets;

        const blocks = segments.map((segment, i) => {
            const passedVariant = overlaps([variant.Hg38_Start, variant.Hg38_End], [segment.span.start, segment.span.end]) ? variant : null;

            const relevantDonors = Object.values(donors)
                    .filter(site => overlaps([site.start, site.end], [segment.span.start, segment.span.end]));
            const relevantAcceptors = Object.values(acceptors)
                    .filter(site => overlaps([site.start, site.end], [segment.span.start, segment.span.end]));

            if (segment.type === 'exon') {
                return (
                    <Exon key={`exon_${segment.id}`} n={segment.id}
                        variant={passedVariant} donors={relevantDonors} acceptors={relevantAcceptors} CIDomains={CIDomains}
                        txStart={segment.span.start} txEnd={segment.span.end}
                        selectedDomain={this.props.selectedDomain}
                        x={offsets[i].x} width={offsets[i].width} height={segmentDims.zoom.exon.height}
                        highlight={true} zoomed={true} isFlipped={isFlipped}
                    />
                );
            }
            else if (segment.type === 'intron') {
                return (
                    <Intron key={`intron_${i + 1}`} n={segment.id}
                        variant={passedVariant} donors={relevantDonors} acceptors={relevantAcceptors} CIDomains={CIDomains}
                        txStart={segment.span.start} txEnd={segment.span.end}
                        selectedDomain={this.props.selectedDomain}
                        x={offsets[i].x} width={offsets[i].width} height={segmentDims.zoom.intron.height}
                        highlight={true} zoomed={true} isFlipped={isFlipped}
                    />
                );
            }
        });

        return (
            <g transform="translate(0, 70)">
            { blocks }
            </g>
        );
    }
}

const SettingsPanel = React.createClass({
    mixins: [CollapsableMixin],

    getCollapsableDOMNode: function() {
        return ReactDOM.findDOMNode(this.refs.panel);
    },

    getCollapsableDimensionValue: function() {
        return ReactDOM.findDOMNode(this.refs.panel).scrollHeight;
    },

    onHandleToggle: function (e) {
        e.preventDefault();

        // ask our parent to toggle us
        this.props.onToggled();
    },

    render: function() {
        const styles = this.getCollapsableClassSet();

        const settings = (
            <div className="container-fluid">
                <div className="row">
                    <div className="col-md-6">
                        { this.props.generateCIDomainSelectors(this.props.meta) }
                    </div>

                    <div className="col-md-6">
                        <div>
                            <label>
                                <input style={{marginRight: '0.5em'}} type="checkbox" name="drawDonors" checked={this.props.drawDonors} onChange={this.props.toggleDrawing} />
                                <svg className="site-indicator" width={18} height={18}>
                                    <rect className="donor" width={18} height={18} />
                                </svg>
                                Donor Sites
                            </label>
                        </div>

                        <div>
                            <label style={{display: 'inline-block', marginRight: '1em', verticalAlign: 'baseline'}}>
                                <input style={{marginRight: '0.5em'}} type="checkbox" name="drawAcceptors" checked={this.props.drawAcceptors} onChange={this.props.toggleDrawing} />
                                <svg className="site-indicator" width={18} height={18}>
                                    <rect className="acceptor" width={18} height={18} />
                                </svg>
                                Acceptor Sites
                            </label>
                        </div>

                        <div>
                            <label>
                                <input style={{marginRight: '0.5em'}} type="checkbox" name="alternatePalette" checked={this.props.alternatePalette} onChange={this.props.toggleDrawing} />
                                Alternate Palette
                            </label>
                        </div>
                    </div>
                </div>
            </div>
        );

        if (this.props.nonCollapsable) {
            return (
                <div style={{marginBottom: 0, borderTop: 'solid 1px #ccc', paddingTop: '10px'}}>{settings}</div>
            );
        }

        return (
            <div>
                <div style={{marginBottom: 0, borderTop: 'solid 1px #ccc', padding: '10px'}}>
                    <div className={`submitter-header ${this.props.expanded ? 'expanded' : ''}`}>
                        <span onClick={this.onHandleToggle} style={{fontWeight: 'bold', cursor: 'pointer'}}>
                            {
                                this.props.expanded
                                    ? <i className="fa fa-caret-down" aria-hidden="true" />
                                    : <i className="fa fa-caret-right" aria-hidden="true" />
                            }
                            <span>&nbsp;Options</span>
                        </span>
                    </div>
                </div>

                <div ref='panel' className={classNames(styles)}>
                {settings}
                </div>
            </div>
        );
    }
});


class Splicing extends React.Component {
    constructor(props) {
        super(props);

        this.state = {
            drawAcceptors: true,
            drawDonors: true,
            drawCIDomains: new Set(),
            selectedDomain: null,
            alternatePalette: false,
            optionsExpanded: localStorage.getItem("transcriptviz-options-expanded") === 'true'
        };

        this.onOptionsToggled = this.onOptionsToggled.bind(this);

        this.toggleDrawing = this.toggleDrawing.bind(this);
        this.toggleCIDomain = this.toggleCIDomain.bind(this);
        this.selectCIDomain = this.selectCIDomain.bind(this);
        this.generateCIDomainSelectors = this.generateCIDomainSelectors.bind(this);
    }

    toggleDrawing(event) {
        this.setState({
            [event.target.name]: event.target.checked
        });
    }

    onOptionsToggled() {
        this.setState((pstate) => ({
            optionsExpanded: !pstate.optionsExpanded
        }), () => {
            // persist expansion state to localstorage
            localStorage.setItem("transcriptviz-options-expanded", this.state.optionsExpanded ? "true" : "false");

            // reflows the parent after our dimensions change
            this.props.onContentsChanged(this.collapser.getCollapsableDOMNode());
        });
    }

    toggleCIDomain(event) {
        const boxChecked = event.target.checked;
        const CIDomain = event.target.name;

        // toggles a set of CI domains
        // clears the selected domain if it's within this domain and the set is being disabled
        this.setState((pstate) => {
            return {
                drawCIDomains: update(pstate.drawCIDomains, boxChecked ? {$add: [CIDomain]} :  {$remove: [CIDomain]}),
                selectedDomain: (!boxChecked && pstate.selectedDomain && pstate.selectedDomain.startsWith(CIDomain))
                    ? null
                    : pstate.selectedDomain
            };
        });
    }

    selectCIDomain(domain, org) {
        this.setState((pstate) => {
            // toggle displaying any domain if we click on the same one again
            return {
                selectedDomain: pstate.selectedDomain !== domain ? domain : null,
                drawCIDomains: update(pstate.drawCIDomains, {$add: [org]})
            };
        });
    }

    render() {
        let { variant } = this.props;

        const width = 800, info = variantInfo(variant);
        const variantStart = variant.Pos | 0;
        const variantEnd = variantStart + info.changed + info.deleted + info.inserted;

        // --- pre-step: get data, sort and format it so we can process it
        const meta = geneMeta[variant['Gene_Symbol']];
        const exons = _.toPairs(meta.exons).map(([name, span]) => ({id: parseInt(name.substr(4)), span}));
        const variantSpan = [variantStart, variantEnd];
        const isFlipped = meta.strand === '-';

        // get the segments (exons and introns interleaved),
        // with segments that overlap the variant's span marked as highlighted
        const segments = this.buildSegments(exons, variantSpan, isFlipped);

        // if segments is null, it indicates that the variant doesn't overlap this gene at all
        if (!segments) {
            console.log(`Variant ${variantSpan} falls outside of gene`);
            return (
                <h4 style={{textAlign: 'center'}}>Variant is outside of transcript.</h4>
            );
        }

        // if it's fully intronic, we blow up the size of the intron in the zoom mode
        const isFullyIntronic = segments.length === 3 && segments[1].segment.type === 'intron';

        let plural = n => n === 1 ? '' : 's';

        // filter the CIDomains according to which ones we're including
        const filteredCIDomains = _.pickBy(meta.CIDomains, (v, k) => this.state.drawCIDomains.has(k));

        return (
            <div className={`transcript-viz ${this.state.alternatePalette ? 'altpalette' : ''}`}>
                <svg viewBox="-4 0 808 160" preserveAspectRatio="xMidYMid">
                    {/* variant fill definitions, declared here instead of in CSS b/c one of them is a <pattern> */}
                    <defs>
                        <pattern id="diagonalHatch" patternUnits="userSpaceOnUse" width="4" height="4">
                            <rect x="0" y="0" width="4" height="4" className="deleted-fill-rect" />
                            <path d="M-1,1 l2,-2
                           M0,4 l4,-4
                           M3,5 l2,-2"
                                className="deleted-diagonal-lines" />
                        </pattern>
                    </defs>

                    {/* transcript and zoomed-in parts */}
                    <Transcript variant={variant} segments={segments} width={width}
                        donors={this.state.drawDonors ? meta.spliceDonors : {}}
                        acceptors={this.state.drawAcceptors ? meta.spliceAcceptors : {}}
                        CIDomains={filteredCIDomains} selectedDomain={this.state.selectedDomain}
                        isFullyIntronic={false} isFlipped={isFlipped}
                    />

                    <Zoom variant={variant} segments={segments.filter(x => x.highlighted)} width={width}
                        donors={this.state.drawDonors ? meta.spliceDonors : {}}
                        acceptors={this.state.drawAcceptors ? meta.spliceAcceptors : {}}
                        CIDomains={filteredCIDomains} selectedDomain={this.state.selectedDomain}
                        isFullyIntronic={isFullyIntronic} isFlipped={isFlipped}
                    />

                    {/* legend */}
                    <g transform="translate(274,140)">
                        <rect x="0" className="changed" stroke="black" width="20" height="10" />
                        <text x="22" y="10">{ `Substitution (${info.changed} base${plural(info.changed)})` }</text>
                        <rect x="192" className="deleted" stroke="black" width="20" height="10" />
                        <text x="214" y="10">{ `Deletion (${info.deleted || 0} base${plural(info.deleted)})` }</text>
                        <rect x="360" className="inserted" stroke="black" width="20" height="10"/>
                        <text x="382" y="10">{ `Insertion (${info.inserted || 0} base${plural(info.inserted)})` }</text>
                    </g>
                </svg>

                <SettingsPanel
                    ref={(me) => { this.collapser = me; }}
                    key="settings-panel"
                    meta={meta}
                    expanded={this.state.optionsExpanded}
                    drawDonors={this.state.drawDonors}
                    drawAcceptors={this.state.drawAcceptors}
                    alternatePalette={this.state.alternatePalette}
                    generateCIDomainSelectors={this.generateCIDomainSelectors}
                    toggleDrawing={this.toggleDrawing}
                    onToggled={this.onOptionsToggled}
                    nonCollapsable={false}
                />
            </div>
        );
    }

    /**
     * Given a set of exons, creates a set of segments (exons and interleaved introns).
     *
     * Each segment has a highlighted field, indicating if it lies within the region of interest, i.e.
     * the segments that overlap the variant plus some buffer exons/introns depending on whether the variant is intronic.
     * @param exons the set of exons for this gene, given as [{id, span: {start, end}}, ...]
     * @param variantSpan start and end of this variant in genomic coords, given as [start, end]
     * @param isFlipped whether this gene is sense ('+') or antisense ('-')
     * @returns {Array} an array of segments, i.e. [{type: 'intron'/'exon', id, span: {start,end}, highlighted}, ...]
     */
    buildSegments(exons, variantSpan, isFlipped) {
        // disregard strandedness so we can just build some intervals
        const sortedExons = _.sortBy(
            exons.map(exon => ({
                id: exon.id,
                span: {
                    // the +1 makes the starting bound exclusive
                    start: Math.min(exon.span.start, exon.span.end) + 1,
                    end: Math.max(exon.span.start, exon.span.end),
                }
            })),
            (a) => a.span.start
        );

        // sanity check: verify that the variant actually overlaps the gene at all
        // FIXME: should this reject variants that don't lie entirely within this gene? currently we look for just an overlap
        const geneSpan = [sortedExons[0].span.start, sortedExons[sortedExons.length - 1].span.end];
        if (!overlaps(variantSpan, geneSpan)) {
            return null;
        }

        // --- step 1: build list of interleaved exons and introns, aka segments

        // a 'segment' is either an exon or intron
        // here we build a set of these segments, i.e. a series of exons with following introns
        // (we need to use the sorted set so that 'start' and 'end' are consistent regardless of the gene's direction)
        let segments = _.flatten(
            // insert introns between each pair of exons
            pairwise(sortedExons).map(([prevExon, nextExon]) => {
                return [
                    {type: 'exon', ...prevExon},
                    {
                        type: 'intron',
                        id: (Math.abs(prevExon.id + nextExon.id) / 2.0), // kind of a hack, but it allows us to sort easily
                        span: {start: prevExon.span.end + 1, end: nextExon.span.start - 1}
                    }
                ];
            })
        );
        // and stick on the last element
        segments.push({type: 'exon', ...(_.last(sortedExons))});


        // --- step 2: identify the region of interest

        // identify the highlighted boundary
        // (variants are defined inclusively in the starting coordinate and exclusively in the ending one)
        const overlappingSegments = segments
            .map((segment, idx) => ({idx: idx, segment: segment}))
            .filter(({segment}) => overlaps([variantSpan[0], variantSpan[1] - 1], [segment.span.start, segment.span.end]));

        const firstSeg = _.first(overlappingSegments), lastSeg = _.last(overlappingSegments);

        // compute region that we'll be zooming in on in the figure
        // the region of interest is one exon before and one exon after the segments that overlap the variant
        // (if we're in an exon we need to move over by 2 elements to get to an exon, but in an intron it's only 1)
        const precedingExonIndex = Math.max(
            (firstSeg.segment.type === 'intron') ? firstSeg.idx - 1 : firstSeg.idx - 2, 0
        );
        const followingExonIndex = Math.min(
            (lastSeg.segment.type === 'intron') ? lastSeg.idx + 1 : lastSeg.idx + 2, segments.length - 1
        );

        // update the objects to highlight them if they're in the region of interest
        _.forEach(segments, (segment, idx) => {
            segment.highlighted = (idx >= precedingExonIndex && idx <= followingExonIndex);
        });

        // if it's antisense, we got the exons in descending ID
        // we need to flip the collection so that it's ordered from exon 1 to exon N
        if (isFlipped) {
            segments.reverse();
        }

        return segments;
    }

    /**
     * Creates UI elements to control visibility of clinically important (CI) domains.
     *
     * @param meta gene metadata, including the CI domains and their locations
     */
    generateCIDomainSelectors(meta) {
        return _.toPairs(meta.CIDomains).map(([org, orgMeta]) =>
            <div key={org}>
                <label style={{display: 'inline-block', marginRight: '1em'}}>
                    <input style={{marginRight: '0.5em'}} type="checkbox"
                        name={org} checked={this.state.drawCIDomains.has(org)} onChange={this.toggleCIDomain}
                    />
                    <svg className="site-indicator" width={18} height={18}>
                        <rect width={18} height={18} className={`domain-${orgMeta.code}`} />
                    </svg>
                    {orgMeta.label}
                </label>

                <ol className="splicing-domain-list">
                    {
                        orgMeta.domains.map((domain, idx) => {
                            const selected = (this.state.selectedDomain === `${org}_${domain.name}`);

                            return (
                                <li key={idx}>
                                    <a style={{fontWeight: selected ? 'bold' : 'normal'}}
                                        onClick={() => this.selectCIDomain(`${org}_${domain.name}`, org)}>
                                        {domain.name}
                                    </a>
                                </li>
                            );
                        })
                    }
                </ol>
            </div>
        );
    }
}

module.exports = Splicing;
