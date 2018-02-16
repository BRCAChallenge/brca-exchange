'use strict';

var React = require('react'),
    _ = require('lodash');

/*
var brca1Exons = JSON.parse(require('raw!../content/brca1LollipopDomain.json'));
var brca2Exons = JSON.parse(require('raw!../content/brca2LollipopDomain.json'));
*/

import {geneMeta} from './SplicingData.js';

const intronWidth = 40,
    exonFill = "#c1ddf0",
    exonStroke = "rgba(0, 0, 0, 0.5)",
    exonBorderRadius = 4,
    intronFill = "#dfe6ec",
    intronStroke = "rgba(0, 0, 0, 0.5)",
    highlightFill = "#ffffbb",
    highlightStroke = "rgba(0, 0, 0, 0.5)",
    // made up values, how do I get UTR sizes for BRCA1 and BRCA2?
    leaderSize = 75,
    tailSize = 50,
    zoomMargin = 20,
    intronMag = 2; // factor by which the intron for a fully-intronic variant is scaled

const donorFill = '#ce67b3',
    acceptorFill = '#7ba8ff';

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
    let [before, after] = _.map(variant.Genomic_Coordinate_hg38.split(':').pop().split('>'), e => e.length);
    if (before === 1 && after === 1) {
        return { changed: 1 };
    }
    before--; after--; // since CG>G is a removal of one base, not a removal of two and insertion of one
    return {
        changed: Math.min(before, after),
        inserted: after > before ? after - before : 0,
        deleted: before > after ? before - after : 0
    };
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

// ensure that the span lies entirely within the parent
function constrain(regionStart, regionWidth, start, width, minWidth) {
    // ensure the starting pos is no less than the parent's pos
    start = Math.max(regionStart, start);

    // if minWidth is given, ensure our pixel width is never less than that (caveats below)
    // * if this is done after the overshoot check, it could cause this element to overflow into the next region
    // * if this is done before the overshoot check, like it is now, then we can't guarantee a minimum width
    // * the final option is to adjust the start position to keep its minimum width, which is a bit of a white lie...
    if (minWidth) {
        width = Math.max(width, minWidth);
    }

    // ensure our width doesn't exceed the parent's end
    if (start + width > regionStart + regionWidth) {
        const overshootWidth = (start + width) - (regionStart + regionWidth);

        // console.log("would overshoot by ", overshootWidth, "; shortened to fit within ", regionWidth);

        // subtract the difference between our larger end and the parent's end,
        // which should put us at just touching the parent
        // FIXME: is this right? how should we deal with an event that overlaps the end of the region?
        start -= (overshootWidth);
    }

    // FIXME: in large regions, short events near the start or end can get obscured by the border
    // (find an example of this, shouldn't be too hard...)
    return { start, width };
}

// given an array, returns pairs of successive elements; e.g. [1,2,3] produces [[1,2],[2,3]]
function pairwise(seq) {
    return _.zip(_.take(seq, seq.length - 1), _.tail(seq));
}


// --------------------------------------------------------------------------------------------------------------
// --- components
// --------------------------------------------------------------------------------------------------------------

class Variant extends React.Component {
    render() {
        let { variant, x, width, height, txStart, txEnd } = this.props;

        // txStart, txEnd are the parent exon/intron's span in bases
        // width, height is the pixel width, height of the parent element
        // x is the pixel position of the parent exon/intron in the SVG

        // length of the parent exon/intro in bases
        const regionWidth = Math.abs(txEnd - txStart);

        let variantStart, variantX, variantChange, variantChangedWidth, variantDeletedWidth, variantInsertedWidth;

        variantStart = variant.Hg38_Start;
        variantX = x + (width * ((variantStart - txStart) / regionWidth));

        // variantStart is the start of this variant in bases
        // variantX is the pixel position of this variant (assumedly) within the parent exon/intron
        // (if variantX is negative, it's because the variant's start is before this region begins)

        variantChange = variantInfo(variant);

        // compute pixel widths for each event
        variantChangedWidth = variantChange.changed ? width * variantChange.changed / regionWidth : 0;
        variantDeletedWidth = variantChange.deleted ? width * variantChange.deleted / regionWidth : 0;
        variantInsertedWidth = variantChange.inserted ? width * variantChange.inserted / regionWidth : 0;

        // it's quite possible that the event begins before this region and/or ends after it.
        // in those cases it's safe to just constrain the rendered event to this elements' boundaries.
        // the other regions with partial overlaps will take care of rendering the event's full extent

        // describe each event, and constrain each one's rendered representation to lie within the parent
        const events = {
            changed: {
                fill: 'lightgreen',
                span: (variantChange.changed > 0 ? constrain(x, width, variantX, variantChangedWidth, 2) : null)
            },
            deleted: {
                fill: 'url(#diagonalHatch)',
                span: (variantChange.deleted > 0 ? constrain(x, width, variantX + variantChangedWidth, variantDeletedWidth, 2) : null)
            },
            inserted: {
                fill: 'lightblue',
                // FIXME: should insertions be drawn as points, not intervals, since there's no corresponding region in the source to annotate?
                span: (variantChange.inserted > 0 ? constrain(x, width, variantX + variantChangedWidth, variantInsertedWidth, 2) : null)
            },
        };

        return (
            <g>
            {
                // map each event with a valid span to a rect
                _.toPairs(events)
                    .filter((keyEvent) => keyEvent[1].span !== null)
                    .map(([key, event]) =>
                        (<rect key={`event_${key}`}
                            x={event.span.start} width={event.span.width} height={height}
                            fill={event.fill} clipPath={this.props.mask && `url(#${this.props.mask})`}
                        />)
                    )
            }
            </g>
        );
    }
}


class Region extends React.Component {
    render() {
        let { region, x, width, height, txStart, txEnd, fill, opacity, mask } = this.props;

        // txStart, txEnd are the parent exon/intron's span in bases
        // width, height is the pixel width, height of the parent element
        // x is the pixel position of the parent exon/intron in the SVG

        // length of the event in bases
        const eventWidth = Math.abs(region.end - region.start);
        const eventMin = Math.min(region.start, region.end), txMin = Math.min(txStart, txEnd);

        // no point drawing an invisible zero-width event
        // TODO: perhaps we can support these to be drawn as lines? maybe a different component is more appropriate
        if (eventWidth <= 0) {
            return null;
        }

        // length of the parent exon/intro in bases
        const regionWidth = Math.abs(txEnd - txStart);

        // pixel position of the event within this region (or possibly outside)
        // (if eventX is negative, it's because the variant's start is before this region begins)
        const eventX = x + (width * ((eventMin - txMin) / regionWidth));

        // compute pixel widths for event and constrain to lie within the parent
        const eventWidthPx = (width * eventWidth) / regionWidth;
        const span = constrain(x, width, eventX, eventWidthPx, 4);

        return (
            <rect x={span.start} width={span.width} height={height} fill={fill} opacity={opacity} clipPath={mask && `url(#${mask})`} />
        );
    }
}


class Exon extends React.Component {
    render() {
        let {
            n, txStart, txEnd, width, height, x, variant, zoomed, highlight, isFlipped,
            donors, acceptors
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
                            <rect x={0} width={width} height={height} rx={exonBorderRadius} ry={exonBorderRadius} />
                        </clipPath>
                    </defs>

                    <rect x={0} width={width} height={height} rx={exonBorderRadius} ry={exonBorderRadius} fill={highlight ? highlightFill : exonFill} stroke={highlight ? highlightStroke : exonStroke} />

                    {
                        variant &&
                        <Variant variant={variant} x={0} width={width} height={height} txStart={txStart} txEnd={txEnd} mask={clipMaskID} />
                    }

                    {
                        donors.map((donorSpan, idx) => (
                            <Region key={`donor_${n}_${idx}`} region={donorSpan}
                                x={0} width={width} height={height} txStart={txStart} txEnd={txEnd}
                                fill={donorFill} opacity={0.5} mask={clipMaskID}
                            />
                        ))
                    }

                    {
                        acceptors.map((acceptorSpan, idx) => (
                            <Region key={`acceptor_${n}_${idx}`} region={acceptorSpan}
                                x={0} width={width} height={height} txStart={txStart} txEnd={txEnd}
                                fill={acceptorFill} opacity={0.5} mask={clipMaskID}
                            />
                        ))
                    }
                </g>
            </g>
        );
    }
}

class Intron extends React.Component {
    render() {
        let {
            txStart, txEnd, x, height, width, variant, highlight, isFlipped,
            donors, acceptors
        } = this.props;

        return (
            <g transform="translate(0, 10)">
                <g transform={isFlipped ? `translate(${x + width}) scale(-1,1)` : `translate(${x})`}>
                    <rect x={0} width={width} height={height} fill={highlight ? highlightFill : intronFill} stroke={highlight ? highlightStroke : intronStroke} />
                    { variant && <Variant variant={variant} x={0} width={width} height={height} txStart={txStart} txEnd={txEnd} /> }

                    {
                        donors.map((donorSpan, idx) => (
                            <Region key={`donor_intron_${idx}`} region={donorSpan}
                                x={0} width={width} height={height} txStart={txStart} txEnd={txEnd}
                                fill={donorFill} opacity={0.5}
                            />
                        ))
                    }

                    {
                        acceptors.map((acceptorSpan, idx) => (
                            <Region key={`acceptor_intron_${idx}`} region={acceptorSpan}
                                x={0} width={width} height={height} txStart={txStart} txEnd={txEnd}
                                fill={acceptorFill} opacity={0.5}
                            />
                        ))
                    }
                </g>
            </g>
        );
    }
}

class Transcript extends React.Component {
    render() {
        let {variant, donors, acceptors, segments, width, isFlipped} = this.props;


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
                        variant={passedVariant} donors={relevantDonors} acceptors={relevantAcceptors}
                        txStart={segment.span.start} txEnd={segment.span.end}
                        x={offsets[i].x} width={offsets[i].width} height={40}
                        highlight={segment.highlighted} zoomed={false} isFlipped={isFlipped}
                    />
                );
            }
            else if (segment.type === 'intron') {
                return (
                    <Intron key={`intron_${segment.id}`} n={segment.id}
                        variant={passedVariant} donors={relevantDonors} acceptors={relevantAcceptors}
                        txStart={segment.span.start} txEnd={segment.span.end}
                        x={offsets[i].x} width={offsets[i].width} height={20}
                        highlight={segment.highlighted} isFlipped={isFlipped}
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
                <line x1={zoomLineStartLeft} y1={50} x2={zoomMargin} y2={95} stroke="#ccc" strokeWidth="3" strokeDasharray="8,3" fill="none" />
                <line x1={zoomLineStartRight} y1={50} x2={width - zoomMargin - 2.5} y2={95} stroke="#ccc" strokeWidth="3" strokeDasharray="8,3" fill="none" />

                <rect x={0} y={53} width={820} height={14} fill="rgba(255,255,255,0.6" />

                <g transform="translate(0, 10)">
                    <rect x={1} y={8} width={scale * leaderSize} height={24} fill={intronFill} stroke={intronStroke} />
                    { blocks }
                    <rect x={tailPos} y={8} width={scale * tailSize} height={24} fill={intronFill} stroke={intronStroke} />
                </g>
            </g>
        );
    }
}

class Zoom extends React.Component {
    render() {
        let { variant, donors, acceptors, segments, width, isFullyIntronic, isFlipped } = this.props;


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
                        variant={passedVariant} donors={relevantDonors} acceptors={relevantAcceptors}
                        txStart={segment.span.start} txEnd={segment.span.end}
                        x={offsets[i].x} width={offsets[i].width} height={80}
                        highlight={true} zoomed={true} isFlipped={isFlipped}
                    />
                );
            }
            else if (segment.type === 'intron') {
                return (
                    <Intron key={`intron_${i + 1}`} n={segment.id}
                        variant={passedVariant} donors={relevantDonors} acceptors={relevantAcceptors}
                        txStart={segment.span.start} txEnd={segment.span.end}
                        x={offsets[i].x} width={offsets[i].width} height={60}
                        highlight={true}  isFlipped={isFlipped}
                    />
                );
            }
        });

        return (
            <g transform="translate(0, 100)">
            { blocks }
            </g>
        );
    }
}

class Splicing extends React.Component {
    constructor(props) {
        super(props);

        this.state = {
            drawAcceptors: true,
            drawDonors: true
        };

        this.toggleDrawing = this.toggleDrawing.bind(this);
    }

    toggleDrawing(event) {
        event.persist();
        console.log(event);
        this.setState({ [event.target.name]: event.target.checked });
    }

    render() {
        let { variant } = this.props;

        let width = 800,
            info = variantInfo(variant),
            variantStart = variant.Hg38_Start,
            variantEnd = variant.Hg38_End;

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
        const siteStyle = {
            display: 'inline-block', verticalAlign: 'text-bottom',
            width: '17px', height: '17px', marginRight: '8px', opacity: 0.5,
            border: 'solid 1px black'
        };

        return (
            <div>
                <svg viewBox="-4 0 808 240" preserveAspectRatio="xMidYMid">
                    <pattern id="diagonalHatch" patternUnits="userSpaceOnUse" width="4" height="4">
                        <rect x="0" y="0" width="4" height="4" fill="#FF8888" />
                        <path d="M-1,1 l2,-2
                           M0,4 l4,-4
                           M3,5 l2,-2"
                            stroke="black" strokeWidth={1} />
                    </pattern>

                    <Transcript variant={variant} segments={segments} width={width}
                        donors={this.state.drawDonors ? meta.spliceDonors : {}}
                        acceptors={this.state.drawAcceptors ? meta.spliceAcceptors : {}}
                        isFlipped={isFlipped} />

                    <Zoom variant={variant} segments={segments.filter(x => x.highlighted)} width={width}
                        donors={this.state.drawDonors ? meta.spliceDonors : {}}
                        acceptors={this.state.drawAcceptors ? meta.spliceAcceptors : {}}
                        isFullyIntronic={isFullyIntronic} isFlipped={isFlipped} />

                    <g transform="translate(274,220)">
                        <rect x="0" fill="lightgreen" stroke="black" width="20" height="10" />
                        <text x="22" y="10">{ `Substitution (${info.changed} base${plural(info.changed)})` }</text>
                        <rect x="192" fill="url(#diagonalHatch)" stroke="black" width="20" height="10" />
                        <text x="214" y="10">{ `Deletion (${info.deleted || 0} base${plural(info.deleted)})` }</text>
                        <rect x="360" fill="#56F" stroke="black" width="20" height="10"/>
                        <text x="382" y="10">{ `Insertion (${info.inserted || 0} base${plural(info.inserted)})` }</text>
                    </g>
                </svg>

                <div style={{padding: '10px'}}>
                    <div>
                        <label>
                            <input style={{marginRight: '0.5em'}} type="checkbox" name="drawDonors" checked={this.state.drawDonors} onChange={this.toggleDrawing} />
                            <span style={{...siteStyle, backgroundColor: donorFill}} />
                            show donor sites
                        </label>
                    </div>

                    <div>
                        <label style={{display: 'inline-block', marginRight: '1em'}}>
                            <input style={{marginRight: '0.5em'}} type="checkbox" name="drawAcceptors" checked={this.state.drawAcceptors} onChange={this.toggleDrawing} />
                            <span style={{...siteStyle, backgroundColor: acceptorFill}} />
                            show acceptor sites
                        </label>
                    </div>
                </div>
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
     * @returns an array of segments, i.e. [{type: 'intron'/'exon', id, span: {start,end}, highlighted}, ...]
     */
    buildSegments(exons, variantSpan, isFlipped) {
        // disregard strandedness so we can just build some intervals
        const sortedExons = _.sortBy(
            exons.map(exon => ({
                id: exon.id,
                span: {
                    start: Math.min(exon.span.start, exon.span.end),
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
        // here we build a set of these segments, i.e. exons with introns between them
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
        const overlappingSegments = segments
            .map((segment, idx) => ({idx: idx, segment: segment}))
            .filter(({segment}) => overlaps(variantSpan, [segment.span.start, segment.span.end]));

        const [firstSeg, lastSeg] = [_.first(overlappingSegments), _.last(overlappingSegments)];

        // compute region that we'll be zooming in on in the figure
        // the region of interest is one exon before and one exon after the segments that overlap the variant
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

        // if it's antisense, we need to flip the collection so that it's ordered from exon 1 to exon N
        if (isFlipped) {
            segments.reverse();
        }

        return segments;
    }
}

module.exports = Splicing;
