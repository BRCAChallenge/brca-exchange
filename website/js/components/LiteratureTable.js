'use strict';

const React = require('react'),
    { Link } = require('react-router'),
    { Grid, Row, Col, ButtonToolbar, Button, DropdownButton, Table } = require('react-bootstrap'),
    backend = require('../backend'),
    _ = require('underscore');

import BetaTag from "./BetaTag";


// Sort first by year, then by PMID (as a proxy for a more specific date of publication?)
const pubsOrdering = function(pub1, pub2) {
    if (isNaN(pub1.year)) {
        return 1;
    }
    if (isNaN(pub2.year)) {
        return -1;
    }
    if (pub2.year - pub1.year !== 0) {
        return pub2.year - pub1.year;
    }
    if (pub2.pmid - pub1.pmid !== 0) {
        return pub2.pmid - pub1.pmid;
    }
    return 1;
};


function formatMatches(matches, count) {
    let ms = matches.slice(0, count);

    return ms.map(match =>
        // we extract entries like <<<some text here>>>, keeping the first < as an indicator that the text is highlighted
        <li>
            <small>...
            {
                _.unescape(match).split(/<<(<.*)>>>/g).map(
                    x => x.startsWith('<')
                        ? <span style={{backgroundColor: 'yellow', borderRadius: 2, border: 'solid 1px #bca723'}}>{x.slice(1)}</span>
                        : x)
            }
            ...</small>
        </li>
    );

    // return ms.map(match => <div><small>... {_.unescape(match)} ...</small></div>);
}

/**
 * Given a string of authors, returns an object with two lists: a full list of parsed-out authors and an abbreviated
 * list of authors. If the number of authors is less than maxBeforeCut, the abbreviated list will be the same as the
 * full list. If it's more, the abbreviated list will contain a maximum of maxAfterCut authors.
 * @param authors the list of semicolon-delimited authors, with each author having a comma-delimited given and first name(s).
 * @param maxBeforeCut the number of authors which triggers the abbreviated list to be shortened
 * @param maxAfterCut the maximum number of authors in the resulting abbreviated list
 * @param abbrevFirstNames whether first name(s) should be shortened to just their first characters, e.g. "Emily Leigh" -> EL
 * @return {{full: *, abbreviated: *, diff: number}}
 */
export function authorsList(authors, maxBeforeCut, maxAfterCut, abbrevFirstNames = false) {
    if (!maxAfterCut) {
        maxAfterCut = maxBeforeCut;
    }

    const authList = authors.split(";").map(x => {
        // for an author like "Rauscher, Emily A" return "Rauscher EA"
        const [given, rest] = x.trim().split(",", 2);

        if (!rest) {
            // sometimes we encounter names of orgs or consortiums, e.g., "ENIGMA"
            // there's no second part there, so we just return the first (aka 'given name') part
            return given.trim();
        }

        if (abbrevFirstNames) {
            return `${given.trim()} ${rest.trim().split(" ").map(x => x.slice(0, 1)).join("")}`;
        }
        else {
            return `${given.trim()}, ${rest.trim()}`;
        }

    });

    const abbrevList = authList.length > maxBeforeCut && authList.length > maxAfterCut
        ? authList.slice(0, maxAfterCut)
        : authList;

    return {
        full: authList,
        abbreviated: abbrevList,
        diff: authList.length - abbrevList.length // if nonzero, indicates caller should add 'et al'
    };
}


class AuthorList extends React.Component {
    constructor(props) {
        super(props);
        this.state = {
            expanded: false
        };
    }

    render() {
        const {authors, maxCount} = this.props;
        const authList = authorsList(authors, maxCount); // authors.split(";");

        return (<div style={{marginTop: '5px'}}>
            {
                this.state.expanded
                    ? `${authList.full.join("; ")}.`
                    : (
                        <span>
                            {authList.abbreviated.join("; ")}
                            ...&nbsp;
                            <a style={{cursor: 'pointer'}} onClick={() => this.setState({ expanded: true})}>
                                {`(${authList.diff} more)`}
                            </a>
                        </span>
                    )
            }
        </div>);
    }
}

AuthorList.propTypes = {
    authors: React.PropTypes.string,
    maxCount: React.PropTypes.number
};

class LiteratureTable extends React.Component {
    constructor(props) {
        super(props);
        this.state = {};
    }

    componentWillMount() {
        if (!this.props.variant) {
            backend.variant(this.props.params.id).subscribe(
                resp => {
                    return this.setState({data: resp.data, error: null});
                },
                () => { this.setState({error: 'Problem connecting to server'}); }
            );
        }

        backend.variantPapers(this.props.variant ? this.props.variant.id : this.props.params.id).subscribe(
            resp => {
                this.setState({papers: resp.data, paperError: null});
            },
            () => { this.setState({paperError: 'Problem retrieving papers'}); }
        );

    }

    toTSV() {
        const headerRow = ["title", "authors", "journal", "year", "keywords", "pmid"].join("\t");
        return this.state.papers ? this.state.papers.map(({title, authors, journal, year, keywords, pmid}) =>
            headerRow + "\n" + [title, authors, journal, year, keywords, pmid].join("\t")).join("\n") : "";
    }

    toCitation() {
        return this.state.papers
            ? this.state.papers.map(({title, authors, journal, year, keywords, pmid}) => {
                const authList = authorsList(authors, 5, 3, true);
                let abbrevList = authList.abbreviated;
                if (authList.diff > 0) {
                    abbrevList.push("et al");
                }

                return `${abbrevList.join(", ")}. (${year}) ${title} ${journal}\nPMID: ${pmid}\nKeywords: ${keywords}\n\n`;
            })
            : "";
    }

    toJSON() {
        return JSON.stringify(this.state.papers ? this.state.papers : []);
    }

    copyTable() {
        let textarea = this.refs.clipboardContent.getDOMNode();
        textarea.value = this.toCitation();
        textarea.select();
        document.execCommand('copy');
    }

    trackLitAccess(elem) {
        const targetURL = elem.target.getAttribute('href');

        if (window.ga) {
            window.ga('send', 'event', {
                eventCategory: 'Literature Outbound Link',
                eventAction: 'click',
                eventLabel: targetURL
            });
        }
    }

    render() {
        if (!this.props.variant && !this.state.data) {
            return (<div />);
        }

        // if no lit results exist, return a 'no lit results' placeholder
        let litResultsExist = false;
        let litRows = (
            <tr>
                <td colSpan={2} style={{textAlign: 'center'}}><i>No literature results for this variant.</i></td>
            </tr>
        );

        if (this.state.papers && this.state.papers.length > 0) {
            litResultsExist = true;

            litRows = this.state.papers.sort(pubsOrdering);

            if (this.props.maxRows) {
                litRows = litRows.slice(0, this.props.maxRows);
            }

            litRows = litRows.map(({title, authors, journal, year, mentions, pmid}) => (
                <tr>
                    <td>
                        <b style={{fontSize: '16px'}}>{title}</b>

                        <AuthorList authors={authors} maxCount={3} />

                        <div className="lit-text-matches">
                            <b>Text Matches:</b>
                            {
                                mentions.length ? (
                                    <ol>
                                        {formatMatches(mentions, 3)}
                                    </ol>
                                ) : null
                            }
                        </div>
                    </td>
                    <td style={{width: '12%'}}>
                        <div className="pmid">PMID: <a onClick={this.trackLitAccess} href={`https://www.ncbi.nlm.nih.gov/pubmed/${pmid}`} target='_blank'>{pmid}</a></div>
                        <div>{year}</div>
                        <div>{journal}</div>
                    </td>
                </tr>
            ));
        }

        let toTSVURL = `data:text/tab-separated-values;charset=utf-8,${encodeURIComponent(this.toTSV())}`;
        let toJSONURL = `data:text/json;charset=utf-8,${encodeURIComponent(this.toJSON())}`;

        let component = (
            (litResultsExist || !this.props.hideEmptyItems) &&
            <div>
                <h4>Literature Search Results: <BetaTag margin="0.25em" verticalALign="top" hoverText="Literature Search is a beta feature, so please beware of erroneous or missing results. You are welcome to contact us about observed errors." /></h4>
                    <Table className='nopointer literature-rows' bordered>
                        <thead>
                            <tr className="active">
                                <th>Publications Mentioning this Variant</th>
                                <th style={{width: '12%'}}>Citation Info</th>
                            </tr>
                        </thead>
                        <tbody>
                        {litRows}
                        </tbody>
                    </Table>

                    { this.state.papers && this.state.papers.length > this.props.maxRows
                        ? ( <div style={{textAlign: "center", marginBottom: '1em'}}>
                                <Link to={`/variant_literature/${this.props.variant.id}`}>
                                    View {`${this.state.papers.length - this.props.maxRows}`} more publications...
                                </Link>
                            </div>
                        ) : null
                    }

                    <div>
                        <em className="pull-left" style={{marginBottom: '1em'}}>
                            To report a false positive, or to include a paper that should be in the list, please <a href="mailto:brca-exchange-contact@genomicsandhealth.org">contact us</a>.
                        </em>

                        <ButtonToolbar className='pull-right'>
                            <Button onClick={this.copyTable.bind(this)}>Copy To Clipboard</Button>
                            <DropdownButton title="Export" className='pull-right'>
                                <li><a href={toTSVURL} download='variant-literature.tsv'>Excel (.tsv format)</a></li>
                                <li><a href={toJSONURL} download='variant-literature.json'>JSON</a></li>
                            </DropdownButton>
                        </ButtonToolbar>
                    </div>

                <textarea ref='clipboardContent' style={{padding: '0', width: '0', height: '0', marginLeft: '-99999999px' }}/>
            </div>
        );

        if (this.props.variant) {
            return component;
        } else {
            return (
                <Grid>
                    <Row>
                        <Col md={12} className="variant-literature-col">
                            <h3>{this.state.data[0]["HGVS_cDNA"]}</h3>
                            {component}
                        </Col>
                    </Row>
                </Grid>
            );
        }
    }
}

export default LiteratureTable;
