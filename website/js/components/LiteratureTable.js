'use strict';

const React = require('react'),
    { Link } = require('react-router'),
    { Grid, Row, Col, ButtonToolbar, Button, DropdownButton, Table } = require('react-bootstrap'),
    backend = require('../backend'),
    _ = require('underscore');

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

/*
function limitAuthorCount(authors, maxCount) {
    let popper = (<Popover>{authors.split(";").map(e => [e, <br />])}</Popover>);
    let auth = authors.split(";");
    return (
        <span>
            { auth.slice(0, maxCount).join(";") }
            { auth.length > maxCount ?
                (<OverlayTrigger placement='bottom' overlay={popper}><span>... (more)</span></OverlayTrigger>)
                : null
            }
        </span>
    );
}
*/

function formatMatches(matches, count) {
    let ms = matches.split("|");
    ms = ms.slice(0, count);
    return ms.map(match => <div><small>... {_.unescape(match)} ...</small></div>);
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
        const authList = authors.split(";");

        return (<div style={{marginTop: '5px'}}>
            {
                this.state.expanded
                    ? `${authList.join(";")}.`
                    : (
                        <span>
                            {authList.slice(0, maxCount).join(";")}
                            ...&nbsp;
                            <a style={{cursor: 'pointer'}} onClick={() => this.setState({ expanded: true})}>(more)</a>
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

    toJSON() {
        return JSON.stringify(this.state.papers ? this.state.papers : []);
    }

    copyTable() {
        let textarea = this.refs.clipboardContent.getDOMNode();
        textarea.value = this.toTSV();
        textarea.select();
        document.execCommand('copy');
    }

    render() {
        if (!this.state.papers) {
            return (<div />);
        }
        if (!this.props.variant && !this.state.data) {
            return (<div />);
        }
        let litRows = this.state.papers.sort(pubsOrdering);
        if (this.props.maxRows) {
            litRows = litRows.slice(0, this.props.maxRows);
        }
        litRows = litRows.map(({title, authors, journal, year, mentions, pmid}) => (
            <tr>
                <td>
                    <b style={{fontSize: '16px'}}>{title}</b>
                    {/*<div>{limitAuthorCount(authors, 3)}</div>*/}

                    <AuthorList authors={authors} maxCount={3} />

                    <div style={{margin: '10px'}}>
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
                    <div className="pmid">PMID: <a href={`https://www.ncbi.nlm.nih.gov/pubmed/${pmid}`} target='_blank'>{pmid}</a></div>
                    <div>{year}</div>
                    <div>{journal}</div>
                </td>
            </tr>
        ));
        let toTSVURL = `data:text/tab-separated-values;charset=utf-8,${encodeURIComponent(this.toTSV())}`;
        let toJSONURL = `data:text/json;charset=utf-8,${encodeURIComponent(this.toJSON())}`;
        let component = (
            <div>
                <h3>Literature Search Results:</h3>
                    <Table className='nopointer literature-rows' bordered>
                        <thead>
                            <tr className="active">
                                <th>Title</th>
                                <th>Publication</th>
                            </tr>
                        </thead>
                        <tbody>
                        {litRows}
                        </tbody>
                    </Table>

                    { this.state.papers.length > this.props.maxRows
                        ? ( <div style={{textAlign: "center"}}>
                                <Link to={`/variant_literature/${this.props.variant.id}`}>
                                    View {`${this.state.papers.length - this.props.maxRows}`} more publications...
                                </Link>
                            </div>
                        ) : null
                    }

                <ButtonToolbar className='pull-right'>
                    <Button onClick={this.copyTable.bind(this)}>Copy To Clipboard</Button>
                    <DropdownButton title="Export">
                        <li><a href={toTSVURL} download='variant-literature.tsv'>tsv</a></li>
                        <li><a href={toJSONURL} download='variant-literature.json'>json</a></li>
                    </DropdownButton>
                </ButtonToolbar>

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
