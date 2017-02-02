/*eslint-env browser */
/*global require: false */
'use strict';

// shims for older browsers
require('babel/polyfill');
require('es5-shim');
require('es5-shim/es5-sham');

require('./favicons');
var React = require('react');
var PureRenderMixin = require('./PureRenderMixin'); // deep-equals version of PRM
var DisclaimerModal = require('./DisclaimerModal');
var RawHTML = require('./RawHTML');
require('bootstrap/dist/css/bootstrap.css');
require('font-awesome-webpack');
require('css/custom.css');
var _ = require('underscore');
var backend = require('./backend');
var {NavBarNew} = require('./NavBarNew');
var Rx = require('rx');
require('rx-dom');
var moment = require('moment');

var brcaLogo = require('./img/BRCA-Exchange-tall-tranparent.png');
var logos = require('./logos');
var slugify = require('./slugify');

var content = require('./content');
var Community = require('./Community');
var {MailingList} = require('./MailingList');

var databaseKey = require('../databaseKey');

var {Grid, Col, Row, Table, Button, Modal} = require('react-bootstrap');

var {VariantTable, ResearchVariantTable, researchModeColumns, columns} = require('./VariantTable');
var {Signup} = require('./Signup');
var {Signin, ResetPassword} = require('./Signin');
var {ConfirmEmail} = require('./ConfirmEmail');
var {ChangePassword} = require('./ChangePassword');
var {Profile} = require('./Profile');
var VariantSearch = require('./VariantSearch');
var {Navigation, State, Route, RouteHandler,
    HistoryLocation, run, DefaultRoute, Link} = require('react-router');
var {Releases, Release} = require('./Releases.js');

var navbarHeight = 70; // XXX This value MUST match the setting in custom.css

var variantPathJoin = row => _.map(databaseKey, k => encodeURIComponent(row[k])).join('@@');

if (typeof console === "undefined") {
    window.console = {
        log: function () {}
    };
}

var Footer = React.createClass({
    mixins: [PureRenderMixin],
    render: function() {
        return (
            <div className="container footer">
                <div className="col-sm-5 left-footer">
                    <ul>
                        <li><a href="/">Home</a></li>
                        <li><a href="/about/history">About</a></li>
                        <li><a href="/variants">Variants</a></li>
                        <li><a href="/help">Help</a></li>
                        <li><a href="/about/api">API</a></li>
                    </ul>
                </div>
                <div className="col-sm-2 logo-footer">
                    <img href="#" src={brcaLogo} alt="brca exchange logo" />
                </div>
                <div className="col-sm-5 right-footer">
                    <ul>
                        <li>
                            <DisclaimerModal text="Disclaimer"/>
                        </li>
                        <li>
                            <a href="mailto:brca-exchange-contact@genomicsandhealth.org?subject=BRCA Exchange website">
                                Contact us
                            </a>
                        </li>
                        <li>
                            <a href="https://github.com/BD2KGenomics/brca-website">
                                Source code
                            </a>
                        </li>
                    </ul>
                </div>
            </div>
        );
    }
});

var Home = React.createClass({
    mixins: [Navigation],
    getInitialState() {
        return {
            index: 0,
            direction: null
        };
    },

    onSearch(value) {
        this.transitionTo('/variants', null, {search: value});
    },
    render: function() {
        var {suggestions} = this.props;
        var logoItems = _.map(logos, ({id, logo, url}) => (
            <Col key={id} lg={4} md={6} xs={12} className="logo-item">
                <a href={url}>
                    <img id={id} src={logo} alt={id + ' logo'} />
                </a>
            </Col>
        ));
        var videoStyles = {
            padding: '56.25% 0 0 0',
            position: 'relative'
        };
        var videoWrapperStyles = {
            height: '100%',
            left: '0',
            position: 'absolute',
            top: '0',
            width: '100%',
        };
        return (
            <Grid id="main-grid" className='home'>
                <Row>
                    <Col smOffset={2} sm={8}>
                        <VariantSearch
                            id='home-search'
                            suggestions={suggestions}
                            onSearch={this.onSearch}/>
                    </Col>
                </Row>
                <Row>
                    <div className="jumbotron">
                        <RawHTML html={content.pages.home} />
                    </div>
                </Row>
                <Row>
                    <Col md={8} mdOffset={2} className="text-center video-container">
                        <iframe src="https://player.vimeo.com/video/199396428" width="640" height="360" frameBorder="0" webkitallowfullscreen mozallowfullscreen allowFullScreen></iframe>
                        <div className="wistia_responsive_padding video-container" style={videoStyles}><div className="wistia_responsive_wrapper" style={videoWrapperStyles}><iframe src="//fast.wistia.net/embed/iframe/dl2zrjowtv?videoFoam=true" title="Wistia video player" allowTransparency="true" frameBorder="0" scrolling="no" className="wistia_embed" name="wistia_embed" allowFullScreen mozallowfullscreen webkitallowfullscreen oallowfullscreen msallowfullscreen width="100%" height="100%"></iframe></div></div>
                        <script src="//fast.wistia.net/assets/external/E-v1.js" async></script>
                    </Col>
                </Row>
                <Row className="logo-block">
                    {logoItems}
                </Row>
            </Grid>
        );
    }
});

var About = React.createClass({
    render: function() {
        var {page} = this.props.params;

        return (
            <Grid id="main-grid" className="main-grid">
                <Row>
                    <Col md={8} mdOffset={2}>
                        <RawHTML html={content.pages[page]} />
                    </Col>
                </Row>
            </Grid>
        );
    }
});

var Help = React.createClass({
    mixins: [State],
    componentDidMount: function () {
        var fragment = slugify(window.location.hash.slice(1));
        if (fragment !== '') {
            setTimeout(function () {
                var el = document.getElementById(fragment);
                if (el) {
                    window.scrollTo(0, el.getBoundingClientRect().top - navbarHeight);
                }
            }, 0);
        }
    },
    render: function () {
        var fragment = slugify(window.location.hash.slice(1));
        var helpContent;
        if (localStorage.getItem("research-mode") === 'true') {
            helpContent = content.pages.helpResearch;
        } else {
            helpContent = content.pages.help;
        }
        return (
            <Grid id="main-grid" className="help">
                {fragment === '' ? null :
                    <style>{`#${fragment} { animation-name: emphasis; animation-duration: 10s; } `}</style>}
                <Row>
                    <Col md={8} mdOffset={2}>
                        <RawHTML ref='content' html={helpContent}/>
                    </Col>
                </Row>
            </Grid>
        );
    }
});

function toNumber(v) {
    return _.isString(v) ? parseInt(v) : v;
}

function databaseParams(paramsIn) {
    var {filter, filterValue, hide, hideSources, excludeSources, orderBy, order, search = '', changeTypes} = paramsIn;
    var numParams = _.mapObject(_.pick(paramsIn, 'page', 'pageLength', 'release'), toNumber);
    var sortBy = {prop: orderBy, order};
    var columnSelection = _.object(hide, _.map(hide, _.constant(false)));
    var sourceSelection = {..._.object(hideSources, _.map(hideSources, _.constant(0))),
                           ..._.object(excludeSources, _.map(excludeSources, _.constant(-1)))};
    var filterValues = _.object(filter, filterValue);
    return {changeTypes, search, sortBy, columnSelection, sourceSelection, filterValues, hide, ...numParams};
}

var transpose = a => _.zip.apply(_, a);

function urlFromDatabase(state) {
    // Need to diff from defaults. The defaults are in DataTable.
    // We could keep the defaults here, or in a different module.
    var {release, changeTypes, columnSelection, filterValues, sourceSelection,
            search, page, pageLength, sortBy: {prop, order}} = state;
    var hide = _.keys(_.pick(columnSelection, v => v === false));
    var hideSources = _.keys(_.pick(sourceSelection, v => v === 0));
    var excludeSources = _.keys(_.pick(sourceSelection, v => v === -1));
    var [filter, filterValue] = transpose(_.pairs(_.pick(filterValues, v => v === true)));
    return _.pick({
        release,
        changeTypes,
        search: search === '' ? null : backend.trimSearchTerm(search),
        filter,
        filterValue,
        page: page === 0 ? null : page,
        pageLength: pageLength === 20 ? null : pageLength,
        orderBy: prop,
        order,
        hideSources: hideSources,
        excludeSources: excludeSources,
        hide: hide.length === 0 ? null : hide
    }, v => v != null);

}

var Database = React.createClass({
    // Note this is not a pure component because of the calls to
    // getQuery().
    mixins: [Navigation, State],
    getInitialState: function () {
        return {
            showModal: false,
        };
    },
    showVariant: function (row) {
          var d3TipDiv = document.getElementsByClassName('d3-tip-selection');
          if (d3TipDiv.length !== 0 && d3TipDiv[0].style.opacity !== '0') {
              d3TipDiv[0].style.opacity = '0';
              d3TipDiv[0].style.pointerEvents = 'none';
          }
          this.transitionTo(`/variant/${variantPathJoin(row)}`);
    },
    showHelp: function (title) {
        var d3TipDiv = document.getElementsByClassName('d3-tip-selection');
        if (d3TipDiv.length !== 0 && d3TipDiv[0].style.opacity !== '0') {
            d3TipDiv[0].style.opacity = '0';
            d3TipDiv[0].style.pointerEvents = 'none';
        }
        this.transitionTo(`/help#${slugify(title)}`);
    },
    componentDidMount: function () {
        var q = this.urlq = new Rx.Subject();
        this.subs = q.debounce(500).subscribe(this.onChange);
    },
    componentWillUnmount: function () {
        this.subs.dispose();
    },
    // XXX An oddity of the state flow here: we update the url when table settings
    // change, so the page can be bookmarked, and forward/back buttons work. We
    // do it on a timeout so we don't generate history entries for every keystroke,
    // which would be bad for the user. Changing the url causes a re-render, passing
    // in new props, which causes DataTable to overwrite its state with the
    // same state that caused us to update the url. It's a bit circular.
    // It would be less confusing if DataTable did not hold these params in state,
    // but just read them from props, and all updates to the props occurred via
    // transitionTo(). Consider for a later refactor.
    onChange: function (state) {
        if (this.props.show) {
            var d3TipDiv = document.getElementsByClassName('d3-tip-selection');
            if (d3TipDiv.length !== 0 && d3TipDiv[0].style.opacity !== '0') {
                d3TipDiv[0].style.opacity = '0';
                d3TipDiv[0].style.pointerEvents = 'none';
            }
            this.transitionTo('/variants', {}, urlFromDatabase(state));
        }
    },
    toggleMode: function () {
        this.props.toggleMode();
        this.setState({ showModal: false });
    },
    render: function () {
        var {show} = this.props,
            params = databaseParams(this.getQuery());
        // XXX is 'keys' used?
        var table, message;
        if (this.props.mode === 'research_mode') {
            table = (
				<ResearchVariantTable
					ref='table'
					initialState={params}
					{...params}
					fetch={backend.data}
					fetchLollipop={backend.lollipopData}
					url={backend.url}
					onChange={s => this.urlq.onNext(s)}
					onToggleMode={this}
					keys={databaseKey}
					onHeaderClick={this.showHelp}
					onRowClick={this.showVariant}/>);
            message = this.renderMessage(content.pages.variantsResearch);
        } else {
            params.columnSelection = {};
            table = (
				<VariantTable
					ref='table'
					initialState={params}
					{...params}
					fetch={backend.data}
					fetchLollipop={backend.lollipopData}
					url={backend.url}
					onChange={s => this.urlq.onNext(s)}
					onToggleMode={this}
					keys={databaseKey}
					onHeaderClick={this.showHelp}
					onRowClick={this.showVariant}/>);
            message = this.renderMessage(content.pages.variantsDefault);
        }
        return (
            <Grid fluid={true} id="main-grid" style={{display: show ? 'block' : 'none'}}>
                {message}
                {table}
            </Grid>
        );
    },
    renderMessage: function(message) {
        return  (
			<Row>
				<Col sm={10} smOffset={1}  className="alert alert-warning">
					{this.props.mode === 'default' && <img id='enigma-logo' src={require('./img/enigma_logo.png')} />}
					<RawHTML ref='content' html={message}/>
					{this.props.mode === 'research_mode' && <Button className="btn-small" onClick={this.toggleMode}>
						Show Expert Reviewed Data Only
					</Button>}
					{this.props.mode === 'default' &&
					<Button className="btn-small" onClick={() =>this.setState({showModal: true})}>
						Show All Public Data
					</Button>}
					{this.props.mode === 'default' && this.state.showModal &&
					<Modal onRequestHide={() => this.setState({ showModal: false })}>
						<RawHTML html={content.pages.researchWarning}/>
						<Button onClick={() => {this.toggleMode();}}>Yes</Button>
						<Button onClick={() => this.setState({ showModal: false })}>No</Button>
					</Modal>}
				</Col>
			</Row>);
    }
});

var Key = React.createClass({
    render() {
        var {onClick, tableKey} = this.props;
        return (
             <td className='help-target'>
                {tableKey}
                <span className="text-nowrap">
                    <span role='button' onClick={onClick}
                        className='help glyphicon glyphicon-question-sign superscript'/>
                </span>
             </td>
        );
    }
});

// get display name for a given key from VariantTable.js column specification,
// if we are in expert reviewed mode, search expert reviewed names then fall back to
// all data, otherwise go straight to all data. Finally, if key is not found, replace
// _ with space in the key and return that.
function getDisplayName(key) {
    var researchMode = (localStorage.getItem("research-mode") === 'true');
    var displayName;
    if (!researchMode) {
        displayName = columns.find(e => e.prop === key);
        displayName = displayName && displayName.title;
    }
    // we are not in expert reviewed more, or key wasn't found in expert reviewed columns
    if (displayName === undefined) {
        displayName = researchModeColumns.find(e => e.prop === key);
        displayName = displayName && displayName.title;
    }
    // key was not found at all
    if (displayName === undefined) {
        displayName = key.replace(/_/g, " ");
    }
    return displayName;
}

// test for the various forms of blank fields
function isEmptyField(value) {
    var v = value.trim();
    return v === '' || v === '-' || v === 'None';
}

var VariantDetail = React.createClass({
    mixins: [Navigation],
    showHelp: function (title) {
        this.transitionTo(`/help#${slugify(title)}`);
    },
    getInitialState: () => ({}),
    componentWillMount: function () {
        backend.variant(this.props.params.id).subscribe(
            resp => {
                return this.setState({data: resp.data, error: null});
            },
            () => { this.setState({error: 'Problem connecting to server'}); });
    },
    onChildToggleMode: function() {
        this.props.toggleMode();
	 this.forceUpdate();
    },
    reformatDate: function(date) { //handles single dates or comma separated dates
        var dates = date.split(',');
        return dates.map(function(date) {
            return moment.utc(new Date(date)).format("DD MMMM YYYY");
        }).join();
    },
    render: function () {
        var {data, error} = this.state;
        if (!data) {
            return <div></div>;
        }

        var variant = data[0],
            release = variant["Data_Release"],
            cols;
        if (localStorage.getItem("research-mode") === 'true') {
            cols = researchModeColumns;
        } else {
            cols = columns;
        }
        var rows = _.map(cols, ({prop, title}) => {
            var rowItem;
            if (prop === "Protein_Change") {
                title = "Abbreviated AA Change";
            }
            if (variant[prop] != null) {
                if (prop === "Gene_Symbol") {
                    rowItem = <i>{variant[prop]}</i>;
                }
                else if (prop === "URL_ENIGMA") {
                    if (variant[prop].length) {
                        rowItem = <a target="_blank" href={variant[prop]}>link to multifactorial analysis</a>;
					}
                } else if (prop === "Assertion_method_citation_ENIGMA") {
                    rowItem = <a target="_blank" href="https://enigmaconsortium.org/library/general-documents/">Enigma Rules version Mar 26, 2015</a>;
// this will be used in All Data display
/*                } else if (prop == "Source_URL") {
                    var url_count = 0;
                    rowItem = _.map(variant[prop].split(','), url => (url.length != 0) && (<span><a key={"Source_URL"+(url_count++)} target="_blank" href={url}>link to multifactorial analysis ({url_count})</a><br /></span>));
  */
                } else if (prop === "Source_URL") {
                    if (variant[prop].startsWith("http://hci-exlovd.hci.utah.edu")) {
                        rowItem = <a target="_blank" href={variant[prop].split(',')[0]}>link to multifactorial analysis</a>;
					}
                } else if (prop === "Comment_on_clinical_significance_ENIGMA" || prop === "Clinical_significance_citations_ENIGMA") {
                    var pubmed = "http://ncbi.nlm.nih.gov/pubmed/";
                    rowItem = _.map(variant[prop].split(/PMID:? ?([0-9]+)/), piece =>
                        (/^[0-9]+$/.test(piece)) ? <a target="_blank" href={pubmed + piece}>PMID: {piece}</a> : piece );
                } else if (prop === "HGVS_cDNA") {
                    rowItem = variant[prop].split(":")[1];
                } else if (prop === "HGVS_Protein") {
                    rowItem = variant[prop].split(":")[1];
                } else if (prop === "Date_last_evaluated_ENIGMA" && !isEmptyField(variant[prop])) {
                    rowItem = moment(variant[prop], "MM/DD/YYYY").format("DD MMMM YYYY");
                } else {
                    rowItem = variant[prop];
                }
            } else if (prop === "HGVS_Protein_ID" && variant["HGVS_Protein"] != null) {
                rowItem = variant["HGVS_Protein"].split(":")[0];
            }
            return (
				<tr key={prop}>
					<Key tableKey={title} columns={cols} onClick={() => this.showHelp(title)}/>
					<td><span className="row-wrap">{rowItem}</span></td>
				</tr>);
        });

        var versionRows = [];
        // which keys represent comma seperated lists? for filtering out re-ordered lists
        // TODO: list keys now handled on the server, remove this code?
        var listKeys = [
            "Pathogenicity_all",
            "Submitter_ClinVar",
            "Method_ClinVar",
            "Source",
            "Date_Last_Updated_ClinVar",
            "Source_URL",
            "SCV_ClinVar",
            "Clinical_Significance_ClinVar",
            "Allele_Origin_ClinVar"
        ];

        // keys that contain date values that need reformatting for the ui
        var dateKeys = [
            "Date_Last_Updated_ClinVar",
            "Date_last_evaluated_ENIGMA"
        ];

        // In research_mode, only show research_mode changes.
        var relevantFieldsToDisplayChanges = cols.map(function(col) {
            return col.prop;
        });

        for (var i = 0; i < data.length; i++) {
            let version = data[i],
                changes = [],
                release = version["Data_Release"],
                highlightRow = false;

            // if this is not the oldest version of the variant, diff them
            // TODO: fix this diff following diff updates in pipeline.
            if (i < data.length - 1) {
                let previous = data[i + 1];
                if (version["Pathogenicity_expert"] !== previous["Pathogenicity_expert"]) {
                    highlightRow = true;
                }
                for (var key in version) {
                    if (relevantFieldsToDisplayChanges.indexOf(key) === -1) {
                        // Do not display changes for these fields.
                        continue;
                    } else if (isEmptyField(version[key].toString() && isEmptyField(previous[key]).toString())) {
                        // If both are empty, there is no change.
                        continue;
                    }

                    if (_.contains(dateKeys, key) && !isEmptyField(version[key])) {
                        version[key] = this.reformatDate(version[key]);
                        previous[key] = this.reformatDate(previous[key]);
                        if (version[key] === previous[key]) {
                            continue;
                        }
                    }

                    if (!_.contains(["Data_Release", "Change_Type", "id", "Synonyms"], key) && version[key] !== previous[key]) {
                        let versionDisplay = isEmptyField(version[key].toString()) ? <span className='empty'></span> : version[key].toString();
                        if (isEmptyField(previous[key].toString())) {
                            changes.push(
                                <span>
                                    <strong>{ getDisplayName(key) }: </strong>
                                    <span className='label label-success'><span className='glyphicon glyphicon-star'></span> New</span>
                                    &nbsp;{ versionDisplay }
                                </span>, <br />
                            );
                        }
                        else if (_.contains(listKeys, key)) {
                            let delimiter = key === "Pathogenicity_all" ? ';' : ',';
                            let trimmedVersion = _.map(version[key].split(delimiter), elem => elem.replace(/_/g, " ").trim());
                            let trimmedPrevious = _.map(previous[key].split(delimiter), elem => elem.replace(/_/g, " ").trim());

                            if (key === "Pathogenicity_all") {
                                let sortSublists = function(elem) {
                                    var submitter = elem.match(/\([a-zA-Z]*\)/)[0];
                                    return elem.replace(submitter, '').trim().split(',').sort().join(',');
                                };
                                trimmedVersion = _.map(trimmedVersion, sortSublists);
                                trimmedPrevious = _.map(trimmedPrevious, sortSublists);
                            }

                            let added = _.map(_.difference(trimmedVersion, trimmedPrevious), elem => `+${elem}`);
                            let deleted = _.map(_.difference(trimmedPrevious, trimmedVersion), elem => `-${elem}`);

                            if (added.length || deleted.length) {
                                changes.push(
                                    <span>
                                        <strong>{ getDisplayName(key) }: </strong> <br />
                                        { added.join(', ') }{ !!(added.length && deleted.length) && ', '}{ deleted.join(', ') }
                                    </span>, <br />
                                );
                            }
                        } else {
                            // exLOVD citation format changed, heuristic for matching: first word (i.e. first author) same -> ignore
                            if (key === "Literature_source_exLOVD" &&
                                version[key].trim().split(' ')[0] === previous[key].trim().split(' ')[0]) {
                                continue;
                            }
                            // ENIGMA rules URL's were previously broken, ignore the change from broken to correct links
                            else if (key === "Assertion_method_citation_ENIGMA") {
                                continue;
                            }
                            // Text change that actually means the same thing
                            else if (key === "Pathogenicity_expert" && version[key] === "Not Yet Reviewed" && previous[key] === "Not Yet Classified") {
                                continue;
                            }
                            changes.push(
                                <span>
                                    <strong>{ getDisplayName(key) }: </strong>
                                    {previous[key].toString()} <span className="glyphicon glyphicon-arrow-right"></span> {versionDisplay}
                                </span>, <br />
                            );
                        }
                    }
                }
            }
            versionRows.push(
                <tr className={highlightRow ? 'danger' : ''}>
                    <td><Link to={`/release/${release.id}`}>{moment(release.date, "YYYY-MM-DDTHH:mm:ss").format("DD MMMM YYYY")}</Link></td>
                    <td>{version["Pathogenicity_expert"]}</td>
                    <td>{changes}</td>
                </tr>
            );
        }

        return (error ? <p>{error}</p> :
            <Grid>
                <Row>
                    <Col md={8} mdOffset={2}>
                        <div className='text-center Variant-detail-title'>
                            <h3>Variant Detail</h3>
                            {variant['Change_Type'] === 'deleted' &&
                                (<p className='deleted text-left'>
                                    Note: This variant has been removed from the BRCA Exchange. For reasons on why this variant was removed please see the <Link to={`/release/${release.id}`}>release notes</Link>.
                                </p>)
                            }
                        </div>
                    </Col>
                </Row>
                <Row>
                    <Col md={8} mdOffset={2}>
                        <Table striped bordered>
                            <tbody>
                                {rows}
                            </tbody>
                        </Table>
                    </Col>
                </Row>
                <Row>
                    <Col md={8} mdOffset={2}>
                        <h3>{variant["HGVS_cDNA"]}</h3>
                        <h4>Previous Versions of this Variant:</h4>
                        <Table className='variant-history' bordered>
                            <thead>
                                <tr className='active'>
                                    <th>Release Date</th>
                                    <th>Clinical Significance</th>
                                    <th>Changes</th>
                                </tr>
                            </thead>
                            <tbody>
                                {versionRows}
                            </tbody>
                        </Table>
                        <p style={{display: this.props.mode === "research_mode" ? 'none' : 'block' }}>There may be additional changes to this variant, click "Show All Public Data on this Variant" to see these changes.</p>
                    </Col>
                </Row>
                <Row>
                    <Col md={8} mdOffset={2}>
                        <DisclaimerModal buttonModal onToggleMode={this.onChildToggleMode} text="Show All Public Data on this Variant"/>
                    </Col>
                </Row>
            </Grid>
        );
    }
});

// XXX implement in server
//var dontSuggest = [
//    'Assertion_method_citation',
//    'URL'
//];

//var flatmap = (coll, fn) => _.flatten(_.map(coll, fn), true);
//var minSuggestion = 3; // minimum length of string to use in autocomplete
//var rowWords = row => flatmap(_.values(_.omit(row, dontSuggest)),
//        v => v.toLowerCase().split(/\s+/));

// Pull out interesting strings from the data, for use in
// auto-completion.
//function getSuggestions(data) {
//    return _.uniq(flatmap(data, row =>
//                _.filter(rowWords(row), w => w.length >= minSuggestion)).sort(),
//            true);
//}

var Application = React.createClass({
    mixins: [State],
    onChildToggleMode: function() {
        this.toggleMode();
    },
    getInitialState: function () {
        return {
            mode: (localStorage.getItem("research-mode") === 'true') ? 'research_mode' : 'default',
        };
    },
    toggleMode: function () {
        if (this.state.mode === 'research_mode') {
            localStorage.setItem('research-mode', false);
            this.setState({mode: 'default'});
        } else {
            localStorage.setItem('research-mode', true);
            this.setState({mode: 'research_mode'});
        }
    },
    render: function () {
        var path = this.getPath().slice(1);
        return (
            <div>
                <NavBarNew path={path} mode={this.state.mode} />
                <RouteHandler toggleMode={this.onChildToggleMode} mode={this.state.mode} />
                <Database
                    mode={this.state.mode}
                    toggleMode={this.onChildToggleMode}
                    show={path.indexOf('variants') === 0} />
                <Footer />
            </div>
        );
    }
});

var routes = (
    <Route handler={Application}>
        <DefaultRoute handler={Home}/>
        <Route path='about/:page' handler={About}/>
        <Route path='help' handler={Help}/>
        <Route path='community' handler={Community}/>
        <Route path='signup' handler={Signup}/>
        <Route path='signin' handler={Signin}/>
        <Route path='mailinglist' handler={MailingList}/>
        <Route path='reset_password' handler={ResetPassword}/>
        <Route path='profile' handler={Profile}/>
        <Route path='confirm/:activationCode' handler={ConfirmEmail}/>
        <Route path='reset/:resetToken' handler={ChangePassword}/>
        <Route path='variants' />
        <Route path='variant/:id' handler={VariantDetail}/>
        <Route path='releases' handler={Releases}/>
        <Route path='release/:id' handler={Release}/>
    </Route>
);

var main = document.getElementById('main');

run(routes, HistoryLocation, (Root) => {
  React.render(<Root/>, main);
});
