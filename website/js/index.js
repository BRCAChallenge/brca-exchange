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
require('css/bootstrap-xlgrid.css'); // adds xl, xxl, xxxl grid sizes to bootstrap 3
require('css/custom.css');
var _ = require('underscore');
var backend = require('./backend');
var {NavBarNew} = require('./NavBarNew');
var Rx = require('rx');
require('rx-dom');
var moment = require('moment');

// faisal: includes for masonry/isotope
var Isotope = require('isotope-layout');
require('isotope-packery');


var brcaLogo = require('./img/BRCA-Exchange-tall-tranparent.png');
var logos = require('./logos');
var slugify = require('./slugify');

var content = require('./content');
var Community = require('./Community');
var FactSheet = require('./FactSheet');
var {MailingList} = require('./MailingList');

var databaseKey = require('../databaseKey');

var {Grid, Col, Row, Table, Button, Modal, Panel, Glyphicon} = require('react-bootstrap');

/* FAISAL: added 'groups' collection that specifies how to map columns to higher-level groups */
var {VariantTable, ResearchVariantTable, researchModeColumns, columns, researchModeGroups, expertModeGroups} = require('./VariantTable');
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
            direction: null,
            showModal: false
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
                        <Button bsStyle="primary" className="center-block video-button" onClick={() => this.setState({ showModal: true })}>
                            <Glyphicon glyph="play-circle" />&nbsp;&nbsp;Video Overview
                        </Button>
                    </div>
                </Row>
                <Row className="logo-block">
                    {logoItems}
                </Row>
                {this.state.showModal && <Modal bsSize="large" onRequestHide={() => this.setState({ showModal: false })}>
                    <iframe className="vimeo-video" src="https://player.vimeo.com/video/199396428" className="vimeo-video" frameBorder="0" webkitallowfullscreen mozallowfullscreen allowFullScreen></iframe>
                </Modal>}
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
    // TODO: diff from defaults
    // Use REACT tools in chrome to see flow of state/props. Make sure to build
    // URL depending on mode -- i.e. if mode is default (expert), only hide non expert columns
    // and deactivate all filters. If state filters/column selections are the same as defaults and
    // if mode is research mode, check local storage to set filters/column selection.
    let {release, changeTypes, columnSelection, filterValues, sourceSelection,
         search, page, pageLength, mode, sortBy: {prop, order}} = state;
    if (mode !== "default") {
        // Default mode (expert portal) has static columns/sources.
        var hide = _.keys(_.pick(columnSelection, v => v === false));
        var hideSources = _.keys(_.pick(sourceSelection, v => v === 0));
        var excludeSources = _.keys(_.pick(sourceSelection, v => v === -1));
    } else {
        hide = '';
        hideSources = '';
        excludeSources = '';
    }
    // Remove empty values.
    clean(filterValues);
    let [filter, filterValue] = transpose(_.pairs(filterValues, v => (v !== null && v !== undefined && v !== '')));
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
    }, v => (!isEmptyVal(v)));

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
            if (this.state.showModal !== true) {
                // Don't change url if modal is open -- user is still deciding whether to change modes.
                this.transitionTo('/variants', {}, urlFromDatabase(state));
            }
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
					onRowClick={this.showVariant}
                    mode={this.props.mode}/>);
            message = this.renderMessage(content.pages.variantsResearch);
        } else {
            // Always reset column and source selections to default in expert mode.
            params.columnSelection = {};
            params.sourceSelection = {};
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
					onRowClick={this.showVariant}
                    mode={this.props.mode}/>);
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
					{this.props.mode === 'default' && <img id='enigma-logo' src={require('./img/enigma_logo.jpeg')} />}
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

const KeyInline = React.createClass({
    render() {
        const {onClick, tableKey} = this.props;
        return (
            <td className='help-target'>
                <span className="help-target-inline" onClick={onClick}>{tableKey}</span>
            </td>
        );
    }
});

const GroupHelpButton = React.createClass({
    render() {
        const {onClick} = this.props;
        return (
            <span role='button' onClick={onClick} aria-label="Help"
                className='panel-help-btn glyphicon glyphicon-question-sign'
            />
        );
    }
});

// get display name for a given key from VariantTable.js column specification,
// if we are in expert reviewed mode, search expert reviewed names then fall back to
// all data, otherwise go straight to all data. Finally, if key is not found, replace
// _ with space in the key and return that.
function getDisplayName(key) {
    const researchMode = (localStorage.getItem("research-mode") === 'true');
    let displayName;
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
    if (Array.isArray(value)) {
        value = value[0];
    }

    if (value === null || (typeof value === 'undefined')) {
        return true;
    }

    var v = value.trim();
    return v === '' || v === '-' || v === 'None';
}

function isEmptyVal(val) {
    if ((typeof val === 'string' || val instanceof String) && val.trim() === '') {
            return true;
        } else if (val === null || val === undefined) {
            return true;
        } else {
            return false;
        }
}

function isEmptyDiff(value) {
    return value === null || value.length < 1;
}

// attempts to parse the given date string using a variety of formats,
// returning the formatted result as something like '08 September 2016'.
// just returns the input if every pattern fails to match
function normalizeDateFieldDisplay(value) {
    // extend this if there are more formats in the future
    const formats = ["MM/DD/YYYY", "YYYY-MM-DD"];

    for (let i = 0; i < formats.length; i++) {
        const q = moment(value, formats[i]);

        if (q.isValid()) {
            return q.format("DD MMMM YYYY");
        }
    }

    return value;
}

// replaces commas with comma-spaces to wrap long lines better, removes blank entries from comma-delimited lists,
// and normalizes blank/null values to a single hyphen
function normalizedFieldDisplay(value) {
    if (value) {
        // replace any number of underscores with spaces
        // make sure commas, if present, wrap
        value = value
            .split(/_+/).join(" ")
            .split(",")
            .map(x => x.trim())
            .filter(x => x && x !== '-')
            .join(", ");

        // ensure that blank entries are always normalized to hyphens
        if (value.trim() === "") {
            value = "-";
        }
    }
    else {
        // similar to above, normalize blank entries to a hyphen
        value = "-";
    }

    return value;
}

function clean(obj) {
    // Removes all empty values from object.
    var propNames = Object.getOwnPropertyNames(obj);
    for (var i = 0; i < propNames.length; i++) {
        let propName = propNames[i];
        let val = obj[propName];
        if (isEmptyVal(val)) {
            delete obj[propName];
        }
    }
}

const IsoGrid = React.createClass({
    displayName: 'IsoGrid',

    // Wrapper to layout child elements passed in
    render: function () {
        const children = this.props.children;
        return (
            <div className="isogrid">
            {children}
            </div>
        );
    },

    // When the DOM is rendered, let Masonry know what's changed
    componentDidUpdate: function() {
        if (this.masonry) {
            this.masonry.reloadItems();
            this.masonry.arrange();
        }
    },

    // Set up Masonry
    componentDidMount: function() {
        if(!this.masonry) {
            // i suppose we're doing this so the children exist when we create it?
            this.masonry = new Isotope('.isogrid', {
                layoutMode: 'packery',
                itemSelector: '.isogrid-item',
                packery: {
                    columnWidth: '.isogrid-sizer',
                    gutter: 0
                }
            });
        }
        else {
            this.masonry.reloadItems();
            this.masonry.arrange();
        }
    }
});

var VariantDetail = React.createClass({
    mixins: [Navigation],
    showHelp: function (title) {
        this.transitionTo(`/help#${slugify(title)}`);
    },
    getInitialState: () => ({
        hideEmptyItems: (localStorage.getItem("hide-empties") === 'true')
    }),
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
    reformatDate: function(date) { //handles single dates or an array of dates
        if (isEmptyField(date)) {
            return date;
        }
        if (!Array.isArray(date)) {
            date = date.split(',');
        }
        return date.map(function(d) {
            return moment.utc(new Date(d)).format("DD MMMM YYYY");
        }).join();
    },
    pathogenicityChanged: function(pathogenicityDiff) {
        return (pathogenicityDiff.added || pathogenicityDiff.removed) ? true : false;
    },
    setEmptyRowVisibility: function(hideEmptyItems) {
        localStorage.setItem('hide-empties', hideEmptyItems);

        this.setState({
            hideEmptyItems: hideEmptyItems
        });
    },
    truncateData: function(field) {
        const fieldsToTruncate = ["Genomic_Coordinate_hg38", "Genomic_Coordinate_hg37", "Genomic_Coordinate_hg36"];
        if (fieldsToTruncate.indexOf(field) > -1) {
            return true;
        } else {
            return false;
        }
    },
    onChangeGroupVisibility(groupTitle, event) {
        // stop the page from scrolling to the top (due to navigating to the fragment '#')
        event.preventDefault();

        // the event target is actually the span *inside* the 'a' tag, but we need to check the 'a' tag for the
        // collapsed state
        const collapsingElem = event.target.parentElement;

        // FIXME: there must be a better way to get at the panel's state than reading the class
        // maybe we'll subclass Panel and let it handle its own visibility persistence

        // this looks silly, but at the time this method is called the item is starting to transition,
        // so its state when it's done will be the opposite of what it currently is
        const isCollapsed = !(collapsingElem.getAttribute("class") === "collapsed");
        localStorage.setItem("collapse-group_" + groupTitle, isCollapsed);

        // defer re-layout until the state change has completed
        const me = this;

        setTimeout(() => {
            // this forces a re-render after a group has expanded/collapsed, fixing the layout
            // note that 300ms just happens to be the duration of the expand/collapse animation
            // it'd be better to run the re-layout whenever the animation ends
            me.forceUpdate();
        }, 300);
    },
    generateDiffRows: function(cols, data) {
        var diffRows = [];

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
            let version = data[i];
            let diff = version.Diff;
            let release = version.Data_Release;
            let highlightRow = false;
            var diffHTML = [];

            if (diff !== null) {
                for (var j = 0; j < diff.length; j++) {
                    let fieldDiff = diff[j];
                    let fieldName = fieldDiff.field;
                    var added;
                    var removed;

                    if (fieldName === "Pathogenicity_expert") {
                        highlightRow = this.pathogenicityChanged(fieldDiff);
                    }

                    if (!_.contains(relevantFieldsToDisplayChanges, fieldName)) {
                        continue;
                    }

                    if (_.contains(dateKeys, fieldName)) {
                        added = this.reformatDate(fieldDiff.added);
                        removed = this.reformatDate(fieldDiff.removed);
                    } else if (fieldDiff.field_type === "list") {
                        added = _.map(fieldDiff.added, elem => elem.replace(/_/g, " ").trim());
                        removed = _.map(fieldDiff.removed, elem => elem.replace(/_/g, " ").trim());
                    } else {
                        added = fieldDiff.added.trim();
                        removed = fieldDiff.removed.trim();
                    }

                    if (added !== null || removed !== null) {
                        if (isEmptyField(removed)) {
                            diffHTML.push(
                                <span>
                                    <strong>{ getDisplayName(fieldName) }: </strong>
                                    <span className='label label-success'><span className='glyphicon glyphicon-star'></span> New</span>
                                    &nbsp;{added}
                                </span>, <br />
                            );
                        } else if (fieldDiff.field_type === "list") {
                            diffHTML.push(
                                <span>
                                    <strong>{ getDisplayName(fieldName) }: </strong> <br />
                                    { !isEmptyDiff(added) && `+${added}` }{ !!(!isEmptyDiff(added) && !isEmptyDiff(removed)) && ', '}{ !isEmptyDiff(removed) && `-${removed}` }
                                </span>, <br />
                            );
                        } else if (fieldDiff.field_type === "individual") {
                            diffHTML.push(
                                <span>
                                    <strong>{ getDisplayName(fieldName) }: </strong>
                                    {removed} <span className="glyphicon glyphicon-arrow-right"></span> {added}
                                </span>, <br />
                            );
                        }
                    }
                }
            }

            diffRows.push(
                <tr className={highlightRow ? 'danger' : ''}>
                    <td><Link to={`/release/${release.id}`}>{moment(release.date, "YYYY-MM-DDTHH:mm:ss").format("DD MMMM YYYY")}</Link></td>
                    <td>{version["Pathogenicity_expert"]}</td>
                    <td>{diffHTML}</td>
                </tr>
            );
        }

        return diffRows;
    },
    render: function () {
        const {data, error} = this.state;
        if (!data) {
            return <div />;
        }

        let variant = data[0],
            release = variant["Data_Release"],
            cols,
            groups;

        if (localStorage.getItem("research-mode") === 'true') {
            cols = researchModeColumns;
            groups = researchModeGroups;
        } else {
            cols = columns;
            groups = expertModeGroups;
        }

        // FAISAL: rather than directly map cols, we create a higher-level groups structure
        // the higher-level groups structure maps a subset of columns to that group
        let groupsEmpty = 0;
        let totalRowsEmpty = 0;

        const groupTables = _.map(groups, ({ groupTitle, innerCols }) => {
            let rowsEmpty = 0;

            // remove the BIC classification and importance fields unless the classification is 1 or 5
            if (groupTitle === 'Clinical Significance (BIC)') {
                const bicClass = variant['Clinical_classification_BIC'];

                if (bicClass !== 'Class 1' && bicClass !== 'Class 5') {
                    innerCols = innerCols.filter(x => x.prop !== 'Clinical_classification_BIC' && x.prop !== 'Clinical_importance_BIC');
                }
            }

            // now map the group's columns to a list of row objects
            const rows = _.map(innerCols, ({prop, title}) => {
                let rowItem;

                if (prop === "Protein_Change") {
                    title = "Abbreviated AA Change";
                }

                if (variant[prop] !== null) {
                    if (prop === "Gene_Symbol") {
                        rowItem = <i>{variant[prop]}</i>;
                    } else if (prop === "URL_ENIGMA") {
                        if (variant[prop].length) {
                            rowItem = <a target="_blank" href={variant[prop]}>link to multifactorial analysis</a>;
                        }
                    } else if (prop === "Assertion_method_citation_ENIGMA") {
                        rowItem = <a target="_blank" href="https://enigmaconsortium.org/library/general-documents/">Enigma Rules version Mar 26, 2015</a>;
                    } else if (prop === "Source_URL") {
                        if (variant[prop].startsWith("http://hci-exlovd.hci.utah.edu")) {
                            rowItem = <a target="_blank" href={variant[prop].split(',')[0]}>link to multifactorial analysis</a>;
                        }
                    } else if (prop === "Comment_on_clinical_significance_ENIGMA" || prop === "Clinical_significance_citations_ENIGMA") {
                        const pubmed = "http://ncbi.nlm.nih.gov/pubmed/";
                        rowItem = _.map(variant[prop].split(/PMID:? ?([0-9]+)/), piece =>
                            (/^[0-9]+$/.test(piece)) ? <a target="_blank" href={pubmed + piece}>PMID: {piece}</a> : piece );
                    } else if (prop === "HGVS_cDNA") {
                        rowItem = variant[prop].split(":")[1];
                    } else if (prop === "HGVS_Protein") {
                        rowItem = variant[prop].split(":")[1];
                    } else if (prop === "Date_last_evaluated_ENIGMA" && !isEmptyField(variant[prop])) {
                        // try a variety of formats until one works, or just display the value if not?
                        rowItem = normalizeDateFieldDisplay(variant[prop]);
                    } else {
                        rowItem = normalizedFieldDisplay(variant[prop]);
                    }
                } else if (prop === "HGVS_Protein_ID" && variant["HGVS_Protein"] !== null) {
                    rowItem = variant["HGVS_Protein"].split(":")[0];
                }

                const isEmptyValue = isEmptyField(variant[prop]);

                if (isEmptyValue) {
                    rowsEmpty += 1;
                    rowItem = '-';
                }

                totalRowsEmpty += rowsEmpty;

                return (
                    <tr key={prop} className={ (isEmptyValue && this.state.hideEmptyItems) ? "variantfield-empty" : "" }>
                        <KeyInline tableKey={title} onClick={() => this.showHelp(title)}/>
                        <td><span className={ this.truncateData(prop) ? "row-value-truncated" : "row-value" }>{rowItem}</span></td>
                    </tr>
                );
            });

            // check if all our rows are empty, in which case our group should be flagged as empty
            const allEmpty = rowsEmpty >= rows.length;

            if (allEmpty) {
                groupsEmpty += 1;
            }

            const header = (
                <h3>{groupTitle} <GroupHelpButton group={groupTitle} onClick={() => this.showHelp(groupTitle)} /></h3>
            );

            return (
                <div key={`group_collection-${groupTitle}`} className={ allEmpty && this.state.hideEmptyItems ? "group-empty" : "" }>
                    <Panel
                        header={header}
                        onSelect={(event) => this.onChangeGroupVisibility(groupTitle, event)}
                        collapsable={true}
                        defaultExpanded={localStorage.getItem("collapse-group_" + groupTitle) !== "true"}>
                        <Table>
                            <tbody>
                                {rows}
                            </tbody>
                        </Table>
                    </Panel>
                </div>
            );
        });

        const diffRows = this.generateDiffRows(cols, data);

        return (error ? <p>{error}</p> :
            <Grid>
                <Row>
                    <Col xs={4} sm={4} smOffset={4} md={4} mdOffset={4} className="vcenterblock">
                        <div className='text-center Variant-detail-title'>
                            <h3>Variant Detail</h3>
                            {variant['Change_Type'] === 'deleted' &&
                                (<p className='deleted text-left'>
                                    Note: This variant has been removed from the BRCA Exchange. For reasons on why this variant was removed please see the <Link to={`/release/${release.id}`}>release notes</Link>.
                                </p>)
                            }
                        </div>
                    </Col>
                    <Col xs={8} sm={4} md={4} className="vcenterblock">
                        <div className="Variant-detail-headerbar">
                            <Button
                                onClick={this.setEmptyRowVisibility.bind(this, !this.state.hideEmptyItems)}
                                bsStyle={!this.state.hideEmptyItems ? "primary" : "default"}>
                                { this.state.hideEmptyItems ?
                                    <span>show empty items</span> :
                                    <span>hide empty items</span>
                                }
                            </Button>
                        </div>
                    </Col>
                </Row>

                <Row>
                    <div className="container-fluid variant-details-body">
                        { (groupTables.length < 3) ?
                            <IsoGrid>
                                <div className={`isogrid-sizer col-xs-12 col-md-${12 / groupTables.length}`}/>
                                {
                                    groupTables.map((x, i) => {
                                        return (
                                            <Col key={"group_col-" + i} xs={12} md={12 / groupTables.length}
                                                className="variant-detail-group isogrid-item">
                                                {x}
                                            </Col>
                                        );
                                    })
                                }
                            </IsoGrid>
                            :
                            <IsoGrid>
                                <div className="isogrid-sizer col-xs-12 col-md-6 col-lg-6 col-xl-4"/>
                                {
                                    // we're mapping each group into a column so we can horizontally stack them
                                    groupTables.map((x, i) => {
                                        return (
                                            <Col key={"group_col-" + i} xs={12} md={6} lg={6}
                                                className="variant-detail-group isogrid-item col-xl-4">
                                                {x}
                                            </Col>
                                        );
                                    })
                                }
                            </IsoGrid>
                        }
                    </div>
                </Row>

                <Row>
                    <Col md={12} className="variant-history-col">
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
                                {diffRows}
                            </tbody>
                        </Table>
                        <p style={{display: this.props.mode === "research_mode" ? 'none' : 'block' }}>There may be additional changes to this variant, click "Show All Public Data on this Variant" to see these changes.</p>
                    </Col>
                </Row>
                <Row>
                    <Col md={12} mdOffset={0}>
                        <DisclaimerModal buttonModal onToggleMode={this.onChildToggleMode} text="Show All Public Data on this Variant"/>
                    </Col>
                </Row>
            </Grid>
        );
    }
});

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
        <Route path='factsheet' handler={FactSheet}/>
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
