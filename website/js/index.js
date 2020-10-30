/*eslint-env browser */
/*global require: false */
'use strict';


import SourceReportsTile from "./components/SourceReportsTile";
import AlleleFrequenciesTile from "./components/AlleleFrequenciesTile";
import LiteratureTable from "./components/LiteratureTable";
import SilicoPredTile from "./components/insilicopred/SilicoPredTile";
import FunctionalAssayTile from "./components/functionalassay/FunctionalAssayTile";
import MupitStructure from './MupitStructure';

// shims for older browsers
require('babel/polyfill');
require('es5-shim');
require('es5-shim/es5-sham');

require('./favicons');
var React = require('react');
var ReactDOM = require('react-dom');
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
import debounce from 'lodash/debounce';

var brcaLogo = require('./img/BRCA-Exchange-tall-tranparent.png');
var logos = require('./logos');
var slugify = require('./slugify');

import content, {parseTooltips} from './content';
var Community = require('./Community');
var FactSheet = require('./FactSheet');
var WhyDonate = require('./WhyDonate');
var FundraisingDetails = require('./FundraisingDetails');
var Splicing = require('./Splicing');

var databaseKey = require('../databaseKey');
var util = require('./util');

var {Grid, Col, Row, Table, Button, Modal, Panel} = require('react-bootstrap');

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
var Help = require('./Help.js');

var KeyInline = require('./components/KeyInline');
var GroupHelpButton = require('./components/GroupHelpButton');

var variantPathJoin = row => _.map(databaseKey, k => encodeURIComponent(row[k])).join('@@');

if (typeof console === "undefined") {
    window.console = {
        log: function () {}
    };
}

// add react-ga tracking code to track each pageview
import ReactGA from 'react-ga';
if (window.config.analytics) {
    ReactGA.initialize(window.config.analytics, { standardImplementation: true });
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
                        <li><a href="/about/app">Mobile App</a></li>
                        <li><a href="/about/api">API</a></li>
                    </ul>
                </div>
                <div className="col-sm-2 logo-footer">
                    <img href="#" src={brcaLogo} alt="brca exchange logo" />
                </div>
                <div className="col-sm-5 right-footer">
                    <ul>
                        <li>
                            <a href="/whydonate">Donate</a>
                        </li>
                        <li>
                            <DisclaimerModal text="Disclaimer"/>
                        </li>
                        <li>
                            <a href="mailto:brca-exchange-contact@genomicsandhealth.org?subject=BRCA Exchange website">
                                Contact us
                            </a>
                        </li>
                        <li>
                            <a href="https://github.com/BRCAChallenge/brca-exchange">
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
        let currentSupporters = _.filter(logos, function(logo) {
                                    return logo.currentSupporter;
                                });
        let notCurrentSupporterLogos = _.filter(logos, function(logo) {
                             return !logo.currentSupporter;
                         });
        var notCurrentSupporterLogoItems = _.map(notCurrentSupporterLogos, ({id, logo, url}) => (
            <Col key={id} lg={4} md={6} xs={12} className="logo-item">
                <a href={url}>
                    <img id={id} src={logo} alt={id + ' logo'} />
                </a>
            </Col>
        ));
        var currentSupporterLogoItems = _.map(currentSupporters, ({id, logo, url}) => (
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
                            onSearch={this.onSearch}/>
                    </Col>
                </Row>
                <Row>
                    <div className="jumbotron homepage-jumbotron">
                        <RawHTML html={content.pages.home} />
                    </div>
                </Row>

                <Row>
                    <Col lg={4} lgOffset={1} md={8} mdOffset={2} xs={12}>
                        <div className="embed-responsive embed-responsive-16by9">
                            <iframe className="vimeo-video embed-responsive-item" src="https://player.vimeo.com/video/199396428" webkitallowfullscreen mozallowfullscreen allowFullScreen />
                        </div>
                        <div className="homepage-under-image-text-container center-block">
                            <div className="homepage-caption">
                            What is the BRCA Exchange?
                            </div>
                            <div className="homepage-caption subtext">
                            Video produced by <a href="http://www.kindealabs.com/">Kindea Labs</a>
                            </div>
                        </div>
                    </Col>

                    <Col lg={4} lgOffset={2} md={8} mdOffset={2} xs={12}>
                        <div className="embed-responsive embed-responsive-16by9">
                            <iframe className="vimeo-video embed-responsive-item" src="https://player.vimeo.com/video/351028818" webkitallowfullscreen mozallowfullscreen allowFullScreen />
                        </div>
                        <div className="homepage-under-image-text-container center-block">
                            <div className="homepage-caption">
                                How Does the BRCA Exchange Benefit Patients?
                            </div>
                            <div className="homepage-caption subtext">
                                Video produced by the <a href="https://www.ga4gh.org/">Global Alliance for Genomics and Health</a> and the <a href="https://www.broadinstitute.org/">Broad Institute</a>
                            </div>
                        </div>
                    </Col>
                </Row>
                <Row className="logo-block">
                    {notCurrentSupporterLogoItems}
                </Row>
                <Row className="logo-block">
                    <h3 className="logo-header">Currently Supported By:</h3>
                    {currentSupporterLogoItems}
                </Row>
                <Row className="logo-block">
                    <h3 className="logo-header no-margin-bottom">Consider supporting this open-source project by <Link to={'/whydonate'}>donating</Link> today.</h3>
                </Row>
            </Grid>
        );
    }
});

var About = React.createClass({
    render: function() {
        let {page} = this.props.params;
        if (page === "thisSite") {
            let currentSupporters = _.filter(logos, function(logo) {
                                        return logo.currentSupporter;
                                    });
            let notCurrentSupporterLogos = _.filter(logos, function(logo) {
                                 return !logo.currentSupporter;
                             });
            var notCurrentSupporterLogoItems = _.map(notCurrentSupporterLogos, ({id, logo, url}) => (
                <Col key={id} lg={4} md={6} xs={12} className="logo-item">
                    <a href={url}>
                        <img id={id} src={logo} alt={id + ' logo'} />
                    </a>
                </Col>
            ));
            var currentSupporterLogoItems = _.map(currentSupporters, ({id, logo, url}) => (
                <Col key={id} lg={4} md={6} xs={12} className="logo-item">
                    <a href={url}>
                        <img id={id} src={logo} alt={id + ' logo'} />
                    </a>
                </Col>
            ));
            return (
                <Grid id="main-grid" className="main-grid">
                    <Row>
                        <Col smOffset={1} sm={10}>
                            <RawHTML html={content.pages[page]} />
                        </Col>
                    </Row>
                    <Row className="logo-block">
                        {notCurrentSupporterLogoItems}
                    </Row>
                    <Row className="logo-block">
                        <h3 className="logo-header">Currently Supported By:</h3>
                        {currentSupporterLogoItems}
                    </Row>
                </Grid>
            );
        } else {
            return (
                <Grid id="main-grid" className="main-grid">
                    <Row>
                        <Col smOffset={1} sm={10}>
                            <RawHTML html={content.pages[page]} />
                        </Col>
                    </Row>
                </Grid>
            );
        }
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
    // Remove empty values from filterValues.
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
            restoringDefaults: false
        };
    },
    showVariant: function (row, event) {
          var d3TipDiv = document.getElementsByClassName('d3-tip-selection');
          if (d3TipDiv.length !== 0 && d3TipDiv[0].style.opacity !== '0') {
              d3TipDiv[0].style.opacity = '0';
              d3TipDiv[0].style.pointerEvents = 'none';
          }

          if (event.metaKey || event.altKey || event.ctrlKey || event.button === 1) {
              // the user is attempting to open the link in a new window/tab
              // (browsers vary in what they consider the 'special' tab-opening key, so we're trapping for them all)
              window.open(`/variant/${variantPathJoin(row)}`, '_blank');
          }
          else {
              // open it in the current window
              this.transitionTo(`/variant/${variantPathJoin(row)}`);
          }
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
    restoreDefaults: function(callback) {
        this.setState({restoringDefaults: true}, function() {
            this.transitionTo('/variants', null, null);

            // Callback resets filters in DataTable.
            // HACK: wrapped in setTimeout to ensure that it happens
            // after transitionTo is complete.
            setTimeout(callback, 0);
        });
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
            if (!this.state.showModal && !this.state.restoringDefaults) {
                // Don't change url if modal is open -- user is still deciding whether to change modes.
                this.transitionTo('/variants', {}, urlFromDatabase(state));
            } else if (this.state.restoringDefaults) {
                // If restoring defaults, transition to is already being called with different params.
                this.setState({restoringDefaults: false});
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
        if (this.state.restoringDefaults) {
            params.columnSelection = {};
            params.sourceSelection = {};
            params.filterValues = {};
        }
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
                    restoreDefaults={this.restoreDefaults}
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
                    restoreDefaults={this.restoreDefaults}
                    mode={this.props.mode}/>);
            message = this.renderMessage(content.pages.variantsDefault);
        }
        return (
            <Grid id="main-grid" style={{display: show ? 'block' : 'none'}}>
                {message}
                {table}
            </Grid>
        );
    },
    renderMessage: function(message) {
        return  (
			<Row>
				<Col className="jumbotron colorized-jumbo">
					{this.props.mode === 'default' && <img id='enigma-logo' src={require('./img/enigma_logo.jpeg')} />}
					<RawHTML ref='content' html={message}/>
					{this.props.mode === 'research_mode' && <Button className="btn-default" onClick={this.toggleMode}>
						Show Summary Data Only
					</Button>}
					{this.props.mode === 'default' &&
					<Button className="btn-default" onClick={() =>this.setState({showModal: true})}>
						Show Detail View
					</Button>}
					{
					    this.props.mode === 'default' && this.state.showModal &&
                        <Modal show={true} onHide={() => this.setState({ showModal: false })}>
                            <RawHTML html={content.pages.researchWarning}/>
                            <Button onClick={() => {this.toggleMode();}}>Yes</Button>
                            <Button onClick={() => this.setState({ showModal: false })}>No</Button>
                        </Modal>
					}
				</Col>
			</Row>);
    }
});

// get display name for a given key from VariantTable.js column specification,
// if we are in summary view mode, search summary view names then fall back to
// all data, otherwise go straight to all data. Finally, if key is not found, replace
// _ with space in the key and return that.
function getDisplayName(key) {
    const researchMode = (localStorage.getItem("research-mode") === 'true');
    let displayName;
    if (!researchMode) {
        displayName = columns.find(e => e.prop === key);
        displayName = displayName && displayName.title;
    }
    // we are not in summary view anymore, or key wasn't found in summary view columns
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

function isEmptyDiff(value) {
    return value === null || value.length < 1;
}

const IsoGrid = React.createClass({
    displayName: 'IsoGrid',

    // wraps contents in a masonry-managed container
    render: function () {
        const children = this.props.children;
        return (
            <div className="isogrid">
            {children}
            </div>
        );
    },

    // create masonry object to manage children's positions
    componentDidMount: function() {
        if (!this.masonry) {
            this.masonry = new Isotope('.isogrid', {
                layoutMode: 'packery',
                itemSelector: '.isogrid-item',
                packery: {
                    columnWidth: '.isogrid-sizer',
                    gutter: 0
                }
            });
        }
    },

    relayout: function(fullRefresh) {
        if (!this.masonry) {
            return;
        }

        if (fullRefresh) {
            this.masonry.reloadItems();
        }

        this.masonry.arrange();
    }
});

var VariantDetail = React.createClass({
    mixins: [Navigation, State],
    showHelp: function (event, title) {
        event.preventDefault();

        this.transitionTo(`/help#${slugify(title)}`);
    },
    getInitialState: () => ({
        hideEmptyItems: (localStorage.getItem("hide-empties") === 'true'),
        tooltips: parseTooltips(localStorage.getItem("research-mode") === 'true')
    }),
    componentWillMount: function () {
        backend.variant(this.props.params.id).subscribe(
            resp => {
                this.setState({data: resp.data, error: null});
            },
            () => { this.setState({error: 'Problem connecting to server'}); }
        );

        backend.variantReports(this.props.params.id).subscribe(
            resp => {
                // we always want reports grouped by source, so we'll do so centrally here
                const groupedReports = _.groupBy(resp.data, 'Source');

                this.setState({reports: groupedReports, error: null}, () => {
                    this.relayoutGrid();
                });
            }, () => {
                this.setState({reportError: 'Problem retrieving reports'});
                console.warn("Couldn't retrieve reports!");
            }
        );

    },
    componentWillUpdate: function(nextProps, nextState) {
        // reparse the tooltips since they're mode-specific
        if (nextProps.mode !== this.props.mode) {
            this.setState({
                tooltips: parseTooltips(nextProps.mode === 'research_mode')
            });
        }

        // ensure that we're viewing the latest version of the variant, and redirect if we're not,
        // unless the querystring param 'noRedirect=true' is specified, in which case we stay here.
        // (if someone *really* wants to link to the old variant, they may also specify 'noRedirectMsg=true'
        // to silence the warning at the top of the page, too.)
        const { data } = nextState;

        if (data && parseInt(this.props.params.id) !== data[0].id) {
            // if the variant we're requesting isn't the newest version, we need to redirect to it
            const { noRedirect } = this.getQuery();
            if (noRedirect !== "true") {
                this.replaceWith(`/variant/:id`, { id: data[0].id }, { redirectedFrom: this.props.params.id });
            }
        }
    },

    componentDidUpdate: function(prevProps, prevState) {
        if (prevProps.mode !== this.props.mode) {
            // if the mode changed, we have to relayout the page on the next available frame
            setTimeout(() => {
                this.relayoutGrid(true);
            }, 0);
        }

        // update the title only when the variant info is first populated
        if (!prevState.data && this.state.data) {
            const data = this.state.data;
            const variantVersionIdx = data.findIndex(x => x.id === parseInt(this.props.params.id));
            const variant = data[variantVersionIdx];

            document.title = `${variant['HGVS_cDNA'].split(":")[1]} (${variant['Gene_Symbol']}) - BRCA Exchange`;
        }
    },
    pathogenicityChanged: function(pathogenicityDiff) {
        return (pathogenicityDiff.added || pathogenicityDiff.removed) ? true : false;
    },
    setEmptyRowVisibility: function(hideEmptyItems) {
        localStorage.setItem('hide-empties', hideEmptyItems);

        this.setState({
            hideEmptyItems: hideEmptyItems
        }, () => {
            this.relayoutGrid();
        });
    },
    truncateData: function(field) {
        const fieldsToTruncate = ["Genomic_Coordinate_hg38", "Genomic_Coordinate_hg37"];
        if (fieldsToTruncate.indexOf(field) > -1) {
            return true;
        } else {
            return false;
        }
    },
    relayoutGrid: debounce(function(fullRefresh) {
        if (this.isogrid) {
            this.isogrid.relayout(fullRefresh);
        }
    }, 100),
    relayoutOnCollapsed: function(/* collapser */) {
        console.warn("Deprecated relayoutOnCollapsed; replace relayoutOnCollapsed handlers w/direct calls to relayoutGrid() in your collapsing components");

        // migration path:
        // 1. remove references to bootstrap's deprecated CollapsableMixin
        // 2. replace collapsible elements (i.e., ones that implement getCollapsableDOMNode,
        //    getCollapsableDimensionValue) with the Collapse element. remove those functions.
        // 3. instead of references to relayoutOnCollapsed, pass refs to this.relayoutGrid down to
        //    collapsible components.
        // 4. call this.relayoutGrid() in Collapse's onEntered, onExited events.
    },
    onChangeGroupVisibility(groupTitle, event) {
        // stop the page from scrolling to the top (due to navigating to the fragment '#')
        event.preventDefault();

        // if there's no existing key or if it's not true, this group is visible and thus collapsing
        const willBeCollapsed = localStorage.getItem("collapse-group_" + groupTitle) !== "true";
        localStorage.setItem("collapse-group_" + groupTitle, willBeCollapsed ? "true" : "false");

        // this.relayoutOnCollapsed(collapser);
        // this.relayoutGrid();
        // note: anything listening to the localstorage setting above should relayout the grid
        // itself, now that collapsing is handled by Collapse elements.
    },
    determineDiffRowColor: function(highlightRow) {
        return highlightRow ? 'danger' : '';
    },
    getPathogenicity: function(version, isReport) {
        if (isReport) {
            if (version.Source === "ClinVar") {
                return util.getFormattedFieldByProp("Clinical_Significance_ClinVar", version);
            } else {
                return util.getFormattedFieldByProp("Classification_LOVD", version);
            }
        } else {
            // Only concerned about expert pathogenicity for diff
            return util.getFormattedFieldByProp("Pathogenicity_expert", version);
        }
    },
    generateDiffRows: function(cols, data, isReports) {
        var diffRows = [];

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

                    if (_.contains(util.dateKeys, fieldName)) {
                        added = util.reformatDate(fieldDiff.added);
                        removed = util.reformatDate(fieldDiff.removed);
                    } else if (fieldDiff.field_type === "list") {
                        added = _.map(fieldDiff.added, elem => elem.replace(/_/g, " ").trim());
                        removed = _.map(fieldDiff.removed, elem => elem.replace(/_/g, " ").trim());
                    } else {
                        added = fieldDiff.added.trim();
                        removed = fieldDiff.removed.trim();
                    }

                    if (fieldName === "Summary_Evidence_ClinVar" || fieldName === "Description_ClinVar" || fieldName === "Review_Status_ClinVar") {
                        added = fieldDiff.added.replace(/_/g, " ").trim();
                        removed = fieldDiff.removed.replace(/_/g, " ").trim();
                    }

                    if (added !== null || removed !== null) {
                        if (util.isEmptyField(removed)) {
                            diffHTML.push(
                                <span>
                                    <strong>{ getDisplayName(fieldName) }: </strong>
                                    <span className='label label-success'><span className='glyphicon glyphicon-star'></span> New</span>
                                    &nbsp;{`${added}`}
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
                <tr key={i} className={this.determineDiffRowColor(highlightRow)}>
                    <td><Link to={`/release/${release.id}`}>{moment(release.date, "YYYY-MM-DDTHH:mm:ss").format("DD MMMM YYYY")}</Link></td>
                    <td>{this.getPathogenicity(version, isReports)}</td>
                    <td>{diffHTML}</td>
                </tr>
            );
        }

        return diffRows;
    },
    toggleSubmitterGroup: function(sourceName, submitter) {
        this.setState((pstate) => {
            // the key under the state collection for the visibility of this source-submitter group
            const k = `submitter-group-${sourceName}-${submitter}`;

            // if the key doesn't exist, set it to false; otherwise, invert it
            return {
                [k]: !(!pstate.hasOwnProperty(k) || pstate[k])
            };
        });
    },
    render: function () {
        const {data, error} = this.state;
        if (!data) {
            return <div />;
        }

        const variantVersionIdx = data.findIndex(x => x.id === parseInt(this.props.params.id));
        const variant = data[variantVersionIdx];
        const release = variant["Data_Release"];
        let cols, groups;

        const { redirectedFrom, noRedirectMsg } = this.getQuery();
        const redirectedFromVariant = redirectedFrom ? data.find(x => x.id === parseInt(redirectedFrom)) : null;

        if (this.props.mode === 'research_mode') {
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

        const groupTables = _.map(groups, ({ groupTitle, innerCols, reportSource, reportBinding, alleleFrequencies, inSilicoPred, innerGroups }) => {
            let rowsEmpty = 0;

            // if it's a report source (i.e. the key reportSource is defined), then we defer
            // to our custom nested-report-rendering method to generate this entire group
            if (reportSource) {
                // hide this panel if we have no reports, or specifically none for this source
                if (!this.state.reports || !this.state.reports[reportSource]) {
                    return null;
                }

                // also hide this reportSource, but with a warning, if we don't have metadata for that source
                // (we expect to have metadata for the source, since we both request the source and define its meta)
                if (!reportBinding) {
                    console.warn("Source report rendering requested for source with missing metadata: ", reportSource);
                    return null;
                }

                return (
                    <SourceReportsTile
                        key="source-reports-tile"
                        groupTitle={groupTitle}
                        sourceName={reportSource}
                        reportBinding={reportBinding}
                        submissions={this.state.reports[reportSource]}
                        onChangeGroupVisibility={this.onChangeGroupVisibility}
                        relayoutGrid={this.relayoutGrid}
                        hideEmptyItems={this.state.hideEmptyItems}
                        helpSection={reportBinding.helpKey}
                        showHelp={this.showHelp}
                        tooltips={this.state.tooltips}
                    />
                );
            }

            if (alleleFrequencies) {
                return (
                    <AlleleFrequenciesTile
                        key="allele-frequency-tile"
                        alleleFrequencyData={innerGroups}
                        groupTitle={groupTitle}
                        onChangeGroupVisibility={this.onChangeGroupVisibility}
                        relayoutGrid={this.relayoutGrid}
                        hideEmptyItems={this.state.hideEmptyItems}
                        helpSection="allele-frequency-reference-sets"
                        showHelp={this.showHelp}
                        tooltips={this.state.tooltips}
                        variant={variant}
                    />
                );
            }
            if (inSilicoPred) {
                return (
                    <SilicoPredTile
                        key="silico-pred-tile"
                        groupTitle='silico-pred-tile'
                        priors={variant.priors}
                        displayTitle={<span><i>In Silico</i> Prior Prediction (prior to considering other evidence)</span>}
                        Genomic_Coordinate_hg38={variant.Genomic_Coordinate_hg38}
                        onChangeGroupVisibility={this.onChangeGroupVisibility}
                        hideEmptyItems={this.state.hideEmptyItems}
                        relayoutGrid={this.relayoutGrid}
                        helpSection="in-silico-prior-probabilities-of-pathogenicity"
                        showHelp={this.showHelp}
                        synonymous={variant.HGVS_Protein.includes('=')}
                    />
                );
            }

            if (groupTitle === "Functional Assay Results") {
                return (
                    <FunctionalAssayTile
                        key="source-reports-tile"
                        groupTitle='functional-assay-tile'
                        score={variant.Functional_Enrichment_Score_Findlay_BRCA1_Ring_Function_Scores}
                        displayTitle="Functional Assay Results"
                        onChangeGroupVisibility={this.onChangeGroupVisibility}
                        hideEmptyItems={this.state.hideEmptyItems}
                        relayoutGrid={this.relayoutGrid}
                        helpSection="functional-assay-results"
                        showHelp={this.showHelp}
                    />
                );
            }

            // remove the BIC classification and importance fields unless the classification is 1 or 5
            if (groupTitle === 'Clinical Significance (BIC)') {
                const bicClass = variant['Clinical_classification_BIC'];

                if (bicClass !== 'Class 1' && bicClass !== 'Class 5') {
                    innerCols = innerCols.filter(x => x.prop !== 'Clinical_classification_BIC' && x.prop !== 'Clinical_importance_BIC');
                }
            } else if (groupTitle === 'Allele Counts (ExAC minus TCGA)') {
                return false;
            }

            // now map the group's columns to a list of row objects
            const rows = _.map(innerCols, (rowDescriptor) => {
                let {prop, title, noHelpLink} = rowDescriptor;
                let rowItem;

                if (prop === "Protein_Change") {
                    title = "Abbreviated AA Change";
                }

                // get mupit structures if they're available
                if (prop === "Mupit_Structure") {
                    rowItem = <MupitStructure variant={variant} prop={prop} onLoad={() => this.relayoutGrid()} />;

                    /*
                    Don't display mupit structures if they don't have an associated Amino Acid change.
                    Note that there shouldn't be mupit structures for these variants in the first place,
                    but there may be as getAminoAcidCode may change after the database is populated
                    */
                    if (util.getAminoAcidCode(variant["HGVS_Protein"]) === false) {
                        rowsEmpty += 1;
                        rowItem = false;
                    }

                    // mupit structures are not displayed if they're empty
                    if (rowItem === false) {
                        return false;
                    }

                    if (!variant[prop]) {
                        rowsEmpty += 1;
                        return false;
                    }
                } else if (prop === "HGVS_Protein_ID" && variant["HGVS_Protein"] !== null) {
                    let val = variant["HGVS_Protein"].split(":")[0];
                    variant[prop] = val;
                    rowItem = val;
                } else if (variant[prop] !== null) {
                    rowItem = util.getFormattedFieldByProp(prop, variant);
                }

                let isEmptyValue = (rowDescriptor.replace || rowDescriptor.dummy)
                    ? rowItem === false
                    : util.isEmptyField(variant[prop]);

                if (title === "Beacons") {
                    if (variant.Ref.length > 1 || variant.Alt.length > 1) {
                        isEmptyValue = true;
                    } else {
                        let websiteUrl = `https://beacon-network.org/#/search?chrom=${variant.Chr}&pos=${variant.Hg37_Start}&ref=${variant.Ref}&allele=${variant.Alt}&rs=GRCh37`;
                        rowItem = <a target="_blank" href={websiteUrl}>{websiteUrl}</a>;
                        isEmptyValue = false;
                    }
                }

                if (!isEmptyValue && prop === "CA_ID") {
                    let websiteUrl = `http://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/by_canonicalid?canonicalid=${variant.CA_ID}`;
                    rowItem = <a target="_blank" href={websiteUrl}>{variant[prop]}</a>;
                    isEmptyValue = false;
                }

                if (isEmptyValue) {
                    rowsEmpty += 1;
                    rowItem = '-';
                }

                totalRowsEmpty += rowsEmpty;
                return (
                    <tr key={prop} className={ (isEmptyValue && this.state.hideEmptyItems) ? "variantfield-empty" : "" }>
                        { rowDescriptor.tableKey !== false &&
                            (<KeyInline
                                tableKey={title} noHelpLink={noHelpLink}
                                tooltip={this.state.tooltips && prop && this.state.tooltips[slugify(prop)]}
                                onClick={(event) => this.showHelp(event, prop)}
                            />)
                        }
                        <td colSpan={rowDescriptor.tableKey === false ? 2 : null} ><span className={ this.truncateData(prop) ? "row-value-truncated" : "row-value" }>{rowItem}</span></td>
                    </tr>
                );
            });

            // check if all our rows are empty, in which case our group should be flagged as empty
            const allEmpty = rowsEmpty >= rows.length;

            if (allEmpty) {
                groupsEmpty += 1;
            }

            const tileTable = (
                <Table>
                    <tbody>
                        {rows}
                    </tbody>
                </Table>
            );

            const bicTileTable = (
                <div>
                    <div className="tile-disclaimer">
                        <div>
                            Please note that the BIC database is no longer actively curated. A copy of all BIC data has
                            been shared with several other variation databases.
                        </div>
                        <div>
                            Though the National Human Genome Research Institute recommends using ClinVar or BRCA Exchange for
                            updated information on BRCA1 and BRCA2 variants, the BIC database will be maintained to allow
                            historical studies and other uses.
                        </div>
                    </div>

                    <Table>
                        <tbody>
                            {rows}
                        </tbody>
                    </Table>
                </div>
            );

            return (
                <div key={`group_collection-${groupTitle}`} className={ (allEmpty && this.state.hideEmptyItems) || (allEmpty && groupTitle === 'CRAVAT - MuPIT 3D Protein View') ? "group-empty" : "" }>
                    <Panel
                        defaultExpanded={localStorage.getItem("collapse-group_" + groupTitle) !== "true"}
                    >
                        <Panel.Heading>
                            <Panel.Title componentClass="h3">
                                <Panel.Toggle componentClass="a" className="title"
                                    onClick={(event) => this.onChangeGroupVisibility(groupTitle, event)}
                                >
                                    {groupTitle}
                                </Panel.Toggle>
                                <GroupHelpButton onClick={(event) => { this.showHelp(event, groupTitle); return true; }} />
                            </Panel.Title>
                        </Panel.Heading>
                        <Panel.Collapse
                            onEntered={this.relayoutGrid}
                            onExited={this.relayoutGrid}
                        >
                            <Panel.Body>
                            {groupTitle === "Clinical Significance (BIC)" ? bicTileTable : tileTable}
                            </Panel.Body>
                        </Panel.Collapse>
                    </Panel>
                </div>
            );
        });

        // generates variant diff rows
        const diffRows = this.generateDiffRows(cols, data, false);

        // generates report diff rows
        if (this.state.reports !== undefined) {
            let sortedSubmissions = {'ClinVar': {}, 'LOVD': {}};

            // get all versions of clinvar submissions organized by accession number
            if (this.state.reports.hasOwnProperty('ClinVar')) {
                let clinvarSubmissions = this.state.reports.ClinVar;
                for (var i = 0; i < clinvarSubmissions.length; i++) {
                    if (clinvarSubmissions[i].Diff === null || clinvarSubmissions[i].Diff === undefined) {
                        // don't show empty diffs for reports
                        continue;
                    }
                    let key = clinvarSubmissions[i].SCV_ClinVar;
                    if (sortedSubmissions.ClinVar.hasOwnProperty(key)) {
                        sortedSubmissions.ClinVar[key].push(clinvarSubmissions[i]);
                    } else {
                        sortedSubmissions.ClinVar[key] = [clinvarSubmissions[i]];
                    }
                }
            }


            // get all versions of lovd submissions organized by submission id
            if (this.state.reports.hasOwnProperty('LOVD')) {
                let lovdSubmissions = this.state.reports.LOVD;
                for (var j = 0; j < lovdSubmissions.length; j++) {
                    if (lovdSubmissions[j].Diff === null || lovdSubmissions[j].Diff === undefined) {
                        // don't show empty diffs for reports
                        continue;
                    }
                    let key = lovdSubmissions[j].Submission_ID_LOVD;
                    if (sortedSubmissions.LOVD.hasOwnProperty(key)) {
                        sortedSubmissions.LOVD[key].push(lovdSubmissions[j]);
                    } else {
                        sortedSubmissions.LOVD[key] = [lovdSubmissions[j]];
                    }
                }
            }

            var clinvarDiffRows = _.map(sortedSubmissions.ClinVar, function(submissions) {
                let newestSubmission = submissions ? submissions[0] : '';
                let oldestSubmission = submissions ? submissions[submissions.length - 1] : '';
                const significance = util.sentenceCase(util.getFormattedFieldByProp("Clinical_Significance_ClinVar", newestSubmission)
                .replace(/(variant of unknown significance|uncertain significance)/i, 'VUS'));
                const submitter = util.abbreviatedSubmitter(util.getFormattedFieldByProp("Submitter_ClinVar", newestSubmission));
                return (
                    <Row>
                        <Col md={12} className="variant-history-col">
                            <h3>ClinVar Submission: {newestSubmission["SCV_ClinVar"]} ({submitter}; {significance})</h3>
                            <h4>Previous Versions of this Submission (since {util.reformatDate(oldestSubmission.Data_Release.date)}):</h4>
                            <Table className='variant-history nopointer' responsive bordered>
                                <thead>
                                    <tr className='active'>
                                        <th>Release Date</th>
                                        <th>Clinical Significance</th>
                                        <th>Changes</th>
                                    </tr>
                                </thead>
                                <tbody>
                                    {this.generateDiffRows(cols, submissions, true)}
                                </tbody>
                            </Table>
                            <p style={{display: this.props.mode === "research_mode" ? 'none' : 'block' }}>There may be additional changes to this variant, click "Show Detail View for this Variant" to see these changes.</p>
                        </Col>
                    </Row>
                );
            }, this);

            var lovdDiffRows = _.map(sortedSubmissions.LOVD, function(submissions) {
                let newestSubmission = submissions ? submissions[0] : '';
                let oldestSubmission = submissions ? submissions[submissions.length - 1] : '';
                const significance = util.sentenceCase(util.getFormattedFieldByProp("Classification_LOVD", newestSubmission)
                .replace(/(variant of unknown significance|uncertain significance)/i, 'VUS'));
                const submitter = util.abbreviatedSubmitter(util.getFormattedFieldByProp("Submitters_LOVD", newestSubmission));
                return (
                    <Row>
                        <Col md={12} className="variant-history-col">
                            <h3>LOVD Submission: {newestSubmission["DBID_LOVD"]} ({submitter}; {significance})</h3>
                            <h4>Previous Versions of this Submission (since {util.normalizeDateFieldDisplay(oldestSubmission.Data_Release.date)}):</h4>
                            <Table className='variant-history nopointer' responsive bordered>
                                <thead>
                                    <tr className='active'>
                                        <th>Release Date</th>
                                        <th>Clinical Significance</th>
                                        <th>Changes</th>
                                    </tr>
                                </thead>
                                <tbody>
                                    {this.generateDiffRows(cols, submissions, true)}
                                </tbody>
                            </Table>
                            <p style={{display: this.props.mode === "research_mode" ? 'none' : 'block' }}>There may be additional changes to this variant, click "Show Detail View for this Variant" to see these changes.</p>
                        </Col>
                    </Row>
                );
            }, this);
        }

        const tileSizeClasses = groupTables.length < 3
            ? `col-xs-12 col-md-${12 / groupTables.length}`
            : `col-xs-12 col-md-6 col-lg-6 col-xl-4`;
        const splicingTileSizeClassse = 'col-xs-12 col-md-12 col-lg-12 col-xl-8';

        return (error ? <p>{error}</p> :
            <Grid>
                <Row>
                    {
                        (this.props.mode !== "research_mode")
                            ? (
                                /* display new header w/genomic coordinates in expert-reviewed mode */
                                <span>
                                    <Col md={2}>
                                        <h3>Variant Details</h3>
                                    </Col>
                                    <Col md={8} className="vcenterblock">
                                        <div className='text-center Variant-detail-title' style={{textAlign: 'center'}}>
                                            <h1 style={{marginTop: 30}}>{variant.Genomic_HGVS_38 ? variant.Genomic_HGVS_38 : variant.Genomic_Coordinate_hg38}</h1>
                                            <div><i>or</i></div>
                                            <h3 style={{marginTop: 10}}>
                                                {variant['Reference_Sequence']}(<i>{variant['Gene_Symbol']}</i>){`:${variant['HGVS_cDNA'].split(":")[1]}`}

                                                {
                                                    (variant['HGVS_Protein'] && variant['HGVS_Protein'] !== "None") &&
                                                        " " + variant['HGVS_Protein'].split(":")[1]
                                                }
                                            </h3>
                                        </div>
                                    </Col>
                                </span>
                            )
                            : (
                                /* display old header if we're in research mode */
                                <Col xs={4} sm={4} smOffset={4} md={4} mdOffset={4} className="vcenterblock">
                                    <div className='text-center Variant-detail-title'>
                                        <h3>Variant Details</h3>
                                    </div>
                                </Col>
                            )
                    }
                    <Col md={2} className={(this.props.mode !== "research_mode") ? "vlowerblock" : "vcenterblock"}>
                        <div className="Variant-detail-headerbar">
                            <Button
                                onClick={this.setEmptyRowVisibility.bind(this, !this.state.hideEmptyItems)}
                                bsStyle={"default"}>
                                { this.state.hideEmptyItems ?
                                    <span>Show Empty Items</span> :
                                    <span>Hide Empty Items</span>
                                }
                            </Button>
                        </div>
                    </Col>

                    {variant['Change_Type'] === 'deleted' &&
                        (<Col xs={12} className="vcenterblock">
                            <p className='variant-message deleted-variant-message'>
                            Note: This variant has been removed from the BRCA Exchange. For reasons on why this variant was removed please see the <Link to={`/release/${release.id}`}>release notes</Link>.
                            </p>
                        </Col>)
                    }

                    {
                        (variantVersionIdx !== 0 && noRedirectMsg !== "true") && (
                          <Col xs={12} classname="vcenterblock">
                              <div className="variant-message outdated-variant-message panel panel-danger">
                                  <div className="panel-body panel-danger">
                                      <h3 style={{marginTop: 0}}>There is new data available on this variant.</h3>

                                      The data below is from {util.reformatDate(variant.Data_Release.date)} (Release {variant.Data_Release.name}). <a href={`/variant/${data[0].id}`}>Click here for updated data on this variant.</a>
                                  </div>
                              </div>
                          </Col>
                        )
                    }

                    {
                        redirectedFromVariant && (
                          <Col xs={12} classname="vcenterblock">
                              <div className="variant-message redirected-variant-msg panel panel-primary">
                                  <div className="panel-body panel-primary">
                                      <h3 style={{marginTop: 0}}>You are viewing the most recent data on this variant.</h3>

                                      The variant url you requested only has data up to {util.reformatDate(redirectedFromVariant.Data_Release.date)}. You have been automatically redirected to the newest data.<br />
                                      <a href={`/variant/${redirectedFrom}?noRedirect=true`}>
                                          Click here to view variant data from {util.reformatDate(redirectedFromVariant.Data_Release.date)} (Release {redirectedFromVariant.Data_Release.name}).
                                      </a>
                                  </div>
                              </div>
                          </Col>
                        )
                    }
                </Row>

                <Row>
                    <div className="container-fluid variant-details-body">
                        <IsoGrid ref={ (me) => { this.isogrid = me; } }>
                            <div className={`isogrid-sizer ${tileSizeClasses}`} />

                            {
                                // show the splicing vis if we're in research mode
                                this.props.mode === "research_mode" && (
                                    <Col key="splicing_vis"
                                        className={`variant-detail-group isogrid-item ${splicingTileSizeClassse}`}>
                                        <Panel
                                            collapsible={true}
                                            defaultExpanded={localStorage.getItem("collapse-group_transcript-visualization") !== "true"}
                                        >
                                            <Panel.Heading>
                                                <Panel.Title componentClass="h3">
                                                    <Panel.Toggle componentClass="a" className="title"
                                                        onClick={(event) => this.onChangeGroupVisibility("transcript-visualization", event)}
                                                    >
                                                        {`${variant['Gene_Symbol']} ${variant['HGVS_cDNA']} Transcript Visualization`}
                                                    </Panel.Toggle>
                                                    <GroupHelpButton group={"transcript-visualization"} onClick={(event) => { this.showHelp(event, "transcript-visualization"); return true; }} />
                                                </Panel.Title>
                                            </Panel.Heading>
                                            <Panel.Collapse
                                                onEntered={this.relayoutGrid}
                                                onExited={this.relayoutGrid}
                                            >
                                                <Panel.Body>
                                                    <Splicing variant={variant}
                                                        relayoutGrid={this.relayoutGrid}
                                                    />
                                                </Panel.Body>
                                            </Panel.Collapse>
                                        </Panel>
                                    </Col>
                                )
                            }

                            {
                                // we're mapping each group into a column so we can horizontally stack them
                                groupTables.map((x, i) => {
                                    return (
                                        <Col key={"group_col-" + i}
                                            className={`variant-detail-group isogrid-item ${tileSizeClasses}`}>
                                            {x}
                                        </Col>
                                    );
                                })
                            }
                        </IsoGrid>
                    </div>
                </Row>

                <Row>
                    <Col md={12} className="variant-history-col">
                        <h3>
                            {variant['Reference_Sequence']}(<i>{variant['Gene_Symbol']}</i>){`:${variant['HGVS_cDNA'].split(":")[1]}`}
                            {
                                (variant['HGVS_Protein'] && variant['HGVS_Protein'] !== "None") &&
                                " " + variant['HGVS_Protein'].split(":")[1]
                            }
                        </h3>
                    </Col>
                </Row>

                { this.props.mode === "research_mode" && (
                    <Row>
                        <Col md={12} className="variant-literature-col">
                            <LiteratureTable maxRows={10} variant={variant} hideEmptyItems={this.state.hideEmptyItems} />
                        </Col>
                    </Row>
                    )
                }

                <Row>
                    <Col md={12} className="variant-history-col">
                        <h4>Variant History:</h4>
                        <p>Variant nomenclature may change between releases, please review submission history below for further details.</p>
                        <Table className='variant-history nopointer' responsive bordered>
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
                        <p style={{display: this.props.mode === "research_mode" ? 'none' : 'block' }}>There may be additional changes to this variant, as well as changes to corresponding submissions. Click "Show Detail View for this Variant" to see these changes.</p>
                    </Col>
                </Row>

                {this.props.mode === "research_mode" ? clinvarDiffRows : ''}
                {this.props.mode === "research_mode" ? lovdDiffRows : ''}

                <Row>
                    <Col md={12} mdOffset={0}>
                        <DisclaimerModal buttonModal onToggleMode={this.props.toggleMode} text="Show Detail View for this Variant"/>
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
    componentDidUpdate() {
        let localStorageMode = (localStorage.getItem("research-mode") === "true") ? "research_mode" : "default";
        if (localStorageMode !== this.state.mode) {
            this.setMode();
        }
    },
    setMode: function () {
        this.setState({mode: (localStorage.getItem("research-mode") === 'true') ? 'research_mode' : 'default'});
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
        const path = this.getPath().slice(1);

        // logs the full path, including the hash, to google analytics (if analytics is enabled)
        if (window.config.analytics) {
            const fullHref = window.location.href;
            const origin = window.location.origin;
            const fullPathWithHash = fullHref.startsWith(origin) ? fullHref.slice(origin.length) : fullHref;
            ReactGA.ga('send', 'pageview', fullPathWithHash);
        }

        return (
            <div>
                <NavBarNew path={path} mode={this.state.mode} toggleMode={this.toggleMode}/>
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
        <Route path='whydonate' handler={WhyDonate}/>
        <Route path='fundraisingdetails' handler={FundraisingDetails}/>
        <Route path='help' handler={Help}/>
        <Route path='community' handler={Community}/>
        <Route path='signup' handler={Signup}/>
        <Route path='signin' handler={Signin}/>
        <Route path='reset_password' handler={ResetPassword}/>
        <Route path='profile' handler={Profile}/>
        <Route path='confirm/:activationCode' handler={ConfirmEmail}/>
        <Route path='reset/:resetToken' handler={ChangePassword}/>
        <Route path='variants' />
        <Route path='variant/:id' handler={VariantDetail}/>
        <Route path='variant_literature/:id' handler={LiteratureTable}/>
        <Route path='releases' handler={Releases}/>
        <Route path='release/:id' handler={Release}/>
    </Route>
);

var main = document.getElementById('main');

run(routes, HistoryLocation, (Root) => {
  ReactDOM.render(<Root/>, main);
});
