/*eslint-env browser */
/*global require: false */
'use strict';

import SourceReportsTile from "./components/SourceReportsTile";
import AlleleFrequenciesTile from "./components/AlleleFrequenciesTile";
import LiteratureTable from "./components/LiteratureTable";
import SilicoPredTile from "./components/insilicopred/SilicoPredTile";
import FunctionalAssayTile from "./components/functionalassay/FunctionalAssayTile";
import ComputationalPredictionTile from "./components/computationalprediction/ComputationalPredictionTile";
import ProvisionalEvidenceTile from "./components/ProvisionalEvidenceTile";
import MupitStructure from './MupitStructure';

require('./favicons');
var React = require('react');
import { createRoot } from 'react-dom/client';
//var PureRenderMixin = require('./PureRenderMixin'); // deep-equals version of PRM
import DisclaimerModal from './DisclaimerModal';
import RawHTML from './RawHTML';

// Keep your existing CSS includes
import 'bootstrap/dist/css/bootstrap.min.css';
import 'font-awesome/css/font-awesome.min.css';
//import 'css/bootstrap-xlgrid.css'; // adds xl, xxl, xxxl grid sizes to bootstrap 3
import 'css/custom.css';

var _ = require('underscore');
var backend = require('./backend');
import NavBarNew from './NavBarNew';
// RxJS 6+ imports
import { Subject } from 'rxjs';
import { debounceTime } from 'rxjs/operators';
var moment = require('moment');
import DonationBar from './components/DonationBar';

// masonry/isotope
var Isotope = require('isotope-layout');
require('isotope-packery');
import debounce from 'lodash/debounce';

var brcaLogo = require('./img/BRCA-Exchange-tall-tranparent.png');
var logos = require('./logos');
var slugify = require('./slugify');

import content, {parseTooltips} from './content';
import Community from './Community';
import FactSheet from './FactSheet';
import WhyDonate from './WhyDonate';
import FundraisingDetails from './FundraisingDetails';
//import {Splicing} from'./Splicing';
import { withRouter, BrowserRouter, Switch, Route, Link } from 'react-router-dom';

//var databaseKey = require('../databaseKey');
//var util = require('./util');

// React-Bootstrap v2+ imports (Bootstrap 5)
import { Container as Grid, Col, Row, Table, Button, Modal, Card, Collapse } from 'react-bootstrap';

/* FAISAL: added 'groups' collection that specifies how to map columns to higher-level groups */
//var {VariantTable, ResearchVariantTable, researchModeColumns, columns, researchModeGroups, expertModeGroups} = require('./VariantTable');
//var {Signup} = require('./Signup');
//var {Signin, ResetPassword} = require('./Signin');
//var {ConfirmEmail} = require('./ConfirmEmail');
//var {ChangePassword} = require('./ChangePassword');
//var {Profile} = require('./Profile');
import VariantSearch from './VariantSearch';
import { Releases, Release } from './Releases.js';
import Help from './Help.js';

//var KeyInline = require('./components/KeyInline');
//var GroupHelpButton = require('./components/GroupHelpButton');

//var variantPathJoin = row => _.map(databaseKey, k => encodeURIComponent(row[k])).join('@@');

if (typeof console === "undefined") {
    window.console = {
        log: function () {}
    };
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

class Footer extends React.PureComponent {
    render() {
        return (
            <div className="container footer">
		<div className = "row">
                    <div className="col-sm-5 left-footer">
                    <ul>
                        <li><a href="/">Home</a></li>
                        <li><a href="/about/history">About</a></li>
                        <li><a href="/variants">Variants</a></li>
                        <li><a href="/about/api">API</a></li>
                        <li><a href="https://brcaexchange.org/blog">Blog</a></li>
                    </ul>
                    </div>
                    <div className="col-sm-2 logo-footer">
                        <img src={brcaLogo} alt="brca exchange logo" />
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
                            <a href="https://github.com/BRCAChallenge/brca-exchange">
                                Source code
                            </a>
                        </li>
                        <li>
                            <a href="mailto:brcaexchange@gmail.com?subject=BRCA Exchange website">
                                Contact us
                            </a>
                        </li>
                        <li><a href="/help">Help</a></li>
                    </ul>
                    </div>
		</div>
            </div>
        );
    }
}

class HomeRaw extends React.Component {
    constructor(props) {
        super(props);
	this.state = {
            index: 0,
            direction: null,
            showModal: false
        };
	this.onSearch = this.onSearch.bind(this);
    }
    onSearch(value) {
	const query = value ? `?search=${encodeURIComponent(value)}` : '';
        this.props.history.push(`/variants${query}`);
    }
    render() {
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
                    <Col sm={{ span: 8, offset: 2 }}>
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
                    <Col lg={{ span: 4, offset: 1 }} md={{ span: 8, offset: 2 }} xs={12}>
                        <div className="ratio ratio-16x9">
                            <iframe className="vimeo-video" src="https://player.vimeo.com/video/199396428" allowFullScreen />
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

                    <Col lg={{ span: 4, offset: 2 }} md={{ span: 8, offset: 2 }} xs={12}>
                        <div className="ratio ratio-16x9">
                            <iframe className="vimeo-video" src="https://player.vimeo.com/video/351028818" allowFullScreen />
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
                    {currentSupporterLogoItems.length ? (<h3 className="logo-header">Currently Supported By:</h3>) : '' }
                    {currentSupporterLogoItems}
                </Row>
                <Row className="logo-block">
                    <h3 className="logo-header no-margin-bottom">Consider supporting this open-source project by <Link to={'/whydonate'}>donating</Link> today.</h3>
                </Row>
            </Grid>
        );
    }
}
const Home = withRouter(HomeRaw);

class About extends React.Component {
    render() {
        const { page } =
		(this.props.match && this.props.match.params) ||
 	      	this.props.params ||
		{};
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
                        <Col sm={{ span: 10, offset: 1 }}>
                            <RawHTML html={content.pages[page]} />
                        </Col>
                    </Row>
                    <Row className="logo-block">
                        {notCurrentSupporterLogoItems}
                    </Row>
                    <Row className="logo-block">
                        {currentSupporterLogoItems.length ? (<h3 className="logo-header">Currently Supported By:</h3>) : '' }
                        {currentSupporterLogoItems}
                    </Row>
                </Grid>
            );
        } else {
            return (
                <Grid id="main-grid" className="main-grid">
                    <Row>
                        <Col sm={{ span: 10, offset: 1 }}>
                            <RawHTML html={content.pages[page]} />
                        </Col>
                    </Row>
                </Grid>
            );
        }
    }
}
/*
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
        // Updated RxJS 6+ syntax
        var q = this.urlq = new Subject();
        this.subs = q.pipe(debounceTime(500)).subscribe(this.onChange);
    },
    componentWillUnmount: function () {
        this.subs.unsubscribe();
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
    onChange: function (state) {
        if (this.props.show) {
            var d3TipDiv = document.getElementsByClassName('d3-tip-selection');
            if (d3TipDiv.length !== 0 && d3TipDiv[0].style.opacity !== '0') {
                d3TipDiv[0].style.opacity = '0';
                d3TipDiv[0].style.pointerEvents = 'none';
            }
            if (!this.state.showModal && !this.state.restoringDefaults) {
                this.transitionTo('/variants', {}, urlFromDatabase(state));
            } else if (this.state.restoringDefaults) {
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
					url={backend.url}
					onChange={s => this.urlq.next(s)}
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
					url={backend.url}
					onChange={s => this.urlq.next(s)}
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
					{this.props.mode === 'research_mode' && <Button variant="secondary" onClick={this.toggleMode}>
						Show Summary Data Only
					</Button>}
					{this.props.mode === 'default' &&
					<Button variant="secondary" onClick={() =>this.setState({showModal: true})}>
						Show Detail View
					</Button>}
					{
					    this.props.mode === 'default' && this.state.showModal &&
                        <Modal show={true} onHide={() => this.setState({ showModal: false })}>
                            <Modal.Body>
                                <RawHTML html={content.pages.researchWarning}/>
                            </Modal.Body>
                            <Modal.Footer>
                                <Button variant="primary" onClick={() => {this.toggleMode();}}>Yes</Button>
                                <Button variant="secondary" onClick={() => this.setState({ showModal: false })}>No</Button>
                            </Modal.Footer>
                        </Modal>
					}
				</Col>
			</Row>);
    }
});

// get display name for a given key from VariantTable.js column specification
function getDisplayName(key) {
    const researchMode = (localStorage.getItem("research-mode") === 'true');
    let displayName;
    if (!researchMode) {
        displayName = columns.find(e => e.prop === key);
        displayName = displayName && displayName.title;
    }
    if (displayName === undefined) {
        displayName = researchModeColumns.find(e => e.prop === key);
        displayName = displayName && displayName.title;
    }
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
    render: function () {
        const children = this.props.children;
        return (
            <div className="isogrid">
            {children}
            </div>
        );
    },
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

// Helpers for new Card/Collapse behavior
function isOpenFromStorage(key) {
    // legacy semantics: localStorage "true" means collapsed
    return localStorage.getItem(key) !== "true";
}

var VariantDetail = React.createClass({
    mixins: [Navigation, State],
    showHelp: function (event, title) {
        event.preventDefault();
        this.transitionTo(`/help#${slugify(title)}`);
    },
    getInitialState: () => ({
        hideEmptyItems: (localStorage.getItem("hide-empties") === 'true'),
        tooltips: parseTooltips(localStorage.getItem("research-mode") === 'true'),
        // track open/closed state for cards (keyed by localStorage key)
        openGroups: {}
    }),
    componentWillMount: function () {
        backend.variant(this.props.params.id).subscribe(
            resp => {
                if (resp.hasOwnProperty('redirect') && resp.redirect === true) {
                    this.transitionTo('/variants', null, {search: resp.data});
                } else {
                    this.setState({data: resp.data, error: null});
                }
            },
            () => { this.setState({error: 'Problem connecting to server'}); }
        );

        backend.variantReports(this.props.params.id).subscribe(
            resp => {
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
        if (nextProps.mode !== this.props.mode) {
            this.setState({
                tooltips: parseTooltips(nextProps.mode === 'research_mode')
            });
        }

        const { data } = nextState;
        if (data && this.props.params.id !== data[0].CA_ID) {
            const { noRedirect } = this.getQuery();
            if (noRedirect !== "true") {
                const redirectToID = data[0].CA_ID  || data[0].id ;
                if (parseInt(this.props.params.id) === data[0].id) {
                    this.replaceWith(`/variant/:id`, { id: redirectToID }, {});
                }
                else {
                    this.replaceWith(`/variant/:id`, { id: redirectToID }, { redirectedFrom: this.props.params.id });
                }
            }
        }
    },
    componentDidUpdate: function(prevProps, prevState) {
        if (prevProps.mode !== this.props.mode) {
            setTimeout(() => {
                this.relayoutGrid(true);
            }, 0);
        }
        if (!prevState.data && this.state.data) {
            const data = this.state.data;
            const variantVersionIdx = data.findIndex(x => x.id === parseInt(this.props.params.id));
            const variant = data[variantVersionIdx] || data[0];
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
    relayoutOnCollapsed: function(/* collapser *//*) {
        console.warn("Deprecated relayoutOnCollapsed; replace relayoutOnCollapsed handlers w/direct calls to relayoutGrid() in your collapsing components");
    },
    // legacy API kept for callers; now just flips storage and local openGroups
    onChangeGroupVisibility(groupTitle, event) {
        event.preventDefault();
        const key = "collapse-group_" + groupTitle;
        const willBeCollapsed = localStorage.getItem(key) !== "true";
        localStorage.setItem(key, willBeCollapsed ? "true" : "false");
        const willBeOpen = !willBeCollapsed;
        this.setState({ openGroups: { ...this.state.openGroups, [key]: willBeOpen } });
    },
    isGroupOpenLS(key) {
        // local state wins, otherwise read from storage (for initial render)
        if (this.state.openGroups.hasOwnProperty(key)) { return !!this.state.openGroups[key]; }
        return isOpenFromStorage(key);
    },
    toggleCard(key) {
        const nowOpen = !this.isGroupOpenLS(key);
        // store "true" when collapsed (legacy semantics)
        localStorage.setItem(key, nowOpen ? "false" : "true");
        this.setState({ openGroups: { ...this.state.openGroups, [key]: nowOpen } });
    },
    determineDiffRowColor: function(highlightRow) {
        return highlightRow ? 'table-danger' : '';
    },
    getPathogenicity: function(version, isReport) {
        if (isReport) {
            if (version.Source === "ClinVar") {
                return util.getFormattedFieldByProp("Clinical_Significance_ClinVar", version);
            } else {
                return util.getFormattedFieldByProp("Classification_LOVD", version);
            }
        } else {
            return util.getFormattedFieldByProp("Pathogenicity_expert", version);
        }
    },
    generateDiffRows: function(cols, data, isReports) {
        var diffRows = [];
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
                                    <span className='badge bg-success'><span className='glyphicon glyphicon-star'></span> New</span>
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
            const k = `submitter-group-${sourceName}-${submitter}`;
            return {
                [k]: !(!pstate.hasOwnProperty(k) || pstate[k])
            };
        });
    },
    // render for VariantDetail
    render: function () {
        const {data, error} = this.state;
        if (!data) {
            return <div />;
        }

        const variantVersionIdx = data.findIndex(x => x.id === parseInt(this.props.params.id));
        const variant = data[variantVersionIdx] || data[0];
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

        let groupsEmpty = 0;
        let totalRowsEmpty = 0;

        const groupTables = _.map(groups, ({ groupTitle, innerCols, reportSource, reportBinding, alleleFrequencies, inSilicoPred, innerGroups }) => {
            let rowsEmpty = 0;

            if (reportSource) {
                if (!this.state.reports || !this.state.reports[reportSource]) {
                    return null;
                }
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
                const results = {
                    Biwas: variant.Result_Biwas_ENIGMA_BRCA12_Functional_Assays,
                    Bouwman1: variant.Result_Bouwman1_ENIGMA_BRCA12_Functional_Assays,
                    Bouwman2: variant.Result_Bouwman2_ENIGMA_BRCA12_Functional_Assays,
                    Fernandes: variant.Result_Fernandes_ENIGMA_BRCA12_Functional_Assays,
                    Findlay: variant.Result_Findlay_ENIGMA_BRCA12_Functional_Assays,
                    Ikegami: variant.Result_Ikegami_ENIGMA_BRCA12_Functional_Assays,
                    Mesman: variant.Result_Mesman_ENIGMA_BRCA12_Functional_Assays,
                    Petitalot: variant.Result_Petitalot_ENIGMA_BRCA12_Functional_Assays,
                    Richardson: variant.Result_Richardson_ENIGMA_BRCA12_Functional_Assays,
                    Starita: variant.Result_Starita_ENIGMA_BRCA12_Functional_Assays
                };

                return (
                    <FunctionalAssayTile
                        key="source-reports-tile"
                        groupTitle='functional-assay-tile'
                        results={results}
                        displayTitle="Functional Assay Results"
                        onChangeGroupVisibility={this.onChangeGroupVisibility}
                        hideEmptyItems={this.state.hideEmptyItems}
                        relayoutGrid={this.relayoutGrid}
                        helpSection="functional-assay-results"
                        showHelp={this.showHelp}
                        tooltips={this.state.tooltips}
                        variant={variant}
                        innerGroups={innerGroups}
                    />
                );
            }

            if (groupTitle === "Computational Predictions") {
                return (
                    <ComputationalPredictionTile
                        key="source-reports-tile"
                        groupTitle='functional-assay-tile'
                        displayTitle="Computational Predictions"
                        onChangeGroupVisibility={this.onChangeGroupVisibility}
                        hideEmptyItems={this.state.hideEmptyItems}
                        relayoutGrid={this.relayoutGrid}
                        helpSection="computational-predictions"
                        showHelp={this.showHelp}
                        tooltips={this.state.tooltips}
                        variant={variant}
                        innerGroups={innerGroups}
                    />
                );
            }

            if (groupTitle === "ACMG Variant Evidence Codes, Provisional Assignment") {
                return (
                    <ProvisionalEvidenceTile
                        groupTitle={groupTitle}
                        onChangeGroupVisibility={this.onChangeGroupVisibility}
                        relayoutGrid={this.relayoutGrid}
                        helpSection="acmg-variant-evidence-codes-provisional-assignment"
                        showHelp={this.showHelp}
                        variant={variant}
                        innerGroups={innerGroups}
                    />
                );
            }

            // standard 2-column table tile
            const rows = _.map(innerCols, (rowDescriptor) => {
                let {prop, title, noHelpLink} = rowDescriptor;
                let rowItem;

                if (prop === "Protein_Change") {
                    title = "Abbreviated AA Change";
                }

                if (prop === "Mupit_Structure") {
                    rowItem = <MupitStructure variant={variant} prop={prop} onLoad={() => this.relayoutGrid()} />;
                    if (util.getAminoAcidCode(variant["HGVS_Protein"]) === false) {
                        rowsEmpty += 1;
                        rowItem = false;
                    }
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

            const storageKey = "collapse-group_" + groupTitle;
            const isOpen = this.isGroupOpenLS(storageKey);

            return (
                <div key={`group_collection-${groupTitle}`} className={ (allEmpty && this.state.hideEmptyItems) || (allEmpty && groupTitle === 'CRAVAT - MuPIT 3D Protein View') ? "group-empty" : "" }>
                    <Card>
                        <Card.Header
                            role="button"
                            aria-expanded={isOpen}
                            className="d-flex justify-content-between align-items-center"
                            onClick={(event) => { event.preventDefault(); this.toggleCard(storageKey); this.onChangeGroupVisibility(groupTitle, event); }}
                        >
                            <span className="title">{groupTitle}</span>
                            <span className="d-flex align-items-center">
                                <GroupHelpButton onClick={(event) => { this.showHelp(event, groupTitle); return true; }} />
                                <i className={`fa fa-chevron-${isOpen ? 'down' : 'right'} ms-2`} aria-hidden="true" />
                            </span>
                        </Card.Header>
                        <Collapse in={isOpen} onEntered={this.relayoutGrid} onExited={this.relayoutGrid}>
                            <div>
                                <Card.Body>
                                    {tileTable}
                                </Card.Body>
                            </div>
                        </Collapse>
                    </Card>
                </div>
            );
        });

        // generates variant diff rows
        const diffRows = this.generateDiffRows(cols, data, false);

        // generates report diff rows
        if (this.state.reports !== undefined) {
            let sortedSubmissions = {'ClinVar': {}, 'LOVD': {}};

            if (this.state.reports.hasOwnProperty('ClinVar')) {
                let clinvarSubmissions = this.state.reports.ClinVar;
                for (var i = 0; i < clinvarSubmissions.length; i++) {
                    if (clinvarSubmissions[i].Diff === null || clinvarSubmissions[i].Diff === undefined) {
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

            if (this.state.reports.hasOwnProperty('LOVD')) {
                let lovdSubmissions = this.state.reports.LOVD;
                for (var j = 0; j < lovdSubmissions.length; j++) {
                    if (lovdSubmissions[j].Diff === null || lovdSubmissions[j].Diff === undefined) {
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
                                    <tr className='table-active'>
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
                                    <tr className='table-active'>
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

        // splicing panel open state
        const splicingKey = "collapse-group_transcript-visualization";
        const splicingOpen = this.isGroupOpenLS(splicingKey);

        return (error ? <p>{error}</p> :
            <Grid>
                <Row>
                    {
                        (this.props.mode !== "research_mode")
                            ? (
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
                                <Col xs={4} sm={{ span: 4, offset: 4 }} md={{ span: 4, offset: 4 }} className="vcenterblock">
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
                                variant="secondary">
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
                        (variant.id !== data[0].id && noRedirectMsg !== "true") && (
                          <Col xs={12} classname="vcenterblock">
                              <div className="variant-message outdated-variant-message alert alert-danger">
                                  <h3 style={{marginTop: 0}}>There is new data available on this variant.</h3>
                                  The data below is from {util.reformatDate(variant.Data_Release.date)} (Release {variant.Data_Release.name}). <a href={`/variant/${data[0].CA_ID}`}>Click here for updated data on this variant.</a>
                              </div>
                          </Col>
                        )
                    }

                    {
                        redirectedFromVariant && (
                          <Col xs={12} classname="vcenterblock">
                              <div className="variant-message redirected-variant-msg alert alert-primary">
                                  <h3 style={{marginTop: 0}}>You are viewing the most recent data on this variant.</h3>
                                  The variant url you requested only has data up to {util.reformatDate(redirectedFromVariant.Data_Release.date)}. You have been automatically redirected to the newest data.<br />
                                  <a href={`/variant/${redirectedFrom}?noRedirect=true`}>
                                      Click here to view variant data from {util.reformatDate(redirectedFromVariant.Data_Release.date)} (Release {redirectedFromVariant.Data_Release.name}).
                                  </a>
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
                                this.props.mode === "research_mode" && (
                                    <Col key="splicing_vis"
                                        className={`variant-detail-group isogrid-item ${splicingTileSizeClassse}`}>
                                        <Card>
                                            <Card.Header
                                                role="button"
                                                aria-expanded={splicingOpen}
                                                className="d-flex justify-content-between align-items-center"
                                                onClick={(e) => { e.preventDefault(); this.toggleCard(splicingKey); this.onChangeGroupVisibility("transcript-visualization", e); }}
                                            >
                                                <span className="title">{`${variant['Gene_Symbol']} ${variant['HGVS_cDNA']} Transcript Visualization`}</span>
                                                <span className="d-flex align-items-center">
                                                    <GroupHelpButton group={"transcript-visualization"} onClick={(event) => { this.showHelp(event, "transcript-visualization"); return true; }} />
                                                    <i className={`fa fa-chevron-${splicingOpen ? 'down' : 'right'} ms-2`} aria-hidden="true" />
                                                </span>
                                            </Card.Header>
                                            <Collapse in={splicingOpen} onEntered={this.relayoutGrid} onExited={this.relayoutGrid}>
                                                <div>
                                                    <Card.Body>
                                                        <Splicing variant={variant} relayoutGrid={this.relayoutGrid} />
                                                    </Card.Body>
                                                </div>
                                            </Collapse>
                                        </Card>
                                    </Col>
                                )
                            }

                            {
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
                                <tr className='table-active'>
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
                    <Col md={{ span: 12, offset: 0 }}>
                        <DisclaimerModal buttonModal onToggleMode={this.props.toggleMode} text="Show Detail View for this Variant"/>
                    </Col>
                </Row>
            </Grid>
        );
    }
});
*/
class Application extends React.Component {
    constructor(props){
	super(props);
	this.state = {
	    mode: localStorage.getItem('research-mode') === 'true' ? 'research_mode' : 'default',
	};
    }

	onChildToggleMode = () => {
        this.toggleMode();
    };

    componentDidUpdate() {
        const localStorageMode = localStorage.getItem("research-mode") === "true" ? "research_mode" : "default";
        if (localStorageMode !== this.state.mode) {
            this.setMode();
        }
    }

    setMode = () => {
        this.setState({mode: localStorage.getItem("research-mode") === 'true' ? 'research_mode' : 'default'});
    };

    toggleMode = () => {
        if (this.state.mode === 'research_mode') {
            localStorage.setItem('research-mode', false);
            this.setState({mode: 'default'});
        } else {
            localStorage.setItem('research-mode', true);
            this.setState({mode: 'research_mode'});
        }
    };

    render() {
        const path = (this.props.location && this.props.location.pathname ? this.props.location.pathname : '/').slice(1);
        return (
            <div>
                <NavBarNew path={path} mode={this.state.mode} toggleMode={this.toggleMode}/>
                <DonationBar />
                {/* If children are rendered here via a parent Route, pass props along */}
         	{this.props.children &&
           		React.cloneElement(this.props.children, {
             			toggleMode: this.onChildToggleMode,
             			mode: this.state.mode,
           	})}
		{path.indexOf('variants') === 0 && (
		<Database
                    mode={this.state.mode}
                    toggleMode={this.onChildToggleMode}
                    show={path.indexOf('variants') === 0} /> 
		)}
                <Footer />
            </div>
        );
    }
}

const routes = (
    <Switch>
        <Route exact path='/' component={Home}/>
	<Route path='/about/:page' component={About}/>
	<Route path='/community' component={Community}/>
	<Route path='/help' component={Help}/>
	<Route path='/factsheet' component={FactSheet}/>
	<Route path='/releases' component={Releases}/>
        <Route path='/release/:id' component={Release}/>
        <Route path='/whydonate' component={WhyDonate}/>
        <Route path='/fundraisingdetails' component={FundraisingDetails}/>
	{/*
        <Route path='signup' component={Signup}/>
        <Route path='signin' component={Signin}/>
        <Route path='reset_password' component={ResetPassword}/>
        <Route path='profile' component={Profile}/>
        <Route path='confirm/:activationCode' component={ConfirmEmail}/>
        <Route path='reset/:resetToken' component={ChangePassword}/>
        // TODO: wire this to your variants page component
        // <Route path='variants' />
        <Route path='variant/:id' component={VariantDetail}/>
        <Route path='variant_literature/:id' component={LiteratureTable}/>
	*/}
    </Switch>
);

const container = document.getElementById('main');
const root = createRoot(container);
root.render(
    <BrowserRouter>
	<Application>
		{routes}
	</Application>
    </BrowserRouter>
);
