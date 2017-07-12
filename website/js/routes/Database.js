import React from 'react';
import { Navigation, State } from 'react-router';
import Rx from 'rx';
import _ from 'underscore';
import { Grid, Col, Row, Button, Modal } from 'react-bootstrap';
import backend from '../backend';
import { VariantTable, ResearchVariantTable } from '../partials/VariantTable';
import databaseKey from '../../databaseKey';
import content from '../data/content';
import RawHTML from '../helpers/RawHTML';
import slugify from '../helpers/slugify';

var variantPathJoin = row => _.map(databaseKey, k => encodeURIComponent(row[k])).join('@@');

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
                    {this.props.mode === 'default' && <img id='enigma-logo' src={require('../img/enigma_logo.jpeg')} />}
                    <RawHTML ref='content' html={message}/>
                    {this.props.mode === 'research_mode' && <Button className="btn-default" onClick={this.toggleMode}>
                        Show Expert Reviewed Data Only
                    </Button>}
                    {this.props.mode === 'default' &&
                    <Button className="btn-default" onClick={() =>this.setState({showModal: true})}>
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

export default Database;
