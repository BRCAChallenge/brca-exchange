/*global module: false, require: false, URL: false, Blob: false */
'use strict';

import React from 'react';
// RxJS imports
import { Subject } from 'rxjs';
import { map, debounceTime, switchMap } from 'rxjs/operators';

// TODO: Uncomment when react-data-components-brcaex is updated/replaced
// import { Table, Pagination } from 'react-data-components-brcaex';

import { Button, Row, Col, Table as BootstrapTable } from 'react-bootstrap';
import VariantSearch from './VariantSearch';
import SelectField from './SelectField';
import _ from 'underscore';
import cx from 'classnames';
import hgvs from './hgvs';
import PureRenderMixin from './PureRenderMixin'; // deep-equals version of PRM
import util from './util';

var filterDisplay = v => v == null ? 'Any' : v;
var filterAny = v => v === 'Any' ? null : v;
var addAny = opts => ['Any', ...opts];

var pluralize = (n, s) => n === 1 ? s : s + 's';

var merge = (...args) => _.extend({}, ...args);


function setPages({data, count, deletedCount, synonyms, releaseName}, pageLength) { //eslint-disable-line camelcase
    return {
        data,
        count,
        deletedCount,
        synonyms,
        releaseName,
        totalPages: Math.ceil(count / pageLength)
    };
}

// TODO: Uncomment when react-data-components-brcaex is updated/replaced
// This is the original FastTable that wraps the Table from react-data-components-brcaex
/*
// Wrap Table with a version having PureRenderMixin
var FastTable = React.createClass({
    mixins: [PureRenderMixin],
    normalizeGenomicCoordinate: function(field, hgvs, variantData) {
        if (!util.isEmptyField(hgvs)) {
            variantData[field] = hgvs;
        } else {
            if (variantData[field].length > 35) {
                variantData[field] = variantData[field].substring(0, 35) + "...";
            }
        }
        return variantData;
    },
    determineGenomicCoordinates: function(variantData) {
        let field = "Genomic_Coordinate_hg37";
        let hgvs = variantData.Genomic_HGVS_37;
        variantData = this.normalizeGenomicCoordinate(field, hgvs, variantData);

        field = "Genomic_Coordinate_hg38";
        hgvs = variantData.Genomic_HGVS_38;
        variantData = this.normalizeGenomicCoordinate(field, hgvs, variantData);
        return variantData;
    },
    render: function () {
        var {dataArray, ...props} = this.props;
        if (dataArray.length > 0) {
            dataArray = _.map(dataArray, function(variantData) {
                return this.determineGenomicCoordinates(variantData);
            }, this);
        }
        return <Table {...props} dataArray={dataArray}/>;
    }
});
*/

// TEMPORARY: Simple table component until react-data-components-brcaex is available
var FastTable = React.createClass({
    mixins: [PureRenderMixin],
    normalizeGenomicCoordinate: function(field, hgvsValue, variantData) {
        if (!util.isEmptyField(hgvsValue)) {
            variantData[field] = hgvsValue;
        } else {
            if (variantData[field] && variantData[field].length > 35) {
                variantData[field] = variantData[field].substring(0, 35) + "...";
            }
        }
        return variantData;
    },
    determineGenomicCoordinates: function(variantData) {
        let field = "Genomic_Coordinate_hg37";
        let hgvsValue = variantData.Genomic_HGVS_37;
        variantData = this.normalizeGenomicCoordinate(field, hgvsValue, variantData);

        field = "Genomic_Coordinate_hg38";
        hgvsValue = variantData.Genomic_HGVS_38;
        variantData = this.normalizeGenomicCoordinate(field, hgvsValue, variantData);
        return variantData;
    },
    onSort: function(column) {
        if (this.props.onSort) {
            const currentSort = this.props.sortBy;
            const newOrder = (currentSort && currentSort.prop === column.prop && currentSort.order === 'asc') ? 'desc' : 'asc';
            this.props.onSort({prop: column.prop, order: newOrder});
        }
    },
    render: function () {
        var {dataArray, columns, className, buildRowOptions, onRowClick, buildHeader, sortBy} = this.props;
        
        if (dataArray.length > 0) {
            dataArray = _.map(dataArray, function(variantData) {
                return this.determineGenomicCoordinates(variantData);
            }, this);
        }

        return (
            <BootstrapTable className={className} striped bordered hover>
                <thead>
                    <tr>
                        {columns.map((column, idx) => {
                            const isSorted = sortBy && sortBy.prop === column.prop;
                            const sortIcon = isSorted ? (sortBy.order === 'asc' ? ' ▲' : ' ▼') : '';
                            
                            return (
                                <th 
                                    key={idx}
                                    onClick={() => this.onSort(column)}
                                    style={{cursor: 'pointer'}}
                                >
                                    {buildHeader ? buildHeader(column) : column.title}
                                    {sortIcon}
                                </th>
                            );
                        })}
                    </tr>
                </thead>
                <tbody>
                    {dataArray.length === 0 ? (
                        <tr>
                            <td colSpan={columns.length} style={{textAlign: 'center', padding: '20px'}}>
                                No variants found
                            </td>
                        </tr>
                    ) : (
                        dataArray.map((row, rowIdx) => {
                            const rowOptions = buildRowOptions ? buildRowOptions(row) : {};
                            return (
                                <tr 
                                    key={rowIdx} 
                                    onClick={(e) => onRowClick && onRowClick(row, e)}
                                    style={{cursor: onRowClick ? 'pointer' : 'default'}}
                                    {...rowOptions}
                                >
                                    {columns.map((column, colIdx) => (
                                        <td key={colIdx}>
                                            {column.render ? column.render(row[column.prop], row) : row[column.prop]}
                                        </td>
                                    ))}
                                </tr>
                            );
                        })
                    )}
                </tbody>
            </BootstrapTable>
        );
    }
});

// Merge new state (e.g. initialState) with existing state,
// deep-merging fields that are objects.
function mergeState(state, newState) {
    var {columnSelection, sourceSelection, filterValues, ...otherProps} = newState,
        cs = {...state.columnSelection, ...columnSelection},
        ss = {...state.sourceSelection, ...sourceSelection},
        fv = {...state.filterValues, ...filterValues};
    return {...state, columnSelection: cs, sourceSelection: ss, filterValues: fv, ...otherProps};
}

var DataTable = React.createClass({
    shouldComponentUpdate: function (nextProps, nextState) {
        return (
            this.state.filtersOpen !== nextState.filtersOpen ||
            this.state.columnSelectorsOpen !== nextState.columnSelectorsOpen ||
            this.state.page !== nextState.page ||
            this.state.count !== nextState.count ||
            this.state.pageLength !== nextState.pageLength ||
            this.props.search !== nextProps.search ||
            !_.isEqual(this.state.sortBy, nextState.sortBy) ||
            !_.isEqual(this.props.sourceSelection, nextProps.sourceSelection) ||
            !_.isEqual(this.props.columnSelection, nextProps.columnSelection) ||
            !_.isEqual(_.sortBy(this.props.hide), _.sortBy(nextProps.hide)) ||
            !_.isEqual(this.state.filterValues, nextState.filterValues) ||
            !_.isEqual(this.state.filterColumns, nextState.filterColumns) ||
            !_.isEqual(_.map(this.state.data, r => r.id), _.map(nextState.data, r=> r.id))
        );
    },
    componentWillMount: function () {
        var q = this.fetchq = new Subject();
        this.subs = q.pipe(
            map(this.props.fetch),
            debounceTime(100),
            switchMap(obs => obs)
        ).subscribe(
            resp => this.setState(setPages(resp, this.state.pageLength)), // set data, count, totalPages
            () => this.setState({error: 'Problem connecting to server'}));
    },
    componentWillUnmount: function () {
        window.removeEventListener('resize', this.handleResize);
        this.subs.unsubscribe();
    },
    componentDidMount: function () {
        window.addEventListener('resize', this.handleResize);
        this.fetch(this.state);
    },
    getInitialState: function () {
        let filterValues = JSON.parse(localStorage.getItem('filterValues'));
        if (filterValues === null || filterValues === undefined) {
            filterValues = {};
        }
        return mergeState({
            data: [],
            filtersOpen: false,
            filterValues: filterValues,
            columnSelectorsOpen: false,
            search: '',
            mode: this.props.mode,
            columnSelection: this.props.columnSelection,
            sourceSelection: this.props.sourceSelection,
            pageLength: 20,
            page: 0,
            totalPages: 0,
            windowWidth: window.innerWidth
        }, this.props.initialState);
    },
    componentWillReceiveProps: function(newProps) {
        var newState = mergeState(this.state, newProps.initialState);
        newState.sourceSelection = newProps.sourceSelection;
        newState.columnSelection = newProps.columnSelection;
        this.setState(newState);
        if (newProps.changeInProgress === true) {
            // if a change is made midway through the component
            // tree, changes must propagate back up to the top before
            // a new request is made
            this.propagateChanges(newState);
        } else {
            this.fetch(newState);
        }
    },
    handleResize: function() {
        this.setState({windowWidth: window.innerWidth});
    },
    setFilters: function (obj) {
        let {filterValues} = this.state;
        let newFilterValues = merge(filterValues, obj);

        localStorage.setItem('filterValues', JSON.stringify(newFilterValues));

        this.propagateChanges({
          filterValues: newFilterValues,
          page: 0
        });
    },
    createDownload: function () {
        var {release, changeTypes, search, sortBy, filterValues, columnSelection, sourceSelection} = this.state;
        return this.props.url(merge({
            format: 'tsv',
            release,
            changeTypes,
            pageLength: null,
            page: null,
            sortBy,
            search,
            searchColumn: _.keys(_.pick(columnSelection, v => v)),
            include: _.keys(_.pick(sourceSelection, v => v === 1)),
            exclude: _.keys(_.pick(sourceSelection, v => v === -1)),
            filterValues}, hgvs.filters(search, filterValues)));
    },
    fetch: function (state) {
        var {pageLength, search, page, sortBy,
            filterValues, columnSelection, sourceSelection,
            release, changeTypes, showDeleted, mode} = state;
        this.fetchq.next(merge({
            release,
            changeTypes,
            showDeleted,
            pageLength,
            page,
            sortBy,
            search,
            mode,
            searchColumn: _.keys(_.pick(columnSelection, v => v)),
            include: _.keys(_.pick(sourceSelection, v => v === 1)),
            exclude: _.keys(_.pick(sourceSelection, v => v === -1)),
            filterValues}, hgvs.filters(search, filterValues)));
    },
    propagateChanges: function (opts) {
        // sends changes up to the top of the tree so the URL
        // can be updated. the update also triggers a new request
        // to the server in this.componentWillReceiveProps
        var newState = mergeState(this.state, opts);
        this.setState(newState);
        this.props.onChange(newState);
    },
    toggleFilters: function () {
        this.setState({filtersOpen: !this.state.filtersOpen});
    },
    toggleColumnSelectors: function() {
        this.setState({columnSelectorsOpen: !this.state.columnSelectorsOpen});
    },
    showDeleted: function () {
        this.propagateChanges({showDeleted: true});
    },
    onChangePage: function (pageNumber) {
        this.propagateChanges({page: pageNumber});
    },
    onSort: function(sortBy) {
        this.propagateChanges({sortBy});
    },
    onPageLengthChange: function(txt) {
        var length = parseInt(txt),
            {page, pageLength} = this.state,
            newPage = Math.floor((page * pageLength) / length);

        this.propagateChanges({page: newPage, pageLength: length});
    },
    restoreDefaults: function() {
        // Clears local storage, resets filters/columns/sources/releases/changetypes,
        // and resets url to /variants (parent method uses a flag to ensure query params remain empty).
        delete localStorage.columnSelection;
        delete localStorage.filterValues;
        delete localStorage.sourceSelection;
        let that = this;
        let defaults = {filterValues: {},
                        release: undefined,
                        changeTypes: undefined,
                        showDeleted: undefined
                       };
        if (this.state.mode === 'default') {
            this.props.expertVariantTableRestoreDefaults(function() {
                that.setState(defaults, function() {that.propagateChanges(defaults);});
            });
        } else {
            this.props.researchVariantTableRestoreDefaults(function() {
                that.setState(defaults, function() {that.propagateChanges(defaults);});
            });
        }
    },
    render: function () {
        var {release, changeTypes, filterValues, filtersOpen, columnSelectorsOpen, search, data, columnSelection,
            page, totalPages, count, synonyms, error} = this.state;
        var {columns, filterColumns, className, columnSelectors, filters, downloadButton, mode} = this.props;
        var renderColumns = _.filter(columns, c => columnSelection[c.prop]);
        var filterFormEls = _.map(filterColumns, ({name, prop, values}) =>
            <SelectField onChange={v => this.setFilters({[prop]: filterAny(v)})}
                         key={prop} label={`${name} is: `} value={filterDisplay(filterValues[prop])}
                         options={addAny(values)}/>);
        // assumes added / changed are lumped together
        var changeString;
        if (changeTypes) {
            if (changeTypes.includes('new')) {
                changeString = "added";
            } else if (changeTypes.includes('added_information')) {
                changeString = "with new or changed information";
            } else if (changeTypes.includes('added_classification')) {
                changeString = "with new or changed classification";
            } else if (changeTypes.includes('deleted')) {
                changeString = "deleted";
            }
        }
        let releaseName = this.state.releaseName;
        var deletedCount = this.state.deletedCount;
        var deletedVariantsNote = '';
        if (deletedCount) {
            let pl = deletedCount !== 1;
            deletedVariantsNote = (<p>
                There {pl ? 'are' : 'is'} {deletedCount} deleted variant{pl ? 's' : ''} that match{pl ? '' : 'es'} your search.
                Click <a href="#" onClick={this.showDeleted}>here</a> to view {pl ? 'these' : 'this'} deleted variant{pl ? 's' : ''}.
            </p>);
        }
        return (error ? <p>{error}</p> :
            <div className={this.props.className}>
            <div id="filters" className="container-fluid">
                <Row id="show-hide" className="btm-buffer">
                    <Col sm={12}>
                        <Button className="btn-default rgt-buffer"
                                onClick={this.toggleFilters}>{(filtersOpen ? 'Hide' : 'Show' ) + ' Filters'}
                        </Button>
                        {mode === "research_mode" && <Button className="btn-default rgt-buffer"
                                onClick={this.toggleColumnSelectors}>{(columnSelectorsOpen ? 'Hide' : 'Show' ) + ' Column Selectors'}
                        </Button>}
                        <Button className="btn-default rgt-buffer"
                                onClick={this.restoreDefaults}>Restore Defaults
                        </Button>
                    </Col>
                </Row>
                <Row id="filters">
                    <Col sm={12}>
                        {filtersOpen && <div className='form-inline'>{filterFormEls}{filters}</div>}
                        {columnSelectorsOpen && mode === "research_mode" && <div className='form-inline'>
                            {columnSelectors}
                        </div>}
                    </Col>
                </Row>
                <Row id="download" className="btm-buffer">
                    <Col sm={8} lg={10}>
                        <div className='form-inline'>
                            <div className='form-group'>
                                <label className='alert-danger matched-variant-count'>
                                {
                                    `${count} matching ${pluralize(count, 'variant')}
                                    ${changeString ? changeString : ''}
                                    ${release ? 'in release ' + releaseName : ''}
                                    ${synonyms ? 'of which ' + synonyms + ' matched on synonyms' : ''}`
                                    .replace(/[\n\t]/g, '')
                                    // using string interpolation prevents react from creating one span tag per string.
                                    // the replace() removes the extra whitespace used to make the code legible
                                }
                                </label>
                                {downloadButton(this.createDownload)}
                            </div>
                        </div>
                        { count === 0 && deletedCount !== 0 &&
                          <div>{deletedVariantsNote}</div> }
                    </Col>
                    <Col sm={4} lg={2}>
                        <div className='form-inline pull-right-sm'>
                            <SelectField
                                label="Page size:"
                                value={this.state.pageLength}
                                options={this.props.pageLengthOptions}
                                onChange={this.onPageLengthChange}
                            />
                        </div>
                    </Col>
                </Row>
                <Row id='variant-search-row' className="btm-buffer">
                    <Col sm={6}>
                        <VariantSearch
                            id='variants-search'
                            value={search}
                            release={release}
                            onChange={v => {
                                // reset the page number to zero on new searches
                                this.propagateChanges({search: v, page: 0});
                            }}
                        />
                    </Col>
                    <Col sm={6}>
                        {/* TODO: Uncomment when react-data-components-brcaex is updated/replaced */}
                        {/*
                        <Pagination
                            className="pagination pull-right-sm"
                            currentPage={page}
                            totalPages={totalPages}
                            onChangePage={this.onChangePage} />
                        */}
                        
                        {/* TEMPORARY: Basic pagination buttons until Pagination component is available */}
                        <div className="pagination pull-right-sm" style={{display: 'inline-block', float: 'right'}}>
                            <Button 
                                variant="secondary"
                                size="sm"
                                disabled={page === 0} 
                                onClick={() => this.onChangePage(page - 1)}
                                style={{marginRight: '5px'}}
                            >
                                Previous
                            </Button>
                            <span style={{padding: '0 10px', display: 'inline-block', lineHeight: '31px'}}>
                                Page {page + 1} of {totalPages}
                            </span>
                            <Button 
                                variant="secondary"
                                size="sm"
                                disabled={page >= totalPages - 1} 
                                onClick={() => this.onChangePage(page + 1)}
                            >
                                Next
                            </Button>
                        </div>
                    </Col>
                </Row>
                <Row>
                    <Col id="data-table-container" sm={12}>
                        <div className="table-responsive">
                            <FastTable
                                className={cx(className, "table table-hover table-bordered table-grayheader")}
                                dataArray={data}
                                columns={renderColumns}
                                keys={this.props.keys}
                                buildRowOptions={this.props.buildRowOptions}
                                onRowClick={this.props.onRowClick}
                                buildHeader={this.props.buildHeader}
                                sortBy={this.state.sortBy}
                                onSort={this.onSort} />
                        </div>
                    </Col>
                </Row>
                <Row>
                    <Col className="text-right" sm={12}>
                        { count !== 0 && deletedVariantsNote }
                    </Col>
                </Row>
                </div>
            </div>
        );
    }
});

export default DataTable;
