/*global module: false, require: false, URL: false, Blob: false */
'use strict';

var React = require('react');
var Rx = require('rx');
require('rx/dist/rx.time');
var {Table, Pagination} = require('react-data-components-bd2k');
var {Button, Row, Col} = require('react-bootstrap');
var VariantSearch = require('./VariantSearch');
var SelectField = require('./SelectField');
var _ = require('underscore');
var cx = require('classnames');
var hgvs = require('./hgvs');
var PureRenderMixin = require('./PureRenderMixin'); // deep-equals version of PRM

var filterDisplay = v => v == null ? 'Any' : v;
var filterAny = v => v === 'Any' ? null : v;
var addAny = opts => ['Any', ...opts];

var pluralize = (n, s) => n === 1 ? s : s + 's';

var merge = (...args) => _.extend({}, ...args);

var Lollipop = require('./d3Lollipop');

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

// Wrap Table with a version having PureRenderMixin
var FastTable = React.createClass({
    mixins: [PureRenderMixin],
    truncateGenomicCoordinates: function(variantData) {
        const genomicCoordinateFields = ["Genomic_Coordinate_hg36", "Genomic_Coordinate_hg37", "Genomic_Coordinate_hg38"];
        genomicCoordinateFields.forEach(function(field) {
            if (variantData[field].length > 35) {
                variantData[field] = variantData[field].substring(0, 35) + "...";
            } else {
                return variantData[field];
            }
        });
        return variantData;
    },
    render: function () {
        var {dataArray, ...props} = this.props;
        if (dataArray.length > 0) {
            dataArray = _.map(dataArray, function(variantData) {
                return this.truncateGenomicCoordinates(variantData);
            }, this);
        }
        return <Table {...props} dataArray={dataArray}/>;
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
            this.state.lollipopOpen !== nextState.lollipopOpen ||
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
        var q = this.fetchq = new Rx.Subject();
        this.subs = q.map(this.props.fetch).debounce(100).switchLatest().subscribe(
            resp => this.setState(setPages(resp, this.state.pageLength)), // set data, count, totalPages
            () => this.setState({error: 'Problem connecting to server'}));
    },
    componentWillUnmount: function () {
        window.removeEventListener('resize', this.handleResize);
        this.subs.dispose();
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
            lollipopOpen: false,
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
        this.setStateFetch(newState);
    },
    handleResize: function() {
        this.setState({windowWidth: window.innerWidth});
    },
    setFilters: function (obj) {
        let {filterValues} = this.state;
        let newFilterValues = merge(filterValues, obj);

        localStorage.setItem('filterValues', JSON.stringify(newFilterValues));

        this.setStateFetch({
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
    lollipopOpts: function () {
        var {search, filterValues, sourceSelection} = this.state;
        return merge({
            search,
            include: _.keys(_.pick(sourceSelection, v => v === 1)),
            exclude: _.keys(_.pick(sourceSelection, v => v === -1)),
            filterValues
        }, hgvs.filters(search, filterValues));
    },
    fetch: function (state) {
        var {pageLength, search, page, sortBy,
            filterValues, columnSelection, sourceSelection,
            release, changeTypes, showDeleted, mode} = state;
        this.fetchq.onNext(merge({
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
    // helper function that sets state, fetches new data,
    // and updates url.
    setStateFetch: function (opts) {
        var newState = mergeState(this.state, opts);
        this.setState(newState);
        this.fetch(newState);
        this.props.onChange(newState);
    },
    toggleLollipop: function () {
        this.setState({lollipopOpen: !this.state.lollipopOpen});
    },
    toggleFilters: function () {
        this.setState({filtersOpen: !this.state.filtersOpen});
    },
    toggleColumnSelectors: function() {
        this.setState({columnSelectorsOpen: !this.state.columnSelectorsOpen});
    },
    showDeleted: function () {
        this.setStateFetch({showDeleted: true});
    },
    onChangePage: function (pageNumber) {
        this.setStateFetch({page: pageNumber});
    },
    onSort: function(sortBy) {
        this.setStateFetch({sortBy});
    },
    onPageLengthChange: function(txt) {
        var length = parseInt(txt),
            {page, pageLength} = this.state,
            newPage = Math.floor((page * pageLength) / length);

        this.setStateFetch({page: newPage, pageLength: length});
    },
    restoreDefaults: function() {
        // Clears local storage, resets filters/columns/sources/releases/changetypes,
        // and resets url to /variants (parent method uses a flag to ensure query params remain empty).
        delete localStorage.columnSelection;
        delete localStorage.filterValues;
        delete localStorage.sourceSelection;
        let that = this;
        if (this.state.mode === 'default') {
            this.props.expertVariantTableRestoreDefaults(function() {
                that.setState({filterValues: {},
                               release: undefined,
                               changeTypes: undefined,
                               showDeleted: undefined
                             });
            });
        } else {
            this.props.researchVariantTableRestoreDefaults(function() {
                that.setState({filterValues: {},
                               release: undefined,
                               changeTypes: undefined,
                               showDeleted: undefined
                             });
            });
        }
    },
    render: function () {
        var {release, changeTypes, filterValues, filtersOpen, columnSelectorsOpen, lollipopOpen, search, data, columnSelection,
            page, totalPages, count, synonyms, error} = this.state;
        var {columns, filterColumns, className, columnSelectors, filters, downloadButton, lollipopButton, mode} = this.props;
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
                        {lollipopButton(this.toggleLollipop, lollipopOpen)}
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
                <Row id="lollipop-chart">
                    <Col sm={12}>
                        {lollipopOpen && this.state.windowWidth > 991 &&
                        <Lollipop fetch={this.props.fetchLollipop} opts={this.lollipopOpts()} onHeaderClick={this.props.onHeaderClick} onRowClick={this.props.onRowClick}/> }
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
                                this.setStateFetch({search: v, page: 0});
                            }}
                        />
                    </Col>
                    <Col sm={6}>
                        <Pagination
                            className="pagination pull-right-sm"
                            currentPage={page}
                            totalPages={totalPages}
                            onChangePage={this.onChangePage} />
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

module.exports = DataTable;
