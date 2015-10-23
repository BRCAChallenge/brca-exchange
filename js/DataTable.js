/*global module: false, require: false, URL: false, Blob: false */
'use strict';

var React = require('react');
var {Table, Pagination, DataMixin} = require('react-data-components-bd2k');
var {Button, Row, Col} = require('react-bootstrap');
var VariantSearch = require('./VariantSearch');
var SelectField = require('./SelectField');
var ColumnCheckbox = require('./ColumnCheckbox');
var _ = require('underscore');
var cx = require('classnames');
var hgvs = require('./hgvs');
var PureRenderMixin = require('./PureRenderMixin'); // deep-equals version of PRM

var filterDisplay = v => v == null ? 'Any' : v;
var filterAny = v => v === 'Any' ? null : v;
var addAny = opts => ['Any', ...opts];

var pluralize = (n, s) => n === 1 ? s : s + 's';

var merge = (...args) => _.extend({}, ...args);

var DataTable = React.createClass({
	mixins: [DataMixin, PureRenderMixin],
	// This differs from the default DataMixin onFilter handler in that multiple filters
	// can be set at once. This is important in order to avoid redundant sorting when
	// setting multiple filters. Also, the onFilter code changes the filter state in-place,
	// which violates the react API contract wherein state is treated as immutable.
	setFilters: function (obj) {
		var {filterValues, sortBy} = this.state,
			{initialData, filter, filters, sort} = this.props,
			newFilterValues = merge(filterValues, obj),
			data = sort(sortBy, filter(filters, newFilterValues, initialData));

		this.setState({
		  data: data,
		  filterValues: newFilterValues,
		  currentPage: 0
		});
	},
	componentWillReceiveProps: function(newProps) {
		var {search} = newProps;
		if (search !== this.props.search) {
			this.setState({search: search});
			this.setFilters(hgvs.filters(search));
		}
	},
	createDownload: function (ev) {
		// XXX This is a bit horrible. In order to build the tsv lazily (on
		// button click, instead of on every table update), we catch the
		// mousedown event and modify the href on the anchor element, behind
		// the back of react. I don't believe this will cause any problems, but
		// it's something to be aware of if react starts doing something
		// strange.  Also needs to be tested cross-browser. We should not offer
		// download on browsers that don't allow client-driven download.
		var data = this.state.data,
			keys = _.keys(data[0]),
			tsvRows = _.map(data, obj => _.map(keys, k => obj[k]).join('\t')).join('\n'), // use os-specific line endings?
			tsv = keys.join('\t') + '\n' + tsvRows;
		ev.target.href = URL.createObjectURL(new Blob([tsv], { type: 'text/tsv' }));
	},
	getInitialState: function () {
		return {filtersOpen: false, search: this.props.search, renderColumns: this.selectColumns()};
	},
	toggleFilters: function () {
		this.setState({filtersOpen: !this.state.filtersOpen});
	},
    toggleColumns: function (title) {
        this.props.columnSelection[title].selectVal = !this.props.columnSelection[title].selectVal;
        this.setState({renderColumns: this.selectColumns()});
    },
    selectColumns () {
        var columnObject = this.props.origionalColumns;
        var newColObject = [];
        for (var i = 0; i < columnObject.length; i++) {
            var title = columnObject[i].prop;
            if (this.props.columnSelection[title].selectVal == true) {
                newColObject.push(columnObject[i]);
            }
        }
        return newColObject;
    },
	render: function () {
		var {filtersOpen, filterValues, search} = this.state,
			{origionalColumns, columnSelection, filterColumns, suggestions, className} = this.props,
			page = this.buildPage(),
			filterFormEls = _.map(filterColumns, ({name, prop, values}) =>
				<SelectField onChange={v => this.setFilters({[prop]: filterAny(v)})}
					key={prop} label={`${name} is: `} value={filterDisplay(filterValues[prop])} options={addAny(values)}/>),
			filterFormCols = _.map(origionalColumns, ({title, prop}) =>
				<ColumnCheckbox onChange={v => this.toggleColumns(prop)} key={prop} label={prop} title={title} initialCheck={columnSelection}/>);

		return (
			<div className={this.props.className}>
				<Row style={{marginBottom: '2px'}}>
					<Col sm={12}>
						<Button bsSize='xsmall' onClick={this.toggleFilters}>{(filtersOpen ? 'Hide' : 'Show' ) + ' Filters'}</Button>
						{filtersOpen && <div className='form-inline'>{filterFormEls}</div>}
                        {filtersOpen && <div className='form-inline'>
                                            <label className='control-label' style={{marginRight: '1em'}}>
                                                Column Selection
                                                {filterFormCols}
                                            </label>
                                        </div>}
					</Col>
				</Row>
				<Row style={{marginBottom: '2px'}}>
					<Col sm={6}>
						<div className='form-inline'>
							<div className='form-group'>
								<label className='control-label'
										style={{marginRight: '1em'}}>
									{this.state.data.length} matching {pluralize(this.state.data.length, 'variant')}
								</label>
								<Button download="variants.tsv" href="#" onMouseDown={this.createDownload}>Download</Button>
							</div>
						</div>
					</Col>
					<Col sm={3} smOffset={3}>
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
				<Row style={{marginBottom: '2px'}}>
					<Col sm={5}>
						<VariantSearch
							id='variants-search'
							suggestions={suggestions}
							value={search}
							/* XXX should debounce this so we don't sort so often */
							onChange={v => {
								this.setState({search: v});
								this.setFilters(hgvs.filters(v));
								this.props.onChange(v);
							}}
						/>
					</Col>
					<Col sm={6} smOffset={1}>
						<Pagination
							className="pagination pull-right-sm"
							currentPage={page.currentPage}
							totalPages={page.totalPages}
							onChangePage={this.onChangePage} />
					</Col>
				</Row>
				<Row>
					<Col className="table-responsive" sm={12}>
						<Table
							className={cx(className, "table table-hover table-bordered table-condensed")}
							dataArray={page.data}
							columns={this.state.renderColumns}
							keys={this.props.keys}
							buildRowOptions={this.props.buildRowOptions}
							buildHeader={this.props.buildHeader}
							sortBy={this.state.sortBy}
							onSort={this.onSort} />
					</Col>
				</Row>
			</div>
		);
	}
});

module.exports = DataTable;
