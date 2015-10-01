/*global module: false, require: false, URL: false, Blob: false */
'use strict';

var React = require('react');
var {Table, Pagination, DataMixin} = require('react-data-components-bd2k');
var {Button, Row, Col} = require('react-bootstrap');
var VariantSearch = require('./VariantSearch');
var SelectField = require('./SelectField');
var _ = require('underscore');
var cx = require('classnames');

var filterDisplay = v => v == null ? 'Any' : v;
var filterAny = v => v === 'Any' ? null : v;
var addAny = opts => ['Any', ...opts];

var pluralize = (n, s) => n === 1 ? s : s + 's';

var DataTable = React.createClass({
	mixins: [DataMixin],
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
		return {filtersOpen: false};
	},
	toggleFilters: function () {
		this.setState({filtersOpen: !this.state.filtersOpen});
	},
	render: function () {
		var {filtersOpen, filterValues} = this.state,
			{columns, filterColumns, suggestions, className} = this.props,
			page = this.buildPage(),
			filterFormEls = _.map(filterColumns, ({name, prop, values}) =>
				<SelectField onChange={v => this.onFilter(prop, filterAny(v))}
					key={prop} label={`${name} is: `} value={filterDisplay(filterValues[prop])} options={addAny(values)}/>);

		return (
			<div className={this.props.className}>
				<Row style={{marginBottom: '2px'}}>
					<Col sm={12}>
						<Button bsSize='xsmall' onClick={this.toggleFilters}>{(filtersOpen ? 'Hide' : 'Show' ) + ' Filters'}</Button>
						{filtersOpen && <div className='form-inline'>{filterFormEls}</div>}
					</Col>
				</Row>
				<Row style={{marginBottom: '2px'}}>
					<Col sm={5}>
						<VariantSearch
							id='variants-search'
							suggestions={suggestions}
							value={this.state.filterValues.visibleSearch}
							onChange={this.onFilter.bind(this, 'visibleSearch')}
						/>
					</Col>
					<Col sm={3} smOffset={4}>
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
					<Col sm={6}>
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
							columns={columns}
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
