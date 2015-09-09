/*global module: false, require: false, console: false */
'use strict';

var React = require('react');
var PureRenderMixin = require('./PureRenderMixin');
var DataTable = require('react-data-components-bd2k').DataTable;
require('react-data-components-bd2k/css/table-twbs.css');

function buildHeader(onClick, title) {
	return (
		<span>
			{title}
			<span onClick={ev => {ev.stopPropagation(); onClick(title); }}
				className='help glyphicon glyphicon-question-sign superscript'/>
		</span>
	);
}

var columns = [
	{title: 'Gene', prop: 'Gene symbol'},
	{title: '  HGVS  ', prop: 'HGVS'},
	{title: 'Pathogenicity', prop: 'Clinical significance'},
	{title: 'Allele origin', prop: 'Allele origin'},
	{title: 'CVA', prop: 'ClinVarAccession'}
];

var VariantTable = React.createClass({
	mixins: [PureRenderMixin],
	getData: function () {
		return this.refs.table.state.data;
	},
	render: function () {
		var {data, onHeaderClick, onRowClick, ...opts} = this.props;
		return (
			<DataTable
				ref='table'
				{...opts}
				buildRowOptions={r => ({title: 'click for details', onClick: () => onRowClick(r.id)})}
				buildHeader={title => buildHeader(onHeaderClick, title)}
				columns={columns}
				initialData={data}
				initialPageLength={10}
                initialSortBy={{ title: 'Gene', prop: 'Gene', order: 'descending' }}
                pageLengthOptions={[ 10, 50, 100 ]}
                keys={['id']}
            />
		);
	}
});


module.exports = VariantTable;
