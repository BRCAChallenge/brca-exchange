/*global module: false, require: false, console: false */
'use strict';

var React = require('react');
var DataTable = require('react-data-components-bd2k').DataTable;
require('react-data-components-bd2k/css/table-twbs.css');

function buildHeader(onClick, title) {
	return (
		<span>
			{title}
			<span onClick={() => onClick(title)}
				className='help glyphicon glyphicon-question-sign superscript'/>
		</span>
	);
}

var VariantTable = React.createClass({
	render: function () {
		var {data, onHeaderClick, onRowClick} = this.props,
			columns = [
				{title: 'Gene', prop: 'Gene symbol'},
				{title: '  HGVS  ', prop: 'HGVS'},
				{title: 'Pathogenicity', prop: 'Clinical significance'},
				{title: 'Allele origin', prop: 'Allele origin'},
				{title: 'CVA', prop: 'ClinVarAccession'}
			];
		return (
			<DataTable
				buildRowOptions={r => ({onClick: () => onRowClick(r.id)})}
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
