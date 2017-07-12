// A table for variants.
//
// The intent here was to split the generic table code
// in DataTable from the variant domain knowledge, which
// would be here. That division has broken down due to
// peculiarities of react-data-partials DataMixin, and
// time pressure. Knowledge about variants is in both files.
// This needs to be revisited.

/*global module: false, require: false, window: false */

var React = require('react');
var {State} = require('react-router');

var _ = require('underscore');
var {Col, Panel, Button, Input} = require('react-bootstrap');
var ColumnCheckbox = require('./ColumnCheckbox');

var {getDefaultExpertColumns, getDefaultResearchColumns, getAllSources} = require('../data/VariantTableDefaults');
var PureRenderMixin = require('../helpers/PureRenderMixin');
var DataTable = require('./DataTable');

require('react-data-components-bd2k/css/table-twbs.css');

var {
    filterColumns,
    subColumns,
    columns,
    researchModeColumns
} = require("../data/columns");

/*eslint-enable camelcase */

// Work-around to allow the user to select text in the table. The browser does not distinguish between
// click and drag: if mouseup and mousedown occur on the same element, a click event is fired even if
// the events occur at very different locations. That makes it hard to select text. This workaround
// defeats the click event if text has been selected.
//
// XXX getSelection().isCollapsed is not available on all platforms. On those platforms we
// will always return false (no selection), so the row click will fire. This makes it hard
// for the user to select text in the table. A better solution would be to add a polyfill for
// getSelection and isCollapsed. There are a few available, though they are much larger than
// what we require:
// https://github.com/Modernizr/Modernizr/wiki/HTML5-Cross-Browser-Polyfills#dom-range-and-selection
// We might want to write a minimal isCollapsed that will use whichever DOM method is available. We could
// also add feature detection, and modify the UI if the feature is not available. E.g. we could style the
// text areas to look clickable instead of selectable..
var hasSelection = () => !(window.getSelection && window.getSelection().isCollapsed);

function buildHeader(onClick, title) {
    return (
        <span>
            {title}
            <span onClick={ev => {ev.stopPropagation(); onClick(title); }}
                className='help glyphicon glyphicon-question-sign superscript'/>
        </span>
    );
}

var Table = React.createClass({
    mixins: [PureRenderMixin],
    render: function () {
        // Expert portal always shows all sources and default columns
        var {data, onHeaderClick, onRowClick, mode, columnSelection, sourceSelection, ...opts} = this.props;
        if (mode === "default") {
            // Always show default columns/sources in expert mode.
            columnSelection = _.object(_.map(columns,
                            c => _.contains(getDefaultExpertColumns(), c.prop) ? [c.prop, true] : [c.prop, false])
                        );
            sourceSelection = getAllSources();
        }
        return (
            <DataTable
                ref='table'
                className='row-clickable data-table table-grayheader'
                {...opts}
                columnSelection={columnSelection}
                sourceSelection={sourceSelection}
                buildRowOptions={r => ({className: r['Change_Type_id'] === 2 ? 'warning data-table-row' : 'data-table-row', title: 'click for details', onClick: () => hasSelection() ? null : onRowClick(r)})}
                buildHeader={title => buildHeader(onHeaderClick, title)}
                onRowClick={onRowClick}
                onHeaderClick={onHeaderClick}
                filterColumns={filterColumns}
                initialData={data}
                initialPageLength={20}
                initialSortBy={{prop: 'Gene_Symbol', order: 'descending'}}
                pageLengthOptions={[ 20, 50, 100 ]}
                mode={mode}/>
        );
    }
});

var ResearchVariantTableSupplier = function (Component) {
    return React.createClass({
        mixins: [State, PureRenderMixin],

        getInitialState: function () {
            /*
             Selections take the following order of priority:
             1. Query params in URL
             2. Local Storage
             3. Default settings

             To accomplish this, first set all column selections to true and all filters off.
             Then, update selections according to query params if present. If no query params
             are present, update according to local storage if present. If neither query params
             nor local storage specify changes, use default settings.
             */

            // Start with all columns set to true and data showing from all sources.
            var selectedColumns = _.object(_.map(this.getColumns(),
                c => [c.prop, true])
            );
            var selectedSources = getAllSources();

            // Get query params.
            const urlParams = this.getQuery();
            const useQueryParams = urlParams.hasOwnProperty("hide") || urlParams.hasOwnProperty("hideSources");

            if (useQueryParams) {
                // If query params are present, use them for settings.
                if (urlParams.hasOwnProperty("hide")) {
                    const columnsToHide = urlParams.hide;
                    for (let i = 0; i < columnsToHide.length; i++) {
                        selectedColumns[columnsToHide[i]] = false;
                    }
                }
                if (urlParams.hasOwnProperty("hideSources")) {
                    const sourcesToHide = urlParams.hideSources;
                    for (let i = 0; i < sourcesToHide.length; i++) {
                        selectedSources[sourcesToHide[i]] = 0;
                    }
                }
            } else {
                // If no query params are present, check local storage.
                const lsSelectedColumns = JSON.parse(localStorage.getItem('columnSelection'));
                if (lsSelectedColumns !== null && lsSelectedColumns !== undefined) {
                    selectedColumns = lsSelectedColumns;
                } else {
                    // If no query params and no local storage, use default settings.
                    selectedColumns = this.getDefaultColumnSelections();
                }
                const lsSelectedSources = JSON.parse(localStorage.getItem('sourceSelection'));
                if (lsSelectedSources !== null && lsSelectedSources !== undefined) {
                    selectedSources = lsSelectedSources;
                }
            }

            return {
                sourceSelection: selectedSources,
                columnSelection: selectedColumns
            };
        },
        toggleColumns: function (prop) {
            let {columnSelection} = this.state,
                val = columnSelection[prop],
                cs = {...columnSelection, [prop]: !val};
            localStorage.setItem('columnSelection', JSON.stringify(cs));
            this.setState({columnSelection: cs});
        },
        setSource: function (prop, event) {
            // this function uses 1, 0 and -1 to accommodate excluding sources as well as not-including them
            // currently only uses 1 and 0 because exclusion is not being used
            let {sourceSelection} = this.state;
            let value = event.target.checked ? 1 : 0;
            let ss = {...sourceSelection, [prop]: value};
            localStorage.setItem('sourceSelection', JSON.stringify(ss));
            this.setState({sourceSelection: ss});
        },
        filterFormCols: function (subColList, columnSelection) {
            return _.map(subColList, ({title, prop}) =>
                <ColumnCheckbox onChange={() => this.toggleColumns(prop)} key={prop} label={prop} title={title}
                    initialCheck={columnSelection}/>);
        },
        onChangeSubcolVisibility(subColTitle, event) {
            // stop the page from scrolling to the top (due to navigating to the fragment '#')
            event.preventDefault();

            const collapsingElem = event.target;

            // FIXME: there must be a better way to get at the panel's state than reading the class
            // maybe we'll subclass Panel and let it handle its own visibility persistence

            const isCollapsed = (collapsingElem.getAttribute("class") === "collapsed");
            localStorage.setItem("collapse-subcol_" + subColTitle, !isCollapsed);
        },
        getColumnSelectors() {
            var filterFormSubCols = _.map(subColumns, ({subColTitle, subColList}) =>
                <Col sm={6} md={4} key={subColTitle}>
                    <Panel
                        header={subColTitle}
                        collapsable={true}
                        defaultExpanded={localStorage.getItem("collapse-subcol_" + subColTitle) !== "true"}
                        onSelect={(event) => this.onChangeSubcolVisibility(subColTitle, event)}>
                        {this.filterFormCols(subColList, this.state.columnSelection)}
                    </Panel>
                </Col>
            );
            return (<label className='control-label'>
                <Panel header="Column Selection">
                    {filterFormSubCols}
                </Panel>
            </label>);
        },
        getFilters: function () {
            var sourceCheckboxes = _.map(this.state.sourceSelection, (value, name) =>
                <Col sm={6} md={3} key={name}>
                    <Input type="checkbox"
                        onChange={v => this.setSource(name, v)}
                        label={name.substring(11).replace(/_/g, " ")} // eg "Variant_in_1000_Genomes" => "1000 Genomes"
                        checked={value > 0}/>
                </Col>
            );
            return (<label className='control-label source-filters'>
                <Panel className="top-buffer" header="Source Selection">
                    {sourceCheckboxes}
                </Panel>
            </label>);
        },
        getDownloadButton: function (callback) {
            return (<Button className="btn-default rgt-buffer" download="variants.tsv"
                href={callback()}>Download</Button>);
        },
        getLollipopButton: function (callback, isOpen) {
            return (<Button id="lollipop-chart-toggle" className="btn-default rgt-buffer"
                onClick={callback}>{(isOpen ? 'Hide' : 'Show' ) + ' Lollipop Chart'}</Button>);
        },
        getColumns: function () {
            return researchModeColumns;
        },
        getDefaultColumnSelections: function () {
            return _.object(_.map(researchModeColumns,
                c => _.contains(getDefaultResearchColumns(), c.prop) ? [c.prop, true] : [c.prop, false])
            );
        },
        getDefaultSourceSelections: function () {
            return getAllSources();
        },
        researchVariantTableRestoreDefaults: function (callback) {
            const columnSelection = this.getDefaultColumnSelections();
            const sourceSelection = this.getDefaultSourceSelections();
            this.setState({
                    columnSelection: columnSelection,
                    sourceSelection: sourceSelection
                },
                function () {
                    this.props.restoreDefaults(callback);
                });
        },
        render: function () {
            return (
                <Component
                    {...this.props}
                    researchVariantTableRestoreDefaults={this.researchVariantTableRestoreDefaults}
                    columns={this.getColumns()}
                    columnSelectors={this.getColumnSelectors()}
                    filters={this.getFilters()}
                    sourceSelection={this.state.sourceSelection}
                    columnSelection={this.state.columnSelection}
                    downloadButton={this.getDownloadButton}
                    lollipopButton={this.getLollipopButton}/>
            );
        }
    });
};

var VariantTableSupplier = function (Component) {
    return React.createClass({
        mixins: [PureRenderMixin],
        getColumns: function () {
            return columns;
        },
        expertVariantTableRestoreDefaults: function (callback) {
            this.props.restoreDefaults(callback);
        },
        render: function () {
            let expertColumns = _.object(_.map(this.getColumns(),
                c => _.contains(getDefaultExpertColumns(), c.prop) ? [c.prop, true] : [c.prop, false])
            );
            // Expert portal always shows all sources
            let sourceSelection = getAllSources();
            return (
                <Component
                    {...this.props}
                    columns={this.getColumns()}
                    columnSelection={expertColumns}
                    sourceSelection={sourceSelection}
                    expertVariantTableRestoreDefaults={this.expertVariantTableRestoreDefaults}
                    downloadButton={() => null}
                    lollipopButton={() => null}/>
            );
        }
    });
};


module.exports = {
    VariantTable: VariantTableSupplier(Table),
    ResearchVariantTable: ResearchVariantTableSupplier(Table)
};
