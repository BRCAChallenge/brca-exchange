/*eslint-env browser */
/*global require: false */
'use strict';

// shims for older browsers
require('babel/polyfill');
require('es5-shim');
require('es5-shim/es5-sham');

require('./favicons');
var React = require('react');
var PureRenderMixin = require('./PureRenderMixin'); // deep-equals version of PRM
var DisclaimerModal = require('./DisclaimerModal');
var RawHTML = require('./RawHTML');
require('bootstrap/dist/css/bootstrap.css');
require('font-awesome-webpack');
require('css/custom.css');
var _ = require('underscore');
var backend = require('./backend');
var {NavBarNew} = require('./NavBarNew');
var Rx = require('rx');
require('rx-dom');

var brcaLogo = require('./img/BRCA-Exchange-tall-tranparent.png');
var logos = require('./logos');
var slugify = require('./slugify');

var content = require('./content');
var Community = require('./Community');

var databaseKey = require('../databaseKey');

var {Grid, Col, Row, Table} = require('react-bootstrap');

var {VariantTable, ResearchVariantTable, research_mode_columns, columns} = require('./VariantTable');
var {Signup} = require('./Signup');
var {Signin} = require('./Signin');
var {Profile} = require('./Profile');
var VariantSearch = require('./VariantSearch');
var {Navigation, State, Route, RouteHandler,
    HistoryLocation, run, DefaultRoute} = require('react-router');

var navbarHeight = 70; // XXX This value MUST match the setting in custom.css

var variantPathJoin = row => _.map(databaseKey, k => encodeURIComponent(row[k])).join('@@');
var variantPathSplit = id => _.object(databaseKey, _.map(id.split(/@@/), decodeURIComponent));

if (typeof console === "undefined") {
    window.console = {
        log: function () {}
    };
}

var Footer = React.createClass({
    mixins: [PureRenderMixin],
    render: function() {
        return (
            <div className="container footer">
                <div className="col-sm-5 left-footer">
                    <ul>
                        <li><a href="/home">Home</a></li>
                        <li><a href="/about/history">About</a></li>
                        <li><a href="/variants">Variants</a></li>
                        <li><a href="/help">Help</a></li>
                    </ul>
                </div>
                <div className="col-sm-2 logo-footer">
                    <img href="#" src={brcaLogo} alt="brca exchange logo" />
                </div>
                <div className="col-sm-5 right-footer">
                    <ul>
                        <li>
                            <DisclaimerModal text="Disclaimer"/>
                        </li>
                        <li>
                            <a href="mailto:brca-exchange-contact@genomicsandhealth.org?subject=BRCA Exchange website">
                                Contact us
                            </a>
                        </li>
                        <li>
                            <a href="https://github.com/BD2KGenomics/brca-website">
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
            direction: null
        };
    },

    onSearch(value) {
        this.transitionTo('/variants', null, {search: value});
    },
    render: function() {
        var {suggestions} = this.props;
        var logoItems = _.map(logos, ({id, logo, url}) => (
            <li key={id}><a href={url}>
                <img id={id} src={logo} alt={id + ' logo'} />
            </a></li>
        ));
        return (
            <Grid id="main-grid" className='home'>
                <Row>
                    <Col smOffset={3} sm={6}>
                        <VariantSearch
                            id='home-search'
                            suggestions={suggestions}
                            onSearch={this.onSearch}/>
                    </Col>
                </Row>
                <Row>
                    <div className="jumbotron">
                        <RawHTML html={content.pages.home} />
                    </div>
                </Row>
                <Row>
                    <div>
                        <ul>
                            <div className="container logo-block">
                                {logoItems}
                            </div>
                        </ul>
                    </div>
                </Row>
            </Grid>
        );
    }
});

var About = React.createClass({
    render: function() {
        var {page} = this.props.params;

        return (
            <Grid id="main-grid" className="main-grid">
                <Row>
                    <Col md={8} mdOffset={2}>
                        <RawHTML html={content.pages[page]} />
                    </Col>
                </Row>
            </Grid>
        );
    }
});

var Help = React.createClass({
    mixins: [State],
    componentDidMount: function () {
        var fragment = slugify(window.location.hash.slice(1));
        if (fragment !== '') {
            setTimeout(function () {
                var el = document.getElementById(fragment);
                if (el) {
                    window.scrollTo(0, el.getBoundingClientRect().top - navbarHeight);
                }
            }, 0);
        }
    },
    render: function () {
        var fragment = slugify(window.location.hash.slice(1));
        var helpContent;
        if (localStorage.getItem("research-mode") === 'true') {
            helpContent = content.pages.help_research;
        } else {
            helpContent = content.pages.help;
        }
        return (
            <Grid id="main-grid" className="help">
                {fragment === '' ? null :
                    <style>{`#${fragment} { animation-name: emphasis; animation-duration: 10s; } `}</style>}
                <Row>
                    <Col md={8} mdOffset={2}>
                        <RawHTML ref='content' html={helpContent}/>
                    </Col>
                </Row>
            </Grid>
        );
    }
});

function toNumber(v) {
    return _.isString(v) ? parseInt(v) : v;
}

function databaseParams(paramsIn) {
    var {filter, filterValue, hide, hideSources, excludeSources, orderBy, order, search = ''} = paramsIn;
    var numParams = _.mapObject(_.pick(paramsIn, 'page', 'pageLength'), toNumber);
    var sortBy = {prop: orderBy, order};
    var columnSelection = _.object(hide, _.map(hide, _.constant(false)));
    var sourceSelection = {..._.object(hideSources, _.map(hideSources, _.constant(0))),
                           ..._.object(excludeSources, _.map(excludeSources, _.constant(-1)))};
    var filterValues = _.object(filter, filterValue);
    return {search, sortBy, columnSelection, sourceSelection, filterValues, hide, ...numParams};
}

var transpose = a => _.zip.apply(_, a);

function urlFromDatabase(state) {
    // Need to diff from defaults. The defaults are in DataTable.
    // We could keep the defaults here, or in a different module.
    var {columnSelection, filterValues, sourceSelection,
            search, page, pageLength, sortBy: {prop, order}} = state;
    var hide = _.keys(_.pick(columnSelection, v => v == false));
    var hideSources = _.keys(_.pick(sourceSelection, v => v == 0));
    var excludeSources = _.keys(_.pick(sourceSelection, v => v == -1));
    var [filter, filterValue] = transpose(_.pairs(_.pick(filterValues, v => v == true)));
    return _.pick({
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
    }, v => v != null);

}

var Database = React.createClass({
    // Note this is not a pure component because of the calls to
    // getQuery().
    mixins: [Navigation, State],
	showVariant: function (row) {
        row.Genomic_Coordinate_hg38 = backend.trimSearchTerm(row.Genomic_Coordinate_hg38);
        var d3TipDiv = document.getElementsByClassName('d3-tip-selection');
        if (d3TipDiv.length != 0 && d3TipDiv[0].style.opacity != '0') {
            d3TipDiv[0].style.opacity='0';
            d3TipDiv[0].style.pointerEvents='none';
        }
        this.transitionTo(`/variant/${variantPathJoin(row)}`);
	},
    showHelp: function (title) {
        var d3TipDiv = document.getElementsByClassName('d3-tip-selection');
        if (d3TipDiv.length != 0 && d3TipDiv[0].style.opacity != '0') {
            d3TipDiv[0].style.opacity='0';
            d3TipDiv[0].style.pointerEvents='none';
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
            if (d3TipDiv.length != 0 && d3TipDiv[0].style.opacity != '0') {
                d3TipDiv[0].style.opacity='0';
                d3TipDiv[0].style.pointerEvents='none';
            }
            this.transitionTo('/variants', {}, urlFromDatabase(state));
        }
    },
    render: function () {
        var {show} = this.props,
            params = databaseParams(this.getQuery());
        // XXX is 'keys' used?
        var table, message;
        if (this.props.mode === 'research_mode') {
            table = <ResearchVariantTable
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
                onRowClick={this.showVariant}/>
            message = this.renderMessage(content.pages.variantsResearch)
        } else {
            table = <VariantTable
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
                onRowClick={this.showVariant}/>
            message = this.renderMessage(content.pages.variantsDefault)
        }
        return (
            <Grid fluid={true} id="main-grid" style={{display: show ? 'block' : 'none'}}>
                {message}
                {table}
            </Grid>
        );
    },
    renderMessage: function(message) {
        return  <Row>
                    <Col sm={10} smOffset={1}  className="alert alert-warning">
                        <RawHTML ref='content' html={message}/>
                    </Col>
                </Row>
    }
});

var Key = React.createClass({
    render() {
        var {onClick, tableKey} = this.props;
        return (
             <td className='help-target'>
                {tableKey}
                <span className="text-nowrap">
                    <span role='button' onClick={onClick}
                        className='help glyphicon glyphicon-question-sign superscript'/>
                </span>
             </td>
        );
    }
});

var VariantDetail = React.createClass({
    mixins: [Navigation],
    showHelp: function (title) {
        this.transitionTo(`/help#${slugify(title)}`);
    },
    componentWillMount: function () {
        backend.data({
            filterValues: variantPathSplit(this.props.params.id),
            pageLength: 1
        }).subscribe(
            resp => {
                return this.setState({data: resp.data[0], error: null})
            },
            this.setState({error: 'Problem connecting to server'}));
    },
    render: function () {
        var {data: variant = {}, error} = this.state;
        var cols;
        if (localStorage.getItem("research-mode") === 'true') {
            cols = research_mode_columns;
        } else {
            cols = columns;
        }
        var rows = _.map(cols, ({prop, title}) => {
            var row_item;
            if ((prop == "URL_ENIGMA" || prop == "Assertion_method_citation_ENIGMA") && variant[prop] != null) {
                row_item = <a href={variant[prop]}>{variant[prop]}</a>
            } else if (prop == "Source_URL" && variant[prop] != null) {
                row_item = _.map(variant[prop].split(','), url => <a href={url}>{url}</a>)
            } else {
                row_item = variant[prop]
            }
            return <tr key={prop}>
                <Key tableKey={title} columns={cols} onClick={() => this.showHelp(title)}/>
                <td><span className="row-wrap">{row_item}</span></td>
            </tr>
        });


        return (error ? <p>{error}</p> :
            <Grid>
                <Row>
                    <div className='text-center Variant-detail-title'>
                        <h3>Variant Detail</h3>
                    </div>
                </Row>
                <Row>
                    <Col md={8} mdOffset={2}>
                        <Table striped bordered>
                            <tbody>
                                {rows}
                            </tbody>
                        </Table>
                    </Col>
                </Row>
                <Row>
                    <Col md={8} mdOffset={2}>
                        <DisclaimerModal buttonModal onToggleMode={this} text="Research information on this variant"/>
                    </Col>
                </Row>
            </Grid>
        );
    }
});

// XXX implement in server
//var dontSuggest = [
//    'Assertion_method_citation',
//    'URL'
//];

//var flatmap = (coll, fn) => _.flatten(_.map(coll, fn), true);
//var minSuggestion = 3; // minimum length of string to use in autocomplete
//var rowWords = row => flatmap(_.values(_.omit(row, dontSuggest)),
//        v => v.toLowerCase().split(/\s+/));

// Pull out interesting strings from the data, for use in
// auto-completion.
//function getSuggestions(data) {
//    return _.uniq(flatmap(data, row =>
//                _.filter(rowWords(row), w => w.length >= minSuggestion)).sort(),
//            true);
//}

var Application = React.createClass({
    mixins: [State],
    onChildToggleMode: function() {
        this.toggleMode();
    },
    getInitialState: function () {
        return {
            mode: (localStorage.getItem("research-mode") === 'true') ? 'research_mode' : 'default',
        }
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
        var path = this.getPath().slice(1);
        return (
            <div>
                <NavBarNew path={path} mode={this.state.mode} toggleMode={this.onChildToggleMode}/>
                <RouteHandler />
                <Database
                    mode={this.state.mode}
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
        <Route path='help' handler={Help}/>
        <Route path='community' handler={Community}/>
        <Route path='signup' handler={Signup}/>
        <Route path='signin' handler={Signin}/>
        <Route path='profile' handler={Profile}/>
        <Route path='variants' />
        <Route path='variant/:id' handler={VariantDetail}/>
    </Route>
);

var main = document.getElementById('main');

run(routes, HistoryLocation, (Root) => {
  React.render(<Root/>, main);
});
