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
require('bootstrap/dist/css/bootstrap.css');
require('font-awesome-webpack');
require('css/custom.css');
var _ = require('underscore');
var backend = require('./backend');
var Rx = require('rx');
require('rx-dom');

var brcaLogo = require('./img/BRCA-Exchange-tall-tranparent.png');
var betaBanner = require('./img/Beta_Banner.png');
var logos = require('./logos');
var slugify = require('./slugify');

var content = require('./content');

var databaseKey = require('../databaseKey');

var {Grid, Col, Row, Navbar, Nav, Table,
    DropdownButton, MenuItem, Modal, Button} = require('react-bootstrap');


var VariantTable = require('./VariantTable');
var VariantSearch = require('./VariantSearch');
var {Navigation, State, Link, Route, RouteHandler,
    HistoryLocation, run, DefaultRoute} = require('react-router');

var navbarHeight = 70; // XXX This value MUST match the setting in custom.css

var D3LollipopSVG = require('./D3LollipopSVG');

var variantPathJoin = row => _.map(databaseKey, k => encodeURIComponent(row[k])).join('@@');
var variantPathSplit = id => _.object(databaseKey, _.map(id.split(/@@/), decodeURIComponent));

if (typeof console === "undefined") {
    window.console = {
        log: function () {}
    };
}

var RawHTML = React.createClass({
    render: function() {
        var {html, ...otherProps} = this.props;
        return (
            <div className='markdown' {...otherProps} dangerouslySetInnerHTML={{__html: html}} />
        );
    }
});

var NavLink = React.createClass({
    render: function () {
        var {children, ...otherProps} = this.props;
        return (
            <li>
                <Link {...otherProps} role='button'>
                    {children}
                </Link>
            </li>
        );
    }
});

var NavBarNew = React.createClass({
    close: function () {
        this.refs.about.setState({open: false});
    },
    shouldComponentUpdate: function (nextProps) {
        // Only rerender if path has change, ignoring query.
        return this.props.path.split(/\?/)[0] !== nextProps.path.split(/\?/)[0];
    },
    activePath: function(path, tab) {
        var navPath = (path === "") ? "home" : path.split("/")[0];
        return ((navPath === tab) ? "active" : "");
    },
    render: function () {
        var {path} = this.props;
        var brand = (
            <a className="navbar-brand" href="http://brcaexchange.org">
                <span>
                    <b className="BRCA">BRCA</b>
                    <span className="exchange"> Exchange</span>
                    <img src={betaBanner} alt="beta banner" />
                </span>
            </a>);
        return (
            <div className="navbar-container">
            <Navbar fixedTop brand={brand} toggleNavKey={0}>
                <Nav eventKey={0} navbar right>
                    <NavLink to='/'>Home</NavLink>
                    <DropdownButton className={this.activePath(path, "about")} ref='about' title='About'>
                        <NavLink onClick={this.close} to='/about/history'>
                            History of the BRCA Exchange
                        </NavLink>
                        <NavLink onClick={this.close} to='/about/variation'>
                            BRCA1, BRCA2, and Cancer
                        </NavLink>
                        <NavLink onClick={this.close} to='/about/lollipop'>
                            DNA Variant BRCA Lollipop Plots
                        </NavLink>
                        <NavLink onClick={this.close} to='/about/thisSite'>
                            This Site
                        </NavLink>
                    </DropdownButton>
                    <NavLink to='/variants'>Variants</NavLink>
                    <NavLink to='/help'>Help</NavLink>
                </Nav>
            </Navbar>
            </div>
        );
    }
});

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
                        <li><DisclaimerModal /></li>
                        <li><a href="mailto:brca-exchange-contact@genomicsandhealth.org?subject=BRCA Exchange website">contact us</a></li>
                        <li>
                            <a href="https://github.com/BD2KGenomics/brca-website">
                                source code
                            </a>
                        </li>
                    </ul>
                </div>
            </div>
        );
    }
});

var DisclaimerModal = React.createClass({
    getInitialState() {
        return { showModal: false };
    },
    close() {
        this.setState({ showModal: false });
    },
    open() {
        this.setState({ showModal: true });
    },
    render() {
        return (
            <div style={{display: "inline"}}>
                <a onClick={this.open}>disclaimer</a>
                {this.state.showModal ?
                    <Modal onHide={this.close}>
                        <RawHTML html={content.pages.disclaimer} />
                        <div className = "close-button">
                            <Button onClick={this.close}>close</Button>
                        </div>
                    </Modal> : null }
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
            <Grid className='home'>
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
            <Grid>
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
    render: function() {
        var fragment = slugify(window.location.hash.slice(1));
        return (
            <Grid className="help">
                {fragment === '' ? null :
                    <style>{`#${fragment} { animation-name: emphasis; animation-duration: 10s; } `}</style>}
                <Row>
                    <Col md={8} mdOffset={2}>
                        <RawHTML ref='content' html={content.pages.help} />
                    </Col>
                </Row>
            </Grid>
        );
    }
});

// wrap scalars in array.
function toArray(v) {
    return _.isArray(v) ? v : [v];
}

function toNumber(v) {
    return _.isString(v) ? parseInt(v) : v;
}

function databaseParams(paramsIn) {
    var {filter, filterValue, hide} = _.mapObject(
            _.pick(paramsIn, 'hide', 'filter', 'filterValue'), toArray),
        numParams = _.mapObject(_.pick(paramsIn, 'page', 'pageLength'),
                toNumber),
        {orderBy, order, search = ''} = _.pick(paramsIn, 'search', 'orderBy', 'order'),
        sortBy = {prop: orderBy, order},
        columnSelection = _.object(hide, _.map(hide, _.constant(false))),
        filterValues = _.object(filter, filterValue);

    return {search, sortBy, columnSelection, filterValues, hide, ...numParams};
}

var transpose = a => _.zip.apply(_, a);

function urlFromDatabase(state) {
    // Need to diff from defaults. The defaults are in DataTable.
    // We could keep the defaults here, or in a different module.
    var {columnSelection, filterValues,
            search, page, pageLength, sortBy: {prop, order}} = state,
        hide = _.keys(_.pick(columnSelection, v => !v)),
        [filter, filterValue] = transpose(_.pairs(_.pick(filterValues, v => v)));
    return _.pick({
        search: search === '' ? null : search,
        filter,
        filterValue,
        page: page === 0 ? null : page,
        pageLength: pageLength === 20 ? null : pageLength,
        orderBy: prop,
        order,
        hide: hide.length === 0 ? null : hide
    }, v => v != null);

}

var Database = React.createClass({
    // Note this is not a pure component because of the calls to
    // getQuery().
    mixins: [Navigation, State],
    showVariant: function (row) {
        this.transitionTo(`/variant/${variantPathJoin(row)}`);
    },
    showHelp: function (title) {
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
            this.transitionTo('/variants', {}, urlFromDatabase(state));
        }
    },
    render: function () {
        var {show} = this.props,
            params = databaseParams(this.getQuery());
        // XXX is 'keys' used?
        return (
            <Grid style={{display: show ? 'block' : 'none'}}>
                <VariantTable
                    ref='table'
                    initialState={params}
                    {...params}
                    fetch={backend.data}
                    url={backend.url}
                    onChange={s => this.urlq.onNext(s)}
                    suggestions={[]}
                    keys={databaseKey}
                    onHeaderClick={this.showHelp}
                    onRowClick={this.showVariant}/>
            </Grid>
        );
    }
});

var _toSpace = s => s.replace(/_/g, ' ');

var Key = React.createClass({
    render() {
        var {onClick, tableKey} = this.props,
            words = tableKey.replace(/_/g, ' ').split(' ');
        return (
             <td className='help-target'>
                {words.slice(0, words.length - 1).join(' ')}
                {' '}
                <span className="text-nowrap">
                    {words[words.length - 1]}
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
        }).take(1).subscribe(
            resp => this.setState({data: resp.data[0], error: null}),
            this.setState({error: 'Problem connecting to server'}));
    },
    render: function() {
        var {data: variant = {}, error} = this.state;

        variant = _.omit(variant, ['__HEADER__']);
        var rows = _.map(variant, (v, k) =>
             <tr key={k}>
                <Key tableKey={k} onClick={() => this.showHelp(_toSpace(k))} />
                <td>{v}</td>
             </tr>);


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
    render: function () {
        var path = this.getPath().slice(1);
        return (
            <div>
                <NavBarNew path={path} />
                <RouteHandler />
                <Database
                    show={path.indexOf('variants') === 0} />
                <Footer />
            </div>
        );
    }
});

var D3StaticLollipop = React.createClass({
    render: function () {
        return (
            <D3LollipopSVG brcakey={this.props.brcakey}/>
        );
    }
});

var Lollipop = React.createClass({
    mixins: [Navigation],
    showHelp: function (title) {
        this.transitionTo(`/help#${slugify(title)}`);
    },
    getInitialState: function() {
        return {brcakey: "BRCA1"};
    },
    onSelect: function(key) {
        this.setState({brcakey: key});
    },
    render: function () {
        return (
            <Grid>
                <Row>
                    <Col md={8} mdOffset={4}>
                        <h1 id="brca-dna-variant-lollipop">{this.state.brcakey} Lollipop Chart</h1>
                    </Col>
                </Row>
                <div>
                    <DropdownButton onSelect={this.onSelect} title="Select Gene" id="bg-vertical-dropdown-1">
                        <MenuItem eventKey="BRCA1">BRCA1</MenuItem>
                        <MenuItem eventKey="BRCA2">BRCA2</MenuItem>
                    </DropdownButton>
                    <span onClick={() => this.showHelp('Lollipop Plots')}
                        className='help glyphicon glyphicon-question-sign superscript'/>
                    <D3StaticLollipop key={this.state.brcakey} brcakey={this.state.brcakey} id='brcaLollipop' ref='d3svgBrca'/>
                </div>
            </Grid>
        );
    }
});

var routes = (
    <Route handler={Application}>
        <DefaultRoute handler={Home}/>
        <Route path='about/lollipop' handler={Lollipop}/>
        <Route path='about/:page' handler={About}/>
        <Route path='help' handler={Help}/>
        <Route path='variants' />
        <Route path='variant/:id' handler={VariantDetail}/>
    </Route>
);

var main = document.getElementById('main');

run(routes, HistoryLocation, (Root) => {
  React.render(<Root/>, main);
});
