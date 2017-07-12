import React from 'react';
import { Navigation, Link } from 'react-router';
import _ from 'underscore';
import moment from 'moment';
import { Grid, Col, Row, Table, Button, Panel } from 'react-bootstrap';

// faisal: includes for masonry/isotope
import Isotope from 'isotope-layout';

import 'isotope-packery';
import backend from '../backend';
import slugify from '../helpers/slugify';
import { columns, researchModeColumns, researchModeGroups, expertModeGroups } from '../data/columns';
import DisclaimerModal from '../partials/DisclaimerModal';

// replaces commas with comma-spaces to wrap long lines better, removes blank entries from comma-delimited lists,
// and normalizes blank/null values to a single hyphen
function normalizedFieldDisplay(value) {
    if (value) {
        // replace any number of underscores with spaces
        // make sure commas, if present, wrap
        value = value
            .split(/_+/).join(" ")
            .split(",")
            .map(x => x.trim())
            .filter(x => x && x !== '-')
            .join(", ");

        // ensure that blank entries are always normalized to hyphens
        if (value.trim() === "") {
            value = "-";
        }
    }
    else {
        // similar to above, normalize blank entries to a hyphen
        value = "-";
    }

    return value;
}

const KeyInline = React.createClass({
    render() {
        const {onClick, tableKey} = this.props;
        return (
            <td className='help-target'>
                <span className="help-target-inline" onClick={onClick}>{tableKey}</span>
            </td>
        );
    }
});

const GroupHelpButton = React.createClass({
    render() {
        const {onClick} = this.props;
        return (
            <span role='button' onClick={onClick} aria-label="Help"
                className='panel-help-btn glyphicon glyphicon-question-sign'
            />
        );
    }
});

// attempts to parse the given date string using a variety of formats,
// returning the formatted result as something like '08 September 2016'.
// just returns the input if every pattern fails to match
function normalizeDateFieldDisplay(value) {
    // extend this if there are more formats in the future
    const formats = ["MM/DD/YYYY", "YYYY-MM-DD"];

    for (let i = 0; i < formats.length; i++) {
        const q = moment(value, formats[i]);

        if (q.isValid()) {
            return q.format("DD MMMM YYYY");
        }
    }

    return value;
}

// get display name for a given key from VariantTable.js column specification,
// if we are in expert reviewed mode, search expert reviewed names then fall back to
// all data, otherwise go straight to all data. Finally, if key is not found, replace
// _ with space in the key and return that.
function getDisplayName(key) {
    const researchMode = (localStorage.getItem("research-mode") === 'true');
    let displayName;
    if (!researchMode) {
        displayName = columns.find(e => e.prop === key);
        displayName = displayName && displayName.title;
    }
    // we are not in expert reviewed more, or key wasn't found in expert reviewed columns
    if (displayName === undefined) {
        displayName = researchModeColumns.find(e => e.prop === key);
        displayName = displayName && displayName.title;
    }
    // key was not found at all
    if (displayName === undefined) {
        displayName = key.replace(/_/g, " ");
    }
    return displayName;
}

// test for the various forms of blank fields
function isEmptyField(value) {
    if (Array.isArray(value)) {
        value = value[0];
    }

    if (value === null || (typeof value === 'undefined')) {
        return true;
    }

    var v = value.trim();
    return v === '' || v === '-' || v === 'None';
}

function isEmptyDiff(value) {
    return value === null || value.length < 1;
}

const IsoGrid = React.createClass({
    displayName: 'IsoGrid',

    // Wrapper to layout child elements passed in
    render: function () {
        const children = this.props.children;
        return (
            <div className="isogrid">
                {children}
            </div>
        );
    },

    // When the DOM is rendered, let Masonry know what's changed
    componentDidUpdate: function() {
        if (this.masonry) {
            this.masonry.reloadItems();
            this.masonry.arrange();
        }
    },

    // Set up Masonry
    componentDidMount: function() {
        if(!this.masonry) {
            // i suppose we're doing this so the children exist when we create it?
            this.masonry = new Isotope('.isogrid', {
                layoutMode: 'packery',
                itemSelector: '.isogrid-item',
                packery: {
                    columnWidth: '.isogrid-sizer',
                    gutter: 0
                }
            });
        }
        else {
            this.masonry.reloadItems();
            this.masonry.arrange();
        }
    }
});

var VariantDetail = React.createClass({
    mixins: [Navigation],
    showHelp: function (event, title) {
        event.preventDefault();

        this.transitionTo(`/help#${slugify(title)}`);
    },
    getInitialState: () => ({
        hideEmptyItems: (localStorage.getItem("hide-empties") === 'true')
    }),
    componentWillMount: function () {
        backend.variant(this.props.params.id).subscribe(
            resp => {
                return this.setState({data: resp.data, error: null});
            },
            () => { this.setState({error: 'Problem connecting to server'}); });
    },
    onChildToggleMode: function() {
        this.props.toggleMode();
        this.forceUpdate();
    },
    reformatDate: function(date) { //handles single dates or an array of dates
        if (isEmptyField(date)) {
            return date;
        }
        if (!Array.isArray(date)) {
            date = date.split(',');
        }
        return date.map(function(d) {
            return moment.utc(new Date(d)).format("DD MMMM YYYY");
        }).join();
    },
    pathogenicityChanged: function(pathogenicityDiff) {
        return (pathogenicityDiff.added || pathogenicityDiff.removed) ? true : false;
    },
    setEmptyRowVisibility: function(hideEmptyItems) {
        localStorage.setItem('hide-empties', hideEmptyItems);

        this.setState({
            hideEmptyItems: hideEmptyItems
        });
    },
    truncateData: function(field) {
        const fieldsToTruncate = ["Genomic_Coordinate_hg38", "Genomic_Coordinate_hg37", "Genomic_Coordinate_hg36"];
        if (fieldsToTruncate.indexOf(field) > -1) {
            return true;
        } else {
            return false;
        }
    },
    onChangeGroupVisibility(groupTitle, event) {
        // stop the page from scrolling to the top (due to navigating to the fragment '#')
        event.preventDefault();

        // the event target is actually the span *inside* the 'a' tag, but we need to check the 'a' tag for the
        // collapsed state
        const collapsingElemParent = event.target.parentElement;

        let willBeCollapsed = true;

        collapsingElemParent.childNodes.forEach(function(child) {
            // FIXME: there must be a better way to get at the panel's state than reading the class
            // Maybe we'll subclass Panel and let it handle its own visibility persistence.
            if (child.getAttribute("class") === "collapsed") {
                // if it's already collapsed, this method should expand it
                willBeCollapsed = false;
            }
        });
        localStorage.setItem("collapse-group_" + groupTitle, willBeCollapsed);

        // defer re-layout until the state change has completed
        const me = this;

        setTimeout(() => {
            // this forces a re-render after a group has expanded/collapsed, fixing the layout
            // note that 300ms just happens to be the duration of the expand/collapse animation
            // it'd be better to run the re-layout whenever the animation ends
            me.forceUpdate();
        }, 300);
    },
    generateDiffRows: function(cols, data) {
        var diffRows = [];

        // keys that contain date values that need reformatting for the ui
        var dateKeys = [
            "Date_Last_Updated_ClinVar",
            "Date_last_evaluated_ENIGMA"
        ];

        // In research_mode, only show research_mode changes.
        var relevantFieldsToDisplayChanges = cols.map(function(col) {
            return col.prop;
        });

        for (var i = 0; i < data.length; i++) {
            let version = data[i];
            let diff = version.Diff;
            let release = version.Data_Release;
            let highlightRow = false;
            var diffHTML = [];

            if (diff !== null) {
                for (var j = 0; j < diff.length; j++) {
                    let fieldDiff = diff[j];
                    let fieldName = fieldDiff.field;
                    var added;
                    var removed;

                    if (fieldName === "Pathogenicity_expert") {
                        highlightRow = this.pathogenicityChanged(fieldDiff);
                    }

                    if (!_.contains(relevantFieldsToDisplayChanges, fieldName)) {
                        continue;
                    }

                    if (_.contains(dateKeys, fieldName)) {
                        added = this.reformatDate(fieldDiff.added);
                        removed = this.reformatDate(fieldDiff.removed);
                    } else if (fieldDiff.field_type === "list") {
                        added = _.map(fieldDiff.added, elem => elem.replace(/_/g, " ").trim());
                        removed = _.map(fieldDiff.removed, elem => elem.replace(/_/g, " ").trim());
                    } else {
                        added = fieldDiff.added.trim();
                        removed = fieldDiff.removed.trim();
                    }

                    if (added !== null || removed !== null) {
                        if (isEmptyField(removed)) {
                            diffHTML.push(
                                <span>
                                    <strong>{ getDisplayName(fieldName) }: </strong>
                                    <span className='label label-success'><span className='glyphicon glyphicon-star' /> New</span>
                                    &nbsp;{`${added}`}
                                </span>, <br />
                            );
                        } else if (fieldDiff.field_type === "list") {
                            diffHTML.push(
                                <span>
                                    <strong>{ getDisplayName(fieldName) }: </strong> <br />
                                    { !isEmptyDiff(added) && `+${added}` }{ !!(!isEmptyDiff(added) && !isEmptyDiff(removed)) && ', '}{ !isEmptyDiff(removed) && `-${removed}` }
                                </span>, <br />
                            );
                        } else if (fieldDiff.field_type === "individual") {
                            diffHTML.push(
                                <span>
                                    <strong>{ getDisplayName(fieldName) }: </strong>
                                    {removed} <span className="glyphicon glyphicon-arrow-right" /> {added}
                                </span>, <br />
                            );
                        }
                    }
                }
            }

            diffRows.push(
                <tr className={highlightRow ? 'danger' : ''}>
                    <td><Link to={`/release/${release.id}`}>{moment(release.date, "YYYY-MM-DDTHH:mm:ss").format("DD MMMM YYYY")}</Link></td>
                    <td>{version["Pathogenicity_expert"]}</td>
                    <td>{diffHTML}</td>
                </tr>
            );
        }

        return diffRows;
    },
    render: function () {
        const {data, error} = this.state;
        if (!data) {
            return <div />;
        }

        let variant = data[0],
            release = variant["Data_Release"],
            cols,
            groups;

        if (localStorage.getItem("research-mode") === 'true') {
            cols = researchModeColumns;
            groups = researchModeGroups;
        } else {
            cols = columns;
            groups = expertModeGroups;
        }

        // FAISAL: rather than directly map cols, we create a higher-level groups structure
        // the higher-level groups structure maps a subset of columns to that group

        const groupTables = _.map(groups, ({ groupTitle, innerCols }) => {
            let rowsEmpty = 0;

            // remove the BIC classification and importance fields unless the classification is 1 or 5
            if (groupTitle === 'Clinical Significance (BIC)') {
                const bicClass = variant['Clinical_classification_BIC'];

                if (bicClass !== 'Class 1' && bicClass !== 'Class 5') {
                    innerCols = innerCols.filter(x => x.prop !== 'Clinical_classification_BIC' && x.prop !== 'Clinical_importance_BIC');
                }
            }

            // now map the group's columns to a list of row objects
            const rows = _.map(innerCols, ({prop, title}) => {
                let rowItem;

                if (prop === "Protein_Change") {
                    title = "Abbreviated AA Change";
                }

                if (variant[prop] !== null) {
                    if (prop === "Gene_Symbol") {
                        rowItem = <i>{variant[prop]}</i>;
                    } else if (prop === "URL_ENIGMA") {
                        if (variant[prop].length) {
                            rowItem = <a target="_blank" href={variant[prop]}>link to multifactorial analysis</a>;
                        }
                    } else if (prop === "Assertion_method_citation_ENIGMA") {
                        rowItem = <a target="_blank" href="https://enigmaconsortium.org/library/general-documents/">Enigma Rules version Mar 26, 2015</a>;
                    } else if (prop === "Source_URL") {
                        if (variant[prop].startsWith("http://hci-exlovd.hci.utah.edu")) {
                            rowItem = <a target="_blank" href={variant[prop].split(',')[0]}>link to multifactorial analysis</a>;
                        }
                    } else if (prop === "Comment_on_clinical_significance_ENIGMA" || prop === "Clinical_significance_citations_ENIGMA") {
                        const pubmed = "http://ncbi.nlm.nih.gov/pubmed/";
                        rowItem = _.map(variant[prop].split(/PMID:? ?([0-9]+)/), piece =>
                            (/^[0-9]+$/.test(piece)) ? <a target="_blank" href={pubmed + piece}>PMID: {piece}</a> : piece );
                    } else if (prop === "HGVS_cDNA") {
                        rowItem = variant[prop].split(":")[1];
                    } else if (prop === "HGVS_Protein") {
                        rowItem = variant[prop].split(":")[1];
                    } else if (prop === "Date_last_evaluated_ENIGMA" && !isEmptyField(variant[prop])) {
                        // try a variety of formats until one works, or just display the value if not?
                        rowItem = normalizeDateFieldDisplay(variant[prop]);
                    } else {
                        rowItem = normalizedFieldDisplay(variant[prop]);
                    }
                } else if (prop === "HGVS_Protein_ID" && variant["HGVS_Protein"] !== null) {
                    rowItem = variant["HGVS_Protein"].split(":")[0];
                }

                const isEmptyValue = isEmptyField(variant[prop]);

                if (isEmptyValue) {
                    rowsEmpty += 1;
                    rowItem = '-';
                }

                return (
                    <tr key={prop} className={ (isEmptyValue && this.state.hideEmptyItems) ? "variantfield-empty" : "" }>
                        <KeyInline tableKey={title} onClick={(event) => this.showHelp(event, title)}/>
                        <td><span className={ this.truncateData(prop) ? "row-value-truncated" : "row-value" }>{rowItem}</span></td>
                    </tr>
                );
            });

            // check if all our rows are empty, in which case our group should be flagged as empty
            const allEmpty = rowsEmpty >= rows.length;

            const header = (
                <h3>
                    <a href="#" onClick={(event) => this.onChangeGroupVisibility(groupTitle, event)}>{groupTitle}</a>
                    <GroupHelpButton group={groupTitle} onClick={(event) => { this.showHelp(event, groupTitle); return true; }} />
                </h3>
            );

            return (
                <div key={`group_collection-${groupTitle}`} className={ allEmpty && this.state.hideEmptyItems ? "group-empty" : "" }>
                    <Panel
                        header={header}
                        collapsable={true}
                        defaultExpanded={localStorage.getItem("collapse-group_" + groupTitle) !== "true"}>
                        <Table>
                            <tbody>
                                {rows}
                            </tbody>
                        </Table>
                    </Panel>
                </div>
            );
        });

        const diffRows = this.generateDiffRows(cols, data);

        return (error ? <p>{error}</p> :
                <Grid>
                    <Row>
                        <Col xs={4} sm={4} smOffset={4} md={4} mdOffset={4} className="vcenterblock">
                            <div className='text-center Variant-detail-title'>
                                <h3>Variant Detail</h3>
                            </div>
                        </Col>
                        <Col xs={8} sm={4} md={4} className="vcenterblock">
                            <div className="Variant-detail-headerbar">
                                <Button
                                    onClick={this.setEmptyRowVisibility.bind(this, !this.state.hideEmptyItems)}
                                    bsStyle={"default"}>
                                    { this.state.hideEmptyItems ?
                                        <span>Show Empty Items</span> :
                                        <span>Hide Empty Items</span>
                                    }
                                </Button>
                            </div>
                        </Col>
                        {variant['Change_Type'] === 'deleted' &&
                        (<Col xs={12} className="vcenterblock">
                            <p className='deleted-variant-message'>
                                Note: This variant has been removed from the BRCA Exchange. For reasons on why this variant was removed please see the <Link to={`/release/${release.id}`}>release notes</Link>.
                            </p>
                        </Col>)
                        }
                    </Row>

                    <Row>
                        <div className="container-fluid variant-details-body">
                            { (groupTables.length < 3) ?
                                <IsoGrid>
                                    <div className={`isogrid-sizer col-xs-12 col-md-${12 / groupTables.length}`}/>
                                    {
                                        groupTables.map((x, i) => {
                                            return (
                                                <Col key={"group_col-" + i} xs={12} md={12 / groupTables.length}
                                                    className="variant-detail-group isogrid-item">
                                                    {x}
                                                </Col>
                                            );
                                        })
                                    }
                                </IsoGrid>
                                :
                                <IsoGrid>
                                    <div className="isogrid-sizer col-xs-12 col-md-6 col-lg-6 col-xl-4"/>
                                    {
                                        // we're mapping each group into a column so we can horizontally stack them
                                        groupTables.map((x, i) => {
                                            return (
                                                <Col key={"group_col-" + i} xs={12} md={6} lg={6}
                                                    className="variant-detail-group isogrid-item col-xl-4">
                                                    {x}
                                                </Col>
                                            );
                                        })
                                    }
                                </IsoGrid>
                            }
                        </div>
                    </Row>

                    <Row>
                        <Col md={12} className="variant-history-col">
                            <h3>{variant["HGVS_cDNA"]}</h3>
                            <h4>Previous Versions of this Variant:</h4>
                            <Table className='variant-history nopointer' bordered>
                                <thead>
                                    <tr className='active'>
                                        <th>Release Date</th>
                                        <th>Clinical Significance</th>
                                        <th>Changes</th>
                                    </tr>
                                </thead>
                                <tbody>
                                    {diffRows}
                                </tbody>
                            </Table>
                            <p style={{display: this.props.mode === "research_mode" ? 'none' : 'block' }}>There may be additional changes to this variant, click "Show All Public Data on this Variant" to see these changes.</p>
                        </Col>
                    </Row>
                    <Row>
                        <Col md={12} mdOffset={0}>
                            <DisclaimerModal buttonModal onToggleMode={this.onChildToggleMode} text="Show All Public Data on this Variant"/>
                        </Col>
                    </Row>
                </Grid>
        );
    }
});

export default VariantDetail;
