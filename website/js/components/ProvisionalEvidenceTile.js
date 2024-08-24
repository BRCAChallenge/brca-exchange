'use strict';
import util from '../util';
import React from 'react';
import {Table} from "react-bootstrap";
import * as _ from 'lodash';
import classNames from "classnames";
import CollapsibleTile from "./collapsible/CollapsibleTile";
import CollapsibleSection from "./collapsible/CollapsibleSection";

// ACMG Variant Evidence Codes, Provisional Assignment
export default class ProvisionalEvidenceTile extends React.Component {

    // For the sub-tile headers that have a result value displayed in the header,
    // Add css so that the value floats to the right
    generateHeader(result) {
        return (
            <div className="func-assay-extras">
                <span className="func-assay-result" style={{float: 'right', paddingRight: '10px'}}>{result}</span>
                <div style={{clear: 'both'}}></div>
            </div>
        );
    }

    // Given a group source (ie groupname), data (list of dicts containing title:, prop:), and variant object
    // generate table rows to be displayed in that group's section
    // retrieving the requested prop values from the variant
    // and applying the corresponding help tooltips
    getRowsAndDetermineIfEmpty(source, data, variant) {
        const rows = _.map(data, (rowDescriptor) => {
            //let {prop, title, noHelpLink} = rowDescriptor;
            let {prop, title} = rowDescriptor;

            const rowItem = util.getFormattedFieldByProp(prop, variant);
            let rowClasses = classNames({
                'variantfield-empty': false, // placeholder until this supports hiding when empty
                // 'variantfield-empty': (isEmptyValue && this.props.hideEmptyItems),
            });
            return (
                <tr key={prop} className={rowClasses}>
                    <td>{title}</td>
                    <td><span className={"row-value" }>{rowItem}</span></td>
                </tr>
            );
        });
        return rows;
    }

    render() {
        const {variant, innerGroups} = this.props;
        let allEmpty = true;

        // innerGroup are: Population Frequency, Computational Prediction
        let sections = _.map(innerGroups, (group) => {
            let groupSource = group.source;
            let groupData = group.data;
            let renderedRows = this.getRowsAndDetermineIfEmpty(groupSource, groupData, variant);

            // TEMPORARY - hide ComputationalPrediction subsection until data is populated
            if (groupSource === "Computational Prediction") {
                return ( <div id={groupSource}/> );
            }

            // Attempt to retrieve the prop name of the provisional code so we can put the value in the subsection header
            let extraHeaderProp = false ;
            for ( const dataItem of groupData ) {
                if (dataItem.title === "Provisionally Assigned" ) {
                    extraHeaderProp = dataItem.prop ;
                }
            }
            // If we found a prop, pull it from the variant. If not, leave as a dash.
            let extraHeaderValue = extraHeaderProp ? util.getFormattedFieldByProp(extraHeaderProp, variant) : "-";

            return ( <CollapsibleSection
                    fieldName={groupSource}
                    id={groupSource} // id needed so that when you collapse one the others don't also collapse
                    extraHeaderItems={this.generateHeader(extraHeaderValue)}
                    twoColumnExtraHeader={true}
                    >
                        <Table>
                            <tbody>
                                { renderedRows }
                            </tbody>
                        </Table>
                    </CollapsibleSection>);
        });

        return (
            <CollapsibleTile allEmpty={allEmpty} {...this.props}>
                <div className="tile-disclaimer">
                    The ClinGen <a href="https://cspec.genome.network/cspec/ui/svi/affiliation/50087">
                        ENIGMA Variant Curation Expert Panel (VCEP) rules (Version 1.1.0, dated 2023-11-22)
                    </a> define how the evidence on a variant can contribute to the variant being curated as
                    Pathogenic or Benign, following the <a href="https://pubmed.ncbi.nlm.nih.gov/25741868/">
                        ACMG/AMP standards and guidelines
                    </a>. We have evaluated the following categories of evidence
                    for each variant against the VCEP rules, to evaluate which evidence codes the variant meets.
                    <b> These are provisional evidence code assignments. In general, they have not yet been reviewed by the
                    VCEP biocuration team.</b>  If the variant has been curated by the VCEP, then the
                    Variant Curation Expert Panel tile will contain information on the evidence used for curation.
                </div>
                {sections}
            </CollapsibleTile>
        );
    }
};

