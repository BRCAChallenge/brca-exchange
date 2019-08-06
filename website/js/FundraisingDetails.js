/*eslint-env browser */
'use strict';

const React = require('react');
const {State} = require('react-router');
const RawHTML = require('./RawHTML');
const {Grid, Row, Col} = require('react-bootstrap');
const content = require('./content');

const FundraisingDetails = React.createClass({
    mixins: [State],


    render() {
        return (
            <Grid id="main-grid" className="help-page fundraising-details">
                <Row>
                    <Col smOffset={1} sm={10}>
                        <RawHTML html={content.pages.fundraisingDetails} />
                        <p className="small margin-top-forty"><sup>1 </sup>“<a href="https://www.cdc.gov/mmwr/volumes/66/ss/ss6615a1.htm">
                            BRCA Genetic Testing and Receipt of Preventive Interventions Among Women Aged 18–64 Years with
                            Employer-Sponsored Health Insurance in Nonmetropolitan and Metropolitan Areas — United States,
                            2009–2014</a>", Center for Disease Control and Prevention.
                        </p>
                        <p className="small"><sup>2 </sup>
                            “<a href="https://www.jax.org/education-and-learning/clinical-and-continuing-education/cancer-resources/hereditary-breast-and-ovarian-cancer-syndrome-factsheet">
                                Hereditary Breast and Ovarian Cancer Syndrome
                            </a>
                            ", The Jackson Laboratory.
                        </p>
                        <p className="small"><sup>3 </sup>“<a href="https://gnomad.broadinstitute.org">
                            The Genome Aggregation Database (gnomAD)</a>", The Broad Institute.
                        </p>
                        <p className="small"><sup>4 </sup><a href="https://www.census.gov/popclock/">
                            https://www.census.gov/popclock/</a>
                        </p>
                    </Col>
                </Row>
            </Grid>
        );
    }
});

module.exports = FundraisingDetails;
