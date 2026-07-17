'use strict';

const React = require('react');
import RawHTML from './RawHTML';
import  { Container, Row, Col } from 'react-bootstrap';
const content = require('./content');

class FundraisingDetails extends React.PureComponent {
    render() {
        return (
            <Container id="main-grid" className="help-page fundraising-details">
                <Row>
                    <Col sm={{ span: 10, offset: 1}}>
                        <RawHTML html={content.pages.fundraisingDetails} />
                        <p className="small margin-top-forty"><sup>1 </sup>&quot;<a href="https://www.cdc.gov/mmwr/volumes/66/ss/ss6615a1.htm">
                            BRCA Genetic Testing and Receipt of Preventive Interventions Among Women Aged 18–64 Years with
                            Employer-Sponsored Health Insurance in Nonmetropolitan and Metropolitan Areas — United States,
                            2009–2014</a>&quot;, Center for Disease Control and Prevention.
                        </p>
                        <p className="small"><sup>2 </sup>
                            &quot;<a href="https://www.jax.org/education-and-learning/clinical-and-continuing-education/cancer-resources/hereditary-breast-and-ovarian-cancer-syndrome-factsheet">
                                Hereditary Breast and Ovarian Cancer Syndrome
                            </a>
                            &quot;, The Jackson Laboratory.
                        </p>
                        <p className="small"><sup>3 </sup>&quot;<a href="https://gnomad.broadinstitute.org">
                            The Genome Aggregation Database (gnomAD)</a>&quot;, The Broad Institute.
                        </p>
                        <p className="small"><sup>4 </sup><a href="https://www.census.gov/popclock/">
                            https://www.census.gov/popclock/</a>
                        </p>
                    </Col>
                </Row>
            </Container>
        );
    }
}

export default FundraisingDetails;
