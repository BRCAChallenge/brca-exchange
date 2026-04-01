/*eslint-env browser */
'use strict';

import React from 'react';
import RawHTML from './RawHTML';
import { Container as Grid, Row, Col, Button} from 'react-bootstrap';
import { Link } from 'react-router-dom';
const content = require('./content');

class WhyDonate extends React.Component {
    render() {
        return (
            <Grid id="main-grid" className="help-page">
                <Row>
                    <Col smOffset={1} sm={10}>
                        <RawHTML html={content.pages.whyDonate} />
                        <h3 className="centered margin-top-forty">Will you help us improve genetic data and variant interpretations worldwide by donating today?</h3>
                        <Button variant="primary" className="d-block mx-auto donate-button" onClick={()=> window.open("https://secure.ucsc.edu/s/1069/bp18/interior.aspx?sid=1069&gid=1001&pgid=780&cid=1749&dids=1004", "_blank")}>
                            Donate Now
                        </Button>
                        <p className="centered small">
                            Your donations are tax-deductible. Donations to support BRCA Exchange are made through
                            the UCSC Foundation (Tax ID: 23-7394590), a registered 501(c)3.
                        </p>
                        <p className="small"><sup>1 </sup>"<a href="https://www.cdc.gov/mmwr/volumes/66/ss/ss6615a1.htm">
                            BRCA Genetic Testing and Receipt of Preventive Interventions Among Women Aged 18–64 Years with
                            Employer-Sponsored Health Insurance in Nonmetropolitan and Metropolitan Areas — United States, 2009–2014</a>”,
                            Center for Disease Control and Prevention.
                        </p>
                        <p className="small"><sup>2 </sup>
                            The BRCA Exchange maintains compliance with the best bioethics standards available
                            by collaborating with experts in international regulatory and privacy laws.
                            This work is facilitated by the Global Alliance for Genomics and Health Regulatory
                            and Ethics Workstream.
                        </p>
                        <p className="small"><sup>3 </sup>
                            For more information, please visit our <Link to={`/fundraisingdetails`}>Fundraising Details Page</Link>.

                        </p>
                    </Col>
                </Row>
            </Grid>
        );
    }
}

export default WhyDonate;
