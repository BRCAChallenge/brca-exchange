/*eslint-env browser */
'use strict';

const React = require('react');
const {Grid, Row, Col} = require('react-bootstrap');
const {Link} = require('react-router');

const DonateDrive = React.createClass({

    render() {
        return (
            <Grid id="main-grid" className="help-page">
                <Row>
                    <Col smOffset={1} sm={10}>
                        <img id='ngbocw-logo' className="centered margin-top-forty" src={require('./img/national_hereditary_cancer_week.jpg')} />
                        <h3 className="margin-top-forty margin-bottom-thirty">National Hereditary Breast and Ovarian Cancer Week is an opportunity to help in the fight against these hereditary cancers. Help us make a difference!</h3>
                        <p>
                            <a href="https://www.cancer.net/research-and-advocacy/cancer-awareness-dates" target="_blank">National HBOC Week</a> marks the transition between Ovarian Cancer
                            Awareness Month in September and Breast Cancer Awareness Month in October, recognizing the link between these cancers. An estimated 10% of breast cancers and 20% of ovarian cancers arise from heritable genetic risk,
                            most often from variation in the BRCA1 and BRCA2 genes; this variation also increases the risk of additional cancers including prostate,
                            pancreatic and male breast cancer<sup><a href="https://www.cancer.net/cancer-types/hereditary-breast-and-ovarian-cancer" target="_blank">1</a></sup>.  Knowing your risk can help you avoid these cancers.  A growing number of “previvors” were born
                            predisposed to these cancers but have not had the disease, thanks in many cases to increased screening or risk-reducing surgeries.
                            We celebrate the courage of these individuals on September 30, National Previvor Day.
                        </p>
                        <p>
                            We at BRCA Exchange are passionate about helping individuals better understand and manage their heritable risk of cancer.
                            In partnership with the ENIGMA Consortium, the international authority on BRCA variant interpretation, we are gathering the evidence
                            that scientists need to evaluate how genetic variation affects individuals’ disease risk.  Here are ways you can help.
                        </p>
                        <ol className="margin-ordered-list">
                            <li><Link to={`/whydonate`}>Donate!</Link>  Every dollar raised will help two or more families better understand their cancer risk.  We suggest donations of $40, to celebrate the 40K variants now available on BRCA Exchange.</li>
                            <li><Link to={`/community`}>Join our community!</Link>  By joining the BRCA Exchange community, you will indicate your support for this work, which will help us raise funds in the future.</li>
                            <li>Get involved!  By partnering with advocacy groups such as <a href="https://www.facingourrisk.org/index.php" target="_blank">FORCE</a> or <a href="https://lightcollective.org/" target="_blank">The Light Collective</a>, you can help push for changes in access to care, privacy, and many other areas that benefit cancer patients and previvors everywhere.</li>
                        </ol>
                    </Col>
                </Row>
            </Grid>
        );
    }
});

module.exports = DonateDrive;
