'use strict';

import React from 'react';
import RawHTML from './RawHTML';
import { Container as Grid, Row, Col} from 'react-bootstrap';
import { resources } from './content';

class Resources extends React.Component {
    render() {
	return (
            <Grid id="main-grid" className="help-page" style={{paddingTop: '40px'}}>
                {resources.map((section, sectionIdx) => (
                    <div key={sectionIdx}>
                        <Row>
                            <Col sm={{ span: 10, offset: 1 }}>
                                <h3>{section.section}</h3>
                            </Col>
                        </Row>
                        {section.tiles.map((tile, tileIdx) => (
                            <Row key={tileIdx}>
                                <Col sm={{ span: 10, offset: 1 }}>
                                    <p>{tile.name}</p>
                                    <RawHTML html={tile.contents} />
                                </Col>
                            </Row>
                        ))}
			<Row>
    			    <Col sm={{ span: 10, offset: 1 }}>
        			<hr />
    			    </Col>
			</Row>
			<Row>
			   <Col sm={{ span: 10, offset: 1 }}>
				<p>Don’t see your resource here? <a href="mailto:brcaexchange@gmail.com?subject=BRCA Exchange website">Contact us</a>!</p>
			   </Col>
			</Row>
                    </div>
                ))}
            </Grid>
        );
    }
}

export default Resources;
