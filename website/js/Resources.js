'use strict';

import React from 'react';
import RawHTML from './RawHTML';
import { Container as Grid, Row, Col} from 'react-bootstrap';
import { resources } from './content';

class Resources extends React.Component {
    render() {
	return (
            <Grid id="main-grid" className="help-page">
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
                                    <h4>{tile.name}</h4>
                                    <RawHTML html={tile.contents} />
                                </Col>
                            </Row>
                        ))}
                    </div>
                ))}
            </Grid>
        );
    }
}

export default Resources;
