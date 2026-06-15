'use strict';

import React from 'react';
import { Container as Grid, Row, Alert } from 'react-bootstrap';
import { withRouter } from 'react-router-dom';
import $ from 'jquery';
import config from './config';

class ConfirmEmailInner extends React.Component {
  state = { success: null, loading: true };

  activate = (result) => {
    this.setState({ success: !!(result && result.success), loading: false });
  };

  componentDidMount() {
    const { activationCode } = (this.props.match && this.props.match.params) || {};
    const url = `${config.backend_url}/accounts/confirm/${activationCode}/`;
    $.get(url, this.activate);
  }

  render() {
    if (this.state.loading) return <Grid id="main-grid"><Row>Loading...</Row></Grid>;
    return (
      <Grid id="main-grid">
        {this.state.success ? (
          <Row id="message">
            <Alert variant="success">Thanks for confirming your email.</Alert>
          </Row>
        ) : (
          <Row id="message">
            <Alert variant="danger">Sorry, this activation link is not valid.</Alert>
          </Row>
        )}
      </Grid>
    );
  }
}

const ConfirmEmail = withRouter(ConfirmEmailInner);
export { ConfirmEmail };
export default ConfirmEmail;
