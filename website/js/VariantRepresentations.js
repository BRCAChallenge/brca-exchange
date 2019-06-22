'use strict';

var React = require('react');
var backend = require('./backend');


var VariantRepresentations = React.createClass({
    getInitialState: () => ({ data: {} }),
    componentWillMount: function() {
        console.log('here');
        backend.variantReps().subscribe(
            resp => this.setState(resp),
            () => this.setState({error: 'Problem connecting to server'}));
    },
    render: function () {
        // Ensure releases are in descending order
        let representations = JSON.stringify(this.state.data);
        return (
            <div>
                <h1 style={{ textAlign: 'center'}}>Variant Representations</h1>
                <div style= {{ padding: '20px'}}>{representations}</div>
            </div>
        );
    }
});

module.exports = ({
    VariantRepresentations: VariantRepresentations,
});
