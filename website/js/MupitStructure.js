'use strict';

var React = require('react');
var content = require('./content');
var _ = require('underscore');


var mupitStructure = function(variant, prop) {
    return <MupitStructure variant={variant} prop={prop} />;
};

class MupitStructure extends React.Component {
    constructor(props) {
        super(props);
    }
    render() {
        let {variant, prop} = this.props;
        let mupitStructure = _.find(content.mupitStructures, function(structure) {return structure.name === variant[prop].name;});
        debugger;
        let mupitUrl = mupitStructure.url + "&gm=chr" + variant.Chr + ":" + variant.Pos + "&altaa=" + variant.Alt.toLowerCase();
        return (
            <div className="mupit-structure">
                <a href={mupitUrl} target="_blank">
                    <img src={mupitStructure.image}></img>
                </a>
            </div>
        );
    }
}

module.exports = mupitStructure;
