'use strict';

var React = require('react');
var content = require('./content');
var _ = require('underscore');
var util = require('./util');


var mupitStructure = function(variant, prop) {
    return <MupitStructure variant={variant} prop={prop} />;
};

class MupitStructure extends React.Component {
    constructor(props) {
        super(props);
    }
    render() {
        let {variant, prop} = this.props;
        let aminoAcidCode = util.getAminoAcidCode(variant.HGVS_Protein);
        let mupitStructure = _.find(content.mupitStructures, function(structure) {return structure.name === variant[prop].name;});
        let mupitUrl = mupitStructure.url + "&gm=chr" + variant.Chr + ":" + variant.Pos + "&altaa=" + aminoAcidCode;
        return (
            <div className="mupit-structure">
                <h5>{variant.Genomic_Coordinate_hg38} in the 3D structure of {mupitStructure.humanReadableName}</h5>
                <a href={mupitUrl} target="_blank">
                    <img src={mupitStructure.image}></img>
                </a>
                <p>Click on this thumbnail to open MuPIT viewer from the CRAVAT project, and view the variant in a 3D protein structure context.</p>
            </div>
        );
    }
}

module.exports = mupitStructure;
