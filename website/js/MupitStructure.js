'use strict';

const React = require('react');
const content = require('./content');
const _ = require('underscore');
const util = require('./util');

export default class MupitStructure extends React.Component {
    constructor(props) {
        super(props);
    }
    render() {
        const {variant, prop} = this.props;
        const aminoAcidCode = util.getAminoAcidCode(variant.HGVS_Protein);
        const mupitStructure = _.find(content.mupitStructures, function(structure) {return structure.name === variant[prop].name;});
        const mupitUrl = mupitStructure.url + "&gm=chr" + variant.Chr + ":" + variant.Pos + "&altaa=" + aminoAcidCode;

        return (
            <div className="mupit-structure">
                <h5>{variant.Genomic_Coordinate_hg38} in the 3D structure of {mupitStructure.humanReadableName}</h5>
                <a href={mupitUrl} target="_blank">
                    <img src={mupitStructure.image} onLoad={this.props.onLoad} />
                </a>
                <p>Click on this thumbnail to open MuPIT viewer from the CRAVAT project, and view the variant in a 3D protein structure context.</p>
            </div>
        );
    }
}
