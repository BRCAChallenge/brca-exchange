import _ from 'underscore';

const PureRenderMixin = {
    shouldComponentUpdate: function (nextProps, nextState) {
        return !_.isEqual(this.props, nextProps) ||
            !_.isEqual(this.state, nextState);
    }
};

// deep equality check for pure rendering
export default PureRenderMixin;
