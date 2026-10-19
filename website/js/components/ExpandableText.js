'use strict';

import React from "react";
import { splitPartsIntoSentences, splitPartsIntoWords } from '../util';

// ExpandableText truncates long content (e.g. a multi-sentence free-text
// field) down to a short preview, with a "Click for More" affordance that
// expands to the full content, and a "Click for Less" affordance at the end
// of the expanded content that collapses it back down.
//
// `content` may be a plain string, a single React node, or an array mixing
// strings and nodes (e.g. the array getFormattedFieldByProp() produces for
// fields with embedded PMID links) -- in all cases it's rendered unchanged
// when it's short enough that truncation isn't needed.
//
// This is intentionally generic so it can be reused for any tile text field
// that may run long (e.g. 'Comment on Clinical Significance' in the ENIGMA
// tile, or 'Summary Evidence' in the ClinVar tile), not just the ones it's
// first applied to.
export default class ExpandableText extends React.PureComponent {
    state = {
        expanded: !!this.props.defaultExpanded
    };

    toggleExpanded = (event) => {
        if (event) {
            event.preventDefault();
        }
        this.setState((prevState) => ({ expanded: !prevState.expanded }), () => {
            if (this.props.relayoutGrid) {
                this.props.relayoutGrid();
            }
        });
    };

    renderToggle(label) {
        return (
            <span
                role="button"
                tabIndex={0}
                className="expandable-text-toggle"
                onClick={this.toggleExpanded}
                onKeyDown={(e) => {
                    if (e.key === 'Enter' || e.key === ' ') {
                        this.toggleExpanded(e);
                    }
                }}
            >
                {label}
            </span>
        );
    }

    render() {
        const {
            content,
            limit = 3,
            mode = "sentences",
            moreLabel = "Click for More",
            lessLabel = "Click for Less",
            className
        } = this.props;
        const { expanded } = this.state;

        const groups = mode === "words"
            ? splitPartsIntoWords(content)
            : splitPartsIntoSentences(content);

        // nothing to truncate: just render the content as-is
        if (groups.length <= limit) {
            return <span className={className}>{content}</span>;
        }

        if (expanded) {
            return (
                <span className={className}>
                    {groups.map((group, idx) => <React.Fragment key={idx}>{group}</React.Fragment>)}
                    {' '}
                    {this.renderToggle(lessLabel)}
                </span>
            );
        }

        const visibleGroups = groups.slice(0, limit);

        return (
            <span className={className}>
                {visibleGroups.map((group, idx) => <React.Fragment key={idx}>{group}</React.Fragment>)}
                {'\u2026 '}
                {this.renderToggle(moreLabel)}
            </span>
        );
    }
}
