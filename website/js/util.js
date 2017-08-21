'use strict';

module.exports.isEmptyField = function isEmptyField(value) {
    if (Array.isArray(value)) {
        value = value[0];
    }

    if (value === null || (typeof value === 'undefined')) {
        return true;
    }

    var v = value.trim();
    return v === '' || v === '-' || v === 'None';
};

