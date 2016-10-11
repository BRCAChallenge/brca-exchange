'use strict';

var months = ["January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"];
var dateFormat = (date) => { let d = new Date(date); return `${months[d.getMonth()]} ${d.getDate()}, ${d.getFullYear()}`; };

module.exports = {
    dateFormat
};
