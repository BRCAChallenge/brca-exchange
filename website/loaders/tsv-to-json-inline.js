// loaders/tsv-to-json-inline.js
// Transforms TSV → compact JSON string so asset/resource will emit .json.
// Uses dsv-loader underneath to parse TSV (add it to your deps).
const dsv = require('d3-dsv'); // or 'd3-dsv'
const { getOptions } = require('loader-utils');

module.exports = function tsvToJsonInline(source) {
  // If you previously used key[]=... selection in the query, you can parse it here:
  const opts = getOptions(this) || {};
  // Parse TSV
  const rows = dsv.tsvParse(source);
  // TODO: If you need to enforce uniqueness on databaseKey, do it here.
  const json = JSON.stringify(rows);
  // Return the raw JSON string so asset/resource writes it out as a file.
  return json;
};

// Tell webpack this loader expects raw text input
module.exports.raw = true;

