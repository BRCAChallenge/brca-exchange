// eslint.config.js
const js = require("@eslint/js");
const react = require("eslint-plugin-react");
const babelParser = require("@babel/eslint-parser");
const globals = require("globals");

module.exports = [
  js.configs.recommended,
  react.configs.flat.recommended,
  react.configs.flat["jsx-runtime"],
  {
    files: ["js/**/*.js"],
    languageOptions: {
      ...react.configs.flat.recommended.languageOptions,
      parser: babelParser,
      ecmaVersion: "latest",
      sourceType: "module",
      parserOptions: {
        requireConfigFile: false,
        babelOptions: {
          presets: [
            ["@babel/preset-env", { modules: false }],
            ["@babel/preset-react", { runtime: "automatic" }],
          ],
        },
      },
      globals: {
        ...globals.browser,
        ...globals.node,
        ...globals.es6,
        google: false,
        grecaptcha: false,
      },
    },
    rules: {
      "camelcase": 2,
      "dot-notation": [2, { "allowKeywords": false, "allowPattern": "^[a-zA-Z]+(_[a-zA-Z]+)+$" }],
      "eol-last": 2,
      "eqeqeq": [2, "allow-null"],
      "curly": ["error", "multi-line"],
      "key-spacing": 2,
      "new-cap": [2, { "capIsNew": false }],
      "no-dupe-keys": 2,
      "no-extra-bind": 2,
      "no-extra-boolean-cast": 2,
      "no-trailing-spaces": 2,
      "no-use-before-define": 2,
      "space-infix-ops": 2,
      "space-before-blocks": 2,
      "yoda": 0,
      "no-extra-parens": [2, "functions"],
      "comma-spacing": 2,
      "no-undef": 2,
      "no-redeclare": [2, { "builtinGlobals": true }],
      "no-underscore-dangle": 0,
      "no-mixed-spaces-and-tabs": 0,
      "no-multi-spaces": 0,
      "no-console": 0,
      "quotes": 0,
      "no-shadow": 0,
      "no-unused-vars": [2, { "vars": "local", "ignoreRestSiblings": true, "args": "after-used", "argsIgnorePattern": "^_" }],
      "react/jsx-no-undef": 1,
      "react/no-did-mount-set-state": 1,
      "react/no-did-update-set-state": 1,
      "react/no-unknown-property": 1,
      "react/self-closing-comp": 1,
      "react/jsx-wrap-multilines": 1,
      "react/prop-types": 0,
      "semi": [2, "always"],
    },
  },
  {
    files: ["test/**/*.js"],
    languageOptions: {
      ...react.configs.flat.recommended.languageOptions,
      parser: babelParser,
      ecmaVersion: "latest",
      sourceType: "module",
      parserOptions: {
        requireConfigFile: false,
        babelOptions: {
          presets: [
            ["@babel/preset-env", { modules: false }],
            ["@babel/preset-react", { runtime: "automatic" }],
          ],
        },
      },
      globals: {
        ...globals.browser,
        ...globals.node,
        ...globals.es6,
        ...globals.mocha,
      },
    },
    rules: {
      "eqeqeq": [2, "allow-null"],
      "curly": ["error", "multi-line"],
      "no-undef": 2,
      "no-redeclare": [2, { "builtinGlobals": false }],
      "no-unused-vars": [2, { "vars": "local", "ignoreRestSiblings": true, "args": "after-used", "argsIgnorePattern": "^_" }],
      "semi": [2, "always"],
    },
  },
];
