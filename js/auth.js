var config = require('./config');
var $ = require('jquery');

module.exports = {
    login: function (username, pass, cb) {
        this.logout();
        // Todo: handle the expiration of the tokens
        //if (localStorage.token) {
        //    if (cb) cb(true);
        //    return
        //}
        this.getToken(username, pass, (res) => {
            if (res.authenticated) {
                localStorage.token = res.token;
                if (cb) cb(true)
            } else {
                if (cb) cb(false)
            }
        })
    },

    logout: function () {
        delete localStorage.token;
    },

    loggedIn: function () {
        return !!localStorage.token
    },

    token: function() {
        return localStorage.token
    },

    getToken: function (username, password, cb) {
        var url = config.backend_url + '/accounts/token-auth/';
        $.ajax({
            url: url,
            data: {
                email: username,
                password: password
            },
            dataType: 'json',
            crossDomain: true,
            method: 'POST',
            success: function (res) {
                cb({
                    authenticated: true,
                    token: res.token
                })
            },
            error: function () {
                cb({
                    authenticated: false
                })
            }
        });
    }
};
