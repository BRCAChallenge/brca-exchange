def index():
    return dict()

def data():
    q_source = (db.brca_variant.Variant_Source.upper().contains([request.vars.source.upper()]))
    order_by = getattr(db.brca_variant, request.vars.order_by)
    page_size = int(request.vars.page_size)
    page_num = int(request.vars.page_num)
    limit_by = (0 + page_size * (page_num -1), page_size * page_num + 1)
    search_term = request.vars.search_term.upper()
    search_columns = request.vars.search_column.split(",")
    # this part saves all search queries in various provided columns in a search query string
    # which then get converted to a python expression through eval()
    search_query_dict = {}
    for i in range(len(search_columns)):
        query = (getattr(db.brca_variant, search_columns[i]).upper().contains([search_term]))
        search_query_dict["search_q{0}".format(i)]=query
    search_query_string = ""
    for key in search_query_dict.keys():
        search_query_string = search_query_string + "search_query_dict[\"{0}\"]".format(str(key)) + "|"
    search_query_string = eval(search_query_string[0:-1])

    brca_data = db(search_query_string & q_source).select(orderby=order_by, limitby=limit_by)
    return dict(brca_data=brca_data)

def user():
    """
    exposes:
    http://..../[app]/default/user/login
    http://..../[app]/default/user/logout
    http://..../[app]/default/user/register
    http://..../[app]/default/user/profile
    http://..../[app]/default/user/retrieve_password
    http://..../[app]/default/user/change_password
    http://..../[app]/default/user/manage_users (requires membership in
    http://..../[app]/default/user/bulk_register
    use @auth.requires_login()
        @auth.requires_membership('group name')
        @auth.requires_permission('read','table name',record_id)
    to decorate functions that need access control
    """
    return dict(form=auth())


@cache.action()
def download():
    """
    allows downloading of uploaded files
    http://..../[app]/default/download/[filename]
    """
    return response.download(request, db)


def call():
    """
    exposes services. for example:
    http://..../[app]/default/call/jsonrpc
    decorate with @services.jsonrpc the functions to expose
    supports xml, json, xmlrpc, jsonrpc, amfrpc, rss, csv
    """
    return service()


