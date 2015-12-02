def tolist(x):
    return x if isinstance(x, list) else [x]

def index():
    if request.vars.format == 'tsv':
        response.view = 'default/data.tsv'
    else:
        response.view = 'default/data.json'
    response.headers['Access-Control-Allow-Origin'] = '*'

    query = (db.brca_variant.Variant_Source.upper().contains([request.vars.source.upper()]))
    direction = request.vars.direction
    order_by = getattr(db.brca_variant, request.vars.order_by)
    order_by_dir = ~order_by if direction == 'descending' else order_by
    if request.vars.page_size:
        page_size = int(request.vars.page_size)
        page_num = int(request.vars.page_num)
        limit_by = (page_size * page_num, page_size * (page_num + 1))
    else:
        limit_by = None

    if request.vars.search_term:
        search_term = request.vars.search_term.upper()
        search_columns = tolist(request.vars.search_column)
        term = lambda c: getattr(db.brca_variant, c).upper().contains([search_term])

        subquery = term(search_columns[0])
        for c in search_columns[1:]:
            subquery |= term(c)

        query &= subquery

    if request.vars.filter:
        filters = zip(tolist(request.vars.filter), tolist(request.vars.filterValue))
        for (p, v) in filters:
            query &= getattr(db.brca_variant, p) == v

    brca_data = db(query).select(orderby=order_by_dir, limitby=limit_by)
    count = db(query).count()

    return dict(data=brca_data, count=count)
