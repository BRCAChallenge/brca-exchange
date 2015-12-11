import cStringIO
import traceback

import csv
from itertools import chain, starmap, groupby, imap

# A iterator future. Evaluates an expression returning an iterator when
# the 1st value is requested. This was an attempt to defer execution of
# the query, however it fails because the db is closed after the action
# function returns. Consequently, we still have a delay of several seconds
# before the browser gets a packet back from web2py. This is horrible.
#class Deferred:
#    def __init__(self, thunk):
#        self.thunk = thunk
#        self.iter = None
#
#    def __iter__(self):
#        return self
#
#    def next(self):
#        if self.iter == None:
#            self.iter = self.thunk()
#        return self.iter.next()

class Keyfn():
    def __init__(self, m):
        self.i = -1
        self.m = m
    def __call__(self, x):
        self.i = self.i + 1
        return self.i / self.m

def chunk(i, rows):
    s = cStringIO.StringIO()
    w = csv.writer(s, dialect='excel-tab')
    for r in rows:
        w.writerow(r)
    buf = s.getvalue()
    s.close()
    return buf

def rowToList(colnames):
    def fn(row):
        return map(lambda c: row[c], colnames)
    return fn

# chunked tsv iterator. Not sure if we need this for HTTP,
# but w/o this we'd be creating a csv writer and StringIO for every
# row, which seems crazy. Maybe there's a better csv API?
def totsv(table, cols, data):
    def header():
        yield cols

    keyfn = Keyfn(10)
    toList = rowToList(map(lambda c: table + '.' + c, cols))

    groups = groupby(chain(header(), imap(toList, data)), keyfn)
    return starmap(chunk, groups)

def tolist(x):
    return x if isinstance(x, list) else [x]

def index():
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

    data = db(query).select(orderby=order_by_dir, limitby=limit_by)

    if request.vars.format == 'tsv':
        colnames = [f for f in db.brca_variant.fields]
        response.headers['Content-Type'] = 'application/vnd.ms-excel'
        # Content-Disposition messes up Content-Type
#        response.headers['Content-Disposition'] = 'attachment;filename="variants.tsv"'
#        data = Deferred(lambda: db(query).select(orderby=order_by_dir, limitby=limit_by))
        return totsv('brca_variant', colnames, data)

    response.view = 'default/data.json'
    count = db(query).count()
    return dict(data=data, count=count)
