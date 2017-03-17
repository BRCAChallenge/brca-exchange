#!/usr/bin/env python

# SOAP client for the Mutalyzer web service in Python using the
# suds library.  This is based on the sample script client-suds.py,
# but rather than calling checkSyntax, calls runMutalyzer and prints
# out the text of the reference and alt allele strings respectively,
# with flanking bases.
#
# See https://www.mutalyzer.nl/webservices
#
# Usage:
#   python ref-alt-string-from-mutalyzer.py 'NM_002001.2:c.1del'
#

import sys
from suds.client import Client

URL = 'https://mutalyzer.nl/services/?wsdl'

if len(sys.argv) < 2:
    print 'Please provide a variant'
    sys.exit(1)

c = Client(URL, cache=None)
o = c.service

print 'Checking ' + sys.argv[1] + ' ...'

#r = o.checkSyntax(sys.argv[1])
mut = o.runMutalyzer(sys.argv[1])
rawVariants = mut.rawVariants
assert(len(rawVariants) == 1)
assert(len(rawVariants[0]) == 1)
variantStrings = rawVariants[0][0].visualisation.split("\n")
altString = variantStrings[0]
refString = variantStrings[1]
print "Ref allele string:", refString
print "Alt allele string:", altString

if r.messages:
    for m in r.messages.SoapMessage:
        print 'Message (%s): %s' % (m.errorcode, m.message)
