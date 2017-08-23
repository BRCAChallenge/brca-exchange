#!/usr/bin/env python2.7
# Melissa Meredith (meredith705)

# Runs Mutalyzer and returns a Mutalyzerfile with binary flags for intronic and 'out of range' variants
# In the event of truncation it returns the locaiton of a stop codon

import sys
import re
import ast
from suds.client import Client

# should be called as an argument when running on sh script
#variant = "NM_007294.3:c.1112delC"

#def mutalyzer(self, file, HGVS_list):
def mutalyzer(self, file, HGVS_list):
	#needed to run mutalyzer
	URL = 'https://mutalyzer.nl/services/?wsdl'
	c = Client(URL, cache=None)
	o = c.service

### runs Mutalyzer creating many fields: mutalyzer.nl/soap-api#op.runMutalyzer
	mm = o.runMutalyzer(variant)

	#prints message if varaint is INTRONIC or "out of range"
	# outputs a binary flag for each situation (0/1)
	if mm.messages:
	    for m in mm.messages.SoapMessage:
#	        print 'Error %s: %s\n' % (m.errorcode, m.message)      
	        if (m.errorcode == "ENOINTRON"):
	        	intronic = 1
	        	inTranscript = 0
	        	NewstopCodon = 0
	        elif (m.errorcode == "ERANGE"):
	        	inTranscript = 0
	        	intronic = 0
	        	NewstopCodon = 0
	elif mm.newProtein:
	# if no error messages then variant is in transcript
		intronic = 0
		inTranscript = 1

		#if mm.newProtein == "mm.newProtein":
		if len(mm.newProtein) < len(mm.origProtein) and len(mm.newCDS) < len(mm.origCDS) and len(mm.mutatedMRNA) < len(mm.origMRNA): 
		    orig_stop = mm.origProtein.index('*')
		    new_stop = mm.newProtein.index('*')
		    if orig_stop != new_stop:
		    	NewstopCodon = new_stop
		else:
			NewstopCodon = 0

	
	
	file.write("%d\t%d\t%d\n" % (intronic, inTranscript, NewstopCodon))
	print intronic, "\t", inTranscript, "\t",  NewstopCodon
	return(intronic, inTranscript, NewstopCodon)
	

# HGVS nomenclature will be passed in here by BrcaParse.py, for now it's a list
HGVS_list = ['NM_007294.3:c.5467+8G>A', 'NM_007294.3:c.2783G>A', 'NM_007294.3:c.1931G>T', 'NM_007294.3:c.441+5A>G', 'NM_007294.3:c.-2943G>A', 'NM_007294.3:c.134+2161G>A', 'NM_007294.3:c.1508delA']

#print "intron\tinTscrip\tNewStop" 
file = open('Mutalyzerfile', 'w')
file.write("intronic\tinTranscript\tNewstopCodon\n")

for variant in HGVS_list:
	mutalyzer(variant, file, HGVS_list)
file.close()

	
