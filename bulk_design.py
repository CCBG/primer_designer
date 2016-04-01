#!/usr/bin/python
import re
import os
import sys
import pipeliners as pipe

infile = sys.argv[1]

# Run designer function 
def run_designer(testname, chro, npos):
	
	if "-" in npos:
		
		npos = re.sub(r"(.*)-(.*)", r"\1 \2", npos)
		#print npos
		cmd = "primer_designer_region -c " + chro + " -r " + npos + " -o " + testname
		#print "region"
		pipe.system_call("Running design_primers", cmd )
	
	else:
		cmd = "primer_designer_region -c " + chro + " -p " + npos + " -o " + testname
		#print "single base"
		pipe.system_call("Running design_primers", cmd )

# Main Loop
with open( infile , "rU" ) as plist:
	
	for line in plist:
		
		line = line.strip("\n")
		line = line.strip(" ")		
		#print line
		
		testname = re.sub(r"^(.*?)\t.*", r"\1", line)
		chro     = re.sub(r".*\t(.*):.*",r"\1",line)
		npos     = re.sub(r".*\t.*:(.*?)",r"\1",line)
		
		#print testname
		#print chro		
		#print npos

		run_designer(testname, chro, npos)
