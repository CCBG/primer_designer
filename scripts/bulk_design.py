#!/usr/bin/python
# 
# 
# 
# 
# Kim Brugger (11 Jul 2016), contact: kim@brugger.dk

import sys
import pprint
pp = pprint.PrettyPrinter(indent=4)

import os
import sys
import re
#import pipeliners as pipe

import subprocess
def system_call( cmd ):

    try:
        subprocess.check_call(cmd, shell=True)

    except subprocess.CalledProcessError as scall:
        print "System call '%s' failed. Exit code: %s Error output: %s" %( cmd,  scall.returncode, scall.output)
        exit()



# Define system call

if ( len(sys.argv ) == 1 ):
	print "USAGE: %s input-file"
	print "input-file one tab seperated output-name and region of interest pr line"
	exit()

infile = sys.argv[1]
if ( infile.find(".txt") == -1):
	print "Infile '%s' does not look like a txt file" % infile
	exit()

outfile = re.sub(r'.txt', '.zip', infile)

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
		if ( line.find("\t") == -1):
			print "'%s' does not contain a 'tab'";
			next

		(testname, region) = line.split("\t")
		if ( region.find(":") == -1 ):
			print "'%s' does not contain a ':'";
			next

		(chrom, npos) = region.split(":")
	        cmd = "/software/packages/primer_designer/primer_designer_region.py -c " + chrom + " -p " + npos + " -o " + testname
#		print cmd
		system_call(cmd )


system_call( "zip -j %s *pdf" % outfile )

print "Output file: %s" % outfile
