import re
import os
import sys
import pipeliners as pipe

# Define system call

infile = sys.argv[1]

with open( infile , "rU" ) as plist:

	for line in plist:
		line = line.strip("\n")

		testname = re.sub(r"(.*)\t.*", r"\1", line)
		chro = re.sub(r".*\t(.*):.*",r"\1",line)
		npos = re.sub(r".*\t.*:(.*)",r"\1",line)
                (testname, chro_pos) = line.split('\t')
                (chro, npos) = chro_pos.split(':')
	        cmd = "primer_designer -c " + chro + " -p " + npos + " -o " + testname
		print cmd
#		print pipe.system_call("Running design_primers", cmd )
