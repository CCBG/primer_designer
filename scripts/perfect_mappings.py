#!/usr/bin/python
# 
# 
# 
# 
# Kim Brugger (15 Sep 2015), contact: kim@brugger.dk

import sys
import os
import pprint
pp = pprint.PrettyPrinter(indent=4)

import shlex
import subprocess
import re

# This is a bit hacky, but it is an easy way to point to the local primer_design modules
path = os.path.dirname(os.path.realpath(__file__))
sys.path.append( path + "/../")

import primer_design.core
import primer_design.config
import primer_design.primer3
import primer_design.blast
import primer_design.analysis



def get_and_parse_arguments():
    primer_design.core.verbose_print( "get_and_parse_arguments", 4)

    import argparse

    parser = argparse.ArgumentParser()

    group = parser.add_mutually_exclusive_group( required = True )
    group.add_argument('-f', '--fasta' )

    args = parser.parse_args()


    return args







#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
#    Main loop !
#
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def main():

    args = get_and_parse_arguments()

#    print "Hello world! \n"

    mappings = primer_design.blast.map_sequences( args.fasta )

#    pp.pprint( mappings )

    for seq in mappings:
        if ( len( mappings[ seq ][ 'CHR' ]  ) > 1 ):
            continue

        print "%s\t%s\t%d\t%s" % ( seq, 
                                   mappings[ seq ][ 'CHR' ][0],
                                   mappings[ seq ][ 'POS' ][0],
                                   mappings[ seq ][ 'STRAND' ][0])



    primer_design.core.delete_tmp_files()


if __name__ == '__main__':
    main()


