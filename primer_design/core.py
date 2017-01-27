"""
  General low level core functions that should be living in a sep module

"""

import sys
import shlex
import subprocess
import re

sys.path.append("/software/lib/python2.7/site-packages/pysam-0.7.5-py2.7-linux-x86_64.egg")
import pysam


verbose_level    =  20

tmp_files = []

import config


if ( __name__ == '__main__'):
    sys.stderr.write( "This is a module and not meant to be run as a stand alone program\n" )
    exit( 1 )




def delete_tmp_files():
    for tmp_filename in tmp_files:
        verbose_print( "deleting tmp file: %s " % tmp_filename, 2)
#        os.remove( tmp_filename )


def add_new_tmp_file( new_file ):
    global tmp_files

    tmp_files.append( new_file )



def set_verbose_level( new_level ):
    
    assert new_level.isdigit(), "Submitted verbose level is not a integer"

    global verbose_level

    verbose_level = new_level



def verbose_print( msg, level ):
    """
    Printing function mainly for development where the developer can
    set the level of debug info s/he wants
    """

    if ( level <= verbose_level ):
        sys.stderr.write( msg+ "\n" )



def revDNA( string ):
    """

    reverses a DNA string. No error catching for non-valid bases or mixed bases.
 
    input: DNA string
    output: reverse DNA string


    Kim Brugger (17 Aug 2016)
    """

    verbose_print( "revDNA", 4)
    rev_bases = { 'A': 'T', 
                  'a': 'T', 
                  'C': 'G',
                  'c': 'G',
                  'G': 'C',
                  'g': 'C',
                  'T': 'A',
                  'T': 'A',
                  '-': '-'}

    rev_str = len(string)*[None]
    for i in range(0, len(string)):
        rev_str[ len(string) - i - 1] =  rev_bases[ string[ i ]]

    return "".join( rev_str )


def align_primers_to_seq( seq, all_primers):
    """

    Find where the primers aligns in a sequence. In this case the region of interest.
    All primers are matched on both strands.

    
    Kim Brugger (17 Aug 2016)
    
    """

    verbose_print( "align_primers_to_seq", 2)

    primers = []
    rev_primers = []
    for primer_set in all_primers:
        ( name, primer) = primer_set
#        print primer
        primers.append(  [name, primer] )
        rev_primers.append( [name, revDNA( primer )])

    mappings = []

    for i in range(0, len(seq)):
        for primer_set in primers:
            (name, primer) = primer_set
            primer_len = len(primer)
            if ("".join(seq[i:i+primer_len]) == primer ):
                mappings.append( [name, primer, i, 0] )

        for primer_set in rev_primers:
            (name, primer) = primer_set
            primer_len = len(primer)
            if ("".join(seq[i:i+primer_len]) == primer ):
                mappings.append( [name, primer, i, 1] )

    return mappings



def fetch_region( chrom, start, end ):
    """
    Extracts a DNA region from the reference genome
    
    input: chromosome (string)
           start position ( int ) 
           end position ( int )

    output: sequence on one line.

    Kim Brugger (25 Jan 2017)
    """

    verbose_print( "fetch_region", 2)

    cmd = "%s faidx %s  %s:%d-%d " % ( config.SAMTOOLS, config.REFERENCE, chrom, start, end )
    args = shlex.split( cmd )
    p = subprocess.Popen(args, stdout=subprocess.PIPE)

    output = p.communicate()
    sequence = ""
    for line in ( output[0].split("\n")):
        if ( re.match('>', line)):
            continue 

        sequence += line
    return sequence 

def fetch_known_SNPs( tabix_file, chrom, start, end):
    """
    Extracts known SNPs from a dbsnp file

    input: location to dbSNP vcf file, 
           chromosom, 
           start pos, 
           end pos 
    
    output: list of list of variants.
    
    Kim Brugger (17 Aug 2016)
    """
    verbose_print( "fetch_known_SNPs", 2)

    cmd = "%s %s  %s:%d-%d " % ( config.TABIX, tabix_file, chrom, start, end )

    verbose_print( cmd, 3)

    args = shlex.split( cmd )
    p = subprocess.Popen(args, stdout=subprocess.PIPE)

    output = p.communicate()
    vars = []

    for line in ( output[0].split("\n")):
        vars.append( line.split("\t"));

    return vars 
