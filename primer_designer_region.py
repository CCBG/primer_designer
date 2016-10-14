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

from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfbase.pdfmetrics import stringWidth 
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import A4
width, height = A4 #keep for later

sys.path.append("/software/lib/python2.7/site-packages/pysam-0.7.5-py2.7-linux-x86_64.egg")
import pysam


TMP_FILES = []

colours = [[255,   0,   0], # red
           [  0, 255,   0], # green
           [  0,   0, 255], # blue
           [255,   0, 255], # Pink
           [0,   255, 255],
           [255, 255,   0], #Dark Green
           [100, 255, 100]] # Yellow, crap!

# External tools we use 
SAMTOOLS   = '/software/bin/samtools ';
REFERENCE  = '/refs/human_1kg/human_g1k_v37.fasta '
TABIX      = '/software/bin/tabix '
SMALT      = '/software/bin/smalt-0.7.6 '
PRIMER3    = '/software/bin/primer3_core '

# default parameters, can be overridden by the user (perhaps)
TARGET_LEAD        = 50
FLANK              = 500 # amount of bp on either side of the target to design primers in
NR_PRIMERS         = 4   # Nr of primers to report back
ALLOWED_MISMATCHES = 4   # Maximum nr of errors when mapping a primer back to the reference
MAX_MAPPINGS       = 5   # Less or equal number of mappings to the reference what chromosomes mapped to are named 
MAX_PRODUCT        = 500

VERBOSE    =  20
VERSION    =  '2.0-beta1'

def verbose_print( msg, level ):
    if ( level <= VERBOSE ):
        print msg

#
# Extracts a DNA region from the reference genome
#
# input: Chromosome, start position, end position
#
# output: sequence on one line.
#
# Kim Brugger (17 Aug 2016), contact: kbr@brugger.dk
def fetch_region( chrom, start, end ):

    verbose_print( "fetch_region", 2)

    cmd = "%s faidx %s  %s:%d-%d " % ( SAMTOOLS, REFERENCE, chrom, start, end )
    args = shlex.split( cmd )
    p = subprocess.Popen(args, stdout=subprocess.PIPE)

    output = p.communicate()
    sequence = ""
    for line in ( output[0].split("\n")):
        if ( re.match('>', line)):
            continue 

        sequence += line
    return sequence 

#
# Extracts known SNPs from a dbsnp file
#
# input: dbSNP-file, chromosom, start pos, end pos 
#
# output: list of list of variants.
#
# Kim Brugger (17 Aug 2016)
def fetch_known_SNPs( tabix_file, chrom, start, end):
    verbose_print( "fetch_known_SNPs", 2)

    cmd = "%s %s  %s:%d-%d " % ( TABIX, tabix_file, chrom, start, end )

    verbose_print( cmd, 3)

    args = shlex.split( cmd )
    p = subprocess.Popen(args, stdout=subprocess.PIPE)

    output = p.communicate()
    vars = []

    for line in ( output[0].split("\n")):
        vars.append( line.split("\t"));

    return vars 


#
# Generate the file required to run primer3
#
# input: sequence-id, sequence, output-file (optional)
#
# output: None
#
# Kim Brugger (17 Aug 2016)
def write_primer3_file(seq_id, seq, primer3_file = ""):
    verbose_print( "write_primer3_file", 2)

    upstream_end      = 1
    downstream_start =  2
    for i in range(0, len(seq)):
            if seq[ i ] == '[':
                upstream_end = i
            elif seq[ i ] == ']':
                downstream_start = i


    template = '''SEQUENCE_ID={seq_id}
SEQUENCE_TEMPLATE={seq}
SEQUENCE_TARGET={flank},{len}
PRIMER_FIRST_BASE_INDEX=1
PRIMER_TASK=generic
PRIMER_MIN_THREE_PRIME_DISTANCE=3
PRIMER_EXPLAIN_FLAG=1
PRIMER_MAX_LIBRARY_MISPRIMING=12.00
PRIMER_PAIR_MAX_LIBRARY_MISPRIMING=20.00
PRIMER_PRODUCT_SIZE_RANGE=400-700
PRIMER_NUM_RETURN=5
PRIMER_MAX_END_STABILITY=9.0
PRIMER_MAX_SELF_ANY_TH=45.00
PRIMER_MAX_SELF_END_TH=35.00
PRIMER_PAIR_MAX_COMPL_ANY_TH=45.00
PRIMER_PAIR_MAX_COMPL_END_TH=35.00
PRIMER_MAX_HAIRPIN_TH=24.00
PRIMER_MAX_TEMPLATE_MISPRIMING_TH=40.00
PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH=70.00
PRIMER_MIN_SIZE=18
PRIMER_OPT_SIZE=20
PRIMER_MAX_SIZE=25
PRIMER_MIN_TM=58.0
PRIMER_OPT_TM=60.0
PRIMER_MAX_TM=62.0
PRIMER_PAIR_MAX_DIFF_TM=5.0
PRIMER_TM_FORMULA=1
PRIMER_SALT_MONOVALENT=50.0
PRIMER_SALT_CORRECTIONS=1
PRIMER_SALT_DIVALENT=1.5
PRIMER_DNTP_CONC=0.6
PRIMER_DNA_CONC=50.0
PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=1
PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT=0
PRIMER_LOWERCASE_MASKING=0
PRIMER_MIN_GC=30.0
PRIMER_MAX_GC=70.0
PRIMER_MAX_NS_ACCEPTED=0
PRIMER_MAX_POLY_X=4
PRIMER_OUTSIDE_PENALTY=0
PRIMER_GC_CLAMP=0
PRIMER_LIBERAL_BASE=1
PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS=0
PRIMER_PICK_ANYWAY=1
PRIMER_WT_TM_LT=1.0
PRIMER_WT_TM_GT=1.0
PRIMER_WT_SIZE_LT=1.0
PRIMER_WT_SIZE_GT=1.0
PRIMER_WT_GC_PERCENT_LT=0.0
PRIMER_WT_GC_PERCENT_GT=0.0
PRIMER_WT_SELF_ANY_TH=0.0
PRIMER_WT_SELF_END_TH=0.0
PRIMER_WT_HAIRPIN_TH=0.0
PRIMER_WT_NUM_NS=0.0
PRIMER_WT_LIBRARY_MISPRIMING=0.0
PRIMER_WT_SEQ_QUAL=0.0
PRIMER_WT_END_QUAL=0.0
PRIMER_WT_POS_PENALTY=0.0
PRIMER_WT_END_STABILITY=0.0
PRIMER_WT_TEMPLATE_MISPRIMING_TH=0.0
PRIMER_PAIR_WT_PRODUCT_SIZE_LT=0.0
PRIMER_PAIR_WT_PRODUCT_SIZE_GT=0.0
PRIMER_PAIR_WT_PRODUCT_TM_LT=0.0
PRIMER_PAIR_WT_PRODUCT_TM_GT=0.0
PRIMER_PAIR_WT_DIFF_TM=0.0
PRIMER_PAIR_WT_COMPL_ANY_TH=0.0
PRIMER_PAIR_WT_COMPL_END_TH=0.0
PRIMER_PAIR_WT_LIBRARY_MISPRIMING=0.0
PRIMER_PAIR_WT_PR_PENALTY=1.0
PRIMER_PAIR_WT_IO_PENALTY=0.0
PRIMER_PAIR_WT_TEMPLATE_MISPRIMING=0.0
PRIMER_INTERNAL_WT_SIZE_LT=1.0
PRIMER_INTERNAL_WT_END_QUAL=0.0
PRIMER_INTERNAL_MAX_SELF_END=12.00
PRIMER_QUALITY_RANGE_MIN=0
PRIMER_PAIR_MAX_COMPL_END=3.00
PRIMER_PRODUCT_MAX_TM=1000000.0
PRIMER_INTERNAL_MAX_SIZE=27
PRIMER_INTERNAL_WT_SELF_ANY=0.0
PRIMER_INTERNAL_MAX_POLY_X=5
PRIMER_INTERNAL_WT_SIZE_GT=1.0
PRIMER_SEQUENCING_ACCURACY=20
PRIMER_INTERNAL_WT_TM_GT=1.0
PRIMER_INTERNAL_WT_LIBRARY_MISHYB=0.0
PRIMER_INTERNAL_MAX_GC=80.0
PRIMER_PAIR_WT_COMPL_ANY=0.0
PRIMER_PICK_INTERNAL_OLIGO=0
PRIMER_MAX_SELF_END=3.00
PRIMER_QUALITY_RANGE_MAX=100
PRIMER_INTERNAL_DNTP_CONC=0.0
PRIMER_INTERNAL_MIN_SIZE=18
PRIMER_INTERNAL_MIN_QUALITY=0
PRIMER_SEQUENCING_INTERVAL=250
PRIMER_INTERNAL_SALT_DIVALENT=1.5
PRIMER_MAX_SELF_ANY=8.00
PRIMER_INTERNAL_WT_SEQ_QUAL=0.0
PRIMER_PAIR_WT_COMPL_END=0.0
PRIMER_INTERNAL_OPT_TM=60.0
PRIMER_SEQUENCING_SPACING=500
PRIMER_INTERNAL_MAX_SELF_ANY=12.00
PRIMER_MIN_END_QUALITY=0
PRIMER_INTERNAL_MIN_TM=57.0
PRIMER_PAIR_MAX_COMPL_ANY=8.00
PRIMER_SEQUENCING_LEAD=50
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_INTERNAL_OPT_SIZE=20
PRIMER_WT_TEMPLATE_MISPRIMING=0.0
PRIMER_MAX_END_GC=5
PRIMER_MIN_QUALITY=0
PRIMER_INTERNAL_MAX_LIBRARY_MISHYB=12.00
PRIMER_INTERNAL_WT_GC_PERCENT_GT=0.0
PRIMER_INTERNAL_MAX_NS_ACCEPTED=0
PRIMER_WT_SELF_ANY=0.0
PRIMER_MAX_TEMPLATE_MISPRIMING=12.00
PRIMER_INTERNAL_WT_NUM_NS=0.0
PRIMER_INTERNAL_WT_SELF_END=0.0
PRIMER_PRODUCT_OPT_SIZE=0
PRIMER_PRODUCT_OPT_TM=0.0
PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING=24.00
PRIMER_INSIDE_PENALTY=-1.0
PRIMER_INTERNAL_MIN_GC=20.0
PRIMER_PRODUCT_MIN_TM=-1000000.0
PRIMER_INTERNAL_SALT_MONOVALENT=50.0
PRIMER_WT_SELF_END=0.0
PRIMER_INTERNAL_DNA_CONC=50.0
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_INTERNAL_MAX_TM=63.0
PRIMER_INTERNAL_WT_GC_PERCENT_LT=0.0
PRIMER_INTERNAL_WT_TM_LT=1.0
PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/software/bin/primer3_config/
='''

    details = {
        "seq_id" : seq_id,
        "seq"    : seq,
        "flank"  : upstream_end,
        "len"    : downstream_start - upstream_end + 1
    }

 
    with open( primer3_file, 'w+') as outfile:
        outfile.write(template.format(**details))

    outfile.close()


#
# runs primer3 and returns the output as a dict
#
# input: sequence ID, sequence, primer3_file (optional)
# output: returns a dict with the results.
#
# Kim Brugger (17 Aug 2016)
def run_primer3( seq_id, seq, primer3_file = ""):
    verbose_print( "run_primer3", 2)

    if ( primer3_file == "" ):
        primer3_file = seq_id + ".primer3"


    primer3_file = re.sub(":", "_", primer3_file)


    TMP_FILES.append( primer3_file )

    write_primer3_file(seq_id, seq, primer3_file)


    cmd = PRIMER3 + " < " + primer3_file


    verbose_print( cmd, 3)
    args = shlex.split( cmd )
    process = subprocess.Popen(cmd,
                               stdout=subprocess.PIPE,
                               shell=True)
  
    output = process.communicate()
    output_dict = dict()
  
    for line in output[0].split("\n"):

        if (line == '='):
            break

        key, value = line.split("=")

        output_dict[ key ] = value

    return output_dict



#
# Maps primers to the reference using smalt to check where they map.
# 
# input: region-ID, primer3_dict
#
# output: dict of primers and their nr of mappings to the genome
#
# Kim Brugger (17 Aug 2016)
def map_primers( region_id, primer3_dict):
    verbose_print( "check_primers", 2)

    res = dict()

    primers_file = region_id + "_p3seq.fasta"

    TMP_FILES.append( primers_file )

    with open( primers_file, 'w+') as primerfasta:

        for i in range(0, 5):

            primer_left_id = "PRIMER_LEFT_%d_SEQUENCE" % i
            primer_right_id = "PRIMER_RIGHT_%d_SEQUENCE" % i

            if ( primer_right_id not in primer3_dict or primer_left_id not in primer3_dict):
                continue

            primerfasta.write(">%s\n%s\n" % (primer_left_id, primer3_dict[ primer_left_id ]))
            primerfasta.write(">%s\n%s\n" % (primer_right_id, primer3_dict[ primer_right_id ]))

            # Store the primer seqs in the res dict as we will needthem later in the program
            left_id = "LEFT_%d" % i
            res[ left_id ] = {}
            res[ left_id ][ 'SEQ' ] = primer3_dict[ primer_left_id ]
            res[ left_id ][ 'TM' ]  = float(primer3_dict[ "PRIMER_LEFT_%d_TM" % i ])
            res[ left_id ][ 'GC_PERCENT' ]  = float(primer3_dict[ "PRIMER_LEFT_%d_GC_PERCENT" % i ])
            res[ left_id ][ 'CHR' ] = []
            res[ left_id ][ 'POS' ] = []

            right_id = "RIGHT_%d" % i
            res[ right_id ] = {}
            res[ right_id ][ 'SEQ' ] = primer3_dict[ primer_right_id ]
            res[ right_id ][ 'TM' ]  = float(primer3_dict[ "PRIMER_RIGHT_%d_TM" % i ])
            res[ right_id ][ 'GC_PERCENT' ]  = float(primer3_dict[ "PRIMER_RIGHT_%d_GC_PERCENT" % i ])
            res[ right_id ][ 'CHR' ] = []
            res[ right_id ][ 'POS' ] = []


        primerfasta.close()

    smalt_results = region_id + ".smalt"
    TMP_FILES.append( smalt_results )

    cmd = SMALT + " map  -d -1 /refs/human_1kg/human_g1k_v37 " + primers_file + "> " + smalt_results + " 2> /dev/null"

    verbose_print( cmd, 4)
    subprocess.call(cmd, shell=True)


#    pp.pprint( res )

    match = {}
    smalt_report = []
    query_region = []


    with open(smalt_results, 'rU') as smalt_output:

        for line in smalt_output:
    
            if (line.startswith("@")):
               continue
       
            line = line.rstrip("\n")
            fields = line.split("\t")

            match_id        = fields[ 0 ] 
            match_id        = re.sub(r'PRIMER_', '', match_id)
            match_id        = re.sub(r'_SEQUENCE', '', match_id)

            match_chr       = fields[ 2 ]
            match_pos       = fields[ 3 ]
            match_flag      = fields[ 1 ]
            
            match_seq       = res[ match_id ][ 'SEQ' ]
            match_length    = len( match_seq )
            match_bp_agree = int(re.sub("AS:i:", '', fields[  12 ]))

            if (match_length <= match_bp_agree + ALLOWED_MISMATCHES):
                res[ match_id ][ 'CHR' ].append( match_chr )
                res[ match_id ][ 'POS' ].append( match_pos )



    # See if the primers map uniquely or not.
    for primer  in  res.keys() :

        res[ primer ]['MAPPING_SUMMARY'] = 'unique mapping'

        nr_of_chromosomes = len(set(res[ primer ][ 'CHR' ]))
        nr_of_mappings    = len( res[ primer ][ 'POS' ])

        if (nr_of_mappings > 1 and nr_of_mappings <= MAX_MAPPINGS ):
            res[ primer ]['MAPPING_SUMMARY'] = '%d mappings' % nr_of_mappings
            res[ primer ][ 'MAPPING_SUMMARY' ] += " to chromosomes: " + ",".join(set(res[ primer ][ 'CHR' ]))


        elif (nr_of_mappings >= MAX_MAPPINGS ):
            res[ primer ][ 'MAPPING_SUMMARY' ] = '%d mappings' % nr_of_mappings
            res[ primer ][ 'MAPPING_SUMMARY' ] += " on %d chromosomes" % nr_of_chromosomes

    return res



def pick_best_primers( primer_data, chromo, start_pos, end_pos ):
    verbose_print("pick_best_primers", 2)

    # Finally we are getting to the crux of this whole ordeal: Pick the best primers.
    # It will be done using the following rules:
    # Unique MAPPING primers are best
    # primers closest to the region of interest
    # Primers generating odd products are eliminated.


    # First group the primers according to the region.

    (closest_fwd, dist_fwd) = (None, None)
    (closest_rev, dist_rev) = (None, None)

    for primer in primer_data:

        if ( primer == 'FULLSEQ'):
            continue

        if ( primer_data[ primer][ 'MAPPING_SUMMARY' ] != 'unique mapping'):
            verbose_print( "Non unique mapping ( %s )" % primer, 5)
            continue

        if ( not primer_data[ primer ][ 'CHR' ] or primer_data[ primer ][ 'CHR' ][ 0 ] != chromo ):
            verbose_print( "No mapping or Unique mapping to different chromosome (%s). Should never happen! " % primer, 5)
            continue

        if  (primer.find( 'LEFT' ) >= 0):
            
            primer_dist = start_pos - int (primer_data[ primer ][ 'POS' ][ 0 ]) + 1

            if ( primer_dist < 0 ):
                verbose_print("Primer %s downstream of region ! ( %d [%d, %d])" % (primer, primer_dist, start_pos, int (primer_data[ primer ][ 'POS' ][ 0 ])), 5)
                continue

            if ( dist_fwd is None or primer_dist < dist_fwd):
                dist_fwd    = primer_dist
                closest_fwd = primer
                continue

        
        elif( primer.find( 'RIGHT' ) >= 0):

            primer_dist =  int (primer_data[ primer ][ 'POS' ][ 0 ]) - end_pos + 1

            if ( primer_dist < 0 ):
                verbose_print( "Primer %s uptream of region ! (%d)" % (primer, primer_dist), 5)
                continue

            if ( dist_rev is None or primer_dist < dist_rev ):
                dist_rev    = primer_dist
                closest_rev = primer
                continue

    

    return closest_fwd, closest_rev



#
# reverses a DNA string. No error catching for non-valid bases or mixed bases.
#
# input: DNA string
# output: reverse DNA string
#
#
# Kim Brugger (17 Aug 2016)
def revDNA( string ):
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


#
# Find where the primers aligns in a sequence. In this case the region of interest.
# All primers are matched on both strands.
#
#
# Kim Brugger (17 Aug 2016)
def align_primers_to_seq( seq, all_primers):
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



def get_and_parse_arguments():
    verbose_print( "get_and_parse_arguments", 4)

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('-c', '--chrom')
    
    # NICK ADDED RANGE
    
    group = parser.add_mutually_exclusive_group( required = True )
    group.add_argument('-p', '--pos')
    group.add_argument('-r', '--range', nargs = 2)

    parser.add_argument('-o', '--output')
    parser.add_argument('-f', '--flank')
    parser.add_argument('-t', '--text_output')

    args = parser.parse_args()



    if ( args.range ):
        args.target_start, args.target_end = [int(x) for x in args.range]
    else:
        args.target_start = int( args.pos )
        args.target_end   = int( args.pos )

    # how many bp on either side are to be used for designing primers in
    args.target_flank = args.flank or FLANK
    args.target_flank = int( FLANK )
    
    if ( args.target_start == args.target_end):
        args.target_id = "%s:%d" % ( args.chrom, args.target_start)
        args.filename  = "%s_%d" % (args.chrom, args.target_start)
    else:
        args.target_id = "%s:%d-%d" % ( args.chrom, args.target_start, args.target_end)
        args.filename  = "%s_%d_%d" % (args.chrom, args.target_start, args.target_end)


    if ( args.output ):
        args.filename = args.output + "_" + args.filename
        
    args.filename = args.filename.rstrip(".pdf")
    args.filename = args.filename.rstrip(".txt")



    return args

#
#  Tag up a sequence with the target-region, common SNPs found in the region
#  furthermore add a flanking sequence to the target to ensure that the primers are 
#  long enough away from the target to generate sequene. Normally this is about 50bp 
#  that is lost during sanger sequencing
#
#
# Kim Brugger (17 Aug 2016)

def markup_sequence( target_sequence, target_chrom, target_start, target_end, target_flank):
    verbose_print( "tag_sequence", 2)

    sequence = list(target_sequence)
    tags = [" "] * len( target_sequence )
    # Our target base
    start = target_flank
    end = len(tags) - target_flank

    for x in range(0, len(tags)):
        if (x >= start and x <= end):
            tags[x] = '*'

        if x in range(start - TARGET_LEAD, start) or x in range(end, end + TARGET_LEAD):
            tags[x] = '-'
    #return tags
        

    verbose_print( "::::: %d - %d, %d" % (target_start, target_end, target_flank), 5)
    # Tag the target sequence, regardless of the nature of the variant only one (1) base is tagged as the target.
    sequence[ target_flank - TARGET_LEAD ] = ' [' + target_sequence [ target_flank - TARGET_LEAD ]
    sequence[(- target_flank + TARGET_LEAD) -1 ] = sequence [ (- target_flank + TARGET_LEAD) -1 ] + '] '

#    return sequence
#   exit ()

    dbSNPs = fetch_known_SNPs( '/software/dev/projects/local_dbsnp/annots-rsIDs-dbSNPv144.20150605.tab.gz', 
                               target_chrom, target_start - target_flank, target_end + target_flank )

    masked_positions = []

    for dbSNP in (dbSNPs):
        #print dbSNP
        #exit()
        if ( len( dbSNP ) < 6):
            continue
        
        #unpack dbSNP entry. 
        (snp_chrom, snp_pos, snp_id, snp_ref, snp_alt, common, vld, caf) = dbSNP

        snp_pos = int( snp_pos )
        if (snp_pos >= target_start - TARGET_LEAD and snp_pos <= target_end + TARGET_LEAD ):
            continue

        if ( common == '1'):
            #        pp.pprint( dbSNP )
            mask_pos = snp_pos - (target_start - target_flank) 

        # Sometimes the SNP we are looking at is also in dbsnp, If this is the case we dont tag it twice
            
        # In the odd case we are looking at common deletion, mask the whole region. Normally this will just be one base
            for i in range(0, len(snp_ref)):
                if (len(sequence)<=  mask_pos + i):
                    break

                # If already masked skip masking it again.
                if ( re.search('<', sequence[ mask_pos + i  ]) or 
                     re.search('\[',sequence[ mask_pos + i  ])):
                    continue

                sequence[ mask_pos + i  ] = ' <' + sequence[ mask_pos + i ] + '> '
                tags[  mask_pos + i  ] = 'X'

            else:
#        target_sequence[ snp_pos - pos + 1  ] = ' {' + target_sequence[ snp_pos - pos + 1  ] + '} '
                pass





    sequence =  "".join( sequence )
    sequence = re.sub(' ', '', sequence)


    #print tags
    #exit()
    tagged_string = "".join( tags )


    return ( sequence, tagged_string  )


#
# Check if the region where we want to put a primer is alreay 
# occupied by another primer.
#
# Input: string, start-pos, end-pos
#
# Output: Boolean. True if they do clash
#
# Kim Brugger (18 Aug 2016)
def check_if_primer_clash( mapped_list, start, end):
    
    verbose_print( "check_if_primers_clash", 4)
    for i in range (start, end):
        if (mapped_list[ i ]  != " "):
            return True 

    return False



#
# Extract the primer seqs from the primer dict
# 
# input: primer dict
#
# output: list of lists containing primer name and sequence
#
# Kim Brugger (19 Aug 2016)
def get_primer_seqs( primers ):
    primer_seqs = []
    for primer_id in sorted(primers):
        primer_seqs.append([primer_id, primers[ primer_id ][ 'SEQ']])

    return primer_seqs

#
# Aligns a set of primers to a target string. 
# Tag them up eg: >>>>> RIGHT_3 >>>>> so the string with the primer tag align to the target sequence
# Collapse primer_tag_string where possigle, ie two primers does not overlap
# And finally assign colours to each of the primer pairs.
#
# input target sequence and primers
#
# output: list of string with mapped primers, list of list with colours to apply
#
def make_mapped_primer_strings( target_sequence, primers):

    # extract the primer sequences
    primer_seqs = get_primer_seqs( primers )
    
    verbose_print( "primer_make_mapped_strings", 2)
    # Align the primers to the target sequence to see where they align
    mappings = align_primers_to_seq(target_sequence, primer_seqs )
    
    # data struct for the align primer tagging, initiate with an list the length of the region we are looking at
    mapped_strings = []
    mapped_strings.append( [" "]*len( target_sequence ) )
    # And then assign colour, by default everything is sat to -1 --> nothing in the PDF generation
    mapped_colours = []
    mapped_colours.append( [ -1 ]*len( target_sequence ) )

    #
    for mapping in mappings:
#        pp.pprint( mapping )
        ( name, primer, primer_pos, strand ) = mapping

        primer_nr = int(re.sub(r'.*_(\d)',r'\1' , name))

        mapping_index = 0

        # Calculate the number of > (or <) for the primer tag. Length of primer-seq - primer-name - 2 (spacing)
        arrows = (len(primer)-len(name)-2)/2
        arrow_type = ">"
        # Flip the type if on the nagative strand
        if ( strand == 1 ):
            primer = revDNA( primer )
            arrow_type = "<"

        # make the tag
        tag = arrow_type*arrows + " " + name + " " + arrow_type*arrows
        # If one to short ammend that. What a hack ;-)
        if ( len(tag) < len(primer)):
            tag += arrow_type
            
        primer = tag 

        mapping_index = 0
        while( 1 ):
            if (check_if_primer_clash ( mapped_strings[ mapping_index ], primer_pos, primer_pos + len( primer ))):

                mapping_index += 1
                
                # if they do clash create new lists that the mapping_index now will point to.
                if (len( mapped_strings ) >= mapping_index ):
                    mapped_strings.append([" "]*len( target_sequence ))
                    mapped_colours.append([ -1 ]*len( target_sequence ))

            else:
                break

        # Add the primer to the string mapping_index points to
        for i in range(0, len( primer)):
            mapped_strings[ mapping_index ] [ primer_pos + i ] = primer[ i ]
            mapped_colours[ mapping_index ] [ primer_pos + i ] = primer_nr


    # Make the lists into a strings
    for i in range(0, len(mapped_strings)):
        mapped_strings[ i ] = "".join( mapped_strings[ i ])

    return mapped_strings, mapped_colours


#
# makes a pretty text output of target region and primers.
#
# input: target snp/target string, mapped primers, position of first base of the reference region
#
# output: string with nicely formatted data
#
def pretty_print_mappings( target_sequence, tagged_string, primer_strings, base1):

    lines = []
    
    for i in range(0, len(target_sequence), 80):


        lines.append("%-9d  %s" % ( base1+ i, target_sequence[i: i+80]))
        lines.append("           " + tagged_string[i: i+80])

        for primer_string in ( primer_strings ):

            line =  "           " + primer_string[i: i+80]
            if ( re.match(r'^ *$', line)):
                continue
            lines.append( line )

        lines.append( "" )

    return lines


#
# makes a pretty text output of target region and primers in a pdf canvas
#
# input: canvas, 
#        offset (where to start writing) 
#        target sequence
#        target snp/target string, 
#        mapped primers,
#        colours to apply to the primer arrows
#        position of first base of the reference region
#
# output: string with nicely formatted data
#
def pretty_pdf_mappings(c, top_offset,  target_sequence, tagged_string, primer_strings, primer_colours, base1):
    verbose_print("pretty_pdf_mappings", 2)

    lines = []
    
    for i in range(0, len(target_sequence), 80):


#        c.drawString(40, top_offset, "%-9d  %s" % ( base1+ i, target_sequence[i: i+80]))
#        top_offset -= 8
#        c.drawString(40, top_offset, "           " + tagged_string[i: i+80])
#        top_offset -= 8

        p_line = "%-9d  %s" % ( base1+ i, target_sequence[i: i+80])
        m_line = "           " + tagged_string[i: i+80]


        x_offset = 40
        for k in range(0, len(p_line)):

            if (m_line[ k ] == "X"):
                c.setFillColorRGB(255,0,0)
            elif (m_line[ k ] == "*"):
                c.setFillColorRGB(0,190,0)

            c.drawString(x_offset , top_offset, p_line[k])
            x_offset += stringWidth(" ", 'mono', 8)
            c.setFillColorRGB(0,0,0)


        top_offset -= 8

        m_line = re.sub(r'X', r' ', m_line)

        if (re.search(r'\*', m_line) or re.search(r'\-', m_line)):
            c.setFillColorRGB(0,190,0)
            c.drawString(40 , top_offset, m_line)
            c.setFillColorRGB(0,0,0)

            top_offset -= 8


#        if (m_line.search(r'\*', name) :


        for j in range(0, len(primer_strings)):
            primer_string = primer_strings[ j ]
            primer_colour = primer_colours[ j ]

            line = primer_string[i: i+80]
            if ( re.match(r'^ *$', line)):
                continue

            x_offset = 40 + stringWidth(" ", 'mono', 8)*11

            for k in range(i, i+80):
                if ( k > len(target_sequence) - 1):
                    break

                if ( primer_colour[k] >= 0 ):

                    c.setFillColorRGB(colours[ primer_colour[k] ][0], 
                                      colours[ primer_colour[k] ][1],
                                      colours[ primer_colour[k] ][2])

                c.drawString(x_offset , top_offset, primer_string[k])
                x_offset += stringWidth(" ", 'mono', 8)
                c.setFillColorRGB(0,0,0)


            top_offset -= 8

        top_offset -= 8

    return top_offset


#
# makes a pretty text output of primer informatino
#
# input: 
#        offset (where to start writing) 
#        target sequence
#        target snp/target string, 
#        mapped primers,
#        colours to apply to the primer arrows
#        position of first base of the reference region
#
# output: string with nicely formatted data
#

#
# makes a pretty text output of primer informatino
#
# input: 
#        pdf canvas
#        offset (where to start writing) 
#        target start
#        target end
#        primers
#        best fwd primer
#        best rev primer
#
# output: string with nicely formatted data
#
def pretty_print_primer_data(primer_dict, target_chrom, target_start, target_end ):

    lines = []


    verbose_print( "extract_passed_primer_seqs", 3)

    lines.append( "\n" )
    lines.append( "\n" )
    lines.append( "\n" )


    lines.append( "_-=-"*15 +"_" )


    if ( target_start == target_end ):
        lines.append( " Primer design report for chr: %s position: %d" % (target_chrom, target_start))
    else:
        lines.append( " Primer design report for chr: %s range: %d-%d" % (target_chrom, target_start, target_end))

    lines.append( "_-=-"*15 +"_")
    lines.append( "\n")

#    lines.append( "\t".join(['ID', '%GC', 'TM', 'Sequence']))
    lines.append( "ID         %GC    TM     Primer sequence           Mapping(s)    ")
    lines.append( "_-=-"*15 +"_")

    primer_seqs = []
    for primer in sorted(primer_dict):
        if ( primer == 'FULLSEQ'):
            continue

        name = primer
        name = re.sub(r'PRIMER_', '', name)
        name = re.sub(r'_SEQUENCE', '', name)


        lines.append( "%-10s %.2f  %.2f  %-25s %s" % (name, 
                                                     primer_dict[ name ][ "GC_PERCENT"], 
                                                     primer_dict[ name ][ "TM"],
                                                     primer_dict[ name ][  "SEQ"], 
                                                     primer_dict[ name ][ 'MAPPING_SUMMARY' ]))


    lines.append( "" )
    lines.append( "-="*46 )
    lines.append( "" )

    return lines



#
# makes a pretty text output of primer informatino
#
# input: 
#        pdf canvas
#        offset (where to start writing) 
#        target start
#        target end
#        primers
#        best fwd primer
#        best rev primer
#
# output: string with nicely formatted data
#
def pretty_pdf_primer_data(c, y_offset, target_chrom, target_start, target_end, passed_primers, fwd_primer=None, rev_primer=None ):
    verbose_print("pretty_pdf_primer_data", 2)

#    c.drawString(40 , y_offset, "_-=-"*15 +"_" )
    c.line(40,y_offset,width - 40  ,y_offset+2)
    y_offset -= 8


    if ( target_start == target_end ):
        c.drawString(40 , y_offset, "Primer design report for chr: %s position: %d" % (target_chrom, target_start))
    else:
        c.drawString(40 , y_offset, "Primer design report for chr: %s range: %d-%d" % (target_chrom, 
                                                                                       target_start, target_end))
            
    y_offset -= 8
    c.line(40,y_offset,width - 40 ,y_offset+2)
    y_offset -= 16

    c.drawString(40 , y_offset, "ID         %GC    TM     Primer sequence           Best primer  Mapping(s)    ")
    y_offset -= 8
    c.line(40,y_offset,width - 40 ,y_offset+2)
    y_offset -= 8




    primer_seqs = []
    for primer in sorted(passed_primers):
        if ( primer == 'FULLSEQ'):
            continue


        name = primer
        name = re.sub(r'PRIMER_', '', name)
        name = re.sub(r'_SEQUENCE', '', name)

        
        if (name == "RIGHT_0"):
            y_offset -= 8


#        lines.append( "\t".join([name, 
#                         primer3_results[ "PRIMER_" + name + "_GC_PERCENT"], 
#                         primer3_results[ "PRIMER_" + name + "_TM"],
#                         primer3_results[ "PRIMER_" + name + "_SEQUENCE"], 
#                         passed_primers[ "PRIMER_" + name + "_SEQUENCE" ][ 'MAPPING_SUMMARY' ]]))

        picked_primer = ' '

        if ( fwd_primer == primer or rev_primer == primer):
            picked_primer = 'Y'


#        c.drawString(40 , y_offset, "%-10s %.2f  %.2f  %-25s %s            %s" % ("", 
#                                                     float(primer3_results[ "PRIMER_" + name + "_GC_PERCENT"]), 
#                                                     float(primer3_results[ "PRIMER_" + name + "_TM"]),
#                                                     primer3_results[ "PRIMER_" + name + "_SEQUENCE"], picked_primer,
#                                                     passed_primers[ "PRIMER_" + name + "_SEQUENCE" ][ 'MAPPING_SUMMARY' ]))



        c.drawString(40 , y_offset, "%-10s %.2f  %.2f  %-25s %s            %s" % ("", 
                                    passed_primers[ name ][ "GC_PERCENT"], 
                                    passed_primers[ name ][ "TM"],
                                    passed_primers[ name ][ "SEQ"], picked_primer,
                                    passed_primers[ name ][ 'MAPPING_SUMMARY' ]))


        primer_nr = re.sub(r'.*_(\d)',r'\1' , name)

#        pp.pprint( primer_nr )
#        print pp.pprint(colours[ int( primer_nr ) ])
        

        c.setFillColorRGB(colours[ int( primer_nr ) ][0], 
                          colours[ int( primer_nr ) ][1],
                          colours[ int( primer_nr ) ][2])

        c.drawString(40 , y_offset, name)
        c.setFillColorRGB(0,0,0)
        y_offset -= 8


            



    y_offset -= 8
    c.line(40,y_offset,width - 40 ,y_offset+2)
    y_offset -= 8
    y_offset -= 8


    return y_offset

def pretty_primer_data(outfile, target_chrom, target_start, target_end,  passed_primers, fwd_primer=None, rev_primer=None  ):
    verbose_print("pretty_primer_data", 2)


    fh = open( outfile, 'w')


    lines = []

    if ( target_start == target_end ):
        lines.append("Primer design report for chr: %s position: %d" % (target_chrom, target_start))
    else:
        lines.append("Primer design report for chr: %s range: %d-%d" % (target_chrom, target_start, target_end ))

    lines.append("ID\t%GC\tTM\tPrimer sequence\tBest primer\tMapping(s)")

#    pp.pprint( passed_primers )

    for primer in sorted(passed_primers):
        if ( primer == 'FULLSEQ'):
            continue

        name = primer
        name = re.sub(r'PRIMER_', '', name)
        name = re.sub(r'_SEQUENCE', '', name)

        
        picked_primer = ''

        if ( fwd_primer == primer or rev_primer == primer):
            picked_primer = 'Y'

        lines.append("\t".join([name, 
                              "%.2f" % passed_primers[ name ][ "GC_PERCENT"], 
                              "%.2f" % passed_primers[ name ][ "TM"],
                              passed_primers[ name ][ "SEQ"], picked_primer,
                              passed_primers[ name ][ 'MAPPING_SUMMARY' ]]))




        # lines.append("\t".join([name, 
        #                       primer3_results[ "PRIMER_" + name + "_GC_PERCENT"], 
        #                       primer3_results[ "PRIMER_" + name + "_TM"],
        #                       primer3_results[ "PRIMER_" + name + "_SEQUENCE"], picked_primer,
        #                       passed_primers[ "PRIMER_" + name + "_SEQUENCE" ][ 'MAPPING_SUMMARY' ]]))




    lines.append("\n")
    fh.write("\n".join( lines ))
    fh.close()



def method_blurb():

    lines = []

    lines.append('primer-designer version: ' + VERSION + ' using dbSNP 144 for SNP checking, and human reference GRCh37.')

    lines.append('Common SNP annotation: A common SNP is one that has at least one 1000Genomes population with a minor ')
    lines.append('allele of frequency >= 1% and for which 2 or more founders contribute to that minor allele frequency.')

    return lines


def pretty_pdf_method(c, top_offset):
    verbose_print("pretty_pdf_method", 3)

    lines = method_blurb()

    top_offset = 80

    for line in lines:
        c.drawString(40, top_offset, line)
        top_offset -= 8



def map_target_region_to_exon(chrom, target_start, target_end):

    gene = {}
    (coding_start, coding_end) = (target_start, target_end)
    (coding_start, coding_end) = (None, None)

    tbx = pysam.Tabixfile("gencode.v19.bed.gz")
    for row in tbx.fetch(chrom, target_start - 10, target_end+10):
        (exon_chrom, exon_start, exon_end, exon_gene, exon_nr, exon_transcript) = row.split("\t")
        
        exon_start = int( exon_start )
        exon_end   = int( exon_end   )

        if ( coding_start is  None or  coding_start > exon_start ):
            coding_start = exon_start

        if ( coding_end is None or coding_end < exon_end ):
            coding_end = exon_end
        
#            print (str(row))

#    print "New region is: %s:%d-%d" % ( chrom, coding_start, coding_end)
    return ( chrom, coding_start, coding_end )






#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
#    Main loop !
#
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def main():

    args = get_and_parse_arguments()

    target_chrom    = args.chrom
    target_start    = args.target_start
    target_end      = args.target_end
    target_flank    = args.target_flank
    target_id       = args.target_id
    target_filename = args.filename

    ( exon_chrom, exon_start, exon_end) = map_target_region_to_exon( target_chrom, target_start, target_end)

    if ( exon_start is not None and  exon_end is not None  and exon_end - exon_start + 1 < MAX_PRODUCT):
        verbose_print("Changing region to an exonic region (%s:%d-%d)" % (exon_chrom, exon_start, exon_end), 1)
        target_start = exon_start
        target_end   = exon_end

#    exit()


    target_sequence                   =   fetch_region( target_chrom, target_start - target_flank, target_end + target_flank     )
    (tagged_sequence, tagged_region)  =   markup_sequence(target_sequence, target_chrom, target_start, target_end, target_flank )

    primer3_results  = run_primer3( target_id , tagged_sequence, "%s_%d.primer3" % ( target_chrom, target_start))
    primer_dict   = map_primers( target_id, primer3_results)
#    pp.pprint( passed_primers )

    best_fwd_primer, best_rev_primer = pick_best_primers(primer_dict, target_chrom, target_start, target_end)
    (mapped_primer_strings, mapped_primer_colours) = make_mapped_primer_strings( target_sequence, primer_dict)

#    pp.pprint( mapped_primer_strings )
#    pp.pprint( mapped_primer_colours )


    pretty_primer_data("%s.txt" % target_filename, target_chrom, target_start, target_end, primer_dict, best_fwd_primer, best_rev_primer )

    if (  args.text_output ):
        lines  = pretty_print_primer_data( primer_dict, target_chrom, target_start, target_end )
        lines += pretty_print_mappings( target_sequence, tagged_region, mapped_primer_strings, target_start - FLANK)
        print "\n".join(lines)
        exit()

    else:
        target_filename = "%s.pdf" % target_filename
        c = canvas.Canvas( target_filename , pagesize=A4)
    
        width, height = A4 #keep for later
        
        font = TTFont('mono', '/usr/share/fonts/truetype/ttf-liberation/LiberationMono-Regular.ttf')
        pdfmetrics.registerFont( font )
        c.setFont('mono', 7)
        top_offset = pretty_pdf_primer_data(c, height - 30, target_chrom, target_start, target_end, primer_dict, best_fwd_primer, best_rev_primer )
            
        c.setFont('mono', 8)
        pretty_pdf_mappings(c, top_offset,  target_sequence, tagged_region, mapped_primer_strings, mapped_primer_colours, target_start - FLANK)

        pretty_pdf_method(c, top_offset)

        c.showPage()
        c.save()

        for tmp_filename in     TMP_FILES:
            print "deleting tmp file: %s " % tmp_filename
            os.remove( tmp_filename )



if __name__ == '__main__':
    main()


