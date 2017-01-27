# 
# Module for designing primers. This is a split out of the initial
# program as some of the functinality will be used in the django part
# of the project
# 
# 
# Kim Brugger (25 Jan 2017), contact: kim@brugger.dk

import sys
import re
import os
import pprint
pp = pprint.PrettyPrinter(indent=4)

import shlex
import subprocess


import primer_design.core as core

sys.path.append("/software/lib/python2.7/site-packages/pysam-0.7.5-py2.7-linux-x86_64.egg")
import pysam


colours = [[255,   0,   0], # red
           [  0, 255,   0], # green
           [  0,   0, 255], # blue
           [255,   0, 255], # Pink
           [0,   255, 255],
           [255, 255,   0], #Dark Green
           [100, 255, 100]] # Yellow, crap!

# External tools we use 
SMALT      = '/software/bin/smalt-0.7.6 '
PRIMER3    = '/software/bin/primer3_core '



# default parameters, can be overridden by the user (perhaps)
TARGET_LEAD          = 50
FLANK                = 500 # amount of bp on either side of the target to design primers in
NR_PRIMERS           = 4   # Nr of primers to report back
ALLOWED_MISMATCHES   = 4   # Maximum nr of errors when mapping a primer back to the reference
MAX_MAPPINGS         = 5   # Less or equal number of mappings to the reference what chromosomes mapped to are named 
MAX_PRODUCT          = 1200
MAX_3_PRIME_MISMATCH = 2

VERSION    =  '2.0-beta1'


import primer_design.core






def markup_sequence( target_sequence, target_chrom, target_start, target_end, target_flank):
    """

    Tag up a sequence with the target-region, common SNPs found in the region
    furthermore add a flanking sequence to the target to ensure that the primers are 
    long enough away from the target to generate sequene. Normally this is about 50bp 
    that is lost during sanger sequencing


    Kim Brugger (17 Aug 2016)

    """


    core.verbose_print( "tag_sequence", 2)

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
        

    core.verbose_print( "::::: %d - %d, %d" % (target_start, target_end, target_flank), 5)
    # Tag the target sequence, regardless of the nature of the variant only one (1) base is tagged as the target.
    sequence[ target_flank - TARGET_LEAD ] = ' [' + target_sequence [ target_flank - TARGET_LEAD ]
    sequence[(- target_flank + TARGET_LEAD) -1 ] = sequence [ (- target_flank + TARGET_LEAD) -1 ] + '] '

#    return sequence
#   exit ()

#    dbSNPs = fetch_known_SNPs( '/software/dev/projects/local_dbsnp/annots-rsIDs-dbSNPv144.20150605.tab.gz', 
#                               target_chrom, target_start - target_flank, target_end + target_flank )


    dbSNPs = core.fetch_known_SNPs( '/refs/human_1kg/annots-rsIDs-dbSNPv144.20150605.tab.gz', target_chrom, 
                                                  target_start - target_flank, 
                                                  target_end + target_flank )

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
#                sequence[ mask_pos + i  ] = 
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


def map_target_region_to_exon(chrom, target_start, target_end):
    """
    function to see if the target is in an exonic region, if it is
    expand the region so it includes the whole region

    input: Chrom, target_start, target_end

    return: adjusted chrom, start, end

    """

    gene = {}
    (coding_start, coding_end) = (target_start, target_end)
    (coding_start, coding_end) = (None, None)
    path = os.path.dirname(os.path.realpath(__file__))
    tbx = pysam.Tabixfile(path + "/gencode.v19.bed.gz")
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



def map_primers_smalt( region_id, primer3_dict, target_chrom, target_start, target_end):
    """

    Maps primers to the reference using smalt to check where they map.
 
    input: region-ID, primer3_dict
    
    output: dict of primers and their nr of mappings to the genome
    
    Kim Brugger (17 Aug 2016)
    """    

    core.verbose_print( "check_primers", 2)

    res = dict()

    primers_file = region_id + "_p3seq.fasta"

    core.add_new_tmp_file( primers_file )

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
    core.add_new_tmp_file( smalt_results )

    cmd = SMALT + " map  -d -1 /refs/human_1kg/human_g1k_v37 " + primers_file + "> " + smalt_results + " 2> /dev/null"

    core.verbose_print( cmd, 4)
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

            match_chrom     = fields[ 2 ]
            match_pos       = int(fields[ 3 ])
            match_flag      = fields[ 1 ]
            
            match_seq       = res[ match_id ][ 'SEQ' ]
            match_length    = len( match_seq )
            match_bp_agree = int(re.sub("AS:i:", '', fields[  12 ]))

            if (match_length <= match_bp_agree + ALLOWED_MISMATCHES):
                res[ match_id ][ 'CHR' ].append( match_chrom )
                res[ match_id ][ 'POS' ].append( match_pos )


            if ( target_chrom == match_chrom and 
                 target_start < match_pos and
                 target_end  > match_pos ):
                res[ match_id ][ 'TARGET_POS'   ] = match_pos
                res[ match_id ][ 'TARGET_CHROM' ] = match_chrom



    # See if the primers map uniquely or not.
    # If a primer does not map to the region of interest (low complexity etc) remove it.                
    for primer  in  res.keys() :

        if 'TARGET_CHROM' not in res[ primer ]:
            del res[ primer ]
            continue

        res[ primer ]['MAPPING_SUMMARY'] = 'unique mapping'

        nr_of_chromosomes = len(set(res[ primer ][ 'CHR' ]))
        nr_of_mappings    = len( res[ primer ][ 'POS' ])

        if (nr_of_mappings > 1 and nr_of_mappings <= MAX_MAPPINGS ):
            res[ primer ]['MAPPING_SUMMARY'] = '%d mappings' % nr_of_mappings
            res[ primer ][ 'MAPPING_SUMMARY' ] += " to chromosome: " + ",".join(set(res[ primer ][ 'CHR' ]))


        elif (nr_of_mappings >= MAX_MAPPINGS ):
            res[ primer ][ 'MAPPING_SUMMARY' ] = '%d mappings' % nr_of_mappings
            res[ primer ][ 'MAPPING_SUMMARY' ] += " on %d chromosomes" % nr_of_chromosomes

    return res





def check_if_primer_clash( mapped_list, start, end):
    """
    Check if the region where we want to put a primer is alreay 
    occupied by another primer.

    Input: string, start-pos, end-pos
    
    Output: Boolean. True if they do clash
    
    Kim Brugger (18 Aug 2016)
    """
    
    core.verbose_print( "check_if_primers_clash", 4)
    for i in range (start, end):
        if (mapped_list[ i ]  != " "):
            return True 

    return False



def get_primer_seqs( primers ):
    """

    Extract the primer seqs from the primer dict
    
    input: primer dict
    
    output: list of lists containing primer name and sequence
    
    Kim Brugger (19 Aug 2016)
    """

    primer_seqs = []
    for primer_id in sorted(primers):
        primer_seqs.append([primer_id, primers[ primer_id ][ 'SEQ']])

    return primer_seqs



def make_mapped_primer_strings( target_sequence, primer_dict):
    """

    Aligns a set of primers to a target string. 
    Tag them up eg: >>>>> RIGHT_3 >>>>> so the string with the primer tag align to the target sequence
    Collapse primer_tag_string where possigle, ie two primers does not overlap
    And finally assign colours to each of the primer pairs.
    
    input target sequence and primer dict
    
    output: list of string with mapped primers, list of list with colours to apply
    
    """

    # extract the primer sequences
    primer_seqs = get_primer_seqs( primer_dict )
    
    core.verbose_print( "primer_make_mapped_strings", 2)
    # Align the primers to the target sequence to see where they align
    mappings = core.align_primers_to_seq(target_sequence, primer_seqs )
    
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
            primer = core.revDNA( primer )
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


def pretty_print_mappings( target_sequence, tagged_string, primer_strings, base1):
    """
    makes a pretty text output of target region and primers.
    
    input: target snp/target string, mapped primers, position of first base of the reference region
    
    output: string with nicely formatted data
    
    """

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


def pick_best_primers( primer_data, chromo, start_pos, end_pos ):
    """
     Finally we are getting to the crux of this whole ordeal: Pick the best primers.
     It will be done using the following rules:
     Unique MAPPING primers are best
     primers closest to the region of interest
     Primers generating odd products are eliminated.

    """

    core.verbose_print("pick_best_primers", 2)

    pp.pprint( primer_data )



    # First group the primers according to the region.

    (closest_fwd, dist_fwd) = (None, None)
    (closest_rev, dist_rev) = (None, None)

    for primer in primer_data:

        if ( primer == 'FULLSEQ'):
            continue

        if ( primer_data[ primer][ 'MAPPING_SUMMARY' ] != 'unique mapping'):
            core.verbose_print( "Non unique mapping ( %s )" % primer, 5)
            continue

        if ( not primer_data[ primer ][ 'CHR' ] or primer_data[ primer ][ 'CHR' ][ 0 ] != chromo ):
            core.verbose_print( "No mapping or Unique mapping to different chromosome (%s). Should never happen! " % primer, 5)
            pp.pprint( primer_data[ primer ] )
            continue

        if  (primer.find( 'LEFT' ) >= 0):
            
            primer_dist = start_pos - int (primer_data[ primer ][ 'POS' ][ 0 ]) + 1

            if ( primer_dist < 0 ):
                core.verbose_print("Primer %s downstream of region ! ( %d [%d, %d])" % (primer, primer_dist, start_pos, int (primer_data[ primer ][ 'POS' ][ 0 ])), 5)
                continue

            if ( dist_fwd is None or primer_dist < dist_fwd):
                dist_fwd    = primer_dist
                closest_fwd = primer
                continue

        
        elif( primer.find( 'RIGHT' ) >= 0):

            primer_dist =  int (primer_data[ primer ][ 'POS' ][ 0 ]) - end_pos + 1

            if ( primer_dist < 0 ):
                core.verbose_print( "Primer %s uptream of region ! (%d)" % (primer, primer_dist), 5)
                continue

            if ( dist_rev is None or primer_dist < dist_rev ):
                dist_rev    = primer_dist
                closest_rev = primer
                continue

    

    return closest_fwd, closest_rev




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


    core.verbose_print( "extract_passed_primer_seqs", 3)

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
    lines.append( "ID         POS             %GC    TM     Primer sequence      Best primer     Mapping(s)    ")
    lines.append( "_-=-"*15 +"_")

    primer_seqs = []
    for primer in sorted(primer_dict):
        if ( primer == 'FULLSEQ'):
            continue

        name = primer
        name = re.sub(r'PRIMER_', '', name)
        name = re.sub(r'_SEQUENCE', '', name)

        lines.append( "%-10s %-14s  %.2f  %.2f  %-25s %s" % (name,
                                                            "%s:%d" % (primer_dict[ name ][ "TARGET_CHROM"], primer_dict[ name ][ "TARGET_POS"]),
                                                            primer_dict[ name ][ "GC_PERCENT"], 
                                                            primer_dict[ name ][ "TM"],
                                                            primer_dict[ name ][  "SEQ"], 
                                                            primer_dict[ name ][ 'MAPPING_SUMMARY' ]))


    lines.append( "" )
    lines.append( "-="*46 )
    lines.append( "" )

    return lines



def pretty_primer_data(outfile, target_chrom, target_start, target_end,  primer_dict, fwd_primer=None, rev_primer=None  ):
    core.verbose_print("pretty_primer_data", 2)

    fh = None

    if outfile:
        fh = open( outfile, 'w')


    lines = []

    if ( target_start == target_end ):
        lines.append("Primer design report for chr: %s position: %d" % (target_chrom, target_start))
    else:
        lines.append("Primer design report for chr: %s range: %d-%d" % (target_chrom, target_start, target_end ))

    lines.append("ID\tPosition\t%GC\tTM\tPrimer sequence\tBest primer\tMapping(s)")

#    pp.pprint( passed_primers )

    for primer in sorted( primer_dict ):
        if ( primer == 'FULLSEQ'):
            continue

        name = primer
        name = re.sub(r'PRIMER_', '', name)
        name = re.sub(r'_SEQUENCE', '', name)

        
        picked_primer = ''

        if ( fwd_primer == primer or rev_primer == primer):
            picked_primer = 'Y'

        lines.append("\t".join([name, 
                              "%s:%d" % (primer_dict[ name ][ "TARGET_CHROM"], primer_dict[ name ][ "TARGET_POS"]),
                              "%.2f" % primer_dict[ name ][ "GC_PERCENT"], 
                              "%.2f" % primer_dict[ name ][ "TM"],
                              primer_dict[ name ][ "SEQ"], picked_primer,
                              primer_dict[ name ][ 'MAPPING_SUMMARY' ]]))




        # lines.append("\t".join([name, 
        #                       primer3_results[ "PRIMER_" + name + "_GC_PERCENT"], 
        #                       primer3_results[ "PRIMER_" + name + "_TM"],
        #                       primer3_results[ "PRIMER_" + name + "_SEQUENCE"], picked_primer,
        #                       passed_primers[ "PRIMER_" + name + "_SEQUENCE" ][ 'MAPPING_SUMMARY' ]]))




    lines.append("\n")
    if ( fh ):
        fh.write("\n".join( lines ))
        fh.close()
    else:
        print "\n".join( lines )



def method_blurb():

    lines = []

    lines.append('primer-designer version: ' + VERSION + ' using dbSNP 144 for SNP checking, and human reference GRCh37.')

    lines.append('Common SNP annotation: A common SNP is one that has at least one 1000Genomes population with a minor ')
    lines.append('allele of frequency >= 1% and for which 2 or more founders contribute to that minor allele frequency.')

    return lines













