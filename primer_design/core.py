"""
  General low level core functions that should be living in a sep module

"""

import sys
import shlex
import subprocess
import re
import os


sys.path.append("/software/lib/python2.7/site-packages/pysam-0.7.5-py2.7-linux-x86_64.egg")
import pysam


verbose_level    =  20

tmp_files = []

import config


if ( __name__ == '__main__'):
    sys.stderr.write( "This is a module and not meant to be run as a stand alone program\n" )
    exit( 1 )






def print_to_file_or_stdout( lines, outfile=None):

    
    if outfile:
        fh = open( outfile, 'w')
        fh.write("\n".join( lines ))
        fh.close()

    else:
        print "\n".join( lines )




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



def markup_sequence( target_sequence, target_chrom, target_start, target_end, target_flank):
    """

    Tag up a sequence with the target-region, common SNPs found in the region
    furthermore add a flanking sequence to the target to ensure that the primers are 
    long enough away from the target to generate sequene. Normally this is about 50bp 
    that is lost during sanger sequencing


    Kim Brugger (17 Aug 2016)

    """


    verbose_print( "tag_sequence", 2)

    sequence = list(target_sequence)
    tags = [" "] * len( target_sequence )
    # Our target base
    start = target_flank
    end = len(tags) - target_flank

    for x in range(0, len(tags)):
        if (x >= start and x <= end):
            tags[x] = '*'

        if x in range(start - config.TARGET_LEAD, start) or x in range(end, end + config.TARGET_LEAD):
            tags[x] = '-'
    #return tags
        

    verbose_print( "::::: %d - %d, %d" % (target_start, target_end, target_flank), 5)
    # Tag the target sequence, regardless of the nature of the variant only one (1) base is tagged as the target.
    sequence[ target_flank - config.TARGET_LEAD ] = ' [' + target_sequence [ target_flank - config.TARGET_LEAD ]
    sequence[(- target_flank + config.TARGET_LEAD) -1 ] = sequence [ (- target_flank + config.TARGET_LEAD) -1 ] + '] '

#    return sequence
#   exit ()

#    dbSNPs = fetch_known_SNPs( '/software/dev/projects/local_dbsnp/annots-rsIDs-dbSNPv144.20150605.tab.gz', 
#                               target_chrom, target_start - target_flank, target_end + target_flank )


    dbSNPs = fetch_known_SNPs( '/refs/human_1kg/annots-rsIDs-dbSNPv144.20150605.tab.gz', target_chrom, 
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
        if (snp_pos >= target_start - config.TARGET_LEAD and snp_pos <= target_end + config.TARGET_LEAD ):
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
    tbx = pysam.Tabixfile(path + "/../gencode.v19.bed.gz")
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






def check_if_primer_clash( mapped_list, start, end):
    """
    Check if the region where we want to put a primer is alreay 
    occupied by another primer.

    Input: string, start-pos, end-pos
    
    Output: Boolean. True if they do clash
    
    Kim Brugger (18 Aug 2016)
    """
    
    verbose_print( "check_if_primers_clash", 4)
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

