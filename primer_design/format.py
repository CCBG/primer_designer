

import re

import config
import core



"""
  General low level core functions that should be living in a sep module

"""


if ( __name__ == '__main__'):
    sys.stderr.write( "This is a module and not meant to be run as a stand alone program\n" )
    exit( 1 )


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
    primer_seqs = core.get_primer_seqs( primer_dict )
    
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




def pretty_mappings( target_sequence, tagged_string, primer_strings, base1):
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



def pretty_primer_data(primer_dict, target_chrom, target_start, target_end ):

    """
    makes a pretty text output of primer informatino

    input: 
          pdf canvas
          offset (where to start writing) 
          target start
          target end
          primers
          best fwd primer
          best rev primer

    output: string with nicely formatted data
    """

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



def tabbed_primer_data(target_chrom, target_start, target_end,  primer_dict, fwd_primer=None, rev_primer=None  ):
    """
    tabbed primer output, mainly to be used for the website of this project

    """

    core.verbose_print("pretty_primer_data", 2)


    lines = []

    if ( target_start == target_end ):
        lines.append("Primer design report for chr: %s position: %d" % (target_chrom, target_start))
    else:
        lines.append("Primer design report for chr: {chrom} range: {start}-{end}".format(chrom=target_chrom, start=target_start, end=target_end ))

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

    lines.append("\n")

    return lines







def method_blurb():
    """
    The medthod section for the reports.
    """

    lines = []

    lines.append('primer-designer version: ' + config.VERSION + ' using dbSNP 144 for SNP checking, and human reference GRCh37.')

    lines.append('Common SNP annotation: A common SNP is one that has at least one 1000Genomes population with a minor ')
    lines.append('allele of frequency >= 1% and for which 2 or more founders contribute to that minor allele frequency.')

    return lines
