
import re

import config
import core


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
def pretty_primer_data(primer_dict, target_chrom, target_start, target_end ):

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

    lines = []

    lines.append('primer-designer version: ' + config.VERSION + ' using dbSNP 144 for SNP checking, and human reference GRCh37.')

    lines.append('Common SNP annotation: A common SNP is one that has at least one 1000Genomes population with a minor ')
    lines.append('allele of frequency >= 1% and for which 2 or more founders contribute to that minor allele frequency.')

    return lines
