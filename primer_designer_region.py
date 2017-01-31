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

from reportlab.lib.pagesizes import A4

import primer_design.core
import primer_design.config
import primer_design.primer3
import primer_design.blast
import primer_design.printing
import primer_design.format
import primer_design.analysis
import primer_design.pdf


def get_and_parse_arguments():
    primer_design.core.verbose_print( "get_and_parse_arguments", 4)

    import argparse

    parser = argparse.ArgumentParser()

    group = parser.add_mutually_exclusive_group( required = True )
    group.add_argument('-r', '--range' )

    parser.add_argument('-o', '--output')
    parser.add_argument('-t', '--text_output', action="store_true")
    parser.add_argument('-w', '--website_output', action="store_true")

    args = parser.parse_args()

    pp.pprint(args.range)
    assert re.match(r'(\d+|[XY]):\d+$', args.range ) or re.match(r'(\d+|[XY]):\d+-\d+$', args.range), "Range needs to be in the following format: Chr:Pos or Chr:Pos-Pos"
    (args.chrom, region) = args.range.split(":")
    
    if ( re.match(r'\d+-\d+', region)):
        args.target_start, args.target_end = map( int, (region.split("-")))
    else:
        print region
        args.target_start, args.target_end = map( int, [region, region])
        
    # how many bp on either side are to be used for designing primers in
    
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







#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
#    Main loop !
#
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def main():

    args = get_and_parse_arguments()
#    pp.pprint( args )
#    exit()

    target_chrom    = args.chrom
    target_start    = args.target_start
    target_end      = args.target_end
    target_flank    = primer_design.config.FLANK
    target_id       = args.target_id
    target_filename = args.filename

    ( exon_chrom, exon_start, exon_end) = primer_design.core.map_target_region_to_exon( target_chrom, target_start, target_end)

    if ( exon_start is not None and  exon_end is not None  and exon_end - exon_start + 1 < primer_design.config.MAX_PRODUCT):
        primer_design.core.verbose_print("Changing region to an exonic region (%s:%d-%d)" % (exon_chrom, exon_start, exon_end), 1)
        target_start = exon_start
        target_end   = exon_end

#    exit()


    target_sequence                   =   primer_design.core.fetch_region( target_chrom, target_start - target_flank, target_end + target_flank     )
    (tagged_sequence, tagged_region)  =   primer_design.primer3.markup_sequence(target_sequence, target_chrom, target_start, target_end, target_flank )

    primer3_results  = primer_design.primer3.run( target_id , tagged_sequence, "%s_%d.primer3" % ( target_chrom, target_start))
#    primer_dict   = map_primers_smalt( target_id, primer3_results, target_chrom, target_start - target_flank, target_end+ target_flank)
    primer_dict   = primer_design.blast.map_primers( target_id, primer3_results, target_chrom, target_start - target_flank, target_end+ target_flank)
#    pp.pprint( primer_dict )
    


#    best_fwd_primer, best_rev_primer = primer_design.analysis.pick_best_primers(primer_dict, target_chrom, target_start, target_end)
    (mapped_primer_strings, mapped_primer_colours) = primer_design.format.make_mapped_primer_strings( target_sequence, primer_dict)

    ePCR_products = primer_design.analysis.digital_PCR( primer_dict )
    (best_fwd_primer, best_rev_primer, size) = primer_design.analysis.best_product( ePCR_products )
    multiple_mappings = primer_design.analysis.multiple_products( ePCR_products )

    for primer in multiple_mappings:
        primer_dict[ primer][ 'MULTIPLE_MAPPINGS' ] = ", ".join( multiple_mappings[ primer ])


#    pp.pprint( mapped_primer_strings )
#    pp.pprint( mapped_primer_colours )



    if (  args.website_output ):

        lines = primer_design.format.tabbed_primer_data(target_chrom, target_start, target_end, primer_dict, best_fwd_primer, best_rev_primer )

        primer_design.core.print_to_file_or_stdout( lines )

    elif (  args.text_output ):
        lines  = primer_design.format.pretty_primer_data( primer_dict, target_chrom, target_start, target_end )
        lines += primer_design.format.pretty_mappings( target_sequence, tagged_region, mapped_primer_strings, target_start - primer_design.config.FLANK)
        primer_design.format.pretty_primer_data("%s.txt" % target_filename, target_chrom, target_start, target_end, primer_dict, best_fwd_primer, best_rev_primer )
        print "\n".join(lines)
        exit()

    else:
        target_filename = "%s.pdf" % target_filename

        primer_design.pdf.create_canvas( target_filename )
    
        width, height = A4 #keep for later

        print height
        primer_design.pdf.set_top_offset( height - 30 )        
        primer_design.pdf.primer_data(target_chrom, target_start, target_end, primer_dict, best_fwd_primer, best_rev_primer )
            
        primer_design.pdf.mappings(target_sequence, tagged_region, mapped_primer_strings, mapped_primer_colours, target_start - primer_design.config.FLANK)

        primer_design.pdf.method()

        primer_design.pdf.save()

    primer_design.core.delete_tmp_files()


if __name__ == '__main__':
    main()


