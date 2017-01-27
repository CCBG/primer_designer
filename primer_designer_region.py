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


import primer_design_module as pd
import primer_design.core
import primer_design.config
import primer_design.primer3
import primer_design.blast
import primer_design.printing
import primer_design.format



def pretty_pdf_mappings(c, top_offset,  target_sequence, tagged_string, primer_strings, primer_colours, base1):
    """

    makes a pretty text output of target region and primers in a pdf canvas
    
    input: canvas, 
           offset (where to start writing) 
           target sequence
           target snp/target string, 
           mapped primers,
           colours to apply to the primer arrows
           position of first base of the reference region
           
           output: string with nicely formatted data
           
    """

    primer_design.core.verbose_print("pretty_pdf_mappings", 2)

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

                    c.setFillColorRGB(pd.colours[ primer_colour[k] ][0], 
                                      pd.colours[ primer_colour[k] ][1],
                                      pd.colours[ primer_colour[k] ][2])

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
def pretty_pdf_primer_data(c, y_offset, target_chrom, target_start, target_end, primer_dict, fwd_primer=None, rev_primer=None ):
    primer_design.core.verbose_print("pretty_pdf_primer_data", 2)

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

    c.drawString(40 , y_offset, "ID         POS             %GC    TM     Primer sequence           Mapping(s)    ")
    y_offset -= 8
    c.line(40,y_offset,width - 40 ,y_offset+2)
    y_offset -= 8




    primer_seqs = []
    for primer in sorted(primer_dict):
        if ( primer == 'FULLSEQ'):
            continue


        name = primer
        name = re.sub(r'PRIMER_', '', name)
        name = re.sub(r'_SEQUENCE', '', name)

        
        if (name == "RIGHT_0"):
            y_offset -= 8

        picked_primer = ' '

        if ( fwd_primer == primer or rev_primer == primer):
            picked_primer = 'Y'

        c.drawString(40 , y_offset, "%-10s %-14s  %.2f  %.2f  %-25s %s            %s" % ("", 
                                    "%s:%d" % (primer_dict[ name ][ "TARGET_CHROM"], primer_dict[ name ][ "TARGET_POS"]),
                                    primer_dict[ name ][ "GC_PERCENT"], 
                                    primer_dict[ name ][ "TM"],
                                    primer_dict[ name ][ "SEQ"], picked_primer,
                                    primer_dict[ name ][ 'MAPPING_SUMMARY' ]))



        primer_nr = re.sub(r'.*_(\d)',r'\1' , name)

#        pp.pprint( primer_nr )
#        print pp.pprint(colours[ int( primer_nr ) ])
        

        c.setFillColorRGB(pd.colours[ int( primer_nr ) ][0], 
                          pd.colours[ int( primer_nr ) ][1],
                          pd.colours[ int( primer_nr ) ][2])

        c.drawString(40 , y_offset, name)
        c.setFillColorRGB(0,0,0)
        y_offset -= 8


            



    y_offset -= 8
    c.line(40,y_offset,width - 40 ,y_offset+2)
    y_offset -= 8
    y_offset -= 8


    return y_offset

def pretty_pdf_method(c, top_offset):
    primer_design.core.verbose_print("pretty_pdf_method", 3)

    lines = primer_design.format.method_blurb()

    top_offset = 80

    for line in lines:
        c.drawString(40, top_offset, line)
        top_offset -= 8





def get_and_parse_arguments():
    primer_design.core.verbose_print( "get_and_parse_arguments", 4)

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('-c', '--chrom')
    
    # NICK ADDED RANGE
    
    group = parser.add_mutually_exclusive_group( required = True )
    group.add_argument('-p', '--pos')
    group.add_argument('-r', '--range', nargs = 2)

    parser.add_argument('-o', '--output')
    parser.add_argument('-f', '--flank')
    parser.add_argument('-t', '--text_output', action="store_true")
    parser.add_argument('-w', '--website_output', action="store_true")

    args = parser.parse_args()



    if ( args.range ):
        args.target_start, args.target_end = [int(x) for x in args.range]
    else:
        args.target_start = int( args.pos )
        args.target_end   = int( args.pos )

    # how many bp on either side are to be used for designing primers in
    args.target_flank = args.flank or primer_design.config.FLANK
    args.target_flank = int( primer_design.config.FLANK )
    
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

    target_chrom    = args.chrom
    target_start    = args.target_start
    target_end      = args.target_end
    target_flank    = args.target_flank
    target_id       = args.target_id
    target_filename = args.filename

    ( exon_chrom, exon_start, exon_end) = pd.map_target_region_to_exon( target_chrom, target_start, target_end)

    if ( exon_start is not None and  exon_end is not None  and exon_end - exon_start + 1 < primer_design.config.MAX_PRODUCT):
        primer_design.core.verbose_print("Changing region to an exonic region (%s:%d-%d)" % (exon_chrom, exon_start, exon_end), 1)
        target_start = exon_start
        target_end   = exon_end

#    exit()


    target_sequence                   =   primer_design.core.fetch_region( target_chrom, target_start - target_flank, target_end + target_flank     )
    (tagged_sequence, tagged_region)  =   pd.markup_sequence(target_sequence, target_chrom, target_start, target_end, target_flank )

    primer3_results  = primer_design.primer3.run( target_id , tagged_sequence, "%s_%d.primer3" % ( target_chrom, target_start))
#    primer_dict   = map_primers_smalt( target_id, primer3_results, target_chrom, target_start - target_flank, target_end+ target_flank)
    primer_dict   = primer_design.blast.map_primers( target_id, primer3_results, target_chrom, target_start - target_flank, target_end+ target_flank)
#    pp.pprint( primer_dict )
    

    best_fwd_primer, best_rev_primer = pd.pick_best_primers(primer_dict, target_chrom, target_start, target_end)
    (mapped_primer_strings, mapped_primer_colours) = pd.make_mapped_primer_strings( target_sequence, primer_dict)

#    pp.pprint( mapped_primer_strings )
#    pp.pprint( mapped_primer_colours )



    if (  args.website_output ):

        lines = primer_design.format.tabbed_primer_data(target_chrom, target_start, target_end, primer_dict, best_fwd_primer, best_rev_primer )

        primer_design.core.print_to_file_or_stdout( lines )

    elif (  args.text_output ):
        lines  = pd.pretty_print_primer_data( primer_dict, target_chrom, target_start, target_end )
        lines += pd.pretty_print_mappings( target_sequence, tagged_region, mapped_primer_strings, target_start - primer_design.config.FLANK)
        pd.pretty_primer_data("%s.txt" % target_filename, target_chrom, target_start, target_end, primer_dict, best_fwd_primer, best_rev_primer )
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
        pretty_pdf_mappings(c, top_offset,  target_sequence, tagged_region, mapped_primer_strings, mapped_primer_colours, target_start - primer_design.config.FLANK)

        pretty_pdf_method(c, top_offset)

        c.showPage()
        c.save()

    primer_design.core.delete_tmp_files()


if __name__ == '__main__':
    main()


