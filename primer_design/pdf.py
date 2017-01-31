"""
   configuration of location of 3rd party software and databases and design paramters.

"""


from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfbase.pdfmetrics import stringWidth 
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import A4
width, height = A4 #keep for later

import core
import format


if ( __name__ == '__main__'):
    sys.stderr.write( "This is a module and not meant to be run as a stand alone program\n" )
    exit( 1 )


c = None # the canvas we will be creating on.
top_offset = 0

def create_canvas( filename, pagesize=A4 ):

    global c
    c = canvas.Canvas( filename , pagesize=A4)

    return c


def save_canvas():

    global c

    c.showPage()
    c.save()



def set_top_offset( new_offset ):
    assert new_offset.isdigit(), "Int not provided as arguement"

    global top_offset
    top_offset = new_offset
    

colours = [[255,   0,   0], # red
           [  0, 255,   0], # green
           [  0,   0, 255], # blue
           [255,   0, 255], # Pink
           [0,   255, 255],
           [255, 255,   0], #Dark Green
           [100, 255, 100]] # Yellow, crap!


def mappings(target_sequence, tagged_string, primer_strings, primer_colours, base1):
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

    core.verbose_print("pdf.mappings", 2)

    global top_offset
    c.setFont('mono', 8)

    lines = []
    
    for i in range(0, len(target_sequence), 80):

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


def primer_data(target_chrom, target_start, target_end, primer_dict, fwd_primer=None, rev_primer=None ):

    global top_offset


    font = TTFont('mono', '/usr/share/fonts/truetype/ttf-liberation/LiberationMono-Regular.ttf')
    pdfmetrics.registerFont( font )
    c.setFont('mono', 7)

    core.verbose_print("pretty_pdf_primer_data", 2)

#    c.drawString(40 , top_offset, "_-=-"*15 +"_" )
    c.line(40,top_offset,width - 40  ,top_offset+2)
    top_offset -= 8


    if ( target_start == target_end ):
        c.drawString(40 , top_offset, "Primer design report for chr: %s position: %d" % (target_chrom, target_start))
    else:
        c.drawString(40 , top_offset, "Primer design report for chr: %s range: %d-%d" % (target_chrom, 
                                                                                       target_start, target_end))
            
    top_offset -= 8
    c.line(40,top_offset,width - 40 ,top_offset+2)
    top_offset -= 16

    c.drawString(40 , top_offset, "ID         POS             %GC    TM     Primer sequence     Best primer  Bad pairs    ")
    top_offset -= 8
    c.line(40,top_offset,width - 40 ,top_offset+2)
    top_offset -= 8




    primer_seqs = []
    for primer in sorted(primer_dict):
        if ( primer == 'FULLSEQ'):
            continue


        name = primer
        name = re.sub(r'PRIMER_', '', name)
        name = re.sub(r'_SEQUENCE', '', name)

        
        if (name == "RIGHT_0"):
            top_offset -= 8

        picked_primer = ' '

        if ( fwd_primer == primer or rev_primer == primer):
            picked_primer = 'Y'

        c.drawString(40 , top_offset, "%-10s %-14s  %.2f  %.2f  %-25s %s      %s" % ("", 
                                    "%s:%d" % (primer_dict[ name ][ "TARGET_CHROM"], primer_dict[ name ][ "TARGET_POS"]),
                                    primer_dict[ name ][ "GC_PERCENT"], 
                                    primer_dict[ name ][ "TM"],
                                    primer_dict[ name ][ "SEQ"], picked_primer,
                                    primer_dict[ name ].get( 'MULTIPLE_MAPPINGS', '')))
#                                    primer_dict[ name ][ 'MAPPING_SUMMARY' ]))



        primer_nr = re.sub(r'.*_(\d)',r'\1' , name)

#        pp.pprint( primer_nr )
#        print pp.pprint(colours[ int( primer_nr ) ])
        

        c.setFillColorRGB(colours[ int( primer_nr ) ][0], 
                          colours[ int( primer_nr ) ][1],
                          colours[ int( primer_nr ) ][2])

        c.drawString(40 , top_offset, name)
        c.setFillColorRGB(0,0,0)
        top_offset -= 8


            



    top_offset -= 8
    c.line(40,top_offset,width - 40 ,top_offset+2)
    top_offset -= 8
    top_offset -= 8


    return top_offset

def method(c, top_offset):
    core.verbose_print("pretty_pdf_method", 3)

    lines = format.method_blurb()

    top_offset = 80

    for line in lines:
        c.drawString(40, top_offset, line)
        top_offset -= 8
