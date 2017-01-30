"""

Various functions to evaluate the mappings of primers to get the best primer-pair.

"""
import pprint as pp

import core


if ( __name__ == '__main__'):
    sys.stderr.write( "This is a module and not meant to be run as a stand alone program\n" )
    exit( 1 )




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



def digital_PCR( primer_mappings ):
    """
    Makes a "digital" PCR by looking at the mappings of primers and
    predict which will produce products, and more important multiple
    products
    """

    primer_names = sorted(primer_mappings.keys())
    nr_primer_names = len( primer_names )

    mappings = {}

    products = {}

    for i in range(0, nr_primer_names):
        primer1 = primer_names[ i ] 

        if ( primer1 == 'FULLSEQ'):
            continue

        if ( not re.search(r'LEFT', primer1 )):
            continue

        mappings[ primer1 ] = {}
        products[ primer1 ] = {}

        for j in range(0, nr_primer_names):
            primer2 = primer_names[ j ] 

            if ( primer2 == 'FULLSEQ'):
                continue

            if ( not re.search(r'RIGHT', primer2 )):
                continue


            mappings[ primer1 ][ primer2 ] = []
            products[ primer1 ][ primer2 ] = []


#            print " -- %s vs %s"  % (primer1, primer2)

            for chr_index1 in range(0, len(primer_mappings[ primer1 ][ 'CHR' ])):
                for chr_index2 in range(0, len(primer_mappings[ primer2 ][ 'CHR' ])):

                    chr1 = primer_mappings[ primer1 ][ 'CHR' ][ chr_index1 ]
                    chr2 = primer_mappings[ primer2 ][ 'CHR' ][ chr_index2 ]


                    pos1 = int( primer_mappings[ primer1 ][ 'POS' ][ chr_index1 ] )
                    pos2 = int( primer_mappings[ primer2 ][ 'POS' ][ chr_index2 ] )

                    strand1 = primer_mappings[ primer1 ][ 'STRAND' ][ chr_index1 ]
                    strand2 = primer_mappings[ primer2 ][ 'STRAND' ][ chr_index2 ]

                    # The primers map to different chromosomes
                    if ( chr1 != chr2 ):
                        continue

                    # the primer are on the same strand. 
                    if ( strand1 == strand2 ):
                        continue


                    if ( pos1 < pos2 and strand1 != 'plus' and strand2 != 'minus'):
                        continue
                    elif( pos1 > pos2 and strand1 != 'minus' and strand2 != 'plus'):
                        continue


                    product_size = ( pos2 - pos1 )
                    if ( product_size < 0  or product_size > MAX_PRODUCT_SIZE):
                        continue


                    print "%s -- %s %s:%d:%s -> %s:%d:%s ==>> %d bp" %( primer1, primer2, chr1, pos1, strand1, chr2, pos2, strand2, product_size)

                    mappings[ primer1 ][ primer2 ].append( product_size )
                    products[ primer1 ][ primer2 ].append( {'chr' : chr1, 'start_pos': pos1, 'end_pos': pos2, 'size': product_size} )




    pp.pprint( products )
#    return products
#    exit()

    longest_product = 0
    longest_product_primer_pairs = ()
    for primer1 in mappings.keys():
        for primer2 in mappings[ primer1 ].keys():
        
            if ( len( mappings[ primer1 ][ primer2 ]) == 0 ):
                print "No usable pcr product from %s and %s" % ( primer1, primer2 )
                continue
            elif ( len( mappings[ primer1 ][ primer2 ]) > 1 ):
                print "multiple pcr products from %s and %s" % ( primer1, primer2 )
                continue

            print "%s + %s > %s bp" % (primer1, primer2,  mappings[ primer1 ][ primer2 ][ 0 ])

            if ( mappings[ primer1 ][ primer2 ][0] > longest_product ):
                longest_product = mappings[ primer1 ][ primer2 ][ 0 ]
                longest_product_primer_pairs = ( primer1, primer2 )
#                print "%s > %s (post)" % (longest_product,  mappings[ primer1 ][ primer2 ][ 0 ])



    print "\n\nLongest product (%d bp) comes from the %s and %s primer pair" % (longest_product, 
                                                                            longest_product_primer_pairs[0], 
                                                                            longest_product_primer_pairs[1])

    
