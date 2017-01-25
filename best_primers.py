#!/usr/bin/python
# 
# 
# 
# 
# Kim Brugger (25 Jan 2017), contact: kim@brugger.dk

import sys
import pprint
pp = pprint.PrettyPrinter(indent=4)
VERBOSE    =  20

MAX_PRODUCT_SIZE   = 2000
MIN_PRODUCT_SIZE   = 120

import re

def verbose_print( msg, level ):
    if ( level <= VERBOSE ):
        sys.stderr.write( msg+ "\n" )




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
            pp.pprint( primer_data[ primer ] )
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



def digital_PCR( primer_mappings ):

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

    

primer_data = {   'LEFT_0': {   'CHR': ['9', '17', '17', '12', '11', '11'],
                  'GC_PERCENT': 60.0,
                  'MAPPING_SUMMARY': '6 mappings on 4 chromosomes',
                  'POS': [   129455306,
                             79131437,
                             71189902,
                             41496819,
                             115091323,
                             129480268],
                  'SEQ': 'CAGGTGTCAACAGAGGGGAC',
                  'STRAND': [   'plus',
                                'minus',
                                'minus',
                                'plus',
                                'minus',
                                'plus'],
                  'TARGET_CHROM': '9',
                  'TARGET_POS': 129455306,
                  'TM': 59.965},
    'LEFT_1': {   'CHR': [   '9',
                             '9',
                             '9',
                             '17',
                             '17',
                             '11',
                             '11',
                             '10',
                             'X',
                             '15',
                             '14',
                             '3',
                             '1'],
                  'GC_PERCENT': 55.0,
                  'MAPPING_SUMMARY': '13 mappings on 9 chromosomes',
                  'POS': [   129455129,
                             136520510,
                             123213546,
                             3419080,
                             3970709,
                             55472574,
                             115008159,
                             28383150,
                             46303535,
                             59270367,
                             105593336,
                             131435105,
                             196026254],
                  'SEQ': 'ACCTGTGTGAGTTCCTGTGC',
                  'STRAND': [   'plus',
                                'minus',
                                'minus',
                                'minus',
                                'plus',
                                'plus',
                                'plus',
                                'minus',
                                'minus',
                                'minus',
                                'minus',
                                'plus',
                                'minus'],
                  'TARGET_CHROM': '9',
                  'TARGET_POS': 129455129,
                  'TM': 60.179},
    'LEFT_2': {   'CHR': [   '9',
                             '14',
                             '6',
                             '5',
                             '3',
                             '17',
                             '16',
                             '15',
                             '11',
                             '2',
                             '1',
                             '1',
                             '1'],
                  'GC_PERCENT': 63.158,
                  'MAPPING_SUMMARY': '13 mappings on 11 chromosomes',
                  'POS': [   129455311,
                             82323837,
                             32997706,
                             77256444,
                             137057871,
                             27523350,
                             27737897,
                             50137329,
                             47708880,
                             33258342,
                             47699543,
                             177862031,
                             235325197],
                  'SEQ': 'GTCAACAGAGGGGACAGGC',
                  'STRAND': [   'plus',
                                'plus',
                                'minus',
                                'plus',
                                'minus',
                                'plus',
                                'minus',
                                'minus',
                                'plus',
                                'minus',
                                'minus',
                                'plus',
                                'minus'],
                  'TARGET_CHROM': '9',
                  'TARGET_POS': 129455311,
                  'TM': 60.003},
    'LEFT_3': {   'CHR': [   '9',
                             '21',
                             '21',
                             '2',
                             '2',
                             '1',
                             '15',
                             '13',
                             '12',
                             '8',
                             '5',
                             '3'],
                  'GC_PERCENT': 55.0,
                  'MAPPING_SUMMARY': '12 mappings on 10 chromosomes',
                  'POS': [   129455065,
                             43873234,
                             43874606,
                             43492074,
                             8299938,
                             166323994,
                             40719530,
                             45866321,
                             27632328,
                             134281606,
                             102031698,
                             22821564],
                  'SEQ': 'TCCTGAAAGTGTGTGCAGCC',
                  'STRAND': [   'plus',
                                'minus',
                                'plus',
                                'minus',
                                'minus',
                                'minus',
                                'minus',
                                'minus',
                                'plus',
                                'minus',
                                'minus',
                                'minus'],
                  'TARGET_CHROM': '9',
                  'TARGET_POS': 129455065,
                  'TM': 60.819},
    'LEFT_4': {   'CHR': ['9', '22', '2', '1'],
                  'GC_PERCENT': 57.895,
                  'MAPPING_SUMMARY': '4 mappings to chromosome: 9,1,2,22',
                  'POS': [129455302, 51076726, 240037323, 20211934],
                  'SEQ': 'ACGGCAGGTGTCAACAGAG',
                  'STRAND': ['plus', 'minus', 'plus', 'plus'],
                  'TARGET_CHROM': '9',
                  'TARGET_POS': 129455302,
                  'TM': 59.93},
    'RIGHT_0': {   'CHR': [   '9',
                              '1',
                              '1',
                              '1',
                              '19',
                              '19',
                              '17',
                              '17',
                              '16',
                              '13',
                              '7',
                              '4',
                              '3',
                              '2'],
                   'GC_PERCENT': 60.0,
                   'MAPPING_SUMMARY': '14 mappings on 10 chromosomes',
                   'POS': [   129455703,
                              41996756,
                              63492662,
                              235258368,
                              3678436,
                              18404993,
                              7341246,
                              25747557,
                              75264729,
                              40115311,
                              23280024,
                              16618029,
                              177680505,
                              232245672],
                   'SEQ': 'GTGCGGAGAGATGGAGGATG',
                   'STRAND': [   'minus',
                                 'plus',
                                 'minus',
                                 'minus',
                                 'plus',
                                 'plus',
                                 'minus',
                                 'plus',
                                 'minus',
                                 'minus',
                                 'minus',
                                 'plus',
                                 'minus',
                                 'minus'],
                   'TARGET_CHROM': '9',
                   'TARGET_POS': 129455703,
                   'TM': 59.967},
    'RIGHT_1': {   'CHR': [   '9',
                              '9',
                              '9',
                              '10',
                              '10',
                              '4',
                              '4',
                              '12',
                              '11',
                              '6',
                              '5',
                              '3',
                              '2',
                              '2'],
                   'GC_PERCENT': 55.0,
                   'MAPPING_SUMMARY': '14 mappings on 9 chromosomes',
                   'POS': [   129455874,
                              97133586,
                              99779108,
                              75555227,
                              122954297,
                              33599871,
                              20717865,
                              2411792,
                              123659962,
                              69907659,
                              59748904,
                              129301523,
                              15901772,
                              76440074],
                   'SEQ': 'GGTGGCCTCTTACCTTTGCT',
                   'STRAND': [   'minus',
                                 'plus',
                                 'minus',
                                 'plus',
                                 'plus',
                                 'plus',
                                 'minus',
                                 'plus',
                                 'minus',
                                 'minus',
                                 'plus',
                                 'minus',
                                 'minus',
                                 'minus'],
                   'TARGET_CHROM': '9',
                   'TARGET_POS': 129455874,
                   'TM': 59.962},
    'RIGHT_2': {   'CHR': ['9', '12', '18', '4', 'X', '6', '5', '2'],
                   'GC_PERCENT': 60.0,
                   'MAPPING_SUMMARY': '8 mappings on 8 chromosomes',
                   'POS': [   129455807,
                              27446513,
                              19078760,
                              54771932,
                              152833167,
                              16426884,
                              146337796,
                              29760559],
                   'SEQ': 'CTCAGCTGCCAGTGTCTCTC',
                   'STRAND': [   'minus',
                                 'minus',
                                 'minus',
                                 'minus',
                                 'minus',
                                 'minus',
                                 'plus',
                                 'minus'],
                   'TARGET_CHROM': '9',
                   'TARGET_POS': 129455807,
                   'TM': 60.109},
    'RIGHT_3': {   'CHR': ['9', '9', '15', '8', '4'],
                   'GC_PERCENT': 60.0,
                   'MAPPING_SUMMARY': '5 mappings to chromosome: 9,8,15,4',
                   'POS': [129455734, 130400726, 92275519, 27586036, 2910649],
                   'SEQ': 'GGAGATGGTGGGTTTAGGGG',
                   'STRAND': ['minus', 'minus', 'plus', 'minus', 'plus'],
                   'TARGET_CHROM': '9',
                   'TARGET_POS': 129455734,
                   'TM': 59.449},
   'RIGHT_4': {   'CHR': [   '9',
                              '9',
                              '9',
                              '9',
                              '13',
                              '13',
                              '3',
                              '3',
                              '3',
                              '3',
                              '3',
                              '1',
                              '1',
                              '1',
                              '1',
                              '1',
                              '1',
                              '22',
                              '22',
                              '22',
                              '20',
                              '16',
                              '16',
                              '16',
                              '16',
                              '15',
                              '15',
                              '15',
                              '14',
                              '14',
                              '14',
                              '14',
                              '12',
                              '11',
                              '11',
                              '10',
                              '10',
                              '8',
                              '8',
                              '8',
                              '8',
                              '7',
                              '7',
                              '7',
                              '6',
                              '6',
                              '6',
                              '6',
                              'Y',
                              'X',
                              'X',
                              '21',
                              '19',
                              '18',
                              '18',
                              '17',
                              '5',
                              '5',
                              '5',
                              '5',
                              '4',
                              '2',
                              '2',
                              '2'],
                   'GC_PERCENT': 55.0,
                   'MAPPING_SUMMARY': '64 mappings on 24 chromosomes',
                   'POS': [   129456019,
                              3453034,
                              101979053,
                              80546688,
                              84132268,
                              22510191,
                              145759075,
                              184401028,
                              66805049,
                              130296557,
                              166867169,
                              49189929,
                              5151140,
                              68770303,
                              90073711,
                              179205534,
                              233404291,
                              33369741,
                              27971982,
                              44280329,
                              38158921,
                              14917119,
                              16315882,
                              18583989,
                              57082031,
                              39786412,
                              58190797,
                              70398763,
                              72768821,
                              100102467,
                              34800836,
                              77953735,
                              110338839,
                              76542920,
                              46329063,
                              105429720,
                              6209943,
                              28244873,
                              108868256,
                              87775046,
                              90691370,
                              31091264,
                              67228655,
                              13579076,
                              28552555,
                              37601684,
                              45993146,
                              111011498,
                              16952885,
                              5811097,
                              93498534,
                              15480514,
                              16551846,
                              5795214,
                              35335340,
                              49411512,
                              35435828,
                              56856724,
                              80390492,
                              173251820,
                              134304661,
                              36331182,
                              85967079,
                              121617924],
                   'SEQ': 'GCCAGCTTCTTCATCTGCAG',
                   'STRAND': [   'minus',
                                 'plus',
                                 'minus',
                                 'plus',
                                 'minus',
                                 'minus',
                                 'plus',
                                 'minus',
                                 'minus',
                                 'plus',
                                 'plus',
                                 'minus',
                                 'minus',
                                 'minus',
                                 'plus',
                                 'plus',
                                 'minus',
                                 'plus',
                                 'minus',
                                 'minus',
                                 'minus',
                                 'minus',
                                 'minus',
                                 'plus',
                                 'plus',
                                 'minus',
                                 'minus',
                                 'minus',
                                 'plus',
                                 'plus',
                                 'minus',
                                 'plus',
                                 'plus',
                                 'minus',
                                 'minus',
                                 'minus',
                                 'plus',
                                 'plus',
                                 'minus',
                                 'minus',
                                 'minus',
                                 'minus',
                                 'plus',
                                 'plus',
                                 'plus',
                                 'plus',
                                 'minus',
                                 'minus',
                                 'minus',
                                 'plus',
                                 'plus',
                                 'plus',
                                 'minus',
                                 'minus',
                                 'plus',
                                 'minus',
                                 'plus',
                                 'minus',
                                 'plus',
                                 'plus',
                                 'plus',
                                 'plus',
                                 'minus',
                                 'plus'],
                   'TARGET_CHROM': '9',
                   'TARGET_POS': 129456019,
                   'TM': 59.266}}

pick_best_primers( primer_data, '9', 129455552, 129455552)
digital_PCR( primer_data )
