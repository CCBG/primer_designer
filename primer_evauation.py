#!/usr/bin/python 
# 
# 
# 
# 
# Kim Brugger (21 Oct 2015), contact: kbr@brugger.dk

import sys
import os
import pprint
pp = pprint.PrettyPrinter(indent=4)

import re

FLANK              = 500
NR_PRIMERS         = 4
ALLOWED_MISMATCHES = 4
MAX_MAPPINGS       = 5
MAX_PRODUCT_SIZE   = 2000
MIN_PRODUCT_SIZE   = 120
 

def check_primers( smalt_results ):


    id_word = "FULLSEQ"
    match = {}
    smalt_report = []
    query_region = []

    res = dict()

    with open(smalt_results, 'rU') as smalt_output:

        for line in smalt_output:
    
            if (line.startswith("@")):
                continue
       
            line = line.rstrip("\n")
            fields = line.split("\t")
#            pp.pprint( fields )
            match[ 'name'           ] = fields[  0  ] #mapping_score
            match[ 'chromosome'     ] = fields[  2  ] #chromosome
            match[ 'pos'            ] = fields[  3  ] #mapping position
            match[ 'length'         ] = len(fields[  9  ]) #mapping_score
            match[ 'length_matched' ] = int(re.sub("AS:i:", '', fields[  12 ])) #mapping_length


            match_id     = fields[ 0 ] 
            match_chr    = fields[ 2 ]
            match_pos    = fields[ 3 ]
            match_length = len(fields[ 9 ])
            match_mathes = int(re.sub("AS:i:", '', fields[  12 ]))

            if ( match_id not in res ):
                res[ match_id ] = dict()
                res[ match_id ][ 'CHR' ] = []
                res[ match_id ][ 'POS' ] = []



            if (match['length'] <= match['length_matched'] + ALLOWED_MISMATCHES):
                res[ match_id ][ 'CHR' ].append( match_chr )
                res[ match_id ][ 'POS' ].append( match_pos )



#            smalt_report.append()
    


    for primer  in  res.keys() :
        if (primer == 'FULLSEQ'):
            continue 

        res[ primer ]['MAPPING_SUMMARY'] = 'unique mapping'


        nr_of_chromosomes = len(set(res[ primer ][ 'CHR' ]))
        nr_of_mappings    = len( res[ primer ][ 'POS' ])

        if (nr_of_mappings > 1 and nr_of_mappings <= MAX_MAPPINGS ):
            res[ primer ]['MAPPING_SUMMARY'] = '%d mappings' % nr_of_mappings
            res[ primer ][ 'MAPPING_SUMMARY' ] += " to chromosomes: " + ",".join(res[ primer ][ 'CHR' ])


        elif (nr_of_mappings >= MAX_MAPPINGS ):
            res[ primer ]['MAPPING_SUMMARY'] = '%d mappings' % nr_of_mappings
            res[ primer ][ 'MAPPING_SUMMARY' ] += " on %d chromosomes" % nr_of_chromosomes



#    pp.pprint( smalt_report)
#    pp.pprint( res )


    return res




def digital_PCR( primer_mappings ):


    primer_names = sorted(primer_mappings.keys())
    nr_primer_names = len( primer_names )

    mappings = {}


    for i in range(0, nr_primer_names):
        primer1 = primer_names[ i ] 

        if ( primer1 == 'FULLSEQ'):
            continue

        if ( not re.search(r'LEFT', primer1 )):
            continue


        mappings[ primer1 ] = {}

        for j in range(0, nr_primer_names):
            primer2 = primer_names[ j ] 

            if ( primer2 == 'FULLSEQ'):
                continue

            if ( not re.search(r'RIGHT', primer2 )):
                continue


            mappings[ primer1 ][ primer2 ] = []


            print " -- %s vs %s"  % (primer1, primer2)


            for chr_index1 in range(0, len(primer_mappings[ primer1 ][ 'CHR' ])):
                for chr_index2 in range(0, len(primer_mappings[ primer2 ][ 'CHR' ])):

                    chr1 = primer_mappings[ primer1 ][ 'CHR' ][ chr_index1 ]
                    chr2 = primer_mappings[ primer2 ][ 'CHR' ][ chr_index2 ]

                    if ( chr1 != chr2 ):
                        continue


                    pos1 = int( primer_mappings[ primer1 ][ 'POS' ][ chr_index1 ] )
                    pos2 = int( primer_mappings[ primer2 ][ 'POS' ][ chr_index2 ] )

                    product_size = ( pos2 - pos1 )

                    if ( product_size > MAX_PRODUCT_SIZE ):
                        continue


#                    if ( product_size < 0 or 
#                         product_size > MAX_PRODUCT_SIZE ):
#                        continue

                    if ( product_size < MIN_PRODUCT_SIZE ):
                        continue

                    print " %s:%d vs %s:%d ==> %d" % ( chr1, pos1, chr2, pos2, product_size )

                    mappings[ primer1 ][ primer2 ].append( product_size )
                                               
        print "\n"

#    pp.pprint( mappings )

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

#            print "%s > %s " % (longest_product,  mappings[ primer1 ][ primer2 ][ 0 ])

            if ( mappings[ primer1 ][ primer2 ][0] > longest_product ):
                longest_product = mappings[ primer1 ][ primer2 ][ 0 ]
                longest_product_primer_pairs = ( primer1, primer2 )
#                print "%s > %s (post)" % (longest_product,  mappings[ primer1 ][ primer2 ][ 0 ])



    print "\n\nLongest product (%d bp) comes from the %s and %s primer pair" % (longest_product, 
                                                                            longest_product_primer_pairs[0], 
                                                                            longest_product_primer_pairs[1])

                




primer_data = check_primers( '8:96259936.smalt' )

#pp.pprint( primer_data )


digital_PCR( primer_data )
