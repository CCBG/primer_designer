#!/usr/bin/python
# 
# 
# 
# 
# Kim Brugger (24 Jan 2017), contact: kim@brugger.dk

import sys
import pprint
pp = pprint.PrettyPrinter(indent=4)

ALLOWED_MISMATCHES = 4   # Maximum nr of errors when mapping a primer back to the reference
MAX_MAPPINGS       = 5   # Less or equal number of mappings to the reference what chromosomes mapped to are named 

target_chrom = '9' 
target_start = 129455552 - 1000
target_end = 129455552 + 1000


import re

res = {   'LEFT_0': {   'CHR': [],
                  'GC_PERCENT': 60.0,
                  'POS': [],
                  'SEQ': 'CAGGTGTCAACAGAGGGGAC',
                  'TM': 59.965},
    'LEFT_1': {   'CHR': [],
                  'GC_PERCENT': 55.0,
                  'POS': [],
                  'SEQ': 'ACCTGTGTGAGTTCCTGTGC',
                  'TM': 60.179},
    'LEFT_2': {   'CHR': [],
                  'GC_PERCENT': 63.158,
                  'POS': [],
                  'SEQ': 'GTCAACAGAGGGGACAGGC',
                  'TM': 60.003},
    'LEFT_3': {   'CHR': [],
                  'GC_PERCENT': 55.0,
                  'POS': [],
                  'SEQ': 'TCCTGAAAGTGTGTGCAGCC',
                  'TM': 60.819},
    'LEFT_4': {   'CHR': [],
                  'GC_PERCENT': 57.895,
                  'POS': [],
                  'SEQ': 'ACGGCAGGTGTCAACAGAG',
                  'TM': 59.93},
    'RIGHT_0': {   'CHR': [],
                   'GC_PERCENT': 60.0,
                   'POS': [],
                   'SEQ': 'GTGCGGAGAGATGGAGGATG',
                   'TM': 59.967},
    'RIGHT_1': {   'CHR': [],
                   'GC_PERCENT': 60.0,
                   'POS': [],
                   'SEQ': 'CTCAGCTGCCAGTGTCTCTC',
                   'TM': 60.109},
    'RIGHT_2': {   'CHR': [],
                   'GC_PERCENT': 55.0,
                   'POS': [],
                   'SEQ': 'GGTGGCCTCTTACCTTTGCT',
                   'TM': 59.962},
    'RIGHT_3': {   'CHR': [],
                   'GC_PERCENT': 60.0,
                   'POS': [],
                   'SEQ': 'GGAGATGGTGGGTTTAGGGG',
                   'TM': 59.449},
    'RIGHT_4': {   'CHR': [],
                   'GC_PERCENT': 55.0,
                   'POS': [],
                   'SEQ': 'GCCAGCTTCTTCATCTGCAG',
                   'TM': 59.266}}



#pp.pprint( res )

blastn_results = "9:129455552.blastn"

with open(blastn_results, 'rU') as blastn_output:

    for line in blastn_output:
        
        if (line.startswith("#")):
            continue

#        print line
        
        line = line.rstrip("\n")
        fields = line.split("\t")
        
        match_id        = fields[ 0 ] 
        match_id        = re.sub(r'PRIMER_', '', match_id)
        match_id        = re.sub(r'_SEQUENCE', '', match_id)
        
        match_chrom     = fields[ 1 ]
        match_pos       = int(fields[ 8 ])
        
        seq             = res[ match_id ][ 'SEQ' ]
        seq_length      = len( seq )
        match_length    = int( fields[ 3 ])
        match_disagree  = int( fields[ 4 ])
        match_gaps      = int( fields[ 5 ])

        if ( match_gaps ):
            continue 

        match_bp_agree  = seq_length - ( seq_length - match_length) - match_disagree 


        if (seq_length   <= match_bp_agree + ALLOWED_MISMATCHES):
            res[ match_id ][ 'CHR' ].append( match_chrom )
            res[ match_id ][ 'POS' ].append( match_pos )

        print " %d -- %d " % ( seq_length, match_bp_agree )

        if ( target_chrom == match_chrom and 
             target_start < match_pos and
             target_end  > match_pos ):
            res[ match_id ][ 'TARGET_POS'   ] = match_pos
            res[ match_id ][ 'TARGET_CHROM' ] = match_chrom
            

pp.pprint( res )

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
        

pp.pprint( res )
