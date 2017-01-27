"""
  blastn related functions

"""

import sys
import shlex
import subprocess
import re

import config
import core


def map_primers( region_id, primer3_dict, target_chrom, target_start, target_end):
    """
    Maps primers to the reference using smalt to check where they map.
    
    input: region-ID, primer3_dict
    
    output: dict of primers and their nr of mappings to the genome

    Kim Brugger (17 Aug 2016)

    """
    core.verbose_print( "check_primers", 2)

    res = dict()

    primers_file = region_id + "_p3seq.fasta"

    core.add_new_tmp_file( primers_file )

    with open( primers_file, 'w+') as primerfasta:

        for i in range(0, 5):

            primer_left_id = "PRIMER_LEFT_%d_SEQUENCE" % i
            primer_right_id = "PRIMER_RIGHT_%d_SEQUENCE" % i

            if ( primer_right_id not in primer3_dict or primer_left_id not in primer3_dict):
                continue

            primerfasta.write(">%s\n%s\n" % (primer_left_id, primer3_dict[ primer_left_id ]))
            primerfasta.write(">%s\n%s\n" % (primer_right_id, primer3_dict[ primer_right_id ]))

            # Store the primer seqs in the res dict as we will needthem later in the program
            left_id = "LEFT_%d" % i
            res[ left_id ] = {}
            res[ left_id ][ 'SEQ' ] = primer3_dict[ primer_left_id ]
            res[ left_id ][ 'TM' ]  = float(primer3_dict[ "PRIMER_LEFT_%d_TM" % i ])
            res[ left_id ][ 'GC_PERCENT' ]  = float(primer3_dict[ "PRIMER_LEFT_%d_GC_PERCENT" % i ])
            res[ left_id ][ 'CHR' ] = []
            res[ left_id ][ 'POS' ] = []
            res[ left_id ][ 'STRAND' ] = []

            right_id = "RIGHT_%d" % i
            res[ right_id ] = {}
            res[ right_id ][ 'SEQ' ] = primer3_dict[ primer_right_id ]
            res[ right_id ][ 'TM' ]  = float(primer3_dict[ "PRIMER_RIGHT_%d_TM" % i ])
            res[ right_id ][ 'GC_PERCENT' ]  = float(primer3_dict[ "PRIMER_RIGHT_%d_GC_PERCENT" % i ])
            res[ right_id ][ 'CHR' ] = []
            res[ right_id ][ 'POS' ] = []
            res[ right_id ][ 'STRAND' ] = []


        primerfasta.close()

    blastn_results = region_id + ".blastn"
#    pp.pprint( res )
    core.add_new_tmp_file( blastn_results )

#    cmd = SMALT + " map  -d -1 /refs/human_1kg/human_g1k_v37 " + primers_file + "> " + smalt_results + " 2> /dev/null"

    cmd = "{blastn} -db {blast_db} -word_size 7 -outfmt 7 -query {query} -out {outfile}".format(blastn   = config.BLASTN, 
                                                                                                blast_db = config.BLAST_DB, 
                                                                                                query    = primers_file,
                                                                                                outfile  = blastn_results)

    core.verbose_print( cmd, 4)
    subprocess.call(cmd, shell=True)


#    pp.pprint( res )

    match = {}
    smalt_report = []
    query_region = []


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
            align_end       = int( fields[ 7 ])
            
            match_disagree  = int( fields[ 4 ])
            match_gaps      = int( fields[ 5 ])
            match_bp_agree  = seq_length - ( seq_length - match_length) - match_disagree 

            match_strand    = "plus"
            if ( int( fields[ 8 ] ) > int( fields[ 9 ])):
                match_strand    = "minus"
                match_pos       = int(fields[ 9 ])


            # Check that the 3' end of the primer looks good
            if ( seq_length > align_end + MAX_3_PRIME_MISMATCH ):
                continue


            if ( match_gaps ):
                continue 


#        print " %d -- %d " % ( seq_length, match_bp_agree )


            if (seq_length <= match_bp_agree + ALLOWED_MISMATCHES):
                res[ match_id ][ 'CHR' ].append( match_chrom )
                res[ match_id ][ 'POS' ].append( match_pos )
                res[ match_id ][ 'STRAND'].append( match_strand )



            if ( target_chrom == match_chrom and 
                 target_start < match_pos and
                 target_end  > match_pos ):
                res[ match_id ][ 'TARGET_POS'   ] = match_pos
                res[ match_id ][ 'TARGET_CHROM' ] = match_chrom
            

#    pp.pprint( res )

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





    return res

