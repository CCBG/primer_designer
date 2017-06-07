"""
  blastn related functions

"""

import sys
import shlex
import subprocess
import re

import config
import core



if ( __name__ == '__main__'):
    sys.stderr.write( "This is a module and not meant to be run as a stand alone program\n" )
    exit( 1 )


def map_primers( region_id, primer3_dict, target_chrom, target_start, target_end, max_hits=None):
    """
    Maps primers to the reference using smalt to check where they map.
    
    input: region-ID, primer3_dict
    
    output: dict of primers and their nr of mappings to the genome

    Kim Brugger (17 Aug 2016)

    """
    core.verbose_print( "check_primers", 2)

    res = dict()

    primers_file = core.tmpfilename(path = '/tmp/', postfix='_p3.fasta')
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

    blastn_results = core.tmpfilename(path = '/tmp/', postfix='.blastn')
#    pp.pprint( res )
    core.add_new_tmp_file( blastn_results )

#    cmd = SMALT + " map  -d -1 /refs/human_1kg/human_g1k_v37 " + primers_file + "> " + smalt_results + " 2> /dev/null"

    cmd = "{blastn} -db {blast_db}   -word_size 7 -outfmt 7 -query {query} -out {outfile}".format(blastn   = config.BLASTN, 
                                                                                                blast_db = config.BLAST_DB, 
                                                                                                query    = primers_file,
                                                                                                outfile  = blastn_results)

    if max_hits is not None:
        cmd += " -max_target_seqs {}".format(  max_hits  )

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
            if ( seq_length > align_end + config.MAX_3_PRIME_MISMATCH ):
                continue


            if ( match_gaps ):
                continue 


#        print " %d -- %d " % ( seq_length, match_bp_agree )


            if (seq_length <= match_bp_agree + config.ALLOWED_MISMATCHES):
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
    
        if (nr_of_mappings > 1 and nr_of_mappings <= config.MAX_MAPPINGS ):
            res[ primer ]['MAPPING_SUMMARY'] = '%d mappings' % nr_of_mappings
            res[ primer ][ 'MAPPING_SUMMARY' ] += " to chromosome: " + ",".join(set(res[ primer ][ 'CHR' ]))
        
        
        elif (nr_of_mappings >= config.MAX_MAPPINGS ):
            res[ primer ][ 'MAPPING_SUMMARY' ] = '%d mappings' % nr_of_mappings
            res[ primer ][ 'MAPPING_SUMMARY' ] += " on %d chromosomes" % nr_of_chromosomes





    return res




def map_sequences( infile = None, perfect_only = 0 ):

    if infile is None:
        print "Infile is reqired!"
        exit()

    tmpfile = core.tmpfilename(path="/tmp/", postfix=".blastn")
    core.add_new_tmp_file( tmpfile )

    cmd = "{blastn} -db {blast_db} -outfmt '7 qacc sacc qlen pident length mismatch gapopen qstart qend sstart send evalue  bitscore'  -word_size 7 -max_hsps 10 -max_target_seqs 4 -query {query} -out {outfile}".format(blastn   = config.BLASTN, 
                                                                                                blast_db = config.BLAST_DB, 
                                                                                                query    = infile,
                                                                                                outfile  = tmpfile)
    print cmd
    subprocess.call(cmd, shell=True)

    QACC = 0
    SACC = 1 
    QLEN =2 
    PIDENT = 3 
    LENGTH = 4 
    MISMATCH = 5 
    GAPOPEN = 6 
    QSTART = 7 
    QEND = 8  
    SSTART = 9 
    SEND  = 10 
    EVALUE = 11
    BITSCORE = 12


    match = {}
    smalt_report = []
    query_region = []

    res = {}
    with open(tmpfile, 'rU') as blastn_output:

        for line in blastn_output:
        
            if (line.startswith("#")):
                continue

            print line
            
            line = line.rstrip("\n")
            fields = line.split("\t")
        
            match_id        = fields[ QACC ] 
        
            match_chrom     = fields[ SACC ]
            match_pos       = int(fields[ SSTART ])
            qstart          = int(fields[ QSTART ])

        
            seq_length      = int( fields[ QLEN ])
            match_length    = int( fields[ LENGTH ])
            match_perc      = float( fields[ PIDENT ])

            mismatch_bases  = int( fields[ MISMATCH ])

            if ( seq_length != match_length ):
                continue

            if ( perfect_only and match_perc != 100.00):
                continue
            
            match_gaps      = int( fields[ GAPOPEN ])

            if ( match_gaps ):
                continue 

            match_strand    = "plus"
            if ( int( fields[ SSTART ] ) > int( fields[ SEND ])):
                match_strand    = "minus"
                match_pos       = int(fields[ SEND ])


            if ( qstart ):
                match_pos - qstart


            if match_id not in res:
                res[ match_id ] = {}
                res[ match_id ][ 'CHR' ]   = []
                res[ match_id ][ 'POS' ]   = []
                res[ match_id ][ 'STRAND'] = []
                res[ match_id ][ 'QUAL'] = []

            res[ match_id ][ 'CHR' ].append( match_chrom )
            res[ match_id ][ 'POS' ].append( match_pos )
            res[ match_id ][ 'STRAND'].append( match_strand )
            res[ match_id ][ 'QUAL'].append( "{}/{}".format( seq_length, seq_length - match_length + mismatch_bases ) )




    return res


    
    
