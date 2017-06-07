"""

  Functions that are primer3 specific.

"""

import re
import shlex
import subprocess

import core
import config


if ( __name__ == '__main__'):
    sys.stderr.write( "This is a module and not meant to be run as a stand alone program\n" )
    exit( 1 )


def markup_sequence( target_sequence, target_chrom, target_start, target_end, target_flank, target_lead=config.TARGET_LEAD):
    """

    Tag up a sequence with the target-region, common SNPs found in the region
    furthermore add a flanking sequence to the target to ensure that the primers are 
    long enough away from the target to generate sequene. Normally this is about 50bp 
    that is lost during sanger sequencing


    Kim Brugger (17 Aug 2016)

    """


    core.verbose_print( "tag_sequence", 2)

    sequence = list(target_sequence)
    tags = [" "] * len( target_sequence )
    # Our target base
    start = target_flank
    end = len(tags) - target_flank

    for x in range(0, len(tags)):
        if (x >= start and x <= end):
            tags[x] = '*'

        if x in range(start - target_lead, start) or x in range(end, end + target_lead):
            tags[x] = '-'
    #return tags
        

    core.verbose_print( "::::: %d - %d, %d" % (target_start, target_end, target_flank), 5)
    # Tag the target sequence, regardless of the nature of the variant only one (1) base is tagged as the target.
    sequence[ target_flank - target_lead ] = ' [' + target_sequence [ target_flank - target_lead ]
    sequence[(- target_flank + target_lead) -1 ] = sequence [ (- target_flank + target_lead) -1 ] + '] '


    dbSNPs = core.fetch_known_SNPs( config.DBSNP_FILE, target_chrom, 
                                    target_start - target_flank, 
                                    target_end + target_flank )

    masked_positions = []

    for dbSNP in (dbSNPs):
        #print dbSNP
        #exit()
        if ( len( dbSNP ) < 6):
            continue
        
        #unpack dbSNP entry. 
        (snp_chrom, snp_pos, snp_id, snp_ref, snp_alt, common, vld, caf) = dbSNP

        snp_pos = int( snp_pos )
        if (snp_pos >= target_start - target_lead and snp_pos <= target_end + target_lead ):
            continue

        if ( common == '1'):
            #        pp.pprint( dbSNP )
            mask_pos = snp_pos - (target_start - target_flank) 

        # Sometimes the SNP we are looking at is also in dbsnp, If this is the case we dont tag it twice
            
        # In the odd case we are looking at common deletion, mask the whole region. Normally this will just be one base
            for i in range(0, len(snp_ref)):
                if (len(sequence)<=  mask_pos + i):
                    break

                # If already masked skip masking it again.
                if ( re.search('<', sequence[ mask_pos + i  ]) or 
                     re.search('\[',sequence[ mask_pos + i  ])):
                    continue

                sequence[ mask_pos + i  ] = ' <' + sequence[ mask_pos + i ] + '> '
#                sequence[ mask_pos + i  ] = 
                tags[  mask_pos + i  ] = 'X'

            else:
#        target_sequence[ snp_pos - pos + 1  ] = ' {' + target_sequence[ snp_pos - pos + 1  ] + '} '
                pass





    sequence =  "".join( sequence )
    sequence = re.sub(' ', '', sequence)


    #print tags
    #exit()
    tagged_string = "".join( tags )


    return ( sequence, tagged_string  )



def write_file(seq_id, seq, primer3_file = "primer3.tmp"):
    """
    Generate the file required to run primer3
    
    input: sequence-id, sequence, output-file (optional)
    
    output: None
    
    Kim Brugger (17 Aug 2016)
    """

    core.verbose_print( "write_file", 2)

    upstream_end      = 1
    downstream_start =  2
    for i in range(0, len(seq)):
            if seq[ i ] == '[':
                upstream_end = i
            elif seq[ i ] == ']':
                downstream_start = i


    template = '''SEQUENCE_ID={seq_id}
SEQUENCE_TEMPLATE={seq}
SEQUENCE_TARGET={flank},{len}
PRIMER_FIRST_BASE_INDEX=1
PRIMER_TASK=generic
PRIMER_MIN_THREE_PRIME_DISTANCE=3
PRIMER_EXPLAIN_FLAG=1
PRIMER_MAX_LIBRARY_MISPRIMING=12.00
PRIMER_PAIR_MAX_LIBRARY_MISPRIMING=20.00
PRIMER_PRODUCT_SIZE_RANGE=200-800
PRIMER_NUM_RETURN=5
PRIMER_MAX_END_STABILITY=9.0
PRIMER_MAX_SELF_ANY_TH=45.00
PRIMER_MAX_SELF_END_TH=35.00
PRIMER_PAIR_MAX_COMPL_ANY_TH=45.00
PRIMER_PAIR_MAX_COMPL_END_TH=35.00
PRIMER_MAX_HAIRPIN_TH=24.00
PRIMER_MAX_TEMPLATE_MISPRIMING_TH=40.00
PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH=70.00
PRIMER_MIN_SIZE=18
PRIMER_OPT_SIZE=20
PRIMER_MAX_SIZE=22
PRIMER_MIN_TM=58.0
PRIMER_OPT_TM=60.0
PRIMER_MAX_TM=62.0
PRIMER_PAIR_MAX_DIFF_TM=5.0
PRIMER_TM_FORMULA=1
PRIMER_SALT_MONOVALENT=50.0
PRIMER_SALT_CORRECTIONS=1
PRIMER_SALT_DIVALENT=1.5
PRIMER_DNTP_CONC=0.6
PRIMER_DNA_CONC=50.0
PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=1
PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT=0
PRIMER_LOWERCASE_MASKING=0
PRIMER_MIN_GC=30.0
PRIMER_MAX_GC=70.0
PRIMER_MAX_NS_ACCEPTED=0
PRIMER_MAX_POLY_X=4
PRIMER_OUTSIDE_PENALTY=0
PRIMER_GC_CLAMP=0
PRIMER_LIBERAL_BASE=1
PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS=0
PRIMER_PICK_ANYWAY=1
PRIMER_WT_TM_LT=1.0
PRIMER_WT_TM_GT=1.0
PRIMER_WT_SIZE_LT=1.0
PRIMER_WT_SIZE_GT=1.0
PRIMER_WT_GC_PERCENT_LT=0.0
PRIMER_WT_GC_PERCENT_GT=0.0
PRIMER_WT_SELF_ANY_TH=0.0
PRIMER_WT_SELF_END_TH=0.0
PRIMER_WT_HAIRPIN_TH=0.0
PRIMER_WT_NUM_NS=0.0
PRIMER_WT_LIBRARY_MISPRIMING=0.0
PRIMER_WT_SEQ_QUAL=0.0
PRIMER_WT_END_QUAL=0.0
PRIMER_WT_POS_PENALTY=0.0
PRIMER_WT_END_STABILITY=0.0
PRIMER_WT_TEMPLATE_MISPRIMING_TH=0.0
PRIMER_PAIR_WT_PRODUCT_SIZE_LT=0.0
PRIMER_PAIR_WT_PRODUCT_SIZE_GT=0.0
PRIMER_PAIR_WT_PRODUCT_TM_LT=0.0
PRIMER_PAIR_WT_PRODUCT_TM_GT=0.0
PRIMER_PAIR_WT_DIFF_TM=0.0
PRIMER_PAIR_WT_COMPL_ANY_TH=0.0
PRIMER_PAIR_WT_COMPL_END_TH=0.0
PRIMER_PAIR_WT_LIBRARY_MISPRIMING=0.0
PRIMER_PAIR_WT_PR_PENALTY=1.0
PRIMER_PAIR_WT_IO_PENALTY=0.0
PRIMER_PAIR_WT_TEMPLATE_MISPRIMING=0.0
PRIMER_INTERNAL_WT_SIZE_LT=1.0
PRIMER_INTERNAL_WT_END_QUAL=0.0
PRIMER_INTERNAL_MAX_SELF_END=12.00
PRIMER_QUALITY_RANGE_MIN=0
PRIMER_PAIR_MAX_COMPL_END=3.00
PRIMER_PRODUCT_MAX_TM=1000000.0
PRIMER_INTERNAL_MAX_SIZE=27
PRIMER_INTERNAL_WT_SELF_ANY=0.0
PRIMER_INTERNAL_MAX_POLY_X=5
PRIMER_INTERNAL_WT_SIZE_GT=1.0
PRIMER_SEQUENCING_ACCURACY=20
PRIMER_INTERNAL_WT_TM_GT=1.0
PRIMER_INTERNAL_WT_LIBRARY_MISHYB=0.0
PRIMER_INTERNAL_MAX_GC=80.0
PRIMER_PAIR_WT_COMPL_ANY=0.0
PRIMER_PICK_INTERNAL_OLIGO=0
PRIMER_MAX_SELF_END=3.00
PRIMER_QUALITY_RANGE_MAX=100
PRIMER_INTERNAL_DNTP_CONC=0.0
PRIMER_INTERNAL_MIN_SIZE=18
PRIMER_INTERNAL_MIN_QUALITY=0
PRIMER_SEQUENCING_INTERVAL=250
PRIMER_INTERNAL_SALT_DIVALENT=1.5
PRIMER_MAX_SELF_ANY=8.00
PRIMER_INTERNAL_WT_SEQ_QUAL=0.0
PRIMER_PAIR_WT_COMPL_END=0.0
PRIMER_INTERNAL_OPT_TM=60.0
PRIMER_SEQUENCING_SPACING=500
PRIMER_INTERNAL_MAX_SELF_ANY=12.00
PRIMER_MIN_END_QUALITY=0
PRIMER_INTERNAL_MIN_TM=57.0
PRIMER_PAIR_MAX_COMPL_ANY=8.00
PRIMER_SEQUENCING_LEAD=50
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_INTERNAL_OPT_SIZE=20
PRIMER_WT_TEMPLATE_MISPRIMING=0.0
PRIMER_MAX_END_GC=5
PRIMER_MIN_QUALITY=0
PRIMER_INTERNAL_MAX_LIBRARY_MISHYB=12.00
PRIMER_INTERNAL_WT_GC_PERCENT_GT=0.0
PRIMER_INTERNAL_MAX_NS_ACCEPTED=0
PRIMER_WT_SELF_ANY=0.0
PRIMER_MAX_TEMPLATE_MISPRIMING=12.00
PRIMER_INTERNAL_WT_NUM_NS=0.0
PRIMER_INTERNAL_WT_SELF_END=0.0
PRIMER_PRODUCT_OPT_SIZE=0
PRIMER_PRODUCT_OPT_TM=0.0
PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING=24.00
PRIMER_INSIDE_PENALTY=-1.0
PRIMER_INTERNAL_MIN_GC=20.0
PRIMER_PRODUCT_MIN_TM=-1000000.0
PRIMER_INTERNAL_SALT_MONOVALENT=50.0
PRIMER_WT_SELF_END=0.0
PRIMER_INTERNAL_DNA_CONC=50.0
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_INTERNAL_MAX_TM=63.0
PRIMER_INTERNAL_WT_GC_PERCENT_LT=0.0
PRIMER_INTERNAL_WT_TM_LT=1.0
PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/software/bin/primer3_config/
='''

    details = {
        "seq_id" : seq_id,
        "seq"    : seq,
        "flank"  : upstream_end,
        "len"    : downstream_start - upstream_end + 1
    }

 
    with open( primer3_file, 'w+') as outfile:
        outfile.write(template.format(**details))

    outfile.close()


def run( seq_id, seq, primer3_file = ""):
    """
    runs primer3 and returns the output as a dict
    
    input: sequence ID, sequence, primer3_file (optional)
    output: returns a dict with the results.
    
    Kim Brugger (17 Aug 2016)
    """

    core.verbose_print( "run_primer3", 2)

    if ( primer3_file is None or primer3_file == "" ):
        primer3_file = seq_id + ".primer3"


    primer3_file = re.sub(":", "_", primer3_file)

    print "PRIMER3 file: ", primer3_file


    core.add_new_tmp_file( primer3_file )

    write_file(seq_id, seq, primer3_file)


    cmd = config.PRIMER3 + " < " + primer3_file


    core.verbose_print( cmd, 3)
    args = shlex.split( cmd )
    process = subprocess.Popen(cmd,
                               stdout=subprocess.PIPE,
                               shell=True)
  
    output = process.communicate()
    output_dict = dict()
  
    for line in output[0].split("\n"):

        if (line == '='):
            break

        key, value = line.split("=")

        output_dict[ key ] = value

    return output_dict


