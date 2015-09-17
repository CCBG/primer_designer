#!/usr/bin/python
# 
# 
# 
# 
# Kim Brugger (15 Sep 2015), contact: kim@brugger.dk

import sys
import pprint
pp = pprint.PrettyPrinter(indent=4)

import shlex
import subprocess
import re


# External tools we use 
SAMTOOLS   = '/software/bin/samtools ';
REFERENCE  = '/refs/human_1kg/human_g1k_v37.fasta '
TABIX      = '/software/bin/tabix '
SMALT      = '/software/bin/smalt-0.7.6 '
PRIMER3    = '/software/bin/primer3_core '

# default parameters, can be overridden by the user (perhaps)
FLANK              = 250
NR_PRIMERS         = 4
ALLOWED_MISMATCHES = 4

VERBOSE    =  3

def verbose_print( msg, level ):
    if ( level <= VERBOSE ):
        print msg

def fetch_region( chr, start, end ):

    verbose_print( "fetch_region", 3)

    cmd = "%s faidx %s  %s:%d-%d " % ( SAMTOOLS, REFERENCE, chr, start, end )
    args = shlex.split( cmd )
    p = subprocess.Popen(args, stdout=subprocess.PIPE)

    output = p.communicate()
    sequence = ""
    for line in ( output[0].split("\n")):
        if ( re.match('>', line)):
            continue 

        sequence += line
    return sequence 


def fetch_known_SNPs( tabix_file, chr, start, end):
    verbose_print( "fetch_known_SNPs", 3)

    cmd = "%s %s  %s:%d-%d " % ( TABIX, tabix_file, chr, start, end )

    args = shlex.split( cmd )
    p = subprocess.Popen(args, stdout=subprocess.PIPE)

    output = p.communicate()
    vars = []

#    pp.pprint( output );

    for line in ( output[0].split("\n")):
        vars.append( line.split("\t"));

    return vars 



def write_primer3_file(seq_id, seq, primer3_file = ""):
    verbose_print( "write_primer3_file", 4)

    template = '''SEQUENCE_ID={seq_id}
SEQUENCE_TEMPLATE={seq}
SEQUENCE_TARGET=251,1
PRIMER_FIRST_BASE_INDEX=1
PRIMER_TASK=generic
PRIMER_MIN_THREE_PRIME_DISTANCE=3
PRIMER_EXPLAIN_FLAG=1
PRIMER_MAX_LIBRARY_MISPRIMING=12.00
PRIMER_PAIR_MAX_LIBRARY_MISPRIMING=20.00
PRIMER_PRODUCT_SIZE_RANGE=150-500
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
PRIMER_MAX_SIZE=25
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
        "seq" : seq
    }



    with open( primer3_file, 'w+') as outfile:
        outfile.write(template.format(**details))

    outfile.close()


# This captures the $_ output of the primer3 engine and returns it as a dictionary.
def run_primer3( seq_id, seq, primer3_file = ""):
    verbose_print( "run_primer3", 3)

    if ( primer3_file == "" ):
        primer3_file = seq_id + ".primer3"

    primer3_file = re.sub(":", "_", primer3_file)


    write_primer3_file(seq_id, seq, primer3_file)


    cmd = "/software/bin/primer3_core -strict_tags < " + primer3_file


#    print cmd
    args = shlex.split( cmd )
    process = subprocess.Popen(cmd,
                               stdout=subprocess.PIPE,
                               shell=True)
  
    output = process.communicate()
    output_dict = dict()
  
    for line in output[0].split("\n"):
    #print line

        if (line == '='):
            break

        key, value = line.split("=")

        output_dict[ key ] = value


  #pp.pprint( output_dict )

    return output_dict



def check_primers( region_id, target_region, primer3_dict, smalt_file = None ):
    verbose_print( "check_primers", 3)

    if ( smalt_file is not None ):

        primers_file = region_id + "_p3seq.fasta"

        with open( primers_file, 'w+') as primerfasta:

            primerfasta.write(">FULLSEQ\n" + target_region + "\n")

            for i in range(0, 5):
                pl_id = "PRIMER_LEFT_%d_SEQUENCE" % i
                pr_id = "PRIMER_RIGHT_%d_SEQUENCE" % i

                if ( pr_id not in primer3_dict or pl_id not in primer3_dict):
                    continue

                primerfasta.write(">%s\n%s\n" % (pl_id, primer3_dict[ pl_id ]))
                primerfasta.write(">%s\n%s\n" % (pr_id, primer3_dict[ pr_id ]))

        primerfasta.close()

        smalt_results = region_id + ".smalt"
        cmd = SMALT + " map  -d -1 /refs/human_1kg/human_g1k_v37 " + primers_file + "> " + smalt_results + " 2> /dev/null"
        cmd = SMALT + " map  -d -1 /refs/human_1kg/human_g1k_v37 " + primers_file + "> " + smalt_results 

#    print cmd
        subprocess.call(cmd, shell=True)


    else:
        smalt_results = smalt_file

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

        if (len( res[ primer ][ 'CHR' ]) > 1 ):
            print primer + " matches more than one place in the genome "
            del (res[ primer ])


#    pp.pprint( smalt_report)
#    pp.pprint( res )


    return res

def revDNA( string ):
    verbose_print( "revDNA", 5)
    rev_bases = { 'A': 'T', 
                  'a': 'T', 
                  'C': 'G',
                  'c': 'G',
                  'G': 'C',
                  'g': 'C',
                  'T': 'A',
                  'T': 'A',
                  '-': '-'}

#    print "STRING FOR REVDNA:: " + string
    rev_str = len(string)*[None]
    for i in range(0, len(string)):
        rev_str[ len(string) - i - 1] =  rev_bases[ string[ i ]]

    return "".join( rev_str )



def align_primers_to_seq( seq, all_primers):
    verbose_print( "align_primers_to_seq", 4)

    primers = []
    rev_primers = []
    for primer in all_primers:
#        print primer
        primers.append(  primer )
        rev_primers.append( revDNA( primer ))

#    pp.pprint( primers )

    mappings = []

    for i in range(0, len(seq)):
        for primer in primers:
            primer_len = len(primer)
            if ("".join(seq[i:i+primer_len]) == primer ):
#                print primer + " Matches at pos " + str( i )

                mappings.append( [primer, i, 0] )

        for primer in rev_primers:
            primer_len = len(primer)
            if ("".join(seq[i:i+primer_len]) == primer ):
                print primer + " Matches at pos " + str( i ) + " minus! "

                mappings.append( [primer, i, 1] )


    return mappings



def get_and_parse_arguments():
    verbose_print( "get_and_parse_arguments", 4)

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('-c', '--chr')
    parser.add_argument('-p', '--pos')
    parser.add_argument('-o', '--output')
    parser.add_argument('-f', '--flank')

    args = parser.parse_args()

    return args



def markup_sequence( chr, pos, flank, sequence):
    verbose_print( "tag_sequence", 3)

    sequence = list(sequence)
    tags = [" "] * len( target_sequence )
    # Our target base
    tags[ FLANK ] = '*'

    # Tag the target sequence, regardless of the nature of the variant only one (1) base is tagged as the target.
    sequence[ flank ] = ' [' + target_sequence [ flank ] + '] '

    dbSNPs = fetch_known_SNPs( '/software/dev/projects/local_dbsnp/test.tab.gz', chr, pos - flank, pos + flank )

    masked_positions = []

    for dbSNP in (dbSNPs):
        if ( len( dbSNP ) < 6):
            continue
    
        (snp_chr, snp_pos, snp_id, snp_ref, snp_alt, common, vld, caf) = dbSNP

        snp_pos = int( snp_pos )
        if ( common == '1'):
            #        pp.pprint( dbSNP )
            mask_pos = snp_pos - (pos - flank) 

        # Sometimes the SNP we are looking at is also in dbsnp, If this is the case we dont tag it twice
            if ( mask_pos == flank ):
                continue
        
        # In the odd case we are looking at common deletion, mask the whole region. Normally this will just be one base
            for i in range(0, len(snp_ref)):

                # If already masked skip masking it again.
                if ( re.search('<', sequence[ mask_pos + i  ]) or 
                     re.search('\[',sequence[ mask_pos + i  ])):
                    continue

                sequence[ mask_pos + i  ] = ' <' + sequence[ mask_pos + i ] + '> '
                tags[  mask_pos + i  ] = 'X'

            else:
#        target_sequence[ snp_pos - pos + 1  ] = ' {' + target_sequence[ snp_pos - pos + 1  ] + '} '
                pass





    sequence =  "".join( sequence )
    sequence = re.sub(' ', '', sequence)



    tagged_string = "".join( tags )


    return ( sequence, tagged_string  )



def check_if_primer_clash( mapped_list, start, end):
    
    verbose_print( "check_if_primers_clash", 3)
    for i in range (start, end):
        if (mapped_list[ i ]  != " "):
            return True 

    return False


def extract_passed_primer_seqs( primer3_results, passed_primers):
    
    verbose_print( "extract_passed_primer_seqs", 3)

    primer_seqs = []
    for primer in passed_primers:
        if ( primer == 'FULLSEQ'):
            continue

        primer_seqs.append( primer3_results[ primer ] )


    return primer_seqs

def make_primer_mapped_strings( target_sequence, passed_primer_seqss):

    verbose_print( "make_primer_mapped_strings", 3)

    

    mappings = align_primers_to_seq(target_sequence, passed_primer_seqs )
    mapped_strings = []
    mapped_strings.append( [" "]*len( target_sequence ) )

    for mapping in mappings:
        ( primer, pos, strand ) = mapping

        mapping_index = 0

        if ( strand == 1 ):
            primer = revDNA( primer )
            primer = "<" + primer + "<"
        else:
            primer = ">" + primer +">"

        pos -= 1

        primer = list( primer )

        scan_index = 0
        while( 1 ):
            if (check_if_primer_clash ( mapped_strings[ scan_index ], pos, pos + len( primer ))):


                print "increasing scan_index by one " + str( scan_index )
                import time
                time.sleep( 1 )
                scan_index += 1

                if (len( mapped_strings ) >= scan_index ):
                    mapped_strings.append([" "]*len( target_sequence ))

            else:
                break


        mapping_index = scan_index 

        for i in range(0, len( primer)):
            mapped_strings[ mapping_index] [ pos + i ] = primer[ i ]



    for i in range(0, len(mapped_strings)):
        mapped_strings[ i ] = "".join( mapped_strings[ i ])

    return mapped_strings


def pretty_print_mappings( target_sequence, tagged_string, primer_strings, base1):
    
    for i in range(0, len(tagged_sequence), 80):

        print "%-9d  %s" % ( base1+ i, target_sequence[i: i+80])
        print "           " + tagged_string[i: i+80]

        for primer_string in ( primer_strings ):
            print "           " + primer_string[i: i+80]

        print ""



#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
#    Main loop !
#
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


args = get_and_parse_arguments()

chr = args.chr
pos = int(args.pos)
FLANK = args.flank or FLANK
FLANK = int( FLANK )

region_id = "%s:%d" % ( chr, pos)

target_sequence                   =   fetch_region( chr, pos - FLANK, pos + FLANK     )
(tagged_sequence, tagged_string)  =   markup_sequence( chr, pos, FLANK, target_sequence  )


primer3_results  = run_primer3( region_id , tagged_sequence, "%s_%d.primer3" % ( chr, pos))
passed_primers   = check_primers( region_id, target_sequence, primer3_results, '3:12393125.smalt')
passed_primer_seqs = extract_passed_primer_seqs( primer3_results, passed_primers )


mapped_primer_strings = make_primer_mapped_strings( target_sequence, passed_primer_seqs)


pretty_print_mappings( target_sequence, tagged_string, mapped_primer_strings, pos - FLANK)

#print tagged_sequence
#print tagged_string

#for mapped_primer_string in mapped_primer_strings:
#    print mapped_primer_string
#pp.pprint( tagged_positions )

#print  target_sequence 
#print  mapped_string 
