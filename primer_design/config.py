"""
   Variables pointing to 3rd party software and databases

"""


VERSION    =  '3.0-beta1'


# External tools and databases used by the program
SMALT      = '/software/bin/smalt-0.7.6 '
PRIMER3    = '/software/bin/primer3_core '
SAMTOOLS   = '/software/bin/samtools ';
REFERENCE  = '/refs/human_1kg/human_g1k_v37.fasta '

PRIMER3    = '/software/bin/primer3_core '
TABIX      = '/software/bin/tabix '
SMALT      = '/software/bin/smalt-0.7.6 '
PRIMER3    = '/software/bin/primer3_core '
BLASTN     = '/software/packages/ncbi-blast-2.5.0+/bin/blastn'
BLAST_DB   = '/refs/human_1kg/human_g1k_v37.fasta'

# data analysis parameters.



# default parameters, can be overridden by the user (perhaps)
TARGET_LEAD          = 50
FLANK                = 500 # amount of bp on either side of the target to design primers in
NR_PRIMERS           = 4   # Nr of primers to report back
ALLOWED_MISMATCHES   = 4   # Maximum nr of errors when mapping a primer back to the reference
MAX_MAPPINGS         = 5   # Less or equal number of mappings to the reference what chromosomes mapped to are named 
MAX_PRODUCT          = 1200
MAX_3_PRIME_MISMATCH = 2

