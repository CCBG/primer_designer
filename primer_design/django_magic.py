"""
  magic functions that makes it possible to transform django data into primer_desin format 

"""

import os
import pprint as pp

import config


def primer_django_to_dicts( primers ):
    '''
    transforms django primer data to primer_design data format

    input: 
        django model data

    return:
        dict of dicts of data

    '''

    primer_dict = {}

    for primer in primers:
        # Just pull data for the untagged primers, otherwise we will get everything twice.
#        if ( primer.tagged == 1):
#            continue

        if not hasattr(primer, 'show_name'):
            continue

        primer_name = primer.show_name
        primer_dict[ primer_name ] = {}

        primer_dict[ primer_name ][ 'PRIMER_ID'  ] = primer.primer_id
        primer_dict[ primer_name ][ 'TARGET_CHR' ] = primer.chrom
        primer_dict[ primer_name ][ 'TARGET_CHR' ] = primer.chrom
        primer_dict[ primer_name ][ 'TARGET_POS' ] = primer.position
        primer_dict[ primer_name ][ 'GC_PERCENT' ] = primer.gc_content
        primer_dict[ primer_name ][ 'SEQ'        ] = primer.sequence
        primer_dict[ primer_name ][ 'TM'         ] = primer.tm

        primer_dict[ primer_name ][ 'CHR' ] = []
        primer_dict[ primer_name ][ 'POS' ] = []
        primer_dict[ primer_name ][ 'STRAND' ] = []


        for mapping in primer.primer_mappings:
            primer_dict[ primer_name ][ 'CHR' ].append( mapping.chrom )
            primer_dict[ primer_name ][ 'POS' ].append( mapping.position )
            strand = "plus"
            if ( mapping.strand == config.MINUS_STRAND ):
                strand = "minus"

            primer_dict[ primer_name ][ 'STRAND' ].append( strand )
        
    
    return primer_dict
