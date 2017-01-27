
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
