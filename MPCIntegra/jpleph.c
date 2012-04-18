/*
*****    jpl planetary and lunar ephemerides    *****     C ver.1.2 *
******************************************************************************
*                                                                            *
*  This program was written in standard fortran-77 and it was manually       *
*  translated to C language by Piotr A. Dybczynski (dybol@phys.amu.edu.pl),  *
*  subsequently revised heavily by Bill J Gray (pluto@gwi.net),  just short  *
*  of a total re-write.                                                      *
*                                                                            *
******************************************************************************
*                 Last modified: July 23, 1997 by PAD                        *
******************************************************************************
21 Apr 2010:  Revised by Bill J. Gray.  The code now determines the kernel
size,  then allocates memory accordingly.  This should 'future-proof' us in
case JPL (or someone else) creates kernels that are larger than the previously
arbitrary MAX_KERNEL_SIZE parameter.  'swap_long' and 'swap_double' have
been replaced with 'swap_32_bit_val' and 'swap_64_bit_val'.  It also works
on 64-bit compiles now.

16 Mar 2001:  Revised by Bill J. Gray.  You can now use binary
ephemerides with either byte order ('big-endian' or 'small-endian');
the code checks to see if the data is in the "wrong" order for the
current platform,  and swaps bytes on-the-fly if needed.  (Yes,  this
                                                           can result in a slowdown... sometimes as much as 1%.  The function is
                                                           so mathematically intensive that the byte-swapping is the least of our
                                                           troubles.)  You can also use DE-200, 403, 404, 405,  or 406 without
recompiling (the constan( ) function now determines which ephemeris is
             in use and its byte order).

Also,  I did some minor optimization of the interp( ) (Chebyshev
                                                       interpolation) function,  resulting in a bit of a speedup.

The code has been modified to be a separately linkable component,  with
details of the implementation encapsulated.
*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>

/**** include variable and type definitions, specific for this C version */

#include "jpleph.h"
#include "jpl_int.h"

#define TRUE 1
#define FALSE 0

double   jpl_get_double( const void *ephem, const int value)
{
    return( *(double *)( (char *)ephem + value));
}

long   jpl_get_long( const void *ephem, const int value)
{
    return( *(int32_t *)( (char *)ephem + value));
}



/*****************************************************************************
 **           jpl_pleph( ephem,et,ntar,ncent,rrd,calc_velocity)              **
 ******************************************************************************
 **                                                                          **
 **    This subroutine reads the jpl planetary ephemeris                     **
 **    and gives the position and velocity of the point 'ntarg'              **
 **    with respect to 'ncent'.                                              **
 **                                                                          **
 **    Calling sequence parameters:                                          **
 **                                                                          **
 **      et = (double) julian ephemeris date at which interpolation          **
 **           is wanted.                                                     **
 **                                                                          **
 **    ntarg = integer number of 'target' point.                             **
 **                                                                          **
 **    ncent = integer number of center point.                               **
 **                                                                          **
 **    The numbering convention for 'ntarg' and 'ncent' is:                  **
 **                                                                          **
 **            1 = mercury           8 = neptune                             **
 **            2 = venus             9 = pluto                               **
 **            3 = earth            10 = moon                                **
 **            4 = mars             11 = sun                                 **
 **            5 = jupiter          12 = solar-system barycenter             **
 **            6 = saturn           13 = earth-moon barycenter               **
 **            7 = uranus           14 = nutations (longitude and obliq)     **
 **                                 15 = librations, if on eph. file         **
 **                                                                          **
 **            (If nutations are wanted, set ntarg = 14.                     **
 **             For librations, set ntarg = 15. set ncent= 0)                **
 **                                                                          **
 **     rrd = output 6-element, double array of position and velocity        **
 **           of point 'ntarg' relative to 'ncent'. The units are au and     **
 **           au/day. For librations the units are radians and radians       **
 **           per day. In the case of nutations the first four words of      **
 **           rrd will be set to nutations and rates, having units of        **
 **           radians and radians/day.                                       **
 **                                                                          **
 **           The option is available to have the units in km and km/sec.    **
 **           for this, set km=TRUE at the begining of the program.          **
 **                                                                          **
 **     calc_velocity = integer flag;  if nonzero,  velocities will be       **
 **           computed,  otherwise not.                                      **
 **                                                                          **
 *****************************************************************************/
int   jpl_pleph( void *ephem, const double et, const int ntarg,
                       const int ncent, double rrd[], const int calc_velocity)
{
    struct jpl_eph_data *eph = (struct jpl_eph_data *)ephem;
    double pv[13][6];/* pv is the position/velocity array
                      NUMBERED FROM ZERO: 0=Mercury,1=Venus,...
                      8=Pluto,9=Moon,10=Sun,11=SSBary,12=EMBary
                      First 10 elements (0-9) are affected by
                      jpl_state(), all are adjusted here.         */
    
    
    int rval = 0, list_val = (calc_velocity ? 2 : 1);
    int i, k, list[12];    /* list is a vector denoting, for which "body"
                            ephemeris values should be calculated by
                            jpl_state():  0=Mercury,1=Venus,2=EMBary,...,
                            8=Pluto,  9=geocentric Moon, 10=nutations in
                            long. & obliq.  11= lunar librations  */
    
    for( i = 0; i < 6; ++i) rrd[i] = 0.0;
    
    if( ntarg == ncent) return( 0);
    
    for( i = 0; i < 12; i++) list[i] = 0;
    
    /*   check for nutation call    */
    
    if( ntarg == 14)
    {
        if( eph->ipt[11][1] > 0) /* there is nutation on ephemeris */
        {
            list[10] = list_val;
            if( jpl_state( ephem, et, list, pv, rrd, 0))
                rval = -1;
        }
        else          /*  no nutations on the ephemeris file  */
            rval = -2;
        return( rval);
    }
    
    /*  check for librations   */
    
    if( ntarg == 15)
    {
        if( eph->ipt[12][1] > 0) /* there are librations on ephemeris file */
        {
            list[11] = list_val;
            if( jpl_state( eph, et, list, pv, rrd, 0))
                rval = -3;
            for( i = 0; i < 6; ++i)
                rrd[i] = pv[10][i]; /* librations */
        }
        else          /*  no librations on the ephemeris file  */
            rval = -4;
        return( rval);
    }
    
    /*  force barycentric output by 'state'     */
    
    /*  set up proper entries in 'list' array for state call     */
    
    for( i = 0; i < 2; i++) /* list[] IS NUMBERED FROM ZERO ! */
    {
        k = ntarg-1;
        if( i == 1) k=ncent-1;            /* same for ntarg & ncent */
        if( k <= 9) list[k] = list_val;   /* Major planets */
        if( k == 9) list[2] = list_val;   /* for moon,  earth state is needed */
        if( k == 2) list[9] = list_val;   /* for earth,  moon state is needed */
        if( k == 12) list[2] = list_val;  /* EMBary state additionaly */
    }
    
    /*   make call to state   */
    
    if( jpl_state( eph, et, list, pv, rrd, 1))
        rval = -5;
    /* Solar System barycentric Sun state goes to pv[10][] */
    if( ntarg == 11 || ncent == 11)
        for( i = 0; i < 6; i++)
            pv[10][i] = eph->pvsun[i];
    
    /* Solar System Barycenter coordinates & velocities equal to zero */
    if( ntarg == 12 || ncent == 12)
        for( i = 0; i < 6; i++)
            pv[11][i] = 0.0;
    
    /* Solar System barycentric EMBary state:  */
    if( ntarg == 13 || ncent == 13)
        for( i = 0; i < 6; i++)
            pv[12][i] = pv[2][i];
    
    /* if moon from earth or earth from moon ..... */
    if( (ntarg*ncent) == 30 && (ntarg+ncent) == 13)
        for( i = 0; i < 6; ++i) pv[2][i]=0.0;
    else
    {
        if( list[2])           /* calculate earth state from EMBary */
            for( i = 0; i < list[2] * 3; ++i)
                pv[2][i] -= pv[9][i]/(1.0+eph->emrat);
        
        if(list[9]) /* calculate Solar System barycentric moon state */
            for( i = 0; i < list[9] * 3; ++i)
                pv[9][i] += pv[2][i];
    }
    
    for( i = 0; i < list_val * 3; ++i)
        rrd[i] = pv[ntarg-1][i] - pv[ncent-1][i];
    
    return( rval);
}

/*****************************************************************************
 **                     interp(buf,t,ncf,ncm,na,ifl,pv)                      **
 ******************************************************************************
 **                                                                          **
 **    this subroutine differentiates and interpolates a                     **
 **    set of chebyshev coefficients to give position and velocity           **
 **                                                                          **
 **    calling sequence parameters:                                          **
 **                                                                          **
 **      input:                                                              **
 **                                                                          **
 **      iinfo   stores certain chunks of interpolation info,  in hopes      **
 **              that if you call again with similar parameters,  the        **
 **              function won't have to re-compute all coefficients/data.    **
 **                                                                          **
 **       coef   1st location of array of d.p. chebyshev coefficients        **
 **              of position                                                 **
 **                                                                          **
 **          t   t[0] is double fractional time in interval covered by       **
 **              coefficients at which interpolation is wanted               **
 **              (0 <= t[0] <= 1).  t[1] is dp length of whole               **
 **              interval in input time units.                               **
 **                                                                          **
 **        ncf   # of coefficients per component                             **
 **                                                                          **
 **        ncm   # of components per set of coefficients                     **
 **                                                                          **
 **         na   # of sets of coefficients in full array                     **
 **              (i.e., # of sub-intervals in full interval)                 **
 **                                                                          **
 **         ifl  integer flag: =1 for positions only                         **
 **                            =2 for pos and vel                            **
 **                                                                          **
 **                                                                          **
 **      output:                                                             **
 **                                                                          **
 **    posvel   interpolated quantities requested.  dimension                **
 **              expected is posvel[ncm*ifl], double precision.              **
 **                                                                          **
 *****************************************************************************/
static void interp( struct interpolation_info *iinfo,
                   const double coef[], const double t[2], const int ncf,
                   const int ncm, const int na, const int ifl, double posvel[])
{
    double dna, dt1, temp, tc, vfac, temp1;
    int l, i, j;
    
    /*  entry point. get correct sub-interval number for this set
     of coefficients and then get normalized chebyshev time
     within that subinterval.                                             */
    
    dna = (double)na;
    modf( t[0], &dt1);
    temp = dna * t[0];
    l = (int)(temp - dt1);
    
    /*  tc is the normalized chebyshev time (-1 <= tc <= 1)    */
    
    tc = 2.0 * (modf( temp, &temp1) + dt1) - 1.0;
    
    /*  check to see whether chebyshev time has changed,
     and compute new polynomial values if it has.
     (the element iinfo->pc[1] is the value of t1[tc] and hence
     contains the value of tc on the previous call.)     */
    
    if(tc != iinfo->pc[1])
    {
        iinfo->np = 2;
        iinfo->nv = 3;
        iinfo->pc[1] = tc;
        iinfo->twot = tc+tc;
    }
    
    /*  be sure that at least 'ncf' polynomials have been evaluated
     and are stored in the array 'iinfo->pc'.    */
    
    if(iinfo->np < ncf)
    {
        double *pc_ptr = iinfo->pc + iinfo->np;
        
        for(i=ncf - iinfo->np; i; i--, pc_ptr++)
            *pc_ptr = iinfo->twot * pc_ptr[-1] - pc_ptr[-2];
        iinfo->np=ncf;
    }
    
    /*  interpolate to get position for each component  */
    
    for( i = 0; i < ncm; ++i)        /* ncm is a number of coordinates */
    {
        const double *coeff_ptr = coef + ncf * (i + l * ncm + 1);
        const double *pc_ptr = iinfo->pc + ncf;
        
        posvel[i]=0.0;
        for( j = ncf; j; j--)
            posvel[i] += (*--pc_ptr) * (*--coeff_ptr);
    }
    
    if(ifl <= 1) return;
    
    /*  if velocity interpolation is wanted, be sure enough
     derivative polynomials have been generated and stored.    */
    
    vfac=(dna+dna)/t[1];
    iinfo->vc[2] = iinfo->twot + iinfo->twot;
    if( iinfo->nv < ncf)
    {
        double *vc_ptr = iinfo->vc + iinfo->nv;
        const double *pc_ptr = iinfo->pc + iinfo->nv - 1;
        
        for( i = ncf - iinfo->nv; i; i--, vc_ptr++, pc_ptr++)
            *vc_ptr = iinfo->twot * vc_ptr[-1] + *pc_ptr + *pc_ptr - vc_ptr[-2];
        iinfo->nv = ncf;
    }
    
    /*  interpolate to get velocity for each component    */
    
    for( i = 0; i < ncm; ++i)
    {
        double tval = 0.;
        const double *coeff_ptr = coef + ncf * (i + l * ncm + 1);
        const double *vc_ptr = iinfo->vc + ncf;
        
        for( j = ncf; j; j--)
            tval += (*--vc_ptr) * (*--coeff_ptr);
        posvel[ i + ncm] = tval * vfac;
    }
    return;
}

/* swap_32_bit_val() and swap_64_bit_val() are used when reading a binary
 ephemeris that was created on a machine with 'opposite' byte order to
 the currently-used machine (signalled by the 'swap_bytes' flag in the
 jpl_eph_data structure).  In such cases,  every double and integer
 value read from the ephemeris must be byte-swapped by these two functions. */

#define SWAP_MACRO( A, B, TEMP)   { TEMP = A;  A = B;  B = TEMP; }

static void swap_32_bit_val( void *ptr)
{
    char *tptr = (char *)ptr, tchar;
    
    SWAP_MACRO( tptr[0], tptr[3], tchar);
    SWAP_MACRO( tptr[1], tptr[2], tchar);
}

static void swap_64_bit_val( void *ptr, long count)
{
    char *tptr = (char *)ptr, tchar;
    
    while( count--)
    {
        SWAP_MACRO( tptr[0], tptr[7], tchar);
        SWAP_MACRO( tptr[1], tptr[6], tchar);
        SWAP_MACRO( tptr[2], tptr[5], tchar);
        SWAP_MACRO( tptr[3], tptr[4], tchar);
        tptr += 8;
    }
}

/*****************************************************************************
 **                        jpl_state(ephem,et2,list,pv,nut,bary)             **
 ******************************************************************************
 ** This subroutine reads and interpolates the jpl planetary ephemeris file  **
 **                                                                          **
 **    Calling sequence parameters:                                          **
 **                                                                          **
 **    Input:                                                                **
 **                                                                          **
 **        et2[] double, 2-element JED epoch at which interpolation          **
 **              is wanted.  Any combination of et2[0]+et2[1] which falls    **
 **              within the time span on the file is a permissible epoch.    **
 **                                                                          **
 **               a. for ease in programming, the user may put the           **
 **                  entire epoch in et2[0] and set et2[1]=0.0               **
 **                                                                          **
 **               b. for maximum interpolation accuracy, set et2[0] =        **
 **                  the most recent midnight at or before interpolation     **
 **                  epoch and set et2[1] = fractional part of a day         **
 **                  elapsed between et2[0] and epoch.                       **
 **                                                                          **
 **               c. as an alternative, it may prove convenient to set       **
 **                  et2[0] = some fixed epoch, such as start of integration,**
 **                  and et2[1] = elapsed interval between then and epoch.   **
 **                                                                          **
 **       list   12-element integer array specifying what interpolation      **
 **              is wanted for each of the "bodies" on the file.             **
 **                                                                          **
 **                        list[i]=0, no interpolation for body i            **
 **                               =1, position only                          **
 **                               =2, position and velocity                  **
 **                                                                          **
 **              the designation of the astronomical bodies by i is:         **
 **                                                                          **
 **                        i = 0: mercury                                    **
 **                          = 1: venus                                      **
 **                          = 2: earth-moon barycenter                      **
 **                          = 3: mars                                       **
 **                          = 4: jupiter                                    **
 **                          = 5: saturn                                     **
 **                          = 6: uranus                                     **
 **                          = 7: neptune                                    **
 **                          = 8: pluto                                      **
 **                          = 9: geocentric moon                            **
 **                          =10: nutations in longitude and obliquity       **
 **                          =11: lunar librations (if on file)              **
 **                                                                          **
 **    output:                                                               **
 **                                                                          **
 **    pv[][6]   double array that will contain requested interpolated       **
 **              quantities.  The body specified by list[i] will have its    **
 **              state in the array starting at pv[i][0]  (on any given      **
 **              call, only those words in 'pv' which are affected by the    **
 **              first 10 'list' entries (and by list(11) if librations are  **
 **              on the file) are set.  The rest of the 'pv' array           **
 **              is untouched.)  The order of components in pv[][] is:       **
 **              pv[][0]=x,....pv[][5]=dz.                                   **
 **                                                                          **
 **              All output vectors are referenced to the earth mean         **
 **              equator and equinox of epoch. The moon state is always      **
 **              geocentric; the other nine states are either heliocentric   **
 **              or solar-system barycentric, depending on the setting of    **
 **              global variables (see below).                               **
 **                                                                          **
 **              Lunar librations, if on file, are put into pv[10][k] if     **
 **              list[11] is 1 or 2.                                         **
 **                                                                          **
 **        nut   dp 4-word array that will contain nutations and rates,      **
 **              depending on the setting of list[10].  the order of         **
 **              quantities in nut is:                                       **
 **                                                                          **
 **                       d psi  (nutation in longitude)                     **
 **                       d epsilon (nutation in obliquity)                  **
 **                       d psi dot                                          **
 **                       d epsilon dot                                      **
 **                                                                          **
 *****************************************************************************/
int   jpl_state( void *ephem, const double et, const int list[12],
                       double pv[][6], double nut[4], const int bary)
{
    struct jpl_eph_data *eph = (struct jpl_eph_data *)ephem;
    int i,j, n_intervals;
    long int nr;
    double prev_midnight, time_of_day;
    double *buf = eph->cache;
    double s,t[2],aufac;
    struct interpolation_info *iinfo = (struct interpolation_info *)eph->iinfo;
    
    
    /*  ********** main entry point **********  */
    
    s=et - 0.5;
    prev_midnight = floor( s);
    time_of_day = s - prev_midnight;
    prev_midnight += 0.5;
    
    /* here prev_midnight contains last midnight before epoch desired (in JED: *.5)
     and time_of_day contains the remaining, fractional part of the epoch    */
    
    /*   error return for epoch out of range  */
    
    if( et < eph->ephem_start || et > eph->ephem_end)
        return( -1);
    
    /*   calculate record # and relative time in interval   */
    
    nr = (long)((prev_midnight-eph->ephem_start)/eph->ephem_step)+2;
    /* add 2 to adjus for the first two records containing header data */
    if( prev_midnight == eph->ephem_end)
        nr--;
    t[0]=( prev_midnight-( (1.0*nr-2.0)*eph->ephem_step+eph->ephem_start) +
          time_of_day )/eph->ephem_step;
    
    /*   read correct record if not in core (static vector buf[])   */
    
    if( nr != eph->curr_cache_loc)
    {
        eph->curr_cache_loc = nr;
        fseek( eph->ifile, nr * eph->recsize, SEEK_SET);
        if( fread( buf, sizeof( double), (size_t)eph->ncoeff, eph->ifile) != (size_t)eph->ncoeff)
            return( -2);
        if( eph->swap_bytes)
            swap_64_bit_val( buf, eph->ncoeff);
    }
    t[1] = eph->ephem_step;
    aufac = 1.0 / eph->au;
    
    for( n_intervals = 1; n_intervals <= 8; n_intervals *= 2)
        for( i = 0; i < 11; i++)
            if( n_intervals == eph->ipt[i][2] && (list[i] || i == 10))
            {
                int flag = ((i == 10) ? 2 : list[i]);
                double *dest = ((i == 10) ? eph->pvsun : pv[i]);
                
                interp( iinfo, &buf[eph->ipt[i][0]-1], t, (int)eph->ipt[i][1], 3,
                       n_intervals, flag, dest);
                /* gotta convert units */
                for( j = 0; j < flag * 3; j++)
                    dest[j] *= aufac;
            }
    
    if( !bary)                             /* gotta correct everybody for */
        for( i = 0; i < 9; i++)            /* the solar system barycenter */
            for( j = 0; j < list[i] * 3; j++)
                pv[i][j] -= eph->pvsun[j];
    
    /*  do nutations if requested (and if on file)    */
    
    if(list[10] > 0 && eph->ipt[11][1] > 0)
        interp( iinfo, &buf[eph->ipt[11][0]-1], t, (int)eph->ipt[11][1], 2,
               (int)eph->ipt[11][2], list[10], nut);
    
    /*  get librations if requested (and if on file)    */
    
    if(list[11] > 0 && eph->ipt[12][1] > 0)
    {
        double pefau[6];
        
        interp( iinfo, &buf[eph->ipt[12][0]-1], t, (int)eph->ipt[12][1], 3,
               (int)eph->ipt[12][2], list[11], pefau);
        for( j = 0; j < 6; ++j) pv[10][j]=pefau[j];
    }
    return( 0);
}

/****************************************************************************
 **    jpl_init_ephemeris( ephemeris_filename, nam, val, n_constants)       **
 *****************************************************************************
 **                                                                         **
 **    this function does the initial prep work for use of binary JPL       **
 **    ephemerides.                                                         **
 **      const char *ephemeris_filename = full path/filename of the binary  **
 **          ephemeris (on the Willmann-Bell CDs,  this is UNIX.200, 405,   **
 **          or 406)
 **      char nam[][6] = array of constant names (max 6 characters each)    **
 **          You can pass nam=NULL if you don't care about the names        **
 **      double *val = array of values of constants                         **
 **          You can pass val=NULL if you don't care about the constants    **
 **      Return value is a pointer to the jpl_eph_data structure            **
 **      NULL is returned if the file isn't opened or memory isn't alloced  **
 ****************************************************************************/

int jpl_init_err_code = 0;

void *   jpl_init_ephemeris( const char *ephemeris_filename, char nam[][6], double *val , int * n_constants )
{
    int i, j;
    long de_version;
    char title[84];
     
    FILE *ifile = fopen( ephemeris_filename, "rb");
    struct jpl_eph_data *rval;
    struct jpl_eph_data temp_data;
    struct interpolation_info *iinfo;
    
    jpl_init_err_code = 0;
    if( !ifile)
    {
        jpl_init_err_code = -1;
        return( NULL);
    }
    
    temp_data.ifile = ifile;
    if( fread( title, 84, 1, ifile) != 1)
        jpl_init_err_code = -3;
    
    fseek( ifile, 2652L, SEEK_SET);      /* skip title & constant name data */
    
   if( fread( &temp_data, JPL_HEADER_SIZE, 1, ifile) != 1)
        jpl_init_err_code = -4;
    
    if( jpl_init_err_code)
    {
        fclose( ifile);
        return( NULL);
    }
    
    de_version = atoi( title + 26);
    
    /* A small piece of trickery:  in the binary file,  data is stored */
    /* for ipt[0...11],  then the ephemeris version,  then the         */
    /* remaining ipt[12] data.  A little switching is required to get  */
    /* the correct order. */
    temp_data.ipt[12][0] = temp_data.ipt[12][1];
    temp_data.ipt[12][1] = temp_data.ipt[12][2];
    temp_data.ipt[12][2] = temp_data.ephemeris_version;
    temp_data.ephemeris_version = de_version;
    
    temp_data.swap_bytes = ( temp_data.ncon < 0 || temp_data.ncon > 65536L);
    if( temp_data.swap_bytes)     /* byte order is wrong for current platform */
    {
        swap_64_bit_val( &temp_data.ephem_start, 1);
        swap_64_bit_val( &temp_data.ephem_end, 1);
        swap_64_bit_val( &temp_data.ephem_step, 1);
        swap_32_bit_val( &temp_data.ncon);
        swap_64_bit_val( &temp_data.au, 1);
        swap_64_bit_val( &temp_data.emrat, 1);
        
        for( j = 0; j < 3; j++)
            for( i = 0; i < 13; i++)
                swap_32_bit_val( &temp_data.ipt[i][j]);
    }
    
    /* Once upon a time,  the kernel size was determined from the */
    /* DE version.  This was not a terrible idea,  except that it */
    /* meant that when the code faced a new version,  it broke.   */
    /* Now we use some logic to compute the kernel size.          */
    temp_data.kernel_size = 4;
    for( i = 0; i < 13; i++)
        temp_data.kernel_size +=
        temp_data.ipt[i][1] * temp_data.ipt[i][2] * ((i == 11) ? 4 : 6);
    temp_data.recsize = temp_data.kernel_size * 4L;
    temp_data.ncoeff = temp_data.kernel_size / 2L;
    
    /* Rather than do three separate allocations,  everything   */
    /* we need is allocated in _one_ chunk,  then parceled out. */
    /* This looks a little weird,  but it does simplify error   */
    /* handling and cleanup.                                    */
    rval = (struct jpl_eph_data *)calloc( sizeof( struct jpl_eph_data)
                                         + sizeof( struct interpolation_info)
                                         + temp_data.recsize, 1);
    if( !rval)
    {
        jpl_init_err_code = -2;
        fclose( ifile);
        return( NULL);
    }
    memcpy( rval, &temp_data, sizeof( struct jpl_eph_data));
    
    /* The 'iinfo' struct is right after the 'jpl_eph_data' struct: */
    iinfo = (struct interpolation_info *)(rval + 1);
    rval->iinfo = (void *)iinfo;
    iinfo->np = 2;
    iinfo->nv = 3;
    iinfo->pc[0] = 1.0;
    iinfo->pc[1] = 0.0;
    iinfo->vc[1] = 1.0;
    rval->curr_cache_loc = -1L;
    
    /* The 'cache' data is right after the 'iinfo' struct: */
    rval->cache = (double *)( iinfo + 1);
    
    if( val)
    {
        fseek( ifile, rval->recsize, SEEK_SET);
        if( fread( val, sizeof( double), (size_t)rval->ncon, ifile)
           != (size_t)rval->ncon)
            jpl_init_err_code = -5;
        if( rval->swap_bytes)     /* gotta swap the constants,  too */
            swap_64_bit_val( val, rval->ncon);
    }
    
    if( nam)
    {
        fseek( ifile, 84L * 3L, SEEK_SET);   /* just after the 3 'title' lines */
        for( i = 0; i < (int)rval->ncon; i++)
            if( fread( nam[i], 6, 1, ifile) != 1)
                jpl_init_err_code = -6;
    }
    *n_constants = rval->ncon;
    
    return( rval);
}

/****************************************************************************
 **    jpl_close_ephemeris( ephem)                                          **
 *****************************************************************************
 **                                                                         **
 **    this function closes files and frees up memory allocated by the      **
 **    jpl_init_ephemeris( ) function.                                      **
 ****************************************************************************/
void   jpl_close_ephemeris( void *ephem)
{
    struct jpl_eph_data *eph = (struct jpl_eph_data *)ephem;
    
    fclose( eph->ifile);
    free( ephem);
}
    
/*************************** THE END ***************************************/

