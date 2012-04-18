//
//  ALPPIntegra.m
//  MPCIntegra
//
//  Created by Michele Bigi on 20/01/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#import "ALPPIntegra.h"
#import "ALMinorPlanet.h"
#import "ALJplEphem.h"
#import "jpleph.h"


#define CUBE_ROOT( X)  (exp( log( X) / 3.))

@implementation ALPPIntegra

int verbose = 0, n_steps_taken = 0, resync_freq = 50;
int asteroid_perturber_number = -1;
long time_in_perturber = 0;
double *position_cache;
unsigned long perturber_mask = PERTURBERS_MERCURY_TO_NEPTUNE;

- (id) initWithMinorPlanet:(ALJplEphem *) ephem:(ALMinorPlanet *) mplanet {
    self = [super init];
    if (self) {
        // Initialization code here.
        // Neqn = 6;
        JPLEphem = ephem;
        minor_planet = mplanet;
    }
    
    return self;
}

- (void)dealloc
{
    [super dealloc];
}

void AL_add_relativistic_accel( double *accel, const double *posnvel)
{
    int i;
    double c = AU_PER_DAY;           /* speed of light in AU per day */
    double r_squared = posnvel[0] * posnvel[0] + posnvel[1] * posnvel[1]
    + posnvel[2] * posnvel[2];
    const double v_squared = posnvel[3] * posnvel[3] + posnvel[4] * posnvel[4]
    + posnvel[5] * posnvel[5];
    const double v_dot_r   = posnvel[0] * posnvel[3] + posnvel[1] * posnvel[4]
    + posnvel[2] * posnvel[5];
    const double r = sqrt( r_squared), r_cubed_c_squared = r_squared * r * c * c;
    const double r_component = (4. * SOLAR_GM / r - v_squared) / r_cubed_c_squared;
    const double v_component = 4. * v_dot_r / r_cubed_c_squared;
    
    for( i = 0; i < 3; i++)
        accel[i] += r_component * posnvel[i] + v_component * posnvel[i + 3];
}

void AL_set_differential_acceleration( const double *posnvel, const double *delta, double *accel)
{
    double p_squared = 0., r_squared = 0.;
    double pfactor, rfactor, posnvel_2[6];
    int i;
    
    for( i = 0; i < 6; i++)
        posnvel_2[i] = posnvel[i] + delta[i];
    for( i = 0; i < 3; i++)
    {
        p_squared += posnvel[i] * posnvel[i];
        r_squared += posnvel_2[i] * posnvel_2[i];
    }
    /* someday,  I'll do it right;  for the nonce,  do it quick: */
    /* SEE: \useless\smalldif.cpp */
    pfactor = 1. / (p_squared * sqrt( p_squared));
    rfactor = 1. / (r_squared * sqrt( r_squared));
    for( i = 0; i < 3; i++)
        accel[i] = pfactor * posnvel[i] - rfactor * posnvel_2[i];
    AL_add_relativistic_accel( accel, posnvel_2);
}

void * AL_jpl_ephemeris;

-(void) AL_ComputePerturber:(int) perturber_no:(double) jd:(double *) perturber_loc {
    static double jd0 = -1;
    double posns[11][6];
    int i;
    
    if( jd0 != jd)
    {
        int list[12];

        double ratio = 1. + [JPLEphem AL_EarthMoonRatio ]; 
        for( i = 0; i < 12; i++)
            list[i] = (i < 10);
        [ JPLEphem JPLState:jd :list : posns ];
        
        for( i = 0; i < 3; ++i)
        {
            posns[2][i] -= posns[9][i] / ratio;
            posns[9][i] += posns[2][i];
        }
        jd0 = jd;
        for( i = 0; i < 10; i++)
        {
            const double sin_obliq = .397777156;
            const double cos_obliq = .917482062;
            double temp = posns[i][1] * cos_obliq + posns[i][2] * sin_obliq;
            
            posns[i][2] = posns[i][2] * cos_obliq - posns[i][1] * sin_obliq;
            posns[i][1] = temp;
        }
    }
    
    for( i = 0; i < 3; ++i)
        perturber_loc[i] = posns[perturber_no - 1][i];  /* - posns[10][i]; */

}

-(double *) AL_MakePositionCache:(double) jd0:(double) stepsize:(int) n_steps {
    
    
    double *rval = (double *)calloc( 2 + n_steps * N_PERTURBERS * 6 * 3, sizeof( double));
    double *tptr = rval+2;
    int i, j;
    
    
    if( !rval)
    {
        printf( "Ran out of memory!\n");
        exit( -1);
    }
    rval[0] = jd0;
    rval[1] = stepsize;
    
    while( n_steps--)
    {
        for( j = 0; j < 6; j++)
        {
            const double avals[6] = { 0., 2. / 9., 1./3., .75, 1., 5./6. };
            
            for( i = 0; i < N_PERTURBERS; i++)
            {
                if( (i < 10) && ((perturber_mask >> i) & 1ul)) {
                    
                    [self AL_ComputePerturber:i+1 :jd0+avals[j] * stepsize :tptr ];
                }
                else {       /* put it far,  far away where it won't do anything: */
                    tptr[0] = tptr[1] = tptr[2] = 1.e+8;
                }
                tptr += 3;
            }
        }
        jd0 += stepsize;
    }
    return( rval );
}

#define MAX_ITERATIONS 7
#define THRESH 1.e-8
#define MIN_THRESH 1.e-15

double AL_near_parabolic( const double ecc_anom, const double e)
{
    const double anom2 = (e > 1. ? ecc_anom * ecc_anom : -ecc_anom * ecc_anom);
    double term = e * anom2 * ecc_anom / 6.;
    double rval = (1. - e) * ecc_anom - term;
    int n = 4;
    
    while( fabs( term) > 1e-15)
    {
        term *= anom2 / (double)(n * (n + 1));
        rval -= term;
        n += 2;
    }
    return( rval);
}

double AL_kepler( const double ecc, double mean_anom)
{
    double curr, err, thresh, offset = 0.;
    double delta_curr = 1.;
    int is_negative = 0, n_iter = 0;
    
    if( !mean_anom)
        return( 0.);
    
    if( ecc < .3)     /* low-eccentricity formula from Meeus,  p. 195 */
    {
        curr = atan2( sin( mean_anom), cos( mean_anom) - ecc);
        /* two correction steps,  and we're done */
        for( n_iter = 2; n_iter; n_iter--)
        {
            err = curr - ecc * sin( curr) - mean_anom;
            curr -= err / (1. - ecc * cos( curr));
        }
        return( curr);
    }
    
    if( ecc < 1.)
        if( mean_anom < -PI || mean_anom > PI)
        {
            double tmod = fmod( mean_anom, PI * 2.);
            
            if( tmod > PI)             /* bring mean anom within -pi to +pi */
                tmod -= 2. * PI;
            else if( tmod < -PI)
                tmod += 2. * PI;
            offset = mean_anom - tmod;
            mean_anom = tmod;
        }
    
    if( mean_anom < 0.)
    {
        mean_anom = -mean_anom;
        is_negative = 1;
    }
    
    curr = mean_anom;
    thresh = THRESH * fabs( 1. - ecc);
    /* Due to roundoff error,  there's no way we can hope to */
    /* get below a certain minimum threshhold anyway:        */
    if( thresh < MIN_THRESH)
        thresh = MIN_THRESH;
    if( (ecc > .8 && mean_anom < PI / 3.) || ecc > 1.)    /* up to 60 degrees */
    {
        double trial = mean_anom / fabs( 1. - ecc);
        
        if( trial * trial > 6. * fabs(1. - ecc))   /* cubic term is dominant */
        {
            if( mean_anom < PI)
                trial = CUBE_ROOT( 6. * mean_anom);
            else        /* hyperbolic w/ 5th & higher-order terms predominant */
                trial = asinh( mean_anom / ecc);
        }
        curr = trial;
    }
    if( ecc > 1. && mean_anom > 4.)    /* hyperbolic, large-mean-anomaly case */
        curr = log( mean_anom);
    
    if( ecc < 1.)
        while( fabs( delta_curr) > thresh)
        {
            if( n_iter++ > MAX_ITERATIONS)
                err = AL_near_parabolic( curr, ecc) - mean_anom;
            else
                err = curr - ecc * sin( curr) - mean_anom;
            delta_curr = -err / (1. - ecc * cos( curr));
            curr += delta_curr;
        }
    else
        while( fabs( delta_curr) > thresh)
        {
            if( n_iter++ > MAX_ITERATIONS)
                err = -AL_near_parabolic( curr, ecc) - mean_anom;
            else
                err = ecc * sinh( curr) - curr - mean_anom;
            delta_curr = -err / (ecc * cosh( curr) - 1.);
            curr += delta_curr;
        }
    return( is_negative ? offset - curr : offset + curr);
}

- (void)  AL_comet_posn_part_ii:(double) t:(double *) loc:(double *) vel {
    double true_anom, r, x, y, r0;
    
    if( [minor_planet Eccentricity ] == 1.)    /* parabolic */
    {
        double g = [minor_planet W0 ] * t * .5;
        
        y = CUBE_ROOT( g + sqrt( g * g + 1.));
        true_anom = 2. * atan( y - 1. / y);
    }
    else           /* got the mean anomaly;  compute eccentric,  then true */
    {
        double ecc_anom;
        
        NSLog(@"Eccentricity: %lf\tMean Anomaly: %lf\n", [minor_planet Eccentricity],[minor_planet MeanAnomaly] );
        ecc_anom = AL_kepler( [minor_planet Eccentricity ], [minor_planet MeanAnomaly ]);
        
        if( [minor_planet Eccentricity ]> 1.)     /* hyperbolic case */
        {
            x = ([minor_planet Eccentricity ] - cosh(ecc_anom));
            y = sinh(ecc_anom);
        }
        else           /* elliptical case */
        {
            x = (cos( ecc_anom) - [minor_planet Eccentricity ]);
            y =  sin( ecc_anom);
        }
        y *= [minor_planet MinorToMajor ];
        true_anom = atan2( y, x);
    }
    
    r0 = [minor_planet Q ] * (1. + [minor_planet Eccentricity ]);
    r = r0 / (1. + [minor_planet Eccentricity ] * cos( true_anom));
    x = r * cos( true_anom);
    y = r * sin( true_anom);
    double * perih_vec = [minor_planet PerihelionVector ];
    double * sideways  = [minor_planet Sideways ];
    loc[0] = perih_vec[0] * x + sideways[0] * y;
    loc[1] = perih_vec[1] * x + sideways[1] * y;
    loc[2] = perih_vec[2] * x + sideways[2] * y;
    loc[3] = r;
    if( vel && ([minor_planet AngularMomentum ]!= 0.))
    {
        double angular_component = [minor_planet AngularMomentum ] / (r * r);
        double radial_component  = [minor_planet Eccentricity ]* sin( true_anom) *
     ( [minor_planet AngularMomentum ] / (r * r0) );
        double x1 = x * radial_component - y * angular_component;
        double y1 = y * radial_component + x * angular_component;
        int i;
        
        for( i = 0; i < 3; i++)
            vel[i] = perih_vec[i] * x1 + sideways[i] * y1;
    }
    
}

- (int)  AL_comet_posn_and_vel:(double) t:(double *) loc:(double *) vel {
    
    
    t -= [minor_planet PerihelionTime ];
    
    if( [minor_planet Eccentricity ] != 1.)    /* not parabolic */
    {
        t /= [minor_planet T0 ];
        if( [minor_planet Eccentricity ] < 1.)     /* elliptical case;  throw out extra orbits */
        {                    /* to fit mean anom between -PI and PI */
            t = fmod( t, PI * 2.);
            if( t < -PI) t += 2. * PI;
            if( t >  PI) t -= 2. * PI;
        }
        [minor_planet MeanAnomaly: t ];
    }	
	[ self AL_comet_posn_part_ii:t :loc :vel ];
    return(0);
}


- (void) AL_computeDerivatives:(double) jd:(double *) delta:(double *) derivs:(double *) posn_data {
    static long double relative_mass[14] = { 1.,
    1.660136795271931e-007,                /* mercury */
    2.447838339664545e-006,                /* venus */
    3.003489596331057e-006,                /* Earth */
    3.227151445053866e-007,                /* Mars */
    0.0009547919384243268,                 /* Jupiter */
    0.0002858859806661309,                 /* saturn */
    4.366244043351564e-005,                /* Uranus */
    5.151389020466116e-005,                /* Neptune */
    7.396449704142013e-009,                /* Pluto */
    3.003489596331057e-006 / EARTH_MOON_RATIO, /* Moon */
    4.7622e-10, 1.0775e-10, 1.3412e-10 };    /* Ceres,  Pallas, Vesta */
    double accel[3];
    double posnvel[6];
    
    int i;
    
    asteroid_perturber_number = 0;
    
    [ self AL_comet_posn_and_vel:jd :posnvel :posnvel+3 ];
    AL_set_differential_acceleration( posnvel, delta, accel);
  
    
    for( i = 0; i < N_PERTURBERS; i++)       /* include perturbers */
        if( (perturber_mask >> i) & 1ul)
        {
            double perturber_loc[3], dfactor;
            long double diff[3], diff_squared = 0.;
            long double radius_squared = 0., rfactor;
            int j;
            
            if( posn_data)
                memcpy( perturber_loc, posn_data + i * 3, 3 * sizeof( double));
            else
                if( i < 10)
                    [ self AL_ComputePerturber:i+1 :jd :perturber_loc ];
                else
                    perturber_loc[0] = perturber_loc[1] = perturber_loc[2] = 1.e+8;
            for( j = 0; j < 3; j++)
            {
                diff[j] = perturber_loc[j] - (posnvel[j] + delta[j]);
                diff_squared += diff[j] * diff[j];
                radius_squared += perturber_loc[j] * perturber_loc[j];
            }
            dfactor = relative_mass[i + 1] / (diff_squared * sqrt( diff_squared));
            rfactor = relative_mass[i + 1] / (radius_squared * sqrt( radius_squared));
            for( j = 0; j < 3; j++)
                accel[j] += diff[j] * dfactor - perturber_loc[j] * rfactor;
            
        }
    
    /* copy in Ceres,  Pallas, Vesta loc if needed: */
    if( posn_data && asteroid_perturber_number >= 0)
        memcpy( posn_data + asteroid_perturber_number * 3, posnvel,
               3 * sizeof( double));
    
    for( i = 0; i < 3; i++)
    {
        derivs[i] = delta[i + 3];
        derivs[i + 3] = SOLAR_GM * accel[i];
    }
    
    NSLog(@"AL_computeDerivatives: %lg   %lg   %lg\n", derivs[0], derivs[1], derivs[2] );
    NSLog(@"AL_computeDerivatives: %lg   %lg   %lg\n", delta[0], delta[1], delta[2] );
}

-(int) AL_TakeStep:(double) jd: (double *) ival: (double *) ovals: (double *)errs :(int) n_vals: (double) step_size {
    

    double *ivals[7];
    double *ivals_p[6];
    double *posn_data = NULL;
    int i, j, k;
    const double bvals[27] = {2. / 9.,
        1. / 12., 1. / 4.,
        69. / 128., -243. / 128., 135. / 64.,
        -17. / 12., 27. / 4., -27. / 5., 16. / 15.,
        65. / 432., -5. / 16., 13 / 16., 4 / 27., 5. / 144.,
        47. / 450., 0., 12 / 25., 32. / 225., 1. / 30., 6. / 25.,
        -1. / 150., 0., .03, -16. / 75., -.05, .24};
    const double *bptr = bvals;
    const double avals[6] = { 0., 2. / 9., 1./3., .75, 1., 5./6. };
    
    
    NSLog(@"AL_TakeStep: jd = %lf\tival = %lf\tovals = %lf\t errs = %f\t n_vals = %i\n", jd, *ival, *ovals, *errs, n_vals);
    ivals[1] = (double *)calloc( 12 * n_vals, sizeof( double ) );
    
    if( !ivals[1])
    {
        printf( "Ran out of memory! (2)\n");
        exit( -1);
    }
    for( i = 0; i < 6; i++)
    {
        ivals[i + 1] = ivals[1] + i * n_vals;
        ivals_p[i] = ivals[1] + (i + 6) * n_vals;
    }
  //  NSLog( @"Position Cache[1] = %lf\n" , position_cache[1] );
    if( fabs( step_size - position_cache[1]) < .000001)
    {
        int cache_loc = (int)floor( (jd - position_cache[0]) / step_size + .5);
        
        posn_data = position_cache + 2 + cache_loc * 6 * N_PERTURBERS * 3;
    }
    
    [ self AL_computeDerivatives:jd : ival:ivals_p[0] :posn_data ];
    
    for( j = 1; j < 7; j++)
    {
        for( i = 0; i < n_vals; i++)
        {
            double tval = 0.;
            
            for( k = 0; k < j; k++)
                tval += bptr[k] * ivals_p[k][i];
            ivals[j][i] = tval * step_size + ival[i];
        }
        bptr += j;
        if( j != 6) {
           
            [ self AL_computeDerivatives:jd+step_size*avals[j]: ivals[j]:ivals_p[j] :posn_data ?
             posn_data + j * N_PERTURBERS * 3 : NULL ];
        }
    }
    
    if( errs)
        for( i = 0; i < n_vals; i++)
        {
            double tval = 0.;
            
            for( k = 0; k < 6; k++)
                tval += bptr[k] * ivals_p[k][i];
            errs[i] = step_size * tval;
        }
    
    memcpy( ovals, ivals[6], n_vals * sizeof( double));
    free( ivals[1]);
    n_steps_taken++;
    return( 0);
}

#define N_VALUES 6
/* i.e.,  a state vector consumes six values: x, y, z, vx, vy, vz */

-(int) AL_FullRKStep:(double *) ivals:(double *) ovals: (double) t0: (double) t1: (double) max_err {
    double step = t1 - t0;
    double errs[N_VALUES], new_vals[N_VALUES];
    int n_chickens = 0;
    
    memcpy( ovals, ivals, N_VALUES * sizeof( double));
    max_err *= max_err;
    while( t0 != t1)
    {
        double err_val = 0.;
        const double chicken_factor = .9;
        int i;
        
        [self AL_TakeStep:t0 :ovals :new_vals :errs :N_VALUES :step ];
        //take_step( t0, elems, ovals, new_vals, errs, N_VALUES, step);
        for( i = 0; i < N_VALUES; i++)
            err_val += errs[i] * errs[i];
        if( err_val < max_err)   /* yeah,  it was a good step */
        {
            memcpy( ovals, new_vals, N_VALUES * sizeof( double));
            t0 += step;
        }
        else
            n_chickens++;
        step *= chicken_factor * exp( log( max_err / err_val) / 5.);
        if( t0 < t1)
            if( t0 + step > t1)
                step = t1 - t0;
        if( t1 < t0)
            if( t0 + step < t1)
                step = t1 - t0;
        /*    if( err_val >= max_err)                             */
        /*       printf( "Chickened out: new step %lf\n", step);  */
    }
    return( n_chickens);
}

double AL_dot_product( const double *a, const double *b)
{
    return( a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}

double AL_vector3_length( const double *vect)
{
    return( sqrt( vect[0] * vect[0] + vect[1] * vect[1] + vect[2] * vect[2]));
}

void AL_vector_cross_product( double *xprod, const double *a, const double *b)
{
    xprod[0] = a[1] * b[2] - a[2] * b[1];
    xprod[1] = a[2] * b[0] - a[0] * b[2];
    xprod[2] = a[0] * b[1] - a[1] * b[0];
}
double acose( const double arg)
{
    if( arg >= 1.)
        return( 0.);
    if( arg <= -1.)
        return( pi);
    return( acos( arg));
}

double asine( const double arg)
{
    if( arg >= 1.)
        return( pi / 2);
    if( arg <= -1.)
        return( -pi / 2.);
    return( asin( arg));
}

-(int) AL_calc_classical_elements:(double *) r:(double) t:(int) ref:(double) gm {
    const double *v = r + 3;
    const double r_dot_v = AL_dot_product( r, v);
    const double dist = AL_vector3_length( r);
    const double v2 = AL_dot_product( v, v);
    const double inv_major_axis = 2. / dist - v2 / gm;
    double h0, n0, tval;
    double h[3], e[3], ecc2;
    double ecc;
    int i;
    
    AL_vector_cross_product( h, r, v);
    n0 = h[0] * h[0] + h[1] * h[1];
    h0 = n0 + h[2] * h[2];
    n0 = sqrt( n0);
    h0 = sqrt( h0);
    
    /* See Danby,  p 204-206,  for much of this: */
    if( ref & 1)
    {
        [ minor_planet LongitudeAscendingNode:atan2( h[0], -h[1]) ];
        [ minor_planet Inclination:asine( n0 / h0) ];
        if( h[2] < 0.)                   /* retrograde orbit */
            [ minor_planet Inclination:PI - [minor_planet Inclination] ];
    }
    AL_vector_cross_product( e, v, h);
    for( i = 0; i < 3; i++)
        e[i] = e[i] / gm - r[i] / dist;
    tval = AL_dot_product( e, h) / h0;     /* "flatten" e vector into the rv */
    for( i = 0; i < 3; i++)             /* plane to avoid roundoff; see   */
        e[i] -= h[i] * tval;             /* above comments                 */
    ecc2 = AL_dot_product( e, e);
    [minor_planet MinorToMajor:sqrt( fabs( 1. - ecc2) ) ];
    [minor_planet Eccentricity:sqrt(ecc2) ];
    ecc = [minor_planet Eccentricity ];
       
    if( !ecc)                     /* for purely circular orbits,  e is */
    {                          /* arbitrary in the orbit plane; choose */
        for( i = 0; i < 3; i++)    /* r normalized                         */
            e[i] = r[i] / dist;
    }
    else                           /* ...and if it's not circular,  */
        for( i = 0; i < 3; i++)     /* normalize e:                  */
            e[i] /= ecc;
    if( inv_major_axis)
    {
        [ minor_planet MajorAxis: 1./inv_major_axis ]; 
        [ minor_planet T0:([minor_planet MajorAxis ] * sqrt( fabs( [minor_planet MajorAxis ]) / gm ))];
    }
    
    if( ecc < .9)
        [minor_planet Q:([minor_planet MajorAxis ] * (1. - ecc ) ) ];
    else        /* at eccentricities near one,  the above suffers  */
    {        /* a loss of precision problem,  and we switch to: */
        const double gm_over_h0 = gm / h0;
        const double perihelion_speed = gm_over_h0 +
        sqrt( gm_over_h0 * gm_over_h0 - inv_major_axis * gm);
        
        [minor_planet Q:( h0 / perihelion_speed ) ];
    }
    
    AL_vector_cross_product( [minor_planet Sideways] , h, e);
    /* At this point,  elem->sideways has length h0.  */
    if( ref & 1)
    {
        const double cos_arg_per = (h[0] * e[1] - h[1] * e[0]) / n0;
        
        if( cos_arg_per < .7 && cos_arg_per > -.7)
            [minor_planet ArgumentPerihelion:acos( cos_arg_per) ];
        else
        {
            const double sin_arg_per = (e[0] * h[0] * h[2] + e[1] * h[1] * h[2] - e[2] * n0 * n0) / (n0 * h0);
            [minor_planet ArgumentPerihelion:fabs( asin( sin_arg_per)) ];
            if( cos_arg_per < 0.)
                [minor_planet ArgumentPerihelion:(PI-[minor_planet ArgumentPerihelion]) ];
                
        }
        if( e[2] < 0.)
            [minor_planet ArgumentPerihelion:(PI + PI - [minor_planet ArgumentPerihelion] )];
    }
    
    if( inv_major_axis > 0.)         /* elliptical case */
    {
        const double r_cos_true_anom = AL_dot_product( r, e);
        const double r_sin_true_anom = AL_dot_product( r, [minor_planet Sideways] ) / h0;
        const double cos_E = r_cos_true_anom * inv_major_axis + ecc;
        const double sin_E = r_sin_true_anom * inv_major_axis / [minor_planet MinorToMajor ];
        const double ecc_anom = atan2( sin_E, cos_E);
        
        [ minor_planet MeanAnomaly:(ecc_anom - ecc * sin( ecc_anom)) ];
        [ minor_planet PerihelionTime:(t - [minor_planet MeanAnomaly ] * [minor_planet T0] ) ];
    }
    else if( inv_major_axis < 0.)         /* hyperbolic case */
    {
        const double z = (1. - dist * inv_major_axis) / ecc;
        double f = log( z + sqrt( z * z - 1.));
        
        if( r_dot_v < 0.)
            f = -f;
        [minor_planet MeanAnomaly:(ecc * sinh(f)-f) ];
        [minor_planet PerihelionTime:(t-[minor_planet MeanAnomaly] * fabs([minor_planet T0] ))];
        h0 = -h0;
    }
    else              /* parabolic case */
    {
        double tau;
        
        tau = sqrt( dist / [minor_planet Q] - 1.);
        if( r_dot_v < 0.)
            tau = -tau;

        [minor_planet W0:( (3. / SQRT_2) / ([minor_planet Q] * sqrt( [minor_planet Q] / gm))) ];
        /*    elem->perih_time = t - tau * (tau * tau / 3. + 1) *                   */
        /*                                      elem->q * sqrt( 2. * elem->q / gm); */
        [minor_planet PerihelionTime:(t - tau * (tau * tau / 3. + 1) * 3. / [minor_planet W0])];
    }
    double * sideways = [minor_planet Sideways ];
    
    [minor_planet PerihelionVector:e ];
    
    for( i = 0; i < 3; i++)
        sideways[i] /= h0;
    [minor_planet AngularMomentum:h0 ];
    return( 0);
}


-(int) AL_IntegrateOrbit:(double) jd_from:(double) jd_to:(double) max_err: (int) n_steps :(double *)pv {
    double delta[6],  posnvel[6], stepsize = (jd_to - jd_from) / (double)n_steps;
    double curr_jd = jd_from;
    int i, j;
    
    NSLog( @"AL_IntegrateOrbit:  0" );
    position_cache = [self AL_MakePositionCache:jd_from :(jd_to-jd_from)/n_steps :n_steps ];
    NSLog( @"AL_IntegrateOrbit:  1" );
    for( i = 0; i < 6; i++)
        delta[i] = 0.;
    NSLog( @"AL_IntegrateOrbit:  2" );
    for( i = 0; i < n_steps; i++)
    {
        double new_delta[6];
        int chickened_out;
     
        chickened_out = [ self AL_FullRKStep:delta : new_delta :curr_jd :curr_jd + stepsize :max_err ];
     //   chickened_out = full_rk_step( elem, delta, new_delta, curr_jd,
     //                                curr_jd + stepsize, max_err);
        memcpy( delta, new_delta, 6 * sizeof( double));
        curr_jd += stepsize;
        if( i && (i % resync_freq == 0 || chickened_out))
        {
            [self AL_comet_posn_and_vel:curr_jd :posnvel :posnvel + 3 ];
            for( j = 0; j < 6; j++)
            {
                posnvel[j] += delta[j];
                delta[j] = 0.;
            }
            [ minor_planet Epoc:curr_jd ];
            [self AL_calc_classical_elements:posnvel :curr_jd :1 :SOLAR_GM ];
        }
    }
    NSLog( @"AL_IntegrateOrbit:  3" );
    [ self AL_comet_posn_and_vel:jd_to :posnvel :posnvel + 3 ];
    for( i = 0; i < 6; i++)
        posnvel[i] += delta[i];
    
    NSLog( @"AL_IntegrateOrbit:  4" );
    for(i=0;i<6;i++)
        pv[i] = posnvel[i];
    NSLog( @"AL_IntegrateOrbit:  5" );
    [ minor_planet Epoc:jd_to ];
    [ self AL_calc_classical_elements:posnvel :jd_to :1 :SOLAR_GM ];
    NSLog( @"AL_IntegrateOrbit:  6" );
    
    
    return(0);
}

@end