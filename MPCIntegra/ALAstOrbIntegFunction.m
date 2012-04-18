//
//  LAAstOrbIntegFunction.m
//  MPCIntegra
//
//  Created by Michele Bigi on 10/07/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import "ALAstOrbIntegFunction.h"
#import "Vec3D.h"
#import "Math3D.h"
#import "costants.h"
//#import "APC_Const.h"
#include "sofa.h"
#include "sofam.h"
#import "gal.h" 

#define EARTH_MOON_RATIO 81.30056
#define GAUSS_K .01720209895

@implementation ALAstOrbIntegFunction

- (id)initWithJPLEphem:(ALJplEphem *) ephem
{
    self = [super init];
    if (self) {
        // Initialization code here.
        Neqn = 6;
        JPLEphem = ephem;
    }
    
    return self;
}

- (void)dealloc
{
    [super dealloc];
}

- (int) AL_getNEQN{
    return  Neqn;
}

- (long double) AL_GM:(int) iPlanet
{
    static const double  tGM[10] = { 
        GM_Sun,                 // Sun
        GM_Sun / 6023600.0,     // Mercury
        GM_Sun /  408523.5,     // Venus
        GM_Sun /  328900.5,     // Earth
        GM_Sun / 3098710.0,     // Mars
        GM_Sun /    1047.355,   // Jupiter
        GM_Sun /    3498.5,     // Saturn
        GM_Sun /   22869.0,     // Uranus
        GM_Sun /   19314.0,     // Neptune
        GM_Sun / 3000000.0      // Pluto
    };
    
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
    
    return relative_mass[iPlanet] * GAUSS_K*GAUSS_K;
}


double norm( double * a )
{
    return  sqrt( a[0]*a[0]+a[1]*a[1]+a[2]*a[2] );
}

void eq2ecl(double t , double a[3][3] ) {
    //static double a[3][3];
    
    double eps = ( 23.43929111-(46.8150+(0.00059-0.001813*t)*t)*t/3600)*PI/180.0;
    //Math3D * eq2ecl = [[ Math3D alloc] initWithRotationX: eps];
    
    gal_rx( eps, a );
    return;
}




double q_norm( double * a )
{
    return  sqrt( a[0]*a[0]+a[1]*a[1]+a[2]*a[2] );
}

void CalcDerivate( long double GM, double * d, long double * a )
{
    long double rt1[3];
    double D = q_norm( d );
    
    rt1[0] = (long double)kGauss * (long double)kGauss *((-(long double)GM ) * (long double )d[0] )  / (long double)(D*D*D);
    rt1[1] = (long double)kGauss * (long double)kGauss *((-(long double)GM ) * (long double )d[1] )  / (long double)(D*D*D);
    rt1[2] = (long double)kGauss * (long double)kGauss *((-(long double)GM ) * (long double )d[2] )  / (long double)(D*D*D);
    a[0] += (double) rt1[0]; a[1] += (double) rt1[1]; a[2] += (double) rt1[2];
}

- (Vec3D *) AL_AccelJPL:( const double )Mjd: (Vec3D *) pr  
{
    int      iPlanet;
    static long double p_GM[14] = { 1.,
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
       
    double r[3] , r_t[6], d_d[3];
    long double a[3];
                  
    double eq2ecl_m[3][3];
    
    double dt = ( Mjd - MJD_2000 ) / 36525.0 ;
    double eps = ( 23.43929111-(46.8150+(0.00059-0.001813*dt)*dt)*dt/3600)*PI/180.0;
    
    double cosE = cos( eps );
    double sinE = sin( eps );
    eq2ecl_m[0][0] = 1; eq2ecl_m[0][1] = 0; eq2ecl_m[0][2] = 0;
    eq2ecl_m[1][0] = 0; eq2ecl_m[1][1] = cosE; eq2ecl_m[1][2] = sinE;
    eq2ecl_m[2][0] = 0; eq2ecl_m[2][1] = -sinE; eq2ecl_m[2][2] = cosE;
    
    r[0] = [pr getElement:0]; r[1] = [pr getElement:1]; r[2] = [pr getElement:2];
       
    a[0] = a[1] = a[2] = 0.0;
    
    CalcDerivate(p_GM[0], r, a );
    //NSLog( @"Accelleration by Sun: %Lg   %Lg   %Lg\n", a[0] , a[1], a[2] );
    
    for ( iPlanet = Mercury; iPlanet <= Moon; iPlanet++ ) {
               
        double * r_planet = [ JPLEphem AL_ElioVector2:Mjd :iPlanet];
        if( iPlanet == Moon ) {
            double * r_earth = [ JPLEphem AL_ElioVector2:Mjd :Earth];
            r_planet[0] += r_earth[0];
            r_planet[1] += r_earth[1];
            r_planet[2] += r_earth[2];
        }
        gal_rxp(eq2ecl_m, r_planet, r_t );
        
        d_d[0] = r_t[0] - r[0];
        d_d[1] = r_t[1] - r[1];
        d_d[2] = r_t[2] - r[2];
        
        CalcDerivate( p_GM[iPlanet] , d_d , a );
        CalcDerivate( p_GM[iPlanet] , r_t , a );
    }
    
    return [ [ Vec3D alloc ] init: (double)a[0] :(double)a[1] :(double)a[2] ];
}

-(void) AL_Run:( double )X: (double *) Y:(double *) dYdX  
{
    
    Vec3D * r = [[ Vec3D alloc] init: Y[ 1 ]: Y[ 2 ]: Y[ 3 ]];
    Vec3D * a = [ self AL_AccelJPL:X :r ];

    dYdX[0] = 0.0;
    dYdX[1] = Y[ 4 ];  // velocity
    dYdX[2] = Y[ 5 ];  // velocity
    dYdX[3] = Y[ 6 ];  // velocity
    dYdX[4] = [ a getElement: 0 ];  // acceleration
    dYdX[5] = [ a getElement:1 ];  // acceleration
    dYdX[6] = [ a getElement:2 ];  // acceleration
     
}


@end
