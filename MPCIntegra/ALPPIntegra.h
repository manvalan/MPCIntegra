//
//  ALPPIntegra.h
//  MPCIntegra
//
//  Created by Michele Bigi on 20/01/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "ALJplEphem.h"
#import "ALMinorPlanet.h"

#ifndef AU_IN_KM
#define AU_IN_KM 1.49597870691e+8
#endif

#ifndef AU_IN_METERS
#define AU_IN_METERS (AU_IN_KM * 1000.)
#endif

#ifndef SPEED_OF_LIGHT
#define SPEED_OF_LIGHT 299792.458
#endif

#ifndef AU_PER_DAY
#define AU_PER_DAY (86400. * SPEED_OF_LIGHT / AU_IN_KM)
#endif

#ifndef JPL_EPHEM_EARTH_MOON_RATIO
#define JPL_EPHEM_EARTH_MOON_RATIO      36
#endif

#ifndef PI
#define PI 3.14159265358979323
#endif

#ifndef GAUSS_K
#define GAUSS_K .01720209895
#endif

#ifndef SOLAR_GM
#define SOLAR_GM (GAUSS_K * GAUSS_K)
#endif

#ifndef PERTURBERS_MERCURY_TO_NEPTUNE
#define PERTURBERS_MERCURY_TO_NEPTUNE 0xff
#define PERTURBERS_PLUTO 0x100
#define PERTURBERS_MOON  0x200
#define PERTURBERS_PLUTO_AND_MOON (PERTURBERS_PLUTO | PERTURBERS_MOON)
#define PERTURBERS_CERES_PALLAS_VESTA 0x1c00
#define N_PERTURBERS 10
#define EARTH_MOON_RATIO 81.30056
#endif

//typedef long double  double;


@interface 

ALPPIntegra : NSObject {
    ALJplEphem * JPLEphem;
    ALMinorPlanet * minor_planet;
}

- (id) initWithMinorPlanet:(ALJplEphem *) ephem:(ALMinorPlanet *) mplanet;

-(int) AL_IntegrateOrbit:(double) jd_from:(double) jd_to:(double) max_err: (int) n_steps:(double *)pv;
-(int) AL_calc_classical_elements:(double *) r:(double) t:(int) ref:(double) gm;
-(int) AL_FullRKStep:(double *) ivals:(double *) ovals: (double) t0: (double) t1: (double) max_err;
-(int) AL_TakeStep:(double) jd: (double *) ival: (double *) ovals: (double *)errs :(int) n_vals: (double) step_size;
-(void) AL_computeDerivatives:(double) jd:(double *) delta:(double *) derivs:(double *) posn_data;
-(int)  AL_comet_posn_and_vel:(double) t:(double *) loc:(double *) vel;
-(void)  AL_comet_posn_part_ii:(double) t:(double *) loc:(double *) vel ;

@end
