//
//  ALMinorPlanet.h
//  MPCIntegra
//
//  Created by Michele Bigi on 05/07/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>

#define ASTORB_RECORD_LEN 268
#define GAUSS_K .01720209895
#define SOLAR_GM (GAUSS_K * GAUSS_K)
#define SQRT_2 1.414213562
#define CALENDAR_GREGORIAN             0
#define CALENDAR_JULIAN                1
#define CALENDAR_JULIAN_GREGORIAN      6
#define JUL_GREG_CALENDAR_EPOCH 1721060L

#ifndef PI
#define PI 3.14159265358979323    
#endif

@interface ALMinorPlanet : NSObject {
@private
    //ELEMENTS  mp_elem;
    
    double perih_time, q, ecc, incl, arg_per, asc_node;
    double epoch,  mean_anomaly, diameter, ceu, dceu;
    /* derived quantities: */
    double lon_per, minor_to_major;
    double perih_vec[3], sideways[3];
    double angular_momentum, major_axis, t0, w0;
    double abs_mag, slope_param;
    int is_asteroid, central_obj;
    double brightness, slope;
    
    long num;
    double jd;
    double ra, dec;
    double magn;
    double dist;
    double earth_distance;
    double elong;
    
    
}

- (double) PerihelionTime;
- (double) Q;
- (double) Eccentricity;
- (double) Inclination;
- (double) ArgumentPerihelion;
- (double) LongitudeAscendingNode;
- (double) MeanAnomaly;
- (double) Diameter;
- (double) Brightness;
- (double) T0;
- (double) Slope;
- (double) W0;
- (double) MinorToMajor;
- (double *) PerihelionVector;
- (double *) Sideways;
- (double) AngularMomentum;
- (double) MajorAxis;

- (void) PerihelionTime:(double) t_val;
- (void) Q:(double) t_val;
- (void) Eccentricity:(double) t_val;
- (void) Inclination:(double) t_val;
- (void) ArgumentPerihelion:(double) t_val;
- (void) LongitudeAscendingNode:(double) t_val;
- (void) MeanAnomaly:(double) t_val;
- (void) Diameter:(double) t_val;
- (void) Brightness:(double) t_val;
- (void) T0:(double) t_val;
- (void) Slope:(double) t_val;
- (void) W0:(double) ttemp; 
- (void) MinorToMajor:(double) ttemp;
- (void) PerihelionVector:(double *) t_val;
- (void) Sideways:(double *) t_val;
- (void) AngularMomentum:(double) t_val;
- (void) MajorAxis:(double) ttemp;

- (id) initWithAstOrb:(NSString *) AstOrbPath : (long) numb;
- (void) ReadAstOrb:(NSString *) AstOrbPath : (long)asteroid_no;
- (void) AL_ExtractAstOrbData:(char *) buff;
- (void) AL_DoRemainingSetup;
- (void) AL_DerivedQuantities:(double) gm;
- (void) AL_SetupOrbitVectors;

- ( double ) getDiameter ;
- ( double ) getDist ;
- ( double ) getEarthDist;
- ( double ) getRA;
- ( double ) getDec;
- ( double ) getDCEU;
- ( double ) getCEU;
- (double) Epoc;
- (void) Epoc:(double) eep;
@end
