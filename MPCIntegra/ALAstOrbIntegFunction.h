//
//  LAAstOrbIntegFunction.h
//  MPCIntegra
//
//  Created by Michele Bigi on 10/07/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "ALJplEphem.h"

#ifndef pi2
#define pi 3.14159265358979324
#define pi2 2.0*pi
#define Rad pi/180.0
#define Deg 180.0/pi
#define Arcs    3600*Deg

#define R_Earth 6378.137
#define R_Sun   696000.0
#define R_Moon  1738.9
#define MJD     2400000.5
#define MJD_2000    (MJD+51544) // DA SISTEMARE NEL CODICE!!! Concettualmente errato!
#define JD_2000     (2400000.5+51544) // DA SISTEMARE NEL CODICE!!! Concettualmente errato!
#define kGauss  0.01720209895 //0.01720209895
#define GM_Sun  (kGauss *kGauss)
#define AU      149597870.0
#define c_light 173.14 // in [AU / d ]
#endif
// Radii of Earth, Sun and Moon
/*const double R_Earth   =   6378.137;     // [km]
const double R_Sun     = 696000.0;       // [km]
const double R_Moon    =   1738.0;       // [km]

const double MJD_J2000 = 51544.5;        // MJD of Epoch J2000.0
const double T_J2000   =  0.0;           // Epoch J2000.0
const double T_B1950   = -0.500002108;   // Epoch B1950

const double kGauss    = 0.01720209895;  // gravitational constant
const double GM_Sun    = kGauss*kGauss;  // [AU^3/d^2]

const double AU        = 149597870.0;    // Astronomical unit [km]

const double c_light   = 173.14;         // speed of light [AU/d]
*/

enum PlanetType { 
    Sun = 0,
    Mercury = 1,
    Venus = 2,
    Earth,
    Mars,
    Jupiter,
    Saturn,
    Uranus,
    Neptune,
    Pluto,
    Moon,
    Ceres,
    Pallas,
    Vesta
};


@interface ALAstOrbIntegFunction : NSObject {
@private
    int     Neqn;              // Number of differential eqns.
    double GM[10];
    ALJplEphem * JPLEphem;
    
}

- (id)initWithJPLEphem:(ALJplEphem *) ephem;
- (Vec3D *) AL_AccelJPL:( const double )Mjd: (Vec3D *) r ;
- (long double) AL_GM:(int) iPlanet;
- (void) AL_Run:( double )X: (double *) Y:(double *) dYdX;
- (int) AL_getNEQN;

@end
