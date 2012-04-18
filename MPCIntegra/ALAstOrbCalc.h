//
//  ALAstOrbCalc.h
//  MPCIntegra
//
//  Created by Michele Bigi on 10/07/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "ALJplEphem.h"
#import "AL_DE.h"
#import "ALAstOrbIntegFunction.h"
#import "Vec3D.h"
#import "Math3D.h"

typedef enum {
    APS_INTEGRATION_DE = 0
} APS_INTEGRATION_TYPE;

static const double   def_eps    = 1.0e-10; // Relative accuracy
static const double   def_abserr = 1.0e-10; // Absolute accuracy

@interface ALAstOrbCalc : NSObject {
@private
  
    
    //APSModuleAstOrbCalc        * pModule;
    const APS_INTEGRATION_TYPE   IntegType;
    //APSDE                      * IntegMethod;
    AL_DE * IntegMethod;
    ALAstOrbIntegFunction * IntegFunction;
    ALJplEphem * jplEphem;
    
    double * Y;
    double ETMjdCurrent;
    double eps;
    double abs_err_val;
}

- (id)initWithPatameter:(double) ObservationEpoch:(Vec3D *) r:(Vec3D *) v: (NSString *) JPLDE405File;
- (int) AL_AOC_Integrate:(double) End;
- (Vec3D *) AL_AOC_GetR;
- (Vec3D *) AL_AOC_GetV;

@end
