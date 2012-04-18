//
//  ALJplEphem.h
//  MPCIntegra
//
//  Created by Michele Bigi on 04/07/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "jpl_int.h"
#import "jpleph.h"
#import "Vec3D.h"
#import "Math3D.h"
//#import "ALPPIntegra.h"

@interface ALJplEphem : NSObject {
@private
    struct jpl_eph_data *eph_data;
    
    
}

- (id)init; 
- (id)initWithFile:(NSString *) file;

- ( struct jpl_eph_data * ) AL_getJplEphData;
- ( void ) AL_setJplEphData:( struct jpl_eph_data * ) data;
- ( struct jpl_eph_data * ) AL_openJplEphemeride: (NSString * ) file;

- (long) AL_EphemerideVersion ;
- (double) AL_EphemerideStart ;
- (double) AL_EphemerideEnd ;
- (double) AL_Stepsize;
- (double) AL_AU ;
- (double) AL_EarthMoonRatio;
- (long) AL_KernelSize;
- (long) AL_RecordSize ;
- (long) AL_NCoeff;
- (long) AL_SwapBytes;
- (void) AL_StateVec:(double) jd : (int) corp: (int) rif : (double *) sv ;
- (void) AL_GeoStateVec:(double) jd : (int) corp: (double *) sv;
- (void) AL_ElioStateVec:(double) jd : (int) corp: (double *) sv;

- (Vec3D *) AL_Vector:(double) jd : (int) corp: (int) rif;
- (Vec3D *) AL_GeoVector:(double) jd : (int) corp ;
- (Vec3D *) AL_ElioVector:(double) jd : (int) corp;
- (double *) AL_ElioVector2:(double) jd : (int) corp;
- (Vec3D *) AL_Vector:(double) jd : (int) corp: (int) rif: (Vec3D * ) Vel;
- (Vec3D *) AL_GeoVector:(double) jd : (int) corp : (Vec3D * ) Vel ;
- (Vec3D *) AL_ElioVector:(double) jd : (int) corp : (Vec3D * ) Vel ;
- (void) JPLState:(double) jd:(int*)list:(  double[][6]) cache;

@end
