//
//  ALJplEphem.m
//  MPCIntegra
//
//  Created by Michele Bigi on 04/07/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import "ALJplEphem.h"


@implementation ALJplEphem

- (id)init
{
    self = [super init];
    if (self) {
        // Initialization code here.
    }
    
    return self;
}

- (id)initWithFile:(NSString *) file
{
    self = [super init];
    if (self) {
        // Initialization code here.
        [ self AL_setJplEphData: [ self AL_openJplEphemeride:file ]];
    }
    
    return self;  
}

- ( struct jpl_eph_data * ) AL_openJplEphemeride: (NSString * ) file
{
    int ncon;
    char costanti_nome[256][6];
    double costanti_val[256];
    
    return  jpl_init_ephemeris( [ file UTF8String ] , costanti_nome , costanti_val , &ncon );
    /*
    
    struct jpl_eph_data * jple = jpl_init_ephemeris( [ file UTF8String ] , costanti_nome , costanti_val , &ncon );

    return  jple;*/
}

- (void)dealloc
{
    jpl_close_ephemeris( eph_data );
    [super dealloc];
}

- (struct jpl_eph_data *) AL_getJplEphData {
    return eph_data;
}

- ( void ) AL_setJplEphData:( struct jpl_eph_data * ) data
{
    eph_data = data;
}

- (long) AL_EphemerideVersion {
    return jpl_get_long( eph_data, JPL_EPHEM_EPHEMERIS_VERSION);
}

- (double) AL_EphemerideStart {
    return jpl_get_double( eph_data, JPL_EPHEM_START_JD);
}

- (double) AL_EphemerideEnd {
    return jpl_get_double( eph_data, JPL_EPHEM_END_JD);
}

- (double) AL_Stepsize {
    return jpl_get_double( eph_data, JPL_EPHEM_STEP);
}


- (double) AL_AU {
    return jpl_get_double( eph_data, JPL_EPHEM_AU_IN_KM);
}

- (double) AL_EarthMoonRatio{
    return  jpl_get_double( eph_data, JPL_EPHEM_EARTH_MOON_RATIO);
}

- (long) AL_KernelSize {
    return jpl_get_long( eph_data, JPL_EPHEM_KERNEL_SIZE);
}

- (long) AL_RecordSize {
    return  jpl_get_long( eph_data, JPL_EPHEM_KERNEL_RECORD_SIZE);
}

- (long) AL_NCoeff {
    return jpl_get_long( eph_data, JPL_EPHEM_KERNEL_NCOEFF);
}

- (long) AL_SwapBytes {
    return jpl_get_long( eph_data, JPL_EPHEM_KERNEL_SWAP_BYTES);
}

- ( void ) AL_StateVec:(double) jd : (int) corp: (int) rif : (double *) sv {
    jpl_pleph( eph_data, jd, rif, corp, sv, 1);
}

- (void) AL_GeoStateVec:(double) jd : (int) corp: (double *) sv {
    [ self AL_StateVec:jd :corp :3 :sv ];
}

- (void) AL_ElioStateVec:(double) jd : (int) corp: (double *) sv {
    [ self AL_StateVec:jd :corp :11 :sv ];
}

- (Vec3D *) AL_Vector:(double) jd : (int) corp: (int) rif {
    double sv[6];
    
    jpl_pleph( eph_data, jd, rif, corp, sv, 0);
    Vec3D * pVec = [[ Vec3D alloc] init: sv[0] : sv[1] : sv[2] ];
    return pVec;
}

- (Vec3D *) AL_GeoVector:(double) jd : (int) corp {
    double sv[6];
    
    jpl_pleph( eph_data, jd, 3, corp, sv, 0);
    Vec3D * pVec = [[ Vec3D alloc] init: sv[0] : sv[1] : sv[2] ];
    return pVec;
}

- (Vec3D *) AL_ElioVector:(double) jd : (int) corp {
    double sv[6];
    
    jpl_pleph( eph_data, jd, 11, corp, sv, 0);
    Vec3D * pVec = [[ Vec3D alloc] init: sv[0] : sv[1] : sv[2] ];
    return pVec;
}

- (double *) AL_ElioVector2:(double) jd : (int) corp {
    static double sv[6];
    
    jpl_pleph( eph_data, jd, 11, corp, sv, 0);
    
    return sv;
}

- (Vec3D *) AL_Vector:(double) jd : (int) corp: (int) rif: (Vec3D * ) Vel {
    double sv[6];
    
    jpl_pleph( eph_data, jd, rif, corp, sv, 1);
    Vec3D * pVec = [[ Vec3D alloc] init: sv[0] : sv[1] : sv[2] ];
    [ Vel set: sv[3] : sv[4] : sv[5] : NO ];
    return pVec;
}

- (Vec3D *) AL_GeoVector:(double) jd : (int) corp : (Vec3D * ) Vel {
    double sv[6];
    
    jpl_pleph( eph_data, jd, 3, corp, sv, 1);
    Vec3D * pVec = [[ Vec3D alloc] init: sv[0] : sv[1] : sv[2] ];
    [ Vel set: sv[3] : sv[4] : sv[5] : NO ];
    return pVec;
}

- (Vec3D *) AL_ElioVector:(double) jd : (int) corp : (Vec3D * ) Vel {
    double sv[6];
    
    jpl_pleph( eph_data, jd, 11, corp, sv, 1);
    Vec3D * pVec = [[ Vec3D alloc] init: sv[0] : sv[1] : sv[2] ];
    [ Vel set: sv[3] : sv[4] : sv[5] : NO ];
    return pVec;
}

- (void) JPLState:(double) jd:(int*)list:(double[][6]) cache {
    
    //double cache_t[12][6];
    int i,j;
    
    jpl_state(eph_data, jd, list, cache, NULL, 0);
    /*
    for(i=0;i<12;i++)
        for(j=0;j<6;j++)
            cache[i][j] = (long double)cache_t[i][j];
    */
    
}


@end
