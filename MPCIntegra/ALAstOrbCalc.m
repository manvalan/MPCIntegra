//
//  ALAstOrbCalc.m
//  MPCIntegra
//
//  Created by Michele Bigi on 10/07/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import "ALAstOrbCalc.h"
#import "ALJplEphem.h"
#import "AL_DE.h"
#import "Math3D.h"
#import "Vec3D.h"
#import "gal.h"
#import "jpl_int.h"
#import "jpleph.h"

@implementation ALAstOrbCalc

- (id)initWithPatameter:(double) ObservationEpoch:(Vec3D *) r:(Vec3D *) v : (NSString *) JPLDE405File
{
    self = [super init];
    if (self) {
        // Initialization code here.
        
        jplEphem = [[ ALJplEphem alloc ] initWithFile: JPLDE405File ];
        IntegFunction = [[ ALAstOrbIntegFunction alloc] initWithJPLEphem:jplEphem ];
        IntegMethod = [[ AL_DE alloc] initDE: IntegFunction ];
        ETMjdCurrent = ObservationEpoch;
     
        Y = malloc( sizeof(double) * ([IntegFunction AL_getNEQN] + 1 ) );
        
        
        Y[1] = [r getElement:0];
        Y[2] = [r getElement:1];
        Y[3] = [r getElement:2];
        Y[4] = [v getElement:0];
        Y[5] = [v getElement:1];
        Y[6] = [v getElement:2];
        Y[0] = 0.0;
    }
    
    return self;
}

- (void)dealloc
{
    [super dealloc];
    free( Y );
    [ IntegFunction dealloc ];
    [ IntegMethod dealloc ];
    [ jplEphem dealloc ];
}

- (int) AL_AOC_Integrate:(double) End
{
    DE_STATE State  = DE_INIT;
    double relerr = def_eps;
    double abserr = def_abserr;
    //int i;
    
    
    do { 
        
       [ IntegMethod AL_DE_Integ:Y :&ETMjdCurrent :End : &relerr : &abserr :&State : YES ];
        NSLog(@"Integrate : %lf\n" , ETMjdCurrent );
        NSLog(@"Y: %lf %lf %lf\n" , Y[1],Y[2], Y[3]);
        //IntegMethod->Integ( Y, ETMjdCurrent, End, relerr, abserr, State );
        
        if( State == DE_INVALID_PARAMS ) { 
            //pModule->ErrorMessage( APS_ASTORBCALC_PARAM );
            return( 1 );
        }
        NSLog(@"State :%i\n" , State );
    } while( State > DE_DONE );
    
    /**************/
    NSLog( @">>> %lf     %lf  %lf  %lf     %lf  %lf  %lf \n", End, Y[1], Y[2], Y[3], Y[4], Y[5], Y[6] );
    
    
    return( 0 );
}

- (Vec3D *) AL_AOC_GetR
{
    return [[Vec3D alloc] init:Y[1] :Y[2]:Y[3] ];
}
 
- (Vec3D *) AL_AOC_GetV
{
        return [[Vec3D alloc] init:Y[4] :Y[5]:Y[6] ];
}


/*

//======================= APSAstOrbCalc ==========================

APSAstOrbCalc :: APSAstOrbCalc( APSSubModule * pAPSSubModule,
                               const APS_INTEGRATION_TYPE aIntegType,
                               const APSAstOrbIntegFunction * IntegFunction,
                               const double ObservationEpoch,
                               const APSVec3d & r,
                               const APSVec3d & v,
                               const double a_eps,
                               const double a_abserr ) :
IntegType( aIntegType ), ETMjdCurrent( ObservationEpoch ),
eps( a_eps ), abs_err_val( a_abserr )
{
   }

APSAstOrbCalc :: ~APSAstOrbCalc( void )
{
    delete [] Y;
    delete IntegMethod;
    delete pModule;
}

int APSAstOrbCalc :: Integrate( const double End )
{
    apsmathlib::DE_STATE State  = apsmathlib::DE_INIT;
    double               relerr = eps;
    double               abserr = abs_err_val;
    
    
    
    do { 
        IntegMethod->Integ( Y, ETMjdCurrent, End, relerr, abserr, State );
        
        if( State == apsmathlib::DE_INVALID_PARAMS ) { 
            pModule->ErrorMessage( APS_ASTORBCALC_PARAM );
            return( 1 );
        }
        
    } while( State < apsmathlib::DE_DONE );
    return( 0 );
}

APSVec3d APSAstOrbCalc :: GetR( void ) const
{
    return( APSVec3d( Y[ 1 ], Y[ 2 ], Y[ 3 ] ) );
}

APSVec3d APSAstOrbCalc :: GetV( void ) const
{
    return( APSVec3d( Y[ 4 ], Y[ 5 ], Y[ 6 ] ) );
}

*/@end
