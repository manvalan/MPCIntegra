//
//  MPCIntegra.m
//  MPCIntegra
//
//  Created by Michele Bigi on 04/07/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import "MPCIntegra.h"
#import "Math3D.h"
#import "ALJplEphem.h"
#import "ALMinorPlanet.h"
#import "ALPPIntegra.h"
#include "sofa.h"
#include "sofam.h"
#import  "gal.h"

@implementation MPCIntegra

- (id)init
{
    self = [super init];
    if (self) {
        // Add your subclass-specific initialization here.
        // If an error occurs here, send a [self release] message and return nil.
    }
    return self;
}

- (NSString *)windowNibName
{
    // Override returning the nib file name of the document
    // If you need to use a subclass of NSWindowController or if your document supports multiple NSWindowControllers, you should remove this method and override -makeWindowControllers instead.
    return @"MPCIntegra";
}

- (void)windowControllerDidLoadNib:(NSWindowController *)aController
{
    [super windowControllerDidLoadNib:aController];
    NSLog(@"inizio" );
    
    // Add any code here that needs to be executed once the windowController has loaded the document's window.
    NSString * file = @"/Users/iMike/Desktop/Occultation/unxp2000N.405.405";
    NSString * astOrbFile = @"/Users/iMike/Desktop/Occultation/astorb.dat";
    NSLog( @"DE405 Ephem file: %@\n", file );
    NSLog( @"AstOrb file:      %@\n", astOrbFile );
    
    long mp_numb = 1;
    double pv[6];
    
    ALJplEphem * Eff = [[ ALJplEphem alloc] initWithFile:file ];
    NSLog( @"MPCIntegra - Check1\n" );
   
    int i;
    for(i=0;i<14;i++) {
        Vec3D * ev = [ Eff AL_ElioVector:2455786.5 : 1 ];
        [ev toString ];
    }
    ALMinorPlanet * mplanet = [[ALMinorPlanet alloc] initWithAstOrb:astOrbFile :mp_numb ];
    
    NSLog( @"MPCIntegra - Check2\n" );
    
    ALPPIntegra * integrator = [[ALPPIntegra alloc] initWithMinorPlanet:Eff :mplanet];
 
    
    NSLog( @"MPCIntegra - Check3\n" );
    double jd_calc = [mplanet Epoc]+1 ;//2455786.5;
    double posnvel[6] = { 2.563074267657057,  -1.437309950582643E+00, -5.171585463325198E-01,
        4.599720288348380E-03,  8.358337188760315E-03, -5.867377166799321E-04};
    
    //[integrator AL_calc_classical_elements:posnvel :  jd_calc  :1 :SOLAR_GM ];
    NSLog( @"Start Epoc: %lf\n",[mplanet Epoc] );
    NSLog( @"Eccentricity            : %lf\n",[mplanet Eccentricity] );
    NSLog( @"Mean Anomaly            : %lf\n",[mplanet MeanAnomaly] );
    NSLog( @"Longitude Ascending Node: %lf\n",[mplanet LongitudeAscendingNode] );
    [integrator AL_IntegrateOrbit : [mplanet Epoc] : jd_calc : 1e-5 : 10: pv ];
    NSLog( @"MPCIntegra - Check4\n" );
    
    NSLog(@"Epoc: %lf\n%lf\t%lf\t%lf\n%lf\t%lf\t%lf\n" , [ mplanet Epoc], pv[0],pv[1],pv[2], pv[3],pv[4],pv[5]);
  
    [ integrator dealloc];
     
    
    [ mplanet dealloc];
    [ Eff dealloc ];
    
}

- (NSData *)dataOfType:(NSString *)typeName error:(NSError **)outError {
    /*
     Insert code here to write your document to data of the specified type. If outError != NULL, ensure that you create and set an appropriate error when returning nil.
    You can also choose to override -fileWrapperOfType:error:, -writeToURL:ofType:error:, or -writeToURL:ofType:forSaveOperation:originalContentsURL:error: instead.
    */
    if (outError) {
        *outError = [NSError errorWithDomain:NSOSStatusErrorDomain code:unimpErr userInfo:NULL];
    }
    return nil;
}

- (BOOL)readFromData:(NSData *)data ofType:(NSString *)typeName error:(NSError **)outError {
    /*
    Insert code here to read your document from the given data of the specified type. If outError != NULL, ensure that you create and set an appropriate error when returning NO.
    You can also choose to override -readFromFileWrapper:ofType:error: or -readFromURL:ofType:error: instead.
    */
    if (outError) {
        *outError = [NSError errorWithDomain:NSOSStatusErrorDomain code:unimpErr userInfo:NULL];
    }
    return YES;
}

@end
