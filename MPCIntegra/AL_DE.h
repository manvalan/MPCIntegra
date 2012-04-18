//
//  AL_DE.h
//  MPCIntegra
//
//  Created by Michele Bigi on 11/07/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "ALAstOrbIntegFunction.h"

#define max(a,b)    ((a>b)?a:b) 
#define min(a,b)    ((a<b)?a:b)

enum eDE_STATE {
    DE_INIT                  = 1,    // Starting step
    DE_DONE                  = 2,    // Successful integration step
    DE_ACCURACY_NOT_ACHIEVED = 3,    // Too stringent accuracy requirements
    DE_TOO_MANY_STEPS        = 4,    // Too many steps required
    DE_STIFF                 = 5,    // Suspect of stiff differential equation
    DE_INVALID_PARAMS        = 6     // Invalid input
};

typedef enum eDE_STATE  DE_STATE;

@interface AL_DE : NSObject {
@private
    ALAstOrbIntegFunction * IntegFunc;
    
    double   *yy,*wt,*p,*yp,*ypout;
    double   **phi;
    double   alpha[13],beta[13],v[13],w[13],psi[13];
    double   sig[14],g[14];
    double   x,h,hold,told,delsgn;
    int      ns,k,kold;
    bool     OldPermit, phase1,start,nornd;
    
    int maxnum;  // Maximum number of steps to take
    
    double umach;
    double twou;   
    double fouru;   
}

-(id)initDE:( ALAstOrbIntegFunction * ) afunc;
-(void) AL_DE_Intrp:(double) x:(double *) y:(double) xout:(double *)yout:(double *)ypout;
-(void) AL_DE_Step:(double *) x:(double *)y:(double *)eps:(bool *) crash;
-(void) AL_DE_Integ:(double *) y:(double *) t:(double) tout:(double *) relerr:(double *) abserr:(DE_STATE *)State:(bool) PermitTOUT;

@end

