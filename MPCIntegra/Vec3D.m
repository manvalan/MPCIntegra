//
//  Vec3D.m
//  iOccult
//
//  Created by Michele Bigi on 18/04/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import "Vec3D.h"
#import "costants.h"

@implementation Vec3D

- (id)init
{
    self = [super init];
    if (self) {
        [ self setElement:0 :0.0]; 
        [ self setElement:1 :0.0]; 
        [ self setElement:2 :0.0]; 
        // Initialization code here.
    }
    
    return self;
}

- (id)init:( double ) e1: ( double ) e2 : (double ) e3 
{
    self = [super init];
    if (self) {
        // Initialization code here.
        m_Vec[0] = e1;
        m_Vec[1] = e2;
        m_Vec[2] = e3;
        m_bPolarValid = FALSE;
        
        [self CalcPolarAngles ]; 
    }
    
    return self;
}

- (id)initPolar:( double ) r: ( double ) theta : ( double )phi
{
    self = [super init];
    if (self) {
        // Initialization code here.
        const double cosE = cos(theta);
        
        m_Vec[0] = r * cos(phi) * cosE;
        m_Vec[1] = r * sin(phi) * cosE;
        m_Vec[2] = r * sin( theta );
        m_r = r ;
        m_phi = phi;
        m_theta = theta;
        m_bPolarValid = TRUE;
        
       // [self CalcPolarAngles ]; 
    }
    
    return self;
}

- (id)initWithVec3d:(Vec3D *) e
{
    self = [super init];
    if (self) {
        // Initialization code here.
        m_Vec[0] = [ e getElement:0 ];
        m_Vec[1] = [ e getElement:1 ];
        m_Vec[2] = [ e getElement:2 ];
        m_bPolarValid = FALSE;
        
        [self CalcPolarAngles ]; 
    }
    
    return self;
}


- (id)init:( double ) e1: ( double ) e2 : (double ) e3 : (bool) polar
{
    self = [super init];
    if (self) {
        // Initialization code here.
        if( !polar ) {
        m_Vec[0] = e1;
        m_Vec[1] = e2;
        m_Vec[2] = e3;
        m_bPolarValid = FALSE;
        
        [self CalcPolarAngles ]; 
        }
        else {
            m_r = e1;
            m_phi = e2;
            m_theta = e3;
            m_bPolarValid = true;
        }
    }
    
    return self;
}

- (void) set:( double ) e1: ( double ) e2 : (double ) e3 : (bool) polar{
    if( !polar ) {
        m_Vec[0] = e1;
        m_Vec[1] = e2;
        m_Vec[2] = e3;
        m_bPolarValid = FALSE;
        
        [self CalcPolarAngles ]; 
    }
    else {
        m_r = e1;
        m_phi = e2;
        m_theta = e3;
        m_bPolarValid = true;
    }
}

- (void)Copy: (Vec3D * ) Right
{
    m_Vec[0] = [Right getVecComponent:0 :NO ];
    m_Vec[1] = [Right getVecComponent:1 :NO ];
    m_Vec[2] = [Right getVecComponent:2 :NO ];
                
    [self CalcPolarAngles ];
}
- (void)dealloc
{
    [super dealloc];
}

- (void) CalcPolarAngles {
    double rhoSqr; 
    rhoSqr = m_Vec[0] * m_Vec[0] + m_Vec[1] * m_Vec[1];
    m_r = sqrt( rhoSqr + m_Vec[2] * m_Vec[2] );
    
    if( ( m_Vec[0] == 0.0 ) && (m_Vec[1] == 0.0 ) )
        m_phi = 0.0;
    else
        m_phi = atan2( m_Vec[1] , m_Vec[0] );
    
    if( m_phi < 0.0 )
        m_phi += 2.0 * PI;
    
    double rho = sqrt( rhoSqr );
    if( ( m_Vec[2] == 0.0 ) && ( rhoSqr == 0.0 )  )
        m_theta = 0.0;
    else
        m_theta = atan2( m_Vec[2] , rho );
    
}

- (void) setElement:(int) i: (double) val {
    m_Vec[i] = val;
}

- (double) getElement:(int) i {
    return m_Vec[i];
}

- ( Vec3D * )  Sum:(Vec3D * ) second {
    double m0, m1, m2;
    
    m0 = m_Vec[0] + [second getVecComponent:0 : FALSE ];
    m1 = m_Vec[1] + [second getVecComponent:1 : FALSE ];
    m2 = m_Vec[2] + [second getVecComponent:2 : FALSE ];
    
    Vec3D * res = [[ Vec3D alloc ] init:m0:m2:m2 ];
    
    return res; 
}

- ( Vec3D * )  Dif:(Vec3D * ) second {
    double m0, m1, m2;
    
    m0 = m_Vec[0] - [second getVecComponent:0 : FALSE ];
    m1 = m_Vec[1] - [second getVecComponent:1 : FALSE ];
    m2 = m_Vec[2] - [second getVecComponent:2 : FALSE ];
    
    Vec3D * res = [[ Vec3D alloc ] init:m0:m2:m2 ];
    
    return res; 
}

- ( Vec3D * ) Prod:(double) val {
    double m0, m1, m2;
    
    m0 = m_Vec[0] * val ;
    m1 = m_Vec[1] * val ;
    m2 = m_Vec[2] * val ;
    
    Vec3D * res = [[ Vec3D alloc ] init:m0:m1:m2 ];
    
    return res; 
    
    
}

- ( double ) ScalarProduct:(Vec3D *) b {
    double m0, m1, m2;
    
    m0 = m_Vec[0] * [b getVecComponent:0 :FALSE];
    m1 = m_Vec[1] * [b getVecComponent:1 :FALSE];
    m2 = m_Vec[2] * [b getVecComponent:2 :FALSE];
    
    return m0+m1+m2; 
}

- (double) Dot:( Vec3D *) b {
    return [self ScalarProduct:b];
}

- ( Vec3D * ) VectProduct:(Vec3D *) b {
    double m0, m1, m2;
    
    m0 = [ self getVecComponent:1:false] * [ b getVecComponent:2 :FALSE] - [ self getVecComponent:2:false] * [b getVecComponent:1 :FALSE];
    m1 = [ self getVecComponent:2:false] * [ b getVecComponent:0 :FALSE] - [ self getVecComponent:0:false] * [b getVecComponent:2 :FALSE];
    m2 = [ self getVecComponent:0:false] * [ b getVecComponent:1 :FALSE] - [ self getVecComponent:1:false] * [b getVecComponent:0 :FALSE];
    
    Vec3D * res = [[ Vec3D alloc ] init:m0:m1:m2 ];
    
    return res;
}

- (double) getVecComponent:(int) comp {
    double rv;
    
    if( !m_bPolarValid )
        rv = m_Vec[ comp ];
    else
        switch (comp) {
            case 0:
                rv = m_r;
                break;
            
            case 1:
                rv = m_phi;
                
            case 2:
                rv = m_theta;
                
            default:
                rv = 0;
                break;
        }
    return rv;
}

- (double) getVecComponent:(int) comp : (bool) polar {
    double rv;
    
    if( !polar )
        rv = m_Vec[ comp ];
    else
        switch (comp) {
            case 0:
                rv = m_r;
                break;
                
            case 1:
                rv = m_phi;
                
            case 2:
                rv = m_theta;
                
            default:
                rv = 0;
                break;
        }
    return rv;
}

- (double) Norm {
    return sqrt( [self getElement:0 ]*[self getElement:0 ] + [self getElement:1 ]*[self getElement:1 ]+ [self getElement:2 ]*[self getElement:2 ]);
}
- (double) PosAngle:(Vec3D *) right
{
    double angle, cosA;
    
    cosA = ( [ self Dot:right ] / ( [self Norm] * [right Norm] ) );
    angle = acos( cosA );
    return angle;
}

- (void) Normalize {
    double Norma;
    
    Norma = [ self Norm ];
    m_Vec[0] /= Norma;
    m_Vec[1] /= Norma;
    m_Vec[2] /= Norma;
}

- (void) toString {
    NSLog( @"%lf   %lf  %lf\n" , [self getElement:0 ], [self getElement:1 ], [self getElement:2 ] );
}

@end
