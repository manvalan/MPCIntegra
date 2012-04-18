//
//  Vec3D.h
//  iOccult
//
//  Created by Michele Bigi on 18/04/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "costants.h"
//#define PI  3.14159265358979324;


@interface Vec3D : NSObject {
    
@private
    double m_Vec[3];
    double m_phi;
    double m_theta;
    double m_r;
    bool   m_bPolarValid;
}

- (id)init;
- (id)init:( double ) e1: ( double ) e2 : (double ) e3 ;
- (id)init:( double ) e1: ( double ) e2 : (double ) e3 : (bool) polar;
- (id)initPolar:( double ) r: ( double ) theta : ( double )phi;
- (id)initWithVec3d:(Vec3D *) e;

- (void) set:( double ) e1: ( double ) e2 : (double ) e3 : (bool) polar;

- (void) CalcPolarAngles;
- (double) getVecComponent:(int) comp;
- (double) getVecComponent:(int) comp : (bool) polar ;

- ( Vec3D * ) Sum:(Vec3D * ) second;
- ( Vec3D * ) Dif:(Vec3D * ) second;
- ( Vec3D * ) Prod:(double) val;
- ( double) ScalarProduct:(Vec3D *) b;
- ( Vec3D * ) VectProduct:(Vec3D *) b;
- (double) Dot:( Vec3D *) b;
- (double) Norm;
- (void) Normalize;
- (double) PosAngle:(Vec3D *) right;
- (void)Copy: (Vec3D * ) Right;

- (void) setElement:(int) i: (double) val;
- (double) getElement:(int) i;

- (void) toString;

@end
