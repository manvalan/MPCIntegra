//
//  Math3D.h
//  iOccult
//
//  Created by Michele Bigi on 19/04/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "costants.h"
#import "Vec3D.h"


@interface Math3D : NSObject {
@private
    double m_Mat[3][3];

}

- (id) init;
- (id) init:(Vec3D *) v1: (Vec3D*) v2: (Vec3D*) v3;
- (id) initWithIdentity;
- (id) initWithRotationX :(double) angle; 
- (id) initWithRotationY :(double) angle;
- (id) initWithRotationZ :(double) angle;
- (id) initEq2Ecl :(double) angle;
- (id) initWithMatrix :(Math3D * ) matr;
- (id) initWithDoubleMatrix:(double [3][3]) val;

- (void) setValue:(int) r: (int) c : (double) val;
- (double) getValue: (int) r: (int) c;
- ( Vec3D * ) getRow:(int) index ;
- ( Vec3D * ) getCol:(int) index ;
- (double) determinant ;
- (void) ScalarProduct:(double) val;
- (Math3D *) MatrixProduct:(Math3D *) b;
- (Math3D *) product:(Math3D *) b;
- (Vec3D *) VectorMatrixProduct: (Vec3D *) b;
- (Math3D *) traspose;
- (Math3D *) MatrixSum:(Math3D *) b;
- (Math3D *) Eq2Ecl:(double) rad;
- (void) toString;

//- (Math3D *) diagonalize;

@end
