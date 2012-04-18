//
//  Math3D.m
//  iOccult
//
//  Created by Michele Bigi on 19/04/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import "costants.h"
#import "Math3D.h"
#import "Vec3D.h"

@implementation Math3D

- (id)init
{
    self = [super init];
    if (self) {
        // Initialization code here.
        int i,j;
        for(i=0;i<3;i++)
            for( j=0;j<3;j++)
                m_Mat[i][j] = 0.0;
        
    }
    
    return self;
}

- (id) initWithMatrix :(Math3D * ) matr {
    self = [super init];
    if (self) {
        // Initialization code here.
        int i,j;
        for(i=0;i<3;i++)
            for( j=0;j<3;j++)
                m_Mat[i][j] = [matr getValue:i :j ];
        
    }
    
    return self;
}

- (id)initWithIdentity
{
    self = [super init];
    if (self) {
        // Initialization code here.
        int i,j;
        for(i=0;i<3;i++)
            for( j=0;j<3;j++)
                if( i != j )
                    m_Mat[i][j] = 0.0;
                else
                    m_Mat[i][j] = 1.0;
        
    }
    
    return self;
}


- (id) init:(Vec3D *) v1: (Vec3D*) v2: (Vec3D*) v3{
    self = [super init];
    if (self) {
        // Initialization code here.
        int i;
        
        for (i = 0; i < 3; i++) {
            m_Mat[i][0]=[v1 getVecComponent:i:false ];
            m_Mat[i][1]=[v2 getVecComponent:i:false ];
            m_Mat[i][2]=[v3 getVecComponent:i:false ];
        }
    }
    
    return self;
}

- (id) initWithRotationX :(double) angle {
    self = [super init];
    if (self) {
        // Initialization code here.
        int i,j;
        for(i=0;i<3;i++)
            for( j=0;j<3;j++)
                if( i != j )
                    m_Mat[i][j] = 0.0;
                else
                    m_Mat[i][j] = 1.0;
        
    }
    double S, C;
    S = sin (angle);
    C = cos (angle);
    
    [ self setValue: 1: 1 :  C ];
    [ self setValue: 1: 2 :  S ];
    [ self setValue: 2: 1 : -S ];
    [ self setValue: 2: 2 :  C ];
    
    return self;
}

- (id) initWithRotationY :(double) angle {
    self = [super init];
    if (self) {
        // Initialization code here.
        int i,j;
        for(i=0;i<3;i++)
            for( j=0;j<3;j++)
                if( i != j )
                    m_Mat[i][j] = 0.0;
                else
                    m_Mat[i][j] = 1.0;
        
    }
    double S, C;
    S = sin (angle);
    C = cos (angle);
    
    [ self setValue: 0: 0 :  C ];
    [ self setValue: 0: 2 : -S ];
    [ self setValue: 2: 0 :  S ];
    [ self setValue: 2: 2 :  C ];
    
    return self;
}

- (id) initWithRotationZ :(double) angle {
    self = [super init];
    if (self) {
        // Initialization code here.
        int i,j;
        for(i=0;i<3;i++)
            for( j=0;j<3;j++)
                if( i != j )
                    m_Mat[i][j] = 0.0;
                else
                    m_Mat[i][j] = 1.0;
        
    }
    double S, C;
    S = sin (angle);
    C = cos (angle);
    
    [ self setValue: 0: 0 :  C ];
    [ self setValue: 0: 1 :  S ];
    [ self setValue: 1: 0 : -S ];
    [ self setValue: 1: 1 :  C ];
    
    return self;
}

- (id) initWithDoubleMatrix:(double [3][3]) val {
    self = [super init];
    if (self) {
        // Initialization code here.
        int i,j;
        for(i=0;i<3;i++)
            for( j=0;j<3;j++)
                m_Mat[i][j] = val[i][j];
        
    }
    
    return self;
}

- (id) initEq2Ecl :(double) angle {
    double eps = ( 23.43929111-(46.8150+(0.00059-0.001813*angle)*angle)*angle/3600)*PI/180.0;
    self = [super init];
    if (self) {
        //        [ self Eq2Ecl:angle];
        
        int i,j;
        for(i=0;i<3;i++)
            for( j=0;j<3;j++)
                if( i != j )
                    m_Mat[i][j] = 0.0;
                else
                    m_Mat[i][j] = 1.0;
    }
    double S, C;
    S = sin (eps);
    C = cos (eps);
    
    [ self setValue: 0: 0 :  C ];
    [ self setValue: 0: 1 :  S ];
    [ self setValue: 1: 0 : -S ];
    [ self setValue: 1: 1 :  C ];
    
    return self;
}


- (void)dealloc
{
    [super dealloc];
}

- (void) setValue:(int) r: (int) c : (double) val {
    m_Mat[r][c] = val;
    
}

- (double) getValue: (int) r: (int) c {
    return m_Mat[r][c];
}


- ( Vec3D * ) getRow:(int) index {
    
    Vec3D * res = [[ Vec3D alloc ] init ];
    [res set: m_Mat[index][0]:m_Mat[index][1]:m_Mat[index][2]:false ];
    
    return res; 
    
}

- ( void ) setRow:(int) index : (Vec3D *) row {
    
    //Vec3D * res = [[ Vec3D alloc ] init ];
    //[res set: m_Mat[index][0]:m_Mat[index][1]:m_Mat[index][2]:false ];
    
    //return res; 
    m_Mat[index][0] = [ row getElement:0 ];
    m_Mat[index][1] = [ row getElement:1 ];
    m_Mat[index][2] = [ row getElement:2 ];
}

- ( Vec3D * ) getCol:(int) index {
    
    Vec3D * res = [[ Vec3D alloc ] init ];
    [res set: m_Mat[0][index]:m_Mat[1][index]:m_Mat[2][index]:false ];
    return res; 
}

- (double) determinant {
    double det;
    
    det = m_Mat[0][0] * m_Mat[1][1] * m_Mat[2][2] + m_Mat[1][0] * m_Mat[2][1] * m_Mat[0][2] + m_Mat[2][2] * m_Mat[1][0] * m_Mat[0][1]; 
    return det;
}

- (void) ScalarProduct:(double) val {
    int i,j;
    
    for(i=0;i<3;i++)
        for(j=0;j<3;j++)
            m_Mat[i][j] *= val;
}

- (Math3D *) MatrixProduct:(Math3D *) b {
    /*
    Math3D * res = [[Math3D alloc] init ];
    int i,j;
    for(i=0;i<3;i++)
        for(j=0;j<3;j++) {
            [res setValue:i :j :[ self getValue:i :j ] * [b getValue:j:i ]];
        }*/
    return [self product:b];
}

- (Math3D *) product:(Math3D *) b {
    Math3D * res = [[Math3D alloc] init ];
    int i,j,k;
    double scalp;
    
    for(i=0;i<3;i++){
        for(j=0;j<3;j++) {
            scalp = 0.0;
            for(k=0;k<3;k++)
                scalp +=  [self getValue:i :k ] * [b getValue:k :j ];
            [res setValue:i :j :scalp ];
        }
    }
    
    return res;
}


- (Vec3D *) VectorMatrixProduct: (Vec3D *) b {
    Vec3D * res = [[Vec3D alloc ] init ];
    int i,j;
    
    for(i=0;i<3;i++)
        for(j=0;j<3;j++) {
            double tr, tt;
            tr = [res getElement:i];
            tt =  [self getValue:i :j] * [ b getElement:j ];
            [ res setElement:i : tr + tt ];
            //[ res setElement:i : [res getElement:i] + [self getValue:i :j] * [ b getElement:j ]];
        }
    return res;
}

- (Math3D *) traspose{
    Math3D * res = [[Math3D alloc] init ];
    
    int i,j;
    for(i=0;i<3;i++)
        for(j=0;j<3;j++)
            [ res setValue:j :i : m_Mat[i][j] ];
    
    return res;
}

- (Math3D *) MatrixSum:(Math3D *) b {
    Math3D * res = [[Math3D alloc] init ];
    int i,j;
    for(i=0;i<3;i++)
        for(j=0;j<3;j++) {
            [res setValue:i :j :m_Mat[i][j] + [b getValue:i:j ]];
        }
    return res;
}

- (int) find_zero:(int) k {
    int k1;
    
    k1 = k;
    if( (k == 2) && ( m_Mat[k][k] == 0.0 ) )
        return k+1;
    if( (k == 2) && ( m_Mat[k][k] != 0.0 ) )
        return k;
    
    do;
    while ( ( m_Mat[k++][k1] == 0.0 ) && (k<3));
    if (m_Mat[k-1][k1] == 0.0 )
        return k;
    return k-1;
}

- (void) pivot:(int) k:(int)k1 {
    Vec3D * rigan = [self getRow:k ];
    Vec3D * rigam = [self getRow:k1];
    
    [self setRow:k  :rigam ];
    [self setRow:k1 :rigan ];
}

- (void) toString{
    int i;
    for(i=0;i<3;i++) {
        Vec3D * v0 = [ self getRow:i ];
        [v0 toString ];
    }
}

- (Math3D *) Eq2Ecl:(double) rad {
    double eps = ( 23.43929111-(46.8150+(0.00059-0.001813*rad)*rad)*rad/3600)*PI/180.0;
    Math3D * eq2ecl = [[ Math3D alloc] initWithRotationX: eps];
    
    return eq2ecl;
}

@end

