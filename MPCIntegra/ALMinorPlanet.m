//
//  ALMinorPlanet.m
//  MPCIntegra
//
//  Created by Michele Bigi on 05/07/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import "ALMinorPlanet.h"
#include <math.h>

@implementation ALMinorPlanet

- (id)init
{
    self = [super init];
    if (self) {
        // Initialization code here.
    }
    
    return self;
}

- (id) initWithAstOrb:(NSString *) AstOrbPath : (long) numb
{
    self = [super init];
    if (self) {
        // Initialization code here.
        [ self ReadAstOrb:AstOrbPath :numb ]; //&mp_elem ];
        //getAstOrbParam( numb, &mp_elem );        
    }
    
    return self;
}

- (void)dealloc
{
    [super dealloc];
}

- (double) Epoc {
    return epoch;
}

- (void) Epoc:(double) eep {
    epoch = eep;
}

- (double) PerihelionTime{
    return perih_time;
}

- (double) Q {
    return q;
}

- (double) Eccentricity{
    return ecc;
}

- (double) Inclination {
    return incl;
}

- (double) ArgumentPerihelion{
    return arg_per;
}

- (double) LongitudeAscendingNode {
    return asc_node;
}

- (double) MeanAnomaly {
    return mean_anomaly;
}

- (double) Diameter {
    return diameter;
}

- (double) MinorToMajor {
    return minor_to_major;
}

- (void) MinorToMajor:(double) ttemp {
    minor_to_major = ttemp;
}

- (double) MajorAxis {
    return major_axis;
}

- (void) MajorAxis:(double) ttemp {
    major_axis = ttemp;
}

- (double) W0 {
    return w0;
}

- (void) W0:(double) ttemp {
    w0 = ttemp;
}

- ( double ) getDiameter {
    return diameter;
}

- ( double ) getDist {
    return dist;
}

- ( double ) getRA {
    return ra;
}

- ( double ) getDec {
    return dec;
}


- ( double ) getEarthDist {
    return earth_distance;
}

- ( double ) getCEU {
    return ceu;
}

- ( double ) getDCEU {
    return dceu;
}

- (double) Brightness {
    return brightness;
}

- (double) T0 {
    return t0;
}

- (double) Slope {
    return slope;
}

- (void) PerihelionTime:(double) t_val {
    perih_time = t_val;
}

- (void) Q:(double) t_val{
   q = t_val;
}

- (void) Eccentricity:(double) t_val {
   ecc = t_val;
}

- (void) Inclination:(double) t_val {
   incl = t_val;
}

- (void) ArgumentPerihelion:(double) t_val {
   arg_per = t_val;
}

- (void) LongitudeAscendingNode:(double) t_val {
   asc_node = t_val;
}

- (void) MeanAnomaly:(double) t_val {
   mean_anomaly = t_val;
}

- (void) Diameter:(double) t_val {
   diameter = t_val;
}

- (void) Brightness:(double) t_val {
   brightness = t_val;
}

- (void) T0:(double) t_val {
   t0 = t_val;
}

- (void) Slope:(double) t_val {
   slope = t_val;
}

- (void) PerihelionVector:(double *) t_val{
    perih_vec[0] = t_val[0];
    perih_vec[1] = t_val[1];
    perih_vec[2] = t_val[2];
}

- (double *) PerihelionVector {
    return perih_vec;
}

- (void) Sideways:(double *) t_val{
    perih_vec[0] = t_val[0];
    perih_vec[1] = t_val[1];
    perih_vec[2] = t_val[2];
}

- (double *) Sideways {
    return perih_vec;
}

- (void) AngularMomentum:(double) t_val {
    angular_momentum = t_val;
}

- (double) AngularMomentum {
    return angular_momentum;
}

- (void) ReadAstOrb:(NSString *) AstOrbPath : (long)asteroid_no  
//getAstOrbParam( long asteroid_no , ELEMENTS * class_elem ) 
{
    char tbuff[300];
    
    FILE *ifile = fopen( [AstOrbPath UTF8String], "rb");
    if( !ifile)
    {
        printf( "Couldn't find 'astorb.dat'\n");
        exit( -1);
    }
    fseek( ifile, (asteroid_no - 1) * ASTORB_RECORD_LEN, SEEK_SET);
    fgets( tbuff, sizeof( tbuff), ifile);
    fclose( ifile);
    
    [ self AL_ExtractAstOrbData:tbuff ];
}

long AL_extract_long( const char *ibuff)
{
    long rval = 0;
    while( *ibuff == ' ')      /* skip leading spaces */
        ibuff++;
    while( *ibuff != ' ')
    {
        if( *ibuff != '.')
            rval = rval * 10 + (long)( *ibuff - '0');
        ibuff++;
    }
    printf( "\n\n" );
    return( rval);
}

- (void) AL_SetupOrbitVectors {
    const double sin_incl = sin( incl), cos_incl = cos( incl);
    double * vec;
    double vec_len;
    double up[3];
    int i;
    
    minor_to_major = sqrt( fabs( 1. - ecc * ecc));
    lon_per = asc_node + atan2( sin( arg_per) * cos_incl, cos( arg_per));
    vec = perih_vec;
    
    vec[0] = cos( lon_per);
    vec[1] = sin( lon_per);
    vec[2] = (sin_incl / cos_incl) * sin( lon_per - asc_node);
    vec_len = sqrt( 1. + vec[2] * vec[2]);
    for( i = 0; i < 3; i++)
        vec[i] /= vec_len;
    /* 'up' is a vector perpendicular to the plane of the orbit */
    up[0] =  sin( asc_node) * sin_incl;
    up[1] = -cos( asc_node) * sin_incl;
    up[2] = cos_incl;
    
    sideways[0] = up[1] * vec[2] - up[2] * vec[1];
    sideways[1] = up[2] * vec[0] - up[0] * vec[2];
    sideways[2] = up[0] * vec[1] - up[1] * vec[0];
}

- (void) AL_DerivedQuantities:(double) gm
{
    if( ecc != 1.)    /* for non-parabolic orbits: */
    {
        major_axis = q / fabs(1. - ecc);
        t0 = major_axis * sqrt( major_axis / gm);
    }
    else
    {
        w0 = (3. / SQRT_2 ) / (q * sqrt( q / gm));
        major_axis = t0 = 0.;
    }
    [ self AL_SetupOrbitVectors ];

}

- (void) AL_DoRemainingSetup 
{
    mean_anomaly *= PI / 180.0 ;
    arg_per      *= PI / 180.0 ;
    asc_node     *= PI / 180.0 ;
    incl         *= PI / 180.0 ;
    q = major_axis * (1. - ecc);
    [ self AL_DerivedQuantities:SOLAR_GM ];
    angular_momentum = sqrt( SOLAR_GM * q * (1. + ecc));
    perih_time = epoch - mean_anomaly * t0;
    is_asteroid = 1;
    central_obj = 0;
}

static void get_jul_greg_year_data( const long year, long *days, char *month_data )
{
    static const char months[13] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 0 };
    
    *days = year * 365L + year / 4L;
    *days += -year / 100L + year / 400L;
          
    memcpy( month_data, months, 13);
    if( !(year % 4))
        if( (year % 100L) || !(year % 400L) )
        {
            month_data[1] = 29;
            (*days)--;
        }
    *days += JUL_GREG_CALENDAR_EPOCH + 1;
}


static int get_calendar_data( const long year, long *days, char *month_data )
{
    int i, rval = 0;
    
    memset( month_data, 0, 13);
    get_jul_greg_year_data( year, days, month_data );
    days[1] = days[0];
    for( i = 0; i < 13; i++)
        days[1] += month_data[i];
    
    return( rval);
}

long     dmy_to_day( const int day, const int month, const long year )
{
    char mdata[13];
    long jd;
    long year_ends[2];
 //   int calendar_to_use = calendar;
    int i, rval;
    
    rval = get_calendar_data( year, year_ends, mdata );
    if( !rval)
    {
        jd = year_ends[0];
        for( i = 0; i < month - 1; i++)
            jd += mdata[i];
        jd += (long)day - 1;
    }
    else
        jd = 0;
    return( jd);
}

double Cal2JED(int m, int d, int y) {
    double time, datnum, a, b, term1;
    double jed;
    
    /* Convert to decimal hours */
    time = 0;
    
    if ((m == 1) || (m == 2)) {
        y = y - 1;
        m = m + 12;
    }
    
    /*
     Test to see if date is in Gregorian calendar.
     Converts date to a number of the form YYYY.MMDD.
     */
    datnum = (double) y + 0.01 * (double) m + 0.0001 * (double) d;
    
    if (datnum >= 1582.1015) {
        /* Gregorian calendar */
        a = (int)(0.01 * (double) y);
        b = 2.0 - a + (int)(0.25 * (double) a);
    } else {
        /* Julian calendar */
        a = 0.0;
        b = 0.0;
    }
    
    if (y < 0) {
        /* Handles negative years */
        term1 = (int)(365.25 * (double) y - 0.75); /* change here */
    } else {
        term1 = (int)(365.25 * (double) y);
    }
    
    jed = term1 + (int)(30.6001 * ((double) m + 1.0)) + (double) d + time / 24.0 + 1720994.5 + (double) b;
    
    return jed;
    
}


- (void) AL_ExtractAstOrbData:(char *) buff
{
    long epoch_jd;
   
    if( strlen( buff) > 267
       && sscanf( buff + 41, "%lf %lf", &abs_mag, &slope_param) == 2)
    {
        brightness   = (double)AL_extract_long( buff + 42 ) / 1.e+4;
        slope        = (double)AL_extract_long( buff + 50 ) / 1.e+3;
        diameter     = (double)AL_extract_long( buff + 59 ) / 1.e+4;
        epoch_jd = atoi( buff + 106);
        mean_anomaly = (double)AL_extract_long( buff + 115) / 1.e+6;
        arg_per      = (double)AL_extract_long( buff + 126) / 1.e+6;
        asc_node     = (double)AL_extract_long( buff + 137) / 1.e+6;
        incl         = (double)AL_extract_long( buff + 148) / 1.e+6;
        ecc          = (double)AL_extract_long( buff + 160) / 1.e+8;
        ceu          = (double)AL_extract_long( buff + 192) / 1.e+6;
        dceu         = (double)AL_extract_long( buff + 200) / 1.e+7;
        major_axis = atof( buff + 169);
        NSLog( @"Epoc 1: %lf\n" , (double) epoch_jd );
        int yy, mm, gg;
        gg = (int)(epoch_jd % 100L);
        mm = (int)((epoch_jd / 100L) % 100L );
        yy = (int)(epoch_jd / 10000L );
                   
        epoch_jd = Cal2JED( mm, gg, yy ); //dmy_to_day( (int)(epoch_jd % 100L), (int)(epoch_jd / 100L) % 100L, epoch_jd / 10000L );
        NSLog( @"Epoc 2: %lf\n" , (double)epoch_jd );
        epoch = (double)epoch_jd - .5;
        [ self AL_DoRemainingSetup ];
        //do_remaining_element_setup( &telem);
        if ( !diameter ) 
            if( brightness ) 
                diameter = pow ( 10.0, 3.52 - 0.2 * brightness );
            
    }
}

@end
