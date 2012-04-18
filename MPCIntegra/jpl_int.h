#include "stdintvc.h"

/* A JPL binary ephemeris header contains five doubles and */
            /* (up to) 41 int32_t integers,  so:                          */
#define JPL_HEADER_SIZE (5 * sizeof( double) + 41 * sizeof( int32_t))

#pragma pack(1)

#define TRUE 1
#define FALSE 0

struct jpl_eph_data {
   double ephem_start, ephem_end, ephem_step;
   int32_t ncon;
   double au;
   double emrat;
   int32_t ipt[13][3];
   int32_t ephemeris_version;
               /* This is the end of the file header.  Following are */
               /* items computed within my code.                     */
   int32_t kernel_size, recsize, ncoeff;
   int32_t swap_bytes;
   int32_t curr_cache_loc;
   double pvsun[6];
   double *cache;
   void *iinfo;
   FILE *ifile;
   };
#pragma pack()

struct interpolation_info
   {
   double pc[18],vc[18], twot;
   int np, nv;
   };

