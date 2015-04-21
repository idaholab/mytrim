#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <cmath>
#include <stdlib.h>

namespace MyTRIM_NS {

inline void v_cross( const double *a1, const double *a2, double *b )
{
    for( int i = 0; i < 3; i++ ) b[i] = a1[(i+1)%3] * a2[(i+2)%3] - a1[(i+2)%3] * a2[(i+1)%3];
}
inline void v_scale( double *a1, double b ) // in=place scale
{
    for( int i = 0; i < 3; i++ ) a1[i] = a1[i] * b;
}
inline double v_dot( const double *a1, const double *a2 )
{
    double b = 0.0;
      for( int i = 0; i < 3; i++ ) b += a1[i] * a2[i];
        return b;
}
inline void v_norm( double *a1, double b = 1.0 ) // in-place normalize to b ( = 1.0 default )
{
    v_scale( a1, b / sqrtf( v_dot( a1, a1) ) );
}

inline double sqr( double a ) { return a*a; }
inline double cub( double a ) { return a*a*a; }

// random numbers
const double drm = double(RAND_MAX)+1.0;
inline double dr250() { return double(rand())/drm; };
inline void r250_init( int s ) { srand(s); }

}

#endif
