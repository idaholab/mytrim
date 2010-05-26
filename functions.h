#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <cmath>

inline void v_cross( const float *a1, const float *a2, float *b )
{
    for( int i = 0; i < 3; i++ ) b[i] = a1[(i+1)%3] * a2[(i+2)%3] - a1[(i+2)%3] * a2[(i+1)%3];
}
inline void v_scale( float *a1, float b ) // in=place scale
{
    for( int i = 0; i < 3; i++ ) a1[i] = a1[i] * b;
}
inline float v_dot( const float *a1, const float *a2 )
{
    float b = 0.0;
      for( int i = 0; i < 3; i++ ) b += a1[i] * a2[i];
        return b;
}
inline void v_norm( float *a1, float b = 1.0 ) // in-place normalize to b ( = 1.0 default )
{
    v_scale( a1, b / sqrtf( v_dot( a1, a1) ) );
}

inline float sqr( float a ) { return a*a; }
inline float cub( float a ) { return a*a*a; }

#endif

