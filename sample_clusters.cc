#include <math.h>
#include "malloc.h"

#include "sample_clusters.h"
#include <r250.h>

#include "functions.h"

sampleClusters::sampleClusters( double x, double y, double z )  : sampleBase( x, y, z)
{ 
  sh = 0;
  cl = 0; cn = 0; cnm = 0;
  for( int i = 0; i < 4; i++ ) c[i] = 0;
}

// look if we are within dr of a cluster
// dr == 0.0 means looking if we are inside the cluster
materialBase* sampleClusters::lookupMaterial( double* pos ) 
{
  int l = lookupCluster( pos, 0.0 );

  if( l == -2 ) return 0;
  if( l == -1 ) return material[0];

  material[1]->tag = l;
  return material[1];
}

// look if we are within dr of a cluster
// dr == 0.0 means looking if we are inside the cluster
int sampleClusters::lookupCluster( double* pos, double dr ) 
{
  double dif[3], r2;
  int k[3], k1[3], k2[3], j[3], l, ks;

  // spatial hash center and span
  // if pos>w || pos<0 drop out early when bc[] = CUT or return material[0] 
  for( int i = 0; i<3; i++ ) 
  {
    k[i] = floor( ( pos[i] * kn[i] ) / w[i] );
    if( pos[i] < 0.0 || pos[i] >= w[i] )
    {
      switch( bc[i] )
      {
        case CUT : return -2;
        case INF : return -1;
        case PBC : k[i] = k[i] % kn[i];
                   if( k[i] < 0 ) k[i] += kn[i];
      }
    }

    // look ks hash boxes around current box
    ks = int( ( cmr + dr ) / kd[i] ) + 1;
    k1[i] = k[i] - ks ;
    k2[i] = k[i] + ks;
    //printf("k12\t%d\t%d\t%d\n", k[i], k1[i], k2[i] );

    if( k1[i] < 0      && bc[i] != PBC ) k1[i] = 0;
    if( k2[i] >= kn[i] && bc[i] != PBC ) k2[i] = kn[i] - 1;
  }

  for( j[0] = k1[0]; j[0] <= k2[0]; j[0]++ )
    for( j[1] = k1[1]; j[1] <= k2[1]; j[1]++ )
      for( j[2] = k1[2]; j[2] <= k2[2]; j[2]++ )
      {
        for( int i = 0; i<3; i++ )
        {
          k[i] = j[i] % kn[i];
          if( k[i] < 0 ) k[i] += kn[i];
        }
        l = k[0] + kn[0] * ( k[1] + kn[1] * k[2] );

        if( sh[l] >= 0 )
        {
          l = sh[l];
          //printf("entering while loop\n" );
          while( l >= 0 ) 
          { 
            for( int i = 0; i<3; i++ )
            {
              dif[i] = pos[i] - c[i][l];
              if( bc[i] == PBC ) dif[i] -= round( dif[i] / w[i] ) *w[i];
            }
            r2 = v_dot( dif, dif );
            //printf(" trying cluster %d, dif=(%f,%f,%f), r2=%f\n", l,dif[0],dif[1],dif[2],r2);
            if( r2 < sqr( c[3][l] + dr ) )
            {
              return l;
            }
            l = cl[l];
          }
        }
      }

  return -1;
}

void sampleClusters::initSpatialhash( int x, int y, int z )
{
  kn[0] = x; kn[1] = y; kn[2] = z;
  sh = (int*)malloc( sizeof(int) * x*y*z );
  clearSpatialHash();
  
  // calculate half the spatial diagonal of a hash block
  sd = 0.0;
  for( int i = 0; i < 3; i++ ) 
  {
    kd[i] = w[i] / double( kn[i] );
    sd += kd[i];
  }
  sd = 0.5 * sqrtf( sd );
  cmr = 0.0;
}

void sampleClusters::clearSpatialHash()
{
  for( int i = 0; i < kn[0]*kn[1]*kn[2]; i++ ) sh[i]=-1;
}

void sampleClusters::reallocClusters( int n )
{
  if( n > cnm )
  {
    cl = (int*)realloc( cl, sizeof(int) * n );
    for( int i = 0; i < 4; i++ )
    {
      c[i] = (double*)realloc( c[i], sizeof(double) * n );
    }
    for( int j = cnm; j < n; j++ ) cl[j] = -1;
    cnm = n;
  }
}

void sampleClusters::clearClusters()
{
  cn = 0;
  clearSpatialHash();
  for( int i = 0; i < cnm; i++ ) cl[i] = -1;
}

void sampleClusters::addCluster( double x, double y, double z, double r )
{
  int k[3], i, l;
  double cb[4];

  if( cn >= cnm ) reallocClusters( cnm + cnm/10 + 10 ); // get 10% more slots plus 10
  c[0][cn] = x; c[1][cn] = y; c[2][cn] = z; c[3][cn] = r;

  for( int i = 0; i < 3; i++ ) 
  {
    k[i] = int( floor( ( c[i][cn] * kn[i] ) / w[i] ) ) % kn[i];
    k[i] = (k[i] < 0 ) ? k[i] + kn[i] : k[i];
  }

  l = k[0] + kn[0] * ( k[1] + kn[1] * k[2] );
  if( sh[l] < 0 ) sh[l] = cn;
  else
  {
    l = sh[l];
    while( cl[l] >= 0 ) { l = cl[l]; };
    cl[l] = cn;
  }
  cl[cn] = -1;

  if( r > cmr ) cmr = r;

  cn++;
}

// add non-overlapping clusters with a minimum surface-surface separation of dr
void sampleClusters::addRandomClusters( int n, double r, double dr )
{
  double npos[3];

  reallocClusters( n + n/10 ); //allocate 10% more for ghost bubbles (get more later if needed)

  for( int i = 0; i < n; i++ )
  {
    while( true )
    {
      for( int j = 0; j < 3; j++ ) npos[j] = dr250() * w[j];

      if( lookupCluster( npos, dr + r ) == -1 )
      {
        addCluster( npos[0], npos[1], npos[2], r );
        break;
      }
      //else fprintf( stderr, "rejected, too close\n" );
    }

    if( i+1 % 10000 == 0 ) fprintf( stderr, " %d (%d) clusters added\n", i+1, cn );
  }
}

