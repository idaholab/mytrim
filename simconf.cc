#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "config.h"
#include "simconf.h"

simconfType *simconf;

simconfType::simconfType( double _alfa )
{
  ed = 25.0; // displacement energy
  alfa = _alfa; // angle of incidence (degrees)
  alpha = alfa * M_PI / 180.0;
  tmin = 1.0; //max impact parameter set by min. transferred energy
  //tmin = 5.0; //max impact parameter set by min. transferred energy
  tau = 0.0;
  da = 3.0; // angular grid for transmitted ions
  cw = 0.001; // channel width 1% of layer thickness

  // output full trajectories
  fullTraj = false;

  // set global ion id to zero (will be incremented for each new projectile)
  id = 0;

  // initialize global statistics
  vacancies_created = 0;
  KP_vacancies = 0.0;

  // read data tables
  read_snuc();
  read_scoef();
}

void simconfType::read_snuc()
{
  char buf[200], fname[500];

  //snprintf( fname, 500, "%s/SNUC03.dat", getenv("HOME") );
  snprintf( fname, 500, "%s/SNUC03.dat", DATA_DIR );
  FILE *sf = fopen( fname, "rt" );
  if( sf == 0 )
  {
    fprintf( stderr, "unable to open %s\n", fname );
    exit(1);
  }
  for( int i = 0; i < 92; i++ )
    for( int j = i; j < 92; j++ )
    {
      fscanf( sf, "%*d %*d %lf %lf %lf %f\n",
        &snuc[j][i][0], &snuc[j][i][1], &snuc[j][i][2], &snuc[j][i][3] );
      for( int n = 0; n < 4; n++ )
        snuc[i][j][n] = snuc[j][i][n];
    }
  fclose( sf );
}

void simconfType::read_scoef()
{
  char buf[2001], fname[500];
  FILE *sf;

  //snprintf( fname, 500, "%s/SCOEF.95A", getenv("HOME") );
  snprintf( fname, 500, "%s/SCOEF.95A", DATA_DIR );
  sf = fopen( fname, "rt" );
  fgets( buf, 2000, sf ); // header
  fgets( buf, 2000, sf ); // header
  for( int i = 0; i < 92; i++ )
  {
    fscanf( sf, "%*d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
      &scoef[i].mm1, &scoef[i].m1, &scoef[i].mnat,
      &scoef[i].rho, &scoef[i].atrho, &scoef[i].vfermi, &scoef[i].heat,
      &pcoef[i][0], &pcoef[i][1], &pcoef[i][2], &pcoef[i][3],
      &pcoef[i][4], &pcoef[i][5], &pcoef[i][6], &pcoef[i][7] );
    //printf( "%f %f\n", scoef[i].mm1, pcoef[i][7] );
  }
  fclose( sf );

  //snprintf( fname, 500, "%s/SLFCTR.dat", getenv("HOME") );
  snprintf( fname, 500, "%s/SLFCTR.dat", DATA_DIR );
  sf = fopen( fname, "rt" );
  fgets( buf, 2000, sf ); // header
  for( int i = 0; i < 92; i++ )
    fscanf( sf, "%*d %lf\n", &scoef[i].lfctr );
  fclose( sf );
}

