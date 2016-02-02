#include <cmath>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
//#include <stdlib.h>
#include <iostream>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "simconf.h"

namespace MyTRIM_NS {
  simconfType *simconf;
}

using namespace MyTRIM_NS;

simconfType::simconfType(Real _alfa)
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
  EelTotal = 0.0;
  EnucTotal = 0.0;

  // read data tables
  read_snuc();
  read_scoef();
}

void
simconfType::read_snuc()
{
  char * data_dir_env = getenv("MYTRIM_DATADIR");
  const unsigned int nbuf = 2000;
  char buffer[nbuf] = {0};
  if (data_dir_env  == 0)
    data_dir_env = strdup(MYTRIM_DATA_DIR);
  FILE *sf;

  snprintf(buffer, nbuf, "%s/%s", data_dir_env, "SNUC03.dat");
  sf = fopen(buffer, "rt");
  if (sf == NULL) fileReadError(buffer);

  for (int i = 0; i < 92; i++)
    for (int j = i; j < 92; j++)
    {
      fscanf(sf, "%*d %*d %lf %lf %lf %lf\n",
        &snuc[j][i][0], &snuc[j][i][1], &snuc[j][i][2], &snuc[j][i][3]);
      for (int n = 0; n < 4; n++)
        snuc[i][j][n] = snuc[j][i][n];
    }
  fclose(sf);
}

void
simconfType::read_scoef()
{
  char * data_dir_env = getenv("MYTRIM_DATADIR");
  const unsigned int nbuf = 2000;
  char buffer[nbuf] = {0};
  if (data_dir_env  == 0)
    data_dir_env = strdup(MYTRIM_DATA_DIR);
  FILE *sf;

  snprintf(buffer, nbuf, "%s/%s", data_dir_env, "SCOEF.95A");
  sf = fopen(buffer, "rt");
  if (sf == NULL) fileReadError(buffer);

  fgets(buffer, nbuf, sf); // header
  fgets(buffer, nbuf, sf); // header
  for (int i = 0; i < 92; i++)
  {
    fscanf(sf, "%*d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
      &scoef[i].mm1, &scoef[i].m1, &scoef[i].mnat,
      &scoef[i].rho, &scoef[i].atrho, &scoef[i].vfermi, &scoef[i].heat,
      &pcoef[i][0], &pcoef[i][1], &pcoef[i][2], &pcoef[i][3],
      &pcoef[i][4], &pcoef[i][5], &pcoef[i][6], &pcoef[i][7]);
  }
  fclose(sf);

  snprintf(buffer, nbuf, "%s/%s", data_dir_env, "SLFCTR.dat");
  sf = fopen(buffer, "rt");
  if (sf == NULL) fileReadError(buffer);

  fgets(buffer, nbuf, sf); // header
  for (int i = 0; i < 92; i++)
    fscanf(sf, "%*d %lf\n", &scoef[i].lfctr);
  fclose(sf);

  snprintf(buffer, nbuf, "%s/%s", data_dir_env, "ELNAME.dat");
  sf = fopen(buffer, "rt");
  if (sf == NULL) fileReadError(buffer);

  for (int i = 0; i < 92; i++)
    fscanf(sf, "%*d %s %s\n", scoef[i].sym, scoef[i].name);
  fclose(sf);
}

void
simconfType::fileReadError(const char * path)
{
#ifdef MYTRIM_ENABLED
    mooseError("Unable to open " << path);
#else
    std::cerr << "Unable to open " << path << std::endl;
    exit(1);
#endif
}
