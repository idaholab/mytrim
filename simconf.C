#include <cmath>
#include <string.h>
#include <stdlib.h>
//#include <stdlib.h>
#include <iostream>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "simconf.h"

namespace MyTRIM_NS {
  SimconfType *simconf;
}

using namespace MyTRIM_NS;

SimconfType::SimconfType(Real _alfa) :
    _data_dir(getenv("MYTRIM_DATADIR") ? getenv("MYTRIM_DATADIR") : MYTRIM_DATA_DIR)
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
  readSnuc();
  readScoef();
}

void
SimconfType::readSnuc()
{
  const unsigned int nbuf = 2000;
  char buffer[nbuf] = {0};
  FILE *sf;

  snprintf(buffer, nbuf, "%s/%s", _data_dir.c_str(), "SNUC03.dat");
  sf = fopen(buffer, "rt");
  if (sf == NULL) fileReadError(buffer);

  for (int i = 0; i < 92; i++)
    for (int j = i; j < 92; j++)
    {
      if (fscanf(sf, "%*d %*d %lf %lf %lf %lf\n",
        &snuc[j][i][0], &snuc[j][i][1], &snuc[j][i][2], &snuc[j][i][3]) != 4)
        fileReadError("contents of SNUC03.dat");

      for (int n = 0; n < 4; n++)
        snuc[i][j][n] = snuc[j][i][n];
    }
  fclose(sf);
}

void
SimconfType::readScoef()
{
  const unsigned int nbuf = 2000;
  char buffer[nbuf] = {0};
  FILE *sf;

  snprintf(buffer, nbuf, "%s/%s", _data_dir.c_str(), "SCOEF.95A");
  sf = fopen(buffer, "rt");
  if (sf == NULL) fileReadError(buffer);

  skipLine(sf); // header
  skipLine(sf); // header
  for (int i = 0; i < 92; i++)
  {
    if (fscanf(sf, "%*d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
      &scoef[i].mm1, &scoef[i].m1, &scoef[i].mnat,
      &scoef[i].rho, &scoef[i].atrho, &scoef[i].vfermi, &scoef[i].heat,
      &pcoef[i][0], &pcoef[i][1], &pcoef[i][2], &pcoef[i][3],
      &pcoef[i][4], &pcoef[i][5], &pcoef[i][6], &pcoef[i][7]) != 15)
      fileReadError("contents of SCOEF.95A");
  }
  fclose(sf);

  snprintf(buffer, nbuf, "%s/%s", _data_dir.c_str(), "SLFCTR.dat");
  sf = fopen(buffer, "rt");
  if (sf == NULL) fileReadError(buffer);

  skipLine(sf); // header
  for (int i = 0; i < 92; i++)
    if (fscanf(sf, "%*d %lf\n", &scoef[i].lfctr) != 1)
      fileReadError("contents of SLFCTR.dat");
  fclose(sf);

  snprintf(buffer, nbuf, "%s/%s", _data_dir.c_str(), "ELNAME.dat");
  sf = fopen(buffer, "rt");
  if (sf == NULL) fileReadError(buffer);

  for (int i = 0; i < 92; i++)
    if (fscanf(sf, "%*d %s %s\n", scoef[i].sym, scoef[i].name) != 2)
      fileReadError("contents of ELNAME.dat");
  fclose(sf);
}

void
SimconfType::fileReadError(const char * path)
{
#ifdef MYTRIM_ENABLED
    mooseError("Error reading " << path);
#else
    std::cerr << "Error reading " << path << std::endl;
    exit(1);
#endif
}

void
SimconfType::skipLine(FILE * sf)
{
  const unsigned int nbuf = 2000;
  char buffer[nbuf] = {0};
  if (!fgets(buffer, nbuf, sf))
  {
    #ifdef MYTRIM_ENABLED
      mooseError("Error reading file");
    #else
      std::cerr << "Error reading file" << std::endl;
      exit(1);
    #endif
  }
}
