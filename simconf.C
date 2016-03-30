#include "simconf.h"
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cmath>
#include <cstdlib>

namespace MyTRIM_NS {
  SimconfType *simconf;
}

using namespace MyTRIM_NS;

SimconfType::SimconfType(Real _alfa) :
    _data_dir(std::getenv("MYTRIM_DATADIR") ? std::getenv("MYTRIM_DATADIR") : MYTRIM_DATA_DIR)
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
  readDataFiles();
}

void
SimconfType::readDataFiles()
{
  std::string fname;
  int dummy;

  fname = _data_dir + "/SNUC03.dat";
  std::ifstream snuc03(fname.c_str());
  if (!snuc03)
    fileReadError(fname);

  for (int i = 0; i < 92; ++i)
    for (int j = i; j < 92; j++)
    {
      snuc03 >> dummy >> dummy
             >> snuc[j][i][0] >> snuc[j][i][1] >> snuc[j][i][2] >> snuc[j][i][3];

      if (!snuc03)
        fileReadError("contents of SNUC03.dat");

      for (int n = 0; n < 4; n++)
        snuc[i][j][n] = snuc[j][i][n];
    }


  fname = _data_dir + "/SCOEF.95A";
  std::ifstream scoef95a(fname.c_str());
  if (!scoef95a)
    fileReadError(fname);

  skipLine(scoef95a); // header
  skipLine(scoef95a); // header
  for (int i = 0; i < 92; ++i)
  {
    scoef95a >> dummy
             >> scoef[i].mm1 >> scoef[i].m1 >> scoef[i].mnat
             >> scoef[i].rho >> scoef[i].atrho >> scoef[i].vfermi >> scoef[i].heat
             >> pcoef[i][0] >> pcoef[i][1] >> pcoef[i][2] >> pcoef[i][3]
             >> pcoef[i][4] >> pcoef[i][5] >> pcoef[i][6] >> pcoef[i][7];

    if (!scoef95a)
      fileReadError("contents of SCOEF.95A");
  }


  fname = _data_dir + "/SLFCTR.dat";
  std::ifstream slfctr(fname.c_str());
  if (!slfctr)
    fileReadError(fname);

  skipLine(slfctr); // header
  for (int i = 0; i < 92; ++i)
    if (!(slfctr >> dummy >> scoef[i].lfctr))
      fileReadError("contents of SLFCTR.dat");


  fname = _data_dir + "/ELNAME.dat";
  std::ifstream elname(fname.c_str());
  if (!elname)
    fileReadError(fname);

  for (int i = 0; i < 92; ++i)
    if (!(elname >> dummy >> scoef[i].sym >> scoef[i].name))
      fileReadError("contents of ELNAME.dat");
}

void
SimconfType::fileReadError(const std::string & path)
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

void
SimconfType::skipLine(std::ifstream & sf)
{
  std::string s;
  if (!std::getline(sf, s))
  {
    #ifdef MYTRIM_ENABLED
      mooseError("Error reading file");
    #else
      std::cerr << "Error reading file" << std::endl;
      exit(1);
    #endif
  }
}
