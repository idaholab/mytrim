/***************************************************************************
 *   Copyright (C) 2008 by Daniel Schwen   *
 *   daniel@schwen.de   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <json/json.h>

#include <stdio.h>
#include <stdlib.h>
#include <queue>
#include <iostream>
#include <string>
#include <limits>

#include "simconf.h"
#include "element.h"
#include "material.h"
#include "sample_layers.h"
#include "ion.h"
#include "trim.h"
#include "invert.h"
#include "functions.h"

// TRIM modules
#include "include/TrimVacCount.h"
#include "include/TrimVacEnergyCount.h"

using namespace MyTRIM_NS;

#define mytrimError(msg)                                 \
  do {                                                   \
    std::cerr << "ERROR: " << msg << '\n' << std::flush; \
    return 1;                                            \
  } while(0)

int main()
{
  // open the input
  Json::Value json_root;
  std::cin >> json_root;

  if (json_root["mytrim"].isObject())
    json_root = json_root["mytrim"];
  else
    mytrimError("No 'mytrim' top level block found in input");

  // initialize global parameter structure and read data tables from file
  SimconfType * simconf = new SimconfType;
  simconf->fullTraj = false;
  simconf->tmin = 0.2;

  //
  // process the "options" block
  //
  if (json_root["options"].isObject())
  {
    // random seed
    int seed;
    if (json_root["options"]["seed"].isNumeric())
    {
      seed = json_root["options"]["seed"].asInt();
      std::cerr << "Using provided seed " << seed << '\n';
    }
    else
    {
      // seed randomnumber generator from system entropy pool
      FILE *urand = fopen("/dev/random", "r");
      if (fread(&seed, sizeof(int), 1, urand) != 1)
        mytrimError("Unable to access /dev/random");
      fclose(urand);
    }
    simconf->seed(seed < 0 ? -seed : seed); // random generator goes haywire with neg. seed

    double scale;
    if (json_root["options"]["scale"].isNumeric())
    {
      scale = json_root["options"]["scale"].asDouble();
      std::cerr << "Using provided length scale " << scale << '\n';
      simconf->setLengthScale(scale);
    }
  }

  //
  // process sample block
  //
  if (!json_root["sample"].isObject())
    mytrimError("Must specify a 'sample' block in the input file");

  if (!json_root["sample"]["layers"].isArray())
    mytrimError("sample.layers must be an array");

  // sample layers
  Json::Value json_layers = json_root["sample"]["layers"];
  unsigned int nlayers = json_layers.size();
  std::cerr << "Building " << nlayers << " layers\n";

  // calculate total thickness
  Real thickness = 0.0;
  for (unsigned int i = 0; i < nlayers; ++i)
    if (json_layers[i]["thickness"].isNumeric())
      thickness += json_layers[i]["thickness"].asDouble();
    else
      mytrimError("No 'thickness' found for layer " << i);

  // initialize sample structure
  SampleLayers *sample = new SampleLayers(thickness, 100.0, 100.0);

  // set up TRIM module
  TrimBase * trim;

  //
  // process output block
  //
  if (!json_root["output"].isObject())
    mytrimError("Must specify an 'output' block in the input file");

  if (!json_root["output"]["type"].isString())
    mytrimError("output.type must be a string");
  const std::string output_type = json_root["output"]["type"].asString();

  // construct TRIM object according to output type
  if (output_type == "vaccount")
    trim = new TrimVacCount(simconf, sample);
  else if (output_type == "vacenergycount")
    trim = new TrimVacEnergyCount(simconf, sample);
  else
    mytrimError("Unknown output type " << output_type);

  if (json_root["output"]["base"].isString())
    trim->setBaseName(json_root["output"]["base"].asString());


  MaterialBase * material;
  Element element;
  for (unsigned int i = 0; i < nlayers; ++i)
  {
    if (!json_layers[i]["rho"].isNumeric())
      mytrimError("Missing 'rho' in layer " << i);
    Real lrho   = json_layers[i]["rho"].asDouble();

    if (!json_layers[i]["elements"].isArray())
      mytrimError("Missing 'elements' in layer " << i);
    Json::Value json_elem = json_layers[i]["elements"];

    material = new MaterialBase(simconf, lrho); // rho

    for (unsigned int j = 0; j < json_elem.size(); ++j)
    {
      if (!json_elem[j]["Z"].isNumeric())
        mytrimError("Missing 'Z' in element " << j << " in layer " << i);
      element._Z = json_elem[j]["Z"].asInt();

      if (!json_elem[j]["mass"].isNumeric())
        mytrimError("Missing 'mass' in element " << j << " in layer " << i);
      element._m = json_elem[j]["mass"].asDouble();

      if (!json_elem[j]["fraction"].isNumeric())
        mytrimError("Missing 'fraction' in element " << j << " in layer " << i);
      element._t = json_elem[j]["fraction"].asDouble();

      material->_element.push_back(element);
    }

    material->prepare(); // all elements added
    sample->material.push_back(material); // add material to sample
    sample->layerThickness.push_back(json_layers[i]["thickness"].asDouble());
  }

  //
  // process ion block
  //
  if (!json_root["ion"].isObject())
    mytrimError("Must specify an 'ion' block in the input file");

  IonBase * pkaTemplate = new IonBase;

  if (!json_root["ion"]["Z"].isNumeric())
    mytrimError("Missing 'Z' in ion block");
  pkaTemplate->_Z = json_root["ion"]["Z"].asInt();

  if (!json_root["ion"]["mass"].isNumeric())
    mytrimError("Missing 'mass' in ion block");
  pkaTemplate->_m = json_root["ion"]["mass"].asDouble();

  if (!json_root["ion"]["energy"].isNumeric())
    mytrimError("Missing 'energy' in ion block");
  pkaTemplate->_E = json_root["ion"]["energy"].asDouble();

  if (!json_root["ion"]["number"].isNumeric())
    mytrimError("Missing 'number' in ion block");
  unsigned int npka = json_root["ion"]["number"].asInt();


  // start output
  trim->startOutput();

  // create a FIFO for recoils
  std::queue<IonBase*> recoils;

  IonBase *pka;
  Point start(0.0, sample->w[1] / 2.0, sample->w[2] / 2.0);

  for (unsigned int n = 0; n < npka; n++)
  {
    if (n % 100 == 0)
      std::cerr << "pka #" << n+1 << '\n';

    pka = new IonBase(pkaTemplate);
    pka->_gen = 0; // generation (0 = PKA)
    pka->_id = simconf->_id++;

    pka->_dir = Point(1.0, 0.0, 0.0);
    pka->_pos = start;

    pka->setEf();
    recoils.push(pka);

    while (!recoils.empty())
    {
      pka = recoils.front();
      recoils.pop();
      sample->averages(pka);

      // do ion analysis/processing BEFORE the cascade here
      // opos = pka->_pos;

      // follow this ion's trajectory and store recoils
      //if (pka->_Z == 29 || pka->_Z == Z)
      //if (pka->_Z == 29 || pka->_Z == Z)
      trim->trim(pka, recoils);

      // do ion analysis/processing AFTER the cascade here
      // if (pka->_Z != Z)
      // {
      //   sum_r2 += (opos - pka->_pos).norm_sq();
      //   nrec++;
      // }
      //
      // // pka is O or Ag
      // //if (pka->_Z == 29 && pka->_pos(0) >= 500.0)
      // if (pka->_Z == 29)
      // {
      //   // output
      //   printf("RP %f %d %d\n", pka->_pos(0), n,  pka->_gen);
      // }

      // done with this recoil
      delete pka;
    }
  }

  // stop output
  trim->stopOutput();

  std::cerr << "Vacancies/ion: " << Real(simconf->vacancies_created)/Real(npka) << '\n';

  return EXIT_SUCCESS;
}
