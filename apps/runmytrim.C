/***************************************************************************
 *   Copyright (C) 2016 by Daniel Schwen   *
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
#include <list>
#include <iostream>
#include <string>
#include <limits>
#include <thread>

#include "simconf.h"
#include "element.h"
#include "material.h"
#include "sample_layers.h"
#include "ion.h"
#include "invert.h"
#include "functions.h"

// TRIM modules
#include "include/ThreadedTrimBase.h"
#include "include/TrimRange.h"
#include "include/TrimVacCount.h"
#include "include/TrimVacEnergyCount.h"

using namespace MyTRIM_NS;

#define mytrimError(msg)                                                                           \
  do                                                                                               \
  {                                                                                                \
    std::cerr << "ERROR: " << msg << '\n' << std::flush;                                           \
    return 1;                                                                                      \
  } while (0)

// thread data
struct ThreadData
{
  SimconfType _simconf;
  SampleLayers * _sample;
  ThreadedTrimBase * _trim;
  std::queue<IonBase *> _recoils;
  std::list<IonBase *> _pka;
};
std::vector<ThreadData> thread_data;

// MC loop threads
void
computeThread(int tid)
{
  auto & td = thread_data[tid];

  for (auto pka : td._pka)
  {
    td._simconf.seed(pka->_seed);
    td._recoils.push(pka);

    while (!td._recoils.empty())
    {
      IonBase * recoil = td._recoils.front();
      td._recoils.pop();
      td._sample->averages(recoil);

      td._trim->trim(recoil, td._recoils);

      // done with this recoil
      delete recoil;
    }
  }
}

int
main(int argc, char * argv[])
{
  // error out if any cli args have been passed in
  if (argc > 1)
    mytrimError("Please supply the input file via stdin (e.g. ./runmytrim < input.json`)");

  // open the input
  Json::Value json_root;
  std::cin >> json_root;

  // number of threads
  unsigned int nthreads = 1;

  if (json_root["mytrim"].isObject())
    json_root = json_root["mytrim"];
  else
    mytrimError("No 'mytrim' top level block found in input");

  // first look at the thread number option
  if (json_root["options"].isObject() && json_root["options"]["threads"].isNumeric())
  {
    nthreads = json_root["options"]["threads"].asInt();
    std::cerr << "Using " << nthreads << " threads\n";
  }

  // set up thread data vector
  thread_data.resize(nthreads);
  for (auto & td : thread_data)
  {
    td._simconf.fullTraj = false;
    td._simconf.tmin = 0.2;
  }

  //
  // process the "options" block
  //
  int master_seed;
  if (json_root["options"].isObject())
  {
    // random seed
    if (json_root["options"]["seed"].isNumeric())
    {
      master_seed = json_root["options"]["seed"].asInt();
      std::cerr << "Using provided master seed " << master_seed << '\n';
    }
    else
    {
      // seed randomnumber generator from system entropy pool
      FILE * urand = fopen("/dev/random", "r");
      if (fread(&master_seed, sizeof(int), 1, urand) != 1)
        mytrimError("Unable to access /dev/random");
      fclose(urand);
    }

    double scale;
    if (json_root["options"]["scale"].isNumeric())
    {
      scale = json_root["options"]["scale"].asDouble();
      std::cerr << "Using provided length scale " << scale << '\n';
      for (auto & td : thread_data)
        td._simconf.setLengthScale(scale);
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
  for (auto & td : thread_data)
    td._sample = new SampleLayers(thickness, 100.0, 100.0);

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
    for (auto & td : thread_data)
      td._trim = new TrimVacCount(&(td._simconf), td._sample);
  else if (output_type == "vacenergycount")
    for (auto & td : thread_data)
      td._trim = new TrimVacEnergyCount(&(td._simconf), td._sample);
  else if (output_type == "range")
    for (auto & td : thread_data)
      td._trim = new TrimRange(&(td._simconf), td._sample);
  else
    mytrimError("Unknown output type " << output_type);

  if (json_root["output"]["base"].isString())
    for (auto & td : thread_data)
      td._trim->setBaseName(json_root["output"]["base"].asString());

  if (json_root["output"]["primaries_only"].isBool())
    for (auto & td : thread_data)
      td._trim->_primaries_only = json_root["output"]["primaries_only"].asBool();

  MaterialBase * material;
  Element element;
  for (auto & td : thread_data)
    for (unsigned int i = 0; i < nlayers; ++i)
    {
      if (!json_layers[i]["rho"].isNumeric())
        mytrimError("Missing 'rho' in layer " << i);
      Real lrho = json_layers[i]["rho"].asDouble();

      if (!json_layers[i]["elements"].isArray())
        mytrimError("Missing 'elements' in layer " << i);
      Json::Value json_elem = json_layers[i]["elements"];

      material = new MaterialBase(&(td._simconf), lrho); // rho

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

        if (json_elem[j]["edisp"].isNumeric())
          element._Edisp = json_elem[j]["edisp"].asDouble();

        if (json_elem[j]["elbind"].isNumeric())
          element._Elbind = json_elem[j]["elbind"].asDouble();

        material->_element.push_back(element);
      }

      material->prepare();                      // all elements added
      td._sample->material.push_back(material); // add material to sample
      td._sample->layerThickness.push_back(json_layers[i]["thickness"].asDouble());
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

  if (json_root["ion"]["final_energy"].isNumeric())
    pkaTemplate->_Ef = json_root["ion"]["final_energy"].asDouble();

  // use the rng in thread 0 to generate the pka seeds from the master seed
  thread_data[0]._simconf.seed(master_seed);

  // fill PKA queues
  for (unsigned int n = 0; n < npka; n++)
  {
    auto & td = thread_data[n % nthreads];
    IonBase * pka = new IonBase(pkaTemplate);
    pka->_gen = 0; // generation (0 = PKA)
    pka->_dir = Point(1.0, 0.0, 0.0);

    Point start(0.0, td._sample->w[1] / 2.0, td._sample->w[2] / 2.0);
    pka->_pos = start;

    // pka->setEf();
    pka->_seed = thread_data[0]._simconf.irand();
    td._pka.push_back(pka);
  }

  // launch threads
  std::vector<std::thread> thread(nthreads);
  for (unsigned int i = 0; i < nthreads; ++i)
    thread[i] = std::thread(computeThread, i);

  // wait for thread completion
  for (unsigned int i = 0; i < nthreads; ++i)
    thread[i].join();

  // join trim data into thread 0
  for (unsigned int i = 1; i < nthreads; ++i)
  {
    thread_data[0]._trim->threadJoin(*thread_data[i]._trim);

    thread_data[0]._simconf.vacancies_created += thread_data[i]._simconf.vacancies_created;
    thread_data[0]._simconf.EelTotal += thread_data[i]._simconf.EelTotal;
    thread_data[0]._simconf.EnucTotal += thread_data[i]._simconf.EnucTotal;
  }

  // write output files
  thread_data[0]._trim->writeOutput();

  // summary
  std::cerr << "Vacancies/ion: " << Real(thread_data[0]._simconf.vacancies_created) / Real(npka)
            << '\n';

  return EXIT_SUCCESS;
}
