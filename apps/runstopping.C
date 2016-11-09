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

  if (json_root["stopping"].isObject())
    json_root = json_root["stopping"];
  else
    mytrimError("No 'stopping' top level block found in input");

  // initialize global parameter structure and read data tables from file
  SimconfType * simconf = new SimconfType;
  simconf->fullTraj = false;
  simconf->tmin = 0.2;

  //
  // process material block
  //
  if (!json_root["material"].isObject())
    mytrimError("Must specify a 'material' block in the input file");


  MaterialBase * material;
  Element element;

  if (!json_root["material"]["rho"].isNumeric())
    mytrimError("Missing 'rho'");
  const Real rho   = json_root["material"]["rho"].asDouble();

  if (!json_root["material"]["elements"].isArray())
    mytrimError("Missing 'elements' in material");
  Json::Value json_elem = json_root["material"]["elements"];

  material = new MaterialBase(simconf, rho);

  for (unsigned int j = 0; j < json_elem.size(); ++j)
  {
    if (!json_elem[j]["Z"].isNumeric())
      mytrimError("Missing 'Z' in element " << j);
    element._Z = json_elem[j]["Z"].asInt();

    if (!json_elem[j]["mass"].isNumeric())
      mytrimError("Missing 'mass' in element " << j);
    element._m = json_elem[j]["mass"].asDouble();

    if (!json_elem[j]["fraction"].isNumeric())
      mytrimError("Missing 'fraction' in element " << j);
    element._t = json_elem[j]["fraction"].asDouble();

    material->_element.push_back(element);
  }

  material->prepare(); // all elements added

  //
  // process ion block
  //
  if (!json_root["ion"].isObject())
    mytrimError("Must specify an 'ion' block in the input file");

  IonBase * pka = new IonBase;

  if (!json_root["ion"]["Z"].isNumeric())
    mytrimError("Missing 'Z' in ion block");
  pka->_Z = json_root["ion"]["Z"].asInt();

  if (!json_root["ion"]["mass"].isNumeric())
    mytrimError("Missing 'mass' in ion block");
  pka->_m = json_root["ion"]["mass"].asDouble();

  // construct list of energies (can be specified as single number, array, or range descriptor)
  std::vector<Real> energies;
  if (json_root["ion"]["energy"].isNumeric())
    energies.push_back(json_root["ion"]["energy"].asDouble());
  else if (json_root["ion"]["energy"].isArray())
  {
    Json::Value json_energy = json_root["ion"]["energy"];
    for (unsigned int j = 0; j < json_energy.size(); ++j)
      energies.push_back(json_energy[j].asDouble());
  }
  else if (json_root["ion"]["energy"].isObject())
  {
    Json::Value json_energy = json_root["ion"]["energy"];

    if (!json_energy["begin"].isNumeric())
      mytrimError("Missing 'begin' in energy block");
    const Real ebegin = json_energy["begin"].asDouble();

    if (!json_energy["end"].isNumeric())
      mytrimError("Missing 'end' in energy block");
    const Real eend = json_energy["end"].asDouble();

    const bool hasStep = json_energy["step"].isNumeric();
    const bool hasMult = json_energy["mult"].isNumeric();
    if (hasStep == hasMult)
      mytrimError("Specify either 'step' or 'mult' energy block");
    if (hasStep)
    {
      const Real estep = json_energy["step"].asDouble();
      for (Real E = ebegin; E <= eend; E += estep)
        energies.push_back(E);
    }
    if (hasMult)
    {
      const Real emult = json_energy["mult"].asDouble();
      if (emult <= 1.0)
        mytrimError("'mult' must be larger than 1.0");
      for (Real E = ebegin; E <= eend; E *= emult)
        energies.push_back(E);
    }
  }
  else
    mytrimError("Missing or invalid 'energy' in ion block");

  for (unsigned int n = 0; n < energies.size(); ++n)
  {
    pka->_E = energies[n];
    std::cout << pka->_E << ' ' << material->getrstop(pka) << '\n';
  }

  return EXIT_SUCCESS;
}
