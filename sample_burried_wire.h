#ifndef SAMPLE_BURRIED_WIRE_H
#define SAMPLE_BURRIED_WIRE_H

#include "sample_wire.h"

namespace MyTRIM_NS {

struct sampleBurriedWire : sampleWire {
  sampleBurriedWire( double x = 10000.0, double y = 10000.0, double z = 10000.0 );

  virtual materialBase* lookupMaterial( double* pos );
};

}

#endif
