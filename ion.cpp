#include <math.h>

#include "ion.h"

ionBase::ionBase()
{ 
  ef = 5.0; // final energy
  t = 0.0; //clock
}

void ionBase::set_ef()
{ 
  ef = fmax( 5.0, 0.00001 * e ); // final energy
}
