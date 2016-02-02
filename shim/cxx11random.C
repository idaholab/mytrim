#include <random>
#include "cxx11random.h"

namespace MyTRIM_NS {

static std::mt19937 * cxx11random_gen = NULL;
static std::uniform_real_distribution<double> cxx11random_dis_Real(0, 1);
static std::uniform_int_distribution<int> cxx11random_dis_int(0, 65535);

/**** Function: r250_init
        Description: initializes r250 random number generator. ****/
void r250_init(int seed)
{
  delete cxx11random_gen;
  cxx11random_gen = new std::mt19937(seed);
}

/**** Function: r250 Description: returns a random unsigned integer k
                                uniformly distributed in the interval 0 <= k < 65536.  ****/
unsigned int r250()
{
  return cxx11random_dis_int(*cxx11random_gen);
}

/**** Function: dr250
                Description: returns a random double z in range 0 <= z < 1.  ****/
double dr250()
{
  return cxx11random_dis_Real(*cxx11random_gen);
}

}
