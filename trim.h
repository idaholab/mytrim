#ifndef TRIM_H
#define TRIM_H 1

#include <vector>
#include <queue>

using namespace std;

#include "material.h"
#include "sample.h"

class trimBase {
public:
  void trim( ionBase *pka, queue<ionBase*> &recoils );
  trimBase( sampleBase *sample_ ) : sample( sample_ ) {};

protected:
  sampleBase *sample;
  ionBase *pka, *recoil;
  materialBase *material;
  elementBase *element;
  queue<ionBase*> *recoil_queue_ptr;
  bool terminate;

  // by default only follow recoils with E > 12eV
  virtual bool followRecoil() {
    // TODO: find a better place for this!
    /*if( pka->md > 0 )
      recoil->md = pka->md +1;
    else
      recoil->md = 0;
    */
    return recoil->e > 12.0;
  };
  virtual void vacancyCreation();
  virtual void checkPKAState() {};
};

// //
// // Do a breadth first rather than depth first recoil simulation
// //
// class trimBreadthFirst : public trimBase {
// public:
//   trimBreadthFirst( sampleBase *sample_ ) : trimBase( sample_ ) {};
// protected:
//   virtual bool followRecoil() {
//     recoil_queue_ptr->push(pka);
//     terminate = true;
//     return recoil->e > 12.0;
//   };
// };


//
// Only follow the primary knock ons
//
class trimPrimaries : public trimBase {
public:
  trimPrimaries( sampleBase *sample_ ) : trimBase( sample_ ) {};
protected:
  virtual bool followRecoil() { return false; };
};


//
// Only follow the first generation of recoils
//
class trimRecoils : public trimBase {
  public:
    trimRecoils( sampleBase *sample_ ) : trimBase( sample_ ) {};
  protected:
    virtual bool followRecoil() { return ( recoil->gen <= 2 ); };
};


//
// store a history of all recoils
//
class trimHistory : public trimBase {
public:
  trimHistory( sampleBase *sample_ ) : trimBase( sample_ ) {};
  vector<double> pos_hist[3];
protected:
  virtual bool followRecoil()
  {
    pos_hist[0].push_back(pka->pos[0]);
    pos_hist[1].push_back(pka->pos[1]);
    pos_hist[2].push_back(pka->pos[2]);
    return recoil->e > 10.0;
  };
};


//
// Full damage cascade, plus log vaccancy creation
//
class trimVacLog : public trimBase {
public:
  trimVacLog( sampleBase *sample_, FILE *vacfile_ ) : vacfile(vacfile_), trimBase( sample_ ) {};
protected:
  FILE *vacfile;
  virtual void vacancyCreation()
  {
    // both atoms have enough energy to leave the site
    // log the original recoil atom position as a vacancy site
    fprintf( vacfile, "%f %f %f %d\n", recoil->pos[0], recoil->pos[1], recoil->pos[2], recoil->z1 );
  };
};

//
// Full damage cascade, plus map vaccancy creation
//
class trimVacMap : public trimBase {
  static const int mx = 20, my = 20;
public:
  int vmap[mx][my][3];
  trimVacMap( sampleBase *sample_, int z1_, int z2_, int z3_ = -1 ) : trimBase( sample_ ), z1(z1_), z2(z2_), z3(z3_)
  {
    for( int e = 0; e < 3; e++ )
      for( int x = 0; x < mx; x++ )
        for( int y = 0; y < my; y++ )
          vmap[x][y][e] = 0;
  };
protected:
  int z1, z2, z3;
  virtual void vacancyCreation()
  {
    // both atoms have enough energy to leave the site
    int x, y;

    x = ( ( recoil->pos[0] * mx ) / sample->w[0] );
    y = ( ( recoil->pos[1] * my ) / sample->w[1] );
    x -= int(x/mx) * mx;
    y -= int(y/my) * my;

    // keep track of vaccancies for the two constituents
    if( recoil->z1 == z1 ) vmap[x][y][0]++;
    else if( recoil->z1 == z2 ) vmap[x][y][1]++;
    else if( recoil->z1 == z3 ) vmap[x][y][2]++;
  };
};


//
// Output all phonon energy losses
//
class trimPhononOut : public trimBase {
public:
  trimPhononOut( sampleBase *sample_, FILE *phonfile_ ) : phonfile(phonfile_), trimBase( sample_ ) {};
protected:
  FILE *phonfile;
  virtual bool followRecoil()
  {
    // dissipate lattice binding energy where recoil branches off
    double Edep = recoil->e + element->Elbind;
    fprintf( phonfile, "%f %f %f %f %d %d %e 0\n", Edep, recoil->pos[0], recoil->pos[1], recoil->pos[2], pka->z1, pka->id, pka->t );
    return (recoil->e > 10.0);
  };
};

#endif
