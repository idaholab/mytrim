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

  // by default only follow recoils with E > 12eV
  virtual bool spawnRecoil() { return recoil->e > 12.0; };
  virtual void vacancyCreation();
};

//
// Only follow the primary knock ons
//
class trimPrimaries : public trimBase {
public:
  trimPrimaries( sampleBase *sample_ ) : trimBase( sample_ ) {};
protected:
  virtual bool spawnRecoil() { return false; };
};

//
// Only follow the first generation of recoils
//
class trimRecoils : public trimBase {
  public:
    trimRecoils( sampleBase *sample_ ) : trimBase( sample_ ) {};
  protected:
    virtual bool spawnRecoil() { return ( recoil->gen <= 2 ); };
};


//
// store a history of all recoils
//
class trimHistory : public trimBase {
public:
  trimHistory( sampleBase *sample_ ) : trimBase( sample_ ) {};
  vector<double> pos_hist[3];
protected:
  virtual bool spawnRecoil() 
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
  virtual bool spawnRecoil() 
  {
    // both atoms must have enough energy to leave the site
    if( recoil->e > 10 )
      fprintf( vacfile, "%f %f %f %d\n", recoil->pos[0], recoil->pos[1], recoil->pos[2], recoil->z1 );
    
    // if pka->e <= pka->ef the resulting interstitial will be at the exact coordinates 
    // of the vaccancy and immediate annihilate

    return recoil->e > 10.0; 
  };
};

//
// Full damage cascade, plus map vaccancy creation
//
class trimVacMap : public trimBase {
  static const int mx = 20, my = 20;
public:
  int vmap[mx][my][3];
  trimVacMap( sampleBase *sample_, int z1_, int z2_ ) : trimBase( sample_ ), z1(z1_), z2(z2_) 
  {
    for( int e = 0; e < 3; e++ )
      for( int x = 0; x < mx; x++ )
        for( int y = 0; y < my; y++ )
          vmap[x][y][e] = 0;
  };
protected:
  int z1, z2;
  virtual bool spawnRecoil() 
  {
    // both atoms must have enough energy to leave the site
    int x, y;
    if( recoil->e > 10 )
    {
      x = ( ( recoil->pos[0] * mx ) / sample->w[0] );
      y = ( ( recoil->pos[1] * my ) / sample->w[1] );
      x -= int(x/mx) * mx;
      y -= int(y/my) * my;

      // keep track of vaccancies for the two constituents
      if( recoil->z1 == z1 ) vmap[x][y][0]++;
      else if( recoil->z1 == z2 ) vmap[x][y][1]++;
      else vmap[x][y][2]++; // this should never happen...
    }
    // if pka->e <= pka->ef the resulting interstitial will be at the exact coordinates 
    // of the vaccancy and immediate annihilate

    return recoil->e > 10.0; 
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
  virtual bool spawnRecoil() 
  { 
    if( recoil->e > 10.0 ) return true;
    else
    {
      fprintf( phonfile, "%f %f %f %f %d %d %e\n", recoil->e, recoil->pos[0], recoil->pos[1], recoil->pos[2], pka->z1, pka->id, pka->t );  
      return false;
    }
  };
};

#endif
