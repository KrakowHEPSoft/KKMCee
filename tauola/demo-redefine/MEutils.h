/*
   These are C++ versions of functions from tauola.f
   Needed for BaBar pi- pi0 current
*/
#include <complex>
using std::complex;

extern "C" {
  // determines whether actual decay is of tau+ or tau-

  // for the calculation of "AXIAL TYPE"  PI-VECTOR  PIA in tau decay matrix element,
  // we must know whether the actual decay is of tau+ or of tau-
  // for that purpose information  from Tauola common blocks must be used.
  // The signtau(int SIGN) should migrate from pipi0.c file it does not depend on the 
  // Tau decay mode. 

extern void taugetsign_(float SIGN1);
}


// clvec will be shifted to private(?) section of utilities class
void clvec(complex<float> *HJ, const float *PN, float *PIV);
// claxi will be shifted to private(?) section of utilities class
void claxi(complex<float> *HJ, const float *PN, float *PIA);
// clnut will be shifted to private(?) section of utilities class
void clnut(complex<float> *HJ, float &B, float *HV);

// tauDecay_ME will be shifted to public(?) section of utilities class
void tauDecay_ME(const int &ifcosCabibo, const float *pt, const float *pn, complex<float> *hadcur, float &amplit, float *hv);

// bwig will be shifted to current specific utilities class
complex<float> bwig(float s, float m, float g);
// fpik will be shifted to current specific utilities class
complex<float> fpik(float W);
