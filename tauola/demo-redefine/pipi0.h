/*
   These are C++ versions of functions from tauola.f
   Needed for BaBar pi- pi0 current
*/
#include <complex>
using std::complex;


void pipi0_curr(const float *pc, const float *pn, complex<float> *hadcur);
void pipi0_ME  (const float *pt, const float *pn, const float *pim1, const float *pim2, float &amplit, float *hv);
