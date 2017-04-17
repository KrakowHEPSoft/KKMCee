#ifndef HisNorm_H
#define HisNorm_H
//////////////////////////////////////////////////////////////////////////////
//                                                                          //
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;


// This works for 1-dim histograms
void HisNorm0( long   Nevt, double Xsav, TH1 *Hst);

/////////////////////////////////////////////////////////////////////
// This works for 1-dim histograms
void HisNorm1(TH1D *NorHst, TH1 *Hst);

/////////////////////////////////////////////////////////////////////
void HisNorm1M(TH1D *NorHst, TH1D *Hst, Float_t msize, Int_t mcolor, Int_t mark);

/////////////////////////////////////////////////////////////////////
// This works for 2-dim histograms
void HisNorm2(TH1D *NorHst, TH2 *Hst);

///////////////////////////////////////////////////////////////////////////////////
void ProjX1(TH2D *Scat, TH1D *&HstProjX);

///////////////////////////////////////////////////////////////////////////////////
void ProjY1(TH2D *Scat, TH1D *&HstProjY);

///////////////////////////////////////////////////////////////////////////////////
void ProjV(TH2D *Scat, TH1D *&hxTot, TH1D *&hxAfb, int NbMax);

///////////////////////////////////////////////////////////////////////////////////
void ProjC(TH2D *Scat, TH1D *&hTot, TH1D *&hAsy, int NbMax);

///////////////////////////////////////////////////////////////////////////////////
void MakeCumul(TH1D *hst1, TH1D *&hcum1);

  ///////////////////////////////////////////////////////////////////////////////////
void MakeAFB(TH1D *hAll, TH1D *&hAFB);

////////////////////////////////////////////////////////////////////////////
#endif
