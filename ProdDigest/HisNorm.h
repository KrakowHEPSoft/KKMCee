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

#include "TH2.h"

double sqr( const Double_t x );

// calculates no. of equiv. WT=1 events
double MCequiv(TH1D *hst1);

// This works for 1-dim histograms
void HisNorm0( long   Nevt, double Xsav, TH1 *Hst);

// Renormalizes 1-dim histograms to Xsav
void HisNorm0(double Xsav, TH1 *Hst);

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
TH1D *HstProjC(TString title, TH2D *Scat, int NbMax);

///////////////////////////////////////////////////////////////////////////////////
TH1D *HstProjCA(TString title, TH2D *Scat, int NbMax);

///////////////////////////////////////////////////////////////////////////////////
void MakeCumul(TH1D *hst1, TH1D *&hcum1);

TH1D *HstCumul(TString title, TH1D *hst1);

TH1D *HstProjV(TString title, TH2D *&Scat, int NbMax);

TH1D *HstProjAv(TString title, TH2D *&Scat, int NbMax);

TH1D *HstProjA(TString title, TH2D *&Scat, int NbMax);

TH1D *HstProjF(TString title, TH2D *&Scat, int NbMax);

TH1D *HstEvent(TString title, TH1D *hst1, double IntLumi);

///////////////////////////////////////////////////////////////////////////////////
void MakeAFB(TH1D *hAll, TH1D *&hAFB);

TH1D *HstDiff(TString title, TH1D *HST1, TH1D *HST2, Int_t kolor);

TH1D *HstAddi(TString title, TH1D *HST1, TH1D *HST2, Int_t kolor);

TH1D *HstRatio(TString title, TH1D *HST1, TH1D *HST2, Int_t kolor);

TH1D *HstRatioSc(TString title, TH1D *HST1, TH1D *HST2, Double_t fact);

TH1D *HstTildeAFB(TString title, TH1D *HST1, TH1D *HST2);

TH1D *HstAFB(TString title, TH1D *HST1, TH1D *HST2);

TH1D *HstAFB3(TString title, TH1D *HST1, TH1D *HST2, TH1D *HST3);

TH1D *HstAFB4(TString title, TH1D *HST21F, TH1D *HST21, TH1D *HST2F, TH1D *HST2);

TH1D *HstAFB2cl(TString title, TH1D *HST1F, TH1D *HST1);

TH1D *HstAFB4cl(TString title, TH1D *HST21F, TH1D *HST21, TH1D *HST2F, TH1D *HST2);

///////////////////////////////////////////////////////////////////////////////////
void PlInitialize(FILE *ltx, int lint);

void PlTable2(int Ncol, TH1D *iHst[], FILE *ltex, Char_t *Capt[], Char_t Mcapt[] , const char *chr1, int k1,int k2,int dk);

void PlEnd(FILE *ltex);

////////////////////////////////////////////////////////////////////////////
#endif
