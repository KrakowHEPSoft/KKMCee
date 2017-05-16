#ifndef TKfig_H
#define TKfig_H

// C++ headers
using namespace std;
#include<stdlib.h>
#include<fstream>
#include <iostream>
#include<iomanip>
#include<math.h>

// ROOT headers
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"

#include <TClass.h>
#include <TFile.h>
#include <TKey.h>
///////////////////////////////////////////////////////////////////////////////
///   Librarry of auxiliary procedures for digesting MC results
///   mainly functions with analytical results and tools for histograms
///////////////////////////////////////////////////////////////////////////////

class TKfig {

 public:
  TString m_PlotSelector;

  TKfig();
  virtual ~TKfig(){};                         // Destructor
  virtual void HistNorm(TH1D *NorHst, TH1D *Hst);
  virtual void HistNorm2(TH1D *NorHst, TH2D *Hst);
  virtual void HistNormalize(TString normname, TFile *file);

  virtual TH1D *DefHist1D(TString, TString,TH1D*, const char*, const char*);
  virtual TH2D *DefHist2D(TString, TString,TH2D*, const char*, const char*);

  //auxiliary; used by DefHist1D and DefHist2D
  virtual double DistSelector1D( double, TString);
//   {return 0;};///3pkt Gauss!!! w skali log lub lin
  virtual double DistSelector2D( double, double, TString);
//   {return 0;};///3pkt Gauss!!! w skali log lub lin
 
  //other helpful routines (to remove???)
  double sqr( const double x ){ return x*x;};
  Double_t Dilogy(double x);
  void StopM(char* message){
    cout <<"++++ TKfig: "<< message << endl; exit(5);}    //Error message
};
/////////////////////////////////////////////////////////////////////////////
#endif
