#ifndef TRobolFOAM_H
#define TRobolFOAM_H
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                CLASS TRobolFOAM                                             //
//                                                                           //
//     Three-stroke engine for analysis: Initialize-Production-Finalize      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
using namespace std;

#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"

//#include "KKMC.h"

#include "TRobol.h"

#include "TMCgenFOAM.h"

class TRobolFOAM: public TRobol{
// constructors destructors
public:
  TRobolFOAM();
  TRobolFOAM(const char*);     // user constructor
  ~TRobolFOAM();
 public:
  long     m_NevGen;           // event serial number
  long     m_count1;           // auxiliary event counter (debug)
  double   m_xpar[10001];      // complete input of KKMC run
//  KKMC     *KKMC_generator;    // MC event generator
  TMCgenFOAM *LibSem;            // Distributions for Foam
  //TFoam    *MC_Gen3;
  //TFoam    *MC_Gen5;
  TFOAM    *MC_Gen3;
  TFOAM    *MC_Gen5;
  TRandom  *PseRan;

  double   m_Xsav3;
  double   m_Xsav5;
  // ============== Histograms follow =================================
  TH1D   *HST_FOAM_NORMA3;    // special histo with KKMC normalization & xpar
  TH1D   *HST_FOAM_NORMA5;    // special histo with KKMC normalization & xpar
  //
  TH1D   *hst_weight3;
  TH1D   *hst_weight5;

  TH1D   *HST_xx_Ceex2n;
  TH1D   *HST_xx_Ceex2;
  //
  TH2D   *SCA_xc_Ceex2;    // Main CEEX2 KKMC , ISR+FSR
  TH2D   *SCA_xc_Ceex2n;   // IFI  off
  //
  TH2D   *SCT_xc_Ceex2;    // Rescricted v-range
  TH2D   *SCT_xc_Ceex2n;
//-------------------------
 public:
// Methods
  void Initialize(long &NevTot);
  void Production(long &iEvent);
  void KKMC_NORMA();
  void Finalize();
  void ReaData(char DiskFile[], int imax, double xpar[]);
////////////////////////////////////////////////////////////////////////////
                      ClassDef(TRobolFOAM,2)
  ////////////////////////////////////////////////////////////////////////////
};
////////////////////////////////////////////////////////////////////////////
#endif
