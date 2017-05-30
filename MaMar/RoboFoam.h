#ifndef ROBOFOAM_H
#define ROBOFOAM_H
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                CLASS RoboFoam                                             //
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

#include "KKMC.h"

#include "KKfoam.h"

class RoboFoam{
 public:
  long     m_NevGen;           // event serial number
  long     m_count1;           // auxiliary event counter (debug)
  double   m_xpar[10001];      // complete input of KKMC run
  KKMC     *KKMC_generator;    // MC event generator
  KKfoam   *LibSem;            // Distributions for Foam
  TFoam    *MC_Gen3;
  TFoam    *MC_Gen5;
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
  //
 public:
  RoboFoam(){
// my constructor
  }
  ~RoboFoam(){
// my destructor
  }
//-------------------------
 public:
// Methods
  void Initialize(long&);
  void Production(long&);
  void KKMC_NORMA();
  void Finalize();
  void ReaData(const char*, int, double[]);
////////////////////////////////////////////////////////////////////////////
//                      ClassDef(RoboFoam,1)
};
////////////////////////////////////////////////////////////////////////////
#endif
