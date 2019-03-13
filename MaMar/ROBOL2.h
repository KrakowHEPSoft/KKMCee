#ifndef ROBOL2_H
#define ROBOL2_H
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                CLASS ROBOL                                                //
//                                                                           //
//     Three-stroke engine for analysis: Initialize-Production-Finalize      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"

#include "KKMC.h"

//class ROBOL : public TNamed {
class ROBOL2{
 public:
  long  m_NevGen;           // event serial number
  long  m_count1;           // auxiliary event counter (debug)
  double m_xpar[10001];      // complete input of KKMC run
  KKMC   *KKMC_generator;    // MC event generator
  // =============== local mirror of KKMC event =======================
  TLorentzVector m_pbea1,m_pbea2;    // initial beams
  TLorentzVector m_pfer1,m_pfer2;    // final fermions
  int           m_Nphot;            // photon multiplicity
  TLorentzVector m_phot[100];        // photon 4-momenta
  int           m_Npart;            // no of particles in Lund/Pythia common
  PartLund m_Event[4000];            // content of /PYJETS/
  // ============== Histograms follow =================================
  TH1D   *HST_KKMC_NORMA;    // special histo with KKMC normalization & xpar
  //
  TH1D   *hst_weight;
  TH1D   *hst_Mff;
  TH1D   *hst_nPhAll;
  TH1D   *hst_nPhVis;
  //
  TH1D   *hst_LnThPhAll;
  TH1D   *hst_LnThPhVis;
  //
  TH1D   *hst_vtNuCeex2;
  //
  TH1D   *hst_vaNuCeex2;
  TH1D   *hst_vaMuCeex2;
  ///
  TH1D   *hst_vPhotNuel;
  TH1D   *hst_vPhotNumu;
  ///
  TH1D   *hst_vtMuCeex2;
  //================================
  TH1D   *hst_vvNuCeex1;
  TH1D   *hst_vvNuCeex2;
  TH1D   *hst_vvNuCeex12;
  //--------------------------------
  TH1D   *hst_vvMuCeex1;
  TH1D   *hst_vvMuCeex2;
  TH1D   *hst_vvMuCeex12;
  TH1D   *hst_vvMuCeex2ifi;
  ///
  TH1D   *hst_vvMuCeex1n;
  TH1D   *hst_vvMuCeex2n;
  TH1D   *hst_vvMuCeex12n;
  ///===============================
  TH1D   *hst_evNuCeex1;
  TH1D   *hst_evNuCeex2;
  TH1D   *hst_evNuCeex12;
  ///-------------------------------
  TH1D   *hst_evMuCeex1;
  TH1D   *hst_evMuCeex2;
  TH1D   *hst_evMuCeex12;
  TH1D   *hst_evMuCeex2ifi;
  ///
  TH1D   *hst_evMuCeex1n;
  TH1D   *hst_evMuCeex2n;
  TH1D   *hst_evMuCeex12n;
  ///==============================
  TH1D   *hst_CosPLCeex2;
  ///
  double m_YSum;   // sum of weights
  double m_YSum2;  // sum of weights^2
  //
 public:
  ROBOL2(){
// my constructor
  }
  ~ROBOL2(){
// my destructor
  }
 public:
// Methods
  TH1D *TH1D_UP(const char*, const char*, int , double, double );
  TH2D *TH2D_UP(const char*, const char*, int , double, double, int, double, double);
  void Initialize(long&);
  void Production(long&);
  void KKMC_NORMA();
  void Finalize();
  void PartImport();
  int PartCount(  const int);
  int PartFindAny(const int);
  int PartFindStable(const int);
  void PyPrint(const int );
  void MomPrint( TLorentzVector&);
  void ReaData(const char*, int, double[]);
////////////////////////////////////////////////////////////////////////////
//                      ClassDef(ROBOL,1)
};
////////////////////////////////////////////////////////////////////////////
#endif
