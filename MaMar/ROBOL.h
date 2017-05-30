#ifndef ROBOL_H
#define ROBOL_H
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
class ROBOL{
 public:
  long  m_NevGen;            // event serial number
  long  m_count1;            // auxiliary event counter (debug)
  double m_xpar[10001];      // complete input of KKMC run
  KKMC   *KKMC_generator;    // MC event generator
  // =============== local mirror of KKMC event =======================
  TLorentzVector m_pbea1,m_pbea2;    // initial beams
  TLorentzVector m_pfer1,m_pfer2;    // final fermions
  int           m_Nphot;            // photon multiplicity
  TLorentzVector m_phot[100];        // photon 4-momenta
  int           m_Npart;            // no of particles in Lund/Pythia common
  PartLund m_Event[4000];            // content of /PYJETS/
  TLorentzVector m_pMu1, m_pMu2;     // muon 4-momenta
  // ============== Histograms follow =================================
  TH1D   *HST_KKMC_NORMA;    // special histo with KKMC normalization & xpar
  //
  TH1D   *hst_weight;
  TH1D   *hst_Mff;
  TH1D   *hst_Q2kloe;
  TH1D   *hst_nPhAll;
  TH1D   *hst_nPhVis;
  TH1D   *hst_vTrueMain;
  TH1D   *hst_vTrueCeex2;
  TH1D   *hst_vAlepCeex2;
  TH1D   *hst_vXGenCeex2;
  TH1D   *hst_Cost1Ceex2;
  TH1D   *hst_CosPLCeex2;
  TH1D   *hst_CosPRCeex2;
  TH1D   *hst_CosPREex2;
  //
  TH2D   *sca_vTcPR_Ceex2;
  TH2D   *sca_vTcPR_Ceex2n;
  TH2D   *sca_vTcPR_Eex2;
  TH2D   *sca_vXcPR_Ceex2;
  TH2D   *sca_vXcPR_Eex2;
  //
  TH2D   *sct_vTcPR_Ceex2;    // vvtrue<02
  TH2D   *sct_vTcPR_Ceex2n;   // vvtrue<02
  //
  TH2D   *sct_vAcPR_Ceex2;    // Main CEEX2 KKMC , ISR+FSR+IFI
  TH2D   *sct_vAcPR_Ceex2n;   // IFI  off
  TH2D   *sct_vKcPL_Ceex2;    // vv from Karlud (pure ISR) thetaPL
  TH2D   *sct_vTcPL_Ceex2;    // vv bare muons
  TH2D   *sct_vAcPL_Ceex2;    // Main CEEX2 KKMC , ISR+FSR
    //
  TH1D   *hst_vACeex2;
  TH1D   *hst_vACeex21F;
  TH1D   *hst_vACeex21B;

  //
  double m_YSum;   // sum of weights
  double m_YSum2;  // sum of weights^2
  //
 public:
  ROBOL(){
// my constructor
  }
  ~ROBOL(){
// my destructor
  }
 public:
// Methods
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
