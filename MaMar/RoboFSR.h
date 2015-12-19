#ifndef ROBOLFSR_H
#define ROBOLFSR_H
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                CLASS RoboFSR                                                //
//                                                                           //
//     Three-stroke engine for analysis: Initialize-Production-Finalize      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"

#include "KKMC.h"

//class RoboFSR : public TNamed {
class RoboFSR{
 public:
  long   m_NevGen;           // event serial number
  long   m_count1;           // auxiliary event counter (debug)
  double m_xpar[10001];      // complete input of KKMC run
  KKMC   *KKMC_generator;    // MC event generator
  // =============== local mirror of KKMC event =======================
  TLorentzVector m_pbea1,m_pbea2;    // initial beams
  TLorentzVector m_pfer1,m_pfer2;    // final fermions
  long           m_Nphot;            // photon multiplicity
  TLorentzVector m_phot[100];        // photon 4-momenta
  long           m_Npart;            // no of particles in Lund/Pythia common
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
  //
  TH1D   *hst_s1Ceex2;
  TH1D   *hst_svk;
  //
  TH1D   *hst_M100mu;
  //
  double m_YSum;   // sum of weights
  double m_YSum2;  // sum of weights^2
  //
 public:
  RoboFSR(){
// my constructor
  }
  ~RoboFSR(){
// my destructor
  }
 public:
// Methods
  void Initialize(long &NevTot);
  void Production(long &iEvent);
  void KKMC_NORMA();
  void Finalize();
  void PartImport();
  long PartCount(  const long);
  long PartFindAny(const long);
  long PartFindStable(const long);
  void PyPrint(const int );
  void MomPrint( TLorentzVector&);
  void ReaData(char DiskFile[], int imax, double xpar[]);
////////////////////////////////////////////////////////////////////////////
//                      ClassDef(RoboFSR,1)
};
////////////////////////////////////////////////////////////////////////////
#endif
