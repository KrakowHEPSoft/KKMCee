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

#include "KKMC.h"

//class ROBOL : public TNamed {
class ROBOL{
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
  TLorentzVector m_PiPl, m_PiMn;     // pion 4-momenta
  // ============== Histograms follow =================================
  TH1D   *HST_KKMC_NORMA;    // special histo with KKMC normalization & xpar
  //
  TH1D   *hst_weight;
  TH1D   *hst_Mff;
  TH1D   *hst_Q2kloe;
  TH1D   *hst_nPhAll;
  TH1D   *hst_nPhTrg1;
  TH1D   *hst_nPhTrg2;
  TH1D   *hst_thPhAll;
  TH1D   *hst_thPhTrg1;
  TH1D   *hst_Q2All;
  TH1D   *hst_Q2Trig0;
  TH1D   *hst_Q2Trig1;
  TH1D   *hst_Q2Trig2;
  TH1D   *hst_Q2Trig2b;
  TH1D   *hst_vvMuTrg0;
  TH1D   *hst_vvMuTrg1;
  TH1D   *hst_Q2MuTrg0;
  TH1D   *hst_Q2MuTrg1;
  TH1D   *hst_Q2dif1;
  TH1D   *hst_Q2dif2;
  TH1D   *hst_Q2dif3;
  TH1D   *hst_PiCosthTrg1;
  TH1D   *hst_PiCosthTrg2;
  TH1D   *hst_PiPtTrg1;
  TH1D   *hst_PiPtTrg2;
  TH1D   *hst_PivvY;
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
//                      ClassDef(ROBOL,1)
};
////////////////////////////////////////////////////////////////////////////
#endif
