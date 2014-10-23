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
  TH1D   *hst_vTrueMain;
  TH1D   *hst_vTrueMu;
  TH1D   *hst_vTrueCeex2;
  ///
  TH1D   *hst_vPhotMain;
  TH1D   *hst_vPhotCeex1;
  TH1D   *hst_vPhotCeex2;
  TH1D   *hst_vPhotCeex12;
  ///
  TH1D   *hst_vPhotNuel;
  TH1D   *hst_vPhotNumu;
  ///
  TH1D   *hst_mPhotCeex1;
  TH1D   *hst_mPhotCeex2;
  TH1D   *hst_mPhotCeex12;
  ///
  TH1D   *hst_nPhotCeex1;
  TH1D   *hst_nPhotCeex2;
  TH1D   *hst_nPhotCeex12;
  TH1D   *hst_lPhotCeex1;
  TH1D   *hst_lPhotCeex2;
  TH1D   *hst_lPhotCeex12;
  ///
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
