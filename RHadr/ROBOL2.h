#ifndef ROBOL2_H
#define ROBOL2_H
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                CLASS ROBOL2                                               //
//                                                                           //
//     Three-stroke engine for analysis: Initialize-Production-Finalize      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "KKMC.h"

//class ROBOL2 : public TNamed {
class ROBOL2{
 public:
  KKMC   *KKMC_generator;    // MC event generator
  TFile  *m_DiskFile;        // DiskFile with all histos
  long   m_NevGen;           // event serial number
  long   m_iBook;            // booking status
  long   m_count1;           // auxiliary event counter (debug)
  double m_xpar[10001];      // complete input of KKMC run
  // =============== local mirror of KKMC event =======================
  TLorentzVector m_pbea1,m_pbea2;    // initial beams
  TLorentzVector m_pfer1,m_pfer2;    // final fermions
  long           m_Nphot;            // photon multiplicity
  TLorentzVector m_phot[100];        // photon 4-momenta
  long           m_Npart;            // no of particles in Lund/Pythia common
  PartLund m_Event[4000];            // content of /PYJETS/
  TLorentzVector m_PiPl, m_PiMn;     // pion 4-momenta
  // ============== Histograms follow =================================
  //
  TH1D   *Hst_Q2hadA;
  TH1D   *Hst_Q2piA;
  TH1D   *Hst_Q2piB;
  TH1D   *Hst_Q2muA;
  TH1D   *Hst_piCosA;
  TH1D   *Hst_piCosB;
  TH1D   *Hst_phCosA;
  TH1D   *Hst_phCosB;
  TH1D   *Hst_Vmu01;
  TH1D   *Hst_Vmu71;
  TH1D   *Hst_Vmu72;
  TH1D   *Hst_Vmu73;
  TH1D   *Hst_Vmu74;
  TH1D   *Hst_Vmu10;
  TH1D   *Hst_Vmu11;
  TH1D   *Hst_Vmu20;
  TH1D   *Hst_Vmu21;
  TH1D   *Hst_Vmu22;
  TH1D   *Hst_Vmu202;
  TH1D   *Hst_Vmu203;
  TH1D   *Hst_Vmu252;
  TH1D   *Hst_Vmu253;
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
  TH1D *HistoUP(const char*, const char*, int, double, double);
  void Initialize(KKMC*, TFile*,long &NevTot, int);
  void Production(long &iEvent);
  void Finalize();
  void PartImport();
  long PartCount(  const long);
  long PartFindAny(const long);
  long PartFindStable(const long);
  void PyPrint(const int );
  void MomPrint( TLorentzVector&);
  void ReaData(char DiskFile[], int imax, double xpar[]);
////////////////////////////////////////////////////////////////////////////
//                      ClassDef(ROBOL2,1)
};
////////////////////////////////////////////////////////////////////////////
#endif
