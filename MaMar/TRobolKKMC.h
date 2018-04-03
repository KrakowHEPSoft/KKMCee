#ifndef TRobolKKMC_H
#define TRobolKKMC_H
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

/// OUR headers
#include "TRobol.h"
#include "TMCgenKKMC.h"

class TRobolKKMC : public TRobol
{
///--- data members
 public:
  long   m_NevGen;             // event serial number
  long   m_count1;             // auxiliary event counter (debug)
  double m_xpar[10001];        // complete input of KKMC run
//  KKMC   *KKMC_generator;    // goes to base class
  // =============== local mirror of KKMC event =======================
  TLorentzVector m_pbea1;      //! initial beams
  TLorentzVector m_pbea2;      //! initial beams
  TLorentzVector m_pfer1;      //! final fermions
  TLorentzVector m_pfer2;      //! final fermions
//  int           m_Nphot;      //! photon multiplicity
  int            m_Nphot;      //! photon multiplicity
  TLorentzVector m_phot[100];  //! photon 4-momenta
  int           m_Npart;       //! no of particles in Lund/Pythia common
  //TPartLund m_Event[4000];     //! content of /PYJETS/
  TLorentzVector m_pMu1;       //! muon 4-momenta
  TLorentzVector m_pMu2;       //! muon 4-momenta
  // ============== Histograms follow =================================
  //TH1D   *HST_KKMC_NORMA;    // goes to base class
  //
  TH1D   *hst_weight;          //!  No sstreamer!!!
//  TH1D   *hst_Mff;             //!  No sstreamer!!!
  TH1D   *hst_Q2kloe;          //!  No sstreamer!!!
  TH1D   *hst_nPhAll;          //!  No sstreamer!!!
  TH1D   *hst_nPhVis;          //!  No sstreamer!!!
  TH1D   *hst_vTrueMain;       //!  No sstreamer!!!
  TH1D   *hst_vTrueCeex2;      //!
  TH1D   *hst_vAlepCeex2;      //!
  TH1D   *hst_vXGenCeex2;      //!
  TH1D   *hst_Cost1Ceex2;      //!
  TH1D   *hst_CosPLCeex2;      //!
  TH1D   *hst_CosPRCeex2;      //!
  TH1D   *hst_CosPREex2;       //!
  //
  TH2D   *sca_vTcPR_Ceex2;     //!  No sstreamer!!!
  TH2D   *sca_vTcPR_Ceex2n;    //!  No sstreamer!!!
  TH2D   *sca_vTcPR_Eex2;      //!

  TH2D   *sca_vXcPR_Ceex2;     //!
  TH2D   *sca_vXcPR_Eex2;      //!

  TH2D   *sca_vTcPL_Ceex0;     //!  No sstreamer!!!
  TH2D   *sca_vTcPL_Ceex0n;    //!  No sstreamer!!!
  TH2D   *sca_vTcPL_Ceex2;     //!  No sstreamer!!!
  TH2D   *sca_vTcPL_Ceex2n;    //!  No sstreamer!!!
  TH2D   *sca_vTcPL_Eex2;      //!  No sstreamer!!!

  TH2D   *sca_vTvA_Eex2;       //!  No sstreamer!!!
  TH2D   *sca_vKvA_Eex2;       //!  No sstreamer!!!
  /////////////////////////////////////////////////////////////////////
  TH2D   *sct_vTcPR_Ceex2;    //! vvtrue<02 No sstreamer!!!
  TH2D   *sct_vTcPR_Ceex2n;   //! vvtrue<02 No sstreamer!!!
  TH2D   *sct_vTcPR_EEX2;     //! vvtrue<02 No sstreamer!!!
  //
  TH2D   *sct_vAcPR_Ceex2;    //! Main CEEX2 KKMC , ISR+FSR+IFI
  TH2D   *sct_vAcPR_Ceex2n;   //! IFI  off
  TH2D   *sct_vKcPL_Ceex2;    //! vv from Karlud (pure ISR) thetaPL

  TH2D   *sct_vTcPL_Ceex2;    //! vv bare muons
  TH2D   *sct_vTcPL_Ceex2n;   //! vv bare muons
  TH2D   *sct_vTcPL_Ceex1;    //! vv bare muons
  TH2D   *sct_vTcPL_Ceex1n;   //! vv bare muons
  TH2D   *sct_vTcPL_Ceex0;    //! vv bare muons
  TH2D   *sct_vTcPL_Ceex0n;   //! vv bare muons

  TH2D   *sct_vAcPL_Ceex2;    //! Main CEEX2 KKMC , ISR+FSR
  /////////////// NEW!!!!
  // CEEX series
  TH1D   *hst_vT_Ceex1;       //!  No streamer!!!
  TH1D   *hst_vT_Ceex2;       //!  No streamer!!!
  TH1D   *hst_vT_Ceex21;      //!  No streamer!!!
  TH1D   *hst_vT_Ceex1_F;     //!  No streamer!!!
  TH1D   *hst_vT_Ceex2_F;     //!  No streamer!!!
  TH1D   *hst_vT_Ceex21_F;    //!  No streamer!!!
  // CEEXn series
  TH1D   *hst_vT_Ceex1n;       //!  No streamer!!!
  TH1D   *hst_vT_Ceex2n;       //!  No streamer!!!
  TH1D   *hst_vT_Ceex21n;      //!  No streamer!!!
  TH1D   *hst_vT_Ceex1n_F;     //!  No streamer!!!
  TH1D   *hst_vT_Ceex2n_F;     //!  No streamer!!!
  TH1D   *hst_vT_Ceex21n_F;    //!  No streamer!!!
  // EEX series
  TH1D   *hst_vT_EEX1;        //!  No streamer!!!
  TH1D   *hst_vT_EEX2;        //!  No streamer!!!
  TH1D   *hst_vT_EEX3;        //!  No streamer!!!
  TH1D   *hst_vT_EEX21;       //!  No streamer!!!
  TH1D   *hst_vT_EEX32;       //!  No streamer!!!
  TH1D   *hst_vT_EEX1_F;      //!  No streamer!!!
  TH1D   *hst_vT_EEX2_F;      //!  No streamer!!!
  TH1D   *hst_vT_EEX3_F;      //!  No streamer!!!
  TH1D   *hst_vT_EEX21_F;     //!  No streamer!!!
  TH1D   *hst_vT_EEX32_F;     //!  No streamer!!!
  // Older quasi-realistic
  TH1D   *hst_vACeex1;        //!  No streamer!!!
  TH1D   *hst_vACeex2;        //!  No streamer!!!
  TH1D   *hst_vACeex1F;       //!  No streamer!!!
  TH1D   *hst_vACeex2F;       //!  No streamer!!!
  TH1D   *hst_vACeex21;       //!  No streamer!!!
  TH1D   *hst_vACeex21F;      //!  No streamer!!!
  TH1D   *hst_vACeex21B;      //!  No streamer!!!
 //]]]]]]]]]]]]
  //
  TH1D   *hst_vT_Ceex2_cPLr90;  //!  No streamer!!!
  //
  TH1D   *hst_vT_Ceex2_xcPL;    //!  No streamer!!!
  TH1D   *hst_vT_Ceex2n_xcPL;   //!  No streamer!!!
  //
  TH1D   *hst_vT_Ceex2_cPL_forw;      //!  No streamer!!!
  TH1D   *hst_vT_Ceex2_cPLr90_forw;  //!  No streamer!!!
//
  TH1D   *hst_Alf0CutB_Ceex2;         //!  No streamer!!!
  TH1D   *hst_Alf0CutB_Ceex2F;        //!  No streamer!!!
  TH1D   *hst_Alf2CutB_Ceex2;         //!  No streamer!!!
  TH1D   *hst_Alf2CutB_Ceex2F;        //!  No streamer!!!
  TH1D   *hst_Al20CutB_Ceex2;         //!  No streamer!!!
  TH1D   *hst_Al20CutB_Ceex2F;        //!  No streamer!!!
  //
  TH1D   *hst_Alf0CutA_Ceex2;         //!  No streamer!!!
  TH1D   *hst_Alf0CutA_Ceex2F;        //!  No streamer!!!
  TH1D   *hst_Alf2CutA_Ceex2;         //!  No streamer!!!
  TH1D   *hst_Alf2CutA_Ceex2F;        //!  No streamer!!!
  TH1D   *hst_Al20CutA_Ceex2;         //!  No streamer!!!
  TH1D   *hst_Al20CutA_Ceex2F;        //!  No streamer!!!
  //
  TH1D   *hst_Alf2CutB_Ceex2n;         //!  No streamer!!!
  TH1D   *hst_Alf2CutB_Ceex2nF;        //!  No streamer!!!
  TH1D   *hst_Al20CutB_Ceex2n;         //!  No streamer!!!
  TH1D   *hst_Al20CutB_Ceex2nF;        //!  No streamer!!!
  //
  TH1D   *hst_Alf2CutA_Ceex2n;         //!  No streamer!!!
  TH1D   *hst_Alf2CutA_Ceex2nF;        //!  No streamer!!!
  TH1D   *hst_Al20CutA_Ceex2n;         //!  No streamer!!!
  TH1D   *hst_Al20CutA_Ceex2nF;        //!  No streamer!!!
  //
  TH1D   *hst_vA_Ceex1;           //!  No streamer!!!
  TH1D   *hst_vA_Ceex2;           //!  No streamer!!!
  TH1D   *hst_vA_Ceex21;          //!  No streamer!!!
  TH1D   *hst_vA_Ceex1_F;         //!  No streamer!!!
  TH1D   *hst_vA_Ceex2_F;         //!  No streamer!!!
  TH1D   *hst_vA_Ceex21_F;        //!  No streamer!!!

  TH1D   *hst_vB_Ceex1;          //!  No streamer!!!
  TH1D   *hst_vB_Ceex2;          //!  No streamer!!!
  TH1D   *hst_vB_Ceex21;         //!  No streamer!!!
  TH1D   *hst_vB_Ceex1_F;         //!  No streamer!!!
  TH1D   *hst_vB_Ceex2_F;         //!  No streamer!!!
  TH1D   *hst_vB_Ceex21_F;        //!  No streamer!!!

  TH1D   *hst_vA_Ceex2i;          //!  No streamer!!!
  TH1D   *hst_vA_Ceex2i_F;        //!  No streamer!!!
  TH1D   *hst_vB_Ceex2i;          //!  No streamer!!!
  TH1D   *hst_vB_Ceex2i_F;        //!  No streamer!!!

  double m_YSum;   // sum of weights
  double m_YSum2;  // sum of weights^2
  double m_DelAlf; // range of AlphaQED (relative)
  double m_vvcut2; // photon energy cut
///////////////////////////////////////////
/// mandatory constructors and destructors
  public:
  TRobolKKMC();                // explicit default constructor for streamer
  TRobolKKMC(const char*);     // user constructor
  virtual ~TRobolKKMC();       // explicit destructor
/// mandatory methods
  virtual void Initialize(ofstream*, TFile*, TFile*);
  virtual void Hbooker();
  virtual void Production(double &);
/*
 public:
  TRobolKKMC(){
// my constructor
  }
  ~TRobolKKMC(){
// my destructor
  }
 public:
// Methods
  void Initialize(int &NevTot);
  void Production(int &iEvent);
  void KKMC_NORMA();
*/
//////////////////////////////////////////////
// Other user methods
  void Finalize();
//  void PartImport();
//  int PartCount(  const int);
//  int PartFindAny(const int);
//  int PartFindStable(const int);
//  void PyPrint(const int );
  void MomPrint( TLorentzVector&);
////////////////////////////////////////////////////////////////////////////
                      ClassDef(TRobolKKMC,1)
};
////////////////////////////////////////////////////////////////////////////
#endif
