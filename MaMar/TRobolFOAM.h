#ifndef TRobolFOAM_H
#define TRobolFOAM_H
///////////////////////////////////////////////////////////////////////////////
///                CLASS TRobol2pB
///     Three-stroke engine for analysis: Initialize-Production-Finalize
///     Daugter class specific to simple multidimensional integration

/// OUR headers
#include "TRobol.h"
#include "TMCgenFOAM.h"

class TRobolFOAM : public TRobol {
/// member functions
  public:
  TRobolFOAM();                // explicit default constructor for streamer
  TRobolFOAM(const char*);     // user constructor
  virtual ~TRobolFOAM();               // explicit destructor
///---
  virtual void Initialize(ofstream*, TFile*, TFile*);
  virtual void Hbooker();
  virtual void Production(double &);
  virtual void Finalize();
/// data members
  double    m_xmin;          // dummy
  double    m_xmax;          // dummy
  double    m_vminL;          // dummy
/// ============== Histograms follow =================================
  TH1D   *hst_weight1;           //! no streamer
  TH1D   *hst_weight3;           //! no streamer
  TH1D   *hst_weight5;           //! no streamer

  TH1D   *HST_xx_Ord1n;          //! no streamer
  TH1D   *HST_xx_Crd1n;          //! no streamer
  TH1D   *HST_xx_Ord1;           //! no streamer
  TH1D   *HST_xx_Crd1;           //! no streamer
  TH1D   *HST_xx_Hrd1;           //! no streamer
  TH1D   *HST_xx_Srd1;           //! no streamer
  TH1D   *HST_xx_Ird1;           //! no streamer

  TH1D   *HST_xx_Ceex2;          //! no streamer
  TH1D   *HST_xx_Ceex2n;         //! no streamer

  TH1D   *HST_tx_Ceex2n;         //! no streamer

  TH2D   *SCA_xc_Ceex2;          //! no streamer
  TH2D   *SCA_xc_Ceex2n;         //! no streamer

  TH2D   *SCA_xc_Ceex0;          //! no streamer
  TH2D   *SCA_xc_Ceex0n;         //! no streamer

  TH2D   *SCT_xc_Ceex2;          //! no streamer
  TH2D   *SCT_xc_Ceex2n;         //! no streamer

  TH2D   *SCT_xc_Ceex0;          //! no streamer
  TH2D   *SCT_xc_Ceex0n;         //! no streamer

  TH2D   *SCT_xc_EEX2;           //! no streamer
  TH2D   *SCT_xc_EEX0;           //! no streamer

  TH2D   *SCT_xc_EEX2i;          //! no streamer
  TH2D   *SCT_xc_EEX2n;          //! no streamer

  TH2D   *SCN_xc_EEX2;           //! no streamer
  TH2D   *SCN_xc_EEX0;           //! no streamer

  TH2D   *SCN_xc_EEX2i;          //! no streamer
  TH2D   *SCN_xc_EEX2n;          //! no streamer
/////////////////////////////////////////////////////////
  TH2D   *SCT_Lxc_EEX2i;         //! no streamer
  TH2D   *SCT_Lxc_EEX2n;         //! no streamer
/////////////////////////////////////////////////////////
  TH1D   *HST5_xx_Ceex2;         //! no streamer
  TH1D   *HST5_xc_Ceex2;         //! no streamer
  TH1D   *HST5_xx_forw_Ceex2;    //! no streamer

  TH1D   *HST5_xx9_Ceex2;        //! no streamer
  TH1D   *HST5_xx9_forw_Ceex2;   //! no streamer
/////////////////////////////////////////////////////////
  TH1D   *HST_cc_EEX2_vmax02;    //! no streamer
  TH1D   *HST_cc_EEX2_vmax002;   //! no streamer
  TH1D   *HST_cc_EEX2_vmax0002;  //! no streamer

  TH1D   *HST_cc_ceex2_vmax02;    //! no streamer
  TH1D   *HST_cc_ceex2_vmax002;   //! no streamer
  TH1D   *HST_cc_ceex2_vmax0002;  //! no streamer

  TH1D   *HST_cs_EEX2_vmax02;    //! no streamer
  TH1D   *HST_cs_EEX2_vmax002;   //! no streamer
  TH1D   *HST_cs_EEX2_vmax0002;  //! no streamer

  TH1D   *HST_cs_ceex2n_vmax02;    //! no streamer
  TH1D   *HST_cs_ceex2n_vmax002;   //! no streamer
  TH1D   *HST_cs_ceex2n_vmax0002;  //! no streamer

  TH1D   *HST_cs_ceex2_vmax02;    //! no streamer
  TH1D   *HST_cs_ceex2_vmax002;   //! no streamer
  TH1D   *HST_cs_ceex2_vmax0002;  //! no streamer

  TH1D   *HST_FOAM_NORMA3;       //! no streamer
  TH1D   *HST_FOAM_NORMA3i;       //! no streamer
  TH1D   *HST_FOAM_NORMA1;       //! no streamer
  TH1D   *HST_FOAM_NORMA2;       //! no streamer
////////////////////////////////////////////////////////////////////////////
  ClassDef(TRobolFOAM,2); // Monte Carlo generator
};
/////////////////////////////////////////////////////////////////////////////
//                End of the class TRobolFOAM                                 //
/////////////////////////////////////////////////////////////////////////////
#endif
