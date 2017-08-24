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
/// ============== Histograms follow =================================
  TH1D   *hst_weight1;           //! no streamer
  TH1D   *hst_weight3;           //! no streamer
  TH1D   *hst_weight5;           //! no streamer

  TH1D   *HST_xx_Ord1;          //! no streamer
  TH1D   *HST_xx_Crd1;          //! no streamer

  TH1D   *HST_xx_Ceex2;          //! no streamer
  TH1D   *HST_xx_Ceex2n;         //! no streamer

  TH1D   *HST_tx_Ceex2n;         //! no streamer

  TH2D   *SCA_xc_Ceex2;          //! no streamer
  TH2D   *SCA_xc_Ceex2n;         //! no streamer

  TH2D   *SCT_xc_Ceex2;          //! no streamer
  TH2D   *SCT_xc_Ceex2n;         //! no streamer
  TH2D   *SCT_xc_EEX2;           //!

  TH1D   *HST_FOAM_NORMA3;       //! no streamer
  TH1D   *HST_FOAM_NORMA1;       //! no streamer
////////////////////////////////////////////////////////////////////////////
  ClassDef(TRobolFOAM,2); // Monte Carlo generator
};
/////////////////////////////////////////////////////////////////////////////
//                End of the class TRobolFOAM                                 //
/////////////////////////////////////////////////////////////////////////////
#endif
