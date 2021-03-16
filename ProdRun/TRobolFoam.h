#ifndef TRobolFoam_H
#define TRobolFoam_H
//                CLASS ROBOL                                                //
#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"

/// OUR headers
#include "TRobol.h"
#include "KKeeFoam.h"

class TRobolFoam : public TRobol
{
///--- data members
 public:
  long   m_NevGen;             // event serial number
  long   m_count1;             // auxiliary event counter (debug)
  double m_xpar[10001];        // complete input of KKMC run
  double m_Mllmin;             //
  double m_Mllmax;             //
// ========== Weight monitoring ===========
//  THwtMon *mon_WtMain;         //!  No streamer!!!
//  THwtMon *mon_WtFoam;         //!  No streamer!!!
////////////////////////////////////////////////////////////////////////////
  TH1D   *hst_weight;          //!  No streamer!!!
  TH1D   *hst6_weight;         //!  No streamer!!!
  //
 TH1D   *Hst_Mll;             //!  No streamer!!!
  TH1D   *Hst_Mll_eex0;        //!  No streamer!!!
  TH1D   *Hst_MllF_eex0;       //!  No streamer!!!
  TH1D   *Hst_Mll_eex2;        //!  No streamer!!!
  TH1D   *Hst_MllF_eex2;       //!  No streamer!!!
  //
  TH1D   *Hst_Mll_ceex2;       //!  No streamer!!!
  TH1D   *Hst_MllF_ceex2;      //!  No streamer!!!
  TH1D   *hst6_Mll_ceex2;       //!  No streamer!!!
  TH1D   *hst6_MllF_ceex2;      //!  No streamer!!!
  //
  TH1D   *Hst_Mll_ceex0;       //!  No streamer!!!
  TH1D   *Hst_MllF_ceex0;      //!  No streamer!!!
  TH1D   *hst6_Mll_ceex0;       //!  No streamer!!!
  TH1D   *hst6_MllF_ceex0;      //!  No streamer!!!
  /////////////////////////////////////////////////////////////////////////////
  /// mandatory constructors and destructors
      public:
      TRobolFoam();                // explicit default constructor for streamer
      TRobolFoam(const char*);     // user constructor
      virtual ~TRobolFoam();       // explicit destructor
  /// mandatory methods
      virtual void Initialize(ofstream*, TFile*, TFile*);
      virtual void Hbooker();
      virtual void Production(double &);
  //////////////////////////////////////////////////////////////////////////////
  // Other user methods
      void Finalize();
  //    void MomPrint( TLorentzVector&);
    ////////////////////////////////////////////////////////////////////////////
                          ClassDef(TRobolFoam,1)
    };
    ////////////////////////////////////////////////////////////////////////////
    #endif
