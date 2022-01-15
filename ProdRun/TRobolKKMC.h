#ifndef TRobolKKMC_H
#define TRobolKKMC_H
//                CLASS ROBOL                                                //
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"

/// OUR headers
#include "TRobol.h"
#include "KKee2f.h"

class TRobolKKMC : public TRobol
{
///--- data members
 public:
 long   m_NevGen;             // event serial number
 long   m_count1;             // auxiliary event counter (debug)
//
 TH1D   *hst_WtMain;          //!  No streamer!!!
 TH1D   *hst_WtMain4;         //!  No streamer!!!
 TH1D   *hst_WtMain8;         //!  No streamer!!!
 TH1D   *hst_WtFoam;          //!  No streamer!!!
 TH1D   *hst_WtCeex2n;        //!  No streamer!!!
//
 TH1D   *hst_vvTrue;          //!  No streamer!!!
 TH1D   *hst_nPhot;           //!  No streamer!!!
 TH1D   *hst_vvBES;           //!  No streamer!!!
 TH1D   *hst_CosTheta;        //!  No streamer!!!
 TH1D   *hst_CosThOve;        //!  No streamer!!!
 //
 TH2D   *sca_vTcPR_Ceex2;     //!  No streamer!!!
 TH2D   *sca_vTcPR_Ceex2n;    //!  No streamer!!!
 TH2D   *sca_vTcPR_Eex0;      //!  No streamer!!!
 TH2D   *sca_vTcPR_Eex2;      //!  No streamer!!!
//
 TH1D   *hst_LnThPhAll;        //!  No streamer!!!
 TH1D   *hst_LnThPhVis;        //!  No streamer!!!
 TH1D   *hst_vtNuCeex2;        //!  No streamer!!!
 TH1D   *hst_vaNuCeex2;        //!  No streamer!!!

 TH1D   *hst_vPhotNuel;        //!  No streamer!!!
 TH1D   *hst_vPhotNumu;        //!  No streamer!!!
 //================================
 TH1D   *hst_vvNuCeex1;        //!  No streamer!!!
 TH1D   *hst_vvNuCeex2;        //!  No streamer!!!
 TH1D   *hst_vvNuCeex12;       //!  No streamer!!!
 //
 TH1D   *hst_nPhAll;           //!  No streamer!!!
 TH1D   *hst_nPhVis;           //!  No streamer!!!

 TH1D   *hst_vaNuElCeex2;           //!  No streamer!!!
 TH1D   *hst_vaNuMuCeex2;           //!  No streamer!!!

 //
 TH2D   *sca_r1r2;            //!  No streamer!!!
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
 void Finalize();
////////////////////////////////////////////////////////////////////////////
                        ClassDef(TRobolKKMC,1)
};
////////////////////////////////////////////////////////////////////////////
#endif
