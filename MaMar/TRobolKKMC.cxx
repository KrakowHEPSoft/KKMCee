///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                CLASS ROBOL                                                //
//                                                                           //
//     Three-stroke engine for analysis: Initialize-Production-Finalize      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TRobolKKMC.h"

# define sw2 setprecision(10) << setw(18)

ClassImp(TRobolKKMC);

///////////////////////////////////////////////////////////////////////////////
//      *************** temporary entries from KKMC ****************
//      SUBROUTINE KarLud_GetVVxx(vv,x1,x2)
extern "C" void  karlud_getvvxx_(double&, double&, double&);
//extern "C" void  pyhepc_(long&);
//extern "C" void  photos_(long&);
//extern "C" void  phoini_();
//extern "C" void  hepevt_setphotosflagtrue_(long&);
//extern "C" void  hepevt_getnhep_(long&);
///////////////////////////////////////////////////////////////////////////////



TRobolKKMC::TRobolKKMC():
  TRobol()
{
  /// This constructor is for ROOT streamers ONLY
  cout<< "@@@@> TRobolKKMC DEFAULT Constructor (for ROOT only) "<<endl;
}


///_____________________________________________________________
TRobolKKMC::TRobolKKMC(const char* Name):
  TRobol(Name)
{
//! Constructor to be used by the user!!!
//! Its important role is to define ALL DEFAULTS.
//! to changed by the user before calling TMCgen::Initialize
  cout<< "@@@@> TRobolKKMC::TRobolFOAM USER Constructor "<<endl;
  m_NevGen=0;
  m_count1=0;

}///TRobolFOAM


///______________________________________________________________________________________
TRobolKKMC::~TRobolKKMC()
{
  //!Explicit destructor
  cout<< "@@@@> TRobolKKMC::TRobolFOAM !!!! DESTRUCTOR !!!! "<<endl;
}///destructor


//______________________________________________________________________________
//////////////////////////////////////////////////////////////
//   Initialize MC generator and analysis programs          //
//////////////////////////////////////////////////////////////
void TRobolKKMC::Initialize(
        ofstream *OutFile, /// Central log-file for messages
        TFile *GenFile,    /// ROOT disk file for CRNG and MC gen.
        TFile *HstFile)    /// ROOT disk file for histograms
{
  cout<< "****> TRobolKKMC::Initialize starts"<<endl;
  TRobol::Initialize(OutFile,GenFile,HstFile);
  ///
  /// Book histograms or read them from the disk
  Hbooker();

//  Mooved to TMCgenKKMC::Initialize
//  dumping KKMC info into normalization histo
/*
  TH1D *HST_KKMC_NORMA = (TH1D*)HstFile->Get("HST_KKMC_NORMA");
  int jmax = ((TMCgenKKMC*)f_MCgen)->m_jmax;
  double CMSene = m_xpar[1];
  cout<<"TRobolKKMC:: CMSene="<<CMSene<<endl; // just for control
  for(int j=1; j<=jmax; j++) {
    HST_KKMC_NORMA->SetBinContent(j,  ((TMCgenFOAM*)f_MCgen)->m_xpar[j]  );    // xpar encoded
  }
  const int jmax =10000;
  ReaData("../../.KK2f_defaults", jmax, m_xpar);  // numbering as in input!!!
  ReaData("./pro.input",         -jmax, m_xpar);  // jmax<0 means no-zeroing
  double ypar[jmax];
  for(int j=0;j<jmax;j++) ypar[j]=m_xpar[j+1];    // ypar has c++ numbering

  NevTot = (long)m_xpar[0];                       // NevTot hidden in xpar[0] !!!
  m_NevGen=0;
  m_count1=0;
  KKMC_generator = new TMCgenKKMC(); // mooved to Start.C
  KKMC_generator->Initialize(ypar);  // done by base class
  cout<<"TRobolKKMC::Initialize:  NevTot = "<<NevTot<<endl;
*/
  cout<< "****> TRobolKKMC::Initialize: finished"<<endl;
}///Initialize

///////////////////////////////////////////////////////////////////////////////
void TRobolKKMC::Hbooker()
{
  ///
  cout<< "****> TRobolFOAM::Hbooker: histogram booking STARTS"<<endl;
  BXOPE(*f_Out);
  BXTXT(*f_Out,"========================================");
  BXTXT(*f_Out,"======    TRobol::Hbooker    ===========");
  BXTXT(*f_Out,"========================================");
  BXCLO(*f_Out);
  f_HstFile->cd();
  //  ************* user histograms  *************
  double CMSene = m_xpar[1];
  int nbin=1000;
  hst_weight  = TH1D_UP("hst_weight" ,  "MC weight",      100, 0.000 , 2.0);
  //hst_Mff     = TH1D_UP("hst_Mff"    ,  "Mass(f-fbar)",  nbin, 0.000 ,CMSene);
  hst_nPhAll  = TH1D_UP("hst_nPhAll" , "No. of photons, all",   8, -0.5 ,7.5);
  hst_nPhVis  = TH1D_UP("hst_nPhVis" , "No. photons, E>10MeV",  8, -0.5 ,7.5);
  int nbv = 50;
  hst_vTrueMain    = TH1D_UP("hst_vTrueMain",   "dSig/dvTrue ", nbv, 0.000 ,1.000);
  hst_vTrueCeex2   = TH1D_UP("hst_vTrueCeex2",  "dSig/dvTrue ", nbv, 0.000 ,1.000);
  hst_vXGenCeex2   = TH1D_UP("hst_vXGenCeex2",  "dSig/dvTrue ", nbv, 0.000 ,1.000);
  hst_vAlepCeex2   = TH1D_UP("hst_vAlepCeex2",  "dSig/dvTrue ", nbv, 0.000 ,1.000);
  int nbc =50;
  hst_Cost1Ceex2= TH1D_UP("hst_Cost1Ceex2",  "dSig/cThet1   ", nbc, -1.000 ,1.000);
  hst_CosPLCeex2= TH1D_UP("hst_CosPLCeex2",  "dSig/cThetPL  ", nbc, -1.000 ,1.000);
  hst_CosPRCeex2= TH1D_UP("hst_CosPRCeex2",  "dSig/cThetPRD ", nbc, -1.000 ,1.000);
  hst_CosPREex2 = TH1D_UP("hst_CosPREex2",   "dSig/cThetPRD ", nbc, -1.000 ,1.000);
  // scatergrams unrestricted vmax<1.0
  sca_vTcPL_Ceex0  = TH2D_UP("sca_vTcPL_Ceex0",    "dSig/dc/dv ", nbv, 0.0 ,1.0, nbc, -1.0 ,1.0);
  sca_vTcPL_Ceex0n = TH2D_UP("sca_vTcPL_Ceex0n",   "dSig/dc/dv ", nbv, 0.0 ,1.0, nbc, -1.0 ,1.0);
  sca_vTcPL_Ceex2  = TH2D_UP("sca_vTcPL_Ceex2",    "dSig/dc/dv ", nbv, 0.0 ,1.0, nbc, -1.0 ,1.0);
  sca_vTcPL_Ceex2n = TH2D_UP("sca_vTcPL_Ceex2n",   "dSig/dc/dv ", nbv, 0.0 ,1.0, nbc, -1.0 ,1.0);
  sca_vTcPL_Eex2   = TH2D_UP("sca_vTcPL_Eex2",     "dSig/dc/dv ", nbv, 0.0 ,1.0, nbc, -1.0 ,1.0);
  sca_vTcPR_Ceex2  = TH2D_UP("sca_vTcPR_Ceex2",  "dSig/dc/dv ", nbv, 0.0 ,1.0, nbc, -1.0 ,1.0);
  sca_vTcPR_Ceex2n = TH2D_UP("sca_vTcPR_Ceex2n",  "dSig/dc/dv ", nbv, 0.0 ,1.0, nbc, -1.0 ,1.0);
  sca_vTcPR_Eex2   = TH2D_UP("sca_vTcPR_Eex2",   "dSig/dc/dv ", nbv, 0.0 ,1.0, nbc, -1.0 ,1.0);
  sca_vXcPR_Ceex2  = TH2D_UP("sca_vXcPR_Ceex2",  "dSig/dc/dv ", nbv, 0.0 ,1.0, nbc, -1.0 ,1.0);
  sca_vXcPR_Eex2   = TH2D_UP("sca_vXcPR_Eex2",   "dSig/dc/dv ", nbv, 0.0 ,1.0, nbc, -1.0 ,1.0);
  // Miscelaneous
  sca_vTvA_Eex2    = TH2D_UP("sca_vTvA_Eex2",   "dSig/dvT/dvA ", 100, 0.0 ,0.2, 100, 0.0 ,0.2);
  sca_vKvA_Eex2    = TH2D_UP("sca_vKvA_Eex2",   "dSig/dvK/dvA ", 100, 0.0 ,0.2, 100, 0.0 ,0.2);
  ///////////////////////////////////////////////////////////////////////////
  //  New bigger scatergrams, restricted vmax<0.2
  int NBv =100; int NBc = 100;
  double vmx2= 0.20;
  sct_vTcPR_Ceex2 = TH2D_UP("sct_vTcPR_Ceex2",  "dSig/dc/dv ", NBv, 0.0 ,vmx2, NBc, -1.0 ,1.0);
  sct_vTcPR_Ceex2n= TH2D_UP("sct_vTcPR_Ceex2n", "dSig/dc/dv ", NBv, 0.0 ,vmx2, NBc, -1.0 ,1.0);
  sct_vTcPR_EEX2  = TH2D_UP("sct_vTcPR_EEX2",   "dSig/dc/dv ", NBv, 0.0 ,vmx2, NBc, -1.0 ,1.0);
  //
  sct_vAcPR_Ceex2= TH2D_UP("sct_vAcPR_Ceex2",  "dSig/dc/dv ", NBv, 0.0 ,vmx2, NBc, -1.0 ,1.0);
  sct_vAcPR_Ceex2n= TH2D_UP("sct_vAcPR_Ceex2n","dSig/dc/dv ", NBv, 0.0 ,vmx2, NBc, -1.0 ,1.0);
  //
  sct_vKcPL_Ceex2= TH2D_UP("sct_vKcPL_Ceex2",  "dSig/dc/dv ", NBv, 0.0 ,vmx2, NBc, -1.0 ,1.0);
  sct_vAcPL_Ceex2= TH2D_UP("sct_vAcPL_Ceex2",  "dSig/dc/dv ", NBv, 0.0 ,vmx2, NBc, -1.0 ,1.0);
///////////////// vTcPL series
  sct_vTcPL_Ceex2 = TH2D_UP("sct_vTcPL_Ceex2",  "dSig/dc/dv ", NBv, 0.0 ,vmx2, NBc, -1.0 ,1.0);
  sct_vTcPL_Ceex2n= TH2D_UP("sct_vTcPL_Ceex2n", "dSig/dc/dv ", NBv, 0.0 ,vmx2, NBc, -1.0 ,1.0);
  sct_vTcPL_Ceex1 = TH2D_UP("sct_vTcPL_Ceex1",  "dSig/dc/dv ", NBv, 0.0 ,vmx2, NBc, -1.0 ,1.0);
  sct_vTcPL_Ceex1n= TH2D_UP("sct_vTcPL_Ceex1n", "dSig/dc/dv ", NBv, 0.0 ,vmx2, NBc, -1.0 ,1.0);
  sct_vTcPL_Ceex0 = TH2D_UP("sct_vTcPL_Ceex0",  "dSig/dc/dv ", NBv, 0.0 ,vmx2, NBc, -1.0 ,1.0);
  sct_vTcPL_Ceex0n= TH2D_UP("sct_vTcPL_Ceex0n", "dSig/dc/dv ", NBv, 0.0 ,vmx2, NBc, -1.0 ,1.0);
  /////////////////////////////////////
  // CEEX series
  hst_vT_Ceex1    = TH1D_UP( "hst_vT_Ceex1",      "dSig/dvTrue ", NBv, 0.000 ,vmx2);
  hst_vT_Ceex2    = TH1D_UP( "hst_vT_Ceex2",      "dSig/dvTrue ", NBv, 0.000 ,vmx2);
  hst_vT_Ceex21   = TH1D_UP( "hst_vT_Ceex21",     "dSig/dvTrue ", NBv, 0.000 ,vmx2);
  hst_vT_Ceex1_F  = TH1D_UP( "hst_vT_Ceex1_F",    "dSig/dvTrue ", NBv, 0.000 ,vmx2);
  hst_vT_Ceex2_F  = TH1D_UP( "hst_vT_Ceex2_F",    "dSig/dvTrue ", NBv, 0.000 ,vmx2);
  hst_vT_Ceex21_F = TH1D_UP( "hst_vT_Ceex21_F",   "dSig/dvTrue ", NBv, 0.000 ,vmx2);
  // CEEXn series
  hst_vT_Ceex1n   = TH1D_UP( "hst_vT_Ceex1n",     "dSig/dvTrue ", NBv, 0.000 ,vmx2);
  hst_vT_Ceex2n   = TH1D_UP( "hst_vT_Ceex2n",     "dSig/dvTrue ", NBv, 0.000 ,vmx2);
  hst_vT_Ceex21n  = TH1D_UP( "hst_vT_Ceex21n",    "dSig/dvTrue ", NBv, 0.000 ,vmx2);
  hst_vT_Ceex1n_F = TH1D_UP( "hst_vT_Ceex1n_F",   "dSig/dvTrue ", NBv, 0.000 ,vmx2);
  hst_vT_Ceex2n_F = TH1D_UP( "hst_vT_Ceex2n_F",   "dSig/dvTrue ", NBv, 0.000 ,vmx2);
  hst_vT_Ceex21n_F= TH1D_UP( "hst_vT_Ceex21n_F",  "dSig/dvTrue ", NBv, 0.000 ,vmx2);
  // EEX series
  hst_vT_EEX1     = TH1D_UP( "hst_vT_EEX1",      "dSig/dvTrue ", NBv, 0.000 ,vmx2);
  hst_vT_EEX2     = TH1D_UP( "hst_vT_EEX2",      "dSig/dvTrue ", NBv, 0.000 ,vmx2);
  hst_vT_EEX3     = TH1D_UP( "hst_vT_EEX3",      "dSig/dvTrue ", NBv, 0.000 ,vmx2);
  hst_vT_EEX21    = TH1D_UP( "hst_vT_EEX21",     "dSig/dvTrue ", NBv, 0.000 ,vmx2);
  hst_vT_EEX32    = TH1D_UP( "hst_vT_EEX32",     "dSig/dvTrue ", NBv, 0.000 ,vmx2);
  hst_vT_EEX1_F   = TH1D_UP( "hst_vT_EEX1_F",    "dSig/dvTrue ", NBv, 0.000 ,vmx2);
  hst_vT_EEX2_F   = TH1D_UP( "hst_vT_EEX2_F",    "dSig/dvTrue ", NBv, 0.000 ,vmx2);
  hst_vT_EEX3_F   = TH1D_UP( "hst_vT_EEX3_F",    "dSig/dvTrue ", NBv, 0.000 ,vmx2);
  hst_vT_EEX21_F  = TH1D_UP( "hst_vT_EEX21_F",   "dSig/dvTrue ", NBv, 0.000 ,vmx2);
  hst_vT_EEX32_F  = TH1D_UP( "hst_vT_EEX32_F",   "dSig/dvTrue ", NBv, 0.000 ,vmx2);
////////////////////////////////////
// for special AFB based on <costheta PL>
  hst_vT_Ceex2_xcPL    = TH1D_UP( "hst_vT_Ceex2_xcPL",   "dSig/dvTrue ", NBv, 0.000 ,vmx2);
  hst_vT_Ceex2n_xcPL   = TH1D_UP( "hst_vT_Ceex2n_xcPL",  "dSig/dvTrue ", NBv, 0.000 ,vmx2);   hst_vT_Ceex2_cPLr90  = TH1D_UP( "hst_vT_Ceex2_cPLr90",  "dSig/dvTrue ", NBv, 0.000 ,vmx2);
// for standard AFB from (F-B)/(F+B) |cost(theta|<1.0
  hst_vT_Ceex2_cPL_forw    = TH1D_UP( "hst_vT_Ceex2_cPL_forw",    "dSig/dvTrue ", NBv, 0.000 ,vmx2);
  hst_vT_Ceex2_cPLr90_forw = TH1D_UP( "hst_vT_Ceex2_cPLr90_forw", "dSig/dvTrue ", NBv, 0.000 ,vmx2);
  // older quasi-realistic
  hst_vACeex1   = TH1D_UP("hst_vACeex1",  "dSig/dvTrue ", NBv, 0.000 ,vmx2);
  hst_vACeex2   = TH1D_UP("hst_vACeex2",  "dSig/dvTrue ", NBv, 0.000 ,vmx2);
  hst_vACeex21F = TH1D_UP("hst_vACeex21F","dSig/dvTrue ", NBv, 0.000 ,vmx2);
  hst_vACeex21B = TH1D_UP("hst_vACeex21B","dSig/dvTrue ", NBv, 0.000 ,vmx2);
  hst_vACeex1F  = TH1D_UP("hst_vACeex1F",  "dSig/dvTrue ", NBv, 0.000 ,vmx2);
  hst_vACeex2F  = TH1D_UP("hst_vACeex2F",  "dSig/dvTrue ", NBv, 0.000 ,vmx2);
  hst_vACeex21  = TH1D_UP("hst_vACeex21",  "dSig/dvTrue ", NBv, 0.000 ,vmx2);
  ////////////////////////////////////
  // New section for alpha dependence
  m_vvcut2 = 0.02;   // photon energy cut
  m_DelAlf = 1e-3;   // relative variation of alphaQED
  int NBalf = 8;
  hst_Al20CutA_Ceex2    = TH1D_UP("hst_Al20CutA_Ceex2",   " dsig/dAlf ", NBalf, -m_DelAlf , m_DelAlf );
  hst_Al20CutA_Ceex2F   = TH1D_UP("hst_Al20CutA_Ceex2F",  " dsig/dAlf ", NBalf, -m_DelAlf , m_DelAlf );
  hst_Alf0CutA_Ceex2    = TH1D_UP("hst_Alf0CutA_Ceex2",   " dsig/dAlf ", NBalf, -m_DelAlf , m_DelAlf );
  hst_Alf0CutA_Ceex2F   = TH1D_UP("hst_Alf0CutA_Ceex2F",  " dsig/dAlf ", NBalf, -m_DelAlf , m_DelAlf );
  hst_Alf2CutA_Ceex2    = TH1D_UP("hst_Alf2CutA_Ceex2",   " dsig/dAlf ", NBalf, -m_DelAlf , m_DelAlf );
  hst_Alf2CutA_Ceex2F   = TH1D_UP("hst_Alf2CutA_Ceex2F",  " dsig/dAlf ", NBalf, -m_DelAlf , m_DelAlf );
  //
  hst_Al20CutB_Ceex2    = TH1D_UP("hst_Al20CutB_Ceex2",   " dsig/dAlf ", NBalf, -m_DelAlf , m_DelAlf );
  hst_Al20CutB_Ceex2F   = TH1D_UP("hst_Al20CutB_Ceex2F",  " dsig/dAlf ", NBalf, -m_DelAlf , m_DelAlf );
  hst_Alf0CutB_Ceex2    = TH1D_UP("hst_Alf0CutB_Ceex2",   " dsig/dAlf ", NBalf, -m_DelAlf , m_DelAlf );
  hst_Alf0CutB_Ceex2F   = TH1D_UP("hst_Alf0CutB_Ceex2F",  " dsig/dAlf ", NBalf, -m_DelAlf , m_DelAlf );
  hst_Alf2CutB_Ceex2    = TH1D_UP("hst_Alf2CutB_Ceex2",   " dsig/dAlf ", NBalf, -m_DelAlf , m_DelAlf );
  hst_Alf2CutB_Ceex2F   = TH1D_UP("hst_Alf2CutB_Ceex2F",  " dsig/dAlf ", NBalf, -m_DelAlf , m_DelAlf );
//
  hst_Al20CutA_Ceex2n    = TH1D_UP("hst_Al20CutA_Ceex2n",   " dsig/dAlf ", NBalf, -m_DelAlf , m_DelAlf );
  hst_Al20CutA_Ceex2nF   = TH1D_UP("hst_Al20CutA_Ceex2nF",  " dsig/dAlf ", NBalf, -m_DelAlf , m_DelAlf );
  hst_Alf2CutA_Ceex2n    = TH1D_UP("hst_Alf2CutA_Ceex2n",   " dsig/dAlf ", NBalf, -m_DelAlf , m_DelAlf );
  hst_Alf2CutA_Ceex2nF   = TH1D_UP("hst_Alf2CutA_Ceex2nF",  " dsig/dAlf ", NBalf, -m_DelAlf , m_DelAlf );
//
  hst_Al20CutB_Ceex2n    = TH1D_UP("hst_Al20CutB_Ceex2n",   " dsig/dAlf ", NBalf, -m_DelAlf , m_DelAlf );
  hst_Al20CutB_Ceex2nF   = TH1D_UP("hst_Al20CutB_Ceex2nF",  " dsig/dAlf ", NBalf, -m_DelAlf , m_DelAlf );
  hst_Alf2CutB_Ceex2n    = TH1D_UP("hst_Alf2CutB_Ceex2n",   " dsig/dAlf ", NBalf, -m_DelAlf , m_DelAlf );
  hst_Alf2CutB_Ceex2nF   = TH1D_UP("hst_Alf2CutB_Ceex2nF",  " dsig/dAlf ", NBalf, -m_DelAlf , m_DelAlf );

/*  mooved to TMCgenKKMC
  //  ************* special histo  *************
  HST_KKMC_NORMA = TH1D_UP("HST_KKMC_NORMA","KKMC normalization &xpar",jmax,0.0,10000.0);
  for(int j=1; j<=jmax; j++)
    HST_KKMC_NORMA->SetBinContent(j,m_xpar[j]);    // xpar encoded
*/

  m_YSum= 0.0;
  m_YSum2=0.0;
}

/*
///////////////////////////////////////////////////////////////////////////////
// !!!!!!!!!!!!!! OBSOLETE !!!!!!!!!!!!
void TRobolKKMC::KKMC_NORMA()
{
  // Transfer normalization Record of KKMC to local histogram.
  // For later use in re-normalizing histostograms
  //
  double XsPrim; long NevPrim;
  KKMC_generator->GetPrimaNorma(XsPrim, NevPrim);
  HST_KKMC_NORMA->SetBinContent(0,XsPrim*NevPrim);
  HST_KKMC_NORMA->SetEntries(NevPrim);
  cout<<" TRobolKKMC::KKMC_NORMA: XsPrim, NevPrim ="<< XsPrim <<"  "<<NevPrim << endl;
}
*/

///////////////////////////////////////////////////////////////////////////////
void TRobolKKMC::Production(double &iEvent)
{
/////////////////////////////////////////////////////////////////////////
//
//   GENERATE AND ANALYZE SINGLE EVENT
//
//
/////////////////////////////////////////////////////////////////////////
  // ****************************************************************
  // ************ Generate event and import it here  ****************
  m_NevGen++;
  TMCgenKKMC *KKMC_generator = (TMCgenKKMC*)f_MCgen;
  KKMC_generator->Generate(); // done in user class
  /// *************************************************
  /// KKMC generation in base class
  /// *************************************************
  /*[[[
  PartImport(); // Import Pythia common block content into local LuPart matrix
  */
  // IMPORT KKMC event and weights
  double WtMain,WtCrude;
  KKMC_generator->GetWt(WtMain,WtCrude);
  KKMC_generator->GetBeams(m_pbea1,m_pbea2);
  KKMC_generator->GetFermions(m_pfer1,m_pfer2);
  KKMC_generator->GetNphot(m_Nphot);                  // photon multiplicity
  TLorentzVector VSumPhot;    // By default all components are initialized by zero.
  long iphot,iphot1;
  for(iphot=0;iphot<m_Nphot;iphot++){
    KKMC_generator->GetPhoton1(iphot+1,m_phot[iphot]);  // photon 4-momenta
    VSumPhot+= m_phot[iphot];
  }
  if(iEvent<10){
    cout<<"-----------------------------------------------------------  "<<iEvent;
    cout<<"  -----------------------------------------------------------"<<endl;
    cout<<"  WtMain= "<< WtMain  << "  WtCrude= "<< WtCrude<<endl;
    cout<<" m_Nphot= "<< m_Nphot<<endl;
    cout<<"VSumPhot= "; MomPrint( VSumPhot );
    KKMC_generator->Print1();
    //KKMC_generator->PyList(2);
    //PyPrint(1);
  }
  // ****************************************************************
  double s  =(m_pbea1+m_pbea2)*(m_pbea1+m_pbea2);
  double s1 =(m_pfer1+m_pfer2)*(m_pfer1+m_pfer2);
  double CMSene = sqrt(s);
  double Mff    = sqrt(s1);
  double vv     = 1-s1/s;
  // ********************************************************************
  // ***   Photon trigger TrigPho is for everybory, all pions, muons etc
  double Pi=4*atan(1.0);
  double phEne,phTheta,phCosth;
  double XEnePho  = 0.010;              // Emin for visible photom
  //**************************************************
  // Loop over photons, just in case
  //**************************************************
  int nph_ene=0;
  for(iphot=0;iphot<m_Nphot;iphot++){
    phEne   = m_phot[iphot].Energy();
    phCosth = m_phot[iphot].CosTheta();
    phTheta = m_phot[iphot].Theta()*180/Pi;
    if(phEne>XEnePho){
      nph_ene++;
    }
  }

  /*[[[[
  //********************************************************************
  // Muon trigger, it is not realy necessary if MC ir run for mu only
  //********************************************************************
  int TrigMu  = 0;
  // find muons, excluding muons from phi decays!!!
  long jMu1 =PartFindStable( 13);    // fortran numbering!!!
  long jMu2 =PartFindStable(-13);    // fortran numbering!!!
  m_pMu1  = m_Event[jMu1-1].fMom;    // fortran numbering!!!
  m_pMu2  = m_Event[jMu2-1].fMom;    // fortran numbering!!!
  long par1=m_Event[jMu1-1].fParent; // fortran numbering!!!
  long par2=m_Event[jMu2-1].fParent; // fortran numbering!!!
  if( (jMu1*jMu1)  && (par1 == par2) && (par1 == 3) ) TrigMu  = 1; // exclude backgr.
  //**************************************************************
  if( TrigMu && (m_count1<17) ){
    m_count1++;
    cout<<"**************************>>> two muons <<<****************************"<<endl;
    KKMC_generator->PyList(2);
  }
  // muons,  vv, Q^2 costheta, etc
  double CosThe1 = m_pMu1.CosTheta();
  double Theta1  = m_pMu1.Theta();
  double E1      = m_pMu1.Energy();
  double CosThe2 = m_pMu2.CosTheta();
  double Theta2  = m_pMu2.Theta();
  double E2      = m_pMu2.Energy();
  */

  double CosThe1 = m_pfer1.CosTheta();
  double Theta1  = m_pfer1.Theta();
  double Phi1    = m_pfer1.Phi();
  double E1      = m_pfer1.Energy();

  double CosThe2 = m_pfer2.CosTheta();
  double Phi2    = m_pfer2.Phi();
  double Theta2  = m_pfer2.Theta();
  double E2      = m_pfer2.Energy();

  double Acol    = fabs(Phi1-(Phi2+3.141594));

  double SinThe1,SinThe2,yy1,yy2,CosThePL,CosPRD,zAleph,s1Aleph,xe1,xe2;
//--------------------------------------------------------------------
//* Various definitions of Theta and s-propagator
//*--------------------------------------------------------------------
//** definition of P.L. B219, 103 (1989)
  CosThePL = ( E1*CosThe1 -E2*CosThe2)/(E1+E2);
//* definition of P.R. D41, 1425 (1990)
  SinThe1 = sqrt(fabs((1-CosThe1)*(1+CosThe1)));
  SinThe2 = sqrt(fabs((1-CosThe2)*(1+CosThe2)));
  yy1 = SinThe2/(SinThe1+SinThe2);
  yy2 = SinThe1/(SinThe1+SinThe2);
  CosPRD = yy1*CosThe1 - yy2*CosThe2;
  xe1 = 2*E1/CMSene;
  xe2 = 2*E2/CMSene;
//*-------------------------------
//* LL formula for s'/s from angles according to ALEPH note 1996
  zAleph =  (sin(Theta1)+sin(Theta2) -fabs(sin(Theta1+Theta2)))
           /(sin(Theta1)+sin(Theta2) +fabs(sin(Theta1+Theta2)));
  s1Aleph  = s*zAleph;
  //
  double vvK,x1,x2;
  karlud_getvvxx_(vvK,x1,x2);    // vv of ISR from MC (illegal in case IFI on)
  double vvF  = 1-(1-vv)/(1-vvK);  // vv of FSR
  // *********************************************************************
  //          Most of histogramming starts here
  // *********************************************************************
  double WtEEX1  = KKMC_generator->GetWtAlter( 72);    //  Second ord. EEX2 O(alf2)
  double WtEEX2  = KKMC_generator->GetWtAlter( 73);    //  Second ord. EEX2 O(alf2)
  double WtEEX3  = KKMC_generator->GetWtAlter( 74);    //  Third order EEX3 O(alf3)
  //cout<< "&&&&&& WtEEX2,3= "<<WtEEX2<<"  "<<WtEEX3<<endl;
  double WtCEEX1 = KKMC_generator->GetWtAlter(202);    //  CEEX Weight O(alf1)
  double WtCEEX1n= KKMC_generator->GetWtAlter(252);    //  CEEX Weight O(alf1)
  double WtCEEX2 = KKMC_generator->GetWtAlter(203);    //  CEEX Weight O(alf2)
  double WtCEEX2n= KKMC_generator->GetWtAlter(253);    //  CEEX Weight O(alf2) IFI off
  double WtCEEX0 = KKMC_generator->GetWtAlter(201);    //  CEEX Weight O(alf0)
  double WtCEEX0n= KKMC_generator->GetWtAlter(251);    //  CEEX Weight O(alf0) IFI off
  //
  double vvA = 1-zAleph;
  hst_nPhAll->Fill(  m_Nphot,WtMain);
  hst_nPhVis->Fill(  nph_ene,WtMain);
  hst_vTrueMain->Fill(   vv, WtMain);
  hst_vTrueCeex2->Fill(  vv, WtCEEX2);          // M(2f) of mun pair
  hst_vXGenCeex2->Fill(  vvK, WtCEEX2);         // M^star from MC (illegal)
  // semirealistic
  hst_vAlepCeex2->Fill(  vvA, WtCEEX2);         // M^star guessed
  /////////////////////////
  // older histos for quasi-experimental variables
  hst_vACeex1->Fill(     vvA, WtCEEX1);         // M^star guessed
  hst_vACeex2->Fill(     vvA, WtCEEX2);         // M^star guessed
  hst_vACeex21->Fill(    vvA, WtCEEX2-WtCEEX1); // M^star guessed
  if( CosPRD > 0) hst_vACeex1F->Fill(  vvA, WtCEEX1);
  if( CosPRD > 0) hst_vACeex2F->Fill(  vvA, WtCEEX2);
  if( CosPRD > 0) hst_vACeex21F->Fill( vvA, WtCEEX2-WtCEEX1); // M^star guessed
  if( CosPRD < 0) hst_vACeex21B->Fill( vvA, WtCEEX2-WtCEEX1); // spurious
  //[[[[[[[[[[[[[[[[[[[[
  // New histos for academic variables
  // CEEX series
  hst_vT_Ceex1->Fill(   vv, WtCEEX1);
  hst_vT_Ceex2->Fill(   vv, WtCEEX2);
  hst_vT_Ceex21->Fill(  vv, WtCEEX2-WtCEEX1);
  if( CosThePL > 0.0) hst_vT_Ceex1_F->Fill(  vv, WtCEEX1);
  if( CosThePL > 0.0) hst_vT_Ceex2_F->Fill(  vv, WtCEEX2);
  if( CosThePL > 0.0) hst_vT_Ceex21_F->Fill( vv, WtCEEX2-WtCEEX1);
  // CEEXn series
  hst_vT_Ceex1n->Fill(   vv, WtCEEX1n);
  hst_vT_Ceex2n->Fill(   vv, WtCEEX2n);
  hst_vT_Ceex21n->Fill(  vv, WtCEEX2n-WtCEEX1n);
  if( CosThePL > 0.0) hst_vT_Ceex1n_F->Fill(  vv, WtCEEX1n);
  if( CosThePL > 0.0) hst_vT_Ceex2n_F->Fill(  vv, WtCEEX2n);
  if( CosThePL > 0.0) hst_vT_Ceex21n_F->Fill( vv, WtCEEX2n-WtCEEX1n);
  // EEX series
  hst_vT_EEX1->Fill(   vv, WtEEX1);
  hst_vT_EEX2->Fill(   vv, WtEEX2);
  hst_vT_EEX3->Fill(   vv, WtEEX3);
  hst_vT_EEX21->Fill(  vv, WtEEX2-WtEEX1);
  hst_vT_EEX32->Fill(  vv, WtEEX3-WtEEX2);
  if( CosThePL > 0.0) hst_vT_EEX1_F->Fill(   vv, WtEEX1);
  if( CosThePL > 0.0) hst_vT_EEX2_F->Fill(   vv, WtEEX2);
  if( CosThePL > 0.0) hst_vT_EEX3_F->Fill(   vv, WtEEX3);
  if( CosThePL > 0.0) hst_vT_EEX21_F->Fill(  vv, WtEEX2-WtEEX1);
  if( CosThePL > 0.0) hst_vT_EEX32_F->Fill(  vv, WtEEX3-WtEEX2);
  //]]]]]]]]]]]]]]]]]]]]
  if(vv<0.9){
    hst_Cost1Ceex2->Fill( CosThe1, WtCEEX2);
    hst_CosPLCeex2->Fill(CosThePL, WtCEEX2);
    hst_CosPRCeex2->Fill(  CosPRD, WtCEEX2);
  }
  hst_CosPREex2->Fill(  CosPRD, WtEEX2);
  // big scatergrams, range vv< 1.0
  sca_vTcPR_Ceex2->Fill(   vv, CosPRD, WtCEEX2);
  sca_vTcPR_Ceex2n->Fill(  vv, CosPRD, WtCEEX2n); // true v, IFI off
  sca_vTcPR_Eex2->Fill(    vv, CosPRD, WtEEX2);

  sca_vTcPL_Ceex2->Fill(   vv, CosThePL, WtCEEX2);
  sca_vTcPL_Ceex2n->Fill(  vv, CosThePL, WtCEEX2n); // true v, IFI off
  sca_vTcPL_Eex2->Fill(    vv, CosThePL,  WtEEX2);   // true v, IFI off

  sca_vTcPL_Ceex0->Fill(   vv, CosThePL, WtCEEX0);
  sca_vTcPL_Ceex0n->Fill(  vv, CosThePL, WtCEEX0n); // true v, IFI off

  sca_vXcPR_Ceex2->Fill(   vvK, CosPRD, WtCEEX2);
  sca_vXcPR_Eex2->Fill(    vvK, CosPRD, WtEEX2);
  //-------------------------------------------
  // New very BIG scaterplots, restricted range vv<0.20
  sct_vTcPR_Ceex2->Fill(   vv, CosPRD, WtCEEX2);  // true v, IFI on
  sct_vTcPR_Ceex2n->Fill(  vv, CosPRD, WtCEEX2n); // true v, IFI off
  sct_vTcPR_EEX2->Fill(    vv, CosPRD, WtEEX2);   // true v, IFI off

  sct_vAcPR_Ceex2->Fill(   vvA, CosPRD,   WtCEEX2);  // Main CEEX2 KKMC , ISR+FSR
  sct_vAcPR_Ceex2n->Fill(  vvA, CosPRD,   WtCEEX2n); // IFI  off
  sct_vKcPL_Ceex2->Fill(   vvK, CosThePL, WtCEEX2);  // vv of Karlud (unphysical, pure ISR) thetaPL
////////////////////
  sct_vTcPL_Ceex2->Fill(    vv, CosThePL, WtCEEX2);  // vv bare muons
  sct_vTcPL_Ceex2n->Fill(   vv, CosThePL, WtCEEX2n); // vv bare muons
  sct_vTcPL_Ceex1->Fill(    vv, CosThePL, WtCEEX1);  // vv bare muons, IFIon
  sct_vTcPL_Ceex1n->Fill(   vv, CosThePL, WtCEEX1n); // vv bare muons, IFIoff
  sct_vTcPL_Ceex0->Fill(    vv, CosThePL, WtCEEX0);  // vv bare muons, IFIon
  sct_vTcPL_Ceex0n->Fill(   vv, CosThePL, WtCEEX0n); // vv bare muons, IFIoff
////////////////////
  sct_vAcPL_Ceex2->Fill(   vvA, CosThePL, WtCEEX2);  // Main CEEX2 KKMC , ISR+FSR
//-------------------------------
// dsigma/dv, any theta, no cut
// specials for AFB from <costheta_PL>
  hst_vT_Ceex2_xcPL->Fill( vv, WtCEEX2*CosThePL);
  hst_vT_Ceex2n_xcPL->Fill(vv, WtCEEX2n*CosThePL);
// AFB from (F-B)/(F+B), cos(theta)>0, cmax=1
  if( CosThePL > 0.0)  hst_vT_Ceex2_cPL_forw->Fill( vv, WtCEEX2);
// AFB from (F-B)/(F+B) with |cos(theta)| < 0.9 cut
  if( fabs(CosThePL) < 0.9 )          hst_vT_Ceex2_cPLr90->Fill(      vv, WtCEEX2);
  if( CosThePL>0.0 && CosThePL < 0.9) hst_vT_Ceex2_cPLr90_forw->Fill( vv, WtCEEX2);
///////////////////////////////////////////////////////////////////////
// Introductory exercises on experimental cut-offs
///////////////////////////////////////////////////////////////////////
// Acoplanarity < 0.35 mrad (i.e., ~3 sigma of the angular resolution).
// Of course, this cut can be relaxed a little bit, but not too much.
//  if( m_Nphot == 1){
  if( Acol < 0.00035){
//  if( Acol < 0.001){
    sca_vTvA_Eex2->Fill(   vv, vvA, WtEEX2); // vv bare vs vALEPH
//  sca_vKvA_Eex2->Fill(  vvK, vvA, WtEEX2); // vv bare vs vALEPH
    sca_vKvA_Eex2->Fill(  vvF, vvA, WtEEX2); // vv bare vs vALEPH
  }// Cuts
///////////////////////////////////////////////////////////////////////
// Miscelaneous
  m_YSum  += WtMain;
  m_YSum2 += WtMain*WtMain;
  hst_weight->Fill(WtMain);              // histogramming
//hst_Mff->Fill(Mff,WtMain);             // histogramming
// debug debug debug debug debug debug debug
  if(iEvent<15){
    cout<<"-----------------------------------------------------------  "<<iEvent;
    cout<<"  -----------------------------------------------------------"<<endl;
    cout<< "vv, 1-zAleph    = "<<     vv<<"  "<<1-zAleph<<endl;
    cout<< "CosThe1,CosThe2 = "<<CosThe1<<"  "<<CosThe2<<endl;
    cout<< "CosThePL,CosPRD = "<<CosThePL<<"  "<<CosPRD<<endl;
  }
///////////////////////////////////////////////////////////////////////
// !!!!!!!!!NEW!!!  NEW!!!  NEW!!!  NEW!!!
///////////////////////////////////////////////////////////////////////

//  double Xmax = hst_Alf0CutA_Ceex2->GetXaxis()->GetXmax();
//  double Xmin = hst_Alf0CutA_Ceex2->GetXaxis()->GetXmin();
//  dalf = Xmin + (Xmax-Xmin)*RN;

  double RN = KKMC_generator->f_RNgen->Rndm();
  double dalf = -m_DelAlf + 2*m_DelAlf*RN;
//
// Modifying alpha_QED
  double xmodif = 1 + dalf;
  bornv_setqedmodif_(xmodif);
// recalculating WC weight
  kk2f_make_wt_();
// Back to normal alpha !!!
  xmodif = 1.0e0;  // Back to normal !!!
  bornv_setqedmodif_( xmodif);
 /////////////////////////////
  double Wt2CEEX1 = KKMC_generator->GetWtAlter(202);    //  CEEX Weight O(alf1)
  double Wt2CEEX1n= KKMC_generator->GetWtAlter(252);    //  CEEX Weight O(alf1)
  double Wt2CEEX2 = KKMC_generator->GetWtAlter(203);    //  CEEX Weight O(alf2)
  double Wt2CEEX2n= KKMC_generator->GetWtAlter(253);    //  CEEX Weight O(alf2) IFI off
  double Wt2EEX1  = KKMC_generator->GetWtAlter( 72);    //  Second ord. EEX2 O(alf2)
  double Wt2EEX2  = KKMC_generator->GetWtAlter( 73);    //  Second ord. EEX2 O(alf2)
  double Wt2EEX3  = KKMC_generator->GetWtAlter( 74);    //  Third order EEX3 O(alf3)
  // (B) Classic cut on "bare" vv from muon pair mass
  if( vv< m_vvcut2 && fabs(CosThePL)<0.9 ){
	                    hst_Al20CutB_Ceex2->Fill( dalf,Wt2CEEX2-WtCEEX2);
	  if( CosThePL>0.0) hst_Al20CutB_Ceex2F->Fill(dalf,Wt2CEEX2-WtCEEX2);
	                    hst_Alf0CutB_Ceex2->Fill( dalf,WtCEEX2);
      if( CosThePL>0.0) hst_Alf0CutB_Ceex2F->Fill(dalf,WtCEEX2);
                        hst_Alf2CutB_Ceex2->Fill( dalf,Wt2CEEX2);
      if( CosThePL>0.0) hst_Alf2CutB_Ceex2F->Fill(dalf,Wt2CEEX2);
  ///////////////////////
                        hst_Al20CutB_Ceex2n->Fill( dalf,Wt2CEEX2n-WtCEEX2n);
      if( CosThePL>0.0) hst_Al20CutB_Ceex2nF->Fill(dalf,Wt2CEEX2n-WtCEEX2n);
                        hst_Alf2CutB_Ceex2n->Fill( dalf,Wt2CEEX2n);
      if( CosThePL>0.0) hst_Alf2CutB_Ceex2nF->Fill(dalf,Wt2CEEX2n);

  }// Cut (B)
  //
  // (A) ALEPHI-like cut on vv_Aleph, acoplanarity and muon energies
  if( vvA< m_vvcut2 && fabs(CosPRD)<0.9  && Acol<0.001 && xe1>0.90 && xe2>0.90){
	                  hst_Al20CutA_Ceex2->Fill( dalf,Wt2CEEX2-WtCEEX2);
	  if( CosPRD>0.0) hst_Al20CutA_Ceex2F->Fill(dalf,Wt2CEEX2-WtCEEX2);
                      hst_Alf0CutA_Ceex2->Fill( dalf,WtCEEX2);
      if( CosPRD>0.0) hst_Alf0CutA_Ceex2F->Fill(dalf,WtCEEX2);
                      hst_Alf2CutA_Ceex2->Fill( dalf,Wt2CEEX2);
      if( CosPRD>0.0) hst_Alf2CutA_Ceex2F->Fill(dalf,Wt2CEEX2);
  ///////////////////////
                      hst_Al20CutA_Ceex2n->Fill( dalf,Wt2CEEX2n-WtCEEX2n);
      if( CosPRD>0.0) hst_Al20CutA_Ceex2nF->Fill(dalf,Wt2CEEX2n-WtCEEX2n);
                      hst_Alf2CutA_Ceex2n->Fill( dalf,Wt2CEEX2n);
      if( CosPRD>0.0) hst_Alf2CutA_Ceex2nF->Fill(dalf,Wt2CEEX2n);
  }// Cut (A)

  if(iEvent<50){
    cout<<"============================================================="<<iEvent;
    cout<<"============================================================="<<endl;
//    cout<< "Wt2CEEX2 = "<<  Wt2CEEX2<<" WtCEEX2 "<<  WtCEEX2<<endl;
    cout<< " dalf= "<<  dalf << "  xmodif= "<<  1+dalf ;
    cout<< " CEEX2 ratio = "<<  Wt2CEEX2/WtCEEX2 ;
    cout<< "  EEX2 ratio = "<<  Wt2EEX2/WtEEX2 <<endl;
  }


} //


///////////////////////////////////////////////////////////////////////////////
void TRobolKKMC::Finalize()
{
//   Finalize MC  run, final printouts, cleaning etc., xcheck of normalization
//   Plotting histograms is done independently using root file
  TMCgenKKMC *KKMC_generator = (TMCgenKKMC*)f_MCgen;
//
  double XsNormPb, XsErroPb;
  KKMC_generator->Finalize();
  XsNormPb =KKMC_generator->m_XsNormPb;
  XsErroPb =KKMC_generator->m_XsErroPb;
  cout << " KKMC: XsNormPb [pb] = "<<  XsNormPb << "  +-  "<< XsErroPb <<endl;
  double xSecPb,xErrPb,xSecNb;
  KKMC_generator->GetXsecMC( xSecPb, xErrPb);
  xSecNb=xSecPb/1000;
  cout << " KKMC: xSecPb   [pb] = "<<  xSecPb << "  +-  "<< xErrPb <<endl;
  cout << " KKMC: xSecNb   [nb] = "<<  xSecNb << "  +-  "<< xErrPb/1000 <<endl;
  // *********************************************************************
  // **** examples of normalizing histogram (will not work on farm)  *****
  /*
  int      nbt  = hst_Mff->GetNbinsX();
  Double_t tmax = hst_Mff->GetXaxis()->GetXmax();
  Double_t tmin = hst_Mff->GetXaxis()->GetXmin();
  Double_t Fact = nbt*XsNormPb/1000/(tmax-tmin)/m_NevGen; // now [nb]
  // **** re-normalized histogram as a clone
  TH1D *hstC_Mff =(TH1D*)hst_Mff->Clone();
  hstC_Mff->SetName("hstC_Mff"); // otherwise you get 2 histograms with same name
  hstC_Mff->SetTitle("dsigma/dQ2 [nb/GeV^2] ");
  hstC_Mff->Scale(Fact);
  // **** alternatively, re-normalized histo can be defined as a new one ****
  TH1D *hstN_Mff =TH1D_UP("hstN_Mff","dsigma/dQ2 [nb/GeV^2]",nbt,tmin,tmax);
  hstN_Mff->Add(hstN_Mff, Fact);
  */
  // *********************************************************************
  // **** alternatively HST_KKMC_NORMA is used at later stage (plotting)
  /*
  long   NevPrim = HST_KKMC_NORMA->GetEntries();
  double XsPrima = HST_KKMC_NORMA->GetBinContent(0)/NevPrim;
  cout << "HST_KKMC_NORMA: XsPrima [nb] = "<< XsPrima << " NevPrim= "<< NevPrim <<endl;
  */
  // *********************************************************************
  // Integrated Xsection canculated on-line
  cout << "///////////////////////////////////////////////////////////////"<<endl;
  double  XsecY = m_YSum/m_NevGen;                                // average weight
  double dXsecY = sqrt((m_YSum2/m_NevGen-XsecY*XsecY)/m_NevGen) ; // dispersion of wt.
  XsecY  *= XsNormPb/1000.0;  //  nanob.
  dXsecY *= XsNormPb/1000.0;  //  nanob.
  cout << " XsecY [nb] = "<< XsecY <<" +- "<< dXsecY <<endl;
  cout << "///////////////////////////////////////////////////////////////"<<endl;
}
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                             UTILITIES                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

/*[[[
///////////////////////////////////////////////////////////////////////////////
void TRobolKKMC::PartImport(){
/////////////////////////////////////////////////////////////////////////
// Import Pythia common block content into local LuPart matrix
/////////////////////////////////////////////////////////////////////////
  m_Npart= ((TMCgenKKMC*)f_MCgen)->GetPyNpart();
  if( (m_Npart<0) || (m_Npart>=4000) ){
    cout<<"++++ TRobolKKMC::Production: STOP m_Npart= "<<m_Npart<<endl;
    exit(5);
  }
  for(long j=0; j<m_Npart;j++){
    ((TMCgenKKMC*)f_MCgen)->GetPyParticle( j, m_Event[j]);  // import one particle
    //m_Event[j].Print(1);
  }
}



///////////////////////////////////////////////////////////////////////////////
long TRobolKKMC::PartCount(const long flav){
  long jCount=0;
  for(long j=0; j<m_Npart;j++)
    if(m_Event[j].fFlafor == flav){
      jCount++;
    }
  return jCount;
}


///////////////////////////////////////////////////////////////////////////////
long TRobolKKMC::PartFindAny(const long flav){
// fortran numbering!!!
  long jPosition=0;
  for(long j=0; j<m_Npart;j++)
    if(m_Event[j].fFlafor == flav){
      jPosition=j+1; break;
    }
  return jPosition;
}
///////////////////////////////////////////////////////////////////////////////
long TRobolKKMC::PartFindStable(const long flav){
// fortran numbering!!!
  long jPosition=0;
  for(long j=0; j<m_Npart;j++)
    if((m_Event[j].fStatus  == 1)&&(m_Event[j].fFlafor == flav)){
      jPosition=j+1; break;
    }
  return jPosition;
}
////////////////////////////////////////////////////////////////////////////////
void TRobolKKMC::PyPrint(const int mode){
//
// PRINT entire Pythia EVENT
//
  TLorentzVector Sum;
  cout<<"  "<<endl;
  cout<<"lser status flavor parent child1 child2";
  cout<<"                Px                Py                Pz               Ene";
  cout<<"              Mass"<<endl;
  for(long j=0; j<m_Npart;j++){
    m_Event[j].Print(mode);
    if(m_Event[j].fStatus == 1) Sum += m_Event[j].fMom;
  }
  cout<<"                  Total 4-momentum --> ";
  MomPrint(Sum);
}
*/

void TRobolKKMC::MomPrint( TLorentzVector &Vect){
//////////////////////////////////////////////////////////////
// printing entire four-vector in one line (with endline)
  for ( int k=0; k < 4 ; k++ )   cout << sw2 << Vect[k];
  cout<<endl;
}
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//           End of Class ROBOL                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
