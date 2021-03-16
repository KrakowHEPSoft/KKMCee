//                CLASS ROBOL                                                //

#include "TRobolFoam.h"

# define sw2 setprecision(10) << setw(18)

ClassImp(TRobolFoam);


TRobolFoam::TRobolFoam():
  TRobol()
{
  /// This constructor is for ROOT streamers ONLY
  cout<< "@@@@> TRobolFoam DEFAULT Constructor (for ROOT only) "<<endl;
}


///_____________________________________________________________
TRobolFoam::TRobolFoam(const char* Name):
  TRobol(Name)
{
//! Constructor to be used by the user!!!
//! Its important role is to define ALL DEFAULTS.
//! to changed by the user before calling TMCgen::Initialize
  cout<< "@@@@> TRobolFoam::TRobolFoam USER Constructor "<<endl;
  m_NevGen=0;
  m_count1=0;
////////////////////////////////////////////////////////////////////////////
// DEFAULT Analysis parameters, can be reset in Start.C !!!!
// Lepton Cuts
  m_Mllmin  = 60.0;
  m_Mllmax  = 150.0;
///////////////////////////////////////////////////////////////////////////////
}///TRobolFOAM


///______________________________________________________________________________________
TRobolFoam::~TRobolFoam()
{
  //!Explicit destructor
  cout<< "@@@@> TRobolFoam::TRobolFoam !!!! DESTRUCTOR !!!! "<<endl;
}///destructor


//______________________________________________________________________________
//////////////////////////////////////////////////////////////
//   Initialize MC generator and analysis programs          //
//////////////////////////////////////////////////////////////
void TRobolFoam::Initialize(
        ofstream *OutFile, /// Central log-file for messages
        TFile *GenFile,    /// ROOT disk file for CRNG and MC gen.
        TFile *HstFile)    /// ROOT disk file for histograms
{
  cout<< "****> TRobolFoam::Initialize starts"<<endl;
////////////////////////////////////////////////////////////////
// Initializing MC generator is done in base class
// f_MCgen and f_RNgen are read from GenFile
// h_TMCgen_NORMA is stored HstFile
// f_MCgen is initialised
  TRobol::Initialize(OutFile,GenFile,HstFile);
//////////////////////////////////////////////
// Book histograms or read them from the disk
  Hbooker();
//////////////////////////////////////////////
  KKeeFoam *MCgen = (KKeeFoam*)f_MCgen;
//
  int jmax = MCgen->maxPar;
// special histogram name "HST_FOAM_NORMA5" set in Start.C ???!!!
  TH1D *h_TMCgen_NORMA = (TH1D*)HstFile->Get("h_TMCgen_NORMA");
  h_TMCgen_NORMA->SetEntries(0);
// storing input data for later use in analysis programs
  for(int j=1; j<=jmax; j++) h_TMCgen_NORMA->SetBinContent(j,  MCgen->m_xpar[j]  );

  f_HstFile->ls();
  f_HstFile->GetListOfKeys()->Print();

  cout<< "****> TRobolFoam::Initialize: finished"<<endl;
}///Initialize



///////////////////////////////////////////////////////////////////////////////
void TRobolFoam::Hbooker()
{
  ///
  cout<< "****> TRobolFoam::Hbooker: histogram booking STARTS"<<endl;
  BXOPE(*f_Out);
  BXTXT(*f_Out,"########################################");
  BXTXT(*f_Out,"###### TRobolFoam::Hbooker  ###########");
  BXTXT(*f_Out,"########################################");
  BXCLO(*f_Out);
  cout<<       "########################################"<<endl;
  cout<<       "###### TRobolFoam::Hbooker  ###########"<<endl;
  cout<<       "########################################"<<endl;
  f_HstFile->cd();
  //
  // *************************************
  // Weight monitoring
  double WtMax = 3.5;
  int NBwt =200;
 // mon_WtFoam = new THwtMon("mon_WtFoam", "WtFoam",  NBwt, WtMax);
//  mon_WtMain = new THwtMon("mon_WtMain", "WtMain",  NBwt, WtMax);
//
  hst_weight   = TH1D_UP("hst_weight" , "MC weight",     NBwt,   0.0 , WtMax);
  hst6_weight  = TH1D_UP("hst6_weight" , "MC weight",    NBwt,   0.0 , WtMax);

  int nBins = 120;
  Hst_Mll   = TH1D_UP("Hst_Mll",   "Dilepton Invariant Mass with IFI;M_{ll};d#sigma/dM_{Z} (nb/GeV)", nBins, m_Mllmin, m_Mllmax);

  Hst_Mll_eex0  = TH1D_UP("Hst_Mll_eex0",   "M_{ll};d#sigma/dM_{Z} (nb/GeV)", nBins, m_Mllmin, m_Mllmax);
  Hst_MllF_eex0 = TH1D_UP("Hst_MllF_eex0",  "AFB(M_{ll})",                    nBins, m_Mllmin, m_Mllmax);

  Hst_Mll_eex2  = TH1D_UP("Hst_Mll_eex2",   "M_{ll};d#sigma/dM_{Z} (nb/GeV)", nBins, m_Mllmin, m_Mllmax);
  Hst_MllF_eex2 = TH1D_UP("Hst_MllF_eex2",  "AFB(M_{ll})",                    nBins, m_Mllmin, m_Mllmax);
////////////////////////////////
  Hst_Mll_ceex2  = TH1D_UP("Hst_Mll_ceex2",   "M_{ll};d#sigma/dM_{Z} (nb/GeV)", nBins, m_Mllmin, m_Mllmax);
  Hst_MllF_ceex2 = TH1D_UP("Hst_MllF_ceex2",  "AFB(M_{ll})",                    nBins, m_Mllmin, m_Mllmax);
//
  hst6_Mll_ceex2  = TH1D_UP("hst6_Mll_ceex2",   "M_{ll};d#sigma/dM_{Z} (nb/GeV)", nBins, m_Mllmin, m_Mllmax);
  hst6_MllF_ceex2 = TH1D_UP("hst6_MllF_ceex2",  "AFB(M_{ll})",                    nBins, m_Mllmin, m_Mllmax);
//////
  Hst_Mll_ceex0  = TH1D_UP("Hst_Mll_ceex0",   "M_{ll};d#sigma/dM_{Z} (nb/GeV)", nBins, m_Mllmin, m_Mllmax);
  Hst_MllF_ceex0 = TH1D_UP("Hst_MllF_ceex0",  "AFB(M_{ll})",                    nBins, m_Mllmin, m_Mllmax);
//
  hst6_Mll_ceex0  = TH1D_UP("hst6_Mll_ceex0",   "M_{ll};d#sigma/dM_{Z} (nb/GeV)", nBins, m_Mllmin, m_Mllmax);
  hst6_MllF_ceex0 = TH1D_UP("hst6_MllF_ceex0",  "AFB(M_{ll})",                    nBins, m_Mllmin, m_Mllmax);

  cout<<       "################Hbooker end##############"<<endl;

}//Hbooker


///////////////////////////////////////////////////////////////////////////////
void TRobolFoam::Production(double &iEvent)
{
/////////////////////////////////////////////////////////////////////////
//
//   GENERATE AND ANALYZE SINGLE EVENT
//
/////////////////////////////////////////////////////////////////////////
  m_NevGen++;
  double xx, vv, CosTheta;
  KKeeFoam *MCgen = (KKeeFoam*)f_MCgen;
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// MC generation using FOAM of the base class
//////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////
  MCgen->m_FoamMode = -7;
  TRobol::Production(iEvent);  // It invokes MCgen->Generate() !!!
/// filling in histos
  double WTfoam;
  MCgen->f_FoamI->GetMCwt(WTfoam);
  xx  = MCgen->m_xx;
  CosTheta = MCgen->m_CosTheta;
  double x1= 1.0 -MCgen->m_y1;
  double x2= 1.0 -MCgen->m_y2;
  // CosTheta knows boost direction of Z
  if(x1<x2) CosTheta *= -1.0; // exploit q/qbar differences in PDF
  double Mll = MCgen->m_CMSene * sqrt( MCgen->m_xx);
  //Mll = MCgen->m_XXXene; // this is for KeyISR=0;

  // monitoring main weight, all events
  //mon_WtFoam->Fill(WTfoam);
  hst_weight->Fill(WTfoam,1.0);

  //////////////////////////////////////////////////
  Hst_Mll->Fill(  Mll, WTfoam);

  double WTborn= WTfoam *MCgen->m_WTset[71]; // EEX0
  Hst_Mll_eex0->Fill(  Mll, WTborn);
  if( CosTheta> 0.0) Hst_MllF_eex0->Fill(  Mll, WTborn);

  double WT73= WTfoam *MCgen->m_WTset[73];  // EEX2
  Hst_Mll_eex2->Fill(  Mll, WT73);
  if( CosTheta> 0.0) Hst_MllF_eex2->Fill(  Mll, WT73);

  double WTceex2= WTfoam *MCgen->m_WTset[253];  // CEEX2 no IFI
  Hst_Mll_ceex2->Fill(  Mll, WTceex2);
  if( CosTheta> 0.0) Hst_MllF_ceex2->Fill(  Mll, WTceex2);

  double WTceex0= WTfoam *MCgen->m_WTset[251];  // ceex0 noIFI
  Hst_Mll_ceex0->Fill(  Mll, WTceex0);
  if( CosTheta> 0.0) Hst_MllF_ceex0->Fill(  Mll, WTceex0);

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
// MC generation using additional FOAM object
/////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
  MCgen->m_FoamMode = -9;
  MCgen->Generate();               //!!!!!!
  MCgen->m_Foam6->GetMCwt(WTfoam);
// event has changed!!!
  hst6_weight->Fill(WTfoam,1.0);
  xx  = MCgen->m_xx;               // event has changed!!!
  CosTheta = MCgen->m_CosTheta;    // event has changed!!!
  Mll = MCgen->m_CMSene * sqrt( MCgen->m_xx);
//
/// filling in histos
  double WTceex2i= WTfoam *MCgen->m_WTset[203];  // CEEX2
//  cout<< "  Production: WTfoam="<<WTfoam<<"  m_WTset[203]="<< MCgen->m_WTset[203]<<" WTceex2i="<<WTceex2i <<endl;
  hst6_Mll_ceex2->Fill(  Mll, WTceex2i);
  if( CosTheta> 0.0) hst6_MllF_ceex2->Fill(  Mll, WTceex2i);

  double WTceex0i= WTfoam *MCgen->m_WTset[201];  // ceex0+IFI
  hst6_Mll_ceex0->Fill(  Mll, WTceex0i);
  if( CosTheta> 0.0) hst6_MllF_ceex0->Fill(  Mll, WTceex0i);

}// Production


///////////////////////////////////////////////////////////////////////////////
void TRobolFoam::Finalize()
{
//   Finalize MC  run, final printouts, cleaning etc., xcheck of normalization
//   Plotting histograms is done independently using root file
BXOPE(*f_Out);
BXTXT(*f_Out,"########################################");
BXTXT(*f_Out,"###### TRobolFoam::Finalize Start #####");
BXTXT(*f_Out,"########################################");
cout<<       "########################################"<<endl;
cout<<       "###### TRobolFoam::Finalize Start #####"<<endl;
cout<<       "########################################"<<endl;
//
/////////////////////////////////////////////////////////////
  Double_t MCresult, MCerror, MCnorm, Errel;
  KKeeFoam *MCgen = (KKeeFoam*)f_MCgen;
  cout << "**************************************************************"<<endl;
  cout << "**************** TRobolFOAM::Finalize  ***********************"<<endl;
  MCgen->Finalize();

  cout << " ###### Some crosschecks: direct evaluation using local <wt> ###### " <<endl;
  double AveWt, ErrAbs;
  //mon_WtFoam->GetAver(AveWt, ErrAbs);
  double XsNormPb = MCgen->f_FoamI->GetPrimary();
  double xSecPb,xErrPb;
  xSecPb   = XsNormPb*AveWt;
  xErrPb   = XsNormPb*ErrAbs;
  cout << " TRobolFoam: Foam xSec [pb] = "<<  xSecPb << "  +-  "<< xErrPb <<endl;
  cout << "**************************************************************"<<endl;
//
}//Finalize







