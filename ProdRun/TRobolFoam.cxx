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
// Name of special histogram name "HST_FOAM_NORMA4" is set in Start.C
  TH1D *h_TMCgen_NORMA = (TH1D*)HstFile->Get("HST_FOAM_NORMA4");
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
  double WtMax = 2.0;
  int NBwt =200;
 // mon_WtFoam = new THwtMon("mon_WtFoam", "WtFoam",  NBwt, WtMax);
//  mon_WtMain = new THwtMon("mon_WtMain", "WtMain",  NBwt, WtMax);
//
  HST_weight4  = TH1D_UP("HST_weight4" , "MC weight",    NBwt,   0.0 , WtMax);
  HST_weight6  = TH1D_UP("HST_weight6" , "MC weight",    NBwt,   0.0 , WtMax);

  int nBins = 100;
//  Hst_Mll   = TH1D_UP("Hst_Mll",   "Dilepton Invariant Mass with IFI;M_{ll};d#sigma/dM_{Z} (nb/GeV)", nBins, m_Mllmin, m_Mllmax);

  HST_vv_eex2 = TH1D_UP("HST_vv_eex2",   "vv", nBins, 0.0, 1.0);

  int nbv =100;
  int nbc = 50;

//  SCA_vTcPR_Ceex2  = TH2D_UP("SCA_vTcPR_Ceex2",   "dSig/dc/dv ", nbv, 0.0 ,1.0, nbc, -1.0 ,1.0);
//  SCA_vTcPR_Ceex2n = TH2D_UP("SCA_vTcPR_Ceex2n",  "dSig/dc/dv ", nbv, 0.0 ,1.0, nbc, -1.0 ,1.0);
  SCA_vTcPR_Eex2   = TH2D_UP("SCA_vTcPR_Eex2",    "dSig/dc/dv ", nbv, 0.0 ,1.0, nbc, -1.0 ,1.0);

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
  KKevent *Event = MCgen->m_Event;
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// MC generation using FOAM of the base class
//////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////
  MCgen->m_FoamMode = -4;
  TRobol::Production(iEvent);  // It invokes MCgen->Generate() !!!
/// filling in histos
  double WTfoam;
  MCgen->f_FoamI->GetMCwt(WTfoam);
  xx  = MCgen->m_xx;
  //
  CosTheta = MCgen->m_CosTheta;
  vv  = MCgen->m_vv;

  HST_weight4->Fill(WTfoam,1.0);

  double WtEEX2 = MCgen->m_WtAlter[73];
//  double WtEEX3 = MCgen->m_WtAlter[74]; // not implemented
  double WtEEX0 = MCgen->m_WtAlter[71];
////  WtEEX2= WtEEX0; // !!!!!!!!!!!!!! DEBUG

  HST_vv_eex2->Fill(xx, WtEEX2);

// big scatergrams, range vv< 1.0
//  SCA_vTcPR_Ceex2->Fill(   vv, CosPRD, WtCEEX2);
//  SCA_vTcPR_Ceex2n->Fill(  vv, CosPRD, WtCEEX2n); // true v, IFI off
  SCA_vTcPR_Eex2->Fill(    xx, CosTheta, WtEEX2);

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
// MC generation using additional FOAM object
/////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
  MCgen->m_FoamMode = -6;
  MCgen->Generate();               //!!!!!!
  MCgen->m_Foam6->GetMCwt(WTfoam);
// event has changed!!!
  HST_weight6->Fill(WTfoam,1.0);

  vv  = Event->m_vv;

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







