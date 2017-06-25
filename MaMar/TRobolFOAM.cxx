#include "TRobolFOAM.h"

////////////////////////////////////////////////////////////////////////////////
/// MC generator class for Two Parton exercises

ClassImp(TRobolFOAM);

TRobolFOAM::TRobolFOAM():
  TRobol()
{
  /// This constructor is for ROOT streamers ONLY
  cout<< "@@@@> TRobolFOAM DEFAULT Constructor (for ROOT only) "<<endl;
}

///_____________________________________________________________
TRobolFOAM::TRobolFOAM(const char* Name):
  TRobol(Name)
{
//! Constructor to be used by the user!!!
//! Its important role is to define ALL DEFAULTS.
//! to changed by the user before calling TMCgen::Initialize
  cout<< "@@@@> TRobolFOAM::TRobolFOAM USER Constructor "<<endl;
  m_xmin    = 0.01;        //! x range
  m_xmax    = 0.99;        //! x range

}///

///______________________________________________________________________________________
TRobolFOAM::~TRobolFOAM()
{
  //!Explicit destructor
  cout<< "@@@@> TRobolFOAM::TRobolFOAM !!!! DESTRUCTOR !!!! "<<endl;
}///destructor

//______________________________________________________________________________
void TRobolFOAM::Initialize(
        ofstream *OutFile, /// Central log-file for messages
        TFile *GenFile,    /// ROOT disk file for CRNG and MC gen.
        TFile *HstFile)    /// ROOT disk file for histograms
{
  cout<< "****> TRobolFOAM::Initialize starts"<<endl;
  TRobol::Initialize(OutFile,GenFile,HstFile);
  ///
  /// Book histograms or read them from the disk
  Hbooker();
  //  dumping KKMC info into normalization histo
  TH1D *HST_FOAM_NORMA5 = (TH1D*)HstFile->Get("HST_FOAM_NORMA5");
  int jmax = ((TMCgenFOAM*)f_MCgen)->m_jmax;
  for(int j=1; j<=jmax; j++) {
    HST_FOAM_NORMA3->SetBinContent(j,  ((TMCgenFOAM*)f_MCgen)->m_xpar[j]  );    // xpar encoded
    HST_FOAM_NORMA5->SetBinContent(j,  ((TMCgenFOAM*)f_MCgen)->m_xpar[j]  );    // xpar encoded
  }
  HST_FOAM_NORMA3->SetEntries(0); // Important!!!
  HST_FOAM_NORMA5->SetEntries(0);

  cout<< "****> TRobolFOAM::Initialize: finished"<<endl;
}///Initialize

///______________________________________________________________________________
void TRobolFOAM::Hbooker()
{
  ///
  cout<< "****> TRobolFOAM::Hbooker: histogram booking STARTS"<<endl;
  BXOPE(*f_Out);
  BXTXT(*f_Out,"========================================");
  BXTXT(*f_Out,"======    TRobol::Hbooker    ===========");
  BXTXT(*f_Out,"========================================");
  BXCLO(*f_Out);
  f_HstFile->cd();

/////////////////////////////////////////////////////////////////////////////////////////
//  ************* user histograms  *************
    hst_weight3 = TH1D_UP("hst_weight3" ,  "MC weight",      100, -1.0, 2.0);
    hst_weight5 = TH1D_UP("hst_weight5" ,  "MC weight",      100, -1.0, 2.0);

    int nbv = 50;
    HST_xx_Ceex2  = TH1D_UP("HST_xx_Ceex2" ,   "dSig/dv",   nbv, 0.0, 1.0);
    HST_xx_Ceex2n = TH1D_UP("HST_xx_Ceex2n" ,  "dSig/dv",   nbv, 0.0, 1.0);
    nbv = 50;
    // scatergrams
    int nbc = 50;
    SCA_xc_Ceex2  = TH2D_UP("SCA_xc_Ceex2",   "dSig/dc/dv ", nbv, 0.0 ,1.0, nbc, -1.0 ,1.0);
    SCA_xc_Ceex2n = TH2D_UP("SCA_xc_Ceex2n",  "dSig/dc/dv ", nbv, 0.0 ,1.0, nbc, -1.0 ,1.0);

    //  New bigger scatergrams, restricted vmax
    int NBv =100; int NBc = 100;
    double vmx2= 0.20;
    SCT_xc_Ceex2=  TH2D_UP("SCT_xc_Ceex2",  "dSig/dc/dv ", NBv, 0.0 ,vmx2, NBc, -1.0 ,1.0);
    SCT_xc_Ceex2n= TH2D_UP("SCT_xc_Ceex2n", "dSig/dc/dv ", NBv, 0.0 ,vmx2, NBc, -1.0 ,1.0);
    SCT_xc_EEX2 =  TH2D_UP("SCT_xc_EEX2",   "dSig/dc/dv ", NBv, 0.0 ,vmx2, NBc, -1.0 ,1.0);

    //************* special normalization histos  *************
    int jmax = ((TMCgenFOAM*)f_MCgen)->m_jmax;
    HST_FOAM_NORMA3 = TH1D_UP("HST_FOAM_NORMA3","Normalization and xpar",jmax,0.0,10000.0);

    // additional histo for testing normalization
    HST_tx_Ceex2n = TH1D_UP("HST_tx_Ceex2n" ,  "dSig/dv",   100, 0.0, 0.2);

///////////////////////////////////////////////////////////////////////////////////////////
  BXOPE(*f_Out);
  BXTXT(*f_Out,"========================================");
  BXTXT(*f_Out,"====== END of TRobol::Hbooker ==========");
  BXTXT(*f_Out,"========================================");
  BXCLO(*f_Out);
  cout<< "****> TRobolFOAM::Hbooker: histogram booking FINISHED"<<endl;
}///Hbooker

///______________________________________________________________________________
void TRobolFOAM::Production(double &iEvent)
{
/////////////////////////////////////////////////////////////
  double wt3, wt5, xx, CosTheta;
/// MC generation in base class, ISR+FSR+IFI event
  TMCgenFOAM *MCgen = (TMCgenFOAM*)f_MCgen;
  MCgen->m_Mode = -5;
  TRobol::Production(iEvent);  // It invokes TMCgen->TMCgenFOAM
  /// filling in histos
  MCgen->f_FoamI->GetMCwt(wt5);
  xx  = MCgen->m_xx;
  CosTheta = MCgen->m_CosTheta;
  hst_weight5->Fill(wt5,1.0);
  SCA_xc_Ceex2->Fill(xx,CosTheta,wt5);
  SCT_xc_Ceex2->Fill(xx,CosTheta,wt5);
//???  HST_FOAM_NORMA5->Fill(-1.0,m_Xsav5);  // fill normalization into underflow

  /////////////////////////////////////////////////////////////
  /// MC generation in user class, ISR+FSR no IFI
  MCgen->m_Mode = -3;
  MCgen->m_Foam3->MakeEvent();         // Additional Foam of the user class
  // filling in histos
  MCgen->m_Foam3->GetMCwt(wt3);
  xx  = MCgen->m_xx;
  CosTheta = MCgen->m_CosTheta;
  hst_weight3->Fill(wt3,1.0);
  HST_xx_Ceex2n->Fill(xx,wt3);
  HST_tx_Ceex2n->Fill(xx,wt3);
  SCA_xc_Ceex2n->Fill(xx,CosTheta,wt3);
  SCT_xc_Ceex2n->Fill(xx,CosTheta,wt3);
  double WTeex2 = wt3 * MCgen->m_WTmodel[2];
  SCT_xc_EEX2->Fill(xx,CosTheta,WTeex2);

  ///  Fill in special normalization histogram
  double Xnorm = MCgen->m_Xsav3;
  HST_FOAM_NORMA3->Fill(-1, Xnorm);      // Normal*Nevtot, new style
  //HST_FOAM_NORMA3->Fill(0.5, Xnorm);
  //HST_FOAM_NORMA3->Fill(1.5, Xnorm);      // Normal*Nevtot
///
}///Production

///______________________________________________________________________________
void TRobolFOAM::Finalize()
{
/////////////////////////////////////////////////////////////
///------------------------
  TMCgenFOAM *MCgen = (TMCgenFOAM*)f_MCgen;
  MCgen->Finalize();

  Double_t MCresult, MCerror, MCnorm, Errel;
  MCgen->m_Foam3->Finalize(MCnorm, Errel);  // Additional Foam of the user class
  MCgen->m_Foam3->GetIntegMC( MCresult, MCerror);  //! get MC integral
  cout << "**************************************************************"<<endl;
  cout << "**************** TRobolFOAM::Finalize  ***********************"<<endl;
  cout << "Directly from FOAM: MCresult= " << MCresult << " +- "<<MCerror <<endl;
  cout << "**************************************************************"<<endl;

}
