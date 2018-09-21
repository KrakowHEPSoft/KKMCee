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
  TMCgenFOAM *MCgen = (TMCgenFOAM*)f_MCgen;
  int jmax = MCgen->m_jmax;
  if( MCgen->m_IsFoam5 == 1) {
    TH1D *HST_FOAM_NORMA5 = (TH1D*)HstFile->Get("HST_FOAM_NORMA5");
    for(int j=1; j<=jmax; j++) HST_FOAM_NORMA5->SetBinContent(j,  MCgen->m_xpar[j]  );    // xpar encoded
    HST_FOAM_NORMA5->SetEntries(0);
  }
  if( MCgen->m_IsFoam3 == 1) {
    for(int j=1; j<=jmax; j++) HST_FOAM_NORMA3->SetBinContent(j,  MCgen->m_xpar[j]  );    // xpar encoded
    HST_FOAM_NORMA3->SetEntries(0); // Important!!!
  }
  if( MCgen->m_IsFoam3i == 1) {
    for(int j=1; j<=jmax; j++) HST_FOAM_NORMA3i->SetBinContent(j,  MCgen->m_xpar[j]  );    // xpar encoded
    HST_FOAM_NORMA3i->SetEntries(0); // Important!!!
  }
  if( MCgen->m_IsFoam1 == 1) {
    for(int j=1; j<=jmax; j++) HST_FOAM_NORMA1->SetBinContent(j,  MCgen->m_xpar[j]  );    // xpar encoded
    HST_FOAM_NORMA1->SetEntries(0);
  }
  if( MCgen->m_IsFoam2 == 1) {
    for(int j=1; j<=jmax; j++) HST_FOAM_NORMA2->SetBinContent(j,  MCgen->m_xpar[j]  );    // xpar encoded
    HST_FOAM_NORMA2->SetEntries(0);
  }
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
    hst_weight1 = TH1D_UP("hst_weight1" ,  "MC weight",      100, -1.5, 2.0);
    hst_weight3 = TH1D_UP("hst_weight3" ,  "MC weight",      100, -1.5, 2.0);
    hst_weight5 = TH1D_UP("hst_weight5" ,  "MC weight",      100, -1.5, 2.0);

    int nbv = 50;

    HST_xx_Ord1n  = TH1D_UP("HST_xx_Ord1n" ,  "dSig/dv",    nbv, 0.0, 1.0);
    HST_xx_Crd1n  = TH1D_UP("HST_xx_Crd1n" ,  "dSig/dv",   nbv, 0.0, 1.0);
    HST_xx_Ord1   = TH1D_UP("HST_xx_Ord1" ,   "dSig/dv",   nbv, 0.0, 1.0);
    HST_xx_Crd1   = TH1D_UP("HST_xx_Crd1" ,   "dSig/dv",   nbv, 0.0, 1.0);
    HST_xx_Hrd1   = TH1D_UP("HST_xx_Hrd1" ,   "dSig/dv",   nbv, 0.0, 1.0);
    HST_xx_Srd1   = TH1D_UP("HST_xx_Srd1" ,   "dSig/dv",   nbv, 0.0, 1.0);
    HST_xx_Ird1   = TH1D_UP("HST_xx_Ird1" ,   "dSig/dv",   nbv, 0.0, 1.0);

    HST_xx_Ceex2  = TH1D_UP("HST_xx_Ceex2" ,   "dSig/dv",   nbv, 0.0, 1.0);
    HST_xx_Ceex2n = TH1D_UP("HST_xx_Ceex2n" ,  "dSig/dv",   nbv, 0.0, 1.0);
    nbv = 50;
// scatergrams
    int nbc = 50;
    SCA_xc_Ceex2  = TH2D_UP("SCA_xc_Ceex2",   "dSig/dc/dv ", nbv, 0.0 ,1.0, nbc, -1.0 ,1.0);
    SCA_xc_Ceex2n = TH2D_UP("SCA_xc_Ceex2n",  "dSig/dc/dv ", nbv, 0.0 ,1.0, nbc, -1.0 ,1.0);
//
    SCA_xc_Ceex0  = TH2D_UP("SCA_xc_Ceex0",   "dSig/dc/dv ", nbv, 0.0 ,1.0, nbc, -1.0 ,1.0);
    SCA_xc_Ceex0n = TH2D_UP("SCA_xc_Ceex0n",  "dSig/dc/dv ", nbv, 0.0 ,1.0, nbc, -1.0 ,1.0);
//
//  New bigger scatergrams, restricted vmax<0.20
    int NBv =100; int NBc = 100;
    double vmx2= 0.20;
    SCT_xc_Ceex2=  TH2D_UP("SCT_xc_Ceex2",  "dSig/dc/dv ", NBv, 0.0 ,vmx2, NBc, -1.0 ,1.0);
    SCT_xc_Ceex2n= TH2D_UP("SCT_xc_Ceex2n", "dSig/dc/dv ", NBv, 0.0 ,vmx2, NBc, -1.0 ,1.0);
//
    SCT_xc_Ceex0=  TH2D_UP("SCT_xc_Ceex0",  "dSig/dc/dv ", NBv, 0.0 ,vmx2, NBc, -1.0 ,1.0);
    SCT_xc_Ceex0n= TH2D_UP("SCT_xc_Ceex0n", "dSig/dc/dv ", NBv, 0.0 ,vmx2, NBc, -1.0 ,1.0);
// EEX series ------------------
    SCT_xc_EEX2 =  TH2D_UP("SCT_xc_EEX2",   "dSig/dc/dv ", NBv, 0.0 ,vmx2, NBc, -1.0 ,1.0);
    SCT_xc_EEX0 =  TH2D_UP("SCT_xc_EEX0",   "dSig/dc/dv ", NBv, 0.0 ,vmx2, NBc, -1.0 ,1.0);
// Adding IFI to ISR, vmax<0.20
    SCT_xc_EEX2i =  TH2D_UP("SCT_xc_EEX2i", "dSig/dc/dv ", NBv, 0.0 ,vmx2, NBc, -1.0 ,1.0);
    SCT_xc_EEX2n =  TH2D_UP("SCT_xc_EEX2n", "dSig/dc/dv ", NBv, 0.0 ,vmx2, NBc, -1.0 ,1.0);
//  The same but restricted vmax<0.02
    double vmx3= 0.02;
    SCN_xc_EEX2  =  TH2D_UP("SCN_xc_EEX2",  "dSig/dc/dv ", NBv, 0.0 ,vmx3, NBc, -1.0 ,1.0);
    SCN_xc_EEX0  =  TH2D_UP("SCN_xc_EEX0",  "dSig/dc/dv ", NBv, 0.0 ,vmx3, NBc, -1.0 ,1.0);
// Adding IFI to ISR, vmax<0.02
    SCN_xc_EEX2i =  TH2D_UP("SCN_xc_EEX2i", "dSig/dc/dv ", NBv, 0.0 ,vmx3, NBc, -1.0 ,1.0);
    SCN_xc_EEX2n =  TH2D_UP("SCN_xc_EEX2n", "dSig/dc/dv ", NBv, 0.0 ,vmx3, NBc, -1.0 ,1.0);
/////////////////////////////////////////////////////////////////////
// scattergram in log10(v)
    m_vminL = 1e-6;
    int NBlv =50, NBlc=100;
    SCT_Lxc_EEX2i =  TH2D_UP("SCT_Lxc_EEX2i", "dSig/dc/dv ", NBlv, log10(m_vminL),-1.0, NBlc, -1.0 ,1.0);
    SCT_Lxc_EEX2n =  TH2D_UP("SCT_Lxc_EEX2n", "dSig/dc/dv ", NBlv, log10(m_vminL),-1.0, NBlc, -1.0 ,1.0);
/////////////////////////////////////////////////////////////////////
    // special histos for AFB from average cos(theta) for IFI on
    HST5_xx_Ceex2  = TH1D_UP("HST5_xx_Ceex2", "dSig/dv",   NBv, 0.0, vmx2);
    HST5_xc_Ceex2  = TH1D_UP("HST5_xc_Ceex2", "c*dSig/dv", NBv, 0.0, vmx2);
    // for AFB=(F-B)/(F+B) with |cos(theta)|<1.0
    HST5_xx_forw_Ceex2  = TH1D_UP("HST5_xx_forw_Ceex2", "dSig/dv",    NBv, 0.0, vmx2);
    // for AFB=(F-B)/(F+B) with |cos(theta)|<0.9
    HST5_xx9_Ceex2      = TH1D_UP("HST5_xx9_Ceex2", "dSig/dv",        NBv, 0.0, vmx2);
    HST5_xx9_forw_Ceex2 = TH1D_UP("HST5_xx9_forw_Ceex2", "dSig/dv",   NBv, 0.0, vmx2);
    //
    HST_cc_EEX2_vmax02   = TH1D_UP("HST_cc_EEX2_vmax02",  "dSig/dc",   NBc, -1, 1);
    HST_cc_EEX2_vmax002  = TH1D_UP("HST_cc_EEX2_vmax002", "dSig/dc",   NBc, -1, 1);
    HST_cc_EEX2_vmax0002 = TH1D_UP("HST_cc_EEX2_vmax0002","dSig/dc",   NBc, -1, 1);
    //
    HST_cc_ceex2_vmax02   = TH1D_UP("HST_cc_ceex2_vmax02",  "dSig/dc",   NBc, -1, 1);
    HST_cc_ceex2_vmax002  = TH1D_UP("HST_cc_ceex2_vmax002", "dSig/dc",   NBc, -1, 1);
    HST_cc_ceex2_vmax0002 = TH1D_UP("HST_cc_ceex2_vmax0002","dSig/dc",   NBc, -1, 1);

    // Testing soft limit
    HST_cs_EEX2_vmax02     = TH1D_UP("HST_cs_EEX2_vmax02",  "dSig/dc",   NBc, -1, 1);
    HST_cs_EEX2_vmax002    = TH1D_UP("HST_cs_EEX2_vmax002", "dSig/dc",   NBc, -1, 1);
    HST_cs_EEX2_vmax0002   = TH1D_UP("HST_cs_EEX2_vmax0002","dSig/dc",   NBc, -1, 1);
    //
    HST_cs_ceex2n_vmax02   = TH1D_UP("HST_cs_ceex2n_vmax02",  "dSig/dc",   NBc, -1, 1);
    HST_cs_ceex2n_vmax002  = TH1D_UP("HST_cs_ceex2n_vmax002", "dSig/dc",   NBc, -1, 1);
    HST_cs_ceex2n_vmax0002 = TH1D_UP("HST_cs_ceex2n_vmax0002","dSig/dc",   NBc, -1, 1);
    //
    HST_cs_ceex2_vmax02   = TH1D_UP("HST_cs_ceex2_vmax02",  "dSig/dc",   NBc, -1, 1);
    HST_cs_ceex2_vmax002  = TH1D_UP("HST_cs_ceex2_vmax002", "dSig/dc",   NBc, -1, 1);
    HST_cs_ceex2_vmax0002 = TH1D_UP("HST_cs_ceex2_vmax0002","dSig/dc",   NBc, -1, 1);

    //************* special normalization histos  *************
    int jmax = ((TMCgenFOAM*)f_MCgen)->m_jmax;
    HST_FOAM_NORMA3 = TH1D_UP("HST_FOAM_NORMA3","Normalization and xpar",jmax,0.0,10000.0);
    HST_FOAM_NORMA3i= TH1D_UP("HST_FOAM_NORMA3i","Normalization and xpar",jmax,0.0,10000.0);
    HST_FOAM_NORMA1 = TH1D_UP("HST_FOAM_NORMA1","Normalization and xpar",jmax,0.0,10000.0);
    HST_FOAM_NORMA2 = TH1D_UP("HST_FOAM_NORMA2","Normalization and xpar",jmax,0.0,10000.0);

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
  double wt1, wt2, wt3, wt5, xx, vv, CosTheta;
  TMCgenFOAM *MCgen = (TMCgenFOAM*)f_MCgen;

/////////////////////////////////////////////////////////////
if( MCgen->m_IsFoam5 == 1) {
  /// MC generation in base class, ISR+FSR+IFI event
  MCgen->m_Mode = -5;
  TRobol::Production(iEvent);  // It invokes MCgen->Generale
  /// filling in histos
  MCgen->f_FoamI->GetMCwt(wt5);
  xx  = MCgen->m_xx;
  CosTheta = MCgen->m_CosTheta;
  hst_weight5->Fill(wt5,1.0);
  SCA_xc_Ceex2->Fill(xx,CosTheta,wt5);
  SCT_xc_Ceex2->Fill(xx,CosTheta,wt5);
  ///
  HST5_xx_Ceex2->Fill(xx,wt5);
  HST5_xc_Ceex2->Fill(xx,wt5*CosTheta);
  // forward
  if( CosTheta> 0.0)                  HST5_xx_forw_Ceex2->Fill(xx,wt5);
  // tot. and forw. for costheta<0.9
  if( fabs(CosTheta) < 0.9 )          HST5_xx9_Ceex2->Fill(xx,wt5);
  if( CosTheta>0.0 && CosTheta < 0.9) HST5_xx9_forw_Ceex2->Fill(xx,wt5);
  //
  double WTceex0 = wt5 * MCgen->m_WTmodel[50];
  SCA_xc_Ceex0->Fill(xx,CosTheta, WTceex0);
  SCT_xc_Ceex0->Fill(xx,CosTheta, WTceex0);
  //
  if( xx < 0.02 )   HST_cc_ceex2_vmax02->Fill(  CosTheta,wt5);
  if( xx < 0.002 )  HST_cc_ceex2_vmax002->Fill( CosTheta,wt5);
  if( xx < 0.0002 ) HST_cc_ceex2_vmax0002->Fill(CosTheta,wt5);

}// m_IsFoam
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
if( MCgen->m_IsFoam3 == 1) {
  /// MC generation in user class, ISR+FSR no IFI
  MCgen->m_Mode = -3;
  MCgen->m_Foam3->MakeEvent();         // Additional Foam of the user class
  // filling in histos
  MCgen->m_Foam3->GetMCwt(wt3);
  double WTeex2, WTeex0, WTceex2n, WTceex0n;
  WTeex2   = wt3;
  WTeex0   = wt3 *MCgen->m_WTmodel[ 3];
  WTceex2n = wt3 *MCgen->m_WTmodel[ 2];
  WTceex0n = wt3 *MCgen->m_WTmodel[52];
//WTceex0n = wt3 *MCgen->m_WTmodel[ 3] *MCgen->m_WTmodel[ 2]; // the same
  xx  = MCgen->m_xx;
  CosTheta = MCgen->m_CosTheta;
  hst_weight3->Fill(wt3,1.0);
  HST_xx_Ceex2n->Fill(xx, WTceex2n);
  HST_tx_Ceex2n->Fill(xx, WTceex2n);
  SCA_xc_Ceex2n->Fill(xx,CosTheta, WTceex2n);
  SCT_xc_Ceex2n->Fill(xx,CosTheta, WTceex2n);
  SCT_xc_EEX2->Fill(xx,CosTheta,WTeex2);
  SCT_xc_EEX0->Fill(xx,CosTheta,WTeex0);
  //
  SCA_xc_Ceex0n->Fill(xx,CosTheta,WTceex0n);
  SCT_xc_Ceex0n->Fill(xx,CosTheta,WTceex0n);
  //
  if( xx < 0.02 )   HST_cc_EEX2_vmax02->Fill(CosTheta,WTeex2);
  if( xx < 0.002 )  HST_cc_EEX2_vmax002->Fill(CosTheta,WTeex2);
  if( xx < 0.0002 ) HST_cc_EEX2_vmax0002->Fill(CosTheta,WTeex2);
  ///  Fill in special normalization histogram
  double Xnorm3 = MCgen->m_Xsav3;
  HST_FOAM_NORMA3->Fill(-1, Xnorm3);      // Normal*Nevtot, new style
}// m_IsFoam
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
if( MCgen->m_IsFoam3i == 1) {
  /// MC generation in user class, ISR+FSR no IFI
  MCgen->m_Mode = -31;
  MCgen->m_Foam3i->MakeEvent();         // Additional Foam of the user class
  xx  = MCgen->m_xx;
  vv  = MCgen->m_vv;
  CosTheta = MCgen->m_CosTheta;
  // filling in histos
  double wtmain;
  MCgen->m_Foam3i->GetMCwt(wtmain);
  double WTeex2i, WTeex0i, WTeex2n, WTeex0n;
  WTeex2i   = wtmain;                       // O(alf2) IFI on
  WTeex0i   = wtmain *MCgen->m_WTmodel[ 3]; // O(alf2) IFI on
  WTeex2n   = wtmain *MCgen->m_WTmodel[ 6]; // O(alf2) IFI off
  WTeex0n   = wtmain *MCgen->m_WTmodel[ 7]; // O(alf2) IFI off
  SCT_xc_EEX2i->Fill(xx,CosTheta,WTeex2i);  // xmax<0.20
  SCT_xc_EEX2n->Fill(xx,CosTheta,WTeex2n);  // xmax<0.20
  //
  SCT_Lxc_EEX2i->Fill(log10(xx +1.001*m_vminL),CosTheta,WTeex2i);  // log scale
  SCT_Lxc_EEX2n->Fill(log10(xx +1.001*m_vminL),CosTheta,WTeex2n);  // log scale
  //
  SCN_xc_EEX2i->Fill(xx,CosTheta,WTeex2i);  // xmax<0.02
  SCN_xc_EEX2n->Fill(xx,CosTheta,WTeex2n);  // xmax<0.02
  ///  Fill in special normalization histogram
  double Xnorm3i = MCgen->m_Xsav3i;
  HST_FOAM_NORMA3i->Fill(-1, Xnorm3i);      // Normal*Nevtot, new style
}// m_IsFoam
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
if( MCgen->m_IsFoam1 == 1) {
/// MC generation in user class, ISR+FSR 1st Order
  MCgen->m_Mode = -1;
  MCgen->m_Foam1->MakeEvent();         // Additional Foam of the user class
  // filling in histos
  MCgen->m_Foam1->GetMCwt(wt1);
  hst_weight1->Fill(wt1,1.0);
  //
  xx  = MCgen->m_vv;
  HST_xx_Ord1n->Fill(xx, wt1);
  HST_xx_Crd1n->Fill(xx, wt1 * MCgen->m_WTmodel[10] );
  HST_xx_Ord1->Fill( xx, wt1 * MCgen->m_WTmodel[11] );
  HST_xx_Crd1->Fill( xx, wt1 * MCgen->m_WTmodel[12] );
  HST_xx_Hrd1->Fill( xx, wt1 * MCgen->m_WTmodel[22] );
  HST_xx_Srd1->Fill( xx, wt1 * MCgen->m_WTmodel[23] );
  HST_xx_Ird1->Fill( xx, wt1 * MCgen->m_WTmodel[24] );

  double Xnorm1 = MCgen->m_Xsav1;
  HST_FOAM_NORMA1->Fill(-1, Xnorm1);      // Normal*Nevtot, new style
}// m_IsFoam

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
if( MCgen->m_IsFoam2 == 1) {
/// MC generation in user class, ISR+FSR+IFI soft limit
  MCgen->m_Mode = -2;
  MCgen->m_Foam2->MakeEvent();         // Additional Foam of the user class
  MCgen->m_Foam2->GetMCwt(wt2);

  CosTheta = MCgen->m_CosTheta;
  //
  HST_cs_EEX2_vmax02->Fill(  CosTheta,wt2);
  HST_cs_EEX2_vmax002->Fill( CosTheta,wt2*MCgen->m_WTmodel[72]);
  HST_cs_EEX2_vmax0002->Fill(CosTheta,wt2*MCgen->m_WTmodel[73]);
  //
  HST_cs_ceex2n_vmax02->Fill(  CosTheta,wt2*MCgen->m_WTmodel[74]);
  HST_cs_ceex2n_vmax002->Fill( CosTheta,wt2*MCgen->m_WTmodel[75]);
  HST_cs_ceex2n_vmax0002->Fill(CosTheta,wt2*MCgen->m_WTmodel[76]);
  //
  HST_cs_ceex2_vmax02->Fill(  CosTheta,wt2*MCgen->m_WTmodel[77]);
  HST_cs_ceex2_vmax002->Fill( CosTheta,wt2*MCgen->m_WTmodel[78]);
  HST_cs_ceex2_vmax0002->Fill(CosTheta,wt2*MCgen->m_WTmodel[79]);

  double Xnorm2 = MCgen->m_Xsav2;
  HST_FOAM_NORMA2->Fill(-1, Xnorm2);   // Normal*Nevtot, new style
}// m_IsFoam

///
}///Production

///______________________________________________________________________________
void TRobolFOAM::Finalize()
{
/////////////////////////////////////////////////////////////
///------------------------

  Double_t MCresult, MCerror, MCnorm, Errel;
  TMCgenFOAM *MCgen = (TMCgenFOAM*)f_MCgen;
  cout << "**************************************************************"<<endl;
  cout << "**************** TRobolFOAM::Finalize  ***********************"<<endl;
  if( MCgen->m_IsFoam5 == 1) {
    MCgen->Finalize();
  }// m_IsFoam
  if( MCgen->m_IsFoam3 == 1) {
    MCgen->m_Foam3->Finalize(MCnorm, Errel);  // Additional Foam of the user class
    MCgen->m_Foam3->GetIntegMC( MCresult, MCerror);  //! get MC integral
  }// m_IsFoam
  if( MCgen->m_IsFoam3i == 1) {
    MCgen->m_Foam3i->Finalize(MCnorm, Errel);  // Additional Foam of the user class
    MCgen->m_Foam3i->GetIntegMC( MCresult, MCerror);  //! get MC integral
  }// m_IsFoam
  if( MCgen->m_IsFoam1 == 1) {
    MCgen->m_Foam1->Finalize(MCnorm, Errel);  // Additional Foam of the user class
    MCgen->m_Foam1->GetIntegMC( MCresult, MCerror);  //! get MC integral
    cout << "Directly from m_Foam1: MCresult= " << MCresult << " +- "<<MCerror <<endl;
 }// m_IsFoam
   cout << "**************************************************************"<<endl;

}
