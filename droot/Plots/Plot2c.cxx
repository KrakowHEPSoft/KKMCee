//////////////////////////////////////////////////////////////////////
//    make Plot2c
//////////////////////////////////////////////////////////////////////
#include <iomanip.h>
#include <fstream.h>
#include <math.h>
#include <TLorentzVector.h>
#include <TLine.h>
#include <TArrow.h>
#include <TLatex.h>
#include "TROOT.h"
#include "TCanvas.h"
#include "TF2.h"
#include "TH2.h"
#include "TGaxis.h"
#include "TApplication.h"
#include "TMarker.h"

// ROOT headers
#include "TROOT.h"
#include "TFile.h"
//=============================================================================
//  ROOT  ROOT ROOT   ROOT  ROOT  ROOT  ROOT  ROOT  ROOT  ROOT   ROOT   ROOT 
//=============================================================================
TFile DiskFileA("../demoC/rmain.root");
//TFile DiskFileA("../demoC/rmain.root.EEX.100M");
TFile DiskFileB("RhoSemi.root","RECREATE","histograms");
//=============================================================================

// Auxiliary procedures for plotting
#include "HisNorm.h"
#include "Marker.h"

Double_t sqr( const Double_t x ){ return x*x;};

///////////////////////////////////////////////////////////////////////////////
//
extern "C" void kk2f_fort_open_( const long&, char*, long);
extern "C" void kk2f_fort_close_(const long&);
//      SUBROUTINE KK2f_GetKeyFSR(KeyFSR)
extern "C" void kk2f_getkeyfsr_( long&);
//      SUBROUTINE KK2f_Initialize(xpar)
extern "C" void kk2f_initialize_(double xpar[]);
//      DOUBLE PRECISION FUNCTION RRes_CorQQ(roots,i,xmi)
extern "C" double rres_corqq_( const double&, const long&, const double&);
//      DOUBLE PRECISION FUNCTION RRes_R_TOT(roots)
extern "C" double rres_r_tot_(const double&);
//      DOUBLE PRECISION FUNCTION RRes_Rqq(kf,roots)
extern "C" double rres_rqq_(const long&, const double&);
//      DOUBLE PRECISION FUNCTION RRes_F_PI_SQ(s)
extern "C" double rres_f_pi_sq_(const double&);
//      DOUBLE PRECISION  FUNCTION F_Pi_Kuehn90_SQ(s)
extern "C" double rres_f_pi_kuehn90_sq_(const double&);
//      DOUBLE PRECISION  FUNCTION F_Pi_Kuehn02_SQ(s)
extern "C" double rres_f_pi_kuehn02_sq_(const double&);
//      SUBROUTINE KKsem_Initialize(xpar_input)
extern "C" void kksem_initialize_(double xpar[]);
//      SUBROUTINE KKsem_VVplot_vec(key,chak,nbin,xmin,xmax,yy)
extern "C" void kksem_vvplot_vec_(const long&, char[5], const long&, const double&, const double&, double[]);
//      SUBROUTINE KKsem_SetKeyFoB(KeyFoB)
extern "C" void kksem_setkeyfob_( const long& );
//      SUBROUTINE KKsem_SetKFfin(KFfin)
extern "C" void kksem_setkffin_( const long& );
//      SUBROUTINE BornV_SetKF(KFferm)
extern "C" void bornv_setkf_( const long& ); // set SINGLE Final State
//      SUBROUTINE BornV_SetIsGenerated(KFferm,IsGenerated)
extern "C" void bornv_setisgenerated_( const long&, const long&); // togle final state flavor generation

///////////////////////////////////////////////////////////////////////////////////
void HistNormalize(){
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  //
  HisNorm(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vvMuTrg0") );
  HisNorm(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vvMuTrg1") );
  //
  HisNorm(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_Q2MuTrg0") );
  HisNorm(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_Q2MuTrg1") );
  //
  HisNorm(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_Q2All") );
  HisNorm(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_Q2Trig0") );
  HisNorm(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_Q2Trig1") );
  HisNorm(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_Q2Trig2") );
  
}

///////////////////////////////////////////////////////////////////////////////////
void KKsem_initialize(){
  //------------------------------------------------------------------------
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  const int jmax =10000;
  double ypar[jmax];
  for(int j=1; j<=jmax; j++)
    ypar[j-1]=HST_KKMC_NORMA->GetBinContent(j);    // xpar encoded
  //[[[[[
  //for(int j=0;j<30;j++)
  //  cout<<j<<"   "<<ypar[j]<<endl;
  //]]]]]
  char *output_file = "./kksem.output";
  long stl2 = strlen(output_file);
  long mout =16;
  kk2f_fort_open_(mout,output_file,stl2);
  kk2f_initialize_(ypar);
  kksem_initialize_(ypar);
  //long kdum; kk2f_getkeyfsr_(kdum); //<-- to avoid linker bug! (if no kk2f_initialize)
}

///////////////////////////////////////////////////////////////////////////////////
void KKsem_make(){
//------------------------------------------------------------------------  
//-------------------------KKsem------------------------------------------  
// Here we produce semianalytical plots using KKsem program, No plotting
//------------------------------------------------------------------------  
  long Key, KeyFob;
  char chak[5]="VRHO2";
  Key = 302;   // O(alf2)
  Key = 304;   // O(alf3) GribovLL
  Key = 303;   // O(alf3)
  Key = 305;   // O(alf3) GribovLL +NLL
  //Key = 662;  // Unexp ???? 
  KeyFob=   10; // BornV_Dizet
  KeyFob=  -10; // KKsem_BornV
  KeyFob=  -11; // BornV_Simple for KeyLib=0
  KeyFob= -100; // KKsem_BornV, without EW, with integration
  long IsGenerated;
  long KF=13;
  //
  Double_t CMSene;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  double s= CMSene*CMSene;
  //
  TH1D *hst_Q2MuTrg0 = (TH1D*)DiskFileA.Get("hst_Q2MuTrg0");
  int     Nb = hst_Q2MuTrg0->GetNbinsX();
  double qmin= hst_Q2MuTrg0->GetXaxis()->GetXmin();
  double qmax= hst_Q2MuTrg0->GetXaxis()->GetXmax();
  cout<<"|||||||||| Nb= "<<Nb<<"  qmin, qmax =  "<<qmin<<"  "<< qmax<<endl;
  double vmin=1-qmax/s;
  double vmax=1-qmin/s;
  cout<<"||||||||||   vmin, vmax =  "<<vmin<<"  "<< vmax<<endl;
  double v,Q2,Jacob;
  double vbin[Nb];
//------------------------------------------------------------------------
//  PiPi Q^2dsig/dQ^2
//------------------------------------------------------------------------
  KF=13; IsGenerated=0; bornv_setisgenerated_( KF, IsGenerated); // exclude muon !!!
  KF=1; kksem_setkffin_(KF); // set one d-quark for special tests only
  Jacob = 1/s; // this is dv/dQ^2
  KeyFob= -100; // KKsem_BornV, DOUBLY KF-summed! WRONG!!!! CosTh-integr.(*10)
  KeyFob=   10; // BornV_Dizet(1,CosTh=0), one m_KFfin only, EW alowed   (*1/3) ????
  KeyFob=  -11; // BornV_Simple(CosTh=0),  one m_KFfin only, KeyLib=0 (R=1/3= NC*(1/3)^2) OK
  KeyFob=  -10; // KKsem_BornV(CosTh=0),   KF-sumed, (R=6/3=2*NC*(1/3)^2+NC*(2/3)^2) OK!!!!
  kksem_setkeyfob_(KeyFob);
  kksem_vvplot_vec_(Key,chak,Nb,vmin,vmax,vbin);
  TH1D *hst_Q2PiOrd3  = new TH1D("hst_Q2PiOrd3","PiPi, Q^2dsig/dQ^2",  Nb, qmin,qmax);
  hst_Q2PiOrd3->Sumw2();
  hst_Q2PiOrd3->SetLineColor(9); // blue
  double dv=(vmax-vmin)/Nb;
  double sigTag=0.0;
  for(int ib=0;ib<Nb;ib++){
    v=vmin+(vmax-vmin)*(ib+0.0)/Nb;
    Q2= s*(1-v);
    cout<<"**** ib= "<<ib<<" v= "<<v<<" Q2= "<<Q2<<" vbin(ib)= "<<vbin[ib]<<endl;
    hst_Q2PiOrd3->SetBinContent(Nb-ib, vbin[ib]*Jacob); // bin content
    hst_Q2PiOrd3->SetBinError(Nb-ib,0.0);
    sigTag += vbin[ib]*dv;
  }
  cout<<"%%%%%%%%%%  sigTag [nb]= "<<sigTag<<endl;
//------------------------------------------------------------------------
//  MuMu Q^2dsig/dQ^2
//------------------------------------------------------------------------  
  Jacob = 1/s; // this is dv/dQ^2
  KF=13; bornv_setkf_( KF ); // set SINGLE Final State
  KF=13; kksem_setkffin_(KF); // set m_KFfin in KKsem
  KeyFob=   10; // BornV_Dizet
  KeyFob=  -11; // BornV_Simple for KeyLib=0
  KeyFob= -100; // KKsem_BornV, without EW, with integration, OK
  KeyFob=  -10; // KKsem_BornV, OK!
  kksem_setkeyfob_(KeyFob);
  kksem_vvplot_vec_(Key,chak,Nb,vmin,vmax,vbin);
  TH1D *hst_Q2MuOrd3  = new TH1D("hst_Q2MuOrd3","MuMu  Q^2dsig/dQ^2",  Nb, qmin,qmax);
  hst_Q2MuOrd3->Sumw2();
  hst_Q2MuOrd3->SetLineColor(8); // green
  for(int ib=0;ib<Nb;ib++){
    v=vmin+(vmax-vmin)*(ib+0.5)/Nb;
    Q2= s*(1-v);
    //cout<<"**** ib= "<<ib<<" v= "<<v<<" vbin(ib)= "<<vbin[ib]<<endl;
    hst_Q2MuOrd3->SetBinContent(Nb-ib, vbin[ib]*Jacob); // bin content
    hst_Q2MuOrd3->SetBinError(Nb-ib,0.0);
  }
//------------------------------------------------------------------------  
//------------------------------------------------------------------------  
}  //KKsem_make

///////////////////////////////////////////////////////////////////////////////////
Double_t Rtot(Double_t *xvar, Double_t *par)
// R total from all quarks
{
  long      KFfin  = par[0];
  Double_t  Mqq  = xvar[0];
  Double_t  Result;
  Result = rres_r_tot_(Mqq);
  return Result;
}
///////////////////////////////////////////////////////////////////////////////////
Double_t Rqq(Double_t *xvar, Double_t *par)
// R from one quark pair
{
  long      KFf  = par[0];
  Double_t  Mqq;
  Double_t  Result =0.0;
  Mqq  = xvar[0]; // Mass in GeV units!!!!
  if(KFf>0 && KFf<6){
    Result = rres_rqq_(KFf,Mqq);
  }else if(KFf==12){
    Result = rres_rqq_(1,Mqq)+rres_rqq_(2,Mqq);
  }else if(KFf==123){
    Result = rres_rqq_(1,Mqq)+rres_rqq_(2,Mqq)+rres_rqq_(3,Mqq);
  }else if(KFf==1234){
    Result = rres_rqq_(1,Mqq)+rres_rqq_(2,Mqq)+rres_rqq_(3,Mqq)+rres_rqq_(4,Mqq);
  }else if(KFf==12345){
    Result = rres_rqq_(1,Mqq)+rres_rqq_(2,Mqq)+rres_rqq_(3,Mqq)
            +rres_rqq_(4,Mqq)+rres_rqq_(5,Mqq);
  }else if(KFf==99){
    Mqq  = sqrt(xvar[0]); //!!!!! s->M translation!!!!
    Result = rres_rqq_(1,Mqq)+rres_rqq_(2,Mqq)+rres_rqq_(3,Mqq);
  }else if(KFf==3001){
    Result = rres_f_pi_sq_(Mqq);          // f_pi Novosibirsk 99
  }else if(KFf==3002){
    Result = rres_f_pi_kuehn90_sq_(Mqq);  // f_pi Kuehn 90
  }else if(KFf==3003){
    Result = rres_f_pi_kuehn02_sq_(Mqq);  // f_pi Kuehn 90
  }else if(KFf==3010){
    Mqq  = sqrt(xvar[0]);
    double mass=0.100;
    long   KF  =1;
    Result =  rres_corqq_(Mqq,KF,mass); // Correction to R(s)
  }else{
    cout<<" STOP in Rqq"<<endl;
    exit(13);
  }
  return Result;
}


///////////////////////////////////////////////////////////////////////////////////
int DefFpi()
{
  cout<<" ========================= DefFpi =========================== "<<endl;
  //------------------------------------------------------------------------
  // parameters for ploting function
  Int_t npar =1;
  Double_t params[9];
  params[0]=1;       // type of quark
  Double_t s1,s2;
  int nbin;
  //------------------------------------------------------------------------
  nbin = 430; s1=0.100; s2=0.960; // rho region KLOE
  //----- R(s) function
  TF1  *funRtot  = new  TF1("funRtot",Rqq, s1,s2,npar);
  params[0]=99;       // type of quark dus
  funRtot->SetParameters(params);
  TH1F *hisRtot  = new TH1F("hisRtot","",  nbin, s1,s2);
  hisRtot->Eval(funRtot,"a");
  //----- Two f_pi function
  TF1  *funFpiM  = new  TF1("funFpiM",Rqq, s1,s2,npar);
  params[0]=3001;       // Fpi^2 of Marteen, Novosibirsk 99
  funFpiM->SetParameters(params);
  TH1F *hisFpiM  = new TH1F("hisFpiM","",  nbin, s1,s2);
  hisFpiM->Eval(funFpiM,"a");
  //-----
  TF1  *funFpiK90  = new  TF1("funFpiK90",Rqq, s1,s2,npar);
  params[0]=3002;       // Fpi^2 Kuehn 1990
  funFpiK90->SetParameters(params);
  TH1F *hisFpiK90  = new TH1F("hisFpiK90","",  nbin, s1,s2);
  hisFpiK90->Eval(funFpiK90,"a");
  //-----
  TF1  *funFpiK02  = new  TF1("funFpiK02",Rqq, s1,s2,npar);
  params[0]=3003;       // Fpi^2 Kuehn 2002
  funFpiK02->SetParameters(params);
  TH1F *hisFpiK02  = new TH1F("hisFpiK02","",  nbin, s1,s2);
  hisFpiK02->Eval(funFpiK02,"a");
  //----- Correction to R(s), for single KF defined in Rqq
  TF1  *funRcoR  = new  TF1("funRcoR",Rqq, s1,s2,npar);
  params[0]=3010;       // KF defined in Rqq
  funRcoR->SetParameters(params);
  TH1F *hisRcoR  = new TH1F("hisRcoR","",  nbin, s1,s2);
  hisRcoR->Eval(funRcoR,"a");
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  return 0;
}


///////////////////////////////////////////////////////////////////////////////////
int ReadHist(TH1D *Hst, char DiskFile[])
{
// Reads histogram from the input file
// and multiplies by the x-value, middle of the bin
  cout<<"===========ReadHist==========ReadHist==============ReadHist============="<<endl;
  ifstream InputFile;
  cout<<"===  "<< DiskFile <<" =========="<<endl;
  InputFile.open(DiskFile);
  int nb = Hst->GetNbinsX();
  cout<<nb<<endl;
  Double_t QQ,QQi,dsig,ddsig;
  for(int i=1; i<nb+1; i++){
    InputFile >>QQi >>dsig >>ddsig;
    QQ = Hst->GetBinCenter(i);
    cout<<i <<"  "<<QQ <<"  "<< dsig<<"  "<< ddsig <<endl;
    Hst->SetBinContent(i,dsig );
    Hst->SetBinError(  i,ddsig);
  }
  InputFile.close();
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  return 0;
}

///////////////////////////////////////////////////////////////////////////////////
int Fig2c()
{
//------------------------------------------------------------------------
  cout<<" ========================= Fig2c =========================== "<<endl;

  TH1D *H_Ste1 = new TH1D("H_Ste1","Q2 distr. Born, cut", 43, 0.10, 0.96);
  TH1D *H_Ste2 = new TH1D("H_Ste2","Q2 distr. Born, all",43, 0.10, 0.96);
  TH1D *H_Ste3 = new TH1D("H_Ste3","Q2 distr. NLO , cut", 43, 0.10, 0.96);
  TH1D *H_Ste4 = new TH1D("H_Ste4","Q2 distr. NLO , all", 43, 0.10, 0.96);
  ReadHist(H_Ste1, "./ds_dqq_phokhara_g_10e6_born_5_21_55_125.dat");
  ReadHist(H_Ste2, "./ds_dqq_phokhara_g_10e6_born_0_180_0_180.dat");
  ReadHist(H_Ste3, "./ds_dqq_phokhara_g_10e6_nlo_5_21_55_125.dat");
  ReadHist(H_Ste4, "./ds_dqq_phokhara_g_10e6_nlo_0_180_0_180.dat");
  H_Ste1->Scale(2.0); // error in trigger corrected
  H_Ste3->Scale(2.0);
  TH1D *H_axel1 = new TH1D("H_axel1","Q2 distr. LO , all", 43, 0.10, 0.96);
  ReadHist(H_axel1, "./axel_nocuts.dat");
  //
  Double_t CMSene;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");

  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  CMSene  *=1000;
  //
  Float_t msize=0.6; // marker size for postscript
  //  R(s) semianalytical
  TH1D *hisRtot = (TH1D*)DiskFileB.Get("hisRtot");
  // 2pi+3pi
  TH1D *hst_Q2All = (TH1D*)DiskFileA.Get("hst_Q2All");
  //  2pi  Q^2 dsig/dQ^2,
  TH1D *hst_Q2Trig0 = (TH1D*)DiskFileA.Get("hst_Q2Trig0");
  TH1D *hst_Q2Trig1 = (TH1D*)DiskFileA.Get("hst_Q2Trig1");
  TH1D *hst_Q2Trig2 = (TH1D*)DiskFileA.Get("hst_Q2Trig2");
  BlackBullet( hst_Q2All, msize);
  RedTriangle( hst_Q2Trig1, msize);  // marker etc
  hst_Q2Trig1->SetLineWidth(2.5);
  hst_Q2All->SetLineWidth(2.5);
  hst_Q2All->SetMinimum(0.0);

  TH1D *hst_Q2MuTrg0 = (TH1D*)DiskFileA.Get("hst_Q2MuTrg0");
  TH1D *hst_Q2MuTrg1 = (TH1D*)DiskFileA.Get("hst_Q2MuTrg1");
  //
  //==================================================================================
  // top caption
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextAlign(23);
  CaptT->SetTextSize(0.050);
  char capt1[100];
  sprintf(capt1,"KKMC: %4.0fMeV", CMSene);
  // bottom caption
  TLatex *CaptB = new TLatex();
  CaptB->SetNDC(); // !!!
  CaptB->SetTextAlign(21);
  CaptB->SetTextSize(0.045);
  //
  TH1D *hst_Q2MuOrd3 = (TH1D*)DiskFileB.Get("hst_Q2MuOrd3");
  TH1D *hst_Q2PiOrd3 = (TH1D*)DiskFileB.Get("hst_Q2PiOrd3");
  //
  ///////////////////////////////////////////////////////////////////////////////
  Float_t  WidPix, HeiPix;
  WidPix = 1000; HeiPix =  900;
  TCanvas *cFig2c = new TCanvas("cFig2c","Fig2c photonic", 70, 20, WidPix,HeiPix);
  cFig2c->SetFillColor(10);
  cFig2c->cd();
  //void Divide(Int_t nx, Int_t ny, Float_t xmargin, Float_t ymargin, Int_t color)
  cFig2c->Divide(2, 2, 0.0, 0);
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  cFig2c->cd(1);
  hst_Q2Trig0->SetStats(0);
  hst_Q2Trig0->SetTitle(0);
  hst_Q2Trig0->DrawCopy("h");   // <-- 2pi
  //
  hst_Q2MuTrg0->SetLineColor(8);  // green
  hst_Q2MuTrg0->DrawCopy("hsame"); // <-- MUON
  //hst_Q2MuOrd3->DrawCopy("hsame");
  // Phokara
  H_Ste2->SetLineColor(6);  // magenta
  H_Ste2->DrawCopy("hsame");
  H_Ste4->SetLineColor(2);  // red
  H_Ste4->DrawCopy("hsame");
  // axel
  H_axel1->SetLineColor(7);  // cyan
  H_axel1->DrawCopy("hsame");

  CaptT->DrawLatex(0.50,0.98, "2#pi KKMC&Phokara, no cuts");
  CaptB->DrawLatex(0.50,0.01,"Q^{2} [GeV]");
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  cFig2c->cd(2);
  TH1D *hisRatio0 = (TH1D*)hst_Q2Trig0->Clone("hisRatio0");
  hisRatio0->Divide(H_Ste4);
  hisRatio0->SetMaximum(1.15);
  hisRatio0->SetMinimum(0.85);
  hisRatio0->SetStats(0);
  hisRatio0->SetTitle(0);
  hisRatio0->SetLineColor(2);     // Red
  hisRatio0->DrawCopy("h");       // NLO
  TH1D *hisRatio1 = (TH1D*)hst_Q2Trig0->Clone("hisRatio1");
  hisRatio1->Divide(H_Ste2);
  hisRatio1->SetLineColor(6);     // Magenta
  hisRatio1->DrawCopy("hsame");   // Born
  TH1D *hisRatio7 = (TH1D*)hst_Q2Trig0->Clone("hisRatio7");
  hisRatio7->Divide(H_axel1);
  hisRatio7->SetLineColor(7);     // cyan
  hisRatio7->DrawCopy("hsame");   // axel
  CaptT->DrawLatex(0.50,0.98, "KKMC/other");
  CaptB->DrawLatex(0.50,0.01,"Q^{2} [GeV]");
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  cFig2c->cd(3);
  hst_Q2Trig2->SetStats(0);
  hst_Q2Trig2->SetTitle(0);
  hst_Q2Trig2->DrawCopy("h");
  H_Ste1->SetLineColor(6);  // Magenta
  H_Ste1->DrawCopy("hsame");// Born
  H_Ste3->SetLineColor(2);  // Red
  H_Ste3->DrawCopy("hsame");// NLO
  CaptT->DrawLatex(0.50,0.98, "2#pi KKMC&Phokara, with cuts");
  CaptB->DrawLatex(0.50,0.01,"Q^{2} [GeV]");
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  cFig2c->cd(4);
  TH1D *hisRatio3 = (TH1D*)hst_Q2Trig2->Clone("hisRatio3");
  hisRatio3->Divide(H_Ste3);
  hisRatio3->SetMaximum(1.15);
  hisRatio3->SetMinimum(0.85);
  hisRatio3->SetLineColor(2);  // Red
  hisRatio3->DrawCopy("h");
  TH1D *hisRatio4 = (TH1D*)hst_Q2Trig2->Clone("hisRatio4");
  hisRatio4->Divide(H_Ste1);
  hisRatio4->SetLineColor(6);     // Magenta
  hisRatio4->DrawCopy("hsame");   // Born
  CaptT->DrawLatex(0.50,0.98, "KKMC/other");
  CaptB->DrawLatex(0.50,0.01,"Q^{2} [GeV]");
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  cFig2c->cd();
  return 0;
}



///////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  //++++++++++++++++++++++++++++++++++++++++
  TApplication theApp("theApp", &argc, argv);
  //++++++++++++++++++++++++++++++++++++++++
  HistNormalize();     // Renormalization of MC histograms
  KKsem_initialize();
  KKsem_make();
  DefFpi();
  Fig2c();
  //++++++++++++++++++++++++++++++++++++++++
  DiskFileA.ls();
  DiskFileB.ls();
  DiskFileB.Write();
  DiskFileB.Close();
  //++++++++++++++++++++++++++++++++++++++++
  theApp.Run();
  //++++++++++++++++++++++++++++++++++++++++
}
