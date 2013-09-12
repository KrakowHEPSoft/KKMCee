//////////////////////////////////////////////////////////////////////
//    make Plot2b
//    almost obsolete
//////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;


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
//TFile DiskFileA("../demoC/rmain.root");
TFile DiskFileA("../demoC/rmain.root.EEX.960M.Ord1"); // does not make much sense!!!
//TFile DiskFileA("../demoC/rmain.root.EEX.500M");
//TFile DiskFileA("../demoC/rmain.root.180M");
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
  HisNorm(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_Q2Trig1") );
  
}

///////////////////////////////////////////////////////////////////////////////////
void KKsem_initialize(){
  //------------------------------------------------------------------------
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  const int jmax =10000;
  double ypar[jmax];
  double Nrun =HST_KKMC_NORMA->GetBinContent(511); // shoul be one for single run
  for(int j=1; j<=jmax; j++){
    ypar[j-1]=HST_KKMC_NORMA->GetBinContent(j)/Nrun; //ypar imported, corrected
    HST_KKMC_NORMA->SetBinContent(j,ypar[j-1]);      //corrected
  }
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
  char chak[6]="VRHO2";
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
//   MuMu  dsigma/dv
//------------------------------------------------------------------------  
  long nbin = 50;
  double  xmin = 0.00;
  double  xmax = 1.00;
  double  bin[nbin];
  bornv_setkf_( KF ); // set SINGLE Final State
  kksem_setkeyfob_(KeyFob);
  kksem_vvplot_vec_(Key,chak,nbin,xmin,xmax,bin);
  //
  TH1D *hst_vvMuOrd3  = new TH1D("hst_vvMuOrd3","MuMu  dsigma/dv ",  nbin, xmin,xmax);
  hst_vvMuOrd3->Sumw2();
  hst_vvMuOrd3->SetLineColor(8); // green
  for(int ib=0;ib<nbin;ib++){
    //cout<<" ib= "<<ib<<"  bin(ib) =  "<<bin[ib]<<endl;
    hst_vvMuOrd3->SetBinContent(ib+1,bin[ib]);
    hst_vvMuOrd3->SetBinError(ib+1,0.0);
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
  Double_t  Mqq  = xvar[0];
  Double_t  Result =0.0;
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
    Mqq  = sqrt(xvar[0]);
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
    //exit(13);
  }
  return Result;
}


///////////////////////////////////////////////////////////////////////////////////
void figFpi()
{
  cout<<" ========================= figFpi =========================== "<<endl;
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
  params[0]=99;       // type of quark duscb
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
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  Float_t  WidPix, HeiPix;
  WidPix = 1100; HeiPix =  600;
  TCanvas *cKloe = new TCanvas("cKloe","Resonances",200,180, WidPix,HeiPix);
  cKloe->Divide(2, 1, 0.0,  0);
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  ///////////////////////////////////////////////////////////////////////////////
  cKloe->SetFillColor(10);
  cKloe->cd(1);
  //
  hisRtot->GetYaxis()->SetTitle("R");
  hisRtot->GetYaxis()->CenterTitle();
  hisRtot->GetYaxis()->SetTitleSize(0.07);
  hisRtot->GetYaxis()->SetTitleOffset(0.5);
  hisRtot->GetXaxis()->SetTitle("Q^{2}  [GeV^{2}]");
  hisRtot->GetXaxis()->CenterTitle();
  hisRtot->GetXaxis()->SetTitleSize(0.05);
  // Kloeal duscb
  //hisRtot->SetMaximum(22.0);
  hisRtot->SetStats(0);
  hisRtot->SetLineColor(1);    // black line
  //hisRtot->SetFillColor(3);    // green fill solid (psi)
  //hisRtot->SetFillColor(4);    // blue fill solid (upsyl)
  //hisRtot->SetFillColor(2);    // red fill solid (rho)
  //
  hisRtot->DrawCopy("");
  //funRcoR->DrawCopy("");
  ///////////////////////////////////////////////////////////////////////////////
  cKloe->cd(2);
  funFpiM->DrawCopy("");
  funFpiK90->SetLineColor(2);    // red line
  funFpiK90->DrawCopy("same");
  funFpiK02->SetLineColor(3);    // red line
  funFpiK02->DrawCopy("same");
  // top caption
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextAlign(23);
  CaptT->SetTextSize(0.050);
  CaptT->DrawLatex(0.50,0.98, "|F_{#pi}(Q^{2})|^{2}");
  ///////////////////////////////////////////////////////////////////////////////
  double s=0.55;
  double fpi1= rres_f_pi_sq_(s);
  double fpi2= rres_f_pi_kuehn90_sq_(s);
  cout<<" fpi_marteen = " <<fpi1<<"  fpi_kuehn = "<<fpi2<<endl;
  ///////////////////////////////////////////////////////////////////////////////
  cKloe->Update();
  cKloe->cd();
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
}


///////////////////////////////////////////////////////////////////////////////////
void Fig2b()
{
//------------------------------------------------------------------------
  cout<<" ========================= Fig2b =========================== "<<endl;
  Double_t CMSene;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");

  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  CMSene  *=1000;
  //
  Float_t msize=0.6; // marker size for postscript
  //  R(s) semianalytical
  TH1D *hisRtot = (TH1D*)DiskFileB.Get("hisRtot");
  //  2pi  Q^2 dsig/dQ^2,
  TH1D *hst_Q2All = (TH1D*)DiskFileA.Get("hst_Q2All");
  TH1D *hst_Q2Trig1 = (TH1D*)DiskFileA.Get("hst_Q2Trig1");
  BlackBullet( hst_Q2All, msize);
  RedTriangle( hst_Q2Trig1, msize);  // marker etc
  //hst_Q2Trig1->SetLineWidth(2.5);
  //hst_Q2All->SetLineWidth(2.5);
  hst_Q2All->SetMinimum(0.0);
  // 2mu   Q^2 dsig/dQ^2,
  TH1D *hst_Q2MuTrg0 = (TH1D*)DiskFileA.Get("hst_Q2MuTrg0");
  TH1D *hst_Q2MuTrg1 = (TH1D*)DiskFileA.Get("hst_Q2MuTrg1");
  //BlackBullet( hst_Q2MuTrg0, msize);
  //RedTriangle( hst_Q2MuTrg1, msize);  // marker etc
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
  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  Float_t  WidPix, HeiPix;
  WidPix = 1200; HeiPix =  600;
  TCanvas *cFig2b1 = new TCanvas("cFig2b1","Fig2b1 photonic", 40, 00, WidPix,HeiPix);
  cFig2b1->SetFillColor(10);
  cFig2b1->Divide(3, 0, 0.0, 0);
  ///////////////////////////////////////////////////////////////////////////////
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  cFig2b1->cd(1);
  hst_Q2All->SetStats(0);
  hst_Q2All->SetTitle(0);
  hst_Q2All->DrawCopy("h");
  hst_Q2MuTrg0->DrawCopy("hsame");   // <-- MUON
  hst_Q2MuOrd3->DrawCopy("hsame");
  hst_Q2PiOrd3->DrawCopy("hsame");
  CaptT->DrawLatex(0.50,0.98, "2#pi+3#pi versus 2 #mu (inclusive) ");
  CaptB->DrawLatex(0.50,0.01,"Q^{2} [GeV]");
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  cFig2b1->cd(2);
  hisRtot->DrawCopy("h");
  hisRtot->SetMinimum(0.0);
  TH1D *hisRatio0 = (TH1D*)hst_Q2All->Clone("hisRatio0");
  hisRatio0->Divide(hst_Q2MuTrg0);    // <-----  TAKE Ratio0 R= pi/mu
  hisRatio0->DrawCopy("hsame");
  CaptT->DrawLatex(0.50,0.98, "The ratio d#sigma_{had} /d#sigma_{#mu}");
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||  
  cFig2b1->cd(3);
  hst_Q2All->SetStats(0);
  hst_Q2All->Divide(hst_Q2PiOrd3);
  hst_Q2All->SetMaximum(1.03);
  hst_Q2All->SetMinimum(0.93);
  hst_Q2All->SetLineColor(9); // blue
  hst_Q2All->DrawCopy("h");
  hst_Q2MuTrg0->Divide(hst_Q2MuOrd3);
  hst_Q2MuTrg0->SetMinimum(0.9);
  hst_Q2MuTrg0->SetLineColor(8); // green
  hst_Q2MuTrg0->DrawCopy("hsame");
  CaptT->DrawLatex(0.50,0.98, "The ratio KKMC/KKsem");
  CaptB->DrawLatex(0.50,0.01,"Q^{2} [GeV]");
  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  WidPix = 800; HeiPix =  600;
  TCanvas *cFig2b2 = new TCanvas("cFig2b2","Fig2b2 photonic", 80, 40, WidPix,HeiPix);
  cFig2b2->SetFillColor(10);
  cFig2b2->Divide(2, 0, 0.0, 0);
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||  
  cFig2b2->cd(1);
  hst_Q2Trig1->SetStats(0);
  hst_Q2Trig1->SetTitle(0);
  hst_Q2Trig1->DrawCopy("h");
  hst_Q2MuTrg1->DrawCopy("hsame");   // <-- MUON
  CaptT->DrawLatex(0.50,0.98, "2#pi photon ang. cut ");
  CaptB->DrawLatex(0.50,0.01,"Q^{2} [GeV]");
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  cFig2b2->cd(2);
  hisRtot->DrawCopy("h");
  TH1D *hisRatio1 = (TH1D*)hst_Q2Trig1->Clone("hisRatio1");
  hisRatio1->Divide(hst_Q2MuTrg1);    // <-----  TAKE RATIO R= pi/mu
  hisRatio1->DrawCopy("hsame");
  CaptT->DrawLatex(0.50,0.98, "The ratio d#sigma_{had} /d#sigma_{#mu}");
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  cFig2b2->cd();
}

///////////////////////////////////////////////////////////////////////////////////
void Fig2mu()
{
//------------------------------------------------------------------------
  cout<<" ========================= Fig2mu =========================== "<<endl;
  // renormalize histograms in nanobarns
  Double_t CMSene;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");

  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  CMSene  *=1000;
  // distr. of Q^2dSig/dQ^2
  TH1D *hst_Q2MuTrg0 = (TH1D*)DiskFileA.Get("hst_Q2MuTrg0");
  TH1D *hst_Q2MuTrg1 = (TH1D*)DiskFileA.Get("hst_Q2MuTrg1");
  // distr. of dSig/dv
  TH1D *hst_vvMuTrg0 = (TH1D*)DiskFileA.Get("hst_vvMuTrg0");
  TH1D *hst_vvMuTrg1 = (TH1D*)DiskFileA.Get("hst_vvMuTrg1");
  //
  TH1D *hst_vvMuOrd3 = (TH1D*)DiskFileB.Get("hst_vvMuOrd3");
  TH1D *hst_Q2MuOrd3 = (TH1D*)DiskFileB.Get("hst_Q2MuOrd3");
//------------------------------------------------------------------------  
  ///////////////////////////////////////////////////////////////////////////////
  Float_t  WidPix, HeiPix;
  WidPix = 1000; HeiPix =  800;
  TCanvas *cFig2mu = new TCanvas("cFig2mu","Fig2b photonic2", 50, 80, WidPix,HeiPix);
  cFig2mu->SetFillColor(10);
  cFig2mu->Divide(2, 2, 0.0, 0);
  //==========plot1==============
  cFig2mu->cd(1);
  gPad->SetLogy(); // !!!!!!
  hst_vvMuTrg0->SetStats(0);
  hst_vvMuTrg0->DrawCopy();
  hst_vvMuOrd3->DrawCopy("same");
  //==========plot2==============
  cFig2mu->cd(2);
  hst_vvMuTrg0->Divide(hst_vvMuOrd3);
  hst_vvMuTrg0->SetMinimum(0.95);
  hst_vvMuTrg0->SetMaximum(1.05);
  hst_vvMuTrg0->DrawCopy();
  //==========plot3==============
  cFig2mu->cd(3);
  hst_Q2MuTrg0->SetStats(0);
  hst_Q2MuTrg0->DrawCopy("");
  hst_Q2MuOrd3->DrawCopy("hsame");
  //==========plot4==============
  cFig2mu->cd(4);
  hst_Q2MuOrd3->Divide(hst_Q2MuTrg0);
  //hst_Q2MuOrd3->SetMinimum(0.95);
  //hst_Q2MuOrd3->SetMaximum(1.05);
  hst_Q2MuOrd3->SetStats(0);
  hst_Q2MuOrd3->DrawCopy("");
  //----------------------------
  cFig2mu->cd();
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
  figFpi();  // <--- cannot comment out this!!!
  Fig2b();
  Fig2mu();
  //++++++++++++++++++++++++++++++++++++++++
  DiskFileA.ls();
  DiskFileB.ls();
  DiskFileB.Write();
  DiskFileB.Close();
  //++++++++++++++++++++++++++++++++++++++++
  theApp.Run();
  return 0;
  //++++++++++++++++++++++++++++++++++++++++
}
