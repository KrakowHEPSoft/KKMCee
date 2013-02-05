//////////////////////////////////////////////////////////////////////
//    make Plot3f
//    muons, FSR ON
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
//TFile DiskFileA("~/sprivate/M.C./KK-all/RHadr/demoC/rmain.root");
//TFile DiskFileA("../demoC/rmain.root");
TFile DiskFileA("../demoC/rmain.root.MuCEEX.FSR.175M"); // last one
//TFile DiskFileA("../demoC/rmain.root.MuCEEX.FSR.150M");
TFile DiskFileB("RhoSemi.root","RECREATE","histograms");
//=============================================================================

// Auxiliary procedures for plotting
#include "HisNorm.h"
#include "Marker.h"

Double_t sqr( const Double_t x ){ return x*x;};

////////////////////////////////////////////////////////////////////////////////
extern "C" void kk2f_fort_open_( const long&, char*, long);
extern "C" void kk2f_fort_close_(const long&);
//      SUBROUTINE KK2f_Initialize(xpar)
extern "C" void kk2f_initialize_(double xpar[]);
//      SUBROUTINE KKsem_Initialize(xpar_input)
extern "C" void kksem_initialize_(double xpar[]);
//      SUBROUTINE BornV_SetKF(KFferm)
extern "C" void bornv_setkf_( const long& ); // set SINGLE Final State
//      SUBROUTINE KKsem_SetKFfin(KFfin)
extern "C" void kksem_setkffin_( const long& );
//      SUBROUTINE KKsem_SetKeyFoB(KeyFoB)
extern "C" void kksem_setkeyfob_( const long& );
//      SUBROUTINE KKsem_VVplot_vec(key,chak,nbin,xmin,xmax,yy)
extern "C" void kksem_vvplot_vec_(const long&, char[5], const long&, const double&, const double&, double[]);
//      SUBROUTINE BornV_SetIsGenerated(KFferm,IsGenerated)
extern "C" void bornv_setisgenerated_( const long&, const long&); // togle final state flavor generation
//      SUBROUTINE BornV_GetIsGenerated(KFferm,IsGenerated)
extern "C" void bornv_getisgenerated_( const long&, const long&); // get final state flavor generation
//////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////
void HistNormalize(){
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  //
  HisNorm(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("Hst_Q2muA") );
  HisNorm(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vvMuTrg0") );
  //
  HisNorm(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("Hst_Vmu01") );
  HisNorm(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("Hst_Vmu71") );
  HisNorm(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("Hst_Vmu72") );
  HisNorm(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("Hst_Vmu73") );
  HisNorm(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("Hst_Vmu74") );

  HisNorm(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("Hst_Vmu10") );
  HisNorm(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("Hst_Vmu11") );
  HisNorm(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("Hst_Vmu20") );
  HisNorm(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("Hst_Vmu21") );
  HisNorm(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("Hst_Vmu22") );
  HisNorm(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("Hst_Vmu202") );
  HisNorm(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("Hst_Vmu203") );
  HisNorm(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("Hst_Vmu252") );
  HisNorm(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("Hst_Vmu253") );
  }


///////////////////////////////////////////////////////////////////////////////////
void KKsem_initialize(){
  //------------------------------------------------------------------------
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  const int jmax =10000;
  double ypar[jmax];
  double Nrun =HST_KKMC_NORMA->GetBinContent(511); // shoul be one for single run
  cout<<"Nrun=   "<<Nrun<<endl;
  for(int j=1; j<=jmax; j++){
    ypar[j-1]=HST_KKMC_NORMA->GetBinContent(j)/Nrun; //ypar imported, corrected
    HST_KKMC_NORMA->SetBinContent(j,ypar[j-1]);      //corrected
  }
  //[[[[[
  //for(int j=0;j<30;j++)
  //  cout<<j<<"<ypar=   "<<ypar[j]<<endl;
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
//-------------------------KKsem------------------------------------------  
// Here we produce semianalytical plots using KKsem program, No plotting
//------------------------------------------------------------------------  
  long KF,Key,KeyX, KeyFob;
  char chak[6]="VRHO2";
  Key = 305;   // O(alf2)
  Key = 302;   // O(alf2)
  Key = 304;   // O(alf3) GribovLL
  Key = 303;   // O(alf3)
  Key = 305;   // O(alf3) GribovLL +NLL
  //----------------------------
  Double_t CMSene;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  double svar= CMSene*CMSene;
  cout<<" CMSene = "<< CMSene<<endl;
  //
  TH1D *Hst_Q2muA = (TH1D*)DiskFileA.Get("Hst_Q2muA");
  int     Nb = Hst_Q2muA->GetNbinsX();
  double qmin= Hst_Q2muA->GetXaxis()->GetXmin();
  double qmax= Hst_Q2muA->GetXaxis()->GetXmax();
  if(qmax > CMSene )  qmax = svar;
  cout<<"|||||||||| Nb= "<<Nb<<"  qmin, qmax =  "<<qmin<<"  "<< qmax<<endl;
  double vmin=1-qmax/svar;
  double vmax=1-qmin/svar;
  cout<<"||||||||||   vmin, vmax =  "<<vmin<<"  "<< vmax<<endl;
  TH1D *hst_Q2line  = new TH1D("hst_Q2line","one",  1, qmin,qmax);
  TH1D *hst_VVline  = new TH1D("hst_VVline","one",  1, 0.0, 1.0);
  double v,Q2,Jacob;
  double vbin[Nb];
//------------------------------------------------------------------------
//  MuMu dsig/dQ^2 ISR+FSR
//------------------------------------------------------------------------  
  cout<<"**** MuMu dsig/dQ^2 **** "<<endl;
  char chakX[6]="XRHO2";
  KeyX = 305302;  // ISR+FSR
  Jacob = 1/svar; // this is dv/dQ^2
  KF=13; bornv_setkf_( KF );  // set SINGLE Final State
  KF=13; kksem_setkffin_(KF); // set m_KFfin in KKsem
  KeyFob=  -10; // KKsem_BornV, OK!
  kksem_setkeyfob_(KeyFob);
  //
  kksem_vvplot_vec_(KeyX,chakX,Nb,vmin,vmax,vbin);
  TH1D *hst_Q2Mu305302  = new TH1D("hst_Q2Mu305302","MuMu  dsig/dQ^2",  Nb, qmin,qmax);
  hst_Q2Mu305302->Sumw2();
  for(int ib=0;ib<Nb;ib++){
    v=vmin+(vmax-vmin)*(ib+0.5)/Nb;
    Q2= svar*(1-v);
    cout<<"+++ ib= "<<ib<<" v= "<<v<<" vbin(ib)= "<<vbin[ib]<<endl;
    hst_Q2Mu305302->SetBinContent(Nb-ib, vbin[ib]*Jacob); // bin content
    hst_Q2Mu305302->SetBinError(Nb-ib,0.0);
  }
//------------------------------------------------------------------------  
// MuMu dsig/dv 
//------------------------------------------------------------------------  
  KeyFob=  -11; // BornV_Simple for KeyLib=0, OK
  KeyFob= -100; // KKsem_BornV, without EW, with integration, OK
  KeyFob=  -10; // KKsem_BornV, OK!
  kksem_setkeyfob_(KeyFob);
  // 74
  TH1D *Hst_Vmu74 = (TH1D*)DiskFileA.Get("Hst_Vmu74");
  Nb =  Hst_Vmu74->GetNbinsX();
  double VVbin[Nb];
  vmin = Hst_Vmu74->GetXaxis()->GetXmin();
  vmax = Hst_Vmu74->GetXaxis()->GetXmax();
  // 305302
  KeyX = 305302;   // O(alf3) GribovLL +NLL
  kksem_vvplot_vec_(KeyX,chakX,Nb,vmin,vmax,VVbin);
  TH1D *hst_vvMu305302  = new TH1D("hst_vvMu305302","MuMu  dsig/dv",  Nb, vmin,vmax);
  hst_vvMu305302->Sumw2();
  for(int ib=0;ib<Nb;ib++){
    cout<<"**** ib= "<<ib<<" v= "<<v<<" VVbin(ib)= "<<VVbin[ib]<<endl;
    hst_vvMu305302->SetBinContent(ib+1, VVbin[ib]); // bin content
    hst_vvMu305302->SetBinError(ib+1,0.0);
 }
  // 301301
  Key = 301301;   // O(alf3) GribovLL +NLL
  kksem_vvplot_vec_(KeyX,chakX,Nb,vmin,vmax,VVbin);
  TH1D *hst_vvMu301301  = new TH1D("hst_vvMu301301","MuMu  dsig/dv",  Nb, vmin,vmax);
  hst_vvMu301301->Sumw2();
  for(int ib=0;ib<Nb;ib++){
    cout<<"**** ib= "<<ib<<" v= "<<v<<" VVbin(ib)= "<<VVbin[ib]<<endl;
    hst_vvMu301301->SetBinContent(ib+1, VVbin[ib]); // bin content
    hst_vvMu301301->SetBinError(ib+1,0.0);
 }
  //  302302
  KeyX = 302302;   // O(alf3) GribovLL +NLL
  kksem_vvplot_vec_(KeyX,chakX,Nb,vmin,vmax,VVbin);
  TH1D *hst_vvMu302302  = new TH1D("hst_vvMu302302","MuMu  dsig/dv",  Nb, vmin,vmax);
  hst_vvMu302302->Sumw2();
  for(int ib=0;ib<Nb;ib++){
    cout<<"**** ib= "<<ib<<" v= "<<v<<" VVbin(ib)= "<<VVbin[ib]<<endl;
    hst_vvMu302302->SetBinContent(ib+1, VVbin[ib]); // bin content
    hst_vvMu302302->SetBinError(ib+1,0.0);
 }
  // 300300
  KeyX = 300300;   // beta0 with NLL
  kksem_vvplot_vec_(KeyX,chakX,Nb,vmin,vmax,VVbin);
  TH1D *hst_vvMu300300  = new TH1D("hst_vvMu300300","MuMu  dsig/dv",  Nb, vmin,vmax);
  hst_vvMu300300->Sumw2();
  for(int ib=0;ib<Nb;ib++){
    cout<<"**** ib= "<<ib<<" v= "<<v<<" VVbin(ib)= "<<VVbin[ib]<<endl;
    hst_vvMu300300->SetBinContent(ib+1, VVbin[ib]); // bin content
    hst_vvMu300300->SetBinError(ib+1,0.0);
 }
  // 400300
  KeyX = 400300;   // beta0 with NLL
  kksem_vvplot_vec_(KeyX,chakX,Nb,vmin,vmax,VVbin);
  TH1D *hst_vvMu400300  = new TH1D("hst_vvMu400300","MuMu  dsig/dv",  Nb, vmin,vmax);
  hst_vvMu400300->Sumw2();
  for(int ib=0;ib<Nb;ib++){
    cout<<"**** ib= "<<ib<<" v= "<<v<<" VVbin(ib)= "<<VVbin[ib]<<endl;
    hst_vvMu400300->SetBinContent(ib+1, VVbin[ib]); // bin content
    hst_vvMu400300->SetBinError(ib+1,0.0);
  }
//------------------------------------------------------------------------  
}  //KKsem_make



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
int Fig3m()
{
//------------------------------------------------------------------------
  cout<<" ========================= Fig3m =========================== "<<endl;

  TH1D *H_Phk1A = new TH1D("H_Phk1A","Q2 distr. Born, incl.", 100, 0.00, 1.0);
  TH1D *H_Phk2A = new TH1D("H_Phk2A","Q2 distr. NLO,  incl.", 100, 0.00, 1.0);
  ReadHist(H_Phk1A, "./phokara_born_1_qq.dat");  // no cut, high stat
  ReadHist(H_Phk2A, "./phokara_nlo_1_qq.dat");   // no cut, high stat
  TH1D *H_Phk1B = new TH1D("H_Phk1B","Q2 distr. Born, cut.", 100, 0.00, 1.0);
  TH1D *H_Phk2B = new TH1D("H_Phk2B","Q2 distr. NLO,  cut.", 100, 0.00, 1.0);
  ReadHist(H_Phk1B, "./phokara_born_2_qq.dat"); // no cut
  ReadHist(H_Phk2B, "./phokara_nlo_2_qq.dat");  // no cut
  //
  TH1D *H_Phk2Amu = new TH1D("H_Phk2Amu","mu Q2 distr. NLO,  incl.", 100, 0.00, 1.0);
  TH1D *H_Phk2Bmu = new TH1D("H_Phk2Bmu","mu Q2 distr. NLO,   cut.", 100, 0.00, 1.0);
  ReadHist(H_Phk2Amu, "./phokhara_mumu_nlo_sel1.dat");  // no cut
  ReadHist(H_Phk2Bmu, "./phokhara_mumu_nlo_sel2.dat");  // cut
  //
  Double_t CMSene;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  // CMSene=xpar(1) stored in NGeISR
  CMSene  = HST_KKMC_NORMA->GetBinContent(1);
  CMSene  *=1000;
  int KeyGPS = HST_KKMC_NORMA->GetBinContent(28);
  //
  Float_t msize=0.6; // marker size for postscript
  TH1D *Hst_Q2muA    = (TH1D*)DiskFileA.Get("Hst_Q2muA");    // no cut
  TH1D *hst_vvMuTrg0 = (TH1D*)DiskFileA.Get("hst_vvMuTrg0"); // no cut
  TH1D *Hst_Vmu01    = (TH1D*)DiskFileA.Get("Hst_Vmu01");    // no cut
  TH1D *Hst_Vmu71    = (TH1D*)DiskFileA.Get("Hst_Vmu71");    // no cut
  TH1D *Hst_Vmu72    = (TH1D*)DiskFileA.Get("Hst_Vmu72");    // no cut
  TH1D *Hst_Vmu73    = (TH1D*)DiskFileA.Get("Hst_Vmu73");    // no cut
  TH1D *Hst_Vmu74    = (TH1D*)DiskFileA.Get("Hst_Vmu74");    // no cut
  TH1D *Hst_Vmu10    = (TH1D*)DiskFileA.Get("Hst_Vmu10");    // no cut
  TH1D *Hst_Vmu11    = (TH1D*)DiskFileA.Get("Hst_Vmu11");    // no cut
  TH1D *Hst_Vmu20    = (TH1D*)DiskFileA.Get("Hst_Vmu20");    // no cut  
  TH1D *Hst_Vmu21    = (TH1D*)DiskFileA.Get("Hst_Vmu21");    // no cut  
  TH1D *Hst_Vmu22    = (TH1D*)DiskFileA.Get("Hst_Vmu22");    // no cut  
  TH1D *Hst_Vmu202   = (TH1D*)DiskFileA.Get("Hst_Vmu202");    // no cut  
  TH1D *Hst_Vmu203   = (TH1D*)DiskFileA.Get("Hst_Vmu203");    // no cut  
  TH1D *Hst_Vmu253   = (TH1D*)DiskFileA.Get("Hst_Vmu253");    // no cut  
  //
  int     Nb = Hst_Q2muA->GetNbinsX();
  cout<<"|||| content of  Hst_Q2muA"<<endl;
  for(int ib=0;ib<Nb;ib++){
    cout<<"|||| ib= "<<ib<<" bin= "<<Hst_Q2muA->GetBinContent(ib+1)
	                 <<" err= "<<Hst_Q2muA->GetBinError(ib+1)<<endl;
 }
  //
  TH1D *hst_Q2MuO2     = (TH1D*)DiskFileB.Get("hst_Q2MuO2");
  TH1D *hst_Q2Mu305302 = (TH1D*)DiskFileB.Get("hst_Q2Mu305302");
  TH1D *hst_vvMu305302 = (TH1D*)DiskFileB.Get("hst_vvMu305302");
  TH1D *hst_vvMu301301 = (TH1D*)DiskFileB.Get("hst_vvMu301301");
  TH1D *hst_vvMu302302 = (TH1D*)DiskFileB.Get("hst_vvMu302302");
  TH1D *hst_vvMu300300 = (TH1D*)DiskFileB.Get("hst_vvMu300300");
  TH1D *hst_vvMu400300 = (TH1D*)DiskFileB.Get("hst_vvMu400300");
  TH1D *hst_Q2line   = (TH1D*)DiskFileB.Get("hst_Q2line");    // no cut
  TH1D *hst_VVline   = (TH1D*)DiskFileB.Get("hst_VVline");    // no cut
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
  // inside caption
  TLatex *CaptC = new TLatex();
  CaptC->SetNDC(); // !!!
  CaptC->SetTextAlign(11);
  CaptC->SetTextSize(0.040);
  //***************************
  TLatex *CaptY = new TLatex();
  CaptY->SetNDC(); CaptY->SetTextAlign(12);
  CaptY->SetTextSize(0.040);
  Float_t verti=0.90;
  Float_t horiz=0.20;
  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  Float_t  WidPix, HeiPix;
  Int_t ndivx =1; 
  WidPix = 600+500*(ndivx-1); HeiPix =  700;
  TCanvas *cFig3f1 = new TCanvas("cFig3f1","Fig3m1 photonic", 20,  0, WidPix,HeiPix);
  cFig3f1->SetFillColor(10);
  cFig3f1->cd();
  //void Divide(Int_t nx, Int_t ny, Float_t xmargin, Float_t ymargin, Int_t color)
  if(ndivx>1) cFig3f1->Divide(ndivx, 0, 0.0, 0);
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  if(ndivx>1) cFig3f1->cd(1);
  hst_Q2line->SetMaximum(1.04);
  hst_Q2line->SetMinimum(0.96);
  hst_Q2line->SetStats(0);
  hst_Q2line->SetTitle(0);
  hst_Q2line->SetBinContent(1,1);
  hst_Q2line->DrawCopy("h");  // horizontal line
  CaptY->SetTextColor(kBlack); verti= 0.85;
  CaptY->DrawLatex(horiz,verti, "(b)   ISR*FSR, mu pair");
  //***
  TH1D *hisRat1 = (TH1D*)hst_Q2Mu305302->Clone("hisRat1");
  hisRat1->Divide(Hst_Q2muA);     // KKsem
  SetBullet( hisRat1, kBlue, 0.7);
  hisRat1->DrawCopy("hsame");
  CaptY->SetTextColor(kBlue); verti-= 0.05;
  CaptY->DrawLatex(horiz,verti, "KKMC.CEEX/KKsem");
  Bullet(kBlue, 1.2, horiz-0.03, verti);
  //***
  CaptB->DrawLatex(0.50,0.01,"Q^{2} [GeV]");
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  if(ndivx>1){
  cFig3f1->cd(2);
  hst_Q2line->SetMaximum(1.02);
  hst_Q2line->SetMinimum(0.98);
  hst_Q2line->SetStats(0);
  hst_Q2line->SetTitle(0);
  hst_Q2line->SetBinContent(1,1);
  hst_Q2line->DrawCopy("h");  // horizontal line
  CaptY->SetTextColor(kBlack); verti= 0.88;
  CaptY->DrawLatex(horiz,verti, "(b) RATIOS");
  //********
  CaptB->DrawLatex(0.50,0.01,"Q^{2} [GeV]");
  }// ndivx
  cFig3f1->cd();
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  return 0;
  ///////////////////////////////////////////////////////////////////////////////
  //Float_t  WidPix, HeiPix;
  ndivx =2;
  WidPix = 600+500*(ndivx-1); HeiPix =  700;
  TCanvas *cFig3f2 = new TCanvas("cFig3f2","Fig3m2 photonic", 40, 20, WidPix,HeiPix);
  cFig3f2->SetFillColor(10);
  cFig3f2->cd();
  cFig3f2->Divide(ndivx, 0, 0.0, 0);
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  cFig3f2->cd(1);
  hst_Q2line->SetMaximum(1.04);
  hst_Q2line->SetMinimum(0.96);
  hst_Q2line->SetStats(0);
  hst_Q2line->SetTitle(0);
  hst_Q2line->SetBinContent(1,1);
  hst_Q2line->DrawCopy("h");  // horizontal line
  CaptY->SetTextColor(kBlack); verti= 0.85;
  CaptY->DrawLatex(horiz,verti, "(a) KKMC/KKsem");
  //***
  TH1D *hisRat74a = (TH1D*)Hst_Vmu74->Clone("hisRat74a");
  hisRat74a->Divide(hst_vvMu305302);
  hisRat74a->SetLineColor(kGreen);
  hisRat74a->DrawCopy("hsame");       // KKMC/KKsem, GREEN
  CaptY->SetTextColor(kGreen); verti-= 0.05;
  CaptY->DrawLatex(horiz,verti, "EEX74/KKsem305302");
  //***
  TH1D *hisRat73a = (TH1D*)Hst_Vmu73->Clone("hisRat73a");
  hisRat73a->Divide(hst_vvMu305302);
  hisRat73a->SetLineColor(kMagenta);
  hisRat73a->DrawCopy("hsame");
  CaptY->SetTextColor(kMagenta); verti-= 0.05;
  CaptY->DrawLatex(horiz,verti, "EEX73/KKsem305302");
  //***
  TH1D *hisRat72a = (TH1D*)Hst_Vmu72->Clone("hisRat72a");
  hisRat72a->Divide(hst_vvMu305302);
  hisRat72a->SetLineColor(kBlue);
  hisRat72a->DrawCopy("hsame");
  CaptY->SetTextColor(kBlue); verti-= 0.05;
  CaptY->DrawLatex(horiz,verti, "EEX72/KKsem305302");
  //***  also CEEX 
  TH1D *hisRat203 = (TH1D*)Hst_Vmu203->Clone("hisRat203");
  hisRat203->Divide(hst_vvMu305302);
  SetTriangle( hisRat203, kBlack, 1.0);
  hisRat203->DrawCopy("hsame");
  CaptY->SetTextColor(kBlack); verti-= 0.05;
  CaptY->DrawLatex(horiz,verti, "CEEX203/KKsem305302");
  Triangle(kBlack, 1.2, horiz-0.03, verti);
  //***  also CEEX no IFI
  TH1D *hisRat253 = (TH1D*)Hst_Vmu253->Clone("hisRat203");
  hisRat253->Divide(hst_vvMu305302);
  SetBox( hisRat253, kRed, 0.8);
  hisRat253->DrawCopy("hsame");
  CaptY->SetTextColor(kRed); verti-= 0.05;
  CaptY->DrawLatex(horiz,verti, "CEEX253/KKsem305302");
  Box(kRed, 1.2, horiz-0.03, verti);
  //
  CaptB->DrawLatex(0.50,0.01,"v=1-Q^{2}/s");
  //
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  if(ndivx>1){
  cFig3f2->cd(2);
  hst_Q2line->SetMaximum(1.04);
  hst_Q2line->SetMinimum(0.96);
  hst_Q2line->SetStats(0);
  hst_Q2line->SetTitle(0);
  hst_Q2line->SetBinContent(1,1);
  hst_Q2line->DrawCopy("h");  // horizontal line
  CaptY->SetTextColor(kBlack); verti= 0.85;
  CaptY->DrawLatex(horiz,verti, "(b) KKMC/KKsem");
  //***
  TH1D *hisRat73b = (TH1D*)Hst_Vmu73->Clone("hisRat73b");
  hisRat73b->Divide(hst_vvMu302302);
  hisRat73b->SetLineColor(kMagenta);
  hisRat73b->DrawCopy("hsame");   // KKMC/KKsem, Magenta
  CaptY->SetTextColor(kMagenta); verti-= 0.05;
  CaptY->DrawLatex(horiz,verti, "EEX73/KKsem302302");
  //***
  TH1D *hisRat72b = (TH1D*)Hst_Vmu72->Clone("hisRat72b");
  hisRat72b->Divide(hst_vvMu301301);
  hisRat72b->SetLineColor(kBlue);
  hisRat72b->DrawCopy("hsame");   // KKMC/KKsem, Blue
  CaptY->SetTextColor(kBlue); verti-= 0.05;
  CaptY->DrawLatex(horiz,verti, "EEX72/KKsem301301");
  //***
  TH1D *hisRat01b = (TH1D*)Hst_Vmu71->Clone("hisRat01b");
  hisRat01b->Divide(hst_vvMu300300);
  hisRat01b->SetLineColor(kRed);
  hisRat01b->DrawCopy("hsame");   // KKMC/KKsem, RED
  CaptY->SetTextColor(kRed); verti-= 0.05;
  CaptY->DrawLatex(horiz,verti, "EEX71/KKsem300300");
  //***
  TH1D *hisRat01g = (TH1D*)Hst_Vmu71->Clone("hisRat01g");
  hisRat01g->Divide(hst_vvMu400300);
  hisRat01g->SetLineColor(kCyan);
  hisRat01g->DrawCopy("hsame");   // KKMC/KKsem, RED
  CaptY->SetTextColor(kCyan); verti-= 0.05;
  CaptY->DrawLatex(horiz,verti, "EEX71/KKsem400300 NLL");
  //
  CaptB->DrawLatex(0.50,0.01,"v=1-Q^{2}/s");
  } //ndivx
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  if(ndivx>2){
  cFig3f2->cd(3);
  gPad->SetLogy(); // !!!!!!
  //
  Hst_Vmu74->SetStats(0);
  Hst_Vmu74->SetTitle(0);
  Hst_Vmu74->DrawCopy("h");

  Hst_Vmu73->SetLineColor(kMagenta);
  Hst_Vmu73->DrawCopy("hsame");       // Magenta

  Hst_Vmu01->SetLineColor(kBlue);     // BLUE
  Hst_Vmu01->DrawCopy("hsame");

  hst_vvMu305302->SetLineColor(kGreen);  // GREEN
  hst_vvMu305302->DrawCopy("hsame");

  hst_vvMu400300->SetLineColor(kRed);    // RED
  hst_vvMu400300->DrawCopy("hsame");
  } //ndivx
  cFig3f2->cd();

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
  Fig3m();
  //++++++++++++++++++++++++++++++++++++++++
  DiskFileA.ls();
  DiskFileB.ls();
  DiskFileB.Write();
  DiskFileB.Close();
  //++++++++++++++++++++++++++++++++++++++++
  theApp.Run();
  //++++++++++++++++++++++++++++++++++++++++
}
