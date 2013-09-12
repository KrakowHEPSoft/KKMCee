//////////////////////////////////////////////////////////////////////
//    make Plot3b
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
//TFile DiskFileA("../demoC/rmain.root.MuCEEX.21M");
TFile DiskFileA("../demoC/rmain.root.EEX.960M.Ord1"); //
//TFile DiskFileA("../demoC/rmain.root.EEX.244M");
//TFile DiskFileA("../demoC/rmain.root.EEX.500M");
//TFile DiskFileA("../demoC/rmain.root.MuEEX.32M");
//TFile DiskFileA("../demoC/rmain.root.180M");
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
  HisNorm(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("Hst_Q2hadA") );
  HisNorm(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("Hst_Q2piA") );
  HisNorm(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("Hst_Q2piB") );
  //
  HisNorm(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("Hst_Q2muA") );
  //
  HisNorm(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("Hst_piCosA") );
  HisNorm(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("Hst_piCosB") );
  HisNorm(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("Hst_phCosA") );
  HisNorm(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("Hst_phCosB") );
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
//-------------------------KKsem------------------------------------------  
// Here we produce semianalytical plots using KKsem program, No plotting
//------------------------------------------------------------------------  
  long KF,Key, KeyFob;
  char chak[6]="VRHO2";
  Key = 302;   // O(alf2)
  Key = 304;   // O(alf3) GribovLL
  Key = 303;   // O(alf3)
  Key = 305;   // O(alf3) GribovLL +NLL
  //Key = 662;  // Unexp ???? 
  //----------------------------
  Double_t CMSene;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  double svar= CMSene*CMSene;
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
  double v,Q2,Jacob;
  double vbin[Nb];
//------------------------------------------------------------------------
//  PiPi dsig/dQ^2
//------------------------------------------------------------------------
  cout<<"**** PiPi dsig/dQ^2 **** "<<endl;
  long IsQuarkGenerated;
  KF=1; bornv_getisgenerated_( KF, IsQuarkGenerated);
  long IsGenerated;
  KF=13; IsGenerated=0; bornv_setisgenerated_( KF, IsGenerated); // exclude muon !!!
  KF=1; kksem_setkffin_(KF); // set one d-quark for special tests only
  Jacob = 1/svar; // this is dv/dQ^2
  KeyFob= -100; // KKsem_BornV, DOUBLY KF-summed! WRONG!!!! CosTh-integr.(*10)
  KeyFob=   10; // BornV_Dizet(1,CosTh=0), one m_KFfin only, EW alowed   (*1/3) ????
  KeyFob=  -11; // BornV_Simple(CosTh=0),  one m_KFfin only, KeyLib=0 (R=1/3= NC*(1/3)^2) OK
  KeyFob=  -10; // KKsem_BornV(CosTh=0),   KF-sumed, (R=6/3=2*NC*(1/3)^2+NC*(2/3)^2) OK!!!!
  kksem_setkeyfob_(KeyFob);
  if(IsQuarkGenerated)  kksem_vvplot_vec_(Key,chak,Nb,vmin,vmax,vbin);
  TH1D *hst_Q2PiOrd3  = new TH1D("hst_Q2PiOrd3","PiPi, dsig/dQ^2",  Nb, qmin,qmax);
  hst_Q2PiOrd3->Sumw2();
  hst_Q2PiOrd3->SetLineColor(9); // blue
  double dv=(vmax-vmin)/Nb;
  double sigTag=0.0;
  for(int ib=0;ib<Nb;ib++){
    v=vmin+(vmax-vmin)*(ib+0.0)/Nb;
    Q2= svar*(1-v);
    cout<<"**** ib= "<<ib<<" v= "<<v<<" Q2= "<<Q2<<" vbin(ib)= "<<vbin[ib]<<endl;
    hst_Q2PiOrd3->SetBinContent(Nb-ib, vbin[ib]*Jacob); // bin content
    hst_Q2PiOrd3->SetBinError(Nb-ib,0.0);
    sigTag += vbin[ib]*dv;
  }
  cout<<"%%%%%%%%%%  sigTag [nb]= "<<sigTag<<endl;
//------------------------------------------------------------------------
//  MuMu dsig/dQ^2
//------------------------------------------------------------------------  
  cout<<"**** MuMu dsig/dQ^2 **** "<<endl;
  Jacob = 1/svar; // this is dv/dQ^2
  KF=13; bornv_setkf_( KF ); // set SINGLE Final State
  KF=13; kksem_setkffin_(KF); // set m_KFfin in KKsem
  KeyFob=   10; // BornV_Dizet, not good
  KeyFob=  -11; // BornV_Simple for KeyLib=0, OK
  KeyFob= -100; // KKsem_BornV, without EW, with integration, OK
  KeyFob=  -10; // KKsem_BornV, OK!
  kksem_setkeyfob_(KeyFob);
  //
  kksem_vvplot_vec_(Key,chak,Nb,vmin,vmax,vbin);
  TH1D *hst_Q2MuOrd3  = new TH1D("hst_Q2MuOrd3","MuMu  dsig/dQ^2",  Nb, qmin,qmax);
  hst_Q2MuOrd3->Sumw2();
  hst_Q2MuOrd3->SetLineColor(8); // green
  for(int ib=0;ib<Nb;ib++){
    v=vmin+(vmax-vmin)*(ib+0.5)/Nb;
    Q2= svar*(1-v);
    cout<<"**** ib= "<<ib<<" v= "<<v<<" vbin(ib)= "<<vbin[ib]<<endl;
    hst_Q2MuOrd3->SetBinContent(Nb-ib, vbin[ib]*Jacob); // bin content
    hst_Q2MuOrd3->SetBinError(Nb-ib,0.0);
  }
  // UNEXP
  Key = 662;  // Unexp old stuff???? 
  Key = 632;  // Unexp O(alf2)
  //Key = 635;  // Unexp O(alf2)
  //Key = 631;  // Unexp O(alf1)
  kksem_vvplot_vec_(Key,chak,Nb,vmin,vmax,vbin);
  TH1D *hst_Q2MuO2  = new TH1D("hst_Q2MuO2","MuMu  dsig/dQ^2",  Nb, qmin,qmax);
  hst_Q2MuO2->Sumw2();
  hst_Q2MuO2->SetLineColor(4); // blue
  for(int ib=0;ib<Nb;ib++){
    v=vmin+(vmax-vmin)*(ib+0.5)/Nb;
    Q2= svar*(1-v);
    cout<<"**** ib= "<<ib<<" v= "<<v<<" vbin(ib)= "<<vbin[ib]<<endl;
    hst_Q2MuO2->SetBinContent(Nb-ib, vbin[ib]*Jacob); // bin content
    hst_Q2MuO2->SetBinError(Nb-ib,0.0);
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
int Fig3b()
{
//------------------------------------------------------------------------
  cout<<" ========================= Fig3b =========================== "<<endl;

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

  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  CMSene  *=1000;
  //
  Float_t msize=0.6; // marker size for postscript
  TH1D *Hst_Q2hadA = (TH1D*)DiskFileA.Get("Hst_Q2hadA");// no cut
  TH1D *Hst_Q2piA  = (TH1D*)DiskFileA.Get("Hst_Q2piA"); // no cut
  TH1D *Hst_Q2muA  = (TH1D*)DiskFileA.Get("Hst_Q2muA"); // no cut
  //
  BlackBullet( Hst_Q2piA, msize);
  Hst_Q2piA->SetLineWidth(2.5);
  //
  TH1D *hst_Q2MuOrd3 = (TH1D*)DiskFileB.Get("hst_Q2MuOrd3");
  TH1D *hst_Q2MuO2   = (TH1D*)DiskFileB.Get("hst_Q2MuO2");
  TH1D *hst_Q2PiOrd3 = (TH1D*)DiskFileB.Get("hst_Q2PiOrd3");
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
  // inside caption
  TLatex *CaptC = new TLatex();
  CaptC->SetNDC(); // !!!
  CaptC->SetTextAlign(11);
  CaptC->SetTextSize(0.040);
  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  Float_t  WidPix, HeiPix;
  WidPix = 1100; HeiPix =  700;
  TCanvas *cFig3b1 = new TCanvas("cFig3b1","Fig3b1 photonic", 40,  0, WidPix,HeiPix);
  cFig3b1->SetFillColor(10);
  cFig3b1->cd();
  //void Divide(Int_t nx, Int_t ny, Float_t xmargin, Float_t ymargin, Int_t color)
  cFig3b1->Divide(2, 0, 0.0, 0);
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  cFig3b1->cd(1);
  gPad->SetLogy(); // !!!!!!
  Hst_Q2muA->SetStats(0);
  Hst_Q2muA->SetTitle(0);
  //
  Hst_Q2muA->SetLineColor(8);      // GREEN
  Hst_Q2muA->SetMinimum(10.0);
  Hst_Q2muA->DrawCopy("h");        // <-- MUON KKMC.EEX
  //
  hst_Q2MuO2->SetLineColor(4);     // BLUE
  hst_Q2MuO2->DrawCopy("hsame");   // <-- MUON KKSEM.Ord2
  //
  H_Phk2Amu->SetLineColor(2);      // RED
  H_Phk2Amu->DrawCopy("hsame");    // phokhara red
  //
  hst_Q2MuOrd3->SetLineColor(1);   // BLACK
  hst_Q2MuOrd3->DrawCopy("hsame"); // <-- MUON KKSEM Exp3
  //
  CaptT->DrawLatex(0.50,0.98, "KKMC muons, NO CUTS");
  CaptB->DrawLatex(0.50,0.01,"Q^{2} [GeV]");
  CaptC->SetTextColor(1);      // BLACK
  CaptC->DrawLatex(0.20,0.85,"KKsem.O3exp O(#alpha^{3})_{exp} ");
  CaptC->SetTextColor(4);      // BLUE
  CaptC->DrawLatex(0.20,0.80,"KKsem.O2 O(#alpha^{2}) ");
  CaptC->SetTextColor(2);      // RED
  CaptC->DrawLatex(0.20,0.75,"PHOKARA.O2 O(#alpha^{2}) ");
  CaptC->SetTextColor(8);      // GREEN
  CaptC->DrawLatex(0.20,0.70,"KKMC.EEX2 O(#alpha^{2})_{exp} ");
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  cFig3b1->cd(2);
  //TH1D *hisRat1 = (TH1D*)hst_Q2MuOrd3->Clone("hisRat1");
  //hisRat1->Divide(Hst_Q2muA);     // KKsem
  TH1D *hisRat1 = (TH1D*)Hst_Q2muA->Clone("hisRat1");
  hisRat1->Divide(hst_Q2MuOrd3);     // KKsem
  hisRat1->SetMaximum(1.05);
  hisRat1->SetMinimum(0.95);
  //
  hisRat1->SetStats(0);
  hisRat1->SetTitle(0);
  hisRat1->SetLineColor(8);     //
  hisRat1->DrawCopy("h");       // KKMC.EEX/KKsemExp3  GREEN
  //
  TH1D *hisRat1mu = (TH1D*)H_Phk2Amu->Clone("hisRat1mu");
  hisRat1mu->Divide(hst_Q2MuOrd3);
  hisRat1mu->SetLineColor(2);   // phokhara/KKsemExp3  RED
  hisRat1mu->DrawCopy("hsame");
  //
  TH1D *hisRat3 = (TH1D*)Hst_Q2muA->Clone("hisRat3");
  hisRat3->Divide(H_Phk2Amu);
  hisRat3->SetLineColor(6);     // KKMC.EEX/Phokhara2  Magenta
  hisRat3->DrawCopy("hsame");
   //
  TH1D *hisRat4 = (TH1D*)hst_Q2MuO2->Clone("hisRat4");
  hisRat4->Divide(hst_Q2MuOrd3);
  hisRat4->SetLineColor(4);     // KKsemO2/KKsemExp3  BLUE
  hisRat4->DrawCopy("hsame");
  //
  CaptT->DrawLatex(0.50,0.98, "Ratios");
  CaptB->DrawLatex(0.50,0.01,"Q^{2} [GeV]");
  CaptC->SetTextColor(8);      // GREEN
  CaptC->DrawLatex(0.15,0.85,"KKMC.EEX2/KKsem.O3exp ");
  CaptC->SetTextColor(6);      // Magenta
  CaptC->DrawLatex(0.15,0.80,"KKMC.EEX2/PHOKARA.O2 ");
  CaptC->SetTextColor(4);      // BLUE
  CaptC->DrawLatex(0.15,0.75,"KKsem.O2/KKsem.O3exp ");
  CaptC->SetTextColor(2);      // RED
  CaptC->DrawLatex(0.15,0.70,"PHOKARA.O2/KKsem.O3exp ");
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  cFig3b1->cd();
  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  //Float_t  WidPix, HeiPix;
  WidPix = 1100; HeiPix =  700;
  TCanvas *cFig3b2 = new TCanvas("cFig3b2","Fig3b2 photonic", 80, 40, WidPix,HeiPix);
  cFig3b2->SetFillColor(10);
  cFig3b2->cd();
  cFig3b2->Divide(2, 0, 0.0, 0);
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  cFig3b2->cd(1);
  gPad->SetLogy(); // !!!!!!
  Hst_Q2piA->SetStats(0);
  Hst_Q2piA->SetTitle(0);
  Hst_Q2piA->DrawCopy("h");       // <-- 2pi
  Hst_Q2piA->DrawCopy("esame");   // <-- 2pi
  //
  Hst_Q2hadA->SetLineColor(6);    // magenta
  Hst_Q2hadA->DrawCopy("esame");  // <-- 2pi+3pi
  //
  hst_Q2PiOrd3->DrawCopy("hsame");

  CaptT->DrawLatex(0.50,0.98, "KKMC&KKsem, NO CUTS");
  CaptB->DrawLatex(0.50,0.01,"Q^{2} [GeV]");
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  cFig3b2->cd(2);
  TH1D *hisRat2 = (TH1D*)Hst_Q2hadA->Clone("hisRat2");
  hisRat2->Divide(hst_Q2PiOrd3);
  hisRat2->SetMaximum(1.05);
  hisRat2->SetMinimum(0.95);
  //
  hisRat2->SetStats(0);
  hisRat2->SetTitle(0);
  hisRat2->SetLineColor(9);     // blue KKsem/KKMC
  hisRat2->DrawCopy("h");       // ratio

  CaptT->DrawLatex(0.50,0.98, "KKMC.EEX/KKsem");
  CaptB->DrawLatex(0.50,0.01,"Q^{2} [GeV]");
  ///////////////////////////////////////////////////////////////////////////////
  cFig3b2->cd();
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
  Fig3b();
  //++++++++++++++++++++++++++++++++++++++++
  DiskFileA.ls();
  DiskFileB.ls();
  DiskFileB.Write();
  DiskFileB.Close();
  //++++++++++++++++++++++++++++++++++++++++
  theApp.Run();
  //++++++++++++++++++++++++++++++++++++++++
}
