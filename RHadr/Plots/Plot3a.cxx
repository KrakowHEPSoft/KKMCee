//////////////////////////////////////////////////////////////////////
//    make Plot3a
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
TFile DiskFileA("../demoC/rmain.root.EEX.960M.Ord1"); //
//TFile DiskFileA("../demoC/rmain.root.EEX.500M");
//TFile DiskFileA("../demoC/rmain.root.180M");
TFile DiskFileB("RhoSemi.root","RECREATE","histograms");
//=============================================================================

// Auxiliary procedures for plotting
#include "HisNorm.h"
#include "Marker.h"

Double_t sqr( const Double_t x ){ return x*x;};


///////////////////////////////////////////////////////////////////////////////////
void HistNormalize(){
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  //
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
int Fig3a()
{
//------------------------------------------------------------------------
  cout<<" ========================= Fig3a =========================== "<<endl;

  TH1D *H_Phk1A = new TH1D("H_Phk1A","Q2 distr. Born, incl.", 100, 0.00, 1.0);
  TH1D *H_Phk2A = new TH1D("H_Phk2A","Q2 distr. Born, incl.", 100, 0.00, 1.0);
  ReadHist(H_Phk1A, "./phokara_born_1_qq.dat");  // no cut, high stat
  ReadHist(H_Phk2A, "./phokara_nlo_1_qq.dat");   // no cut, high stat
  TH1D *H_Phk1B = new TH1D("H_Phk1B","Q2 distr. Born, incl.", 100, 0.00, 1.0);
  TH1D *H_Phk2B = new TH1D("H_Phk2B","Q2 distr. Born, incl.", 100, 0.00, 1.0);
  ReadHist(H_Phk1B, "./phokara_born_2_qq.dat"); // no cut
  ReadHist(H_Phk2B, "./phokara_nlo_2_qq.dat");  // no cut
  //
  Double_t CMSene;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");

  CMSene  = HST_KKMC_NORMA->GetBinContent(1)/HST_KKMC_NORMA->GetBinContent(511);
  CMSene  *=1000;
  //
  Float_t msize=0.6; // marker size for postscript
  TH1D *Hst_Q2piA = (TH1D*)DiskFileA.Get("Hst_Q2piA"); // no cut
  TH1D *Hst_Q2piB = (TH1D*)DiskFileA.Get("Hst_Q2piB"); // with cut
  TH1D *Hst_Q2muA = (TH1D*)DiskFileA.Get("Hst_Q2muA"); // no cut
  Hst_Q2piB->Scale(0.5);
  //
  TH1D *Hst_piCosA = (TH1D*)DiskFileA.Get("Hst_piCosA");
  TH1D *Hst_piCosB = (TH1D*)DiskFileA.Get("Hst_piCosB");
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
  //***************************
  TLatex *CaptY = new TLatex();
  CaptY->SetNDC(); CaptY->SetTextAlign(12);
  CaptY->SetTextSize(0.040);
  Float_t verti=0.85;
  Float_t horiz=0.50;
  TH1D *hst_Q2line  = new TH1D("hst_Q2line","one",  1, 0.0,1.0);
  hst_Q2line->SetBinContent(1,1);
  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  Float_t  WidPix, HeiPix;
  WidPix = 1100; HeiPix =  550;
  TCanvas *cFig3a1 = new TCanvas("cFig3a1","Fig3a1 photonic", 40,  0, WidPix,HeiPix);
  gPad->SetFillColor(10);
  //void Divide(Int_t nx, Int_t ny, Float_t xmargin, Float_t ymargin, Int_t color)
  cFig3a1->Divide(2, 0, 0.0, 0);
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  cFig3a1->cd(1);
  gPad->SetFillColor(10);
  gPad->SetLogy(); // !!!!!!
  Hst_Q2piA->SetStats(0);
  Hst_Q2piA->SetTitle(0);
  Hst_Q2piA->DrawCopy("h");   // <-- 2pi
  //***
  Hst_Q2muA->SetLineColor(kGreen);  // green
  Hst_Q2muA->DrawCopy("hsame"); // <-- MUON
  //*** Phokara
  //H_Phk1A->SetLineColor(kMagenta);  // Born magenta
  //H_Phk1A->DrawCopy("hsame");
  //***
  H_Phk2A->SetLineColor(kBlue);  // NLO  red
  H_Phk2A->DrawCopy("hsame");
  CaptY->DrawLatex(horiz,verti, "(a)");
  //***
  CaptT->DrawLatex(0.50,0.98, "2#pi KKMC&PHOKHARA, NO CUTS");
  CaptB->DrawLatex(0.50,0.01,"Q^{2} [GeV^{2}]");
  //**************************************************************
  cFig3a1->cd(2);
  gPad->SetFillColor(10);
  hst_Q2line->SetMaximum(1.05);
  hst_Q2line->SetMinimum(0.95);
  hst_Q2line->SetStats(0);
  hst_Q2line->SetTitle(0);
  hst_Q2line->SetBinContent(1,1);
  hst_Q2line->DrawCopy("h");  // horizontal line
  CaptY->SetTextColor(kBlack); verti= 0.85;
  CaptY->DrawLatex(horiz-0.05,verti, "(b) KKMC/PHOKHARA");
  /*
  TH1D *hisRat1 = (TH1D*)Hst_Q2piA->Clone("hisRat1");
  hisRat1->Divide(H_Phk1A);     // Born
  SetTriangle( hisRat1, kMagenta, 1.0);
  hisRat1->DrawCopy("hsame");       // Born
  CaptY->SetTextColor(kMagenta); verti-= 0.05;
  CaptY->DrawLatex(horiz,verti, "EEX3/PHOKHARA1");
  Triangle(kMagenta, 1.2, horiz-0.03, verti);
  */
  TH1D *hisRat2 = (TH1D*)Hst_Q2piA->Clone("hisRat2");
  hisRat2->Divide(H_Phk2A);
  SetBullet( hisRat2, kBlue, 0.9);
  hisRat2->DrawCopy("hsame");   // NLO
  CaptY->SetTextColor(kBlue); verti-= 0.05;
  CaptY->DrawLatex(horiz,verti, "EEX3/PHOKHARA2");
  Bullet(kBlue, 1.1, horiz-0.03, verti);
  //***
  CaptT->DrawLatex(0.50,0.98, "NO CUTS");
  CaptB->DrawLatex(0.50,0.01,"Q^{2} [GeV^{2}]");
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  cFig3a1->cd();
  //
  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  //Float_t  WidPix, HeiPix;
  WidPix = 1100; HeiPix =  550;
  TCanvas *cFig3a2 = new TCanvas("cFig3a2","Fig3a2 photonic", 80, 40, WidPix,HeiPix);
  gPad->SetFillColor(10);
  cFig3a2->Divide(2, 0, 0.0, 0);
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  cFig3a2->cd(1);
  gPad->SetFillColor(10);
  gPad->SetLogy(); // !!!!!!
  Hst_Q2piB->SetStats(0);
  Hst_Q2piB->SetTitle(0);
  Hst_Q2piB->DrawCopy("h");   // <-- 2pi
  Hst_Q2piB->DrawCopy("esame");   // <-- 2pi
  // Phokara
  //H_Phk1B->SetLineColor(6);  // magenta
  //H_Phk1B->DrawCopy("hsame");
  H_Phk2B->SetLineColor(2);  // red
  H_Phk2B->DrawCopy("hsame");
  //***
  CaptY->SetTextColor(kBlack);verti=0.85;
  CaptY->DrawLatex(horiz,verti, "(c)");
  CaptT->DrawLatex(0.50,0.98, "2#pi KKMC&PHOKHARA, WITH CUTS");
  CaptB->DrawLatex(0.50,0.01,"Q^{2} [GeV^{2}]");
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  cFig3a2->cd(2);
  gPad->SetFillColor(10);
  hst_Q2line->SetMaximum(1.05);
  hst_Q2line->SetMinimum(0.95);
  hst_Q2line->SetStats(0);
  hst_Q2line->SetTitle(0);
  hst_Q2line->SetBinContent(1,1);
  hst_Q2line->DrawCopy("h");  // horizontal line
  CaptY->SetTextColor(kBlack); verti= 0.85;
  CaptY->DrawLatex(horiz-0.05,verti, "(d) KKMC/PHOKHARA");
  /*
  TH1D *hisRat3 = (TH1D*)Hst_Q2piB->Clone("hisRat3");
  hisRat3->Divide(H_Phk1B);
  SetTriangle( hisRat3, kMagenta, 1.0);
  hisRat3->DrawCopy("hsame");       //
  CaptY->SetTextColor(kMagenta); verti-= 0.05;
  CaptY->DrawLatex(horiz,verti, "EEX3/PHOKHARA1");
  Triangle(kMagenta, 1.2, horiz-0.03, verti);
  */
  TH1D *hisRat4 = (TH1D*)Hst_Q2piB->Clone("hisRat4");
  hisRat4->Divide(H_Phk2B);
  SetBullet( hisRat4, kBlue, 0.9);
  hisRat4->DrawCopy("hsame");   //
  CaptY->SetTextColor(kBlue); verti-= 0.05;
  CaptY->DrawLatex(horiz,verti, "EEX3/PHOKHARA2");
  Bullet(kBlue, 1.1, horiz-0.03, verti);
  //***
  CaptT->DrawLatex(0.50,0.98, "WITH CUTS");
  CaptB->DrawLatex(0.50,0.01,"Q^{2} [GeV^{2}]");
  //***
  cFig3a2->Update();
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  // 
  int nb = Hst_Q2piA->GetNbinsX(); int i;
  Double_t QQ,QQi,dsig,ddsig;
  cout<<nb<<endl;
  for(int i=1; i<nb+1; i++){
    QQ = Hst_Q2piA->GetBinCenter(i);
    cout<<i <<"  "<<QQ <<"  "<< Hst_Q2piA->GetBinContent(i)
                       <<"  "<< Hst_Q2piA->GetBinError(i) <<endl;
  }
  ///////////////////////////////////////////////////////////////////////////////
  return 0;
}



///////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  //++++++++++++++++++++++++++++++++++++++++
  TApplication theApp("theApp", &argc, argv);
  //++++++++++++++++++++++++++++++++++++++++
  HistNormalize();     // Renormalization of MC histograms
  Fig3a();
  //++++++++++++++++++++++++++++++++++++++++
  DiskFileA.ls();
  DiskFileB.ls();
  DiskFileB.Write();
  DiskFileB.Close();
  //++++++++++++++++++++++++++++++++++++++++
  theApp.Run();
  //++++++++++++++++++++++++++++++++++++++++
}
