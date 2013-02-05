
#include "HisNorm.h"
#include "Marker.h"

////////////////////////////////////////////////////////////////////////////
//    root Fig1b.C
////////////////////////////////////////////////////////////////////////////
void Fig1b(){
  gROOT.Reset();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  //TApplication theApp("App", 0, 0);
  //(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((
  //           Connect ROOT filec
  TFile fileA("../demoC/rmain.root");
  //TFile fileA("../demoC/rmain.root.CEEX_50M"); // CEEX
  //TFile fileA("../demoC/rmain.root.50M"); //EEX
  //))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
  // renormalize histograms in nanobarns
  fileA->cd();
  Double_t CMSene;
  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  CMSene  *=1000;
  //
  hst10  = hst_PiPtTrg1;
  hst11  = hst_PiPtTrg2;
  hst20  = hst_PiCosthTrg1;
  hst21  = hst_PiCosthTrg2;
  HisNorm(HST_KKMC_NORMA, hst10);
  HisNorm(HST_KKMC_NORMA, hst11);
  HisNorm(HST_KKMC_NORMA, hst20);
  HisNorm(HST_KKMC_NORMA, hst21);
  //
  Float_t msize=0.8; // marker size for postscript
  RedTriangle( hst10, msize); // marker etc
  BlueBox(     hst11, msize); // marker etc
  RedTriangle( hst20, msize); // marker etc
  BlueBox(     hst21, msize); // marker etc
  hst10->SetLineWidth(2.5);
  hst11->SetLineWidth(2.5);
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
  ///////////////////////////////////////////////////////////////////////////////
  Float_t  WidPix, HeiPix;
  WidPix = 1000; HeiPix =  600;
  TCanvas *cFig1b = new TCanvas("cFig1b","photonic", 20, 50, WidPix,HeiPix);
  cFig1b.SetFillColor(10);
  cFig1b.cd();
  //void Divide(Int_t nx, Int_t ny, Float_t xmargin, Float_t ymargin, Int_t color)
  cFig1b->Divide(2, 1, 0.0, 0);
  ///////////////////////////////////////////////////////////////////////////////
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  cFig1b.cd(1);
  //gPad.SetLogy(); // !!!!!!
  hst10->GetYaxis()->SetTitle("d#sigma/dP^{T}_{#pi} [nb/GeV]");
  //hst10->GetYaxis()->SetTitleSize(0.05);
  hst10->GetYaxis()->SetTitleOffset(1.0);
  hst10->GetYaxis()->CenterTitle();
  hst10.SetMinimum(0.08);
  hst10.DrawCopy("h");
  hst11.DrawCopy("hsame");
  //    CAPTIONS
  CaptT->DrawLatex(0.50,0.98, capt1);
  CaptB->DrawLatex(0.50,0.01,"P^{T}_{#pi} [GeV]");
  //  description of MARKERS as a separate pad in the pad
  MarkPad = new TPad("MarkPad", "MarkPad",0.20,0.60,0.50,0.85);
  MarkPad1(MarkPad);
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  cFig1b.cd(2);
  //gPad.SetLogy(); // !!!!!!
  // POSITIONING
  //Double_t YMin=0;
  //hst20.SetMinimum(YMin);
  //Double_t YMax= 1.1*hst20.GetMaximum();
  //hst20.SetMaximum(YMax);
  //   PLOTTING CURVES
  hst20.SetMinimum(0.03);
  hst20->GetYaxis()->SetTitle("d#sigma/d cos#theta_{#pi} [nb]");
  hst20->GetYaxis()->SetTitleOffset(1.0);
  hst20->GetYaxis()->CenterTitle();
  hst20.DrawCopy("c");
  hst21.DrawCopy("hsame");
  //    CAPTIONS
  CaptT->DrawLatex(0.50,0.98, capt1);
  CaptB->DrawLatex(0.50,0.01,"cos#theta_{#pi}");
  //MarkPad = new TPad("MarkPad", "MarkPad",0.50,0.60,0.80,0.85);
  //MarkPad1(MarkPad);
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  cFig1b.cd();
  //------------------------------------------------------------------------
  //char ans;
  //------------------------------------------------------------------------
  //scanf("%s",&ans); // this freezes menu
  //------------------------------------------------------------------------
  //theApp->Run();
  delete theApp; // prevents form disappearing plots in canvas
};
