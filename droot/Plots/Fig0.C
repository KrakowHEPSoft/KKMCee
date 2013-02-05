
#include "HisNorm.h"
#include "Marker.h"

////////////////////////////////////////////////////////////////////////////
//     TWO METHODS OF RE-NORMALIZing histograms
//     root Fig0.C
////////////////////////////////////////////////////////////////////////////
void Fig0(){
  gROOT.Reset();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  //TApplication theApp("App", 0, 0);
  //(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((
  //           Connect ROOT filec
  TFile fileA("../demoC/rmain.root");
  //))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
  // renormalize histograms in nanobarns
  fileA->cd();
  Double_t CMSene;
  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  CMSene  *=1000;
  //
  hst0=HST_Q2kloe; // Normalized at the end of MC run
  hst1=hst_Q2kloe; // Normalize here, locally
  HisNorm(HST_KKMC_NORMA, hst1);
  //
  Float_t msize=0.4; // marker size for postscript
  BlueBox(     hst0, msize); // marker etc
  BlackBullet( hst1, msize);
  //==================================================================================
  // top caption
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextAlign(23);
  CaptT->SetTextSize(0.050);
  // bottom caption
  TLatex *CaptB = new TLatex();
  CaptB->SetNDC(); // !!!
  CaptB->SetTextAlign(21);
  CaptB->SetTextSize(0.045);
  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFig0 = new TCanvas("cFig0","Q2 distribution",1000,600);
  cFig0.SetFillColor(10);
  cFig0.cd();
  //void Divide(Int_t nx, Int_t ny, Float_t xmargin, Float_t ymargin, Int_t color)
  cFig0->Divide(2, 1, 0.0, 0);
  ///////////////////////////////////////////////////////////////////////////////
  // POSITIONING
  Double_t YMin=0;
  hst0.SetMinimum(YMin);
  hst1.SetMinimum(YMin);
  Double_t YMax= 1.1*hst0.GetMaximum();
  hst0.SetMaximum(YMax);
  hst1.SetMaximum(YMax);
  //Double_t XMax= hst1.GetMaximum();
  //Double_t Xcenter = 0.5*(hst1->GetXaxis()->GetXmax())+0.5*(hst1->GetXaxis()->GetXmin());
  //   PLOTTING CURVES
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  cFig0.cd(1);
  hst0.DrawCopy("p");
  CaptT->DrawLatex(0.50,0.98, "normalized in MC main prog.");
  CaptB->DrawLatex(0.50,0.01,"Q^{2}");
  //  description of MARKERS as a separate pad in the pad
  MarkPad = new TPad("MarkPad", "MarkPad",0.20,0.50,0.50,0.75);
  MarkPad1(MarkPad);
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  cFig0.cd(2);
  hst1.DrawCopy("c");
  //    CAPTIONS
  CaptT->DrawLatex(0.50,0.98, "normalized here, localy");
  CaptB->DrawLatex(0.50,0.01,"Q^{2}");
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  cFig0.cd();
  //------------------------------------------------------------------------
  //char ans;
  //------------------------------------------------------------------------
  //scanf("%s",&ans); // this freezes menu
  //------------------------------------------------------------------------
  //theApp->Run();
  delete theApp; // prevents form disappearing plots in canvas
};
