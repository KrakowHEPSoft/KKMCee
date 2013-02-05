
#include "HisNorm.h"
#include "Marker.h"

////////////////////////////////////////////////////////////////////////////
//     TWO METHODS OF RE-NORMALIZing histograms
//     root Fig0a.C
////////////////////////////////////////////////////////////////////////////
void Fig0a(){
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
  cout<<"CMSene= "<<CMSene<<endl;
  //
  hst0=hst_Qaux0; // Normalized at the end of MC run
  hst1=hst_Qaux1; // Normalize here, locally
  HisNorm(HST_KKMC_NORMA, hst0);
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
  TCanvas *cFig0a = new TCanvas("cFig0a","Q2 distribution",1000,600);
  cFig0a.SetFillColor(10);
  cFig0a.cd();
  //void Divide(Int_t nx, Int_t ny, Float_t xmargin, Float_t ymargin, Int_t color)
  cFig0a->Divide(2, 1, 0.0, 0);
  ///////////////////////////////////////////////////////////////////////////////
  //   PLOTTING CURVES
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  cFig0a.cd(1);
  hst0.SetMinimum(0);
  hst0.DrawCopy("p");
  CaptT->DrawLatex(0.50,0.98, " psi 3.096 resonance.");
  CaptB->DrawLatex(0.50,0.01,"Q [GeV]");
  //  description of MARKERS as a separate pad in the pad
  //MarkPad = new TPad("MarkPad", "MarkPad",0.20,0.50,0.50,0.75);
  //MarkPad1(MarkPad);
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  cFig0a.cd(2);
  hst1.SetMinimum(0);
  hst1.DrawCopy("c");
  //    CAPTIONS
  CaptT->DrawLatex(0.50,0.98, " psi 3.685 resonance ");
  CaptB->DrawLatex(0.50,0.01,"Q [GeV]");
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  cFig0a.cd();
  //------------------------------------------------------------------------
  //char ans;
  //------------------------------------------------------------------------
  //scanf("%s",&ans); // this freezes menu
  //------------------------------------------------------------------------
  //theApp->Run();
  delete theApp; // prevents form disappearing plots in canvas
};
