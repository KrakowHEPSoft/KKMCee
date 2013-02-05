
#include "HisNorm.h"
#include "Marker.h"

////////////////////////////////////////////////////////////////////////////
//   root Fig3.C
////////////////////////////////////////////////////////////////////////////
void Fig3(){
  gROOT.Reset();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  //TApplication theApp("App", 0, 0);
  //(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((
  //           Connect ROOT filec
  TFile fileA("../demoC/rmain.root");
  //TFile fileA("../demoC/rmain.root.CEEX_50M"); // CEEX
  //TFile fileA("../demoC/rmain.root.50M"); // EEX
  //))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
  // renormalize histograms in nanobarns
  fileA->cd();
  Double_t CMSene;
  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  CMSene  *=1000;
  //
  hst00  = hst_Q2Trig2;
  HisNorm(HST_KKMC_NORMA, hst00);
  hst10  = hst_Q2dif1;  // EEX3  -CEEX2, O(alf^2 L), O(alf^3 L^3)
  hst11  = hst_Q2dif3;  // CEEX1 -CEEX2, O(alf^2 L^2), O(alf^3 L^3)
  hst20  = hst_Q2dif2;  // EEX2  -CEEX3, O(alf^3 L^3)
  HisNorm(HST_KKMC_NORMA, hst10);
  HisNorm(HST_KKMC_NORMA, hst11);
  HisNorm(HST_KKMC_NORMA, hst20);
  hst10->Divide(hst00);
  hst11->Divide(hst00);
  hst20->Divide(hst00);
  //
  Float_t msize=0.8; // marker size for postscript
  BlackBullet( hst10, msize); // EEX3  -CEEX2, O(alf^2 L),   O(alf^3 L^3)
  RedTriangle( hst11, msize); // CEEX1 -CEEX2, O(alf^2 L^2), O(alf^2 L)
  BlueBox(     hst20, msize); // EEX2  -CEEX3, O(alf^3 L^3)
  hst10->SetLineWidth(1.5);
  hst11->SetLineWidth(1.5);
  hst20->SetLineWidth(1.5);
  // create zero histogram, one possibility is cloning
  TH1D *hzer = (TH1D*)hst00->Clone();
  hzer->SetName("hzer"); //recommended, otherwise you get 2 histograms 
  hzer->Scale(0.0);
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
  TCanvas *cFig3 = new TCanvas("cFig3","photonic", 20, 50, WidPix,HeiPix);
  cFig3.SetFillColor(10);
  cFig3.cd();
  //void Divide(Int_t nx, Int_t ny, Float_t xmargin, Float_t ymargin, Int_t color)
  cFig3->Divide(2, 1, 0.0, 0);
  ///////////////////////////////////////////////////////////////////////////////
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  cFig3.cd(1);
  //
  hst10.SetMaximum( 0.015);
  hst10.SetMinimum(-0.015);
  hst10.DrawCopy("h");
  hst11.DrawCopy("hsame");
  hzer->DrawCopy("hsame");
  //    CAPTIONS
  CaptT->DrawLatex(0.50,0.98, capt1);
  CaptB->DrawLatex(0.50,0.01,"Q^{2} [GeV]");
  //  description of MARKERS as a separate pad in the pad
  MarkPad = new TPad("MarkPad", "MarkPad",0.55,0.60,0.85,0.85);
  MarkPad3(MarkPad);
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  cFig3.cd(2);
  //
  hst20.SetMaximum( 0.0010);
  hst20.SetMinimum(-0.0005);
  hst20.DrawCopy("h");
  hzer->DrawCopy("hsame");
  //    CAPTIONS
  CaptT->DrawLatex(0.50,0.98, capt1);
  CaptB->DrawLatex(0.50,0.01,"Q^{2} [GeV]");
  MarkPad = new TPad("MarkPad", "MarkPad",0.55,0.60,0.85,0.85);
  MarkPad3(MarkPad);
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  cFig3.cd();
  //------------------------------------------------------------------------
  //char ans;
  //------------------------------------------------------------------------
  //scanf("%s",&ans); // this freezes menu
  //------------------------------------------------------------------------
  //theApp->Run();
  delete theApp; // prevents form disappearing plots in canvas
};
