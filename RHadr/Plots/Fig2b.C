
#include "HisNorm.h"
#include "Marker.h"

////////////////////////////////////////////////////////////////////////////
//   root Fig2b.C
////////////////////////////////////////////////////////////////////////////
void Fig2b(){
  gROOT.Reset();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  //TApplication theApp("App", 0, 0);
  //(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((
  //           Connect ROOT filec
  //TFile fileA("../demoC/rmain.root");
  TFile fileA("../demoC/rmain.root.180M");
  TFile fileB("./RhoSemi.root");
  //TFile fileA("../demoC/rmain.root.CEEX_50M"); // CEEX
  //TFile fileA("../demoC/rmain.root.50M"); // EEX
  //))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
  // renormalize histograms in nanobarns
  fileA->cd();
  Double_t CMSene;
  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  CMSene  *=1000;
  //
  hst10  = hst_Q2All;
  hst11  = hst_Q2Trig1;
  HisNorm(HST_KKMC_NORMA, hst10);
  HisNorm(HST_KKMC_NORMA, hst11);
  //
  hst21  = hst_Q2Trig2;
  hst22  = hst_Q2Trig2b;
  HisNorm(HST_KKMC_NORMA, hst21);
  HisNorm(HST_KKMC_NORMA, hst22);
  //
  hst30  = hst_Q2MuTrg1;
  HisNorm(HST_KKMC_NORMA, hst30);
  // Subtract muon contribution, now not necessary
  //hst11->Add(hst30,-1e0);
  //
  Float_t msize=0.6; // marker size for postscript
  BlackBullet( hst10, msize);
  RedTriangle( hst11, msize); // marker etc
  BlueBox(     hst21, msize); // marker etc
  hst10->SetLineWidth(2.5);
  hst11->SetLineWidth(2.5);
  hst10->SetMinimum(0.0);
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
  TCanvas *cFig2b = new TCanvas("cFig2b","photonic", 20, 50, WidPix,HeiPix);
  cFig2b.SetFillColor(10);
  cFig2b.cd();
  //void Divide(Int_t nx, Int_t ny, Float_t xmargin, Float_t ymargin, Int_t color)
  cFig2b->Divide(2, 1, 0.0, 0);
  ///////////////////////////////////////////////////////////////////////////////
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  cFig2b.cd(1);
  //gPad.SetLogy(); // !!!!!!
  //hst10->GetYaxis()->SetTitle("d#sigma/dQ^{2}");
  //hst10->GetYaxis()->CenterTitle();
  //hst10->GetYaxis()->SetTitleSize(0.07);
  //hst10->GetYaxis()->SetTitleOffset(0.5);
  hst11.DrawCopy("h");
  hst30.DrawCopy("hsame");   // <-- MUON
  //    CAPTIONS
  //CaptT->DrawLatex(0.50,0.98, capt1);
  CaptT->DrawLatex(0.50,0.98, "2#pi versus 2#mu");
  CaptB->DrawLatex(0.50,0.01,"Q^{2} [GeV]");
  //  description of MARKERS as a separate pad in the pad
  MarkPad = new TPad("MarkPad", "MarkPad",0.15,0.60,0.45,0.85);
  MarkPad2(MarkPad);
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  cFig2b.cd(2);
  //   PLOTTING CURVES
  //gPad.SetLogy(); // !!!!!!
  //
  hst11->Divide(hst30);
  hst11->SetMinimum(0.0);
  //
  hst11.DrawCopy("h");
  //
  fileB->cd();
  fileB->ls();
  hisKloe.DrawCopy("hsame");
  //    CAPTIONS
  CaptT->DrawLatex(0.50,0.98, "The ratio d#sigma_{had}/d#sigma_{#mu}");
  CaptB->DrawLatex(0.50,0.01,"Q^{2} [GeV]");
  //MarkPad = new TPad("MarkPad", "MarkPad",0.15,0.60,0.45,0.85);
  //MarkPad2(MarkPad);
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  cFig2b.cd();
  //------------------------------------------------------------------------
  //char ans;
  //------------------------------------------------------------------------
  //scanf("%s",&ans); // this freezes menu
  //------------------------------------------------------------------------
  //theApp->Run();
  delete theApp; // prevents form disappearing plots in canvas
};
