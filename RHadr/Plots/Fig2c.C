
#include "HisNorm.h"
#include "Marker.h"

////////////////////////////////////////////////////////////////////////////
//   root Fig2c.C
////////////////////////////////////////////////////////////////////////////
void Fig2c(){
  gROOT.Reset();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  //TApplication theApp("App", 0, 0);
  //(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((
  //           Connect ROOT filec
  //TFile fileA("../demoC/rmain.root");
  //TFile fileA("../demoC/rmain.root.EEX.100M"); // EEX
  TFile fileA("../demoC/rmain.root.180M"); // EEX
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
  hst30  = hst_Q2MuTrg1;
  HisNorm(HST_KKMC_NORMA, hst30);
  //==================================================================================
  cout<<"================================================================="<<endl;
  ifstream InputFile;
  InputFile.open("./43bins.dat");
  int nb=43;
  Double_t  bin[100],err[100];
  Double_t QQmin=0.10;
  Double_t QQmax=0.96;
  Double_t QQ,QQi,dsig,ddsig,sigtax;
  TH1D *H_axel = new TH1D("H_axel","Q2 distr. from Axel", nb,QQmin,QQmax);
  sigtax=0;
  for(int i=1; i<nb+1; i++){
    InputFile >>QQi >>dsig >>ddsig;
    QQ= QQmin+ (QQmax-QQmin)/nb*(i-0.5);
    bin[i-1]=dsig;
    err[i-1]=ddsig;
    cout<<i <<"  "<<QQ <<"  "<< QQi <<"  "<< dsig<<"  "<< ddsig <<endl;
    H_axel->SetBinContent(i,dsig*QQ);
    H_axel->SetBinError(i,ddsig*QQ);
    sigtax+=dsig;
  }
  InputFile.close();
  sigtax *= (QQmax-QQmin)/nb;
  cout<<"%%%%% sigma tot [nb ] from Axel= " << sigtax <<endl;
  Double_t sigtot=0;
  for(int i=1; i<nb+1; i++){
    QQ= QQmin+ (QQmax-QQmin)/nb*(i-0.5);
    sigtot+= hst21->GetBinContent(i)/QQ;
  }
  sigtot *= (QQmax-QQmin)/nb;
  cout<<"%%%%% sigma tot [nb ] from KKMC= " << sigtot <<endl;
  //==================================================================================
  Float_t msize=0.6; // marker size for postscript
  BlackBullet( hst10, msize);
  //BlueBox(     hst11, msize); // marker etc
  RedTriangle( hst11, msize); // marker etc
  //BlackBullet( hst20, msize);
  BlueBox(     hst21, msize); // marker etc
  //RedTriangle( hst21, msize); // marker etc
  hst10->SetLineWidth(2.5);
  hst11->SetLineWidth(2.5);
  //hst21->SetLineWidth(3.5);
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
  TCanvas *cFig2 = new TCanvas("cFig2","photonic", 20, 50, WidPix,HeiPix);
  cFig2.SetFillColor(10);
  cFig2.cd();
  //void Divide(Int_t nx, Int_t ny, Float_t xmargin, Float_t ymargin, Int_t color)
  cFig2->Divide(2, 1, 0.0, 0);
  ///////////////////////////////////////////////////////////////////////////////
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  cFig2.cd(1);
  hst21->SetMinimum(0.0);
  hst21.DrawCopy("h");
  //hst30.DrawCopy("hsame");   // <-- MUON
  H_axel.DrawCopy("hsame");
  //    CAPTIONS
  CaptT->DrawLatex(0.50,0.98, capt1);
  CaptB->DrawLatex(0.50,0.01,"Q^{2} [GeV]");
  //  description of MARKERS as a separate pad in the pad
  MarkPad = new TPad("MarkPad", "MarkPad",0.15,0.60,0.45,0.85);
  MarkPad2(MarkPad);
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  cFig2.cd(2);
  //   PLOTTING CURVES
  H_axel->Divide(hst21);    // <-----  TAKE RATIO
  H_axel->DrawCopy("h");
  //    CAPTIONS
  sprintf(capt1,"Axel/KKMC");
  CaptT->DrawLatex(0.50,0.98, capt1);
  CaptB->DrawLatex(0.50,0.01,"Q^{2} [GeV]");
  //MarkPad = new TPad("MarkPad", "MarkPad",0.15,0.60,0.45,0.85);
  //MarkPad2(MarkPad);
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  cFig2.cd();
  //------------------------------------------------------------------------
  theApp->Run();
  delete theApp; // prevents form disappearing plots in canvas
};
