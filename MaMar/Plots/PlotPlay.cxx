
//    (make PlotPlay; ./PlotPlay)

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

TFile DiskFileB("PlotPlay.root","RECREATE","histograms");
///////////////////////////////////////////////////////////////////////////////////
//              GLOBAL stuff
///////////////////////////////////////////////////////////////////////////////////
int    kGold=kOrange-3, kBrune=46, kPine=kGreen+3;
//
float  gXcanv = 0, gYcanv = 0;
//

///////////////////////////////////////////////////////////////////////////////////
void FigPrag2()
{
//------------------------------------------------------------------------
  int nbx=5, nby = 5;
  cout<<" ========================= FigPrag2 =========================== "<<endl;
  TH2D *ISRee = new TH2D("ISRee",    " ISRee ", nbx, 0.0 ,5.0, nby, 0.0 ,5.0);
  TH2D *ISRmu = new TH2D("ISRmu",    " ISRmu ", nbx, 0.0 ,5.0, nby, 0.0 ,5.0);

  double alfpi = 1.0/400.0;
  double MZ = 91;
  double m_e  = 0.000501;
  double m_mu = 0.105;
  double Le  = 2*log(MZ/m_e);
  double Lmu = 2*log(MZ/m_mu);
  cout<<"FigPrag2:  Le="<< Le<<endl;
  double A,B;
//
  for(int n=0; n <= 5; n++)
	  for(int j=0; j <= 5; j++)
	  {
		  if(j <= n){
    	  	  A = exp(n*log(alfpi)) * exp(j*log(Le));
    	  	  B = exp(n*log(alfpi)) * exp(j*log(Lmu));
		  }else{
			  A= 1.0e-100;
			  B= 1.0e-100;
		  }
		  ISRee->Fill(0.5+n,0.5+j, A);
		  ISRmu->Fill(0.5+n,0.5+j, B);
		  cout<<"n,j ="<< n <<"  " << j << "  A ="<< A << "  B ="<< B <<endl;
	  }

  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex(); CaptT->SetNDC(); // !!!
  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cPrag2 = new TCanvas("cPrag2","cPrag2", gXcanv,  gYcanv,   1600,  800);
  //                            Name    Title            xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  gXcanv += 50; gYcanv += 50;
  cPrag2->SetFillColor(10);
  cPrag2->Divide( 2,  0);
  //cPrag2->Divide( 2,  0,     0.0,     0.0,   10);
  //              nx, ny, xmargin, ymargin, color
  //////////////////////////////////////////////
  cPrag2->cd(1);
  CaptT->DrawLatex(0.12,0.95,"QED LO, NLO ...");

  //////////////////////////////////////////////
  TString OptSurf;
  //OptSurf="      "; // 2D scatergram, points
  //OptSurf="col   "; // 2D histogram, color
  //OptSurf="colz  "; // 2D kolorowe paski, ze skala
  //OptSurf="surf1 "; // 3D surface color
  OptSurf="lego2 "; // 3D histogram color
  //OptSurf="surf3 "; // 3D histogram, z plotem "na dachu"
  //OptSurf="surf2z"; // 3D kolorowe paski, ze skala
  //OptSurf="surf2 "; // 3D kolorowe paski bez skali
  //OptSurf="surf4 "; // 3D gladka powierchnia
  //-------------------------------------
  //
  gPad->SetLogz(); // !!!!!!
  gPad->SetTheta(11);
  gPad->SetPhi(-117);
  //
  ISRee->SetStats(0);
  ISRee->SetTitle(0);
  double zmax = 1.0;
  //zmax= ISRee->GetMaximum();
  ISRee->SetMaximum(zmax*1.1);
  ISRee->SetMinimum(zmax*1e-7);
  //
  ISRee->GetXaxis()->SetTitle("n");
  ISRee->GetXaxis()->CenterTitle();
  //ISRee->GetXaxis()->SetTitleOffset(1.5);
  ISRee->GetXaxis()->SetTitleSize(0.045);
  ISRee->GetXaxis()->SetNdivisions(5);
  //
  ISRee->GetYaxis()->SetTitle("r");
  ISRee->GetYaxis()->CenterTitle();
  //ISRee->GetYaxis()->SetTitleOffset(1.5);
  ISRee->GetYaxis()->SetTitleSize(0.045);
  ISRee->GetYaxis()->SetNdivisions(5);
  //
  ISRee->DrawCopy(OptSurf);
  //
  CaptT->DrawLatex(0.10,0.96,"QED strength, ISR e^{+}e^{-}");

  ///////////////////////////////////////////////////////
  cPrag2->cd(2);
  gPad->SetLogz(); // !!!!!!
  gPad->SetTheta(11);
  gPad->SetPhi(-117);
  //
  ISRmu->SetStats(0);
  ISRmu->SetTitle(0);
  //zmax= ISRmu->GetMaximum();
  ISRmu->SetMaximum(zmax*1.1);
  ISRmu->SetMinimum(zmax*1e-7);
  //
  ISRmu->GetXaxis()->SetTitle("n");
  ISRmu->GetXaxis()->CenterTitle();
  //ISRmu->GetXaxis()->SetTitleOffset(1.5);
  ISRmu->GetXaxis()->SetTitleSize(0.045);
  ISRmu->GetXaxis()->SetNdivisions(5);
  //
  ISRmu->GetYaxis()->SetTitle("r");
  ISRmu->GetYaxis()->CenterTitle();
  //ISRmu->GetYaxis()->SetTitleOffset(1.5);
  ISRmu->GetYaxis()->SetTitleSize(0.045);
  ISRmu->GetYaxis()->SetNdivisions(5);
  //
  ISRmu->DrawCopy(OptSurf);
  //
  CaptT->DrawLatex(0.10,0.96,"QED strength, FSR #mu^{+}#mu^{-}");
  cPrag2->cd();
  //================================================
  cPrag2->SaveAs("cPrag2.pdf");
  cPrag2->SaveAs("cPrag2.jpg");
  cPrag2->SaveAs("cPrag2.png");
//
}// FigTempl


///////////////////////////////////////////////////////////////////////////////////
void FigFCCee1()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigFCCee1 =========================== "<<endl;
  TH1D *h_lep = new TH1D("h_lep" ,  "LEP exp error.",      10, 0.0, 10.0);
  TH1D *h_qed = new TH1D("h_qed" ,  "LEP EQD uncert",      10, 0.0, 10.0);
  TH1D *h_fcc = new TH1D("h_fcc" ,  "FCCee1 exp. err",      10, 0.0, 10.0);

  h_lep->SetBinContent( 1,  2.2);
  h_lep->SetBinContent( 2,  2.2);
  h_lep->SetBinContent( 3, 25.0);
  h_lep->SetBinContent( 4, 37.0);
  h_lep->SetBinContent( 5,  8.0);
  h_lep->SetBinContent( 6,150.0);
  h_lep->SetBinContent( 7, 53.0);
  h_lep->SetBinContent( 8, 41.0);
  h_lep->SetBinContent( 9, 33.0);
  h_lep->SetBinContent(10, 2000);

  h_qed->SetBinContent( 1,  0.3);
  h_qed->SetBinContent( 2,  0.2);
  h_qed->SetBinContent( 3, 12.0);
  h_qed->SetBinContent( 4, 25.0);
  h_qed->SetBinContent( 5,  6.0);
  h_qed->SetBinContent( 6, 60.0);
  h_qed->SetBinContent( 7, 28.0);
  h_qed->SetBinContent( 8, 12.0);
  h_qed->SetBinContent( 9,  6.0);
  h_qed->SetBinContent(10,  100);

  h_fcc->SetBinContent( 1,  0.1);
  h_fcc->SetBinContent( 2,  0.1);
  h_fcc->SetBinContent( 3,  1.0);
  h_fcc->SetBinContent( 4,  4.0);
  h_fcc->SetBinContent( 5,  1.0);
  h_fcc->SetBinContent( 6,  1.0);
  h_fcc->SetBinContent( 7,  0.5);
  h_fcc->SetBinContent( 8,  0.6);
  h_fcc->SetBinContent( 9,  0.5);
  h_fcc->SetBinContent(10,    1);

  // normalize to LEP
  h_qed->Divide(h_lep);
  h_fcc->Divide(h_lep);
  h_lep->Divide(h_lep);
  // Improvement factor
  TH1D *h_rat = (TH1D*)h_qed->Clone("h_rat");
  h_rat->Divide(h_fcc);

  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex(); CaptT->SetNDC(); // !!!
  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cFCCee1 = new TCanvas("cFCCee1","cFCCee1", gXcanv,  gYcanv,   1200,  800);
  //                            Name    Title            xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  gXcanv += 50; gYcanv += 50;
  cFCCee1->SetFillColor(10);
  //////////////////////////////////////////////
  cFCCee1->cd(1);
  gPad->SetLogy(); // !!!!!!

  TH1D *Hst=h_lep;

  Hst->SetStats(0);
  Hst->SetTitle(0);
  Hst->GetXaxis()->SetTitle("Observable");
  Hst->SetMaximum(1.0);
  Hst->SetMinimum(0.0003);

  Hst->SetLineColor(kBlack);    // blue line
  Hst->SetLineWidth(6);
  Hst->SetFillColor(kBlue);    // blues fill solid
  Hst->DrawCopy("h");

  h_qed->SetLineColor(kBlack);
  h_qed->SetLineWidth(6);
  h_qed->SetFillColor(kGreen);
  h_qed->DrawCopy("hsame");

  CaptT->DrawLatex(0.12,0.95,"Induced QED error in LEP pseudo-observables");
  //-------------------------------------
  //
  cFCCee1->cd();
  cFCCee1->SaveAs("cFCCee1.pdf");
//
}// FigFCCee1


///////////////////////////////////////////////////////////////////////////////////
void FigFCCee2()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigFCCee2 =========================== "<<endl;
  TH1D *h_lep = (TH1D*)DiskFileB.Get("h_lep");
  TH1D *h_qed = (TH1D*)DiskFileB.Get("h_qed");
  TH1D *h_fcc = (TH1D*)DiskFileB.Get("h_fcc");
  TH1D *h_rat = (TH1D*)DiskFileB.Get("h_rat");

  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex(); CaptT->SetNDC(); // !!!
  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cFCCee2 = new TCanvas("cFCCee2","cFCCee2", gXcanv,  gYcanv,   1200,  800);
  //                            Name    Title            xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  gXcanv += 50; gYcanv += 50;
  cFCCee2->SetFillColor(10);
  //////////////////////////////////////////////
  cFCCee2->cd(1);
  gPad->SetLogy(); // !!!!!!

  TH1D *Hst=h_lep;

  Hst->SetStats(0);
  Hst->SetTitle(0);
  Hst->GetXaxis()->SetTitle("Observable");
  Hst->SetMaximum(1.0);
  Hst->SetMinimum(0.0003);

  Hst->SetLineColor(kBlack);    // blue line
  Hst->SetLineWidth(6);
  Hst->SetFillColor(kBlue);    // blues fill solid
  Hst->DrawCopy("h");

  h_qed->SetLineColor(kBlack);
  h_qed->SetLineWidth(6);
  h_qed->SetFillColor(kGreen);
  h_qed->DrawCopy("hsame");

  h_fcc->SetLineColor(kBlack);
  h_fcc->SetLineWidth(6);
  h_fcc->SetFillColor(kRed);
  h_fcc->DrawCopy("hsame");


  CaptT->DrawLatex(0.12,0.95,"Current QED precision vs. FCCee exp. error");
  //-------------------------------------
  //
  cFCCee2->cd();
  cFCCee2->SaveAs("cFCCee2.pdf");
//
}// FigFCCee2


///////////////////////////////////////////////////////////////////////////////////
void FigFCCee3()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigFCCee3 =========================== "<<endl;
  TH1D *h_lep = (TH1D*)DiskFileB.Get("h_lep");
  TH1D *h_qed = (TH1D*)DiskFileB.Get("h_qed");
  TH1D *h_fcc = (TH1D*)DiskFileB.Get("h_fcc");
  TH1D *h_rat = (TH1D*)DiskFileB.Get("h_rat");

  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex(); CaptT->SetNDC(); // !!!
  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cFCCee3 = new TCanvas("cFCCee3","cFCCee3", gXcanv,  gYcanv,   1200,  800);
  //                            Name    Title            xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  gXcanv += 50; gYcanv += 50;
  cFCCee3->SetFillColor(10);
  //////////////////////////////////////////////
  cFCCee3->cd(1);
  gPad->SetLogy(); // !!!!!!

  TH1D *Hst=h_rat;

  Hst->SetStats(0);
  Hst->SetTitle(0);
  Hst->GetXaxis()->SetTitle("Observable");
  //Hst->SetMaximum(1.0);
  //Hst->SetMinimum(0.0003);

  Hst->SetLineColor(kBlack);
  Hst->SetLineWidth(6);
  Hst->SetFillColor(kRed);
  Hst->DrawCopy("h");

  Hst->Scale(3);
  Hst->SetFillColor(0);
  Hst->SetLineStyle(2);
  Hst->DrawCopy("hsame");

  CaptT->DrawLatex(0.12,0.95,"Needed improvement for QED precision at FCCee");
  //-------------------------------------
  //
  cFCCee3->cd();
  cFCCee3->SaveAs("cFCCee3.pdf");
//
}// FigFCCee3


///////////////////////////////////////////////////////////////////////////////////
void FigTempl()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigTempl =========================== "<<endl;
  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex(); CaptT->SetNDC(); // !!!
  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cTempl = new TCanvas("cTempl","cTempl", gXcanv,  gYcanv,   1200,  600);
  //                            Name    Title            xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  gXcanv += 50; gYcanv += 50;
  cTempl->SetFillColor(10);
  cTempl->Divide( 2,  0);
  //cTempl->Divide( 2,  0,     0.0,     0.0,   10);
  //              nx, ny, xmargin, ymargin, color
  //////////////////////////////////////////////
  cTempl->cd(1);
  CaptT->DrawLatex(0.12,0.95,"A_{FB}(v_{max}), ????");
  //-------------------------------------
  cTempl->cd(2);
  //
  cTempl->cd();
//
}// FigTempl


///////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  //++++++++++++++++++++++++++++++++++++++++
  TApplication theApp("theApp", &argc, argv);
  //++++++++++++++++++++++++++++++++++++++++
  //========== PLOTTING ==========
  FigPrag2();
  FigFCCee1();
  FigFCCee2();
  FigFCCee3();
  // Template empty canvas  with 2 figures
  //FigTempl();
  //++++++++++++++++++++++++++++++++++++++++
  DiskFileB.ls();
  DiskFileB.Write();
  DiskFileB.Close();
  //++++++++++++++++++++++++++++++++++++++++
  theApp.Run();
  //++++++++++++++++++++++++++++++++++++++++
}
