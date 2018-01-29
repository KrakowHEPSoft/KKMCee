
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
