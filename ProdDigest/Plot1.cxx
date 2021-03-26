#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

#include <math.h>
#include <TLorentzVector.h>
#include <TLine.h>
#include <TArrow.h>
#include <TLatex.h>
#include <TSystem.h>

#include "TROOT.h"
#include "TCanvas.h"
#include "TF2.h"
//#include "TH2.h"
#include "TGaxis.h"
#include "TApplication.h"
#include "TMarker.h"
#include "TFile.h"

#include "HisNorm.h"
//=============================================================================
//  ROOT  ROOT ROOT   ROOT  ROOT  ROOT  ROOT  ROOT  ROOT  ROOT   ROOT   ROOT
//=============================================================================
TFile *DiskFileA;
TFile *DiskFileF;
TFile *DiskFileB;

TString FileA= "../ProdRun/work1/histo.root";

TString FileF= "../ProdRun/workFoam/histo.root";

///////////////////////////////////////////////////////////////////////////////////
//              GLOBAL stuff
///////////////////////////////////////////////////////////////////////////////////
double gCMSene, gNevTot; // from KKMC run
char   gTextEne[100], gTextNev[100], gTextNev2[100];
int    kGold=kOrange-3, kBrune=46, kPine=kGreen+3;
//
float  gXcanv = 0, gYcanv = 0, gDcanv = 30;
///////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////
void PlotSame(TH1D *&HST, double &ycapt, Int_t kolor, TString opis)
{
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.035);
  HST->SetLineColor(kolor);
  HST->DrawCopy("hsame");      // Magenta
  CaptT->SetTextColor(kolor);
  ycapt += -0.04;
  //double yy=ycapt;
  CaptT->DrawLatex(0.40,ycapt, opis);
}// PlotSame


///////////////////////////////////////////////////////////////////////////////////
void PlotSame2(TH1D *HST, double &xcapt, double &ycapt, Int_t kolor, double xx,  TString label,  TString opis)
{
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.035);
  HST->SetLineColor(kolor);
  HST->DrawCopy("hsame");
  CaptT->SetTextColor(kolor);
  ycapt += -0.04;
  CaptT->DrawLatex(xcapt,ycapt, opis);
  CaptT->DrawLatex(xcapt-0.05,ycapt, label);
  //
  TLatex *CaptS = new TLatex();
  CaptS->SetTextSize(0.040);
  CaptS->SetTextAlign(21);
  CaptS->SetTextColor(kolor);
  int ib = HST->FindBin(xx);
  double yy= HST->GetBinContent(ib);
  CaptS->DrawLatex(xx,yy,label);
}// PlotSame2



///////////////////////////////////////////////////////////////////////////////////
void HistNormalize(){
  //
  cout<<"----------------------------- HistNormalize ------------------------------------"<<endl;
  DiskFileA->ls("");
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA->Get("HST_KKMC_NORMA");
  //
  //HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA->Get("hst_WtMain") );
  //HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA->Get("hst_WtFoam") );
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA->Get("hst_nPhot") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA->Get("hst_CosTheta") );
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA->Get("hst_vvTrue") );
 //
  //HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA->Get("sca_vTcPR_Ceex2") );
  //HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA->Get("sca_vTcPR_Ceex2n") );
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA->Get("sca_vTcPR_Eex2") );

  //[[[[[[[[[[[[[[[[[[[[[
  //cout<<"----------------DiskFileA->GetListOfKeys--------------------"<<endl;
  //DiskFileA->GetListOfKeys()->Print();
  //cout<<"----------------DiskFileA->ShowStreamerInfo-------------------------"<<endl;
  //DiskFileA->ShowStreamerInfo();
  //cout<<"----------------DiskFileF->ls--------------------"<<endl;
  //DiskFileF->ls("");
  //cout<<"----------------DiskFileF->ShowStreamerInfo-------------------------"<<endl;
  //DiskFileF->ShowStreamerInfo();
  //]]]]]]]]]]]]]]]]]]]]]]
//----------------------
  TH1D    *HST_FOAM_NORMA4= (TH1D*)DiskFileF->Get("HST_FOAM_NORMA4");
  HisNorm1(HST_FOAM_NORMA4, (TH1D*)DiskFileF->Get("HST_weight4") );
  HisNorm1(HST_FOAM_NORMA4, (TH1D*)DiskFileF->Get("HST_vv_eex2") );
  HisNorm2(HST_FOAM_NORMA4, (TH2D*)DiskFileF->Get("SCA_vTcPR_Eex2") );
//------------------------------
  TH1D    *HST_FOAM_NORMA6= (TH1D*)DiskFileF->Get("HST_FOAM_NORMA6");
  HisNorm1(HST_FOAM_NORMA6, (TH1D*)DiskFileF->Get("HST_weight6") );

}// HistNormalize


///////////////////////////////////////////////////////////////////////////////////
void ReMakeMChisto(){
	// Some KKMC histos are preprocessed
	//------------------------------------------------------------------------
  cout<<"==================================================================="<<endl;
  cout<<"================ ReMakeMChisto  BEGIN  ============================"<<endl;
  TH1D *HST_KKMC_NORMA  = (TH1D*)DiskFileA->Get("HST_KKMC_NORMA");
  TH1D *HST_FOAM_NORMA4 = (TH1D*)DiskFileF->Get("HST_FOAM_NORMA4");
  //****************************************************************************************
  // Pure MC reprocessing part
  //
  TH2D *sca_vTcPR_Eex2  = (TH2D*)DiskFileA->Get("sca_vTcPR_Eex2");
  //TH2D *sca_vTcPR_Ceex2 = (TH2D*)DiskFileA->Get("sca_vTcPR_Ceex2");
  //TH2D *sca_vTcPR_Ceex2n = (TH2D*)DiskFileA->Get("sca_vTcPR_Ceex2n");

  ///****************************************************************************************
  /// Distributions of v=vTrue with unlimited c=cos(theta)
  //  without cutoff on c=cos(thetaPRD)
  int nbMax=0;   // cosThetaMax = 1.0
  TH1D                   *hTot_vTcPR_Eex2, *hAfb_vTcPR_Eex2;
  ProjV( sca_vTcPR_Eex2,  hTot_vTcPR_Eex2,  hAfb_vTcPR_Eex2, nbMax);  //!!!!
  hTot_vTcPR_Eex2->SetName("hTot_vTcPR_Eex2");
  hAfb_vTcPR_Eex2->SetName("hAfb_vTcPR_Eex2");
  nbMax=0;   // cosThetaMax = 1.0
  /*
  TH1D                    *HTot_vTcPR_Ceex2, *HAfb_vTcPR_Ceex2;
  ProjV( sca_vTcPR_Ceex2,  HTot_vTcPR_Ceex2,  HAfb_vTcPR_Ceex2, nbMax);  //!!!!
  HTot_vTcPR_Ceex2->SetName("HTot_vTcPR_Ceex2");
  HAfb_vTcPR_Ceex2->SetName("HAfb_vTcPR_Ceex2");
  //
  // IFI off
  nbMax=0;   // cosThetaMax = 1.0
  TH1D                    *HTot_vTcPR_Ceex2n, *HAfb_vTcPR_Ceex2n;
  ProjV( sca_vTcPR_Ceex2n,  HTot_vTcPR_Ceex2n,  HAfb_vTcPR_Ceex2n, nbMax);  //!!!!
  HTot_vTcPR_Ceex2n->SetName("HTot_vTcPR_Ceex2n");
  HAfb_vTcPR_Ceex2n->SetName("HAfb_vTcPR_Ceex2n");
  */
  //
  ///****************************************************************************************
  //  dsigma/dv unlimited cos(theta)
  TH1D *hPro_vT_Eex2;
  ProjX1(sca_vTcPR_Eex2, hPro_vT_Eex2);
  hPro_vT_Eex2->SetName("hPro_vT_Eex2");

  /*
  TH1D *Hpro_vT_Ceex2;
  ProjX1(sca_vTcPR_Ceex2, Hpro_vT_Ceex2);
  Hpro_vT_Ceex2->SetName("Hpro_vT_Ceex2");

  //  dsigma/dv unlimited cos(theta)
  TH1D *Hpro_vT_Ceex2n;
  ProjX1(sca_vTcPR_Ceex2, Hpro_vT_Ceex2n);
  Hpro_vT_Ceex2->SetName("Hpro_vT_Ceex2n");
*/
//-------------------------------------------------------------------------------------------
//-----------------------------FOAM----------------------------------------------------------
//-------------------------------------------------------------------------------------------
  TH2D *SCA_vTcPR_Eex2  = (TH2D*)DiskFileF->Get("SCA_vTcPR_Eex2");
  //TH2D *SCA_vTcPR_Ceex2 = (TH2D*)DiskFileF->Get("SCA_vTcPR_Ceex2");
  //TH2D *SCA_vTcPR_Ceex2n = (TH2D*)DiskFileF->Get("SCA_vTcPR_Ceex2n");

  ///****************************************************************************************
  /// Distributions of v=vTrue with unlimited c=cos(theta)
  //  without cutoff on c=cos(thetaPRD)
  nbMax=0;   // cosThetaMax = 1.0
  TH1D                   *HTot_vTcPR_Eex2, *HAfb_vTcPR_Eex2;
  ProjV( SCA_vTcPR_Eex2,  HTot_vTcPR_Eex2,  HAfb_vTcPR_Eex2, nbMax);  //!!!!
  HTot_vTcPR_Eex2->SetName("HTot_vTcPR_Eex2");
  HAfb_vTcPR_Eex2->SetName("HAfb_vTcPR_Eex2");

  ///****************************************************************************************
  //  dsigma/dv unlimited cos(theta)
  TH1D *HPro_vT_Eex2;
  ProjX1(SCA_vTcPR_Eex2, HPro_vT_Eex2);
  HPro_vT_Eex2->SetName("HPro_vT_Eex2");


  cout<<"================ ReMakeMChisto ENDs  ============================="<<endl;
  cout<<"==================================================================="<<endl;
}//RemakeMChisto


///////////////////////////////////////////////////////////////////////////////////
void FigInfo()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigInfo =========================== "<<endl;
  char capt1[100];
  double CMSene=100.0;
  sprintf(capt1,"#sqrt{s} =%4.0fGeV", CMSene);
  //
  TH1D *hst_WtMain    = (TH1D*)DiskFileA->Get("hst_WtMain");
  TH1D *hst_WtFoam    = (TH1D*)DiskFileA->Get("hst_WtFoam");

  TH1D *hst_nPhot     = (TH1D*)DiskFileA->Get("hst_nPhot");
  TH1D *hst_CosTheta  = (TH1D*)DiskFileA->Get("hst_CosTheta");
  //
  TH1D *HST_weight4    = (TH1D*)DiskFileF->Get("HST_weight4");
  TH1D *HST_weight6    = (TH1D*)DiskFileF->Get("HST_weight6");
  //------------------------------------------------------------------------
  //////////////////////////////////////////////
  TLatex *CaptE = new TLatex();
  CaptE->SetNDC(); // !!!
  CaptE->SetTextAlign(23);
//  CaptE->SetTextSize(0.055);
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
//  CaptT->SetTextSize(0.060);
  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigInfo = new TCanvas("cFigInfo","FigInfo: general info ",  gXcanv,  gYcanv,    1000,  800);
  //                            Name    Title               xoff,yoff, WidPix,HeiPix
  gXcanv += gDcanv; gYcanv += gDcanv;
  cFigInfo->SetFillColor(10);
  ////////////////////////////////////////////////////////////////////////////////
  cFigInfo->Divide( 2,  2);
  ///==========plot1==============
  cFigInfo->cd(1);
  gPad->SetLogy(); // !!!!!!
  ///
  CaptT->DrawLatex(0.10,0.95,"weight distribution");
  hst_WtMain->GetXaxis()->SetLabelSize(0.05);
  hst_WtMain->DrawCopy("h");

  hst_WtFoam->SetLineColor(kRed);
  hst_WtFoam->DrawCopy("hsame");

  //==========plot2==============
  cFigInfo->cd(2);
  //gPad->SetLogy(); // !!!!!!
  TH1D *HST = hst_nPhot;
  //HST->SetStats(0);
  //HST->SetMinimum( 1e-3*HST->GetMaximum());
  HST->SetTitle(0);
  HST->SetLineColor(kBlue);
  HST->GetXaxis()->SetLabelSize(0.06);
  HST->DrawCopy("h");
  CaptT->DrawLatex(0.10,0.93,"No. of #gamma's. tagged (red) and All (blue)");
  CaptE->DrawLatex(0.70,0.85, capt1);

  ///==========plot3==============
  cFigInfo->cd(3);
  //-----------------------------
  gPad->SetLogy(); // !!!!!!
  HST = HST_weight4;

  //HST->SetTitle(0);
  //HST->SetStats(0);
  HST->SetLineColor(4);
  //HST->SetMinimum(0);

  HST->DrawCopy("h");
  HST_weight6->DrawCopy("hsame");

  //=================
  cFigInfo->cd();

  cFigInfo->SaveAs("cFigInfo.pdf");
}//FigInfo

///////////////////////////////////////////////////////////////////////////////////
void FigVplot()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigVplot =========================== "<<endl;
  char capt1[100];
  double CMSene=100.0;
  sprintf(capt1,"#sqrt{s} =%4.0fGeV", CMSene);
  //
  TH1D *hst_vvTrue    = (TH1D*)DiskFileA->Get("hst_vvTrue");    //KKMCee
  TH1D *hPro_vT_Eex2  = (TH1D*)DiskFileB->Get("hPro_vT_Eex2");  //KKMCee

  TH1D *HST_vv_eex2   = (TH1D*)DiskFileF->Get("HST_vv_eex2");   //KKeeFoam
  TH1D *HPro_vT_Eex2  = (TH1D*)DiskFileB->Get("HPro_vT_Eex2");  //KKeeFoam

  //------------------------------------------------------------------------
  //////////////////////////////////////////////
  TLatex *CaptE = new TLatex();
  CaptE->SetNDC(); // !!!
  CaptE->SetTextAlign(23);
//  CaptE->SetTextSize(0.055);
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.040);
  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigVplot = new TCanvas("cFigVplot","FigVplot: general info ",  gXcanv,  gYcanv,    1000,  500);
  //                            Name    Title               xoff,yoff, WidPix,HeiPix
  gXcanv += gDcanv; gYcanv += gDcanv;
  cFigVplot->SetFillColor(10);
  ////////////////////////////////////////////////////////////////////////////////
  cFigVplot->Divide( 2,  1);

  //==========plot1==============
  cFigVplot->cd(1);
  gPad->SetLogy(); // !!!!!!
  TH1D *HST; //
  HST = hst_vvTrue;   //
  HST = hPro_vT_Eex2; //

  HST->SetTitle(0);
  HST->SetStats(0);

  HST->GetXaxis()->SetTitle("v=1-#hat{s}/s");
  HST->DrawCopy("h");                  //KKMCee

  CaptT->DrawLatex(0.06,0.95, "#sigma(v_{max}) ");
  double ycapt = 0.80; double xcapt=0.40;
  CaptT->SetTextColor(kBlack); ycapt += -0.04;
  CaptT->DrawLatex(xcapt,ycapt, "e^{+}e^{-} -> #mu^{+} #mu^{-}, 100GeV");

  PlotSame2(hst_vvTrue,    xcapt, ycapt, kRed,      0.10, "(a2)", "  ISR EEX2, KKMCee ");
  PlotSame2(hPro_vT_Eex2,  xcapt, ycapt, kGreen,    0.20, "(b2)", "  ISR EEX2, KKMCee ");
  PlotSame2(HST_vv_eex2,   xcapt, ycapt, kBlue,     0.60, "(f2)", "  ISR EEX2, KKeeFoam ");
  PlotSame2(HPro_vT_Eex2,  xcapt, ycapt, kMagenta,  0.70, "(g2)", "  ISR EEX2, KKeeFoam ");

//==========plot1==============
  cFigVplot->cd(2);
  TH1D *RAT_Hh  = HstRatio("RAT_Hh",  hst_vvTrue,  HST_vv_eex2,  kBlue ); //

  HST =RAT_Hh;
  HST->SetTitle(0);
  HST->SetStats(0);
  HST->SetLineColor(kBlue);
  HST->SetMinimum(1-0.02);
  HST->SetMaximum(1+0.02);
//  HST->SetMinimum(1-0.2);
//  HST->SetMaximum(1+0.2);
  HST->GetXaxis()->SetTitle("v=1-#hat{s}/s");
  HST->DrawCopy("h");

  CaptT->DrawLatex(0.06,0.95, " (a2)/(f2) ");

  TH1D *hOne = (TH1D*)RAT_Hh->Clone("hOne");  // zero line
  for(int i=1; i <= hOne->GetNbinsX() ; i++) { hOne->SetBinContent(i, 1); hOne->SetBinError(i, 0);}
  hOne->SetLineColor(kBlack);
  hOne->DrawCopy("hsame");

  //=================
  cFigVplot->cd();

  cFigVplot->SaveAs("cFigVplot.pdf");
}//FigVplot

///////////////////////////////////////////////////////////////////////////////////
void FigVCplot()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigVCplot =========================== "<<endl;
  char capt1[100];
  double CMSene=100.0;
  sprintf(capt1,"#sqrt{s} =%4.0fGeV", CMSene);
  TH1D *hTot_vTcPR_Eex2 = (TH1D*)DiskFileB->Get("hTot_vTcPR_Eex2");
  TH1D *HTot_vTcPR_Eex2 = (TH1D*)DiskFileB->Get("HTot_vTcPR_Eex2");
  //------------------------------------------------------------------------
  //////////////////////////////////////////////
  TLatex *CaptE = new TLatex();
  CaptE->SetNDC(); // !!!
  CaptE->SetTextAlign(23);
//  CaptE->SetTextSize(0.055);
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.040);
  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigVCplot = new TCanvas("cFigVCplot","FigVCplot: general info ",  gXcanv,  gYcanv,    1000,  500);
  //                            Name    Title               xoff,yoff, WidPix,HeiPix
  gXcanv += gDcanv; gYcanv += gDcanv;
  cFigVCplot->SetFillColor(10);
  ////////////////////////////////////////////////////////////////////////////////
  cFigVCplot->Divide( 2,  1);

  //==========plot1==============
  cFigVCplot->cd(1);
  //gPad->SetLogy(); // !!!!!!
  TH1D *HST; //
  HST = hTot_vTcPR_Eex2; //

  HST->SetTitle(0);
  HST->SetStats(0);
  HST->SetLineColor(kRed);
  HST->SetMinimum(0);
  HST->GetXaxis()->SetTitle("v=1-#hat{s}/s");
  HST->GetYaxis()->SetTitle("  ");

  HST->DrawCopy("h");

  CaptT->DrawLatex(0.06,0.95, "#sigma(v_{max}) ");
  double ycapt = 0.70; double xcapt=0.40;
  CaptT->SetTextColor(kBlack); ycapt += -0.04;
  CaptT->DrawLatex(xcapt,ycapt, "e^{+}e^{-} -> #mu^{+} #mu^{-}, 100GeV");

  PlotSame2(hTot_vTcPR_Eex2,    xcapt, ycapt, kBlack,   0.10, "(a2)", "  ISR EEX2, KKMCee:  ");
  PlotSame2(HTot_vTcPR_Eex2,    xcapt, ycapt, kBlue,    0.40, "(f2)", "  ISR EEX2, KKeeFoam:");

//==========plot1==============
  cFigVCplot->cd(2);
  //  TH1D *RAT_Hh  = HstRatio("RAT_Hh",  hst_Mll_ceex2_B,  hst_Mll_ceex2_A,  kBlack ); //
  TH1D *RAT_Ha  = HstRatio("RAT_Ha",  hTot_vTcPR_Eex2,  HTot_vTcPR_Eex2,  kBlue ); //

  HST =RAT_Ha;
  HST->SetTitle(0);
  HST->SetStats(0);
  HST->SetLineColor(kBlue);
//  HST->SetMinimum(1-0.006);
//  HST->SetMaximum(1+0.006);
  HST->SetMinimum(1-0.01);
  HST->SetMaximum(1+0.01);

  HST->GetXaxis()->SetTitle("v=1-#hat{s}/s");
  HST->GetYaxis()->SetTitle("  ");

  HST->DrawCopy("h");

  CaptT->DrawLatex(0.06,0.95, "(a2)/(f2) ");

  TH1D *hOne = (TH1D*)RAT_Ha->Clone("hOne");  // zero line
  for(int i=1; i <= hOne->GetNbinsX() ; i++) { hOne->SetBinContent(i, 1); hOne->SetBinError(i, 0);}
  hOne->SetLineColor(kBlack);
  hOne->DrawCopy("hsame");

  //=================
  cFigVCplot->cd();

  cFigVCplot->SaveAs("cFigVCplot.pdf");
}//FigVCplot


///////////////////////////////////////////////////////////////////////////////////
void FigAFBvv()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigAFBvv =========================== "<<endl;
  char capt1[100];
  double CMSene=100.0;
  sprintf(capt1,"#sqrt{s} =%4.0fGeV", CMSene);
  TH1D *hAfb_vTcPR_Eex2= (TH1D*)DiskFileB->Get("hAfb_vTcPR_Eex2");
  TH1D *HAfb_vTcPR_Eex2= (TH1D*)DiskFileB->Get("HAfb_vTcPR_Eex2");
  //------------------------------------------------------------------------
  //////////////////////////////////////////////
  TLatex *CaptE = new TLatex();
  CaptE->SetNDC(); // !!!
  CaptE->SetTextAlign(23);
//  CaptE->SetTextSize(0.055);
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.040);
  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigAFBvv = new TCanvas("cFigAFBvv","FigAFBvv: general info ",  gXcanv,  gYcanv,    1000,  500);
  //                            Name    Title               xoff,yoff, WidPix,HeiPix
  gXcanv += gDcanv; gYcanv += gDcanv;
  cFigAFBvv->SetFillColor(10);
  ////////////////////////////////////////////////////////////////////////////////
  cFigAFBvv->Divide( 2,  1);

//==========plot1==============
  cFigAFBvv->cd(1);
  //gPad->SetLogy(); // !!!!!!
  TH1D *HST; //
  HST = hAfb_vTcPR_Eex2; //
  HST->SetTitle(0);
  HST->SetStats(0);
  HST->SetLineColor(kRed);
//  HST->SetMinimum(0);
  HST->GetXaxis()->SetTitle("v");
  HST->DrawCopy("h");

  CaptT->DrawLatex(0.06,0.95, "A_{FB}(v_{max}) ");
  double ycapt = 0.70; double xcapt=0.40;
  CaptT->SetTextColor(kBlack); ycapt += -0.04;
  CaptT->DrawLatex(xcapt,ycapt, "e^{+}e^{-} -> #mu^{+} #mu^{-}, 100GeV");
  PlotSame2(hAfb_vTcPR_Eex2,    xcapt, ycapt, kBlack,   0.30, "(a2)", "  ISR EEX2, KKMCee:  ");
  PlotSame2(HAfb_vTcPR_Eex2,    xcapt, ycapt, kBlue,    0.50, "(f2)", "  ISR EEX2, KKeeFoam:");

//==========plot1==============
  cFigAFBvv->cd(2);
  TH1D *afb_diff_Hh = HstDiff(  "afb_diff_h", hAfb_vTcPR_Eex2, HAfb_vTcPR_Eex2, kBlue); //
  HST = afb_diff_Hh;
  HST->SetTitle(0);
  HST->SetStats(0);
  HST->SetLineColor(kBlack);
  HST->SetMinimum(-0.006);
  HST->SetMaximum(+0.006);
  HST->GetXaxis()->SetTitle("v");
  HST->DrawCopy("h");

  CaptT->DrawLatex(0.06,0.95, "(a2)-(f2)");

  TH1D *hZero0 = (TH1D*)afb_diff_Hh->Clone("hZero0");  // zero line
  for(int i=1; i <= hZero0->GetNbinsX() ; i++) { hZero0->SetBinContent(i, 0); hZero0->SetBinError(i, 0);}
  hZero0->SetLineColor(kBlack);
  hZero0->DrawCopy("hsame");

  //=================
  cFigAFBvv->cd();

  cFigAFBvv->SaveAs("cFigAFBvv.pdf");

}//FigAFBvv



///////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{

  gSystem->Load("../SRCee/.libs/libKKee.so");     // needed ???
  gSystem->Load("../SRCee/.libs/libKKfm.so");     // NEEDED! why?

  DiskFileA = new TFile(FileA);
  DiskFileF = new TFile(FileF);
//
  DiskFileB = new TFile("Plot1.root","RECREATE","histograms");

  //++++++++++++++++++++++++++++++++++++++++
  TApplication theApp("theApp", &argc, argv);
  //++++++++++++++++++++++++++++++++++++++++
  /*
  /////////////////////////////////////////////////////////
  // Reading directly KKMC input (farming)
  int Nodes;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA->Get("HST_KKMC_NORMA");
  Nodes    = HST_KKMC_NORMA->GetBinContent(511);       // No of farm nodes (trick)
  gCMSene  = HST_KKMC_NORMA->GetBinContent(1)/Nodes;   // CMSene=xpar(1), farn adjusted
  gNevTot  = HST_KKMC_NORMA->GetEntries();             // MC statistics from KKMC
  sprintf(gTextEne,"#sqrt{s} =%4.2fGeV", gCMSene);
  sprintf(gTextNev,"KKMC:%10.2e events", gNevTot);
 */
  DiskFileB->cd();
  HistNormalize();     // Renormalization of MC histograms
  ReMakeMChisto();
  //========== PLOTTING ==========
  FigInfo();
  FigVplot();
  FigVCplot();
//  FigAFBvv();
 //++++++++++++++++++++++++++++++++++++++++
  DiskFileA->ls();
//
//  cout<< "CMSene[GeV] = "<< gCMSene<< endl;
//  cout<< "KKMC: No. of farm nodes="<< Nodes  << "  Tot no. of events = "<<gNevTot<< endl;
//
  //++++++++++++++++++++++++++++++++++++++++
  theApp.Run();
  //++++++++++++++++++++++++++++++++++++++++
}
