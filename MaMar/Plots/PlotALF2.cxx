//    make PlotALF-run
//    Plots for new second paper on AFB

// This is only AFB study for several energies,
// renomalizing scattergrams is no necesssary!

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

#include "KKplot.h"
#include "HisNorm.h"

// ROOT headers
#include "TROOT.h"
#include "TFile.h"
//=============================================================================
//  ROOT  ROOT ROOT   ROOT  ROOT  ROOT  ROOT  ROOT  ROOT  ROOT   ROOT   ROOT
//=============================================================================
// Latest from /workKKMC

//TFile DiskFileA88("../workKKMC/histo.root");  // apr.2018
//TFile DiskFileA95("../workKKMC/histo.root");  // apr.2018

//TFile DiskFileA88("../workKKMC/histo.root_88GeV.new");  // recent
//TFile DiskFileA95("../workKKMC/histo.root_95GeV.new");  // recent

TFile DiskFileA88("../workKKMC/histo.root_88GeV_apr_13G");  // apr.2018
TFile DiskFileA95("../workKKMC/histo.root_95GeV_apr_38G");  // apr.2018

// current new histos
TFile DiskFileB("RhoAFB.root","RECREATE","histograms");
///////////////////////////////////////////////////////////////////////////////////
//              GLOBAL stuff
///////////////////////////////////////////////////////////////////////////////////
int    kGold=kOrange-3, kBrune=46, kPine=kGreen+3;
//
float  gXcanv = 0, gYcanv = 0;
//
int    gTogEne = 1;   // 10 GeV and MZ included
//int    gTogEne = 0; // 10 GeV and MZ exluded
//
//int    gTogle  = 0;  // excluding new data files
int    gTogle  = 1;  // including new data files
//
int  gNBmax =45;    // for |cos(theta)|<0.90
double gCosTheta=0.90;       // to be synchronized with gNbMax

//
//int    gNbMax   =0;          // for 100bins, default=0 for gCosTheta = 1.00
//double gCosTheta=1.00;       // to be synchronized with gNbMax
//
int    gNbMax2=0;            // for 50 bins, default=0 for gCosTheta = 1.00
//
//KKplot LibSem("KKplot");
///////////////////////////////////////////////////////////////////////////////////
//Double_t sqr( const Double_t x ){ return x*x;};


///////////////////////////////////////////////////////////////////////////////////
void PlotSame3(TH1D *HST, double &xcapt, double &ycapt, Int_t kolor, int kMarker, int kStyle,  TString opis)
{
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.035);
  HST->SetLineColor(kolor);
  HST->SetMarkerStyle(kMarker);
  HST->SetLineStyle(kStyle);
  HST->DrawCopy("hsame");
  CaptT->SetTextColor(kolor);
  CaptT->DrawLatex(xcapt,ycapt, opis);
//
  TMarker *marker3 = new TMarker(xcapt-0.05, ycapt+0.01, kMarker);
  marker3->SetMarkerSize(1.2);
  marker3->SetNDC();
  marker3->Draw();
//
  TLine *lineM = new TLine(xcapt-0.08, ycapt+0.01, xcapt-0.02, ycapt+0.01);
  lineM->SetLineStyle(kStyle);
  lineM->SetNDC(); // !!!
  lineM->Draw();
//
  ycapt += -0.047;
}// PlotSame3



///////////////////////////////////////////////////////////////////////////////////
void ReMakeHistoMC()
{
//------------------------------------------------------------------------
  cout<<" ========================= ReMakeHistoMC =========================== "<<endl;
/////////////////////////////////////////////////////////////////////////////////////////////////
// CEEX2-CEEX1   v-distributions
/////////////////////////////////////////////////////////////////////////////////////////////////
  TH1D *hst9_vA_Ceex1       = (TH1D*)DiskFileA95.Get("hst_vA_Ceex1");     // total CEEX2
  TH1D *hst9_vA_Ceex2       = (TH1D*)DiskFileA95.Get("hst_vA_Ceex2");     //
  TH1D *hst9_vA_Ceex21      = (TH1D*)DiskFileA95.Get("hst_vA_Ceex21");    //
//
  TH1D *hst9_vA_Ceex1_F     = (TH1D*)DiskFileA95.Get("hst_vA_Ceex1_F");   //
  TH1D *hst9_vA_Ceex2_F     = (TH1D*)DiskFileA95.Get("hst_vA_Ceex2_F");   //
  TH1D *hst9_vA_Ceex21_F    = (TH1D*)DiskFileA95.Get("hst_vA_Ceex21_F");  //
//
  TH1D *hst9_vB_Ceex1       = (TH1D*)DiskFileA95.Get("hst_vB_Ceex1");     //
  TH1D *hst9_vB_Ceex2       = (TH1D*)DiskFileA95.Get("hst_vB_Ceex2");     //
  TH1D *hst9_vB_Ceex21      = (TH1D*)DiskFileA95.Get("hst_vB_Ceex21");    //
//
  TH1D *hst9_vB_Ceex1_F     = (TH1D*)DiskFileA95.Get("hst_vB_Ceex1_F");   //
  TH1D *hst9_vB_Ceex2_F     = (TH1D*)DiskFileA95.Get("hst_vB_Ceex2_F");   //
  TH1D *hst9_vB_Ceex21_F    = (TH1D*)DiskFileA95.Get("hst_vB_Ceex21_F");  //

  TH1D *hst9_vB_Ceex2n      = (TH1D*)DiskFileA95.Get("hst_vB_Ceex2n");  //
  TH1D *hst9_vB_Ceex2n_F    = (TH1D*)DiskFileA95.Get("hst_vB_Ceex2n_F");  //
////////////////////////////////////////////////////////////////////////////////////////////////////
// AFB calculated directly
  TH1D *HAfb9_vB_Ceex2 = HstAFB("HAfb9_vB_Ceex2", hst9_vB_Ceex2_F, hst9_vB_Ceex2);
  TH1D *HAfb9_vB_Ceex1 = HstAFB("HAfb9_vB_Ceex1", hst9_vB_Ceex1_F, hst9_vB_Ceex1); // ???
//
  TH1D *HAfb9_vB_Ceex2n = HstAFB("HAfb9_vB_Ceex2n", hst9_vB_Ceex2n_F, hst9_vB_Ceex2n);

//
  TH1D *HAfb9_vA_Ceex2 = HstAFB("HAfb9_vA_Ceex2", hst9_vA_Ceex2_F, hst9_vA_Ceex2);
  TH1D *HAfb9_vA_Ceex1 = HstAFB("HAfb9_vA_Ceex1", hst9_vA_Ceex1_F, hst9_vA_Ceex1); // ???
  // Exact formula for AFB from weight differences
  TH1D *HAfb9_vB_Ceex21 = HstAFB4( "HAfb9_vB_Ceex21", hst9_vB_Ceex21_F, hst9_vB_Ceex21,
		                                              hst9_vB_Ceex2_F,  hst9_vB_Ceex2 );
  TH1D *HAfb9_vA_Ceex21 = HstAFB4( "HAfb9_vA_Ceex21", hst9_vA_Ceex21_F, hst9_vA_Ceex21,
		                                              hst9_vA_Ceex2_F,  hst9_vA_Ceex2 );
//
/////////////////////////////////////////////////////////////////////////////////////////////////
//  CEEX2-CEEX1 v-distributions
/////////////////////////////////////////////////////////////////////////////////////////////////
  TH1D *hst8_vA_Ceex1       = (TH1D*)DiskFileA88.Get("hst_vA_Ceex1");     // total CEEX2
  TH1D *hst8_vA_Ceex2       = (TH1D*)DiskFileA88.Get("hst_vA_Ceex2");     //
  TH1D *hst8_vA_Ceex21      = (TH1D*)DiskFileA88.Get("hst_vA_Ceex21");    //
//
  TH1D *hst8_vA_Ceex1_F     = (TH1D*)DiskFileA88.Get("hst_vA_Ceex1_F");   //
  TH1D *hst8_vA_Ceex2_F     = (TH1D*)DiskFileA88.Get("hst_vA_Ceex2_F");   //
  TH1D *hst8_vA_Ceex21_F    = (TH1D*)DiskFileA88.Get("hst_vA_Ceex21_F");  //
//
  TH1D *hst8_vB_Ceex1       = (TH1D*)DiskFileA88.Get("hst_vB_Ceex1");     //
  TH1D *hst8_vB_Ceex2       = (TH1D*)DiskFileA88.Get("hst_vB_Ceex2");     //
  TH1D *hst8_vB_Ceex21      = (TH1D*)DiskFileA88.Get("hst_vB_Ceex21");    //
//
  TH1D *hst8_vB_Ceex1_F     = (TH1D*)DiskFileA88.Get("hst_vB_Ceex1_F");   //
  TH1D *hst8_vB_Ceex2_F     = (TH1D*)DiskFileA88.Get("hst_vB_Ceex2_F");   //
  TH1D *hst8_vB_Ceex21_F    = (TH1D*)DiskFileA88.Get("hst_vB_Ceex21_F");  //

  TH1D *hst8_vB_Ceex2n      = (TH1D*)DiskFileA88.Get("hst_vB_Ceex2n");  //
  TH1D *hst8_vB_Ceex2n_F    = (TH1D*)DiskFileA88.Get("hst_vB_Ceex2n_F");  //
////////////////////////////////////////////////////////////////////////////////////////////////////
// AFB calculated directly
  TH1D *HAfb8_vB_Ceex2 = HstAFB("HAfb8_vB_Ceex2", hst8_vB_Ceex2_F, hst8_vB_Ceex2);
  TH1D *HAfb8_vB_Ceex1 = HstAFB("HAfb8_vB_Ceex1", hst8_vB_Ceex1_F, hst8_vB_Ceex1); // ???
//
   TH1D *HAfb8_vB_Ceex2n = HstAFB("HAfb8_vB_Ceex2n", hst8_vB_Ceex2n_F, hst8_vB_Ceex2n);
//
  TH1D *HAfb8_vA_Ceex2 = HstAFB("HAfb8_vA_Ceex2", hst8_vA_Ceex2_F, hst8_vA_Ceex2);
  TH1D *HAfb8_vA_Ceex1 = HstAFB("HAfb8_vA_Ceex1", hst8_vA_Ceex1_F, hst8_vA_Ceex1); // ???
// Exact formula for AFB from weight differences
  TH1D *HAfb8_vB_Ceex21 = HstAFB4( "HAfb8_vB_Ceex21", hst8_vB_Ceex21_F, hst8_vB_Ceex21,
    	                                              hst8_vB_Ceex2_F,  hst8_vB_Ceex2 );
  TH1D *HAfb8_vA_Ceex21 = HstAFB4( "HAfb8_vA_Ceex21", hst8_vA_Ceex21_F, hst8_vA_Ceex21,
  		                                              hst8_vA_Ceex2_F,  hst8_vA_Ceex2 );
////////////////////////////////////////////////////////////////////////////////////////////////////
// IFI
////////////////////////////////////////////////////////////////////////////////////////////////////
  TH1D *hst8_vB_Ceex2i       = (TH1D*)DiskFileA88.Get("hst_vB_Ceex2i");     //
  TH1D *hst8_vB_Ceex2i_F     = (TH1D*)DiskFileA88.Get("hst_vB_Ceex2i_F");   //
  TH1D *hst8_vA_Ceex2i       = (TH1D*)DiskFileA88.Get("hst_vA_Ceex2i");     //
  TH1D *hst8_vA_Ceex2i_F     = (TH1D*)DiskFileA88.Get("hst_vA_Ceex2i_F");   //
//
  TH1D *HAfb8_vA_Ceex2_IFI = HstAFB4( "HAfb8_vA_Ceex2_IFI",hst8_vA_Ceex2i_F, hst8_vA_Ceex2i,
  		                                                   hst8_vA_Ceex2_F,  hst8_vA_Ceex2 );
  TH1D *HAfb8_vB_Ceex2_IFI = HstAFB4( "HAfb8_vB_Ceex2_IFI",hst8_vB_Ceex2i_F, hst8_vB_Ceex2i,
  		                                                   hst8_vB_Ceex2_F,  hst8_vB_Ceex2 );
///////////////
  TH1D *hst8_vB_Ceex1i       = (TH1D*)DiskFileA88.Get("hst_vB_Ceex1i");     //
  TH1D *hst8_vB_Ceex1i_F     = (TH1D*)DiskFileA88.Get("hst_vB_Ceex1i_F");   //
  TH1D *hst8_vA_Ceex1i       = (TH1D*)DiskFileA88.Get("hst_vA_Ceex1i");     //
  TH1D *hst8_vA_Ceex1i_F     = (TH1D*)DiskFileA88.Get("hst_vA_Ceex1i_F");   //
//
  TH1D *HAfb8_vA_Ceex1_IFI = HstAFB4( "HAfb8_vA_Ceex1_IFI",hst8_vA_Ceex1i_F, hst8_vA_Ceex1i,
   	                                                       hst8_vA_Ceex1_F,  hst8_vA_Ceex1 );
  TH1D *HAfb8_vB_Ceex1_IFI = HstAFB4( "HAfb8_vB_Ceex1_IFI",hst8_vB_Ceex1i_F, hst8_vB_Ceex1i,
   	                                                       hst8_vB_Ceex1_F,  hst8_vB_Ceex1 );

////////////////////////////////////////////////////////////////////////////////////////////////////
// IFI
////////////////////////////////////////////////////////////////////////////////////////////////////
  TH1D *hst9_vB_Ceex2i       = (TH1D*)DiskFileA95.Get("hst_vB_Ceex2i");     //
  TH1D *hst9_vB_Ceex2i_F     = (TH1D*)DiskFileA95.Get("hst_vB_Ceex2i_F");   //
  TH1D *hst9_vA_Ceex2i       = (TH1D*)DiskFileA95.Get("hst_vA_Ceex2i");     //
  TH1D *hst9_vA_Ceex2i_F     = (TH1D*)DiskFileA95.Get("hst_vA_Ceex2i_F");   //
//
  TH1D *HAfb9_vA_Ceex2_IFI = HstAFB4( "HAfb9_vA_Ceex2_IFI",hst9_vA_Ceex2i_F, hst9_vA_Ceex2i,
  		                                                   hst9_vA_Ceex2_F,  hst9_vA_Ceex2 );
  TH1D *HAfb9_vB_Ceex2_IFI = HstAFB4( "HAfb9_vB_Ceex2_IFI",hst9_vB_Ceex2i_F, hst9_vB_Ceex2i,
                                                           hst9_vB_Ceex2_F,  hst9_vB_Ceex2 );
//////////////
  TH1D *hst9_vB_Ceex1i       = (TH1D*)DiskFileA95.Get("hst_vB_Ceex1i");     //
  TH1D *hst9_vB_Ceex1i_F     = (TH1D*)DiskFileA95.Get("hst_vB_Ceex1i_F");   //
  TH1D *hst9_vA_Ceex1i       = (TH1D*)DiskFileA95.Get("hst_vA_Ceex1i");     //
  TH1D *hst9_vA_Ceex1i_F     = (TH1D*)DiskFileA95.Get("hst_vA_Ceex1i_F");   //
//
  TH1D *HAfb9_vA_Ceex1_IFI = HstAFB4( "HAfb9_vA_Ceex1_IFI",hst9_vA_Ceex1i_F, hst9_vA_Ceex1i,
  		                                                   hst9_vA_Ceex1_F,  hst9_vA_Ceex1 );
  TH1D *HAfb9_vB_Ceex1_IFI = HstAFB4( "HAfb9_vB_Ceex1_IFI",hst9_vB_Ceex1i_F, hst9_vB_Ceex1i,
                                                           hst9_vB_Ceex1_F,  hst9_vB_Ceex1 );

  cout<<" ===================== End of ReMakeHistoMC ======================== "<<endl;

}// ReMakeHistoMC


///////////////////////////////////////////////////////////////////////////////////
void KKsemMakeHisto(){
  // Here we produce semianalytical plots using KKsem program, No plotting
  //------------------------------------------------------------------------
  cout<<"==================================================================="<<endl;
  cout<<"================ KKsem MakeHisto  BEGIN ============================"<<endl;
  //
  long KF=13; // muon
  long KeyDis, KeyFob;
  char chak[5];
  KeyFob=    0; // With EW (BornV_Dizet) With integration OK!
//------------------------------------------------------------------------
// Tempate is for vmax<0.20 and 100 bins in cos(theta)
  TH1D *hstVtemplate = (TH1D*)DiskFileB.Get("HAfb8_vB_Ceex2");
//------------------------------------------------------------------------
//   MuMu  Sigma(vmax) and AFB(vmax) with ulimited c=cos(theta)
//
  sprintf(chak,"XCHI2");  // ISR*FSR Mff
  KeyDis = 302302;        // ISR*FSR O(alf2)
//
//  sprintf(chak,"VCHI2");  // ISR only Mff
//  KeyDis = 302;   // ISR O(alf2)
  //  KeyDis = 303;   // ISR O(alf3)

  // ************ 95GeV ************
  KKplot LibSem9("LibSem9");
  LibSem9.Initialize(DiskFileA95);  // for non-farm case
  //
  TH1D *vcum_KKsem_95 =(TH1D*)hstVtemplate->Clone("vcum_KKsem_95");
  TH1D *afbv_KKsem_95 =(TH1D*)hstVtemplate->Clone("afbv_KKsem_95");
  LibSem9.VVmake( vcum_KKsem_95, afbv_KKsem_95, KF, chak, KeyDis, KeyFob, gCosTheta);
  // ************  88GeV ************
  KKplot LibSem8("LibSem8");
  LibSem8.Initialize(DiskFileA88);  // for non-farm case
  //
  TH1D *vcum_KKsem_88 =(TH1D*)hstVtemplate->Clone("vcum_KKsem_88");
  TH1D *afbv_KKsem_88 =(TH1D*)hstVtemplate->Clone("afbv_KKsem_88");
  LibSem8.VVmake( vcum_KKsem_88, afbv_KKsem_88, KF, chak, KeyDis, KeyFob, gCosTheta);

  cout<<"================ KKsem MakeHisto ENDs ============================="<<endl;
  cout<<"==================================================================="<<endl;
//------------------------------------------------------------------------
//------------------------------------------------------------------------
}//  KKsemMakeHisto




///////////////////////////////////////////////////////////////////////////////////
void FigExp1()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigExp1 =========================== "<<endl;

  TH1D *HAfb9_vA_Ceex2  = (TH1D*)DiskFileB.Get("HAfb9_vA_Ceex2");
  TH1D *HAfb9_vA_Ceex1  = (TH1D*)DiskFileB.Get("HAfb9_vA_Ceex1");
  TH1D *HAfb9_vA_Ceex21 = (TH1D*)DiskFileB.Get("HAfb9_vA_Ceex21");

  TH1D *HAfb9_vB_Ceex2  = (TH1D*)DiskFileB.Get("HAfb9_vB_Ceex2");
  TH1D *HAfb9_vB_Ceex1  = (TH1D*)DiskFileB.Get("HAfb9_vB_Ceex1");
  TH1D *HAfb9_vB_Ceex21 = (TH1D*)DiskFileB.Get("HAfb9_vB_Ceex21");
  TH1D *HAfb9_vB_Ceex2n = (TH1D*)DiskFileB.Get("HAfb9_vB_Ceex2n");

  TH1D *afbv_KKsem_95  = (TH1D*)DiskFileB.Get("afbv_KKsem_95");

  TLatex *CaptNDC = new TLatex(); CaptNDC->SetNDC(); // !!!
  CaptNDC->SetTextSize(0.037);
  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cExp1 = new TCanvas("cExp1","cExp1", gXcanv,  gYcanv,   1200,  600);
  //                            Name    Title            xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  gXcanv += 50; gYcanv += 50;
  cExp1->SetFillColor(10);
  cExp1->Divide( 2,  0);
  //cExp1->Divide( 2,  0,     0.0,     0.0,   10);
  //              nx, ny, xmargin, ymargin, color
  //////////////////////////////////////////////
  cExp1->cd(1);
  TH1D *HST1;
  HST1 =HAfb9_vB_Ceex2;

  HST1 =HAfb9_vB_Ceex2n;


  HST1->SetStats(0);
  HST1->SetTitle(0);
  HST1->GetXaxis()->SetNdivisions(8);
  HST1->SetMaximum(0.350); HST1->SetMinimum(0.250); // 95GeV
  HST1->DrawCopy("h");

  double ycapt =0.91, xcapt=0.3;
  CaptNDC->DrawLatex(xcapt-0.1,ycapt," #sqrt{s} = 94.3GeV, Ceex2, IFI on");
  ycapt += -0.047;
  PlotSame3(HAfb9_vB_Ceex2,     xcapt, ycapt, kBlue,      25, 1, "(B) Idealized");
  PlotSame3(HAfb9_vA_Ceex2,     xcapt, ycapt, kBlack,     24, 1, "(A) Realistic");

  PlotSame3(HAfb9_vB_Ceex2n,    xcapt, ycapt, kRed,       27, 1, "(B) noIFI vtr");
  PlotSame3(afbv_KKsem_95,     xcapt, ycapt, kGreen,     24, 1, "(X) KKsem");

  CaptNDC->DrawLatex(0.04,0.95,"A_{FB}(v)  ");
  CaptNDC->DrawLatex(0.60,0.02,"v_{max} ");

  ///////////////////////////////////////////////
  cExp1->cd(2);
  TH1D *HST2;
  HST2 =HAfb9_vB_Ceex21;

  HST2->SetStats(0);
  HST2->SetTitle(0);
  HST2->GetXaxis()->SetNdivisions(8);
  HST2->SetMaximum(1.1e-5); HST2->SetMinimum(-1.1e-5);
  HST2->DrawCopy("h");

  ycapt =0.91, xcapt=0.3;
  CaptNDC->DrawLatex(xcapt-0.1,ycapt," #sqrt{s} = 94.3GeV, Ceex2-Ceex1, IFI on");
  ycapt += -0.047;
  PlotSame3(HAfb9_vB_Ceex21,   xcapt, ycapt, kBlue,     25, 1, "(B) Idealized");
  PlotSame3(HAfb9_vA_Ceex21,   xcapt, ycapt, kBlack,    25, 1, "(A) Realistic");

  CaptNDC->DrawLatex(0.04,0.95,"#delta A_{FB}(v)  ");
  CaptNDC->DrawLatex(0.60,0.02,"v_{max} ");

  //
  cExp1->cd();
  //
  cExp1->SaveAs("cExp1.pdf");
//
}// FigExp1


///////////////////////////////////////////////////////////////////////////////////
void FigExp2()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigExp2 =========================== "<<endl;

  TH1D *HAfb8_vA_Ceex2  = (TH1D*)DiskFileB.Get("HAfb8_vA_Ceex2");
  TH1D *HAfb8_vA_Ceex1  = (TH1D*)DiskFileB.Get("HAfb8_vA_Ceex1");
  TH1D *HAfb8_vA_Ceex21 = (TH1D*)DiskFileB.Get("HAfb8_vA_Ceex21");

  TH1D *afbv_KKsem_88  = (TH1D*)DiskFileB.Get("afbv_KKsem_88");

  TH1D *HAfb8_vB_Ceex2  = (TH1D*)DiskFileB.Get("HAfb8_vB_Ceex2");
  TH1D *HAfb8_vB_Ceex1  = (TH1D*)DiskFileB.Get("HAfb8_vB_Ceex1");
  TH1D *HAfb8_vB_Ceex21 = (TH1D*)DiskFileB.Get("HAfb8_vB_Ceex21");
  TH1D *HAfb8_vB_Ceex2n = (TH1D*)DiskFileB.Get("HAfb8_vB_Ceex2n");

  TLatex *CaptNDC = new TLatex(); CaptNDC->SetNDC(); // !!!
  CaptNDC->SetTextSize(0.037);
  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cExp2 = new TCanvas("cExp2","cExp2", gXcanv,  gYcanv,   1200,  600);
  //                            Name    Title            xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  gXcanv += 50; gYcanv += 50;
  cExp2->SetFillColor(10);
  cExp2->Divide( 2,  0);
  //cExp2->Divide( 2,  0,     0.0,     0.0,   10);
  //              nx, ny, xmargin, ymargin, color
  //////////////////////////////////////////////
  cExp2->cd(1);
  TH1D *HST1;
  HST1 =HAfb8_vB_Ceex2;

  HST1->SetStats(0);
  HST1->SetTitle(0);
  HST1->GetXaxis()->SetNdivisions(8);
  HST1->SetMaximum(-0.215); HST1->SetMinimum(-0.290); // 95GeV
  HST1->DrawCopy("h");

  double ycapt =0.91, xcapt=0.3;
  CaptNDC->DrawLatex(xcapt-0.1,ycapt," #sqrt{s} = 87.90GeV, Ceex2, IFI on");
  ycapt += -0.047;
  PlotSame3(HAfb8_vB_Ceex2,     xcapt, ycapt, kBlue,      25, 1, "(B) Idealized");
  PlotSame3(HAfb8_vA_Ceex2,     xcapt, ycapt, kBlack,     24, 1, "(A) Realistic");

  PlotSame3(HAfb8_vB_Ceex2n,    xcapt, ycapt, kRed,       27, 1, "(B) noIFI ");
  PlotSame3(afbv_KKsem_88,      xcapt, ycapt, kGreen,     24, 1, "(X) KKsem");

  CaptNDC->DrawLatex(0.04,0.95,"A_{FB}(v)  ");
  CaptNDC->DrawLatex(0.60,0.02,"v_{max} ");

  ///////////////////////////////////////////////
  cExp2->cd(2);
  TH1D *HST2;
  HST2 =HAfb8_vB_Ceex21;

  HST2->SetStats(0);
  HST2->SetTitle(0);
  HST2->GetXaxis()->SetNdivisions(8);
  HST2->SetMaximum(1.1e-5); HST2->SetMinimum(-1.1e-5);
  HST2->DrawCopy("h");

  ycapt =0.91, xcapt=0.3;
  CaptNDC->DrawLatex(xcapt-0.1,ycapt," #sqrt{s} = 87.90GeV, Ceex2-Ceex1, IFI on");
  ycapt += -0.047;
  PlotSame3(HAfb8_vB_Ceex21,   xcapt, ycapt, kBlue,     25, 1, "(B) Idealized");
  PlotSame3(HAfb8_vA_Ceex21,   xcapt, ycapt, kBlack,    25, 1, "(A) Realistic");

  CaptNDC->DrawLatex(0.04,0.95,"#delta A_{FB}(v)  ");
  CaptNDC->DrawLatex(0.60,0.02,"v_{max} ");

  //
  cExp2->cd();
//
  cExp2->SaveAs("cExp2.pdf");
//
}// FigExp2



///////////////////////////////////////////////////////////////////////////////////
void FigExp3()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigExp3 =========================== "<<endl;

  TH1D *HAfb8_vA_Ceex2  = (TH1D*)DiskFileB.Get("HAfb8_vA_Ceex2");
  TH1D *HAfb8_vA_Ceex21 = (TH1D*)DiskFileB.Get("HAfb8_vA_Ceex21");
  TH1D *HAfb8_vB_Ceex2  = (TH1D*)DiskFileB.Get("HAfb8_vB_Ceex2");
  TH1D *HAfb8_vB_Ceex21 = (TH1D*)DiskFileB.Get("HAfb8_vB_Ceex21");

  TH1D *HAfb9_vA_Ceex2  = (TH1D*)DiskFileB.Get("HAfb9_vA_Ceex2");
  TH1D *HAfb9_vA_Ceex21 = (TH1D*)DiskFileB.Get("HAfb9_vA_Ceex21");
  TH1D *HAfb9_vB_Ceex2  = (TH1D*)DiskFileB.Get("HAfb9_vB_Ceex2");
  TH1D *HAfb9_vB_Ceex21 = (TH1D*)DiskFileB.Get("HAfb9_vB_Ceex21");

  TH1D *HAfbDel_B_Ceex2 = HstDiff( "HAfbDel_B_Ceex2", HAfb9_vB_Ceex2, HAfb8_vB_Ceex2, kBlue);
  TH1D *HAfbDel_A_Ceex2 = HstDiff( "HAfbDel_A_Ceex2", HAfb9_vA_Ceex2, HAfb8_vA_Ceex2, kBlue);

  TH1D *HAfbDel_B_Ceex21 = HstDiff( "HAfbDel_B_Ceex21", HAfb9_vB_Ceex21, HAfb8_vB_Ceex21, kBlue);
  TH1D *HAfbDel_A_Ceex21 = HstDiff( "HAfbDel_A_Ceex21", HAfb9_vA_Ceex21, HAfb8_vA_Ceex21, kBlue);


  TLatex *CaptNDC = new TLatex(); CaptNDC->SetNDC(); // !!!
  CaptNDC->SetTextSize(0.037);
  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cExp3 = new TCanvas("cExp3","cExp3", gXcanv,  gYcanv,   1200,  600);
  //                            Name    Title            xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  gXcanv += 50; gYcanv += 50;
  cExp3->SetFillColor(10);
  cExp3->Divide( 2,  0);
  //cExp3->Divide( 2,  0,     0.0,     0.0,   10);
  //              nx, ny, xmargin, ymargin, color
  //////////////////////////////////////////////
  cExp3->cd(1);
  TH1D *HST1;
  HST1 =HAfbDel_B_Ceex2;

  HST1->SetStats(0);
  HST1->SetTitle(0);
  HST1->GetXaxis()->SetNdivisions(8);
  HST1->DrawCopy("h");

  double ycapt =0.85, xcapt=0.5;
  CaptNDC->DrawLatex(xcapt-0.1,ycapt," Ceex2, IFI on");
  ycapt += -0.047;
  PlotSame3(HAfbDel_B_Ceex2,     xcapt, ycapt, kBlue,      25, 1, "(B) Idealized");
  PlotSame3(HAfbDel_A_Ceex2,     xcapt, ycapt, kBlack,     24, 1, "(A) Realistic");

  CaptNDC->DrawLatex(0.04,0.95,"#Delta A_{FB}(v) = A_{FB}(v,s_{+})- A_{FB}(v,s_{-})");
  CaptNDC->DrawLatex(0.60,0.02,"v_{max} ");

  ///////////////////////////////////////////////
  cExp3->cd(2);
  TH1D *HST2;
  HST2 =HAfbDel_B_Ceex21;

  HST2->SetStats(0);
  HST2->SetTitle(0);
  HST2->GetXaxis()->SetNdivisions(8);
  HST2->SetMaximum(1.1e-5); HST2->SetMinimum(-1.1e-5);
  HST2->DrawCopy("h");

  ycapt =0.85, xcapt=0.5;
  CaptNDC->DrawLatex(xcapt-0.1,ycapt,"  Ceex2-Ceex1, IFI on");
  ycapt += -0.047;
  PlotSame3(HAfbDel_B_Ceex21,   xcapt, ycapt, kBlue,     25, 1, "(B) Idealized");
  PlotSame3(HAfbDel_A_Ceex21,   xcapt, ycapt, kBlack,    25, 1, "(A) Realistic");

  CaptNDC->DrawLatex(0.04,0.95,"#delta #Delta A_{FB}(v)  ");
  CaptNDC->DrawLatex(0.60,0.02,"v_{max} ");
//
  cExp3->cd();
//
  cExp3->SaveAs("cExp3.pdf");
//
}// FigExp3



///////////////////////////////////////////////////////////////////////////////////
void FigIFI2()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigIFI2 =========================== "<<endl;
//////////////////////
  TH1D *HAfb8_vA_Ceex2_IFI  = (TH1D*)DiskFileB.Get("HAfb8_vA_Ceex2_IFI");
  TH1D *HAfb8_vB_Ceex2_IFI  = (TH1D*)DiskFileB.Get("HAfb8_vB_Ceex2_IFI");
//
  TH1D *HAfb9_vA_Ceex2_IFI  = (TH1D*)DiskFileB.Get("HAfb9_vA_Ceex2_IFI");
  TH1D *HAfb9_vB_Ceex2_IFI  = (TH1D*)DiskFileB.Get("HAfb9_vB_Ceex2_IFI");
///////
  TH1D *HAfbDel_A_Ceex2_IFI = HstDiff( "HAfbDel_A_Ceex2_IFI", HAfb9_vA_Ceex2_IFI, HAfb8_vA_Ceex2_IFI, kBlue);
  TH1D *HAfbDel_B_Ceex2_IFI = HstDiff( "HAfbDel_B_Ceex2_IFI", HAfb9_vB_Ceex2_IFI, HAfb8_vB_Ceex2_IFI, kBlue);
/////////////////////
  TH1D *HAfb8_vA_Ceex1_IFI  = (TH1D*)DiskFileB.Get("HAfb8_vA_Ceex1_IFI");
  TH1D *HAfb8_vB_Ceex1_IFI  = (TH1D*)DiskFileB.Get("HAfb8_vB_Ceex1_IFI");
//
  TH1D *HAfb9_vA_Ceex1_IFI  = (TH1D*)DiskFileB.Get("HAfb9_vA_Ceex1_IFI");
  TH1D *HAfb9_vB_Ceex1_IFI  = (TH1D*)DiskFileB.Get("HAfb9_vB_Ceex1_IFI");
///////
  TH1D *HAfbDel_A_Ceex1_IFI = HstDiff( "HAfbDel_A_Ceex1_IFI", HAfb9_vA_Ceex1_IFI, HAfb8_vA_Ceex1_IFI, kBlue);
  TH1D *HAfbDel_B_Ceex1_IFI = HstDiff( "HAfbDel_B_Ceex1_IFI", HAfb9_vB_Ceex1_IFI, HAfb8_vB_Ceex1_IFI, kBlue);
///////
  TH1D *HAfbDel_A_Ceex21_IFI = HstDiff( "HAfbDel_A_Ceex21_IFI", HAfbDel_A_Ceex2_IFI, HAfbDel_A_Ceex1_IFI, kBlue);
  TH1D *HAfbDel_B_Ceex21_IFI = HstDiff( "HAfbDel_B_Ceex21_IFI", HAfbDel_B_Ceex2_IFI, HAfbDel_B_Ceex1_IFI, kBlue);

  TLatex *CaptNDC = new TLatex(); CaptNDC->SetNDC(); // !!!
  CaptNDC->SetTextSize(0.037);
  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cIFI2 = new TCanvas("cIFI2","cIFI2", gXcanv,  gYcanv,   1200,  600);
  //                            Name    Title            xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  gXcanv += 50; gYcanv += 50;
  cIFI2->SetFillColor(10);
  cIFI2->Divide( 2,  0);
  //cIFI2->Divide( 2,  0,     0.0,     0.0,   10);
  //              nx, ny, xmargin, ymargin, color
  //////////////////////////////////////////////
  cIFI2->cd(1);
  TH1D *HST1;
  HST1 = HAfbDel_B_Ceex2_IFI;

  HST1->SetStats(0);
  HST1->SetTitle(0);
  HST1->GetXaxis()->SetNdivisions(8);
  HST1->DrawCopy("h");

  double ycapt =0.35, xcapt=0.3;
  CaptNDC->DrawLatex(xcapt-0.1,ycapt,"IFI  Ceex2,1");
  ycapt += -0.047;
  PlotSame3(HAfbDel_B_Ceex2_IFI,     xcapt, ycapt, kBlack,     26, 1, "(B) Ceex2");
  PlotSame3(HAfbDel_B_Ceex1_IFI,     xcapt, ycapt, kBlack,     27, 1, "(B) Ceex1");
  PlotSame3(HAfbDel_A_Ceex2_IFI,     xcapt, ycapt, kBlue,      24, 1, "(A) Ceex2");
  PlotSame3(HAfbDel_A_Ceex1_IFI,     xcapt, ycapt, kBlue,      25, 1, "(A) Ceex1");

  CaptNDC->DrawLatex(0.04,0.95,"#Delta A_{FB}^{IFI}(v) = A_{FB}^{IFI}(v,s_{+}) -A_{FB}^{IFI}(v,s_{-})");
  CaptNDC->DrawLatex(0.60,0.02,"v_{max} ");

  ///////////////////////////////////////////////
  cIFI2->cd(2);
  TH1D *HST2;
  HST2 = HAfbDel_A_Ceex21_IFI;

  HST2->SetStats(0);
  HST2->SetTitle(0);
  HST2->GetXaxis()->SetNdivisions(8);
  HST2->SetMaximum(1.1e-4); HST2->SetMinimum(-1.1e-4);
  HST2->DrawCopy("h");

  ycapt =0.85, xcapt=0.5;
  CaptNDC->DrawLatex(xcapt-0.1,ycapt,"IFI  Ceex2,1");
  ycapt += -0.047;
  PlotSame3(HAfbDel_B_Ceex21_IFI,   xcapt, ycapt, kBlack,    24, 1, "(B) IFI Ceex2-Ceex1");
  PlotSame3(HAfbDel_A_Ceex21_IFI,   xcapt, ycapt, kBlue,     25, 1, "(A) IFI Ceex2-Ceex1");

  CaptNDC->DrawLatex(0.04,0.95,
		  "#delta #Delta A_{FB}^{IFI}(v) = Delta A_{FB}^{IFI}(v)|_{ceex2} - Delta A_{FB}^{IFI}(v)|_{ceex2} ");
  CaptNDC->DrawLatex(0.60,0.02,"v_{max} ");
//
  cIFI2->cd();
//
  cIFI2->SaveAs("cIFI2.pdf");
//
}// FigIFI2



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
  /////////////////////////////////////////////////////////
//  LibSem.Initialize(DiskFileA95);  // for non-farm case
  /////////////////////////////////////////////////////////
  //HistNormalize();     // Renormalization of MC histograms
  ReMakeHistoMC();       // reprocessing MC histos
  KKsemMakeHisto();
  //========== PLOTTING ==========
  //
  FigExp1(); // 94GeV
  FigExp2(); // 88GeV
  //FigExp3();  // differences
  //
  //FigIFI2();
  // Template empty canvas  with 2 figures
  //FigTempl();
  //++++++++++++++++++++++++++++++++++++++++
  DiskFileA95.ls();
  DiskFileB.ls();
  DiskFileB.Write();
  DiskFileB.Close();
  //++++++++++++++++++++++++++++++++++++++++
  theApp.Run();
  //++++++++++++++++++++++++++++++++++++++++
}


