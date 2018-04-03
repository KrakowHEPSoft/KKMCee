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

TFile DiskFileA88("../workKKMC/histo.root");  // apr.2018
//TFile DiskFileA95("../workKKMC/histo.root");  // apr.2018

//TFile DiskFileA88("../workKKMC/histo.root_88GeV.new");  // jan.2018
TFile DiskFileA95("../workKKMC/histo.root_95GeV.new");  // jan.2018

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
//int  gNBmax =0;    // for |cos(theta)|<1.00

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


////////////////////////////////////////////////////////////////////////////////
double BWR(double E, double M, double G)
{
	double s = E*E;
	return sqr(s)/(sqr(s-M*M)+ sqr(s*G/M)); // runing
//	return sqr(s)/(sqr(s-M*M)+ sqr(G*M));   // non-runing
}

////////////////////////////////////////////////////////////////////////////////
double AFBc(double E1, double MZ, double GamZ, double SinW2, double AlfRun)
{
double GFermi = 1.166397e-5;  // one digit more than in old KKdefaults!
double Pi =3.1415926535897932;
double ga,gv,cG,cZ;
ga = -1.0/2.0;
gv = ga*(1-4*SinW2);
cG = AlfRun;
cZ = MZ*MZ*GFermi/(2*sqrt(2.0)*Pi);
double afb1,x1,y1;
x1  =   cG*cG                                             // GG*(1+c*c)
    + 2*cG*cZ *sqr(gv) *(1-sqr(MZ/E1))*BWR(E1,MZ,GamZ)   // GZ*(1+c*c)
    +   cZ*cZ *sqr(ga*ga+gv*gv)       *BWR(E1,MZ,GamZ);  // ZZ*(1+c*c)
y1  = 2*cG*cZ *sqr(ga) *(1-sqr(MZ/E1))*BWR(E1,MZ,GamZ)   // GZ*2*c
    + cZ*cZ*4* ga*ga*gv*gv            *BWR(E1,MZ,GamZ);  // ZZ*2*c
afb1 = 0.75*y1/x1;
return afb1;
}



///////////////////////////////////////////////////////////////////////////////////
void ReMakeHistoMC()
{
//------------------------------------------------------------------------
  cout<<" ========================= ReMakeHistoMC =========================== "<<endl;
  ////////////// Cut (A) Realistic /////
    TH1D *hst9_Al20CutA_Ceex2    = (TH1D*)DiskFileA95.Get("hst_Al20CutA_Ceex2");     // total CEEX2
    TH1D *hst9_Al20CutA_Ceex2F   = (TH1D*)DiskFileA95.Get("hst_Al20CutA_Ceex2F");    // forw CEEX2
  //
    TH1D *hst9_Alf0CutA_Ceex2    = (TH1D*)DiskFileA95.Get("hst_Alf0CutA_Ceex2");     // total CEEX2
    TH1D *hst9_Alf0CutA_Ceex2F   = (TH1D*)DiskFileA95.Get("hst_Alf0CutA_Ceex2F");    // forw CEEX2
  //
    TH1D *hst9_Alf2CutA_Ceex2    = (TH1D*)DiskFileA95.Get("hst_Alf2CutA_Ceex2");     // total CEEX2
    TH1D *hst9_Alf2CutA_Ceex2F   = (TH1D*)DiskFileA95.Get("hst_Alf2CutA_Ceex2F");    // forw CEEX2
  ////////////// Cut (B ) Idealized /////
    TH1D *hst9_Al20CutB_Ceex2    = (TH1D*)DiskFileA95.Get("hst_Al20CutB_Ceex2");     // total CEEX2
    TH1D *hst9_Al20CutB_Ceex2F   = (TH1D*)DiskFileA95.Get("hst_Al20CutB_Ceex2F");    // forw CEEX2
  //
    TH1D *hst9_Alf0CutB_Ceex2    = (TH1D*)DiskFileA95.Get("hst_Alf0CutB_Ceex2");     // total CEEX2
    TH1D *hst9_Alf0CutB_Ceex2F   = (TH1D*)DiskFileA95.Get("hst_Alf0CutB_Ceex2F");    // forw CEEX2
  //
    TH1D *hst9_Alf2CutB_Ceex2    = (TH1D*)DiskFileA95.Get("hst_Alf2CutB_Ceex2");     // total CEEX2
    TH1D *hst9_Alf2CutB_Ceex2F   = (TH1D*)DiskFileA95.Get("hst_Alf2CutB_Ceex2F");    // forw CEEX2
  /////////////////////////////////////////////////////////////////////////////////////////////////
  // calculating AFB directly
    TH1D *hst9_Afb2_A =  HstAFB2cl("hst9_Afb2_A", hst9_Alf2CutA_Ceex2F , hst9_Alf2CutA_Ceex2);
    TH1D *hst9_Afb0_A =  HstAFB2cl("hst9_Afb0_A", hst9_Alf0CutA_Ceex2F , hst9_Alf0CutA_Ceex2);
  //
    TH1D *hst9_Afb2_B =  HstAFB2cl("hst9_Afb2_B", hst9_Alf2CutB_Ceex2F , hst9_Alf2CutB_Ceex2);
    TH1D *hst9_Afb0_B =  HstAFB2cl("hst9_Afb0_B", hst9_Alf2CutB_Ceex2F , hst9_Alf0CutB_Ceex2);
  //
    TH1D *hst9_Afb20_A = HstDiff( "hst9_Afb20_A", hst9_Afb2_A, hst9_Afb0_A, kBlue);
    TH1D *hst9_Afb20_B = HstDiff( "hst9_Afb20_B", hst9_Afb2_B, hst9_Afb0_B, kRed);
  /////////////////////////////////////////////////////////////////////////////////////////////////
  // calculating AFB from weight differences
    TH1D *HAfb_9A_Ceex2= HstAFB4cl("HAfb_9A_Ceex2", hst9_Al20CutA_Ceex2F, hst9_Al20CutA_Ceex2,
  		                                          hst9_Alf2CutA_Ceex2F, hst9_Alf2CutA_Ceex2);
  //
    TH1D *HAfb_9B_Ceex2= HstAFB4cl("HAfb_9B_Ceex2", hst9_Al20CutB_Ceex2F, hst9_Al20CutB_Ceex2,
  		                                          hst9_Alf2CutB_Ceex2F, hst9_Alf2CutB_Ceex2);
  /////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////// Cut (A) Realistic /////
   TH1D *hst9_Al20CutA_Ceex2n    = (TH1D*)DiskFileA95.Get("hst_Al20CutA_Ceex2n");     // total CEEX2
   TH1D *hst9_Al20CutA_Ceex2nF   = (TH1D*)DiskFileA95.Get("hst_Al20CutA_Ceex2nF");    // forw CEEX2
   TH1D *hst9_Alf2CutA_Ceex2n    = (TH1D*)DiskFileA95.Get("hst_Alf2CutA_Ceex2n");     // total CEEX2
   TH1D *hst9_Alf2CutA_Ceex2nF   = (TH1D*)DiskFileA95.Get("hst_Alf2CutA_Ceex2nF");    // forw CEEX2
  ////////////// Cut (B ) Idealized /////
   TH1D *hst9_Al20CutB_Ceex2n    = (TH1D*)DiskFileA95.Get("hst_Al20CutB_Ceex2n");     // total CEEX2
   TH1D *hst9_Al20CutB_Ceex2nF   = (TH1D*)DiskFileA95.Get("hst_Al20CutB_Ceex2nF");    // forw CEEX2
   TH1D *hst9_Alf2CutB_Ceex2n    = (TH1D*)DiskFileA95.Get("hst_Alf2CutB_Ceex2n");     // total CEEX2
   TH1D *hst9_Alf2CutB_Ceex2nF   = (TH1D*)DiskFileA95.Get("hst_Alf2CutB_Ceex2nF");    // forw CEEX2
  /////////////////////////////////////////////////////////////////////////////////////////////////
  // calculating AFB from weight differences
   TH1D *HAfb_9A_Ceex2n= HstAFB4cl("HAfb_9A_Ceex2n",hst9_Al20CutA_Ceex2nF, hst9_Al20CutA_Ceex2n,
  	                                              hst9_Alf2CutA_Ceex2nF, hst9_Alf2CutA_Ceex2n);
  //
   TH1D *HAfb_9B_Ceex2n= HstAFB4cl("HAfb_9B_Ceex2n",hst9_Al20CutB_Ceex2nF, hst9_Al20CutB_Ceex2n,
                                                    hst9_Alf2CutB_Ceex2nF, hst9_Alf2CutB_Ceex2n);
  /////////////////////////////////////////////////////////////////////////////////////////////////
    double E1 =  87.90e0, E2 =  94.30e0;
  //
    int    nPt  = 8;
    double MZ = 91.187e0;
    double alfinv0 = 137.0359895e0;
    double GammZ = 2.50072032;    // from KKdefaults!
    double SinW2 =  0.23152;      // m_xpar[503]
  //
    double Alf0MZ =0.00775874495;
    Alf0MZ=1/128.886825;  // the same
  /////////////////////////////////////////////////////////////////
    TH1D *hst8_DelAfb = (TH1D*)hst9_Afb0_B->Clone("hst8_DelAfb");
    hst8_DelAfb->Reset();
    TH1D *hst9_DelAfb = (TH1D*)hst9_Afb0_B->Clone("hst9_DelAfb");
    hst9_DelAfb->Reset();
    Double_t Xmax = hst8_DelAfb->GetXaxis()->GetXmax();
    Double_t Xmin = hst8_DelAfb->GetXaxis()->GetXmin();
    double DelAlf = (Xmax-Xmin)/2; // relative variation of alphaQED
  //
    double beta = 20.0/(9.0*3.141593); // RGE beta for alpha_QED for 3 leptons and 5 quarks
    double Y1  = beta* 2*log(MZ/E1),  Y2  = beta* 2*log(E2/MZ);  // linear interpolation
  // Born from analytical formula
    double DelAfb1, DelAfb2, AlfRun1,AlfRun2;
    double AlfRunE1, AlfRunE2;
    for(int ix=1; ix <= nPt; ix++){
      AlfRunE1 = Alf0MZ*(1-Y1*Alf0MZ );
      AlfRunE2 = Alf0MZ*(1+Y2*Alf0MZ );
  //
      AlfRun1 = AlfRunE1*(1 -DelAlf +(ix-0.5)/nPt* 2*DelAlf);
      AlfRun2 = AlfRunE2*(1 -DelAlf +(ix-0.5)/nPt* 2*DelAlf);
      DelAfb1 = AFBc(E1, MZ, GammZ, SinW2, AlfRun1) - AFBc(E1, MZ, GammZ, SinW2, AlfRunE1);
      DelAfb2 = AFBc(E2, MZ, GammZ, SinW2, AlfRun2) - AFBc(E2, MZ, GammZ, SinW2, AlfRunE2) ;
   //
       hst8_DelAfb->SetBinContent(  ix, DelAfb1 );
       hst8_DelAfb->SetBinError(    ix, 0.0 );
       hst9_DelAfb->SetBinContent(  ix, DelAfb2 );
       hst9_DelAfb->SetBinError(    ix, 0.0 );
    }//ix
   /////////////////////////////////////////////////////////////////////////////////////////////////

  cout<<" ===================== End of ReMakeHistoMC ======================== "<<endl;

}// ReMakeHistoMC

///////////////////////////////////////////////////////////////////////////////////
void FigAlf1()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigAlf1 =========================== "<<endl;

  TH1D *hst9_DelAfb   = (TH1D*)DiskFileB.Get("hst9_DelAfb");
  TH1D *HAfb_9A_Ceex2 = (TH1D*)DiskFileB.Get("HAfb_9A_Ceex2");
  TH1D *HAfb_9B_Ceex2 = (TH1D*)DiskFileB.Get("HAfb_9B_Ceex2");
  TH1D *HAfb_9A_Ceex2n= (TH1D*)DiskFileB.Get("HAfb_9A_Ceex2n");
  TH1D *HAfb_9B_Ceex2n= (TH1D*)DiskFileB.Get("HAfb_9A_Ceex2n");

  TLatex *CaptNDC = new TLatex(); CaptNDC->SetNDC(); // !!!
  CaptNDC->SetTextSize(0.037);

  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cAlf1 = new TCanvas("cAlf1","cAlf1", gXcanv,  gYcanv,   1200,  600);
  //                            Name    Title            xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  gXcanv += 50; gYcanv += 50;
  cAlf1->SetFillColor(10);
  cAlf1->Divide( 2,  0);
  //cAlf1->Divide( 2,  0,     0.0,     0.0,   10);
  //              nx, ny, xmargin, ymargin, color
  //////////////////////////////////////////////
  cAlf1->cd(1);
  TH1D *HST1;
  HST1 = hst9_DelAfb;

  HST1->SetStats(0);
  HST1->SetTitle(0);
  HST1->GetXaxis()->SetNdivisions(4);
  HST1->SetLineStyle(3);
  HST1->DrawCopy("h");

  // https://root.cern.ch/doc/master/classTAttLine.html
  // https://root.cern.ch/doc/master/classTAttMarker.html
  double ycapt =0.80, xcapt=0.3;
  CaptNDC->DrawLatex(xcapt,ycapt," #sqrt{s} = 94.3GeV");
  ycapt += -0.047;
  PlotSame3(HAfb_9A_Ceex2, xcapt, ycapt, kMagenta, 24, 1, "(A) Realistic, IFIon");
  PlotSame3(HAfb_9B_Ceex2, xcapt, ycapt, kBlue,    25, 1, "(B) Idealized, IFIon");
  PlotSame3(hst9_DelAfb,   xcapt, ycapt, kPine,     0, 3, "Born");

  CaptNDC->DrawLatex(0.04,0.95,"#delta A_{FB}(#alpha)  ");
  CaptNDC->DrawLatex(0.60,0.02,"#delta #alpha / #alpha ");
  ///////////////////////////////////////////////
  cAlf1->cd(2);
  TH1D *HST2 = hst9_DelAfb;

  HST2->SetStats(0);
  HST2->SetTitle(0);
  HST2->GetXaxis()->SetNdivisions(4);
  HST2->DrawCopy("h");

  ycapt =0.80, xcapt=0.3;
  CaptNDC->DrawLatex(xcapt,ycapt," #sqrt{s} = 94.3GeV");
  ycapt += -0.047;
  PlotSame3(HAfb_9A_Ceex2n, xcapt, ycapt, kMagenta, 24, 1, "(A) Realistic, IFIoff");
  PlotSame3(HAfb_9B_Ceex2n, xcapt, ycapt, kBlue,    25, 1, "(B) Idealized, IFIoff");
  PlotSame3(hst9_DelAfb,    xcapt, ycapt, kPine,     0, 3, "Born");

  CaptNDC->DrawLatex(0.12,0.95,"A_{FB}(#alpha) ");
  CaptNDC->DrawLatex(0.60,0.02,"#delta #alpha / #alpha ");

  //
  cAlf1->cd();
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
  //HistNormalize();     // Renormalization of MC histograms
  ReMakeHistoMC();       // reprocessing MC histos
  //========== PLOTTING ==========
  //
  FigAlf1();
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


