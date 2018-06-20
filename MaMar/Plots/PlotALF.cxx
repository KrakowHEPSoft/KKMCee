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
  /////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////// Cut (A) Realistic, IFI on /////
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
  //  Cut (A) calculating AFB directly, IFI on
    TH1D *hst9_Afb2_A =  HstAFB2cl("hst9_Afb2_A", hst9_Alf2CutA_Ceex2F , hst9_Alf2CutA_Ceex2);
    TH1D *hst9_Afb0_A =  HstAFB2cl("hst9_Afb0_A", hst9_Alf0CutA_Ceex2F , hst9_Alf0CutA_Ceex2);
  //
    TH1D *hst9_Afb2_B =  HstAFB2cl("hst9_Afb2_B", hst9_Alf2CutB_Ceex2F , hst9_Alf2CutB_Ceex2);
    TH1D *hst9_Afb0_B =  HstAFB2cl("hst9_Afb0_B", hst9_Alf2CutB_Ceex2F , hst9_Alf0CutB_Ceex2);
  //
    TH1D *hst9_Afb20_A = HstDiff( "hst9_Afb20_A", hst9_Afb2_A, hst9_Afb0_A, kBlue);
    TH1D *hst9_Afb20_B = HstDiff( "hst9_Afb20_B", hst9_Afb2_B, hst9_Afb0_B, kRed);
  /////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////// Cut (A) Realistic, IFI off /////
   TH1D *hst9_Al20CutA_Ceex2n    = (TH1D*)DiskFileA95.Get("hst_Al20CutA_Ceex2n");     // total CEEX2
   TH1D *hst9_Al20CutA_Ceex2nF   = (TH1D*)DiskFileA95.Get("hst_Al20CutA_Ceex2nF");    // forw CEEX2
   TH1D *hst9_Alf2CutA_Ceex2n    = (TH1D*)DiskFileA95.Get("hst_Alf2CutA_Ceex2n");     // total CEEX2
   TH1D *hst9_Alf2CutA_Ceex2nF   = (TH1D*)DiskFileA95.Get("hst_Alf2CutA_Ceex2nF");    // forw CEEX2
  ////////////// Cut (B ) Idealized, IFI off /////
   TH1D *hst9_Al20CutB_Ceex2n    = (TH1D*)DiskFileA95.Get("hst_Al20CutB_Ceex2n");     // total CEEX2
   TH1D *hst9_Al20CutB_Ceex2nF   = (TH1D*)DiskFileA95.Get("hst_Al20CutB_Ceex2nF");    // forw CEEX2
   TH1D *hst9_Alf2CutB_Ceex2n    = (TH1D*)DiskFileA95.Get("hst_Alf2CutB_Ceex2n");     // total CEEX2
   TH1D *hst9_Alf2CutB_Ceex2nF   = (TH1D*)DiskFileA95.Get("hst_Alf2CutB_Ceex2nF");    // forw CEEX2
  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
  //  Cut (A) calculating AFB from weight differences, IFI on
   TH1D *HAfb_9A_Ceex2= HstAFB4cl("HAfb_9A_Ceex2",  hst9_Al20CutA_Ceex2F, hst9_Al20CutA_Ceex2,
  	                                                hst9_Alf2CutA_Ceex2F, hst9_Alf2CutA_Ceex2);
  //
   TH1D *HAfb_9B_Ceex2= HstAFB4cl("HAfb_9B_Ceex2",  hst9_Al20CutB_Ceex2F, hst9_Al20CutB_Ceex2,
                                                    hst9_Alf2CutB_Ceex2F, hst9_Alf2CutB_Ceex2);
  // calculating AFB from weight differences, IFI off
   TH1D *HAfb_9A_Ceex2n= HstAFB4cl("HAfb_9A_Ceex2n",hst9_Al20CutA_Ceex2nF, hst9_Al20CutA_Ceex2n,
  	                                                hst9_Alf2CutA_Ceex2nF, hst9_Alf2CutA_Ceex2n);
  //
   TH1D *HAfb_9B_Ceex2n= HstAFB4cl("HAfb_9B_Ceex2n",hst9_Al20CutB_Ceex2nF, hst9_Al20CutB_Ceex2n,
                                                    hst9_Alf2CutB_Ceex2nF, hst9_Alf2CutB_Ceex2n);
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
////////////// Cut (A) Realistic, IFI on /////
   TH1D *hst8_Al20CutA_Ceex2    = (TH1D*)DiskFileA88.Get("hst_Al20CutA_Ceex2");     // total CEEX2
   TH1D *hst8_Al20CutA_Ceex2F   = (TH1D*)DiskFileA88.Get("hst_Al20CutA_Ceex2F");    // forw CEEX2
//
   TH1D *hst8_Alf2CutA_Ceex2    = (TH1D*)DiskFileA88.Get("hst_Alf2CutA_Ceex2");     // total CEEX2
   TH1D *hst8_Alf2CutA_Ceex2F   = (TH1D*)DiskFileA88.Get("hst_Alf2CutA_Ceex2F");    // forw CEEX2
////////////// Cut (B ) Idealized /////
   TH1D *hst8_Al20CutB_Ceex2    = (TH1D*)DiskFileA88.Get("hst_Al20CutB_Ceex2");     // total CEEX2
   TH1D *hst8_Al20CutB_Ceex2F   = (TH1D*)DiskFileA88.Get("hst_Al20CutB_Ceex2F");    // forw CEEX2
//
   TH1D *hst8_Alf2CutB_Ceex2    = (TH1D*)DiskFileA88.Get("hst_Alf2CutB_Ceex2");     // total CEEX2
   TH1D *hst8_Alf2CutB_Ceex2F   = (TH1D*)DiskFileA88.Get("hst_Alf2CutB_Ceex2F");    // forw CEEX2
/////////////////////////////////////////////////////////////////////////////////////////////////
////////////// Cut (A) Realistic, IFI off /////
   TH1D *hst8_Al20CutA_Ceex2n    = (TH1D*)DiskFileA88.Get("hst_Al20CutA_Ceex2n");     // total CEEX2
   TH1D *hst8_Al20CutA_Ceex2nF   = (TH1D*)DiskFileA88.Get("hst_Al20CutA_Ceex2nF");    // forw CEEX2
   TH1D *hst8_Alf2CutA_Ceex2n    = (TH1D*)DiskFileA88.Get("hst_Alf2CutA_Ceex2n");     // total CEEX2
   TH1D *hst8_Alf2CutA_Ceex2nF   = (TH1D*)DiskFileA88.Get("hst_Alf2CutA_Ceex2nF");    // forw CEEX2
////////////// Cut (B ) Idealized, IFI off /////
   TH1D *hst8_Al20CutB_Ceex2n    = (TH1D*)DiskFileA88.Get("hst_Al20CutB_Ceex2n");     // total CEEX2
   TH1D *hst8_Al20CutB_Ceex2nF   = (TH1D*)DiskFileA88.Get("hst_Al20CutB_Ceex2nF");    // forw CEEX2
   TH1D *hst8_Alf2CutB_Ceex2n    = (TH1D*)DiskFileA88.Get("hst_Alf2CutB_Ceex2n");     // total CEEX2
   TH1D *hst8_Alf2CutB_Ceex2nF   = (TH1D*)DiskFileA88.Get("hst_Alf2CutB_Ceex2nF");    // forw CEEX2
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//  Cut (A) calculating AFB from weight differences, IFI on
   TH1D *HAfb_8A_Ceex2= HstAFB4cl("HAfb_8A_Ceex2",  hst8_Al20CutA_Ceex2F, hst8_Al20CutA_Ceex2,
   	                                                hst8_Alf2CutA_Ceex2F, hst8_Alf2CutA_Ceex2);
   //
   TH1D *HAfb_8B_Ceex2= HstAFB4cl("HAfb_8B_Ceex2",  hst8_Al20CutB_Ceex2F, hst8_Al20CutB_Ceex2,
                                                    hst8_Alf2CutB_Ceex2F, hst8_Alf2CutB_Ceex2);
// calculating AFB from weight differences, IFI off
   TH1D *HAfb_8A_Ceex2n= HstAFB4cl("HAfb_8A_Ceex2n",hst8_Al20CutA_Ceex2nF, hst8_Al20CutA_Ceex2n,
   	                                                hst8_Alf2CutA_Ceex2nF, hst8_Alf2CutA_Ceex2n);
//
   TH1D *HAfb_8B_Ceex2n= HstAFB4cl("HAfb_8B_Ceex2n",hst8_Al20CutB_Ceex2nF, hst8_Al20CutB_Ceex2n,
                                                    hst8_Alf2CutB_Ceex2nF, hst8_Alf2CutB_Ceex2n);
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
// CEEX2-CEEX1
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
////////////////////////////////////////////////////////////////////////////////////////////////////
// AFB calculated directly
  TH1D *HAfb9_vB_Ceex2 = HstAFB("HAfb9_vB_Ceex2", hst9_vB_Ceex2_F, hst9_vB_Ceex2);
  TH1D *HAfb9_vB_Ceex1 = HstAFB("HAfb9_vB_Ceex1", hst9_vB_Ceex1_F, hst9_vB_Ceex1); // ???
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
//  CEEX2-CEEX1
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
////////////////////////////////////////////////////////////////////////////////////////////////////
// AFB calculated directly
  TH1D *HAfb8_vB_Ceex2 = HstAFB("HAfb8_vB_Ceex2", hst8_vB_Ceex2_F, hst8_vB_Ceex2);
  TH1D *HAfb8_vB_Ceex1 = HstAFB("HAfb8_vB_Ceex1", hst8_vB_Ceex1_F, hst8_vB_Ceex1); // ???
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
  CaptNDC->DrawLatex(0.60,0.02,"#delta#alpha/#alpha ");
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
  CaptNDC->DrawLatex(0.60,0.02,"#delta#alpha/#alpha ");

  //
  cAlf1->cd();
  //
  cAlf1->SaveAs("cAlf1.pdf");
//
}// FigAlf1

///////////////////////////////////////////////////////////////////////////////////
void FigAlf2()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigAlf2 =========================== "<<endl;

  TH1D *HAfb_9A_Ceex2 = (TH1D*)DiskFileB.Get("HAfb_9A_Ceex2");
  TH1D *HAfb_9B_Ceex2 = (TH1D*)DiskFileB.Get("HAfb_9B_Ceex2");
  TH1D *HAfb_9A_Ceex2n= (TH1D*)DiskFileB.Get("HAfb_9A_Ceex2n");
  TH1D *HAfb_9B_Ceex2n= (TH1D*)DiskFileB.Get("HAfb_9A_Ceex2n");
  TH1D *hst9_DelAfb   = (TH1D*)DiskFileB.Get("hst9_DelAfb");

  TH1D *HAfb_9B_Ceex2ifi = HstDiff( "HAfb_9B_Ceex2ifi", HAfb_9B_Ceex2, HAfb_9B_Ceex2n, kBlue);
  HAfb_9B_Ceex2ifi->Scale(10.0);

  TLatex *CaptNDC = new TLatex(); CaptNDC->SetNDC(); // !!!
  CaptNDC->SetTextSize(0.037);
  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cAlf2 = new TCanvas("cAlf2","cAlf2", gXcanv,  gYcanv,   1200,  600);
  //                            Name    Title            xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  gXcanv += 50; gYcanv += 50;
  cAlf2->SetFillColor(10);
  cAlf2->Divide( 2,  0);
  //cAlf2->Divide( 2,  0,     0.0,     0.0,   10);
  //              nx, ny, xmargin, ymargin, color
  //////////////////////////////////////////////
  cAlf2->cd(1);
  TH1D *HST1 = HAfb_9A_Ceex2;
  HST1 =hst9_DelAfb;

  HST1->SetStats(0);
  HST1->SetTitle(0);
  HST1->GetXaxis()->SetNdivisions(8);
  HST1->DrawCopy("h");

  double ycapt =0.80, xcapt=0.3;
  CaptNDC->DrawLatex(xcapt-0.1,ycapt," #sqrt{s} = 94.3GeV, (A) Realistic");
  ycapt += -0.047;
  PlotSame3(hst9_DelAfb,     xcapt, ycapt, kPine,     0, 3, "Born");
  PlotSame3(HAfb_9A_Ceex2,   xcapt, ycapt, kMagenta, 24, 1, "IFI on");
  PlotSame3(HAfb_9A_Ceex2n,  xcapt, ycapt, kBlue,    25, 1, "IFI off");
  PlotSame3(HAfb_9B_Ceex2ifi,xcapt, ycapt, kRed,      0, 5, "IFIx10");

  CaptNDC->DrawLatex(0.04,0.95,"#delta A_{FB}(#alpha)  ");
  CaptNDC->DrawLatex(0.60,0.02,"#delta#alpha/#alpha ");

  ///////////////////////////////////////////////
  cAlf2->cd(2);
  TH1D *HST2 = HAfb_9B_Ceex2;
  HST2 =hst9_DelAfb;

  HST2->SetStats(0);
  HST2->SetTitle(0);
  HST2->GetXaxis()->SetNdivisions(8);
  HST2->DrawCopy("h");

  ycapt =0.80, xcapt=0.3;
  CaptNDC->DrawLatex(xcapt-0.1,ycapt," #sqrt{s} = 94.3GeV, (B) Idealized,");
  ycapt += -0.047;
  PlotSame3(hst9_DelAfb,   xcapt, ycapt, kPine,     0, 3, "Born");
  PlotSame3(HAfb_9B_Ceex2, xcapt, ycapt, kMagenta, 24, 1, "IFI on");
  PlotSame3(HAfb_9B_Ceex2n,xcapt, ycapt, kBlue,    25, 1, "IFI off");

  CaptNDC->DrawLatex(0.04,0.95,"#delta A_{FB}(#alpha)  ");
  CaptNDC->DrawLatex(0.60,0.02,"#delta#alpha/#alpha ");

  //
  cAlf2->cd();
  //
  cAlf2->SaveAs("cAlf2.pdf");
//
}// FigAlf2

///////////////////////////////////////////////////////////////////////////////////
void FigAlf3()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigAlf3 =========================== "<<endl;

  TH1D *HAfb_9A_Ceex2 = (TH1D*)DiskFileB.Get("HAfb_9A_Ceex2");
  TH1D *HAfb_9B_Ceex2 = (TH1D*)DiskFileB.Get("HAfb_9B_Ceex2");
  TH1D *HAfb_9A_Ceex2n= (TH1D*)DiskFileB.Get("HAfb_9A_Ceex2n");
  TH1D *HAfb_9B_Ceex2n= (TH1D*)DiskFileB.Get("HAfb_9A_Ceex2n");

  TH1D *HAfb_9A_Ceex2sym = (TH1D*)HAfb_9A_Ceex2->Clone("HAfb_9A_Ceex2sym");
  HAfb_9A_Ceex2sym->Reset();
  TH1D *HAfb_9A_Ceex2nsym = (TH1D*)HAfb_9A_Ceex2n->Clone("HAfb_9A_Ceex2nsym");
  HAfb_9A_Ceex2nsym->Reset();
  double bin,err;
  int Nb = HAfb_9A_Ceex2sym->GetXaxis()->GetNbins();
  for(int ix=1; ix <= Nb; ix++){
	  bin =     HAfb_9A_Ceex2->GetBinContent(ix)       +HAfb_9A_Ceex2->GetBinContent(Nb+1-ix);
	  err = sqr(HAfb_9A_Ceex2->GetBinError(  ix))  +sqr(HAfb_9A_Ceex2->GetBinError(  Nb+1-ix));
	  HAfb_9A_Ceex2sym->SetBinContent(ix,bin);
	  HAfb_9A_Ceex2sym->SetBinError(  ix,sqrt(err));
	  bin =     HAfb_9A_Ceex2n->GetBinContent(ix)       +HAfb_9A_Ceex2n->GetBinContent(Nb+1-ix);
	  err = sqr(HAfb_9A_Ceex2n->GetBinError(  ix))  +sqr(HAfb_9A_Ceex2n->GetBinError(  Nb+1-ix));
	  HAfb_9A_Ceex2nsym->SetBinContent(ix,bin);
	  HAfb_9A_Ceex2nsym->SetBinError(  ix,sqrt(err));
  }//ix

  TLatex *CaptNDC = new TLatex(); CaptNDC->SetNDC(); // !!!
  CaptNDC->SetTextSize(0.037);
  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cAlf3 = new TCanvas("cAlf3","cAlf3", gXcanv,  gYcanv,   600,  600);
  //                            Name    Title            xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  gXcanv += 50; gYcanv += 50;
  cAlf3->SetFillColor(10);
  //cAlf3->Divide( 2,  0);
  //              nx, ny, xmargin, ymargin, color
  //////////////////////////////////////////////
  cAlf3->cd(1);
  TH1D *HST1 = HAfb_9A_Ceex2sym;

  HST1->SetStats(0);
  HST1->SetTitle(0);
  HST1->GetXaxis()->SetNdivisions(8);
  HST1->DrawCopy("h");

  double ycapt =0.91, xcapt=0.45;
  CaptNDC->DrawLatex(xcapt-0.1,ycapt," #sqrt{s} = 94.3GeV, Test of linearity");
  ycapt += -0.047;
  PlotSame3(HAfb_9A_Ceex2sym,   xcapt, ycapt, kBlack, 24, 1, "(A) Realistic");
  PlotSame3(HAfb_9A_Ceex2nsym,  xcapt, ycapt, kBlue,  25, 2, "(B) Idealized");

  CaptNDC->DrawLatex(0.04,0.95,"#delta A_{FB}(#alpha)  ");
  CaptNDC->DrawLatex(0.60,0.02,"#delta#alpha/#alpha ");
  //
  cAlf3->cd();
  //
  cAlf3->SaveAs("cAlf3.pdf");
//
}// FigAlf3



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

  HST1->SetStats(0);
  HST1->SetTitle(0);
  HST1->GetXaxis()->SetNdivisions(8);
  HST1->DrawCopy("h");

  double ycapt =0.91, xcapt=0.3;
  CaptNDC->DrawLatex(xcapt-0.1,ycapt," #sqrt{s} = 94.3GeV, Ceex2, IFI on");
  ycapt += -0.047;
  PlotSame3(HAfb9_vB_Ceex2,     xcapt, ycapt, kBlue,      25, 1, "(B) Idealized");
  PlotSame3(HAfb9_vA_Ceex2,     xcapt, ycapt, kBlack,     24, 1, "(A) Realistic");

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

  TH1D *HAfb8_vB_Ceex2  = (TH1D*)DiskFileB.Get("HAfb8_vB_Ceex2");
  TH1D *HAfb8_vB_Ceex1  = (TH1D*)DiskFileB.Get("HAfb8_vB_Ceex1");
  TH1D *HAfb8_vB_Ceex21 = (TH1D*)DiskFileB.Get("HAfb8_vB_Ceex21");

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
  HST1->DrawCopy("h");

  double ycapt =0.91, xcapt=0.3;
  CaptNDC->DrawLatex(xcapt-0.1,ycapt," #sqrt{s} = 87.90GeV, Ceex2, IFI on");
  ycapt += -0.047;
  PlotSame3(HAfb8_vB_Ceex2,     xcapt, ycapt, kBlue,      25, 1, "(B) Idealized");
  PlotSame3(HAfb8_vA_Ceex2,     xcapt, ycapt, kBlack,     24, 1, "(A) Realistic");

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
  //HistNormalize();     // Renormalization of MC histograms
  ReMakeHistoMC();       // reprocessing MC histos
  //========== PLOTTING ==========
  //
  FigAlf1();
  FigAlf2();
  FigAlf3();
  //
  FigExp1();
  FigExp2();
  FigExp3();
  //
  FigIFI2();
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


