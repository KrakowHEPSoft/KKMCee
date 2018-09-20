//    make PlotALF2-run
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

//TFile DiskFileA88("../workKKMC/histo.root");  // current
//TFile DiskFileA95("../workKKMC/histo.root");  // current

TFile DiskFileA95("../workKKMC/histo.root_95GeV.sept_200M");  // sep.2018
TFile DiskFileA88("../workKKMC/histo.root_88GeV.sept_290M");  // sep.2018
TFile DiskFileA91("../workKKMC/histo.root_91GeV.sept_180M");  // sep.2018

//TFile DiskFileA88("../workKKMC/histo.root_88GeV.new");  // recent
//TFile DiskFileA95("../workKKMC/histo.root_95GeV.new");  // recent

//TFile DiskFileA88("../workKKMC/histo.root_88GeV_apr_13G");  // apr.2018
//TFile DiskFileA95("../workKKMC/histo.root_95GeV_apr_38G");  // apr.2018

////  ****************** FOAM ******************

//TFile DiskFileF91("../workFoam3/histo.root"); // current
//TFile DiskFileF95("../workFoam3/histo.root"); // current
//TFile DiskFileF88("../workFoam3/histo.root"); // current

TFile DiskFileF95("../workFoam3/histo.root_95GeV_gps_58M");  // gps
TFile DiskFileF88("../workFoam3/histo.root_88GeV_gps_70M");  // gps
TFile DiskFileF91("../workFoam3/histo.root_91GeV_gps_94M");  // gps


//TFile DiskFileF91("../workFoam3/histo.root_91GeV_kksem_180M"); // kksem
//TFile DiskFileF95("../workFoam3/histo.root_95GeV_kksem_80M");  // kksem
//TFile DiskFileF88("../workFoam3/histo.root_88GeV_kksem_110M"); // current


//TFile DiskFileF95("../workFoam3/histo.root_95GeV_new2"); // current
//TFile DiskFileF88("../workFoam3/histo.root_88GeV_new"); // current
//TFile DiskFileF91("../workFoam3/histo.root_91GeV_new"); // current

//  Sept. 2018
//TFile DiskFileF95("../workFoam3/histo.root_95GeV_2G"); //


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
void PlotSame2(TH1D *HST, double &ycapt, Int_t kolor, double xx,  TString label,  TString opis)
{
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.035);
  HST->SetLineColor(kolor);
  HST->DrawCopy("hsame");      // Magenta
  CaptT->SetTextColor(kolor);
  ycapt += -0.047;
  double xcapt = 0.40;
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
void HisReMakeFOAMi()
{
//
  cout<<"==================================================================="<<endl;
  cout<<"================ HisReMakeFOAMi  BEGIN   ============================"<<endl;
//////////////////////////////////////////////////////////////////
  cout<<"  Renormalizing  and reprocessing histograms from FOAM"<<endl;

  TH1D *HST_FOAM_NORMA3i9 = (TH1D*)DiskFileF95.Get("HST_FOAM_NORMA3i");

  TH2D *SCT_xc_EEX2i9   = (TH2D*)DiskFileF95.Get("SCT_xc_EEX2i");    // FOAM small range x<0.20
  TH2D *SCT_xc_EEX2n9   = (TH2D*)DiskFileF95.Get("SCT_xc_EEX2n");    // FOAM small range x<0.20

  TH2D *SCN_xc_EEX2i9   = (TH2D*)DiskFileF95.Get("SCN_xc_EEX2i");    // FOAM small range x<0.02
  TH2D *SCN_xc_EEX2n9   = (TH2D*)DiskFileF95.Get("SCN_xc_EEX2n");    // FOAM small range x<0.02

  TH1D *HST_FOAM_NORMA3i8 = (TH1D*)DiskFileF88.Get("HST_FOAM_NORMA3i");

  TH2D *SCT_xc_EEX2i8   = (TH2D*)DiskFileF88.Get("SCT_xc_EEX2i");    // FOAM small range x<0.20
  TH2D *SCT_xc_EEX2n8   = (TH2D*)DiskFileF88.Get("SCT_xc_EEX2n");    // FOAM small range x<0.20

  TH2D *SCN_xc_EEX2i8   = (TH2D*)DiskFileF88.Get("SCN_xc_EEX2i");    // FOAM small range x<0.02
  TH2D *SCN_xc_EEX2n8   = (TH2D*)DiskFileF88.Get("SCN_xc_EEX2n");    // FOAM small range x<0.02

  TH1D *HST_FOAM_NORMA3iZ = (TH1D*)DiskFileF91.Get("HST_FOAM_NORMA3i");

  TH2D *SCT_xc_EEX2iZ   = (TH2D*)DiskFileF91.Get("SCT_xc_EEX2i");    // FOAM small range x<0.20
  TH2D *SCT_xc_EEX2nZ   = (TH2D*)DiskFileF91.Get("SCT_xc_EEX2n");    // FOAM small range x<0.20

  TH2D *SCN_xc_EEX2iZ   = (TH2D*)DiskFileF91.Get("SCN_xc_EEX2i");    // FOAM small range x<0.02
  TH2D *SCN_xc_EEX2nZ   = (TH2D*)DiskFileF91.Get("SCN_xc_EEX2n");    // FOAM small range x<0.02

/////////////////////////////////////////////////////////////////////////
///////////// For AFB renormalizing not needed
// FOAM3i
//   HisNorm2(HST_FOAM_NORMA3i9, SCT_xc_EEX2i9 );    // normalizing
//   HisNorm2(HST_FOAM_NORMA3i9, SCT_xc_EEX2n9 );    // normalizing

//   HisNorm2(HST_FOAM_NORMA3i9, SCN_xc_EEX2i9 );    // normalizing
//   HisNorm2(HST_FOAM_NORMA3i9, SCN_xc_EEX2n9 );    // normalizing

////////////////////////////////////////////////////////////////////////
// sigma(vmax) and AFB(vmax) from KKfoam scat.
// vmax<0.2 and/or <0.02, 100 bins in ctheta, 100 bins in vv
   int NbMax=45;
//   NbMax=0;
////////////////------------------95GeV-------------------------------////////
   TH1D  *Htot2_xmax_EEX2i9 = HstProjV("Htot2_xmax_EEX2i9",SCT_xc_EEX2i9,NbMax); //  x<0.20
   TH1D  *Hafb2_xmax_EEX2i9 = HstProjA("Hafb2_xmax_EEX2i9",SCT_xc_EEX2i9,NbMax); //  x<0.20
//
   TH1D  *Htot2_xmax_EEX2n9 = HstProjV("Htot2_xmax_EEX2n9",SCT_xc_EEX2n9,NbMax); //  x<0.20
   TH1D  *Hafb2_xmax_EEX2n9 = HstProjA("Hafb2_xmax_EEX2n9",SCT_xc_EEX2n9,NbMax); //  x<0.20
//
   TH1D  *Htot1_xmax_EEX2i9 = HstProjV("Htot1_xmax_EEX2i9",SCN_xc_EEX2i9,NbMax); //  x<0.02
   TH1D  *Hafb1_xmax_EEX2i9 = HstProjA("Hafb1_xmax_EEX2i9",SCN_xc_EEX2i9,NbMax); //  x<0.02
//
   TH1D  *Htot1_xmax_EEX2n9 = HstProjV("Htot1_xmax_EEX2n9",SCN_xc_EEX2n9,NbMax); //  x<0.02
   TH1D  *Hafb1_xmax_EEX2n9 = HstProjA("Hafb1_xmax_EEX2n9",SCN_xc_EEX2n9,NbMax); //  x<0.02
///////////////---------------- 88GeV-------------------------------////////
   TH1D  *Htot2_xmax_EEX2i8 = HstProjV("Htot2_xmax_EEX2i8",SCT_xc_EEX2i8,NbMax); //  x<0.20
   TH1D  *Hafb2_xmax_EEX2i8 = HstProjA("Hafb2_xmax_EEX2i8",SCT_xc_EEX2i8,NbMax); //  x<0.20
//
   TH1D  *Htot2_xmax_EEX2n8 = HstProjV("Htot2_xmax_EEX2n8",SCT_xc_EEX2n8,NbMax); //  x<0.20
   TH1D  *Hafb2_xmax_EEX2n8 = HstProjA("Hafb2_xmax_EEX2n8",SCT_xc_EEX2n8,NbMax); //  x<0.20
//
   TH1D  *Htot1_xmax_EEX2i8 = HstProjV("Htot1_xmax_EEX2i8",SCN_xc_EEX2i8,NbMax); //  x<0.02
   TH1D  *Hafb1_xmax_EEX2i8 = HstProjA("Hafb1_xmax_EEX2i8",SCN_xc_EEX2i8,NbMax); //  x<0.02
//
   TH1D  *Htot1_xmax_EEX2n8 = HstProjV("Htot1_xmax_EEX2n8",SCN_xc_EEX2n8,NbMax); //  x<0.02
   TH1D  *Hafb1_xmax_EEX2n8 = HstProjA("Hafb1_xmax_EEX2n8",SCN_xc_EEX2n8,NbMax); //  x<0.02
   ///////////////---------------- 91GeV-------------------------------////////
   TH1D  *Htot2_xmax_EEX2iZ = HstProjV("Htot2_xmax_EEX2iZ",SCT_xc_EEX2iZ,NbMax); //  x<0.20
   TH1D  *Hafb2_xmax_EEX2iZ = HstProjA("Hafb2_xmax_EEX2iZ",SCT_xc_EEX2iZ,NbMax); //  x<0.20
//
   TH1D  *Htot2_xmax_EEX2nZ = HstProjV("Htot2_xmax_EEX2nZ",SCT_xc_EEX2nZ,NbMax); //  x<0.20
   TH1D  *Hafb2_xmax_EEX2nZ = HstProjA("Hafb2_xmax_EEX2nZ",SCT_xc_EEX2nZ,NbMax); //  x<0.20
//
   TH1D  *Htot1_xmax_EEX2iZ = HstProjV("Htot1_xmax_EEX2iZ",SCN_xc_EEX2iZ,NbMax); //  x<0.02
   TH1D  *Hafb1_xmax_EEX2iZ = HstProjA("Hafb1_xmax_EEX2iZ",SCN_xc_EEX2iZ,NbMax); //  x<0.02
//
   TH1D  *Htot1_xmax_EEX2nZ = HstProjV("Htot1_xmax_EEX2nZ",SCN_xc_EEX2nZ,NbMax); //  x<0.02
   TH1D  *Hafb1_xmax_EEX2nZ = HstProjA("Hafb1_xmax_EEX2nZ",SCN_xc_EEX2nZ,NbMax); //  x<0.02


}// HisReMakeFOAMi


///////////////////////////////////////////////////////////////////////////////////
void HisReMakeKKMCi()
{
//------------------------------------------------------------------------
  cout<<" ========================= HisReMakeKKMCi =========================== "<<endl;
/////////////////////////////////////////////////////////////////////////////////////////////////

//****************************************************************************************
// v<vmax<0.20, c=cos(theta) with 100
//****************************************************************************************
  TH1D *HST_KKMC_NORMA9 = (TH1D*)DiskFileA95.Get("HST_KKMC_NORMA");
  int NbMax=45;

  cout<<"  Renormalizing  and reprocessing histograms from KKMC"<<endl;

  TH2D *sct_vTcPR_Ceex2i9  = (TH2D*)DiskFileA95.Get("sct_vTcPR_Ceex2");
  TH2D *sct_vTcPR_Ceex2n9  = (TH2D*)DiskFileA95.Get("sct_vTcPR_Ceex2n");

  TH2D *sct_vAcPR_Ceex2i9  = (TH2D*)DiskFileA95.Get("sct_vAcPR_Ceex2");
  TH2D *sct_vAcPR_Ceex2n9  = (TH2D*)DiskFileA95.Get("sct_vAcPR_Ceex2n");

  TH2D *sct_vTcPR_Ceex2i8  = (TH2D*)DiskFileA88.Get("sct_vTcPR_Ceex2");
  TH2D *sct_vTcPR_Ceex2n8  = (TH2D*)DiskFileA88.Get("sct_vTcPR_Ceex2n");

  TH2D *sct_vAcPR_Ceex2i8  = (TH2D*)DiskFileA88.Get("sct_vAcPR_Ceex2");
  TH2D *sct_vAcPR_Ceex2n8  = (TH2D*)DiskFileA88.Get("sct_vAcPR_Ceex2n");

//
  TH2D *sct_vTcPR_Ceex2iZ  = (TH2D*)DiskFileA91.Get("sct_vTcPR_Ceex2");
  TH2D *sct_vTcPR_Ceex2nZ  = (TH2D*)DiskFileA91.Get("sct_vTcPR_Ceex2n");

  TH2D *sct_vAcPR_Ceex2iZ  = (TH2D*)DiskFileA91.Get("sct_vAcPR_Ceex2");
  TH2D *sct_vAcPR_Ceex2nZ  = (TH2D*)DiskFileA91.Get("sct_vAcPR_Ceex2n");

//==================
  TH1D *HAfb9_vB_Ceex2   = HstProjA("HAfb9_vB_Ceex2", sct_vTcPR_Ceex2i9,NbMax); //  x<0.20
  TH1D *HAfb9_vB_Ceex2n  = HstProjA("HAfb9_vB_Ceex2n",sct_vTcPR_Ceex2n9,NbMax); //  x<0.20
//
  TH1D *HAfb9_vA_Ceex2   = HstProjA("HAfb9_vA_Ceex2",  sct_vAcPR_Ceex2i9,NbMax); //  x<0.20
  TH1D *HAfb9_vA_Ceex2n  = HstProjA("HAfb9_vA_Ceex2n", sct_vAcPR_Ceex2n9,NbMax); //  x<0.20
//==================
  TH1D *HAfbZ_vB_Ceex2   = HstProjA("HAfbZ_vB_Ceex2", sct_vTcPR_Ceex2iZ,NbMax); //  x<0.20
  TH1D *HAfbZ_vB_Ceex2n  = HstProjA("HAfbZ_vB_Ceex2n",sct_vTcPR_Ceex2nZ,NbMax); //  x<0.20
//
  TH1D *HAfbZ_vA_Ceex2   = HstProjA("HAfbZ_vA_Ceex2",  sct_vAcPR_Ceex2iZ,NbMax); //  x<0.20
  TH1D *HAfbZ_vA_Ceex2n  = HstProjA("HAfbZ_vA_Ceex2n", sct_vAcPR_Ceex2nZ,NbMax); //  x<0.20
//==================
  TH1D *HAfb8_vB_Ceex2   = HstProjA("HAfb8_vB_Ceex2", sct_vTcPR_Ceex2i8,NbMax); //  x<0.20
  TH1D *HAfb8_vB_Ceex2n  = HstProjA("HAfb8_vB_Ceex2n",sct_vTcPR_Ceex2n8,NbMax); //  x<0.20
//
  TH1D *HAfb8_vA_Ceex2   = HstProjA("HAfb8_vA_Ceex2",  sct_vAcPR_Ceex2i8,NbMax); //  x<0.20
  TH1D *HAfb8_vA_Ceex2n  = HstProjA("HAfb8_vA_Ceex2n", sct_vAcPR_Ceex2n8,NbMax); //  x<0.20

/*
  // NEW!!!!! NEW!!!!! NEW!!!!! NEW!!!!! NEW!!!!! NEW!!!!!
  TH2D *sct_CutAcPR_Ceex2  = (TH2D*)DiskFileA95.Get("sct_CutAcPR_Ceex2");
  TH2D *sct_CutAcPR_Ceex2n = (TH2D*)DiskFileA95.Get("sct_CutAcPR_Ceex2n");

  TH1D *HAfb9_vA_Ceex2   = HstProjA("HAfb9_vA_Ceex2",  sct_CutAcPR_Ceex2,NbMax);  //  x<0.20
  TH1D *HAfb9_vA_Ceex2n  = HstProjA("HAfb9_vA_Ceex2n", sct_CutAcPR_Ceex2n,NbMax); //  x<0.20
*/
  cout<<" ===================== End of HisReMakeKKMCi ======================== "<<endl;

}// HisReMakeKKMCi


///////////////////////////////////////////////////////////////////////////////////
void KKsemMakeHisto(){
  // Here we produce semianalytical plots using KKsem program, No plotting
  //------------------------------------------------------------------------
  cout<<"==================================================================="<<endl;
  cout<<"================ KKsem MakeHisto  BEGIN ============================"<<endl;
  // MuMu  Sigma(vmax) and AFB(vmax) with limited c=cos(theta)
  //
  long KF=13; // muon
  long KeyDis, KeyFob;
  char chak[5];
//------------------------------------------------------------------------
// Tempate is for 0<vmax<0.02,  cos(theta) < gCosTheta=0.90
  TH1D *hstVtemplate = (TH1D*)DiskFileB.Get("Hafb1_xmax_EEX2i9");
//------------------------------------------------------------------------
//  sprintf(chak,"XCHI2");  // ISR*FSR Mff
//  KeyDis = 302302;        // ISR*FSR O(alf2)
//
  sprintf(chak,"VCHI2");  // ISR only Mff
  KeyDis = 302;   // ISR O(alf2)
  //  KeyDis = 303;   // ISR O(alf3)

  // ************ 95GeV ************
  KKplot LibSem9("LibSem9");
  LibSem9.Initialize(DiskFileA95);  // for non-farm case
  // Classic test KKMC-KKsem for IFI off
  KeyFob=    0; // With EW (BornV_Dizet) With integration OK!
  TH1D *vcum_KKsem_95 =(TH1D*)hstVtemplate->Clone("vcum_KKsem_95");
  TH1D *afbv_KKsem_95 =(TH1D*)hstVtemplate->Clone("afbv_KKsem_95");
  LibSem9.VVmake( vcum_KKsem_95, afbv_KKsem_95, KF, chak, KeyDis, KeyFob, gCosTheta);
  // New test KKMC-KKsem for IFI off based on KKsem_Born Gmu scheme
  KeyFob=  -100; //
  //KeyFob=  -200; //
  TH1D *vcum_KKsem_95b =(TH1D*)hstVtemplate->Clone("vcum_KKsem_95b");
  TH1D *afbv_KKsem_95b =(TH1D*)hstVtemplate->Clone("afbv_KKsem_95b");
  LibSem9.VVmake( vcum_KKsem_95b, afbv_KKsem_95b, KF, chak, KeyDis, KeyFob, gCosTheta);
  // ************  88GeV ************
  KKplot LibSem8("LibSem8");
  LibSem8.Initialize(DiskFileA88);  // for non-farm case
  // Classic test KKMC-KKsem for IFI off
  KeyFob=    0; // With EW (BornV_Dizet) With integration OK!
  TH1D *vcum_KKsem_88 =(TH1D*)hstVtemplate->Clone("vcum_KKsem_88");
  TH1D *afbv_KKsem_88 =(TH1D*)hstVtemplate->Clone("afbv_KKsem_88");
  LibSem8.VVmake( vcum_KKsem_88, afbv_KKsem_88, KF, chak, KeyDis, KeyFob, gCosTheta);
  KeyFob=  -100; //
  //KeyFob=  -200; //
  TH1D *vcum_KKsem_88b =(TH1D*)hstVtemplate->Clone("vcum_KKsem_88b");
  TH1D *afbv_KKsem_88b =(TH1D*)hstVtemplate->Clone("afbv_KKsem_88b");
  LibSem9.VVmake( vcum_KKsem_88b, afbv_KKsem_88b, KF, chak, KeyDis, KeyFob, gCosTheta);
  // ************  88GeV ************
  KKplot LibSemZ("LibSemZ");
  LibSemZ.Initialize(DiskFileA91);  // for non-farm case
  // Classic test KKMC-KKsem for IFI off
  KeyFob=    0; // With EW (BornV_Dizet) With integration OK!
  TH1D *vcum_KKsem_91 =(TH1D*)hstVtemplate->Clone("vcum_KKsem_91");
  TH1D *afbv_KKsem_91 =(TH1D*)hstVtemplate->Clone("afbv_KKsem_91");
  LibSemZ.VVmake( vcum_KKsem_91, afbv_KKsem_91, KF, chak, KeyDis, KeyFob, gCosTheta);
  KeyFob=  -100; //
  //KeyFob=  -200; //
  TH1D *vcum_KKsem_91b =(TH1D*)hstVtemplate->Clone("vcum_KKsem_91b");
  TH1D *afbv_KKsem_91b =(TH1D*)hstVtemplate->Clone("afbv_KKsem_91b");
  LibSemZ.VVmake( vcum_KKsem_91b, afbv_KKsem_91b, KF, chak, KeyDis, KeyFob, gCosTheta);

  cout<<"================ KKsem MakeHisto ENDs ============================="<<endl;
  cout<<"==================================================================="<<endl;
//------------------------------------------------------------------------
//------------------------------------------------------------------------
}//  KKsemMakeHisto




///////////////////////////////////////////////////////////////////////////////////
void FigExp9w()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigExp9w =========================== "<<endl;

  TH1D *HAfb9_vA_Ceex2  = (TH1D*)DiskFileB.Get("HAfb9_vA_Ceex2");
  TH1D *HAfb9_vA_Ceex2n = (TH1D*)DiskFileB.Get("HAfb9_vA_Ceex2n");

  TH1D *HAfb9_vB_Ceex2  = (TH1D*)DiskFileB.Get("HAfb9_vB_Ceex2");
  TH1D *HAfb9_vB_Ceex2n = (TH1D*)DiskFileB.Get("HAfb9_vB_Ceex2n");

//  TH1D *afbv_KKsem_95  = (TH1D*)DiskFileB.Get("afbv_KKsem_95");
//  TH1D *afbv_KKsem_95b = (TH1D*)DiskFileB.Get("afbv_KKsem_95b");

  TH1D *Hafb2_xmax_EEX2i9 =(TH1D*)DiskFileB.Get("Hafb2_xmax_EEX2i9");
  TH1D *Hafb2_xmax_EEX2n9 =(TH1D*)DiskFileB.Get("Hafb2_xmax_EEX2n9");

  TLatex *CaptNDC = new TLatex(); CaptNDC->SetNDC(); // !!!
  CaptNDC->SetTextSize(0.037);
  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cExp9w = new TCanvas("cExp9w","cExp9w", gXcanv,  gYcanv,   1200,  600);
  //                            Name    Title            xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  gXcanv += 50; gYcanv += 50;
  cExp9w->SetFillColor(10);
  cExp9w->Divide( 2,  0);
  //cExp9w->Divide( 2,  0,     0.0,     0.0,   10);
  //              nx, ny, xmargin, ymargin, color
  //////////////////////////////////////////////
  cExp9w->cd(1);
  TH1D *HST1;
//  HST1 =HAfb9_vB_Ceex2;
//  HST1 =HAfb9_vB_Ceex2n;
  HST1 =Hafb2_xmax_EEX2i9;

  HST1->SetStats(0);
  HST1->SetTitle(0);
  HST1->GetXaxis()->SetNdivisions(8);
//  HST1->SetMaximum(0.350);
//  HST1->SetMinimum(0.250); // 95GeV
  HST1->DrawCopy("h");

  double ycapt =0.91, xcapt=0.3;
  CaptNDC->DrawLatex(xcapt-0.1,ycapt," #sqrt{s} = 94.3GeV");
  ycapt += -0.047;

  PlotSame2(HAfb9_vB_Ceex2,     ycapt, kBlue,      0.02, "(a)", "KKMC Bcut IFIon ");
  PlotSame2(HAfb9_vB_Ceex2n,    ycapt, kBlue,      0.02, "(b)", "KKMC Bcut IFIoff");

  PlotSame2(HAfb9_vA_Ceex2,     ycapt, kBlack,     0.007,"(c)", "KKMC Acut IFIon");
  PlotSame2(HAfb9_vA_Ceex2n,    ycapt, kBlack,     0.08, "(d)", "KKMC Acut IFIoff");

  PlotSame2(Hafb2_xmax_EEX2i9,  ycapt, kRed,       0.15, "(e)", "KKfoam3 IFIon");
  PlotSame2(Hafb2_xmax_EEX2n9,  ycapt, kRed,       0.15, "(f)", "KKfoam3 IFIoff");

//  PlotSame2(afbv_KKsem_95,      ycapt, kGreen,     0.05, "(g)", "KKsem DIZET");
//  PlotSame2(afbv_KKsem_95b,     ycapt, kMagenta,   0.05, "(h)", "KKsem Gmu");

  CaptNDC->DrawLatex(0.04,0.95,"A_{FB}(v)  ");
  CaptNDC->DrawLatex(0.60,0.02,"v_{max} ");

  ///////////////////////////////////////////////
  cExp9w->cd(2);
  TH1D *Hst9vB_IFI   = HstDiff("Hst9vB_IFI",   HAfb9_vB_Ceex2,   HAfb9_vB_Ceex2n,   kBlue);
  TH1D *Hst9vA_IFI   = HstDiff("Hst9vA_IFI",   HAfb9_vA_Ceex2,   HAfb9_vA_Ceex2n,   kBlack);
  TH1D *Hst9vE_IFI   = HstDiff("Hst9vE_IFI",   Hafb2_xmax_EEX2i9,Hafb2_xmax_EEX2n9, kRed);

  TH1D *HST2;
//  HST2 = Hst9vB_IFI;
  HST2 = Hst9vE_IFI;

  HST2->SetStats(0);
  HST2->SetTitle(0);
//  HST2->GetXaxis()->SetNdivisions(8);
//  HST2->SetMaximum(0.07);
//  HST2->SetMinimum(0.00);
  HST2->DrawCopy("h");

  ycapt =0.91, xcapt=0.3;
  CaptNDC->DrawLatex(xcapt-0.1,ycapt," #sqrt{s} = 94.3GeV, IFI");
  ycapt += -0.047;

  PlotSame2(Hst9vB_IFI,       ycapt, kBlue,    0.03, "(a)", "KKMC Bcut IFI ");
  PlotSame2(Hst9vA_IFI,       ycapt, kBlack,   0.05, "(b)", "KKMC Acut IFI ");
  PlotSame2(Hst9vE_IFI,       ycapt, kRed,     0.05, "(c)", "KKFoam3   IFI ");

//  PlotSame3(HAfb9_vB_Ceex21,   xcapt, ycapt, kBlue,     25, 1, "(B) Idealized");
//  PlotSame3(HAfb9_vA_Ceex21,   xcapt, ycapt, kBlack,    25, 1, "(A) Realistic");

  CaptNDC->DrawLatex(0.04,0.95,"A_{FB}^{IFI}(v)  ");
  CaptNDC->DrawLatex(0.60,0.02,"v_{max} ");

  //
  cExp9w->cd();
  //
  cExp9w->SaveAs("cExp9w.pdf");
//
}// FigExp9w



///////////////////////////////////////////////////////////////////////////////////
void FigExp8w()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigExp8w =========================== "<<endl;

  TH1D *HAfb8_vA_Ceex2  = (TH1D*)DiskFileB.Get("HAfb8_vA_Ceex2");
  TH1D *HAfb8_vA_Ceex2n = (TH1D*)DiskFileB.Get("HAfb8_vA_Ceex2n");

  TH1D *HAfb8_vB_Ceex2  = (TH1D*)DiskFileB.Get("HAfb8_vB_Ceex2");
  TH1D *HAfb8_vB_Ceex2n = (TH1D*)DiskFileB.Get("HAfb8_vB_Ceex2n");

  TH1D *afbv_KKsem_88  = (TH1D*)DiskFileB.Get("afbv_KKsem_88");
  TH1D *afbv_KKsem_88b = (TH1D*)DiskFileB.Get("afbv_KKsem_88b");

  TH1D *Hafb2_xmax_EEX2i8 =(TH1D*)DiskFileB.Get("Hafb2_xmax_EEX2i8");
  TH1D *Hafb2_xmax_EEX2n8 =(TH1D*)DiskFileB.Get("Hafb2_xmax_EEX2n8");

  TLatex *CaptNDC = new TLatex(); CaptNDC->SetNDC(); // !!!
  CaptNDC->SetTextSize(0.037);
  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cExp8w = new TCanvas("cExp8w","cExp8w", gXcanv,  gYcanv,   1200,  600);
  //                            Name    Title            xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  gXcanv += 50; gYcanv += 50;
  cExp8w->SetFillColor(10);
  cExp8w->Divide( 2,  0);
  //cExp8w->Divide( 2,  0,     0.0,     0.0,   10);
  //              nx, ny, xmargin, ymargin, color
  //////////////////////////////////////////////
  cExp8w->cd(1);
  TH1D *HST1;
  HST1 =HAfb8_vB_Ceex2;
//  HST1 =HAfb8_vB_Ceex2n;
//  HST1 =Hafb2_xmax_EEX2i8;

  HST1->SetStats(0);
  HST1->SetTitle(0);
  HST1->GetXaxis()->SetNdivisions(8);
  HST1->SetMaximum(-0.230);
  HST1->SetMinimum(-0.300); // 88GeV
  HST1->DrawCopy("h");

  double ycapt =0.91, xcapt=0.3;
  CaptNDC->DrawLatex(xcapt-0.1,ycapt," #sqrt{s} = 87.9GeV");
  ycapt += -0.047;

  PlotSame2(HAfb8_vB_Ceex2,     ycapt, kBlue,      0.05, "(a)", "KKMC Bcut IFIon ");
  PlotSame2(HAfb8_vB_Ceex2n,    ycapt, kBlue,      0.05, "(b)", "KKMC Bcut IFIoff");

  PlotSame2(HAfb8_vA_Ceex2,     ycapt, kBlack,     0.05, "(c)", "KKMC Acut IFIon");
  PlotSame2(HAfb8_vA_Ceex2n,    ycapt, kBlack,     0.09, "(d)", "KKMC Acut IFIoff");

  PlotSame2(Hafb2_xmax_EEX2i8,  ycapt, kRed,       0.05, "(e)", "KKfoam3 IFIon");
  PlotSame2(Hafb2_xmax_EEX2n8,  ycapt, kRed,       0.03, "(f)", "KKfoam3 IFIoff");

//  PlotSame2(afbv_KKsem_88,      ycapt, kGreen,     0.013, "(g)", "KKsem DIZET");
//  PlotSame2(afbv_KKsem_88b,     ycapt, kMagenta,   0.009, "(h)", "KKsem Gmu");

  CaptNDC->DrawLatex(0.04,0.95,"A_{FB}(v)  ");
  CaptNDC->DrawLatex(0.60,0.02,"v_{max} ");

  ///////////////////////////////////////////////
  cExp8w->cd(2);
  TH1D *Hst8vB_IFI   = HstDiff("Hst8vB_IFI",   HAfb8_vB_Ceex2,   HAfb8_vB_Ceex2n,   kBlue);
  TH1D *Hst8vA_IFI   = HstDiff("Hst8vA_IFI",   HAfb8_vA_Ceex2,   HAfb8_vA_Ceex2n,   kBlack);
  TH1D *Hst8vE_IFI   = HstDiff("Hst8vE_IFI",   Hafb2_xmax_EEX2i8,Hafb2_xmax_EEX2n8, kRed);

  TH1D *HST2;
  HST2 = Hst8vB_IFI;

  HST2->SetStats(0);
  HST2->SetTitle(0);
  HST2->GetXaxis()->SetNdivisions(8);
  HST2->SetMaximum(0.05); HST2->SetMinimum(0.00);
  HST2->DrawCopy("h");

  ycapt =0.91, xcapt=0.3;
  CaptNDC->DrawLatex(xcapt-0.1,ycapt," #sqrt{s} = 87.9GeV, IFI");
  ycapt += -0.047;

  PlotSame2(Hst8vB_IFI,       ycapt, kBlue,    0.15, "(a)", "KKMC Bcut IFI ");
  PlotSame2(Hst8vA_IFI,       ycapt, kBlack,   0.10, "(b)", "KKMC Acut IFI ");
  PlotSame2(Hst8vE_IFI,       ycapt, kRed,     0.15, "(c)", "KKFoam3   IFI ");

//  PlotSame3(HAfb8_vB_Ceex21,   xcapt, ycapt, kBlue,     25, 1, "(B) Idealized");
//  PlotSame3(HAfb8_vA_Ceex21,   xcapt, ycapt, kBlack,    25, 1, "(A) Realistic");

  CaptNDC->DrawLatex(0.04,0.95,"A_{FB}^{IFI}(v)  ");
  CaptNDC->DrawLatex(0.60,0.02,"v_{max} ");

  //
  cExp8w->cd();
  //
  cExp8w->SaveAs("cExp8w.pdf");
//
}// FigExp8w



///////////////////////////////////////////////////////////////////////////////////
void FigExpZw()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigExpZw =========================== "<<endl;

  TH1D *HAfbZ_vA_Ceex2  = (TH1D*)DiskFileB.Get("HAfbZ_vA_Ceex2");
  TH1D *HAfbZ_vA_Ceex2n = (TH1D*)DiskFileB.Get("HAfbZ_vA_Ceex2n");

  TH1D *HAfbZ_vB_Ceex2  = (TH1D*)DiskFileB.Get("HAfbZ_vB_Ceex2");
  TH1D *HAfbZ_vB_Ceex2n = (TH1D*)DiskFileB.Get("HAfbZ_vB_Ceex2n");

  TH1D *afbv_KKsem_91  = (TH1D*)DiskFileB.Get("afbv_KKsem_91");
  TH1D *afbv_KKsem_91b = (TH1D*)DiskFileB.Get("afbv_KKsem_91b");

  TH1D *Hafb2_xmax_EEX2iZ =(TH1D*)DiskFileB.Get("Hafb2_xmax_EEX2iZ");
  TH1D *Hafb2_xmax_EEX2nZ =(TH1D*)DiskFileB.Get("Hafb2_xmax_EEX2nZ");

  TLatex *CaptNDC = new TLatex(); CaptNDC->SetNDC(); // !!!
  CaptNDC->SetTextSize(0.037);
  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cExpZw = new TCanvas("cExpZw","cExpZw", gXcanv,  gYcanv,   1200,  600);
  //                            Name    Title            xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  gXcanv += 50; gYcanv += 50;
  cExpZw->SetFillColor(10);
  cExpZw->Divide( 2,  0);
  //cExpZw->Divide( 2,  0,     0.0,     0.0,   10);
  //              nx, ny, xmargin, ymargin, color
  //////////////////////////////////////////////
  cExpZw->cd(1);
  TH1D *HST1;
  HST1 =HAfbZ_vB_Ceex2;
//  HST1 =HAfbZ_vB_Ceex2n;
//  HST1 =Hafb2_xmax_EEX2iZ;

  HST1->SetStats(0);
  HST1->SetTitle(0);
  HST1->GetXaxis()->SetNdivisions(8);
  HST1->SetMaximum( 0.055);
  HST1->SetMinimum( 0.000); // 88GeV
  HST1->DrawCopy("h");

  double ycapt =0.91, xcapt=0.3;
  CaptNDC->DrawLatex(xcapt-0.1,ycapt," #sqrt{s} = 91.20GeV");
  ycapt += -0.047;

  PlotSame2(HAfbZ_vB_Ceex2,     ycapt, kBlue,      0.01, "(a)", "KKMC Bcut IFIon ");
  PlotSame2(HAfbZ_vB_Ceex2n,    ycapt, kBlue,      0.01, "(b)", "KKMC Bcut IFIoff");

  PlotSame2(HAfbZ_vA_Ceex2,     ycapt, kBlack,     0.02, "(c)", "KKMC Acut IFIon");
  PlotSame2(HAfbZ_vA_Ceex2n,    ycapt, kBlack,     0.02, "(d)", "KKMC Acut IFIoff");

  PlotSame2(Hafb2_xmax_EEX2iZ,  ycapt, kRed,       0.05, "(e)", "KKfoam3 IFIon");
  PlotSame2(Hafb2_xmax_EEX2nZ,  ycapt, kRed,       0.05, "(f)", "KKfoam3 IFIoff");

//  PlotSame2(afbv_KKsem_88,      ycapt, kGreen,     0.013, "(g)", "KKsem DIZET");
//  PlotSame2(afbv_KKsem_88b,     ycapt, kMagenta,   0.009, "(h)", "KKsem Gmu");

  CaptNDC->DrawLatex(0.04,0.95,"A_{FB}(v)  ");
  CaptNDC->DrawLatex(0.60,0.02,"v_{max} ");

  ///////////////////////////////////////////////
  cExpZw->cd(2);
  TH1D *HstZvB_IFI   = HstDiff("HstZvB_IFI",   HAfbZ_vB_Ceex2,   HAfbZ_vB_Ceex2n,   kBlue);
  TH1D *HstZvA_IFI   = HstDiff("HstZvA_IFI",   HAfbZ_vA_Ceex2,   HAfbZ_vA_Ceex2n,   kBlack);
  TH1D *HstZvE_IFI   = HstDiff("HstZvE_IFI",   Hafb2_xmax_EEX2iZ,Hafb2_xmax_EEX2nZ, kRed);

  TH1D *HST2;
  HST2 = HstZvB_IFI;

  HST2->SetStats(0);
  HST2->SetTitle(0);
  HST2->GetXaxis()->SetNdivisions(8);
  HST2->SetMaximum(0.035); HST2->SetMinimum(0.00);
  HST2->DrawCopy("h");

  ycapt =0.91, xcapt=0.3;
  CaptNDC->DrawLatex(xcapt-0.1,ycapt," #sqrt{s} = 91.20GeV, IFI");
  ycapt += -0.047;

  PlotSame2(HstZvB_IFI,       ycapt, kBlue,    0.03, "(a)", "KKMC Bcut IFI ");
  PlotSame2(HstZvA_IFI,       ycapt, kBlack,   0.05, "(b)", "KKMC Acut IFI ");
  PlotSame2(HstZvE_IFI,       ycapt, kRed,     0.02, "(c)", "KKFoam3   IFI ");

//  PlotSame3(HAfbZ_vB_Ceex21,   xcapt, ycapt, kBlue,     25, 1, "(B) Idealized");
//  PlotSame3(HAfbZ_vA_Ceex21,   xcapt, ycapt, kBlack,    25, 1, "(A) Realistic");

  CaptNDC->DrawLatex(0.04,0.95,"A_{FB}^{IFI}(v)  ");
  CaptNDC->DrawLatex(0.60,0.02,"v_{max} ");

  //
  cExpZw->cd();
  //
  cExpZw->SaveAs("cExpZw.pdf");
//
}// FigExpZw



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
  HisReMakeFOAMi();
  HisReMakeKKMCi();       // reprocessing MC histos
  KKsemMakeHisto();
  //========== PLOTTING ==========
  //
  FigExp9w(); // 94GeV
  FigExp8w(); // 88GeV
  FigExpZw(); // 91GeV
  //
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


