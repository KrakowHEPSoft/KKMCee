//////////////////////////////////////////////////////////////////////
//    make Plot1
//////////////////////////////////////////////////////////////////////

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

#include "KKsem.h"

// ROOT headers
#include "TROOT.h"
#include "TFile.h"
//=============================================================================
//  ROOT  ROOT ROOT   ROOT  ROOT  ROOT  ROOT  ROOT  ROOT  ROOT   ROOT   ROOT 
//=============================================================================
//TFile DiskFileA("../test0/rmain.root_CEEX.31M"); /// new
//TFile DiskFileA("../test0/rmain.root.2.5M"); // KeyElw=0
//TFile DiskFileA("../test0/rmain.root.6M.EW"); // KeyElw=1
TFile DiskFileA("../test0/rmain.root");
TFile DiskFileB("RhoSemi.root","RECREATE","histograms");
//=============================================================================

Double_t sqr( const Double_t x ){ return x*x;};
// Auxiliary procedures for plotting
#include "HisNorm.h"
#include "Marker.h"

///////////////////////////////////////////////////////////////////////////////////
void HistNormalize(){
  //
  cout<<"----------------------------- HistNormalize ------------------------------------"<<endl;
  DiskFileA.ls("");
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vTrueMain") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vTrueCeex2") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vAlepCeex2") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_vXGenCeex2") );
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_Cost1Ceex2") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_CosPLCeex2") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_CosPRCeex2") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_CosPREex2") );
  //
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sca_vTcPR_Ceex2") );
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sca_vTcPR_Eex2") );
  //
}

///////////////////////////////////////////////////////////////////////////////////
void KKsemMakeHisto(){
	// Here we produce semianalytical plots using KKsem program, No plotting
	// also some MC histos are preprocessed
	//------------------------------------------------------------------------
  cout<<"==================================================================="<<endl;
  cout<<"================ KKsem MakeHisto  BEGIN ============================"<<endl;
  // initilalization of KKsem
  KKsem LibSem;
  LibSem.Initialize(DiskFileA);
  //
  long KF=13; // muon
  long KeyDis, KeyFob;
  char chak[5];
  //KeyDis = 302;   // ISR O(alf2)
  //KeyDis = 304;   // ISR O(alf3) GribovLL
  //KeyDis = 303;   // ISR O(alf3)
  //KeyDis = 305;   // ISR O(alf3) GribovLL +NLL
  //KeyDis = 662;  // Unexp ????
  //KeyDis = 302302;   // ISR*FSR O(alf3)
  //
  KeyFob=   10; // BornV_Dizet, with EW and without integration ???
  KeyFob=  -11; // BornV_Simple, for KeyLib=0, NO EW, NO integration OK
  KeyFob=  -10; // KKsem_BornV, NO EW, NO integration OK!
  KeyFob= -100; // KKsem_BornV, NO EW, WITH integration, OK
  KeyFob=    0; // With EW (BornV_Dizet) With integration OK!
//------------------------------------------------------------------------
//   MuMu  dsigma/dv
//------------------------------------------------------------------------  
  TH1D *hstVtemplate = (TH1D*)DiskFileA.Get("hst_vTrueMain");
  TH1D *hstCtemplate = (TH1D*)DiskFileA.Get("hst_Cost1Ceex2");
  // ISR*FSR
  KeyDis = 302302;        // ISR*FSR O(alf2)
  sprintf(chak,"XRHO2");  // ISR*FSR Mff
  TH1D *vdis_ISR2_FSR2 =(TH1D*)hstVtemplate->Clone("vdis_ISR2_FSR2");
  LibSem.VVplot(vdis_ISR2_FSR2, KF, chak, KeyDis, KeyFob);
  // ISR only
  KeyDis = 303;           // ISR O(alf3)
  sprintf(chak,"VRHO2");  // ISR only
  TH1D *vdis_ISR2 =(TH1D*)hstVtemplate->Clone("vdis_ISR2");
  LibSem.VVplot(vdis_ISR2, KF, chak, KeyDis, KeyFob);
//------------------------------------------------------------------------
//   MuMu  Sigma(vmax)
//------------------------------------------------------------------------  
  // ISR*FSR
  KeyDis = 302302;        // ISR*FSR O(alf2)
  sprintf(chak,"XCHI2");  // ISR*FSR Mff
  TH1D *vcum_ISR2_FSR2 =(TH1D*)hstVtemplate->Clone("vcum_ISR2_FSR2");
  LibSem.VVplot(vcum_ISR2_FSR2, KF, chak, KeyDis, KeyFob);
  // with costhe cut
  kksem_setcrange_(-22.0/25, 22.0/25);
  TH1D *vcum2_ISR2_FSR2 =(TH1D*)hstVtemplate->Clone("vcum2_ISR2_FSR2");
  LibSem.VVplot(vcum2_ISR2_FSR2, KF, chak, KeyDis, KeyFob);
  kksem_setcrange_(-1.0,1.0); // back to normal
//------------------------------------------------------------------------
//   MuMu  dsigma/dCosTheta, KKsem ISR only
//------------------------------------------------------------------------  
  // ISR only
  KeyDis = 303;           // ISR O(alf3)
  sprintf(chak,"VRHO2");  // ISR only
  //KeyDis = 302302;        // ISR*FSR O(alf2)
  //sprintf(chak,"XRHO2");  // ISR*FSR Mff
  KeyFob= -100; // KKsem_BornV, WITH integration, OK
  double vmin=0.0;
  double vmax=0.02;
  TH1D        *cdisKS_ISR2 =(TH1D*)hstCtemplate->Clone("cdisKS_ISR2");
  LibSem.Cplot(cdisKS_ISR2, KF, chak, KeyDis, KeyFob, vmin,vmax);
  //
  KeyFob=    0; // With EW (BornV_Dizet) With integration OK!
  TH1D        *cdisDZ_ISR2 =(TH1D*)hstCtemplate->Clone("cdisDZ_ISR2");
  LibSem.Cplot(cdisDZ_ISR2, KF, chak, KeyDis, KeyFob, vmin,vmax);
  //----------------------------------------------
  // for x-check only, the distribution without Z
  long KeyZet=0;
  kksem_setkeyzet_(KeyZet);
  TH1D *cdis_ISR2_Zoff =(TH1D*)hstCtemplate->Clone("cdis_ISR2_Zoff");
  LibSem.Cplot(cdis_ISR2_Zoff, KF, chak, KeyDis, KeyFob, vmin,vmax);

  //****************************************************************************************
  // Pure MC reprocessing part
  //
  TH2D *sca_vTcPR_Ceex2 = (TH2D*)DiskFileA.Get("sca_vTcPR_Ceex2");
  TH2D *sca_vTcPR_Eex2  = (TH2D*)DiskFileA.Get("sca_vTcPR_Eex2");
  TH2D *sca_vTcPR_Ceex2n = (TH2D*)DiskFileA.Get("sca_vTcPR_Ceex2n");

  int nbMax=0; // this is cut on costheta, 0= no cut
  //nbMax=45;    // vMax = 45/50=0.9
  TH1D                    *Hcth_vTcPR_Ceex2_vmax90, *Hcas_vTcPR_Ceex2_vmax90;
  ProjC( sca_vTcPR_Ceex2,  Hcth_vTcPR_Ceex2_vmax90,  Hcas_vTcPR_Ceex2_vmax90, nbMax);
  Hcth_vTcPR_Ceex2_vmax90->SetName("Hcth_vTcPR_Ceex2_vmax90");
  Hcas_vTcPR_Ceex2_vmax90->SetName("Hcas_vTcPR_Ceex2_vmax90");
  //
  nbMax=15;    // vMax = 15/50=0.3
  TH1D                    *Hcth_vTcPR_Ceex2_vmax30, *Hcas_vTcPR_Ceex2_vmax30;
  ProjC( sca_vTcPR_Ceex2,  Hcth_vTcPR_Ceex2_vmax30,  Hcas_vTcPR_Ceex2_vmax30, nbMax);
  Hcth_vTcPR_Ceex2_vmax30->SetName("Hcth_vTcPR_Ceex2_vmax30");
  Hcas_vTcPR_Ceex2_vmax30->SetName("Hcas_vTcPR_Ceex2_vmax30");
  //
  nbMax=5;     // vMax = 5/50=0.1
  TH1D                    *Hcth_vTcPR_Ceex2_vmax10, *Hcas_vTcPR_Ceex2_vmax10;
  ProjC( sca_vTcPR_Ceex2,  Hcth_vTcPR_Ceex2_vmax10,  Hcas_vTcPR_Ceex2_vmax10, nbMax);
  Hcth_vTcPR_Ceex2_vmax10->SetName("Hcth_vTcPR_Ceex2_vmax10");
  Hcas_vTcPR_Ceex2_vmax10->SetName("Hcas_vTcPR_Ceex2_vmax10");
  //
  nbMax=1;     // vMax = 5/50=0.02
  TH1D                    *Hcth_vTcPR_Ceex2_vmax02, *Hcas_vTcPR_Ceex2_vmax02;
  ProjC( sca_vTcPR_Ceex2,  Hcth_vTcPR_Ceex2_vmax02,  Hcas_vTcPR_Ceex2_vmax02, nbMax);
  Hcth_vTcPR_Ceex2_vmax02->SetName("Hcth_vTcPR_Ceex2_vmax02");
  Hcas_vTcPR_Ceex2_vmax02->SetName("Hcas_vTcPR_Ceex2_vmax02");
  // IFI off
  TH1D                     *Hcth_vTcPR_Ceex2n_vmax02, *Hcas_vTcPR_Ceex2n_vmax02;
  ProjC( sca_vTcPR_Ceex2n,  Hcth_vTcPR_Ceex2n_vmax02,  Hcas_vTcPR_Ceex2n_vmax02, nbMax);
  Hcth_vTcPR_Ceex2n_vmax02->SetName("Hcth_vTcPR_Ceex2n_vmax02");
  Hcas_vTcPR_Ceex2n_vmax02->SetName("Hcas_vTcPR_Ceex2n_vmax02");
  //
  TH1D                  *casyKS_ISR2;
  MakeAFB(cdisKS_ISR2,   casyKS_ISR2);
  casyKS_ISR2->SetName("casyKS_ISR2");
  //
  TH1D                  *casyDZ_ISR2;
  MakeAFB(cdisDZ_ISR2,   casyDZ_ISR2);
  casyDZ_ISR2->SetName("casyDZ_ISR2");
  //  
  cout<<"================ KKsem MakeHisto ENDs ============================="<<endl;
  cout<<"==================================================================="<<endl;
//------------------------------------------------------------------------  
//------------------------------------------------------------------------  
}//  KKsemMakeHisto


///////////////////////////////////////////////////////////////////////////////////
void FigScatA()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigScat =========================== "<<endl;
  // renormalize histograms in nanobarns
  Double_t CMSene;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");

  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  //
  TH2D *sca_vTcPR_Ceex2 = (TH2D*)DiskFileA.Get("sca_vTcPR_Ceex2");
  TH2D *sca_vTcPR_Eex2  = (TH2D*)DiskFileA.Get("sca_vTcPR_Eex2");
  //
  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cScatA = new TCanvas("cScatA","2dim big picture", 70,  50,   1200,  600);
  //                            Name    Title            xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  cScatA->SetFillColor(10);
  cScatA->Divide( 2,  0);
  //cScatA->Divide( 2,  0,     0.0,     0.0,   10);
  //              nx, ny, xmargin, ymargin, color
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
  cScatA->cd(1);
  gPad->SetLogz(); // !!!!!!
  gPad->SetTheta(25);
  gPad->SetPhi( -38);
  //
  double zmax = sca_vTcPR_Ceex2->GetMaximum();
  sca_vTcPR_Ceex2->SetMaximum(zmax*1.2);
  sca_vTcPR_Ceex2->SetMinimum(zmax*1e-3);
  sca_vTcPR_Ceex2->Draw(OptSurf);
  //-------------------------------------
  cScatA->cd(2);
  gPad->SetLogz(); // !!!!!!
  gPad->SetTheta(25);
  gPad->SetPhi( -38);
  sca_vTcPR_Eex2->SetMaximum(zmax*1.2);
  sca_vTcPR_Eex2->SetMinimum(zmax*1e-3);
  sca_vTcPR_Eex2->Draw(OptSurf);
  cScatA->cd();
}


///////////////////////////////////////////////////////////////////////////////////
void FigInfo()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigInfo =========================== "<<endl;
  // renormalize histograms in nanobarns
  Double_t CMSene;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");

  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  //
  TH1D *hst_nPhAll     = (TH1D*)DiskFileA.Get("hst_nPhAll");
  TH1D *hst_nPhVis     = (TH1D*)DiskFileA.Get("hst_nPhVis");
  TH1D *hst_weight     = (TH1D*)DiskFileA.Get("hst_weight");

  TH1D *hst_vTrueCeex2 = (TH1D*)DiskFileA.Get("hst_vTrueCeex2");
  TH1D *hst_vAlepCeex2 = (TH1D*)DiskFileA.Get("hst_vAlepCeex2");
  TH1D *hst_vXGenCeex2 = (TH1D*)DiskFileA.Get("hst_vXGenCeex2");

  TH1D *hst_Cost1Ceex2 = (TH1D*)DiskFileA.Get("hst_Cost1Ceex2");
  TH1D *hst_CosPLCeex2 = (TH1D*)DiskFileA.Get("hst_CosPLCeex2");
  TH1D *hst_CosPRCeex2 = (TH1D*)DiskFileA.Get("hst_CosPRCeex2");
//------------------------------------------------------------------------  
  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigInfo = new TCanvas("cFigInfo","FigInfo: general info ", 50, 80,    1000,  800);
  //                            Name    Title               xoff,yoff, WidPix,HeiPix
  cFigInfo->SetFillColor(10);
  ////////////////////////////////////////////////////////////////////////////////
  cFigInfo->Divide( 2,  2);
  //cFigInfo->Divide( 2,  2,     0.0,     0.0,   10);
  //              nx, ny, xmargin, ymargin, color
  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.04);
  //==========plot1==============
  cFigInfo->cd(1);
  hst_nPhVis->DrawCopy("h");
  hst_nPhAll->SetLineColor(2);
  hst_nPhAll->DrawCopy("hsame");
  //==========plot2==============
  cFigInfo->cd(2);
  hst_weight->DrawCopy("h");
  //==========plot3==============
  cFigInfo->cd(3);
  gPad->SetLogy(); // !!!!!!
  hst_vTrueCeex2->SetStats(0);
  hst_vTrueCeex2->SetTitle(0);
  hst_vTrueCeex2->DrawCopy("h");
  //
  hst_vAlepCeex2->SetLineColor(kRed);
  hst_vAlepCeex2->DrawCopy("hsame");
  //
  hst_vXGenCeex2->SetLineColor(4);
  hst_vXGenCeex2->DrawCopy("hsame");
  CaptT->DrawLatex(0.10,0.95,"d#sigma/dv (Ceex2); Black=Bare, Red=Aleph, Blue=Gener");
  //==========plot4==============
  cFigInfo->cd(4);
  //-----------------------------
  hst_Cost1Ceex2->SetStats(0);
  hst_Cost1Ceex2->SetTitle(0);
  hst_Cost1Ceex2->DrawCopy("h");
  //
  hst_CosPRCeex2->SetLineColor(2);
  hst_CosPRCeex2->DrawCopy("hsame");
  //
  hst_CosPLCeex2->SetLineColor(4);
  hst_CosPLCeex2->DrawCopy("hsame");  
  CaptT->DrawLatex(0.10,0.95,"d#sigma/dc (Ceex2); Black=   #theta_{1}, Red=PRD, Blue=PL");
  //----------------------------
  cFigInfo->cd();
}

///////////////////////////////////////////////////////////////////////////////////
void FigVtest()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigVtest =========================== "<<endl;
  Double_t CMSene;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  //
  TH1D *hst_vTrueCeex2 = (TH1D*)DiskFileA.Get("hst_vTrueCeex2");
  TH1D *hst_vXGenCeex2 = (TH1D*)DiskFileA.Get("hst_vXGenCeex2");
  //
  TH1D *vdis_ISR2      = (TH1D*)DiskFileB.Get("vdis_ISR2");
  TH1D *vdis_ISR2_FSR2 = (TH1D*)DiskFileB.Get("vdis_ISR2_FSR2");
  //
  TH2D *sca_vTcPR_Ceex2 = (TH2D*)DiskFileA.Get("sca_vTcPR_Ceex2");
  TH2D *sca_vTcPR_Eex2  = (TH2D*)DiskFileA.Get("sca_vTcPR_Eex2");
  //------------------------------------------------------------------------
  //****************************************************************************************
  //************************** Developement corner *****************************************
  //TH1D *hst_ProjV = (TH1D*)sca_vTcPR_Ceex2->ProjectionX("hst_projV",1,50,"e");
  TH1D *hst_projV;
  ProjX1(sca_vTcPR_Ceex2, hst_projV);
  hst_projV->SetName("hst_projV");
  //****************************************************************************************
  //****************************************************************************************
  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigVtest = new TCanvas("cFigVtest","FigVtest: photonic2", 50, 50,    1000, 800);
  //                            Name    Title               xoff,yoff, WidPix,HeiPix
  cFigVtest->SetFillColor(10);
  ////////////////////////////////////////////////////////////////////////////////
  cFigVtest->Divide( 2,  2);
  //cFigVtest->Divide( 2,  2,     0.0,     0.0,   10);
  //              nx, ny, xmargin, ymargin, color
  //////////////////////////////////////////////
  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.03);
  //==========plot1==============
  cFigVtest->cd(1);
  gPad->SetLogy(); // !!!!!!
  hst_vTrueCeex2->SetStats(0);
  hst_vTrueCeex2->SetTitle(0);
  hst_vTrueCeex2->DrawCopy("h");  // black
  //
  vdis_ISR2->SetLineColor(4); // blue
  vdis_ISR2->DrawCopy("hsame");
  //
  vdis_ISR2_FSR2->SetLineColor(6); // magenta
  vdis_ISR2_FSR2->DrawCopy("hsame");
  //
  hst_projV->SetLineColor(8); // green
  hst_projV->DrawCopy("same");
  //
  //HstNew->SetLineColor(7); // cyan
  //HstNew->DrawCopy("same");
  //
  CaptT->DrawLatex(0.02,0.95,"d#sigma/dv (Ceex2); Black=Bare, Red=Gener, Blue=ISR, Mag=ISR+FSR");
  //==========plot2==============
  cFigVtest->cd(2);
  hst_vTrueCeex2->Divide(vdis_ISR2_FSR2);
  hst_vTrueCeex2->SetStats(0);
  hst_vTrueCeex2->SetTitle(0);
  hst_vTrueCeex2->SetMinimum(0.85);
  hst_vTrueCeex2->SetMaximum(1.15);
  hst_vTrueCeex2->DrawCopy("h");
  //
  CaptT->DrawLatex(0.02,0.95,"d#sigma/dv (Ceex2); red: Gener/KKsemISR");
  //==========plot3==============
  cFigVtest->cd(3);
  gPad->SetLogy(); // !!!!!!
  hst_vXGenCeex2->SetStats(0);
  hst_vXGenCeex2->SetTitle(0);
  hst_vXGenCeex2->SetLineColor(2); // red
  hst_vXGenCeex2->DrawCopy("h");
  //
  vdis_ISR2->SetLineColor(4); // blue
  vdis_ISR2->DrawCopy("hsame");
  CaptT->DrawLatex(0.02,0.95,"d#sigma/dv (Ceex2); Red=Gener, Blue=ISR");
  //==========plot4==============
  cFigVtest->cd(4);
  hst_vXGenCeex2->Divide(vdis_ISR2);
  hst_vXGenCeex2->SetStats(0);
  hst_vXGenCeex2->SetTitle(0);
  hst_vXGenCeex2->SetMinimum(0.85);
  hst_vXGenCeex2->SetMaximum(1.15);
  hst_vXGenCeex2->DrawCopy("h");  // black
  //----------------------------
  cFigVtest->cd();
  //================================================
}//FigVtest

///////////////////////////////////////////////////////////////////////////////////
void FigCtest()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigCtest =========================== "<<endl;
  TH1D *hst_Cost1Ceex2 = (TH1D*)DiskFileA.Get("hst_Cost1Ceex2");
  TH1D *hst_CosPLCeex2 = (TH1D*)DiskFileA.Get("hst_CosPLCeex2");
  TH1D *hst_CosPRCeex2 = (TH1D*)DiskFileA.Get("hst_CosPRCeex2");
  TH1D *hst_CosPREex2  = (TH1D*)DiskFileA.Get("hst_CosPREex2");
  //
  TH1D *cdisKS_ISR2      = (TH1D*)DiskFileB.Get("cdisKS_ISR2");
  //
  TH2D *sca_vTcPR_Ceex2 = (TH2D*)DiskFileA.Get("sca_vTcPR_Ceex2");
  TH2D *sca_vTcPR_Eex2  = (TH2D*)DiskFileA.Get("sca_vTcPR_Eex2");
 //------------------------------------------------------------------------  
  //****************************************************************************************
  //************************** Developement corner *****************************************
  //TH1D *hst_ProjV = (TH1D*)sca_vTcPR_Ceex2->ProjectionX("hst_projV",1,50,"e");
  TH1D *hst_projC;
  ProjY1(sca_vTcPR_Ceex2, hst_projC);
  hst_projC->SetName("hst_projC");
  //****************************************************************************************
  //****************************************************************************************
   ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigCtest = new TCanvas("cFigCtest","FigCtest: cos(thet) dis.", 30, 70,    1000, 800);
  //                                Name    Title            xoff,yoff, WidPix,HeiPix
  cFigCtest->SetFillColor(10);
  ////////////////////////////////////////////////////////////////////////////////
  cFigCtest->Divide( 2,  2);
  //cFigCtest->Divide( 2,  2,     0.0,     0.0,   10);
  //              nx, ny, xmargin, ymargin, color
  //////////////////////////////////////////////
  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.03);
  //==========plot1==============
  cFigCtest->cd(1);
  hst_CosPLCeex2->SetStats(0);
  hst_CosPLCeex2->SetTitle(0);
  hst_CosPLCeex2->DrawCopy("h");  // black
  //
  hst_CosPRCeex2->SetLineColor(2); // red
  hst_CosPRCeex2->DrawCopy("hsame");
  //
  hst_Cost1Ceex2->SetLineColor(4); // blue
  hst_Cost1Ceex2->DrawCopy("hsame");
  //
  CaptT->DrawLatex(0.02,0.95,"d#sigma/dc (Ceex2); Black=PL, Red=PRD, Blue=Cth1");
  //==========plot2==============
  cFigCtest->cd(2);
  hst_CosPRCeex2->SetLineColor(2); // red
  hst_CosPRCeex2->DrawCopy("h");
  //
  hst_projC->DrawCopy("same");
  //
  CaptT->DrawLatex(0.02,0.95,"d#sigma/dc (Ceex2); Red PRD");
  //==========plot3==============
  cFigCtest->cd(3);
  hst_CosPRCeex2->SetStats(0);
  hst_CosPRCeex2->SetTitle(0);
  hst_CosPRCeex2->SetLineColor(2); // red
  hst_CosPRCeex2->DrawCopy("h");
  //
  cdisKS_ISR2->DrawCopy("hsame");  // black
  //
  CaptT->DrawLatex(0.02,0.95,"d#sigma/dc (Ceex2); Red=PRD, Black=KKsem ");
  //==========plot4==============
  cFigCtest->cd(4);
  hst_CosPREex2->DrawCopy("h");
  //----------------------------
  cFigCtest->cd();
  //================================================
}

///////////////////////////////////////////////////////////////////////////////////
void FigVprod()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigVprod =========================== "<<endl;
  //
  TH2D *sca_vTcPR_Ceex2 = (TH2D*)DiskFileA.Get("sca_vTcPR_Ceex2");
  TH2D *sca_vTcPR_Eex2  = (TH2D*)DiskFileA.Get("sca_vTcPR_Eex2");
  // KKsem
  TH1D *vcum_ISR2_FSR2   = (TH1D*)DiskFileB.Get("vcum_ISR2_FSR2");
  TH1D *vcum2_ISR2_FSR2  = (TH1D*)DiskFileB.Get("vcum2_ISR2_FSR2"); // with ctheta cut
  //****************************************************************************************
  int nbMax=0; // this is cut on costheta, 0= no cut
  //nbMax=22;    // cosThetaMax = 22/25
  TH1D                    *HTot_vTcPR_Ceex2, *HAfb_vTcPR_Ceex2;
  ProjV( sca_vTcPR_Ceex2,  HTot_vTcPR_Ceex2,  HAfb_vTcPR_Ceex2, nbMax);
  HTot_vTcPR_Ceex2->SetName("HTot_vTcPR_Ceex2");
  HAfb_vTcPR_Ceex2->SetName("HAfb_vTcPR_Ceex2");
  nbMax=22;    // cosThetaMax = 22/25
  TH1D                    *HTot2_vTcPR_Ceex2, *HAfb2_vTcPR_Ceex2;
  ProjV( sca_vTcPR_Ceex2,  HTot2_vTcPR_Ceex2,  HAfb2_vTcPR_Ceex2, nbMax);
  HTot2_vTcPR_Ceex2->SetName("HTot2_vTcPR_Ceex2");
  HAfb2_vTcPR_Ceex2->SetName("HAfb2_vTcPR_Ceex2");
  //
  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigVprod = new TCanvas("cFigVprod","FigVprod: photonic2", 70, 20,    1000, 800);
  //                            Name    Title               xoff,yoff, WidPix,HeiPix
  cFigVprod->SetFillColor(10);
  ////////////////////////////////////////////////////////////////////////////////
  cFigVprod->Divide( 2,  2);
  //cFigVprod->Divide( 2,  2,     0.0,     0.0,   10);
  //              nx, ny, xmargin, ymargin, color
  //////////////////////////////////////////////
  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.04);
  TLatex *CaptTb = new TLatex(0.40,0.01,"v_{max}");
  CaptTb->SetNDC(); // !!!
  CaptTb->SetTextSize(0.04);
  //==========plot1==============
  cFigVprod->cd(1);
  //gPad->SetLogy(); // !!!!!!
  HTot_vTcPR_Ceex2->SetStats(0);
  HTot_vTcPR_Ceex2->SetTitle(0);
  HTot_vTcPR_Ceex2->SetLineColor(4); // blue
  HTot_vTcPR_Ceex2->DrawCopy(" ");
  HTot2_vTcPR_Ceex2->SetLineColor(2); // red
  HTot2_vTcPR_Ceex2->DrawCopy("same");
  // KKsem
  vcum_ISR2_FSR2->DrawCopy("hsame");
  vcum2_ISR2_FSR2->DrawCopy("hsame");
  //
  CaptT->DrawLatex(0.02,0.95,"(a) Ceex2: #sigma(v_{max}), Blue=NoCut, Red=Ccut, Black=KKsem ");
  CaptTb->Draw();
  //==========plot2==============
  cFigVprod->cd(2);
  HTot_vTcPR_Ceex2->Divide(vcum_ISR2_FSR2);
  HTot_vTcPR_Ceex2->SetMinimum(0.99);
  HTot_vTcPR_Ceex2->SetMaximum(1.01);
  HTot_vTcPR_Ceex2->DrawCopy(" ");
  HTot2_vTcPR_Ceex2->Divide(vcum2_ISR2_FSR2);
  HTot2_vTcPR_Ceex2->DrawCopy("same");
  CaptT->DrawLatex(0.02,0.95,"(b) Ceex2: #sigma(v_{max}), MC/KKsem, Blue=NoCut, Red=Ccut ");
  CaptTb->Draw();
  //==========plot3==============
  cFigVprod->cd(3);
  HAfb_vTcPR_Ceex2->SetStats(0);
  HAfb_vTcPR_Ceex2->SetTitle(0);
  HAfb_vTcPR_Ceex2->SetLineColor(4); // blue
  HAfb_vTcPR_Ceex2->DrawCopy(" ");
  HAfb2_vTcPR_Ceex2->SetLineColor(2); // red
  HAfb2_vTcPR_Ceex2->DrawCopy("same");
  CaptT->DrawLatex(0.02,0.95,"(c) Ceex2: A_{FB}(v_{max}), Blue=NoCut, Red=22/25"); 
  CaptTb->Draw();
  //----------------------------
  cFigVprod->cd();
  //================================================
}//FigVprod


///////////////////////////////////////////////////////////////////////////////////
void FigCprod()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigCprod =========================== "<<endl;
  // from KKsem
  TH1D *cdisKS_ISR2      = (TH1D*)DiskFileB.Get("cdisKS_ISR2");
  TH1D *cdisDZ_ISR2      = (TH1D*)DiskFileB.Get("cdisDZ_ISR2");
  TH1D *casyKS_ISR2     = (TH1D*)DiskFileB.Get("casyKS_ISR2");
  TH1D *casyDZ_ISR2     = (TH1D*)DiskFileB.Get("casyDZ_ISR2");
  TH1D *cdis_ISR2_Zoff       = (TH1D*)DiskFileB.Get("cdis_ISR2_Zoff"); // no Z, pure QED
//
  TH1D *Hcth_vTcPR_Ceex2_vmax90 = (TH1D*)DiskFileB.Get("Hcth_vTcPR_Ceex2_vmax90");
  TH1D *Hcas_vTcPR_Ceex2_vmax90 = (TH1D*)DiskFileB.Get("Hcas_vTcPR_Ceex2_vmax90");
  TH1D *Hcth_vTcPR_Ceex2_vmax30= (TH1D*)DiskFileB.Get("Hcth_vTcPR_Ceex2_vmax30");
  TH1D *Hcth_vTcPR_Ceex2_vmax10= (TH1D*)DiskFileB.Get("Hcth_vTcPR_Ceex2_vmax10");
  TH1D *Hcth_vTcPR_Ceex2_vmax02= (TH1D*)DiskFileB.Get("Hcth_vTcPR_Ceex2_vmax02");
  TH1D *Hcas_vTcPR_Ceex2_vmax02= (TH1D*)DiskFileB.Get("Hcas_vTcPR_Ceex2_vmax02");
  TH1D *Hcas_vTcPR_Ceex2n_vmax02=(TH1D*)DiskFileB.Get("Hcas_vTcPR_Ceex2n_vmax02");
  TH1D *Hcas_vTcPR_Ceex2_vmax10= (TH1D*)DiskFileB.Get("Hcas_vTcPR_Ceex2_vmax10");
  TH1D *Hcas_vTcPR_Ceex2_vmax30= (TH1D*)DiskFileB.Get("Hcas_vTcPR_Ceex2_vmax30");


  //****************************************************************************************
  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigCprod = new TCanvas("cFigCprod","FigCprod: cos theta production", 50, 50,    1000, 800);
  //                            Name    Title               xoff,yoff, WidPix,HeiPix
  cFigCprod->SetFillColor(10);
  ////////////////////////////////////////////////////////////////////////////////
  cFigCprod->Divide( 2,  2);
  //cFigCprod->Divide( 2,  2,     0.0,     0.0,   10);
  //              nx, ny, xmargin, ymargin, color
  //////////////////////////////////////////////
  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.04);
  //==========plot1==============
  cFigCprod->cd(1);
  TH1D *Hst;
  Hst = cdisKS_ISR2;
  Hst->SetStats(0);
  Hst->SetTitle(0);
  Hst->SetMinimum(0);
  Hst->GetYaxis()->CenterTitle();
  Hst->GetYaxis()->SetTitleSize(0.05);
  Hst->GetYaxis()->SetTitle("d#sigma/dcos(#theta)");
  Hst->GetXaxis()->SetTitleSize(0.04);
  Hst->GetXaxis()->SetTitle("cos(#theta)");
  //
  Hst->SetLineColor(kMagenta);
  Hst->DrawCopy("h");
 // Dizet??
  cdisDZ_ISR2->SetLineColor(kGreen);
  cdisDZ_ISR2->DrawCopy("same");
  //
  Hcth_vTcPR_Ceex2_vmax02->SetLineColor(kBlue);
  Hcth_vTcPR_Ceex2_vmax02->DrawCopy("hsame");
  //
  CaptT->DrawLatex(0.02,0.96,"(a) KKsem EW-on/off Green/Magenta, v=<0.02, ISR only");
  CaptT->DrawLatex(0.02,0.91,"     CEEX2 v=<0.02, blue");
  //==========plot2==============
  cFigCprod->cd(2);
  TH1D *hZero = (TH1D*)cdisKS_ISR2->Clone("hZero"); // zero line
  hZero->Reset();
  Hst = casyKS_ISR2;
  Hst->SetStats(0);
  Hst->SetTitle(0);
  Hst->GetYaxis()->CenterTitle();
  Hst->GetYaxis()->SetTitleSize(0.05);
  Hst->GetYaxis()->SetTitle("AFB(#theta)");
  Hst->GetXaxis()->SetTitleSize(0.04);
  Hst->GetXaxis()->SetTitle("cos(#theta)");
  //
  Hst->SetLineColor(kMagenta); // magenta
  Hst->DrawCopy("h");
  // no Z
  //cdis4_ISR2->SetLineColor(2); // red no Z
  //cdis4_ISR2->DrawCopy("hsame");
  // equivalent operation on histograms
  //cdis_ISR2->Add(cdis_ISR2_Zoff,-1.0); // subtract pure QED
  //cdis_ISR2->Divide(cdis_ISR2_Zoff);   // divide over pure QED
  //cdis_ISR2->SetLineColor(8); // green
  //
  //cdisKS_ISR2->DrawCopy("same");
  //hZero->DrawCopy("hsame");
  // Dizet???
  casyDZ_ISR2->SetLineColor(8); // green
  casyDZ_ISR2->DrawCopy("same");
  //
  Hcas_vTcPR_Ceex2_vmax02->SetLineColor(kBlue);
  Hcas_vTcPR_Ceex2_vmax02->DrawCopy("hsame");
  //
  Hcas_vTcPR_Ceex2n_vmax02->SetLineColor(kCyan);
  Hcas_vTcPR_Ceex2n_vmax02->DrawCopy("hsame");
  //
  hZero->DrawCopy("hsame");
  //
  CaptT->DrawLatex(0.02,0.96,"(b) KKsem, v<0.02, EW on/off=green/magenta ,ISR only, ");
  CaptT->DrawLatex(0.02,0.91,"      CEEX2, v<0.02, IFI on/off= blue/cyan");
  //==========plot3==============
  cFigCprod->cd(3);
  Hst=Hcth_vTcPR_Ceex2_vmax90;
  Hst->SetStats(0);
  Hst->SetTitle(0);
  Hst->GetYaxis()->CenterTitle();
  Hst->GetYaxis()->SetTitleSize(0.05);
  Hst->GetYaxis()->SetTitle("d#sigma/dcos(#theta)");
  Hst->GetXaxis()->SetTitleSize(0.04);
  Hst->GetXaxis()->SetTitle("cos(#theta)");
  //
  Hst->SetLineColor(kBlack);
  Hst->SetMinimum(0);
  Hst->DrawCopy("h");
  //
  Hcth_vTcPR_Ceex2_vmax30->SetLineColor(kRed);
  Hcth_vTcPR_Ceex2_vmax30->DrawCopy("hsame");
  //
  Hcth_vTcPR_Ceex2_vmax10->SetLineColor(kGreen);
  Hcth_vTcPR_Ceex2_vmax10->DrawCopy("hsame");
  //
  Hcth_vTcPR_Ceex2_vmax02->SetLineColor(kBlue);
  Hcth_vTcPR_Ceex2_vmax02->DrawCopy("hsame");
  //cdisKS_ISR2->DrawCopy("hsame");
  CaptT->DrawLatex(0.02,0.95,"(c) Ceex2: Black/Red/Green/Blue= v< 0.90, 0.30, 0.10, 0.02");
  //==========plot4==============
  cFigCprod->cd(4);
  Hst=Hcas_vTcPR_Ceex2_vmax02;
  Hst->SetStats(0);
  Hst->SetTitle(0);
  Hst->GetYaxis()->CenterTitle();
  Hst->GetYaxis()->SetTitleSize(0.05);
  Hst->GetYaxis()->SetTitle("AFB(#theta)");
  Hst->GetXaxis()->SetTitleSize(0.04);
  Hst->GetXaxis()->SetTitle("cos(#theta)");
  Hst->SetLineColor(kBlue);
  Hst->DrawCopy("h");
  // IFI off
  Hcas_vTcPR_Ceex2n_vmax02->SetLineColor(kCyan);
  Hcas_vTcPR_Ceex2n_vmax02->DrawCopy("hsame");
  //
  Hcas_vTcPR_Ceex2_vmax10->SetLineColor(kGreen);
  Hcas_vTcPR_Ceex2_vmax10->DrawCopy("hsame");
  //
  Hcas_vTcPR_Ceex2_vmax30->SetLineColor(kRed);
  Hcas_vTcPR_Ceex2_vmax30->DrawCopy("hsame");
  //
  Hcas_vTcPR_Ceex2_vmax90->SetLineColor(kBlack);
  Hcas_vTcPR_Ceex2_vmax90->DrawCopy("hsame");
 //
  hZero->DrawCopy("hsame");
  CaptT->DrawLatex(0.02,0.96,"(d) CEEX2: Black/Red/Green/Blue =  v< 0.9, 0.3, 0.10, 0.02");
  CaptT->DrawLatex(0.02,0.91,"    CEEX2: IFI off, Cyan v< 0.02");
  //----------------------------
  cFigCprod->cd();
  //================================================
}//FigCprod



///////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  //++++++++++++++++++++++++++++++++++++++++
  TApplication theApp("theApp", &argc, argv);
  //++++++++++++++++++++++++++++++++++++++++
  HistNormalize();     // Renormalization of MC histograms
  KKsemMakeHisto();    // prepave histos for plotting
  //========== PLOTTING ==========
  FigScatA();
  FigInfo();
  FigVtest();
  FigCtest();
  FigVprod();
  FigCprod();
  //++++++++++++++++++++++++++++++++++++++++
  DiskFileA.ls();
  DiskFileB.ls();
  DiskFileB.Write();
  DiskFileB.Close();
  //++++++++++++++++++++++++++++++++++++++++
  theApp.Run();
  //++++++++++++++++++++++++++++++++++++++++
}
