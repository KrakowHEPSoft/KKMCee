//////////////////////////////////////////////////////////////////////
//    make Plot1
//////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
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
#include "TObjString.h"
#include "TFile.h"

#include "HisNorm.h"
#include "KKcol.h"

//=============================================================================
//  ROOT  ROOT ROOT   ROOT  ROOT  ROOT  ROOT  ROOT  ROOT  ROOT   ROOT   ROOT
//=============================================================================
// current
//TFile DiskFileA("../test0/rmain.root");
TFile DiskFileA("../workAFB/rmain.root");
TFile DiskFileB("RhoSemi.root","RECREATE","histograms");
FILE *DiskFileT;

// Interface to KKcol and some extra plotting facilities
KKcol LibSem;

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
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sca_vTcPR_Ceex2n") );
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sca_vTcPR_Eex2") );
  //
}



///////////////////////////////////////////////////////////////////////////////////
void KKsemMakeHisto(){
  // Here we produce semianalytical plots using KKsem program, No plotting
  //------------------------------------------------------------------------
  cout<<"==================================================================="<<endl;
  cout<<"================ KKsem MakeHisto  BEGIN ============================"<<endl;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  double CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR

  // initialization of KKsem
  //KKcol LibSem;
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

  TH1D *hstVtemplate = (TH1D*)DiskFileA.Get("hst_vTrueCeex2");
  TH1D *hstCtemplate = (TH1D*)DiskFileA.Get("hst_Cost1Ceex2");
//------------------------------------------------------------------------
//   MuMu  Sigma(vmax) with ulimited c=cos(theta)
//------------------------------------------------------------------------
// ISR*FSR
  KeyDis = 302302;        // ISR*FSR O(alf2)
  sprintf(chak,"XCHI2");  // ISR*FSR Mff
  TH1D *vcum_ISR2_FSR2 =(TH1D*)hstVtemplate->Clone("vcum_ISR2_FSR2");
  LibSem.VVplot(vcum_ISR2_FSR2, KF, chak, KeyDis, KeyFob);
//-------------------------------------------------
//    AFB(vmax) for unlimited c=cos(theta)
//-------------------------------------------------
  kksem_setcrange_(0, 25.0/25); // forward, unlimited cos(theta)
  TH1D *afbv_ISR2_FSR2 =(TH1D*)hstVtemplate->Clone("afbv_ISR2_FSR2");
  LibSem.VVplot(afbv_ISR2_FSR2, KF, chak, KeyDis, KeyFob);// Forward
  afbv_ISR2_FSR2->Add(afbv_ISR2_FSR2, vcum_ISR2_FSR2, 2.0, -1.0) ; // numerator F-B = 2F-(F+B)
  afbv_ISR2_FSR2->Divide(vcum_ISR2_FSR2);                          // finally (F-B)(F+B)
  //
  cout<<"================ KKsem MakeHisto ENDs ============================="<<endl;
  cout<<"==================================================================="<<endl;
//------------------------------------------------------------------------
//------------------------------------------------------------------------
}//  KKsemMakeHisto

///////////////////////////////////////////////////////////////////////////////////
void ReMakeMChisto(){
	// Some KKMC histos are preprocessed
	//------------------------------------------------------------------------
  cout<<"==================================================================="<<endl;
  cout<<"================ ReMakeMChisto  BEGIN  ============================"<<endl;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  double CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR

  //****************************************************************************************
  // Pure MC reprocessing part
  //
  TH2D *sca_vTcPR_Eex2  = (TH2D*)DiskFileA.Get("sca_vTcPR_Eex2");
  TH2D *sca_vTcPR_Ceex2 = (TH2D*)DiskFileA.Get("sca_vTcPR_Ceex2");
  TH2D *sca_vTcPR_Ceex2n = (TH2D*)DiskFileA.Get("sca_vTcPR_Ceex2n");

  ///****************************************************************************************
  /// Distributions of v=vTrue with unlimited c=cos(theta)
  //  without cutoff on c=cos(thetaPRD)
  int nbMax=0;   // cosThetaMax = 1.0
  TH1D                    *HTot_vTcPR_Eex2, *HAfb_vTcPR_Eex2;
  ProjV( sca_vTcPR_Eex2,  HTot_vTcPR_Eex2,  HAfb_vTcPR_Eex2, nbMax);  //!!!!
  HTot_vTcPR_Eex2->SetName("HTot_vTcPR_Eex2");
  HAfb_vTcPR_Eex2->SetName("HAfb_vTcPR_Eex2");
  nbMax=0;   // cosThetaMax = 1.0
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
  //if( CMSene<91.0 ) HAfb_vTcPR_Ceex2n->Scale(-1);

  cout<<"================ ReMakeMChisto ENDs  ============================="<<endl;
  cout<<"==================================================================="<<endl;
}//RemakeMChisto


///////////////////////////////////////////////////////////////////////////////////
void FigOldBench()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigOldBench =========================== "<<endl;
  //
  Double_t CMSene;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  char TextEne[100]; sprintf(TextEne,"#sqrt{s} =%4.2fGeV", CMSene);
  //
  // KKsem
  TH1D *vcum_ISR2_FSR2   = (TH1D*)DiskFileB.Get("vcum_ISR2_FSR2");
  TH1D *afbv_ISR2_FSR2   = (TH1D*)DiskFileB.Get("afbv_ISR2_FSR2");

  // Distributions of v=vTrue
  // without cutoff on c=cos(thetaPRD)
  TH1D *HTot_vTcPR_Ceex2  = (TH1D*)DiskFileB.Get("HTot_vTcPR_Ceex2");
  TH1D *HAfb_vTcPR_Ceex2  = (TH1D*)DiskFileB.Get("HAfb_vTcPR_Ceex2");

  TH1D *HTot_vTcPR_Ceex2n= (TH1D*)DiskFileB.Get("HTot_vTcPR_Ceex2n");
  TH1D *HAfb_vTcPR_Ceex2n= (TH1D*)DiskFileB.Get("HAfb_vTcPR_Ceex2n");
  //
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.04);
  TLatex *CaptTb = new TLatex(0.40,0.01,"v_{max}");
  CaptTb->SetNDC(); // !!!
  CaptTb->SetTextSize(0.04);
  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigOldBench = new TCanvas("cFigOldBench","FigOldBench", 70, 20,    700, 700);
  //                                   Name    Title     xoff,yoff, WidPix,HeiPix
  cFigOldBench->SetFillColor(10);

  TH1D *HTot_rat_Ceex2 =(TH1D*)HTot_vTcPR_Ceex2->Clone("HTot_rat_Ceex2");
  HTot_rat_Ceex2->Divide(vcum_ISR2_FSR2);

  HTot_rat_Ceex2->SetMinimum(0.975);
  HTot_rat_Ceex2->SetMaximum(1.100);
  HTot_rat_Ceex2->SetStats(0);
  HTot_rat_Ceex2->SetTitle(0);
  HTot_rat_Ceex2->DrawCopy("h");
  //
  TH1D *HTot_rat_Ceex2n = (TH1D*)HTot_vTcPR_Ceex2n->Clone("HTot_rat_Ceex2n");
  HTot_rat_Ceex2n->Divide(vcum_ISR2_FSR2);

  HTot_rat_Ceex2n->SetLineColor(kMagenta);
  HTot_rat_Ceex2n->DrawCopy("hsame");

  CaptT->DrawLatex(0.12,0.85,"Ceex2/KKsem, Blue/Magenta for IFI on/off");
  CaptT->DrawLatex(0.60,0.75,TextEne);
  CaptTb->Draw();
  //----------------------------
  cFigOldBench->cd();
  cFigOldBench->SaveAs("cFigOldBench.pdf");

  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigOldBench2 = new TCanvas("cFigOldBench2","FigOldBench2",100, 50,    700, 700);
  //                                   Name    Title     xoff,yoff, WidPix,HeiPix
  cFigOldBench2->SetFillColor(10);

  TH1D *HAfb_diff_Ceex2 =(TH1D*)HAfb_vTcPR_Ceex2->Clone("HAfb_diff_Ceex2");
  HAfb_diff_Ceex2->Add(HAfb_diff_Ceex2,  afbv_ISR2_FSR2, 1.0, -1.0);

  HAfb_diff_Ceex2->SetLineColor(kRed); // red
  HAfb_diff_Ceex2->SetStats(0);
  HAfb_diff_Ceex2->SetTitle(0);
  HAfb_diff_Ceex2->SetMinimum(-0.02);
  HAfb_diff_Ceex2->SetMaximum( 0.06);
  HAfb_diff_Ceex2->DrawCopy("h");
  //
  TH1D *HAfb_diff_Ceex2n =(TH1D*)HAfb_vTcPR_Ceex2n->Clone("HAfb_diff_Ceex2n");
  HAfb_diff_Ceex2n->Add(HAfb_diff_Ceex2n, afbv_ISR2_FSR2, 1.0, -1.0) ;
  HAfb_diff_Ceex2n->SetLineColor(kGreen); // green
  HAfb_diff_Ceex2n->SetLineWidth(2);
  HAfb_diff_Ceex2n->DrawCopy("hsame");
  //
   //
  CaptT->DrawLatex(0.12,0.85,"A^{KKMC}_{FB}-A^{KKsem}_{FB}, Red/Green = IFI on/off");
  CaptT->DrawLatex(0.60,0.75,TextEne);

  //----------------------------
  cFigOldBench2->cd();
  cFigOldBench2->SaveAs("cFigOldBench2.pdf");
  //================================================
}//FigOldBench


///////////////////////////////////////////////////////////////////////////////////
void TabOldBench()
{
//------------------------------------------------------------------------
  cout<<" ========================= TabOldBench start=========================== "<<endl;
//
  Double_t CMSene;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  char TextEne[100]; sprintf(TextEne,"#sqrt{s} =%4.2fGeV", CMSene);
  //
  // KKsem
  TH1D *vcum_ISR2_FSR2   = (TH1D*)DiskFileB.Get("vcum_ISR2_FSR2");
  TH1D *afbv_ISR2_FSR2   = (TH1D*)DiskFileB.Get("afbv_ISR2_FSR2");

  // Distributions of v=vTrue
  // without cutoff on c=cos(thetaPRD)
  TH1D *HTot_vTcPR_Eex2  = (TH1D*)DiskFileB.Get("HTot_vTcPR_Eex2");
  TH1D *HAfb_vTcPR_Eex2  = (TH1D*)DiskFileB.Get("HAfb_vTcPR_Eex2");

  TH1D *HTot_vTcPR_Ceex2  = (TH1D*)DiskFileB.Get("HTot_vTcPR_Ceex2");
  TH1D *HAfb_vTcPR_Ceex2  = (TH1D*)DiskFileB.Get("HAfb_vTcPR_Ceex2");

  TH1D *HTot_vTcPR_Ceex2n= (TH1D*)DiskFileB.Get("HTot_vTcPR_Ceex2n");
  TH1D *HAfb_vTcPR_Ceex2n= (TH1D*)DiskFileB.Get("HAfb_vTcPR_Ceex2n");

//  Char_t Capt[20][132];

// Column captions
  int nPlt=4;   // KORALZ eliminated
  Char_t *Capt[nPlt+1];
  for( int i=0; i<=nPlt; i++ ) Capt[i]=new char[132];
  strcpy(Capt[0],"{\\color{blue}$v_{\\max}$}");
  strcpy(Capt[1],"{\\color{blue} ${\\cal KK}$sem Refer.}");
  strcpy(Capt[2],"{\\color{blue}${\\cal O}(\\alpha^3)_{\\rm EEX3}$ }");
  strcpy(Capt[3],"{\\color{red}${\\cal O}(\\alpha^2)_{\\rm CEEX}$ intOFF}");
  strcpy(Capt[4],"{\\color{red}${\\cal O}(\\alpha^2)_{\\rm CEEX}$ }");

// formats, not used in PlTable2
//  Char_t fmt[3][10];
//  strcpy(fmt[0],"f10.2");
//  strcpy(fmt[1],"f10.4");
//  strcpy(fmt[2],"f8.4");

// pointers to histograms
  TH1D *iHst[nPlt+1];
  iHst[1]= vcum_ISR2_FSR2;     //KKsem
  iHst[2]= HTot_vTcPR_Eex2;    //EEX2
  iHst[3]= HTot_vTcPR_Ceex2n;  //CEEX2 INT off
  iHst[4]= HTot_vTcPR_Ceex2;   //CEEX2
  iHst[1]->Scale(1e3);    // nano- to pico-barns
  iHst[2]->Scale(1e3);    // nano- to pico-barns
  iHst[3]->Scale(1e3);    // nano- to pico-barns
  iHst[4]->Scale(1e3);    // nano- to pico-barns
// multicolumn caption
  Char_t Mcapt[132];
  strcpy(Mcapt,"{\\color{red}$\\sigma(v_{\\max})$ [pb]}");

///************************************
  DiskFileT = fopen("TabOldBench.txp","w");
//************************************
// Initialization of the latex source file
  LibSem.PlInitialize(DiskFileT, 2);


//  int k1,k2,dk;
//  k1=10; k2=90; dk=20;  //original
//  k1= 5; k2=45; dk=10;
  LibSem.PlTable2( nPlt, iHst, DiskFileT, Capt,  Mcapt, "B", 1, 1, 1); // for 50 bins
  LibSem.PlTable2(-nPlt, iHst, DiskFileT, Capt,  Mcapt, "T", 5,45,10); // for 50 bins
  LibSem.PlTable2(-nPlt, iHst, DiskFileT, Capt,  Mcapt, "T",50,50, 1); // for 50 bins

  iHst[1]= afbv_ISR2_FSR2;     //KKsem
  iHst[2]= HAfb_vTcPR_Eex2;    //EEX2
  iHst[3]= HAfb_vTcPR_Ceex2n;  //CEEX2 INT off
  iHst[4]= HAfb_vTcPR_Ceex2;   //CEEX2

  strcpy(Mcapt,"{\\color{red}$A_{\\rm FB}(v_{\\max})$}");
  LibSem.PlTable2( nPlt, iHst, DiskFileT, Capt,  Mcapt, "T", 1, 1, 1); // for 50 bins
  LibSem.PlTable2(-nPlt, iHst, DiskFileT, Capt,  Mcapt, "T", 5,45,10); // for 50 bins
  LibSem.PlTable2(-nPlt, iHst, DiskFileT, Capt,  Mcapt, "E",50,50, 1); // for 50 bins

// finalizing latex source file
  LibSem.PlEnd(DiskFileT);
//************************************
  fclose(DiskFileT);
//************************************
  cout<<" ========================= TabOldBench end =========================== "<<endl;
}//TabOldBench


////////////////////////////////////////////////////////////////////////////////
//
//  Monte Carlo integration/run using Foam
//
////////////////////////////////////////////////////////////////////////////////
//______________________________________________________________________________
void ISRgener()
{
  cout<<"--- ISRgener started ---"<<endl;

  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  double CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  double vvmax   = HST_KKMC_NORMA->GetBinContent(17);

  // %%% Histogramming %%%
  DiskFileB.cd();
  TH1D *hst_Vtemplate = (TH1D*)DiskFileA.Get("hst_vTrueCeex2");
  TH1D *hst_Ctemplate = (TH1D*)DiskFileA.Get("hst_Cost1Ceex2");
  TH2D *sca_VCtemplate= (TH2D*)DiskFileA.Get("sca_vTcPR_Ceex2");

  TH1D *HST_vv_Ceex2n = (TH1D*)hst_Vtemplate->Clone("HST_vv_Ceex2n");
  HST_vv_Ceex2n->Reset();

  // %%% Configuring integrand for Foam %%%
  Rho4Foam *Rho1= new Rho4Foam("Rho4Foam");
  Rho1->m_CMSene = CMSene;
  Rho1->m_vvmax  = vvmax;
  cout<<"%%%%%   ISRgener:  vvmax= "<< vvmax << " from KKMC" <<endl;
  //------------------------------------------
  TRandom  *PseRan   = new TRandom3();  // Create random number generator
  PseRan->SetSeed(4357);
  // %%% FOAM simulator/integrator %%%
  //------------------------------------------
  TFoam   *MC_fisr    = new TFoam("MC_fisr");   // Create Simulator
  MC_fisr->SetkDim(2);         // No. of dimensions, obligatory!
  MC_fisr->SetnCells( 10000);  // No. of cells, can be omitted, default=2000
  MC_fisr->SetnSampl(100000);  // No. of MC evts/cell in exploration, default=200
  MC_fisr->SetRho(Rho1);       // Set 2-dim distribution, included above
  MC_fisr->SetPseRan(PseRan);  // Set random number generator, mandatory!
  MC_fisr->SetOptRej(0);       // wted events (=0), default wt=1 events (=1)
  MC_fisr->Initialize();       // Initialize simulator, may take time...

  // %%% loop over MC events %%%
  double wt,Mll, wt2, Mka, vv;
  long NevTot = 2000000;  // 2M
  NevTot      = 8000000;  // 8M
  //NevTot =    100000000;  // 100M
  for(long loop=0; loop<NevTot; loop++)
  {
	//----------------------------------------------------------
	/// Generate ISR+FSR event
    MC_fisr->MakeEvent();            // generate MC event
    MC_fisr->GetMCwt(wt2);
    Mll = Rho1->m_Mll;
    vv = 1-sqr(Mll/CMSene);
    HST_vv_Ceex2n->Fill(vv,wt2);
    if( 1000000*(loop/1000000) == loop) cout<<" Nev ="<< loop<< endl;
  }// loop
  //  Renormalizing histograms
  double Xsav, dXsav;
  MC_fisr->GetIntNorm(Xsav,dXsav);
  HisNorm0( NevTot, Xsav, HST_vv_Ceex2n);
  //
  TH1D *HST_vmax_Ceex2n;
  MakeCumul(HST_vv_Ceex2n,HST_vmax_Ceex2n);
  HST_vmax_Ceex2n->SetName("HST_vmax_Ceex2n");
//
  cout<<"--- ISRgener ended ---"<<endl;
}//ISRgener


///////////////////////////////////////////////////////////////////////////////////
void FigVtest()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigVtest =========================== "<<endl;
  Double_t CMSene;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  //
  TH1D *hst_vTrueCeex2 = (TH1D*)DiskFileB.Get("HTot_vTcPR_Ceex2n");
  TH1D *hst_vXGenCeex2 = (TH1D*)DiskFileA.Get("hst_vXGenCeex2");
  //
  TH1D *HST_vv_Ceex2n        = (TH1D*)DiskFileB.Get("HST_vv_Ceex2n");
  TH1D *HST_vmax_Ceex2n      = (TH1D*)DiskFileB.Get("HST_vmax_Ceex2n");
  //
  //TH1D *vdis_ISR2      = (TH1D*)DiskFileB.Get("vdis_ISR2");
  //TH1D *vdis_ISR2_FSR2 = (TH1D*)DiskFileB.Get("vdis_ISR2_FSR2");
  //TH1D *Hpro_vT_Ceex2  = (TH1D*)DiskFileB.Get("Hpro_vT_Ceex2");
  //
  //*****************************************************************************
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
  CaptT->SetTextSize(0.04);
  //==========plot1==============
  cFigVtest->cd(1);
  gPad->SetLogy(); // !!!!!!
  // MC v-true direct
  hst_vTrueCeex2->SetStats(0);
  hst_vTrueCeex2->SetTitle(0);
  hst_vTrueCeex2->SetMinimum(1e-4*hst_vTrueCeex2->GetMaximum());
  hst_vTrueCeex2->DrawCopy("h");  // black

  HST_vmax_Ceex2n->SetLineColor(kGreen); // green
  HST_vmax_Ceex2n->DrawCopy("hsame");
  // MC vtrue from scatergram
//  Hpro_vT_Ceex2->SetLineColor(kGreen); // green
//  Hpro_vT_Ceex2->DrawCopy("same");
  // KKsem ISR+FSR
//  vdis_ISR2_FSR2->SetLineColor(kMagenta); // magenta
//  vdis_ISR2_FSR2->DrawCopy("hsame");
  // KKsem ISR only
//  vdis_ISR2->SetLineColor(kBlue); // blue
//  vdis_ISR2->DrawCopy("hsame");
  //
  //
  CaptT->DrawLatex(0.02,0.95, "d#sigma/dv(ISR+FSR) KKMC CEEX2: Black, Red; KKsem=Magenta");
  CaptT->DrawLatex(0.02,0.91, "           ISR only, KKsem=Blue");
  //==========plot2==============
  cFigVtest->cd(2);
  hst_vTrueCeex2->Divide(HST_vmax_Ceex2n);
//  hst_vTrueCeex2->SetStats(0);
//  hst_vTrueCeex2->SetTitle(0);
  hst_vTrueCeex2->SetMinimum(0.85);
  hst_vTrueCeex2->SetMaximum(1.15);
  hst_vTrueCeex2->SetLineColor(kBlue);
  hst_vTrueCeex2->DrawCopy("h");
  //
  CaptT->DrawLatex(0.02,0.95,"d#sigma/dv(ISR+FSR); KKMC_CEEX2/KKsem");
  //==========plot3==============
  cFigVtest->cd(3);
  gPad->SetLogy(); // !!!!!!
  hst_vXGenCeex2->SetStats(0);
  hst_vXGenCeex2->SetTitle(0);
  hst_vXGenCeex2->SetLineColor(kRed); // red
  hst_vXGenCeex2->DrawCopy("h");
  //
//  vdis_ISR2->SetLineColor(kBlue); // blue
//  vdis_ISR2->DrawCopy("hsame");
  CaptT->DrawLatex(0.02,0.95,"d#sigma/dv(ISR),  KKMC_CEEX2=Red, Blue=KKsem");
  //==========plot4==============
  cFigVtest->cd(4);
//  hst_vXGenCeex2->Divide(vdis_ISR2);
//  hst_vXGenCeex2->SetStats(0);
//  hst_vXGenCeex2->SetTitle(0);
//  hst_vXGenCeex2->SetMinimum(0.85);
//  hst_vXGenCeex2->SetMaximum(1.15);
//  hst_vXGenCeex2->SetLineColor(kRed);
//  hst_vXGenCeex2->DrawCopy("h");  // black
  CaptT->DrawLatex(0.02,0.95,"d#sigma/dv(ISR); KKMC_CEEX2/KKsem");
  //----------------------------
  cFigVtest->cd();
  //================================================
}//FigVtest



///////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  //++++++++++++++++++++++++++++++++++++++++
  TApplication theApp("theApp", &argc, argv);
  //++++++++++++++++++++++++++++++++++++++++

  DiskFileB.cd();
  HistNormalize();     // Renormalization of MC histograms
  KKsemMakeHisto();    // prepare histos from KKsem
  ReMakeMChisto();     // reprocessing MC histos
  //========== PLOTTING ==========
  // New benchmarks KKMC vs. KKcol with Foam integrator
  ISRgener();
  FigVtest();
  // Old benchmarks KKMC vs. KKsem with Gauss integrator
  FigOldBench();
  TabOldBench();

  //++++++++++++++++++++++++++++++++++++++++
  DiskFileA.ls();
  DiskFileB.ls();
  DiskFileB.Write();
  DiskFileB.Close();

  //++++++++++++++++++++++++++++++++++++++++
  theApp.Run();
  //++++++++++++++++++++++++++++++++++++++++
}

