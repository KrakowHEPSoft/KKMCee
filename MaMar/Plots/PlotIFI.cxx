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

// current
//TFile DiskFileA("../test0/rmain.root");
TFile DiskFileA("../workAFB/rmain.root");
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

  TH1D *hstVtemplate = (TH1D*)DiskFileA.Get("hst_vTrueMain");
  TH1D *hstCtemplate = (TH1D*)DiskFileA.Get("hst_Cost1Ceex2");
//------------------------------------------------------------------------
//   MuMu  Sigma(vmax) with limited c=cos(theta)
//------------------------------------------------------------------------
// ISR*FSR
  KeyDis = 302302;        // ISR*FSR O(alf2)
  sprintf(chak,"XCHI2");  // ISR*FSR Mff
  TH1D *vcum_ISR2_FSR2 =(TH1D*)hstVtemplate->Clone("vcum_ISR2_FSR2");
  LibSem.VVplot(vcum_ISR2_FSR2, KF, chak, KeyDis, KeyFob);

//-------------------------------------------------
// and finally AFB(vmax) for ulimited c=cos(theta)
  kksem_setcrange_(0, 25.0/25); // forward
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
	// Here we produce semianalytical plots using KKsem program, No plotting
	// also some MC histos are preprocessed
	//------------------------------------------------------------------------
  cout<<"==================================================================="<<endl;
  cout<<"================ ReMakeMChisto  BEGIN  ============================"<<endl;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  double CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR

  //****************************************************************************************
  // Pure MC reprocessing part
  //
  TH2D *sca_vTcPR_Ceex2 = (TH2D*)DiskFileA.Get("sca_vTcPR_Ceex2");
  TH2D *sca_vTcPR_Eex2  = (TH2D*)DiskFileA.Get("sca_vTcPR_Eex2");
  TH2D *sca_vTcPR_Ceex2n = (TH2D*)DiskFileA.Get("sca_vTcPR_Ceex2n");

  ///****************************************************************************************
  /// Distributions of v=vTrue with limited c=cos(theta)
  //  without cutoff on c=cos(thetaPRD)
  int nbMax=0;   // cosThetaMax = 1.0
  TH1D                    *HTot_vTcPR_Ceex2, *HAfb_vTcPR_Ceex2;
  ProjV( sca_vTcPR_Ceex2,  HTot_vTcPR_Ceex2,  HAfb_vTcPR_Ceex2, nbMax);  //!!!!
  HTot_vTcPR_Ceex2->SetName("HTot_vTcPR_Ceex2");
  HAfb_vTcPR_Ceex2->SetName("HAfb_vTcPR_Ceex2");
  //if( CMSene<91.0 ) HAfb_vTcPR_Ceex2->Scale(-1);
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
void FigVprod()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigVprod =========================== "<<endl;
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
  TCanvas *cFigVprod = new TCanvas("cFigVprod","FigVprod", 70, 20,    700, 700);
  //                                   Name    Title     xoff,yoff, WidPix,HeiPix
  cFigVprod->SetFillColor(10);

  HTot_vTcPR_Ceex2->Divide(vcum_ISR2_FSR2);
  HTot_vTcPR_Ceex2->SetMinimum(0.975);
  HTot_vTcPR_Ceex2->SetMaximum(1.100);
  HTot_vTcPR_Ceex2->DrawCopy("h");
  //
  HTot_vTcPR_Ceex2n->SetLineColor(kMagenta);
  HTot_vTcPR_Ceex2n->Divide(vcum_ISR2_FSR2);
  HTot_vTcPR_Ceex2n->DrawCopy("hsame");

  CaptT->DrawLatex(0.12,0.95,"(b) Ceex2/KKsem, IFIon, Blue=|cos(#theta|<1, Red=|cos(#theta|<0.88 ");
  CaptT->DrawLatex(0.16,0.85,"    Magenta: Ceex2/KKsem for IFIoff, |cos(#theta|<1 ");
  CaptT->DrawLatex(0.60,0.75,TextEne);
  CaptTb->Draw();
  //----------------------------
  cFigVprod->cd();
  cFigVprod->SaveAs("cFigVprod.jpg");

  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigVprod2 = new TCanvas("cFigVprod2","FigVprod2",100, 50,    700, 700);
  //                                   Name    Title     xoff,yoff, WidPix,HeiPix
  cFigVprod2->SetFillColor(10);

  HAfb_vTcPR_Ceex2->Add(HAfb_vTcPR_Ceex2,  afbv_ISR2_FSR2,1.0, -1.0);
  HAfb_vTcPR_Ceex2->SetLineColor(kRed); // red
  HAfb_vTcPR_Ceex2->SetStats(0);
  HAfb_vTcPR_Ceex2->SetMinimum(-0.02);
  HAfb_vTcPR_Ceex2->SetMaximum( 0.06);
  HAfb_vTcPR_Ceex2->DrawCopy("h");
  //
  HAfb_vTcPR_Ceex2n->Add(HAfb_vTcPR_Ceex2n,afbv_ISR2_FSR2,1.0, -1.0) ;
  HAfb_vTcPR_Ceex2n->SetLineColor(kGreen); // green
  HAfb_vTcPR_Ceex2n->SetLineWidth(2);
  HAfb_vTcPR_Ceex2n->DrawCopy("hsame");
  //
  //afbv_ISR2_FSR2->SetLineColor(kBlack);
  //afbv_ISR2_FSR2->DrawCopy("hsame");
  //
  CaptT->DrawLatex(0.12,0.95,"(d) A^{KKMC}_{FB}(|cos(#theta|<0.88), Red/Green=IFIon/off");
  CaptT->DrawLatex(0.17,0.85,"  A^{KKsem}_{FB}(|cos(#theta|<0.88), Black=IFIoff");
  //----------------------------
  cFigVprod2->cd();
  cFigVprod2->SaveAs("cFigVprod2.jpg");
  //================================================
}//FigVprod

///////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  //++++++++++++++++++++++++++++++++++++++++
  TApplication theApp("theApp", &argc, argv);
  //++++++++++++++++++++++++++++++++++++++++
  HistNormalize();     // Renormalization of MC histograms
  KKsemMakeHisto();    // prepare histos from KKsem
  ReMakeMChisto();     // reprocessing MC histos
  //========== PLOTTING ==========

  FigVprod();

  //++++++++++++++++++++++++++++++++++++++++
  DiskFileA.ls();
  DiskFileB.ls();
  DiskFileB.Write();
  DiskFileB.Close();
  //++++++++++++++++++++++++++++++++++++++++
  theApp.Run();
  //++++++++++++++++++++++++++++++++++++++++
}
