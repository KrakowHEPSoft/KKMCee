//////////////////////////////////////////////////////////////////////
//    make Plot1-run
//////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <iomanip>
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
#include "TFile.h"

#include "HisNorm.h"
#include "KKplot.h"

//=============================================================================
//  ROOT  ROOT ROOT   ROOT  ROOT  ROOT  ROOT  ROOT  ROOT  ROOT   ROOT   ROOT 
//=============================================================================
// current
//TFile DiskFileA("../workKKMC/histo.root");

TFile DiskFileA("../workKKMC/histo.root_95GeV_26G");
//TFile DiskFileA("../workKKMC/histo.root_95GeV.4G");
//TFile DiskFileA("../workKKMC/histo.root_10GeV_5.7G"); //

//TFile DiskFileA("../workAFB/rmain.root");
// archive
//TFile DiskFileA("../workAFB/rmain.root_95GeV_100M");
//TFile DiskFileA("../workAFB/rmain.root_88GeV_100M"); // archive
//TFile DiskFileA("../workAFB/rmain.root_10GeV_30M");
//
TFile DiskFileB("RhoSemi.root","RECREATE","histograms");
//=============================================================================

KKplot LibSem("KKplot");

///////////////////////////////////////////////////////////////////////////////////
//              GLOBAL stuff
///////////////////////////////////////////////////////////////////////////////////
double gCMSene, gNevTot, gNevTot2; // from KKMC and KKfoam MC runs (histograms)
char   gTextEne[100], gTextNev[100], gTextNev2[100];
int    kGold=92, kBrune=46, kPine=71;
//
//
float  gXcanv = 10, gYcanv = 10;

///////////////////////////////////////////////////////////////////////////////////
void PlotSame2(TH1D *HST, double &ycapt, Int_t kolor, double xx,  TString label,  TString opis)
{
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.035);
  HST->SetLineColor(kolor);
  HST->DrawCopy("hsame");      // Magenta
  CaptT->SetTextColor(kolor);
  ycapt += -0.04;
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
  cout<<"================ KKsemMakeHisto  BEGIN ============================"<<endl;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  double CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  CMSene        /= HST_KKMC_NORMA->GetBinContent(511); // farm adjusted
  cout<< "KKsemMakeHisto: CMSene  = "<<CMSene<<endl;
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
  cout<<" MuMu  dsigma/dv, unlimited cos(theta)"<<endl;
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
  cout<<" MuMu  Sigma(vmax) with limited c=cos(theta)"<<endl;
//------------------------------------------------------------------------  
  // ISR*FSR
  KeyDis = 302302;        // ISR*FSR O(alf2)
  sprintf(chak,"XCHI2");  // ISR*FSR Mff
  TH1D *vcum_ISR2_FSR2 =(TH1D*)hstVtemplate->Clone("vcum_ISR2_FSR2");
  LibSem.VVplot(vcum_ISR2_FSR2, KF, chak, KeyDis, KeyFob);
  // with costhe cut
  // does it make sense for ISR*FSR????
  kksem_setcrange_(-22.0/25, 22.0/25);
  TH1D *vcum2_ISR2_FSR2 =(TH1D*)hstVtemplate->Clone("vcum2_ISR2_FSR2");
  LibSem.VVplot(vcum2_ISR2_FSR2, KF, chak, KeyDis, KeyFob);
  //-------------------------------------------------
  cout<<" finally AFB(vmax) for limited c=cos(theta)"<<endl;
  // does it make sense for ISR*FSR????
  kksem_setcrange_(0, 22.0/25); // forward
  TH1D *afb2v_ISR2_FSR2 =(TH1D*)hstVtemplate->Clone("afb2v_ISR2_FSR2");
  LibSem.VVplot(afb2v_ISR2_FSR2, KF, chak, KeyDis, KeyFob);// Forward
  afb2v_ISR2_FSR2->Add(afb2v_ISR2_FSR2, vcum2_ISR2_FSR2, 2.0, -1.0) ; // numerator F-B = 2F-(F+B)
  afb2v_ISR2_FSR2->Divide(vcum2_ISR2_FSR2);                           // finally (F-B)(F+B)
  if(CMSene < 91.0 ) afb2v_ISR2_FSR2->Scale(-1);

//------------------------------------------------------------------------
  cout<<"  MuMu  dsigma/dCosTheta, limited v "<<endl;
//------------------------------------------------------------------------  
  kksem_setcrange_(-1.0,1.0); // back to normal
  // ISR only
  KeyDis = 303;           // ISR O(alf3)
  sprintf(chak,"VRHO2");  // ISR only
  //
  KeyDis = 302302;        // ISR*FSR O(alf2)
  sprintf(chak,"XRHO2");  // ISR only

  KeyFob= -100; // KKsem_BornV, WITH integration, OK
  double vmin=0.0;
  double vmax=0.02;
  TH1D        *cdisKS_ISR2 =(TH1D*)hstCtemplate->Clone("cdisKS_ISR2");
  LibSem.Cplot(cdisKS_ISR2, KF, chak, KeyDis, KeyFob, vmin,vmax);
  //
  KeyFob=    0; // With EW (BornV_Dizet) With integration OK!
  TH1D        *cdisDZ_ISR2 =(TH1D*)hstCtemplate->Clone("cdisDZ_ISR2");
  LibSem.Cplot(cdisDZ_ISR2, KF, chak, KeyDis, KeyFob, vmin,vmax);

  //--------------------------------------------------------------------------
  cout<<"  x-check only, the distribution without Z (for 10GeV only)"<<endl;
  long KeyZet=0;
  kksem_setkeyzet_(KeyZet);
  TH1D *cdis_ISR2_Zoff =(TH1D*)hstCtemplate->Clone("cdis_ISR2_Zoff");
  LibSem.Cplot(cdis_ISR2_Zoff, KF, chak, KeyDis, KeyFob, vmin,vmax);

  // ******************** reprocessing plots from KKsem **********************
  cout<<"  *** reprocessing plots from KKsem ****"<<endl;
  TH1D                  *casyKS_ISR2;
  MakeAFB(cdisKS_ISR2,   casyKS_ISR2);
  casyKS_ISR2->SetName("casyKS_ISR2");
  //if(CMSene < 91.0 ) casyKS_ISR2->Scale(-1);
  //
  TH1D                  *casyDZ_ISR2;
  MakeAFB(cdisDZ_ISR2,   casyDZ_ISR2);
  casyDZ_ISR2->SetName("casyDZ_ISR2");
  //if(CMSene < 91.0 ) casyDZ_ISR2->Scale(-1);
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
  CMSene        /= HST_KKMC_NORMA->GetBinContent(511); // farm adjusted

  //****************************************************************************************
  // Pure MC reprocessing part
  //
  TH2D *sca_vTcPR_Ceex2  = (TH2D*)DiskFileA.Get("sca_vTcPR_Ceex2");
  TH2D *sca_vTcPR_Eex2   = (TH2D*)DiskFileA.Get("sca_vTcPR_Eex2");
  TH2D *sca_vTcPR_Ceex2n = (TH2D*)DiskFileA.Get("sca_vTcPR_Ceex2n");

  // *********  distributions in cos(theta) and limited v *************
  int nbMax=0; // it is cut on vv, 0= no cut
  nbMax=45;    // vMax = 45/50=0.9
  TH1D *Hcth_vTcPR_Ceex2_vmax90  = HstProjC( "Hcth_vTcPR_Ceex2_vmax90",  sca_vTcPR_Ceex2,  nbMax);
  TH1D *Hcas_vTcPR_Ceex2_vmax90  = HstProjCA("Hcas_vTcPR_Ceex2_vmax90",  sca_vTcPR_Ceex2,  nbMax);
  //
  nbMax=15;    // vMax = 15/50=0.3
  TH1D *Hcth_vTcPR_Ceex2_vmax30  = HstProjC( "Hcth_vTcPR_Ceex2_vmax30",  sca_vTcPR_Ceex2,  nbMax);
  TH1D *Hcas_vTcPR_Ceex2_vmax30  = HstProjCA("Hcas_vTcPR_Ceex2_vmax30",  sca_vTcPR_Ceex2,  nbMax);
  //
  nbMax=5;     // vMax = 5/50=0.1
  TH1D *Hcth_vTcPR_Ceex2_vmax10  = HstProjC( "Hcth_vTcPR_Ceex2_vmax10",  sca_vTcPR_Ceex2,  nbMax);
  TH1D *Hcas_vTcPR_Ceex2_vmax10  = HstProjCA("Hcas_vTcPR_Ceex2_vmax10",  sca_vTcPR_Ceex2,  nbMax);
  //
  nbMax=1;     // vMax = 5/50=0.02
  TH1D *Hcth_vTcPR_Ceex2_vmax02  = HstProjC( "Hcth_vTcPR_Ceex2_vmax02",  sca_vTcPR_Ceex2,  nbMax);
  TH1D *Hcas_vTcPR_Ceex2_vmax02  = HstProjCA("Hcas_vTcPR_Ceex2_vmax02",  sca_vTcPR_Ceex2,  nbMax);
  //
  TH1D *Hcth_vTcPR_Ceex2n_vmax02 = HstProjC( "Hcth_vTcPR_Ceex2n_vmax02", sca_vTcPR_Ceex2n, nbMax);
  TH1D *Hcas_vTcPR_Ceex2n_vmax02 = HstProjCA("Hcas_vTcPR_Ceex2n_vmax02", sca_vTcPR_Ceex2n, nbMax);

  //  *********** distrib. of cos(theta) unlimited v
  TH1D                   *Hpro_cosPR_Ceex2;
  ProjY1(sca_vTcPR_Ceex2, Hpro_cosPR_Ceex2);
  Hpro_cosPR_Ceex2->SetName("Hpro_cosPR_Ceex2");

  //  *********** distrib. of v unlimited cos(theta)
  TH1D *Hpro_vT_Ceex2;
  ProjX1(sca_vTcPR_Ceex2, Hpro_vT_Ceex2);
  Hpro_vT_Ceex2->SetName("Hpro_vT_Ceex2");

  //************************** Developement corner *****************************************
  // ******* distrib. of v unlimited c, (failed xcheck)
  //TH1D *hst_ProjV = (TH1D*)sca_vTcPR_Ceex2->ProjectionX("hst_projV",1,50,"e");

  ///****************************************************************************************
  /// Distributions of v=vTrue with limited c=cos(theta)
  //  without cutoff on c=cos(thetaPRD)
  nbMax=0;   // cosThetaMax = 1.0
  TH1D  *HTot_vTcPR_Ceex2  = HstProjV("HTot_vTcPR_Ceex2",sca_vTcPR_Ceex2,nbMax);
  TH1D  *HAfb_vTcPR_Ceex2  = HstProjA("HAfb_vTcPR_Ceex2",sca_vTcPR_Ceex2,nbMax);
  nbMax=22;      // cosThetaMax = 22/25 =0.88
  TH1D  *HTot2_vTcPR_Ceex2  = HstProjV("HTot2_vTcPR_Ceex2",sca_vTcPR_Ceex2,nbMax);
  TH1D  *HAfb2_vTcPR_Ceex2  = HstProjA("HAfb2_vTcPR_Ceex2",sca_vTcPR_Ceex2,nbMax);
  nbMax=0;   // cosThetaMax = 1.0
  TH1D  *HTot_vTcPR_Ceex2n  = HstProjV("HTot_vTcPR_Ceex2n",sca_vTcPR_Ceex2n,nbMax);
  TH1D  *HAfb_vTcPR_Ceex2n  = HstProjA("HAfb_vTcPR_Ceex2n",sca_vTcPR_Ceex2n,nbMax);
  nbMax=22;      // cosThetaMax = 22/25 =0.88
  TH1D  *HTot2_vTcPR_Ceex2n  = HstProjV("HTot2_vTcPR_Ceex2n",sca_vTcPR_Ceex2n,nbMax);
  TH1D  *HAfb2_vTcPR_Ceex2n  = HstProjA("HAfb2_vTcPR_Ceex2n",sca_vTcPR_Ceex2n,nbMax);

  cout<<"================ ReMakeMChisto ENDs  ============================="<<endl;
  cout<<"==================================================================="<<endl;
}//RemakeMChisto






///////////////////////////////////////////////////////////////////////////////////
void FigScatA()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigScat =========================== "<<endl;
  // renormalize histograms in nanobarns
  Double_t CMSene;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  CMSene /= HST_KKMC_NORMA->GetBinContent(511); // farm adjusted
  //
  TH2D *sca_vTcPR_Ceex2 = (TH2D*)DiskFileA.Get("sca_vTcPR_Ceex2");
  TH2D *sca_vTcPR_Eex2  = (TH2D*)DiskFileA.Get("sca_vTcPR_Eex2");
  //
  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cScatA = new TCanvas("cScatA","cScatA", gXcanv,  gYcanv,   1000,  600);
  //                            Name    Title            xoff,yoff, WidPix,HeiPix
  gXcanv += 25, gYcanv += 25;
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
  CMSene /= HST_KKMC_NORMA->GetBinContent(511); // farm adjusted
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
  TCanvas *cFigInfo = new TCanvas("cFigInfo","cFigInfo ", gXcanv, gYcanv,    1000,  1000);
  //                            Name    Title               xoff,yoff, WidPix,HeiPix
  gXcanv += 25, gYcanv += 25;
  //
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
  hst_vTrueCeex2->SetMinimum( 1e-4* hst_vTrueCeex2->GetMaximum() );
  hst_vTrueCeex2->DrawCopy("h");
  //
  hst_vAlepCeex2->SetLineColor(kRed);
  hst_vAlepCeex2->DrawCopy("hsame");
  //
  hst_vXGenCeex2->SetLineColor(4);
  hst_vXGenCeex2->DrawCopy("hsame");
  CaptT->DrawLatex(0.10,0.95,"d#sigma/dv (Ceex2); Black=v_{Bare}, Red=v_{Aleph}, Blue=v_{ISR}");
  //==========plot4==============
  cFigInfo->cd(4);
  //-----------------------------
  hst_Cost1Ceex2->SetStats(0);
  hst_Cost1Ceex2->SetTitle(0);
  hst_Cost1Ceex2->SetMinimum(0);
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
  CMSene /= HST_KKMC_NORMA->GetBinContent(511); // farm adjusted
  //
  TH1D *hst_vTrueCeex2 = (TH1D*)DiskFileA.Get("hst_vTrueCeex2");
  TH1D *hst_vXGenCeex2 = (TH1D*)DiskFileA.Get("hst_vXGenCeex2");
  //
  TH1D *vdis_ISR2      = (TH1D*)DiskFileB.Get("vdis_ISR2");
  TH1D *vdis_ISR2_FSR2 = (TH1D*)DiskFileB.Get("vdis_ISR2_FSR2");
  TH1D *Hpro_vT_Ceex2  = (TH1D*)DiskFileB.Get("Hpro_vT_Ceex2");
  //
  //TH2D *sca_vTcPR_Ceex2 = (TH2D*)DiskFileA.Get("sca_vTcPR_Ceex2");
  //TH2D *sca_vTcPR_Eex2  = (TH2D*)DiskFileA.Get("sca_vTcPR_Eex2");
  //------------------------------------------------------------------------
  //****************************************************************************************
  //****************************************************************************************
  //****************************************************************************************
  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigVtest = new TCanvas("cFigVtest","cFigVtest", gXcanv, gXcanv,    1000, 1000);
  //                            Name    Title               xoff,yoff, WidPix,HeiPix
  gXcanv += 25, gYcanv += 25;
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
  // MC vtrue from scatergram
  Hpro_vT_Ceex2->SetLineColor(kGreen); // green
  Hpro_vT_Ceex2->DrawCopy("same");
  // KKsem ISR+FSR
  vdis_ISR2_FSR2->SetLineColor(kMagenta); // magenta
  vdis_ISR2_FSR2->DrawCopy("hsame");
  // KKsem ISR only
  vdis_ISR2->SetLineColor(kBlue); // blue
  vdis_ISR2->DrawCopy("hsame");
  //
  //HstNew->SetLineColor(7); // cyan
  //HstNew->DrawCopy("same");
  //
  CaptT->DrawLatex(0.02,0.95, "d#sigma/dv(ISR+FSR) KKMC CEEX2: Black, Red; KKsem=Magenta");
  CaptT->DrawLatex(0.02,0.91, "           ISR only, KKsem=Blue");
  //==========plot2==============
  cFigVtest->cd(2);
  hst_vTrueCeex2->Divide(vdis_ISR2_FSR2);
  hst_vTrueCeex2->SetStats(0);
  hst_vTrueCeex2->SetTitle(0);
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
  vdis_ISR2->SetLineColor(kBlue); // blue
  vdis_ISR2->DrawCopy("hsame");
  CaptT->DrawLatex(0.02,0.95,"d#sigma/dv(ISR),  KKMC_CEEX2=Red, Blue=KKsem");
  //==========plot4==============
  cFigVtest->cd(4);
  hst_vXGenCeex2->Divide(vdis_ISR2);
  hst_vXGenCeex2->SetStats(0);
  hst_vXGenCeex2->SetTitle(0);
  hst_vXGenCeex2->SetMinimum(0.85);
  hst_vXGenCeex2->SetMaximum(1.15);
  hst_vXGenCeex2->SetLineColor(kRed);
  hst_vXGenCeex2->DrawCopy("h");  // black
  CaptT->DrawLatex(0.02,0.95,"d#sigma/dv(ISR); KKMC_CEEX2/KKsem");
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
  TH1D *cdisKS_ISR2     = (TH1D*)DiskFileB.Get("cdisKS_ISR2");
  TH1D *Hpro_cosPR_Ceex2= (TH1D*)DiskFileB.Get("Hpro_cosPR_Ceex2");
  //------------------------------------------------------------------------
  //****************************************************************************************
  //************************** Developement corner *****************************************
  //TH1D *hst_ProjV = (TH1D*)sca_vTcPR_Ceex2->ProjectionX("hst_projV",1,50,"e");
  //****************************************************************************************
  //****************************************************************************************
   ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigCtest = new TCanvas("cFigCtest","cFigCtest: cos(thet) dis.", gXcanv, gXcanv,    1000, 1000);
  //                                Name    Title            xoff,yoff, WidPix,HeiPix
  gXcanv += 25, gYcanv += 25;

  cFigCtest->SetFillColor(10);
  ////////////////////////////////////////////////////////////////////////////////
  cFigCtest->Divide( 2,  2);
  //cFigCtest->Divide( 2,  2,     0.0,     0.0,   10);
  //              nx, ny, xmargin, ymargin, color
  //////////////////////////////////////////////
  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.04);
  //==========plot1==============
  cFigCtest->cd(1);
  hst_CosPLCeex2->SetStats(0);
  hst_CosPLCeex2->SetTitle(0);
  hst_CosPLCeex2->SetMinimum(0);
  hst_CosPLCeex2->DrawCopy("h");  // black
  //
  hst_CosPRCeex2->SetLineColor(2); // red
  hst_CosPRCeex2->DrawCopy("hsame");
  //
  hst_Cost1Ceex2->SetLineColor(4); // blue
  hst_Cost1Ceex2->DrawCopy("hsame");
  //
  CaptT->DrawLatex(0.02,0.95,"d#sigma/dc(Ceex2); Black=PL, Red=PRD, Blue=Cth1");
  //==========plot2==============
  cFigCtest->cd(2);
  hst_CosPRCeex2->SetLineColor(kRed); // red
  hst_CosPRCeex2->SetMinimum(0);
  hst_CosPRCeex2->SetTitle(0);
  hst_CosPRCeex2->DrawCopy("h");
  //
  Hpro_cosPR_Ceex2->DrawCopy("same");
  //
  CaptT->DrawLatex(0.02,0.95,"d#sigma/dc(Ceex2); Red PRD");
  //==========plot3==============
  cFigCtest->cd(3);
  hst_CosPRCeex2->SetStats(0);
  hst_CosPRCeex2->SetTitle(0);
  hst_CosPRCeex2->SetLineColor(kRed); // red
  hst_CosPRCeex2->SetMinimum(0);
  hst_CosPRCeex2->DrawCopy("h");
  //
  cdisKS_ISR2->SetLineColor(kBlack);
  cdisKS_ISR2->DrawCopy("hsame");  // black
  //
  CaptT->DrawLatex(0.02,0.95,"d#sigma/dc (Ceex2); Red=PRD, Black=KKsem ");
  //==========plot4==============
  cFigCtest->cd(4);
  hst_CosPREex2->SetMinimum(0);
  hst_CosPREex2->SetTitle(0);
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
  Double_t CMSene;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  CMSene /= HST_KKMC_NORMA->GetBinContent(511); // farm adjusted
  char TextEne[100]; sprintf(TextEne,"#sqrt{s} =%4.2fGeV", CMSene);
  //
  // KKsem
  TH1D *vcum_ISR2_FSR2   = (TH1D*)DiskFileB.Get("vcum_ISR2_FSR2");
  TH1D *vcum2_ISR2_FSR2  = (TH1D*)DiskFileB.Get("vcum2_ISR2_FSR2"); // with ctheta cut

  /// Distributions of v=vTrue with limited cos(thetaPRD)
  // without cutoff on c=cos(thetaPRD)
  TH1D *HTot_vTcPR_Ceex2  = (TH1D*)DiskFileB.Get("HTot_vTcPR_Ceex2");
  TH1D *HAfb_vTcPR_Ceex2  = (TH1D*)DiskFileB.Get("HAfb_vTcPR_Ceex2");
  // cosThetaMax = 22/25 =0.88, IFI on
  TH1D *HTot2_vTcPR_Ceex2 = (TH1D*)DiskFileB.Get("HTot2_vTcPR_Ceex2");
  TH1D *HAfb2_vTcPR_Ceex2 = (TH1D*)DiskFileB.Get("HAfb2_vTcPR_Ceex2");
  // cosThetaMax = 22/25 =0.88, IFI off
  TH1D *HTot2_vTcPR_Ceex2n= (TH1D*)DiskFileB.Get("HTot2_vTcPR_Ceex2n");
  TH1D *HAfb2_vTcPR_Ceex2n= (TH1D*)DiskFileB.Get("HAfb2_vTcPR_Ceex2n");
  //
  TH1D *afb2v_ISR2_FSR2   = (TH1D*)DiskFileB.Get("afb2v_ISR2_FSR2");
  //
  TH1D *HTot_vTcPR_Ceex2n = (TH1D*)DiskFileB.Get("HTot_vTcPR_Ceex2n");
  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigVprod = new TCanvas("cFigVprod","FigVprod", gXcanv, gYcanv,    1000, 1000);
  //                                   Name    Title               xoff,yoff, WidPix,HeiPix
  gXcanv += 25, gYcanv += 25;
  cFigVprod->SetFillColor(10);
  ////////////////////////////////////////////////////////////////////////////////
  cFigVprod->Divide( 2,  2);
  //cFigVprod->Divide( 2,  2,     0.0,     0.0,   10);
  //                  nx, ny, xmargin, ymargin, color
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
  HTot_vTcPR_Ceex2->SetLineColor(kBlue); // blue
  HTot_vTcPR_Ceex2->DrawCopy("h");
  HTot2_vTcPR_Ceex2->SetLineColor(kRed); // red
  HTot2_vTcPR_Ceex2->DrawCopy("hsame");
  // KKsem
  vcum_ISR2_FSR2->DrawCopy("hsame");
  vcum2_ISR2_FSR2->DrawCopy("hsame");
  //
  CaptT->DrawLatex(0.06,0.95,"(a) KKMC: #sigma^{IFIon}_{ISR+FSR}(v_{max}), Blue: |cos(#theta|<1, Red: |cos(#theta|<0.88");
  CaptT->DrawLatex(0.12,0.85,"    KKsem (IFI off): Black ");
  CaptTb->Draw();
  //==========plot2==============
  cFigVprod->cd(2);
  HTot_vTcPR_Ceex2->Divide(vcum_ISR2_FSR2);
  HTot_vTcPR_Ceex2->SetMinimum(0.98);
  HTot_vTcPR_Ceex2->SetMaximum(1.02);
  HTot_vTcPR_Ceex2->DrawCopy("h");
  //
  HTot2_vTcPR_Ceex2->Divide(vcum2_ISR2_FSR2);
  HTot2_vTcPR_Ceex2->DrawCopy("hsame");
  //
  HTot_vTcPR_Ceex2n->SetLineColor(kMagenta);
  HTot_vTcPR_Ceex2n->Divide(vcum_ISR2_FSR2);
  HTot_vTcPR_Ceex2n->DrawCopy("hsame");

  CaptT->DrawLatex(0.12,0.95,"(b) Ceex2/KKsem, IFIon, Blue=|cos(#theta|<1, Red=|cos(#theta|<0.88 ");
  CaptT->DrawLatex(0.16,0.85,"    Magenta: Ceex2/KKsem for IFIoff, |cos(#theta|<1 ");
  CaptT->DrawLatex(0.60,0.75,TextEne);
  CaptTb->Draw();
  //==========plot3==============
  cFigVprod->cd(3);
  HAfb_vTcPR_Ceex2->SetStats(0);
  HAfb_vTcPR_Ceex2->SetTitle(0);
  HAfb_vTcPR_Ceex2->SetLineColor(kBlue); // blue
  HAfb2_vTcPR_Ceex2->SetTitle(0);
  HAfb2_vTcPR_Ceex2->SetStats(0);
  HAfb2_vTcPR_Ceex2->SetLineColor(kRed); // red
  //HAfb2_vTcPR_Ceex2->SetMinimum(0.12);
  //HAfb2_vTcPR_Ceex2->SetMaximum(0.32);
  HAfb2_vTcPR_Ceex2->DrawCopy("h");
  //
  HAfb_vTcPR_Ceex2->DrawCopy("hsame");
  CaptT->DrawLatex(0.06,0.95,"(c) A^{KKMC}_{FB}, Blue=|cos(#theta|<1, Red=|cos(#theta|<0.88");
  CaptTb->Draw();
  //==========plot4==============
  cFigVprod->cd(4);

  HAfb2_vTcPR_Ceex2->SetLineColor(kRed); // red
  HAfb2_vTcPR_Ceex2->SetStats(0);
  HAfb2_vTcPR_Ceex2->DrawCopy("h");
  //
  HAfb2_vTcPR_Ceex2n->SetLineColor(kGreen); // green
  HAfb2_vTcPR_Ceex2n->SetLineWidth(2);
  HAfb2_vTcPR_Ceex2n->DrawCopy("hsame");
  //
  afb2v_ISR2_FSR2->SetLineColor(kBlack);
  afb2v_ISR2_FSR2->DrawCopy("hsame");
  //
  CaptT->DrawLatex(0.12,0.95,"(d) A^{KKMC}_{FB}(|cos(#theta|<0.88), Red/Green=IFIon/off");
  CaptT->DrawLatex(0.17,0.85,"  A^{KKsem}_{FB}(|cos(#theta|<0.88), Black=IFIoff");
  //----------------------------
  cFigVprod->cd();
  //================================================
}//FigVprod


///////////////////////////////////////////////////////////////////////////////////
void FigCprod()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigCprod =========================== "<<endl;
  //
  Double_t CMSene;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  CMSene /= HST_KKMC_NORMA->GetBinContent(511); // farm adjusted
  char TextEne[100]; sprintf(TextEne,"#sqrt{s} =%4.2fGeV", CMSene);
  //
  // from KKsem
  TH1D *cdisKS_ISR2      = (TH1D*)DiskFileB.Get("cdisKS_ISR2");
  TH1D *cdisDZ_ISR2      = (TH1D*)DiskFileB.Get("cdisDZ_ISR2");
  TH1D *casyKS_ISR2      = (TH1D*)DiskFileB.Get("casyKS_ISR2");
  TH1D *casyDZ_ISR2      = (TH1D*)DiskFileB.Get("casyDZ_ISR2");
  TH1D *cdis_ISR2_Zoff   = (TH1D*)DiskFileB.Get("cdis_ISR2_Zoff"); // no Z, pure QED
  // from KKMC
  TH1D *Hcth_vTcPR_Ceex2_vmax90 = (TH1D*)DiskFileB.Get("Hcth_vTcPR_Ceex2_vmax90");
  TH1D *Hcas_vTcPR_Ceex2_vmax90 = (TH1D*)DiskFileB.Get("Hcas_vTcPR_Ceex2_vmax90");
  TH1D *Hcth_vTcPR_Ceex2_vmax30 = (TH1D*)DiskFileB.Get("Hcth_vTcPR_Ceex2_vmax30");
  TH1D *Hcth_vTcPR_Ceex2_vmax10 = (TH1D*)DiskFileB.Get("Hcth_vTcPR_Ceex2_vmax10");
  TH1D *Hcth_vTcPR_Ceex2_vmax02 = (TH1D*)DiskFileB.Get("Hcth_vTcPR_Ceex2_vmax02");
  TH1D *Hcth_vTcPR_Ceex2n_vmax02= (TH1D*)DiskFileB.Get("Hcth_vTcPR_Ceex2n_vmax02");
  //
  TH1D *Hcas_vTcPR_Ceex2_vmax02 = (TH1D*)DiskFileB.Get("Hcas_vTcPR_Ceex2_vmax02");
  TH1D *Hcas_vTcPR_Ceex2n_vmax02=(TH1D*)DiskFileB.Get("Hcas_vTcPR_Ceex2n_vmax02");
  TH1D *Hcas_vTcPR_Ceex2_vmax10 = (TH1D*)DiskFileB.Get("Hcas_vTcPR_Ceex2_vmax10");
  TH1D *Hcas_vTcPR_Ceex2_vmax30 = (TH1D*)DiskFileB.Get("Hcas_vTcPR_Ceex2_vmax30");

  //****************************************************************************************
  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigCprod = new TCanvas("cFigCprod","cFigCprod", gXcanv, gYcanv,    1000, 1000);
  //                            Name    Title               xoff,yoff, WidPix,HeiPix
  gXcanv += 25, gYcanv += 25;
  cFigCprod->SetFillColor(10);
  ////////////////////////////////////////////////////////////////////////////////
  //cFigCprod->Divide( 2,  1);  // lower raw temporarily off
  cFigCprod->Divide( 2,  2);
  //cFigCprod->Divide( 2,  2,     0.0,     0.0,   10);
  //                nx, ny, xmargin, ymargin, color
  //////////////////////////////////////////////
  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.04);
  //==========plot1==============
  cFigCprod->cd(1);
  TH1D *Hst, *Hst2;
  Hst = cdisDZ_ISR2; //ISR*FSR EW on (Dizet)
  Hst->SetStats(0);
  Hst->SetTitle(0);
  Hst->GetYaxis()->CenterTitle();
  Hst->GetYaxis()->SetTitleSize(0.05);
  Hst->GetYaxis()->SetTitle("d#sigma/dcos(#theta)");
  Hst->GetXaxis()->SetTitleSize(0.04);
  Hst->GetXaxis()->SetTitle("cos(#theta)");
  Hst->SetMinimum(0);
  Hst->SetLineColor(kGreen);
  Hst->SetLineWidth(2);
  Hst->DrawCopy("h");
  // EW off
  Hst2 = cdisKS_ISR2;  // EW switched off
  Hst2->SetLineColor(kMagenta);
  //Hst2->DrawCopy("hsame");
  // CEEX2 ISR*FSR IFIon
  Hcth_vTcPR_Ceex2_vmax02->SetLineColor(kBlue);
  Hcth_vTcPR_Ceex2_vmax02->DrawCopy("hsame");
  // CEEX2  ISR*FSR IFIoff
  Hcth_vTcPR_Ceex2n_vmax02->SetLineColor(kBlack);
  Hcth_vTcPR_Ceex2n_vmax02->DrawCopy("hsame");
  //
  CaptT->DrawLatex(0.12,0.96," (a) CEEX2: Blue=IFIon, Black=IFIoff, v_{Bare}<0.02, ISR*FSR");
  CaptT->DrawLatex(0.12,0.91,"     KKsem: Green=IFIoff, ISR*FSR");
  //CaptT->DrawLatex(0.12,0.85,"     KKsem: Magenta=IFIoff&ELWoff");
  //==========plot2==============
  cFigCprod->cd(2);
  TH1D *hZero = (TH1D*)cdisKS_ISR2->Clone("hZero"); // zero line
  hZero->Reset();
  Hst = casyDZ_ISR2;  // ELW on
  Hst->SetStats(0);
  Hst->SetTitle(0);
  Hst->GetYaxis()->CenterTitle();
  Hst->GetYaxis()->SetTitleSize(0.05);
  Hst->GetYaxis()->SetTitle("AFB(#theta)");
  Hst->GetXaxis()->SetTitleSize(0.04);
  Hst->GetXaxis()->SetTitle("cos(#theta_{max})");
  //
  //gPad->DrawFrame(0.7, 0.0,   1.0, 0.8);
  gPad->DrawFrame(15./25., 0.0,   1.0, 0.6, " (b); cos(#theta_{max}) ; A_{FB}");
  Hst->SetLineColor(kGreen);
  Hst->SetLineWidth(2);
  Hst->DrawCopy("hsame");
  //
  Hst2=casyKS_ISR2; // ELW off
  Hst2->SetLineColor(kMagenta);
  //Hst2->DrawCopy("hsame");
  //
  Hcas_vTcPR_Ceex2_vmax02->SetLineColor(kBlue);
  Hcas_vTcPR_Ceex2_vmax02->DrawCopy("hsame");
  //
  Hcas_vTcPR_Ceex2n_vmax02->SetLineColor(kBlack); // no IFI
  Hcas_vTcPR_Ceex2n_vmax02->DrawCopy("hsame");   // no IFI
  //
  //hZero->DrawCopy("hsame");
  //
  CaptT->DrawLatex(0.11,0.86," CEEX2: v_{Bare}<0.02, Blue=IFIon, Black=IFIoff, ISR*FSR");
  CaptT->DrawLatex(0.11,0.81," KKsem: Green=IFIoff,  ISR*FSR");
  //CaptT->DrawLatex(0.15,0.77," KKsem: Magenta=IFIoff&ELWoff ");
  CaptT->DrawLatex(0.60,0.75,TextEne);
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
  gPad->DrawFrame(0.0, 0.0,   1.0, 0.6, " (b); cos(#theta_{max}) ; A_{FB}");

  Hst=Hcas_vTcPR_Ceex2_vmax02;
  Hst->SetStats(0);
  Hst->SetTitle(0);
  Hst->GetYaxis()->CenterTitle();
  Hst->GetYaxis()->SetTitleSize(0.05);
  Hst->GetYaxis()->SetTitle("AFB(#theta)");
  Hst->GetXaxis()->SetTitleSize(0.04);
  Hst->GetXaxis()->SetTitle("cos(#theta)");
  Hst->SetLineColor(kBlue);
  Hst->DrawCopy("hsame");
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
  CaptT->DrawLatex(0.12,0.86,"CEEX2: Black/Red/Green/Blue =  v< 0.9, 0.3, 0.10, 0.02");
  CaptT->DrawLatex(0.12,0.81,"CEEX2: IFI off, Cyan v< 0.02");

  //----------------------------
  cFigCprod->cd();
  //================================================
}//FigCprod


///////////////////////////////////////////////////////////////////////////////////
void FigCosThe()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigCosThe =========================== "<<endl;

  TH1D *Hcth_vTcPR_Ceex2_vmax02 = (TH1D*)DiskFileB.Get("Hcth_vTcPR_Ceex2_vmax02");
  TH1D *Hcth_vTcPR_Ceex2n_vmax02= (TH1D*)DiskFileB.Get("Hcth_vTcPR_Ceex2n_vmax02");

  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.04);
 ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigCosThe = new TCanvas("cFigCosThe","cFigCosThe", gXcanv, gYcanv,    600, 600);
  //                                      Name    Title        xoff,yoff, WidPix,HeiPix
  gXcanv += 25, gYcanv += 25;
  cFigCosThe->SetFillColor(10);

  cFigCosThe->cd();

  TH1D *Hst=Hcth_vTcPR_Ceex2_vmax02;
  Hst->SetStats(0);
  Hst->SetTitle(0);
  Hst->GetXaxis()->CenterTitle();
  Hst->GetYaxis()->SetTitleSize(0.04);
  Hst->GetYaxis()->SetTitle("d#sigma/dcos(#theta) [nb]");
  Hst->GetXaxis()->SetTitleSize(0.04);
  Hst->GetXaxis()->SetTitle("cos(#theta)");

  Hst->SetLineColor(kBlue);
  Hst->SetMinimum(0.0);
  Hst->DrawCopy("h");

  double ycapt = 0.80;
  CaptT->DrawLatex(0.40, ycapt,gTextEne);
  PlotSame2(Hcth_vTcPR_Ceex2_vmax02,  ycapt, kBlue,   +0.80, "(a)", "KKMC, IFI on ");
  PlotSame2(Hcth_vTcPR_Ceex2n_vmax02, ycapt, kBlack,  -0.80, "(b)", "KKMC, IFI off ");


  //================================================
}//FigCosThe

///////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  //++++++++++++++++++++++++++++++++++++++++
  TApplication theApp("theApp", &argc, argv);
  //++++++++++++++++++++++++++++++++++++++++
  LibSem.Initialize(DiskFileA);
  /////////////////////////////////////////////////////////
  // Reading directly KKMC input (farming)
  int Nodes, Nodes2;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  Nodes    = HST_KKMC_NORMA->GetBinContent(511);       // No of farm nodes (trick)
  gCMSene  = HST_KKMC_NORMA->GetBinContent(1)/Nodes;   // CMSene=xpar(1), farn adjusted
  gNevTot  = HST_KKMC_NORMA->GetEntries();             // MC statistics from KKMC
  sprintf(gTextEne,"#sqrt{s} =%4.2fGeV", gCMSene);
  sprintf(gTextNev,"KKMC:%10.2e events", gNevTot);
  //
  /*
  TH1D *HST_FOAM_NORMA3 = (TH1D*)DiskFileF.Get("HST_FOAM_NORMA3");
  Nodes2   =  HST_FOAM_NORMA3->GetBinContent(511);    // No of farm nodes (trick)
  double  CMSeneF  = HST_FOAM_NORMA3->GetBinContent(1)/Nodes2; // CMSene=xpar(1)
  if( fabs(gCMSene/CMSeneF-1) >1e-4 ){
	  cout<<" +++++ Wrong input files !!!! KKMC "<< gCMSene <<"GeV and  FOAM "<< CMSeneF<<"GeV"<<endl;
	  exit(19);
  }
  gNevTot2  = HST_FOAM_NORMA3->GetEntries();       // MC statistics from KKMC
  sprintf(gTextNev2,"FOAM:%10.2e events", gNevTot2);
  */
  //
  cout<< "CMSene[GeV] = "<< gCMSene<< endl;
  cout<< "KKMC: No. of farm nodes="<< Nodes  << "  Tot no. of events = "<<gNevTot<< endl;
  //cout<< "FOAM: No. of farm nodes="<< Nodes2 << "  Tot no. of events = "<<gNevTot2<<endl;
//////////////////////////////////////////////////////////////////////////

  HistNormalize();     // Renormalization of MC histograms
  ReMakeMChisto();     // reprocessing MC histos
  //
  KKsemMakeHisto();    // prepare histos from KKsem
  //========== PLOTTING ==========
  //
  //FigScatA();
  //FigInfo();

  FigVtest();  // introduct. tests/calibrations
  FigCtest();  // introduct. tests/calibrations

  FigVprod();
  FigCprod();
  //
  FigCosThe();
  //++++++++++++++++++++++++++++++++++++++++
  DiskFileA.ls();
  DiskFileB.ls();
  DiskFileB.Write();
  DiskFileB.Close();
  //++++++++++++++++++++++++++++++++++++++++
  theApp.Run();
  //++++++++++++++++++++++++++++++++++++++++
}
