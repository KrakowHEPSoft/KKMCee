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
// KKMC current
TFile DiskFileA("../workKKMC/histo.root");
// Jan. 2018
//TFile DiskFileA("../workKKMC/histo.root_88GeV_11G"); // Jan. 2018
//TFile DiskFileA("../workKKMC/histo.root_10GeV_10G"); // Jan. 2018
//
//TFile DiskFileA("../workKKMC/histo.root_95GeV_26G");
//TFile DiskFileA("../workKKMC/histo.root_10GeV_5.8G"); //

//////  *** KKFOAM, bacis xcheck including IFI
//TFile DiskFileF("../workFOAM/histo.root"); // current
// Dec 2017 run
//TFile DiskFileF("../workFOAM/histo.root_88GeV_22G");
TFile DiskFileF("../workFOAM/histo.root_95GeV_23G");
//TFile DiskFileF("../workFOAM/histo.root_95GeV_23G");
//TFile DiskFileF("../workFOAM/histo.root_10GeV_18G");

// Sept. 2017 runs
//TFile DiskFileF("../workFOAM/histo.root_95GeV_57G");  // last
//TFile DiskFileF("../workFOAM/histo.root_88GeV_15G");
//TFile DiskFileF("../workFOAM/histo.root_10GeV_25G");

//////  *** KKFOAM1  Soft limit study
//TFile DiskFileF2("../workFOAM1/histo.root"); // current
// Dec 2017 run
//TFile DiskFileF2("../workFOAM1/histo.root_88GeV_13G"); //
TFile DiskFileF2("../workFOAM1/histo.root_95GeV_14G"); //
//TFile DiskFileF2("../workFOAM1/histo.root_10GeV_14G"); //

//////////////////OBSOLETE Oldies ///////////////////
//TFile DiskFileA("../workAFB/rmain.root");
// archive obsolete!
//TFile DiskFileA("../workAFB/rmain.root_95GeV_100M");
//TFile DiskFileA("../workAFB/rmain.root_88GeV_100M"); // archive
//TFile DiskFileA("../workAFB/rmain.root_10GeV_30M");

TFile DiskFileB("RhoSemi.root","RECREATE","histograms");
//=============================================================================

KKplot LibSem("KKplot");

///////////////////////////////////////////////////////////////////////////////////
//              GLOBAL stuff
///////////////////////////////////////////////////////////////////////////////////
double gCMSene, gNevTot, gNevTot2; // from KKMC and KKfoam MC runs (histograms)
char   gTextEne[100], gTextNev[100], gTextNev2[100];
int    kGold=kOrange-3, kBrune=46, kPine=kGreen+3;
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
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sct_vTcPL_Ceex2") );
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sct_vTcPL_Ceex2n") );

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
void ReMakeKKMC(){
// Some KKMC histos are preprocessed
//------------------------------------------------------------------------
  cout<<"==================================================================="<<endl;
  cout<<"================ ReMakeMChisto  BEGIN  ============================"<<endl;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  double CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  CMSene        /= HST_KKMC_NORMA->GetBinContent(511); // farm adjusted

  //****************************************************************************************
  // Pure MC reprocessing part
  // wide range scaterplots, v<1, 50x50
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
  //=======================================================
  // narower range scaterplots, v<0.2, 100x100
  TH2D *sct_vTcPL_Ceex2  = (TH2D*)DiskFileA.Get("sct_vTcPL_Ceex2");
  TH2D *sct_vTcPL_Ceex2n = (TH2D*)DiskFileA.Get("sct_vTcPL_Ceex2n");
  nbMax=10;     // vMax = 10/100*0.2=0.02
  TH1D *Hcth_vTcPL_Ceex2_vmax02  = HstProjC( "Hcth_vTcPL_Ceex2_vmax02",  sct_vTcPL_Ceex2,  nbMax);
  TH1D *Hcas_vTcPL_Ceex2_vmax02  = HstProjCA("Hcas_vTcPL_Ceex2_vmax02",  sct_vTcPL_Ceex2,  nbMax);
  //
  TH1D *Hcth_vTcPL_Ceex2n_vmax02 = HstProjC( "Hcth_vTcPL_Ceex2n_vmax02", sct_vTcPL_Ceex2n, nbMax);
  TH1D *Hcas_vTcPL_Ceex2n_vmax02 = HstProjCA("Hcas_vTcPL_Ceex2n_vmax02", sct_vTcPL_Ceex2n, nbMax);
  //
  nbMax=1;     // vMax = 1/100*0.2=0.002
  TH1D *Hcth_vTcPL_Ceex2_vmax002  = HstProjC( "Hcth_vTcPL_Ceex2_vmax002",  sct_vTcPL_Ceex2,  nbMax);
  TH1D *Hcth_vTcPL_Ceex2n_vmax002 = HstProjC( "Hcth_vTcPL_Ceex2n_vmax002", sct_vTcPL_Ceex2n, nbMax);

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
//void HisReMakeFoam35(TFile *DiskFileF, int NbMax, int NbMax2){
void HisReMakeFoam35(){

//
  cout<<"==================================================================="<<endl;
  cout<<"================ ReMakeFoam35  BEGIN   ============================"<<endl;
//////////////////////////////////////////////////////////////////
  cout<<"  Renormalizing  and reprocessing histograms from FOAM"<<endl;

  TH1D *HST_FOAM_NORMA3 = (TH1D*)DiskFileF.Get("HST_FOAM_NORMA3");
  TH1D *HST_FOAM_NORMA5 = (TH1D*)DiskFileF.Get("HST_FOAM_NORMA5");

  TH2D *SCT_xc_Ceex2  = (TH2D*)DiskFileF.Get("SCT_xc_Ceex2");   // FOAM5 small range x<0.20
  TH2D *SCT_xc_Ceex2n = (TH2D*)DiskFileF.Get("SCT_xc_Ceex2n");  // FOAM3 small range x<0.20

  HisNorm2(HST_FOAM_NORMA5, SCT_xc_Ceex2 );   // normalizing FOAM5
  HisNorm2(HST_FOAM_NORMA3, SCT_xc_Ceex2n );  // normalizing FOAM3


  int nbMax=10;     // vMax = 10/100*0.2=0.02
  TH1D *Hcth_foam_Ceex2_vmax02  = HstProjC( "Hcth_foam_Ceex2_vmax02",  SCT_xc_Ceex2,  nbMax);
  TH1D *Hcth_foam_Ceex2n_vmax02 = HstProjC( "Hcth_foam_Ceex2n_vmax02", SCT_xc_Ceex2n, nbMax);
  nbMax=1;     // vMax = 1/100*0.2=0.002
  TH1D *Hcth_foam_Ceex2_vmax002  = HstProjC( "Hcth_foam_Ceex2_vmax002",  SCT_xc_Ceex2,  nbMax);
  TH1D *Hcth_foam_Ceex2n_vmax002 = HstProjC( "Hcth_foam_Ceex2n_vmax002", SCT_xc_Ceex2n, nbMax);

  HisNorm1(HST_FOAM_NORMA3, (TH1D*)DiskFileF.Get("HST_cc_EEX2_vmax02") );
  HisNorm1(HST_FOAM_NORMA3, (TH1D*)DiskFileF.Get("HST_cc_EEX2_vmax002") );
  HisNorm1(HST_FOAM_NORMA3, (TH1D*)DiskFileF.Get("HST_cc_EEX2_vmax0002") );


  //////////////////////////////////////////////////////////////////
  // testing soft limit using Foam
  TH1D *HST_FOAM_NORMA2   = (TH1D*)DiskFileF2.Get("HST_FOAM_NORMA2");

  HisNorm1(HST_FOAM_NORMA2, (TH1D*)DiskFileF2.Get("HST_cs_EEX2_vmax02") );
  HisNorm1(HST_FOAM_NORMA2, (TH1D*)DiskFileF2.Get("HST_cs_EEX2_vmax002") );
  HisNorm1(HST_FOAM_NORMA2, (TH1D*)DiskFileF2.Get("HST_cs_EEX2_vmax0002") );

  HisNorm1(HST_FOAM_NORMA2, (TH1D*)DiskFileF2.Get("HST_cs_ceex2n_vmax02") );
  HisNorm1(HST_FOAM_NORMA2, (TH1D*)DiskFileF2.Get("HST_cs_ceex2n_vmax002") );
  HisNorm1(HST_FOAM_NORMA2, (TH1D*)DiskFileF2.Get("HST_cs_ceex2n_vmax0002") );


  HisNorm1(HST_FOAM_NORMA2, (TH1D*)DiskFileF2.Get("HST_cs_ceex2_vmax02") );
  HisNorm1(HST_FOAM_NORMA2, (TH1D*)DiskFileF2.Get("HST_cs_ceex2_vmax002") );
  HisNorm1(HST_FOAM_NORMA2, (TH1D*)DiskFileF2.Get("HST_cs_ceex2_vmax0002") );

  cout<<"================ ReMakeFoam35 ENDs  ==============================="<<endl;
  cout<<"==================================================================="<<endl;
}//ReMakeFoam35



///////////////////////////////////////////////////////////////////////////////////
void FigScatV()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigScat =========================== "<<endl;
  //
  TH2D *sca_vTvA_Eex2  = (TH2D*)DiskFileA.Get("sca_vTvA_Eex2");
  TH2D *sca_vKvA_Eex2  = (TH2D*)DiskFileA.Get("sca_vKvA_Eex2");
  //
  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.04);
  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cScatV = new TCanvas("cScatV","cScatV", gXcanv,  gYcanv,   1600,  800);
  //                            Name    Title            xoff,yoff, WidPix,HeiPix
  gXcanv += 25, gYcanv += 25;
  cScatV->SetFillColor(10);
  cScatV->Divide( 2,  0);
  //cScatV->Divide( 2,  0,     0.0,     0.0,   10);
  //              nx, ny, xmargin, ymargin, color
  //////////////////////////////////////////////
  TString OptSurf;
  //OptSurf="      "; // 2D scatergram, points
  //OptSurf="col   "; // 2D histogram, color
  OptSurf="colz  "; // 2D kolorowe paski, ze skala
  //OptSurf="surf1 "; // 3D surface color
  //OptSurf="lego2 "; // 3D histogram color
  //OptSurf="surf3 "; // 3D histogram, z plotem "na dachu"
  //OptSurf="surf2z"; // 3D kolorowe paski, ze skala
  //OptSurf="surf2 "; // 3D kolorowe paski bez skali
  //OptSurf="surf4 "; // 3D gladka powierchnia
  //-------------------------------------
  cScatV->cd(1);
  gPad->SetLogz(); // !!!!!!
  gPad->SetTheta(25);
  gPad->SetPhi( -38);
  //
  sca_vTvA_Eex2->GetYaxis()->CenterTitle();
  sca_vTvA_Eex2->GetYaxis()->SetTitleOffset(1.4);
  sca_vTvA_Eex2->GetYaxis()->SetTitleSize(0.035);
  sca_vTvA_Eex2->GetYaxis()->SetNdivisions(5);
  sca_vTvA_Eex2->GetYaxis()->SetTitle("vA=1-z_{ALEPH}");
  sca_vTvA_Eex2->GetXaxis()->CenterTitle();
  sca_vTvA_Eex2->GetXaxis()->SetNdivisions(5);
  sca_vTvA_Eex2->GetXaxis()->SetTitle("vT=1-M^{2}/s");
  //
  sca_vTvA_Eex2->Draw(OptSurf);
  //
  CaptT->DrawLatex(0.20,0.75, "KKMC ISR+FSR");
//  CaptT->DrawLatex(0.20,0.75, "KKMC ISR, FSR off");
  //-------------------------------------
  cScatV->cd(2);
  gPad->SetLogz(); // !!!!!!

  gPad->SetTheta(25);
  gPad->SetPhi( -38);
  //
  //
  sca_vKvA_Eex2->GetYaxis()->CenterTitle();
  sca_vKvA_Eex2->GetYaxis()->SetTitleOffset(1.4);
  sca_vKvA_Eex2->GetYaxis()->SetTitleSize(0.035);
  sca_vKvA_Eex2->GetYaxis()->SetNdivisions(5);
  sca_vKvA_Eex2->GetYaxis()->SetTitle("vAleph=1-z_{AEPH}");
  sca_vKvA_Eex2->GetXaxis()->CenterTitle();
  sca_vKvA_Eex2->GetXaxis()->SetNdivisions(5);
  sca_vKvA_Eex2->GetXaxis()->SetTitle("vK=1-M^{2}_{ISR}/s");
  //
  sca_vKvA_Eex2->Draw(OptSurf);
  //
  CaptT->DrawLatex(0.20,0.75, "KKMC ISR+FSR");
//  CaptT->DrawLatex(0.20,0.75, "KKMC ISR, FSR off");

  cScatV->cd();

  cScatV->SaveAs("cScatV.pdf");

}
// FigScatV()



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
  Hpro_vT_Ceex2->SetLineColor(kPine); // green
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
  HAfb2_vTcPR_Ceex2n->SetLineColor(kPine); // green
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
  Hst->SetLineColor(kPine);
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
  Hst->SetLineColor(kPine);
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
  Hcth_vTcPR_Ceex2_vmax10->SetLineColor(kPine);
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
  Hcas_vTcPR_Ceex2_vmax10->SetLineColor(kPine);
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
// Older version with PRD angle
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

  TH1D *Hst=Hcth_vTcPR_Ceex2n_vmax02;

  Hst->SetStats(0);
  Hst->SetTitle(0);
  Hst->GetXaxis()->CenterTitle();
  //Hst->GetYaxis()->SetTitleSize(0.04);
  //Hst->GetYaxis()->SetTitle("d#sigma/dcos(#theta) [nb]");
  Hst->GetXaxis()->SetTitleSize(0.04);
  Hst->GetXaxis()->SetTitle("cos(#theta)");

  Hst->SetLineColor(kBlue);
  Hst->SetMinimum(0.0);
  Hst->DrawCopy("h");

  double ycapt = 0.80;
  CaptT->DrawLatex(0.40, ycapt,gTextEne);
  PlotSame2(Hcth_vTcPR_Ceex2_vmax02,  ycapt, kBlue,   +0.80, "(a)", "KKMC, IFI on ");
  PlotSame2(Hcth_vTcPR_Ceex2n_vmax02, ycapt, kBlack,  -0.80, "(b)", "KKMC, IFI off ");

  CaptT->DrawLatex(0.00,0.96,"d#sigma/d(cos #theta)");

  //================================================
}//FigCosThe


///////////////////////////////////////////////////////////////////////////////////
void FigCosThe2()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigCosThe2 =========================== "<<endl;
 //
  TH1D *Hcth_vTcPL_Ceex2_vmax002 = (TH1D*)DiskFileB.Get("Hcth_vTcPL_Ceex2_vmax02");
  TH1D *Hcth_vTcPL_Ceex2n_vmax002= (TH1D*)DiskFileB.Get("Hcth_vTcPL_Ceex2n_vmax02");
  //
  TH1D *Hcth_foam_Ceex2_vmax002  = (TH1D*)DiskFileB.Get("Hcth_foam_Ceex2_vmax02");
  TH1D *Hcth_foam_Ceex2n_vmax002 = (TH1D*)DiskFileB.Get("Hcth_foam_Ceex2n_vmax02");

  TH1D *HST_cs_ceex2_vmax002     = (TH1D*)DiskFileF2.Get("HST_cs_ceex2_vmax02");

  //////
  TH1D *hOne = (TH1D*)Hcth_foam_Ceex2_vmax002->Clone("hOne");  // zero line
  for(int i=1; i <= hOne->GetNbinsX() ; i++) { hOne->SetBinContent(i, 1); hOne->SetBinError(i, 0);}

////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.04);
 ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigCosThe2 = new TCanvas("cFigCosThe2","cFigCosThe2", gXcanv, gYcanv,    1200, 600);
  //                                      Name    Title        xoff,yoff, WidPix,HeiPix
  gXcanv += 25, gYcanv += 25;
  cFigCosThe2->SetFillColor(10);
  ////////////////////////////////////////////////////////////////////////////////
  cFigCosThe2->Divide( 2,  0);
  //====================plot1========================
  cFigCosThe2->cd(1);
//  TH1D *Hst=Hcth_vTcPR_Ceex2_vmax002;
  TH1D *Hst=Hcth_vTcPL_Ceex2n_vmax002;
  Hst->SetStats(0);
  Hst->SetTitle(0);
  Hst->GetXaxis()->CenterTitle();
  //Hst->GetYaxis()->SetTitleSize(0.04);
  //Hst->GetYaxis()->SetTitle("d#sigma/dcos(#theta) [nb]");
  Hst->GetXaxis()->SetTitleSize(0.04);
  Hst->GetXaxis()->SetTitle("cos #theta");

  Hst->SetLineColor(kBlue);
  Hst->SetMinimum(0.0);
  Hst->DrawCopy("h");

  CaptT->DrawLatex(0.00,0.96,"d#sigma/d(cos #theta)");
  double ycapt = 0.80;
  CaptT->DrawLatex(0.40, ycapt,gTextEne);
  PlotSame2(Hcth_vTcPL_Ceex2_vmax002,  ycapt, kBlue,   +0.80, "(a)", "KKMC,   IFI on ");
  PlotSame2(Hcth_foam_Ceex2_vmax002,   ycapt, kRed,    +0.70, "(b)", "KKfoam, IFI on ");
  PlotSame2(Hcth_vTcPL_Ceex2n_vmax002, ycapt, kBlack,  -0.80, "(c)", "KKMC,   IFI off ");
  PlotSame2(Hcth_foam_Ceex2n_vmax002,  ycapt, kPine,   -0.70, "(d)", "KKfoam, IFI off ");
  //
  //PlotSame2(HST_cs_ceex2_vmax002,      ycapt, kPine,  -0.50, "(s)", "KKfoam2, IFIon ");

  //====================plot2========================
  cFigCosThe2->cd(2);
  TH1D *Hst_ratio1  = HstRatio("Hst_ratio1",   Hcth_vTcPL_Ceex2_vmax002,  Hcth_foam_Ceex2_vmax002, kBlack);
  TH1D *Hst_ratio2  = HstRatio("Hst_ratio2",   Hcth_vTcPL_Ceex2n_vmax002, Hcth_foam_Ceex2n_vmax002,kBlue);
  //
  TH1D *Hst_ratio5  = HstRatio("Hst_ratio5",   HST_cs_ceex2_vmax002,      Hcth_foam_Ceex2_vmax002, kBlack);
  Hst=Hst_ratio1;
  Hst->SetStats(0);
  Hst->SetTitle(0);
  Hst->GetXaxis()->SetTitle("cos #theta");

  Hst->SetMinimum(1.0 -0.01); Hst->SetMaximum(1.0 +0.01);
  Hst->DrawCopy("h");
  CaptT->DrawLatex(0.00,0.96,"Ratios of d#sigma/d(cos #theta)");
  ycapt = 0.85; // starting value, to be decremented below
  CaptT->DrawLatex(0.40,ycapt,gTextEne);  ycapt += -0.04;
  CaptT->DrawLatex(0.40,ycapt,gTextNev);  ycapt += -0.04;
  CaptT->DrawLatex(0.40,ycapt,gTextNev2); ycapt += -0.04;
  CaptT->DrawLatex(0.40,ycapt,"vmax = 0.02");
  //
  PlotSame2(Hst_ratio2,  ycapt, kBlue,     0.60, "(a)", "KKMC/KKfoam3   IFIoff");
  PlotSame2(Hst_ratio1,  ycapt, kBlack,    0.30, "(b)", "KKMC/KKfoam5   IFIon");
  //
  //PlotSame2(Hst_ratio5,  ycapt, kPine,    0.40, "(x)", "KKfoam2/KKfoam5 IFIon");

  hOne->SetLineColor(kBlack);
  hOne->DrawCopy("hsame");

  //================================================
  cFigCosThe2->SaveAs("cFigCosThe2.pdf");

}//FigCosThe2


///////////////////////////////////////////////////////////////////////////////////
void FigCtheSoft()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigCtheSoft =========================== "<<endl;
  // KKMC
  TH1D *Hcth_vTcPL_Ceex2n_vmax02= (TH1D*)DiskFileB.Get("Hcth_vTcPL_Ceex2n_vmax02");
  // KKfoam3 from scat
  TH1D *Hcth_foam3_Ceex2n_vmax02  = (TH1D*)DiskFileB.Get("Hcth_foam_Ceex2n_vmax02");
  TH1D *Hcth_foam3_Ceex2n_vmax002 = (TH1D*)DiskFileB.Get("Hcth_foam_Ceex2n_vmax002");
  // KKfoam3 from hst
  TH1D *HST_cc_EEX2_vmax02   = (TH1D*)DiskFileF.Get("HST_cc_EEX2_vmax02");
  TH1D *HST_cc_EEX2_vmax002  = (TH1D*)DiskFileF.Get("HST_cc_EEX2_vmax002");
  TH1D *HST_cc_EEX2_vmax0002 = (TH1D*)DiskFileF.Get("HST_cc_EEX2_vmax0002");
  // KKfoam2 SOFT!!!
  TH1D *HST_cs_EEX2_vmax02   = (TH1D*)DiskFileF2.Get("HST_cs_EEX2_vmax02");
  TH1D *HST_cs_EEX2_vmax002  = (TH1D*)DiskFileF2.Get("HST_cs_EEX2_vmax002");
  TH1D *HST_cs_EEX2_vmax0002 = (TH1D*)DiskFileF2.Get("HST_cs_EEX2_vmax0002");
  //
  TH1D *HST_cs_ceex2n_vmax02   = (TH1D*)DiskFileF2.Get("HST_cs_ceex2n_vmax02");
  TH1D *HST_cs_ceex2n_vmax002  = (TH1D*)DiskFileF2.Get("HST_cs_ceex2n_vmax002");

  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.04);
 ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigCtheSoft = new TCanvas("cFigCtheSoft","cFigCtheSoft", gXcanv, gYcanv,    1200, 600);
  //                                      Name    Title        xoff,yoff, WidPix,HeiPix
  gXcanv += 25, gYcanv += 25;
  cFigCtheSoft->SetFillColor(10);
  ////////////////////////////////////////////////////////////////////////////////
  cFigCtheSoft->Divide( 2,  0);
  //====================plot1========================
  cFigCtheSoft->cd(1);
  TH1D *Hst=Hcth_vTcPL_Ceex2n_vmax02;
  Hst->SetStats(0);
  Hst->SetTitle(0);
  Hst->GetXaxis()->CenterTitle();
  //Hst->GetYaxis()->SetTitleSize(0.04);
  //Hst->GetYaxis()->SetTitle("d#sigma/dcos(#theta) [nb]");
  Hst->GetXaxis()->SetTitleSize(0.04);
  Hst->GetXaxis()->SetTitle("cos #theta");

  Hst->SetMinimum(0);Hst->SetMaximum(0.8);

  Hst->SetLineColor(kBlue);
  Hst->SetMinimum(0.0);
  Hst->DrawCopy("h");

  CaptT->DrawLatex(0.00,0.96,"d#sigma/d(cos #theta)");
  double ycapt = 0.90;
  CaptT->DrawLatex(0.40, ycapt,gTextEne);
  PlotSame2(Hcth_vTcPL_Ceex2n_vmax02,  ycapt, kBlack,  -0.80, "(a1)", "KKMC:  IFIoff, v<0.02 ");
  //
  PlotSame2(HST_cc_EEX2_vmax02,        ycapt, kRed,   -0.40, "(c1)", "foam3: IFIoff, v<0.02  (hst)");
  PlotSame2(HST_cc_EEX2_vmax002,       ycapt, kRed,   -0.40, "(c2)", "foam3: IFIoff, v<0.002 (hst)");
  PlotSame2(HST_cc_EEX2_vmax0002,      ycapt, kRed,   -0.40, "(c3)", "foam3: IFIoff, v<0.0002(hst)");
  //
  PlotSame2(HST_cs_EEX2_vmax02,       ycapt, kPine,  -0.70, "(d1)", "foam2, IFIoff, v<0.02   soft");
  PlotSame2(HST_cs_EEX2_vmax002,      ycapt, kPine,  -0.70, "(d2)", "foam2, IFIoff, v<0.002  soft");
  PlotSame2(HST_cs_EEX2_vmax0002,     ycapt, kPine,  -0.70, "(d3)", "foam2, IFIoff, v<0.0002 soft");
//
  PlotSame2(Hcth_foam3_Ceex2n_vmax02,  ycapt, kGold,  -0.20, "(f1)", "foam3: ceex2n, v<0.02  (scat)");
  PlotSame2(Hcth_foam3_Ceex2n_vmax002, ycapt, kGold,  -0.20, "(f2)", "foam3: ceex2n, v<0.002 (scat)");
//
  PlotSame2(HST_cs_ceex2n_vmax02,      ycapt, kBlue,  -0.00, "(g1)", "foam2: ceex2n, v<0.02  soft");
  PlotSame2(HST_cs_ceex2n_vmax002,     ycapt, kBlue,  -0.00, "(g2)", "foam2: ceex2n, v<0.002 soft");

  //====================plot2========================
  cFigCtheSoft->cd(2);
  TH1D *Hst_rats1  = HstRatio("Hst_rats1",   HST_cc_EEX2_vmax02,   HST_cs_EEX2_vmax02,  kBlack);
  TH1D *Hst_rats2  = HstRatio("Hst_rats2",   HST_cc_EEX2_vmax002,  HST_cs_EEX2_vmax002, kBlue);
  TH1D *Hst_rats3  = HstRatio("Hst_rats3",   HST_cc_EEX2_vmax0002, HST_cs_EEX2_vmax0002,kPine);

  TH1D *Hst_rats10 = HstRatio("Hst_rats10", Hcth_foam3_Ceex2n_vmax02,  HST_cs_ceex2n_vmax02,  kBlack);
  TH1D *Hst_rats20 = HstRatio("Hst_rats20", Hcth_foam3_Ceex2n_vmax002, HST_cs_ceex2n_vmax002, kBlack);

  Hst=Hst_rats1;
  Hst->SetStats(0);
  Hst->SetTitle(0);
  Hst->SetMinimum(1.0 -0.01); Hst->SetMaximum(1.0 +0.01);

  Hst->DrawCopy("h");
  ycapt = 0.90;
  PlotSame2(Hst_rats1,    ycapt, kBlack,  -0.80, "(e1)", "foam2/3, eex, v<0.02  ");
  PlotSame2(Hst_rats2,    ycapt, kBlue,   -0.80, "(e2)", "foam2/3, eex, v<0.002 ");
  PlotSame2(Hst_rats3,    ycapt, kPine,  -0.80, "(e3)", "foam2/3, eex, v<0.0002");

  PlotSame2(Hst_rats10,   ycapt, kGold,   -0.60, "(v1)", "foam2/3 ceex2n, v<0.02  ");
  PlotSame2(Hst_rats20,   ycapt, kRed,    -0.60, "(v2)", "foam2/3 ceex2n, v<0.002 ");
  //================================================
  cFigCtheSoft->SaveAs("cFigCtheSoft.pdf");

}//FigCtheSoft


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
  TH1D *HST_FOAM_NORMA3 = (TH1D*)DiskFileF.Get("HST_FOAM_NORMA3");
  Nodes2   =  HST_FOAM_NORMA3->GetBinContent(511);    // No of farm nodes (trick)
  double  CMSeneF  = HST_FOAM_NORMA3->GetBinContent(1)/Nodes2; // CMSene=xpar(1)
  if( fabs(gCMSene/CMSeneF-1) >1e-4 ){
	  cout<<" +++++ Wrong input files !!!! KKMC "<< gCMSene <<"GeV and  FOAM "<< CMSeneF<<"GeV"<<endl;
	  exit(19);
  }
  gNevTot2  = HST_FOAM_NORMA3->GetEntries();       // MC statistics from KKMC
  sprintf(gTextNev2,"FOAM:%10.2e events", gNevTot2);

  //
  cout<< "CMSene[GeV] = "<< gCMSene<< endl;
  cout<< "KKMC: No. of farm nodes="<< Nodes  << "  Tot no. of events = "<<gNevTot<< endl;
  //cout<< "FOAM: No. of farm nodes="<< Nodes2 << "  Tot no. of events = "<<gNevTot2<<endl;
//////////////////////////////////////////////////////////////////////////
//
  HistNormalize();     // Renormalization of MC histograms
  ReMakeKKMC();        // reprocessing KKMC histos
  HisReMakeFoam35();   // prepare histos from KKfoam
//  KKsemMakeHisto();    // prepare histos from KKsem
  //========== PLOTTING ==========
  //
  FigScatV();
/*
  FigScatA();
  FigInfo();
  //FigVtest();  // introduct. tests/calibrations
  //FigCtest();  // introduct. tests/calibrations

  FigVprod();
  FigCprod();
  //
  FigCosThe();
  FigCosThe2();
  //
  //FigCtheSoft();

*/
  //++++++++++++++++++++++++++++++++++++++++
  DiskFileA.ls();
  DiskFileB.ls();
  DiskFileB.Write();
  DiskFileB.Close();
  //++++++++++++++++++++++++++++++++++++++++
  theApp.Run();
  //++++++++++++++++++++++++++++++++++++++++
}
