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

FILE *DiskFileTeX;

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
  HisNorm2(HST_FOAM_NORMA4, (TH2D*)DiskFileF->Get("SCA_vTcPR_Eex2") );
//------------------------------
  TH1D    *HST_FOAM_NORMA6= (TH1D*)DiskFileF->Get("HST_FOAM_NORMA6");

}// HistNormalize
//////////////////////////////////////////////////////////////////////////////////
void ReMakeMChisto(){
	// Some KKMC histos are preprocessed
	//------------------------------------------------------------------------
  cout<<"==================================================================="<<endl;
  cout<<"================ ReMakeMChisto  BEGIN  ============================"<<endl;
  TH1D *HST_KKMC_NORMA  = (TH1D*)DiskFileA->Get("HST_KKMC_NORMA");
  TH1D *HST_FOAM_NORMA4 = (TH1D*)DiskFileF->Get("HST_FOAM_NORMA4");
  //****************************************************************************************
  // KKMCee reprocessing part
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
void TableIII()
{
//------------------------------------------------------------------------
  cout<<" ========================= TableIII start=========================== "<<endl;
//
  Double_t CMSene;
//  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
//  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
//  CMSene /= HST_KKMC_NORMA->GetBinContent(511); // farm adjusted
//  char TextEne[100]; sprintf(TextEne,"#sqrt{s} =%4.2fGeV", CMSene);
  //
  // Distributions of sigma(vvmax) without cutoff on c=cos(thetaPRD)
  TH1D *HTot_vTcPR_Eex2  = (TH1D*)DiskFileB->Get("HTot_vTcPR_Eex2"); // KKeeFoam
  TH1D *HAfb_vTcPR_Eex2  = (TH1D*)DiskFileB->Get("HAfb_vTcPR_Eex2"); // KKeeFoam
  //
  // Distributions of sigma(vvmax) without cutoff on c=cos(thetaPRD)
  TH1D *hTot_vTcPR_Eex2  = (TH1D*)DiskFileB->Get("hTot_vTcPR_Eex2"); // KKMC
  TH1D *hAfb_vTcPR_Eex2  = (TH1D*)DiskFileB->Get("hAfb_vTcPR_Eex2"); // KKMC

//  Char_t Capt[20][132];

// Column captions
  int nPlt=2;   //
  Char_t *Capt[nPlt+1];
  for( int i=0; i<=nPlt; i++ ) Capt[i]=new char[132];
  strcpy(Capt[0],"{\\color{blue}$v_{\\max}$}");
  strcpy(Capt[1],"{\\color{blue} KKeeFoam ${\\cal O}(\\alpha^2)$}");
  strcpy(Capt[2],"{\\color{blue} KKMCee ${\\cal O}(\\alpha^2)_{\\rm EEX2}$ }");
//  strcpy(Capt[3],"{\\color{red}${\\cal O}(\\alpha^2)_{\\rm CEEX}$ intOFF}");
//  strcpy(Capt[4],"{\\color{red}${\\cal O}(\\alpha^2)_{\\rm CEEX}$ }");

// formats, not used in PlTable2
//  Char_t fmt[3][10];
//  strcpy(fmt[0],"f10.2");
//  strcpy(fmt[1],"f10.4");
//  strcpy(fmt[2],"f8.4");

// pointers to histograms
  TH1D *iHst[nPlt+1];
  iHst[1]= HTot_vTcPR_Eex2;     //KKeeFoam
  iHst[2]= hTot_vTcPR_Eex2;     //KKMCee EEX2
//  iHst[3]= HTot_vTcPR_Ceex2n;  //CEEX2 INT off
//  iHst[4]= HTot_vTcPR_Ceex2;   //CEEX2
  iHst[1]->Scale(1e3);    // nano- to pico-barns
  iHst[2]->Scale(1e3);    // nano- to pico-barns
//  iHst[3]->Scale(1e3);    // nano- to pico-barns
//  iHst[4]->Scale(1e3);    // nano- to pico-barns
// multicolumn caption
  Char_t Mcapt[132];
  strcpy(Mcapt,"{\\color{red}$\\sigma(v_{\\max})$ [pb]}");

///************************************
  DiskFileTeX = fopen("TabIIIPRD.txp","w");
//************************************
// Initialization of the latex source file
  PlInitialize(DiskFileTeX, 2);


//  int k1,k2,dk;
//  k1=10; k2=90; dk=20;  //original
//  k1= 5; k2=45; dk=10;
  PlTable2( nPlt, iHst, DiskFileTeX, Capt,  Mcapt, "B", 2, 2, 2); // for 100 bins
  PlTable2(-nPlt, iHst, DiskFileTeX, Capt,  Mcapt, "T",10,90,20); // for 100 bins
  PlTable2(-nPlt, iHst, DiskFileTeX, Capt,  Mcapt, "T",99,99, 2); // for 100 bins
  iHst[1]->Scale(1e-3);    // back to nano-barns
  iHst[2]->Scale(1e-3);    //
//  iHst[3]->Scale(1e-3);    //
//  iHst[4]->Scale(1e-3);    //

  iHst[1]= HAfb_vTcPR_Eex2;    //KKeeFoam
  iHst[2]= hAfb_vTcPR_Eex2;    //KKMCee EEX2
  //iHst[3]= HAfb_vTcPR_Ceex2n;  //CEEX2 INT off
  //iHst[4]= HAfb_vTcPR_Ceex2;   //CEEX2

  strcpy(Mcapt,"{\\color{red}$A_{\\rm FB}(v_{\\max})$}");
  PlTable2( nPlt, iHst, DiskFileTeX, Capt,  Mcapt, "T", 2, 2, 2); // for 100 bins
  PlTable2(-nPlt, iHst, DiskFileTeX, Capt,  Mcapt, "T",10,90,20); // for 100 bins
  PlTable2(-nPlt, iHst, DiskFileTeX, Capt,  Mcapt, "E",99,99, 2); // for 100 bins

// finalizing latex source file
  PlEnd(DiskFileTeX);
//************************************
  fclose(DiskFileTeX);
//************************************
  cout<<" ========================= TableIII end =========================== "<<endl;
}//TableIII


///////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{

  gSystem->Load("../SRCee/.libs/libKKee.so");     // needed ???
  gSystem->Load("../SRCee/.libs/libKKfm.so");     // NEEDED! why?

  DiskFileA = new TFile(FileA);
  DiskFileF = new TFile(FileF);
//
  DiskFileB = new TFile("PlotIIIPRD.root","RECREATE","histograms");

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
  TableIII();
 //++++++++++++++++++++++++++++++++++++++++
  DiskFileA->ls();
//
//  cout<< "CMSene[GeV] = "<< gCMSene<< endl;
//  cout<< "KKMC: No. of farm nodes="<< Nodes  << "  Tot no. of events = "<<gNevTot<< endl;
//
  //++++++++++++++++++++++++++++++++++++++++
//  theApp.Run();
  //++++++++++++++++++++++++++++++++++++++++
}
