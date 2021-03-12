//////////////////////////////////////////////////////////////////////
// make TableNev-pdf
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
#include "KKplot.h"

//=============================================================================
//  ROOT  ROOT ROOT   ROOT  ROOT  ROOT  ROOT  ROOT  ROOT  ROOT   ROOT   ROOT
//=============================================================================
//
TFile  DiskFileA("../workKKMC/histo.root");           // KKMCee current

//
TFile DiskFileB("RhoSemi.root","RECREATE","histograms");
FILE *DiskFileTeX;

// Interface to KKplot and some extra plotting facilities
KKplot LibSem("KKplot");

///////////////////////////////////////////////////////////////////////////////////

void HistNormalize(){
  //
  cout<<"----------------------------- HistNormalize ------------------------------------"<<endl;
  DiskFileA.ls("");
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_axib1") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_axib2") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA.Get("hst_axib3") );
  //
  cout<<"----------------------------- HistNormalize ended -------------------------------"<<endl;
  //
}



///////////////////////////////////////////////////////////////////////////////////
void TableNew()
{
//------------------------------------------------------------------------
  cout<<" ========================= TableNew start=========================== "<<endl;
//

  TH1D *hst_axib1= (TH1D*)DiskFileA.Get("hst_axib1");
  TH1D *hst_axib2= (TH1D*)DiskFileA.Get("hst_axib2");
  TH1D *hst_axib3= (TH1D*)DiskFileA.Get("hst_axib3");

//  Char_t Capt[20][132];

// Column captions
  int nPlt=3;   // KORALZ eliminated
  Char_t *Capt[nPlt+1];
  for( int i=0; i<=nPlt; i++ ) Capt[i]=new char[132];
  strcpy(Capt[0],"{\\color{blue}$N$}");
  strcpy(Capt[1],"{\\color{blue} Tot.}");
  strcpy(Capt[2],"{\\color{blue} 1 phot. }");
  strcpy(Capt[3],"{\\color{blue} 1 phot. NO lept.}");
//  strcpy(Capt[4],"{\\color{red}${\\cal O}(\\alpha^2)_{\\rm CEEX}$ }");

// formats, not used in PlTable2
//  Char_t fmt[3][10];
//  strcpy(fmt[0],"f10.2");
//  strcpy(fmt[1],"f10.4");
//  strcpy(fmt[2],"f8.4");

// pointers to histograms
  TH1D *iHst[nPlt+1];
  iHst[1]= hst_axib1;  //
  iHst[2]= hst_axib2;  //
  iHst[3]= hst_axib3;  //
//  iHst[4]= HTot_vTcPR_Ceex2;   //CEEX2
  iHst[1]->Scale(1e3);    // nano- to pico-barns
  iHst[2]->Scale(1e3);    // nano- to pico-barns
  iHst[3]->Scale(1e3);    // nano- to pico-barns
//  iHst[4]->Scale(1e3);    // nano- to pico-barns
// multicolumn caption
  Char_t Mcapt[132];
//  strcpy(Mcapt,"{\\color{red} 161GeV, $\\mu^+\\mu^-\\gamma$, $\\sigma$[pb]}");
  strcpy(Mcapt,"{\\color{red} 91.2GeV, $\\nu\\bar{\\nu}\\gamma$, $\\sigma$[pb]}");

///************************************
  DiskFileTeX = fopen("TableNew.txp","w");
//************************************
// Initialization of the latex source file
  PlInitialize(DiskFileTeX, 2);


//  int k1,k2,dk;
//  k1=10; k2=90; dk=20;  //original
//  k1= 5; k2=45; dk=10;
    PlTable2( nPlt, iHst, DiskFileTeX, Capt,  Mcapt, " ", 1, 1, 1); // for 1 bins
//  PlTable2( nPlt, iHst, DiskFileTeX, Capt,  Mcapt, "B", 1, 1, 1); // for 1 bins
//  PlTable2(-nPlt, iHst, DiskFileTeX, Capt,  Mcapt, "T", 1, 1, 1); // for 1 bins
//  PlTable2(-nPlt, iHst, DiskFileTeX, Capt,  Mcapt, "E", 1, 1, 1); // for 1 bins

// finalizing latex source file
  PlEnd(DiskFileTeX);
//************************************
  fclose(DiskFileTeX);
//************************************
  cout<<" ========================= TableNew end =========================== "<<endl;
}//TableNew



///////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  //++++++++++++++++++++++++++++++++++++++++
  TApplication theApp("theApp", &argc, argv);
  //++++++++++++++++++++++++++++++++++++++++
//  LibSem.Initialize(DiskFileA);
/////////////////////////////////////////////////////////////////////////
  DiskFileB.cd();
  HistNormalize();     // Renormalization of MC histograms
//  ReMakeMChisto();     // reprocessing MC histos
//  KKsemMakeHisto();    // prepare histos from KKsem
  //========== PLOTTING ==========
  // Some comparisons with KKsem
  // Old benchmarks KKMC vs. KKsem with Gauss integrator
  TableNew();
  //++++++++++++++++++++++++++++++++++++++++
  DiskFileA.ls();
  DiskFileB.ls();
  DiskFileB.Write();
  DiskFileB.Close();

  //++++++++++++++++++++++++++++++++++++++++
  //  theApp.Run();  // only for displaying plots
  //++++++++++++++++++++++++++++++++++++++++
}

