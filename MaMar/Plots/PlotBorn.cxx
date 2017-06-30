// make PlotBorn-run
//////////////////////////////////////////////////////////
//line-shape i Afb, with and without box corrections
//    w funkcji sqrt(s), zakres 70 - 150 GeV
//    fixed costheta = 0.0, 0.3, 0.6, 0.9
//    procesy: ee-> mumu, uubar, ddbar
//////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
using namespace std;

#include <math.h>
#include "TROOT.h"
#include "TCanvas.h"
#include "TF2.h"
#include "TH2.h"
#include "TLatex.h"
#include "TApplication.h"
#include "TObjString.h"
#include "TFile.h"

#include "HisNorm.h"
#include "KKplot.h"

KKplot LibSem("kkplot");

//
TFile DiskFileB("RhoSemi.root","RECREATE","histograms");
FILE *DFile;


void Vdef(double v[4], const double v1, const double v2, const double v3, const double v4)
  { // define a 4-vector (avoids initialization warnings)
       v[0] = v1; v[1] = v2; v[2] = v3; v[3] = v4;
  }


void InitBorn(){
cout<<"==========================================================="<<endl;
cout<<"================ InitBorn BEGIN ============================"<<endl;

//////////////////////////////////////////////////////////////
//   Initialize MC generator and analysis programs          //
//////////////////////////////////////////////////////////////
double   m_xpar[10001];      // complete input of KKMC run
const int jmax =10000;
//ReaData("./KK2f_defaults_ERW",    jmax, m_xpar);  // numbering as in input!!!
LibSem.ReaData("../../.KK2f_defaults",     jmax, m_xpar);  // numbering as in input!!!
// User input data
//LibSem.ReaData("../workKKMC/workKKMC_189GeV.input", -jmax, m_xpar);  // jmax<0 means no-zeroing
LibSem.ReaData("../workKKMC/PlotBorn_189GeV.input", -jmax, m_xpar);  // jmax<0 means no-zeroing
double ypar[jmax];
for(int j=0;j<jmax;j++) ypar[j]=m_xpar[j+1];    // ypar has c++ numbering
// KKMC and KKsem initialization
const char *output_file = "./kksem.output";
long stl2 = strlen(output_file);
long mout =16;
//kk2f_fort_open_(mout,output_file,stl2);
//kk2f_initialize_(ypar);
//kksem_initialize_(ypar);
LibSem.Initialize(m_xpar);

//************************************
DFile = fopen("TableBorn.txt","w");
//************************************

cout<<"================ InitBorn END   ============================"<<endl;

}//InitBorn


void TabBorn(){
cout<<"==========================================================="<<endl;
cout<<"================ TabBorn BEGIN ============================"<<endl;
//  ************* user histograms  *************
double svar2, sigma,  CosTheta;
double dSig0, dSig3, dSig6, dSig9;
int m_KFini = 11;
int m_KFf   = 13;
fprintf(DFile," e+e- --> mu+ mu- \n");
fprintf(DFile," d(sigma)/d(cos_theta) [nb], cos_theta = 0.0, 0.3, 0.6, 0.9 \n");
double Emin = 88;
double Emax = 94;
int    nPt  = 10;
double Ene;
//
for( int i=0; i<=nPt; i++ ){
   Ene = Emin +i*((Emax-Emin)/nPt);
   svar2 = Ene*Ene;
   CosTheta=0; bornv_interpogsw_(m_KFf,svar2, CosTheta);
   dSig0 = bornv_dizet_( 1, m_KFini, m_KFf, svar2, CosTheta, 0.0, 0.0, 0.0, 0.0);
   //
   CosTheta=0.3; bornv_interpogsw_(m_KFf,svar2, CosTheta);
   dSig3 = bornv_dizet_( 1, m_KFini, m_KFf, svar2, CosTheta, 0.0, 0.0, 0.0, 0.0);
   //
   CosTheta=0.6; bornv_interpogsw_(m_KFf,svar2, CosTheta);
   dSig6 = bornv_dizet_( 1, m_KFini, m_KFf, svar2, CosTheta, 0.0, 0.0, 0.0, 0.0);
   //
   CosTheta=0.9; bornv_interpogsw_(m_KFf,svar2, CosTheta);
   dSig9 = bornv_dizet_( 1, m_KFini, m_KFf, svar2, CosTheta, 0.0, 0.0, 0.0, 0.0);
   fprintf(DFile,"Ene= %10.5f  dSig_dCos= %12.7f  %12.7f  %12.7f  %12.7f  \n", Ene, dSig0, dSig3, dSig6, dSig9);
   cout<<"Ene="<<Ene; cout<< "  Sig_dCos= "<< dSig0 <<"  "<< dSig3 <<"  "<< dSig6 <<"  "<< dSig9<<endl;
}//i

cout<<"==========================================================="<<endl;

long KeyFob;
KeyFob=   10; // BornV_Dizet, with EW and without integration ???
KeyFob=  -11; // BornV_Simple, for KeyLib=0, NO EW, NO integration OK
kksem_setkeyfob_( KeyFob );
double xBorn;
for( int i=0; i<=nPt; i++ ){
   Ene = Emin +i*((Emax-Emin)/nPt);
   svar2 = Ene*Ene;
   kksem_makeborn_( svar2, xBorn);
   cout<< "Ene= "<<Ene<< "xBorn [nb]= "<<xBorn<<endl;
}

//************************************
  fclose(DFile);
//************************************

cout<<"================ TabBorn END   ============================"<<endl;
}//TabBorn


void FigBorn3(){
cout<<"==========================================================="<<endl;
cout<<"================ FigBorn3 BEGIN ==========================="<<endl;
double Emin = 85;
double Emax = 95;
int    nPt  = 200;
double Ene;
TH1D *hst_sigma3 = new TH1D("hst_sigma3" ,  "sigma(CMSene) [nb]",    nPt, Emin, Emax);
hst_sigma3->Sumw2();

long KeyFob;
KeyFob=   10; // BornV_Dizet, with EW and without integration ???
KeyFob=  -11; // BornV_Simple, for KeyLib=0, NO EW, NO integration OK
KeyFob=    0; // BornV_Dizet, with EW and with integration
kksem_setkeyfob_( KeyFob );
double xBorn, svar2;
for(int ix=1; ix <= nPt; ix++){
   Ene = Emin +(ix-1)*((Emax-Emin)/nPt);
   svar2 = Ene*Ene;
   kksem_makeborn_( svar2, xBorn);
   hst_sigma3->SetBinContent(  ix, xBorn );
   hst_sigma3->SetBinError(    ix, 0.0 );
   cout<< "Ene= "<<Ene<< "xBorn [nb]= "<<xBorn<<endl;
}//ix
//------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////
TCanvas *cFigBorn3 = new TCanvas("cFigBorn3","cFigBorn3: sigma(CMSene) ", 50, 100,    500,  500);
//                                      Name    Title                  xoff,yoff, WidPix,HeiPix
cFigBorn3->SetFillColor(10);
TH1D *Hst1 = hst_sigma3;
cFigBorn3->cd();
Hst1->SetStats(0);
//Hst1->SetTitle(0);
Hst1->SetMinimum(0);
Hst1->SetLineColor(kRed);
Hst1->DrawCopy("h");
cout<<"================ FigBorn3 END   ==========================="<<endl;
}//FigBorn3



void FigBorn1(){
cout<<"==========================================================="<<endl;
cout<<"================ FigBorn1 BEGIN ==========================="<<endl;
double CosTheta = 0; // for some reason 0 no allowed!!!
CosTheta = 0.01;
double Emin = 85;
double Emax = 95;
Emin=10;
int    nPt  = 200;
double CMSene;
//
double m_p1[4];         //!
double m_p2[4];         //!
double m_p3[4];         //!
double m_p4[4];         //!
//
TH1D *hst_sigEEX = new TH1D("hst_sigEEX" ,  "dSigma(CMSene,costhe)",    nPt, Emin, Emax);
hst_sigEEX->Sumw2();
TH1D *hst_sigGPS = new TH1D("hst_sigGPS" ,  "dSigma(CMSene,costhe)",    nPt, Emin, Emax);
hst_sigGPS->Sumw2();

long KeyFob;
KeyFob=   10; // BornV_Dizet, with EW and without integration ???
KeyFob=  -11; // BornV_Simple, for KeyLib=0, NO EW, NO integration OK
KeyFob=    0; // BornV_Dizet, with EW and with integration
kksem_setkeyfob_( KeyFob );
int KFf =13, KFini=11;
double xBorn, svar2;
for(int ix=1; ix <= nPt; ix++){
   CMSene = Emin +(ix-1)*((Emax-Emin)/nPt);
   svar2 = CMSene*CMSene;
   //
   bornv_interpogsw_(KFf,svar2, CosTheta);
   double dSig_EEX = bornv_dizet_( 1, KFini, KFf, svar2, CosTheta, 0.0, 0.0, 0.0, 0.0);
   cout<< "CMSene= "<<CMSene<< "   dSig_EEX= "<<dSig_EEX<<endl;
   //
   // =============== Sigm/dOmega from spin amplitudes ===============
   // Effective 4-momenta, KKMC convention: p={px,py,pz,E)
   double Ene   = CMSene/2;
   double   m_beam  = 0.510999e-3;
   double	m_fin   = 0.1056583;
   double Pmb  = sqrt( (Ene-m_beam)*(Ene+m_beam) ); // modulus
   Vdef(m_p1, 0, 0 , Pmb, Ene);  // beam
   Vdef(m_p2, 0, 0 ,-Pmb, Ene);  // beam
   double Pmf  =sqrt( (Ene-m_fin)*(Ene+m_fin) ); // modulus
   Vdef(m_p3, Pmf*sqrt(1-sqr(CosTheta)), 0 , Pmf*CosTheta,  Ene); // final
   Vdef(m_p4,-Pmf*sqrt(1-sqr(CosTheta)), 0 ,-Pmf*CosTheta,  Ene); // final
   double PX[4] = {0, 0, 0, 2*Ene};
   //***** pure Born of CEEX
   double dSig_GPS;
   gps_bornf_(KFini, KFf ,PX, CosTheta, m_p1,m_beam, m_p2, -m_beam,
                                        m_p3,m_fin,  m_p4, -m_fin,   dSig_GPS);
   /////////////////////////////////////////////////////////////////
    cout<< "CMSene= "<<CMSene<< "  dSig_GPS= "<<dSig_GPS<<"  ratio="<< dSig_GPS/dSig_EEX  <<endl;

    hst_sigEEX->SetBinContent(  ix, dSig_EEX );
    hst_sigEEX->SetBinError(    ix, 0.0 );
    hst_sigGPS->SetBinContent(  ix, dSig_GPS );
    hst_sigGPS->SetBinError(    ix, 0.0 );

}//ix

//------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////
TCanvas *cFigBorn1 = new TCanvas("cFigBorn1","cFigBorn1: sigma(CMSene) ", 250, 200,    900,  600);
//                                      Name    Title                  xoff,yoff, WidPix,HeiPix
cFigBorn1->SetFillColor(10);
TPad*    UpperPad = new TPad("upperPad", "upperPad",
			       .005, .3525, .995, .995);
TPad*    LowerPad = new TPad("lowerPad", "lowerPad",
			       .005, .005, .995, .3575);
UpperPad->Draw();
LowerPad->Draw();
/////////////////////////////////////////////////////
UpperPad->cd();
TH1D *Hst1 = hst_sigEEX;
Hst1->SetStats(0);
//Hst1->SetTitle(0);
Hst1->SetMinimum(0);
Hst1->SetLineColor(kRed);
Hst1->DrawCopy("h");
hst_sigGPS->DrawCopy("hsame");

TLatex *CaptT = new TLatex(); CaptT->SetNDC(); // !!!
char TextCos[100]; sprintf(TextCos,"cos#theta =%6.2f ", CosTheta);
CaptT->DrawLatex(0.60,0.75,TextCos);

////////////////////////////////////////////////////
LowerPad->cd();

TH1D *dSigRatio =(TH1D*)hst_sigGPS->Clone("dSigRatio");
dSigRatio->Divide(dSigRatio, hst_sigEEX); // sigma(xmax), FoamEEX/KKsem IFIoff
dSigRatio->SetLineColor(kRed);                 // FoamEEX2/KKsem IFIoff

dSigRatio->GetYaxis()->SetLabelSize(0.08); // default is 0.04
dSigRatio->GetXaxis()->SetLabelSize(0.08); // default is 0.04

dSigRatio->SetMaximum(1+0.0003);
dSigRatio->SetMinimum(1-0.0003);
dSigRatio->SetStats(0);
dSigRatio->SetTitle(0);
dSigRatio->DrawCopy("h");


TH1D *hOne = (TH1D*)dSigRatio->Clone("hOne");  // unity line
for(int i=1; i <= hOne->GetNbinsX() ; i++) { hOne->SetBinContent(i, 1); hOne->SetBinError(i, 0);}
hOne->SetLineColor(kRed);
hOne->DrawCopy("hsame");


cout<<"================ FigBorn1 END   ==========================="<<endl;
}//FigBorn1






///////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  //++++++++++++++++++++++++++++++++++++++++
  TApplication theApp("theApp", &argc, argv);
  //++++++++++++++++++++++++++++++++++++++++
  //
  DiskFileB.cd();
  //
  InitBorn();
  //
  TabBorn();
  //
  FigBorn3();
  FigBorn1();
  //
  //++++++++++++++++++++++++++++++++++++++++
  DiskFileB.ls();
  DiskFileB.Write();
  DiskFileB.Close();

  //++++++++++++++++++++++++++++++++++++++++
  theApp.Run();
  //++++++++++++++++++++++++++++++++++++++++
}

