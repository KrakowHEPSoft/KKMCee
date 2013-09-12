//////////////////////////////////////////////////////////////////////
//    make SkyLine
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

// ROOT headers
#include "TROOT.h"
#include "TFile.h"


//=============================================================================
//  ROOT  ROOT ROOT   ROOT  ROOT  ROOT  ROOT  ROOT  ROOT  ROOT   ROOT   ROOT 
//=============================================================================
TROOT tuples("rmain","Main program ");
TFile DiskFile("Rplots.root","RECREATE","histograms");
//=============================================================================

Double_t sqr( const Double_t x ){ return x*x;};

///////////////////////////////////////////////////////////////////////////////
//extern "C" void vacpol2001_(const long&, const double&, const double&, double&, double& );

//      DOUBLE PRECISION FUNCTION RRes_R_TOT(roots)
extern "C" double rres_r_tot_(const double&);
//      DOUBLE PRECISION FUNCTION RRes_Rqq(kf,roots)
extern "C" double rres_rqq_(const long&, const double&);
//      DOUBLE PRECISION FUNCTION RRes_F_PI_SQ(s)
extern "C" double rres_f_pi_sq_(const double&);
//      DOUBLE PRECISION  FUNCTION F_Pi_Kuehn_SQ(s)
//extern "C" double f_pi_kuehn02_sq_(const double&);
extern "C" double rres_f_pi_kuehn02_sq_(const double&);

///////////////////////////////////////////////////////////////////////////////////
Double_t Rtot(Double_t *xvar, Double_t *par)
// R total from all quarks
{
  long      KFfin  = par[0];
  Double_t  Mqq  = xvar[0];
  Double_t  Result;
  Result = rres_r_tot_(Mqq);
  return Result;
}
///////////////////////////////////////////////////////////////////////////////////
Double_t Rqq(Double_t *xvar, Double_t *par)
// R from one quark pair
{
  long      KFf  = par[0];
  Double_t  Mqq  = xvar[0];
  Double_t  Result =0.0;
  if(KFf>0 && KFf<6){
    Result = rres_rqq_(KFf,Mqq);
  }else if(KFf==12){
    Result = rres_rqq_(1,Mqq)+rres_rqq_(2,Mqq);
  }else if(KFf==123){
    Result = rres_rqq_(1,Mqq)+rres_rqq_(2,Mqq)+rres_rqq_(3,Mqq);
  }else if(KFf==1234){
    Result = rres_rqq_(1,Mqq)+rres_rqq_(2,Mqq)+rres_rqq_(3,Mqq)+rres_rqq_(4,Mqq);
  }else if(KFf==12345){
    Result = rres_rqq_(1,Mqq)+rres_rqq_(2,Mqq)+rres_rqq_(3,Mqq)
            +rres_rqq_(4,Mqq)+rres_rqq_(5,Mqq);
  }else if(KFf==99){
    Mqq  = sqrt(xvar[0]);
    Result = rres_rqq_(1,Mqq)+rres_rqq_(2,Mqq)+rres_rqq_(3,Mqq);;
  }else if(KFf==3001){
    Result = rres_f_pi_sq_(Mqq);
  }else if(KFf==3002){
    Result = rres_f_pi_kuehn02_sq_(Mqq);
  }
  return Result;
}

///////////////////////////////////////////////////////////////////////////////////
int figRtot()
{
  cout<<" ========================= figRtot =========================== "<<endl;
  //------------------------------------------------------------------------
  // parameters for ploting function
  Int_t npar =1;
  Double_t params[9];
  params[0]=1;       // type of quark
  Double_t mas1,mas2;
  int nbin;
  //------------------------------------------------------------------------
  mas1 =  0.000; mas2 = 12.000; nbin = 12000;
  TF1  *funRtot  = new  TF1("funRtot",Rtot, mas1,mas2,npar);
  funRtot->SetParameters(params);
  TH1F *hisRtot  = new TH1F("hisRtot","",  nbin, mas1,mas2);
  hisRtot->Eval(funRtot,"a");
  //------------------------------------------------------------------------
  mas1 =  0.000; mas2 = 12.000; nbin = 12000;
  TF1  *funRtot2  = new  TF1("funRtot",Rqq, mas1,mas2,npar);
  params[0]=12345;       // type of quark duscb
  funRtot2->SetParameters(params);
  TH1F *hisRtot2  = new TH1F("hisRtot2","",  nbin, mas1,mas2);
  hisRtot2->Eval(funRtot2,"a");
  //------------------------------------------------------------------------
  TF1  *funRdu  = new  TF1("funRdu",Rqq, mas1,mas2,npar);
  params[0]=12;       // type of quark du
  funRdu->SetParameters(params);
  TH1F *hisRdu  = new TH1F("hisRdu","",  nbin, mas1,mas2);
  hisRdu->Eval(funRdu,"a");
  //------------------------------------------------------------------------
  TF1  *funRdus  = new  TF1("funRdus",Rqq, mas1,mas2,npar);
  params[0]=123;       // type of quark du
  funRdus->SetParameters(params);
  TH1F *hisRdus  = new TH1F("hisRdus","",  nbin, mas1,mas2);
  hisRdus->Eval(funRdus,"a");
  //------------------------------------------------------------------------
  TF1  *funRdusc  = new  TF1("funRdusc",Rqq, mas1,mas2,npar);
  params[0]=1234;       // type of quark du
  funRdusc->SetParameters(params);
  TH1F *hisRdusc  = new TH1F("hisRdusc","",  nbin, mas1,mas2);
  hisRdusc->Eval(funRdusc,"a");
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  Float_t  WidPix, HeiPix;
  WidPix = 1200; HeiPix =  800;
  TCanvas *cRtot = new TCanvas("cRtot","Rtot grand view",150,70, WidPix,HeiPix);
  cRtot->SetFillColor(10);
  //cRtot->Divide(1, 2, 0.0,  0);
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  ///////////////////////////////////////////////////////////////////////////////
  cRtot->cd();
  hisRtot->GetYaxis()->SetTitle("R");
  hisRtot->GetYaxis()->CenterTitle();
  hisRtot->GetYaxis()->SetTitleSize(0.07);
  hisRtot->GetYaxis()->SetTitleOffset(0.5);
  hisRtot->GetXaxis()->SetTitle("#sqrt{s}  [GeV]");
  hisRtot->GetXaxis()->CenterTitle();
  hisRtot->GetXaxis()->SetTitleSize(0.05);
  // Rtotal duscb
  hisRtot->SetMaximum(10.0);
  hisRtot->SetStats(0);
  hisRtot->SetLineColor(4);    // blue line
  hisRtot->SetFillColor(4);    // blues fill solid
  hisRtot->DrawCopy();
  // dusc
  hisRdusc->SetLineColor(3);   // green line
  hisRdusc->SetFillColor(3);   // green fill solid
  hisRdusc->DrawCopy("same");
  // dus
  hisRdus->SetLineColor(2);    // red line
  hisRdus->SetFillColor(2);    // red fill solid
  hisRdus->DrawCopy("same");
  // du
  hisRdu->SetFillStyle(3004);  // transparent "////"  fill
  hisRdu->SetFillStyle(3010);  // transparent brickwall fill
  hisRdu->SetLineColor(1);     // black line
  hisRdu->SetFillColor(1);     // black fill solid
  hisRdu->DrawCopy("same");
  ///////////////////////////////////////////////////////////////////////////////
  cRtot->cd();
  cRtot->Update();
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
}

///////////////////////////////////////////////////////////////////////////////////
int figRres()
{
  cout<<" ========================= RERONANCES =========================== "<<endl;
  //------------------------------------------------------------------------
  // parameters for ploting function
  Int_t npar =1;
  Double_t params[9];
  params[0]=1;       // type of quark
  Double_t mas1,mas2;
  int nbin;
  //------------------------------------------------------------------------
  mas1 =  3.000; mas2 = 4.700;  nbin = 10000; // ccbar threshold
  mas1 =  9.200; mas2 = 11.500; nbin = 10000; // bbbar threshold
  mas1 =  0.600; mas2 = 1.200;  nbin = 10000; // rho region
  TF1  *funRres  = new  TF1("funRres",Rqq, mas1,mas2,npar);
  params[0]=12345;       // type of quark duscb
  funRres->SetParameters(params);
  TH1F *hisRres  = new TH1F("hisRres","",  nbin, mas1,mas2);
  hisRres->Eval(funRres,"a");
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  Float_t  WidPix, HeiPix;
  WidPix = 1100; HeiPix =  800;
  TCanvas *cRres = new TCanvas("cRres","Resonances",170,50, WidPix,HeiPix);
  //cRres->Divide(1, 2, 0.0,  0);
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  ///////////////////////////////////////////////////////////////////////////////
  cRres->SetFillColor(10);
  cRres->cd();
  //
  cRres->SetLogy(); // logarithmic scale
  //
  hisRres->GetYaxis()->SetTitle("R");
  hisRres->GetYaxis()->CenterTitle();
  hisRres->GetYaxis()->SetTitleSize(0.07);
  hisRres->GetYaxis()->SetTitleOffset(0.5);
  hisRres->GetXaxis()->SetTitle("#sqrt{s}  [GeV]");
  hisRres->GetXaxis()->CenterTitle();
  hisRres->GetXaxis()->SetTitleSize(0.05);
  // Rresal duscb
  //hisRres->SetMaximum(22.0);
  hisRres->SetStats(0);
  hisRres->SetLineColor(1);    // black line
  hisRres->SetFillColor(3);    // green fill solid (psi)
  hisRres->SetFillColor(4);    // blue fill solid (upsyl)
  hisRres->SetFillColor(2);    // red fill solid (rho)
  hisRres->DrawCopy();
  ///////////////////////////////////////////////////////////////////////////////
  cRres->Update();
  cRres->cd();
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
}


///////////////////////////////////////////////////////////////////////////////////
int figKloe()
{
  cout<<" ========================= RERONANCES =========================== "<<endl;
  //------------------------------------------------------------------------
  // parameters for ploting function
  Int_t npar =1;
  Double_t params[9];
  params[0]=1;       // type of quark
  Double_t s1,s2;
  int nbin;
  //------------------------------------------------------------------------
  nbin = 430; s1=0.100; s2=0.960; // rho region KLOE
  //-----
  TF1  *funKloe  = new  TF1("funKloe",Rqq, s1,s2,npar);
  params[0]=99;       // type of quark duscb
  funKloe->SetParameters(params);
  TH1F *hisKloe  = new TH1F("hisKloe","",  nbin, s1,s2);
  hisKloe->Eval(funKloe,"a");
  //-----
  TF1  *funFpiM  = new  TF1("funFpiM",Rqq, s1,s2,npar);
  params[0]=3001;       // Fpi^2 of marteen
  funFpiM->SetParameters(params);
  TH1F *hisFpiM  = new TH1F("hisFpiM","",  nbin, s1,s2);
  hisFpiM->Eval(funFpiM,"a");
  //-----
  TF1  *funFpiK  = new  TF1("funFpiK",Rqq, s1,s2,npar);
  params[0]=3002;       // Fpi^2 of marteen
  funFpiK->SetParameters(params);
  TH1F *hisFpiK  = new TH1F("hisFpiK","",  nbin, s1,s2);
  hisFpiK->Eval(funFpiK,"a");
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  Float_t  WidPix, HeiPix;
  WidPix = 800; HeiPix =  600;
  TCanvas *cKloe = new TCanvas("cKloe","Resonances",200,80, WidPix,HeiPix);
  //cKloe->Divide(1, 2, 0.0,  0);
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  ///////////////////////////////////////////////////////////////////////////////
  cKloe->SetFillColor(10);
  cKloe->cd();
  //
  //cKloe->SetLogy(); // logarithmic scale
  //
  hisKloe->GetYaxis()->SetTitle("R");
  hisKloe->GetYaxis()->CenterTitle();
  hisKloe->GetYaxis()->SetTitleSize(0.07);
  hisKloe->GetYaxis()->SetTitleOffset(0.5);
  hisKloe->GetXaxis()->SetTitle("Q^{2}  [GeV^{2}]");
  hisKloe->GetXaxis()->CenterTitle();
  hisKloe->GetXaxis()->SetTitleSize(0.05);
  // Kloeal duscb
  //hisKloe->SetMaximum(22.0);
  hisKloe->SetStats(0);
  hisKloe->SetLineColor(1);    // black line
  //hisKloe->SetFillColor(3);    // green fill solid (psi)
  //hisKloe->SetFillColor(4);    // blue fill solid (upsyl)
  //hisKloe->SetFillColor(2);    // red fill solid (rho)
  //
  hisKloe->DrawCopy();
  //-------------------------------------------
  //          Test of F_pi formfactor 
  //-------------------------------------------
  //funFpiM->DrawCopy("");
  //funFpiK->SetLineColor(2);    // red line
  //funFpiK->DrawCopy("same");
  ///////////////////////////////////////////////////////////////////////////////
  //double s=0.55;
  //double fpi1= rres_f_pi_sq_(s);
  //double fpi2= rres_f_pi_kuehn_sq_(s);
  //cout<<" fpi_marteen = " <<fpi1<<"  fpi_kuehn = "<<fpi2<<endl;
  ///////////////////////////////////////////////////////////////////////////////
  cKloe->Update();
  cKloe->cd();
  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
}


///////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  //++++++++++++++++++++++++++++++++++++++++
  TApplication theApp("theApp", &argc, argv);
  //++++++++++++++++++++++++++++++++++++++++
  figRtot();
  figRres();
  figKloe();
  //++++++++++++++++++++++++++++++++++++++++
  DiskFile.Write();
  DiskFile.ls();
  DiskFile.Close();
  //++++++++++++++++++++++++++++++++++++++++
  theApp.Run();
  //++++++++++++++++++++++++++++++++++++++++
}
