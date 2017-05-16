using namespace std;
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <TROOT.h>
#include <TClass.h>
#include <TFile.h>
#include <TH2.h>
#include <TKey.h>

/////////////////////////////////////////////////////////////////////
void HistNorm(TH1D *NorHst, TH1D *Hst){
  // normalize histogram in nanobarns
  Hst->ls();
  Double_t Nevt = NorHst->GetBinContent(2);
  Double_t Xsav = NorHst->GetBinContent(1)/Nevt; // NANOBARNS
  //
  int      nbt  = Hst->GetNbinsX();
  Double_t tmax = Hst->GetXaxis()->GetXmax();
  Double_t tmin = Hst->GetXaxis()->GetXmin();
  Double_t Fact = nbt*Xsav/(tmax-tmin)/Nevt;
  cout<<"HistNorm: Xsav = "<<Xsav<<"  Nevt =  "<<"  Fact = "<<Fact<<endl;
  Hst->Scale(Fact);
}
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
void HistNorm2(TH1D *NorHst, TH2D *Hst){
  // normalize histogram in nanobarns
  //Hst->ls();
  Double_t Nevt = NorHst->GetBinContent(2);
  Double_t Xsav = NorHst->GetBinContent(1)/Nevt; // NANOBARNS
  //
  int      nbtx = Hst->GetNbinsX();
  Double_t xmax = Hst->GetXaxis()->GetXmax();
  Double_t xmin = Hst->GetXaxis()->GetXmin();
  int      nbty = Hst->GetNbinsY();
  Double_t ymax = Hst->GetYaxis()->GetXmax();
  Double_t ymin = Hst->GetYaxis()->GetXmin();
  Double_t Fact = Xsav*nbtx*nbty/(xmax-xmin)/(ymax-ymin)/Nevt;
  cout<<"HistNorm2: Xsav = "<<Xsav<<"  Nevt =  "<<Nevt<<"  Fact = "<<Fact<<endl;
  Hst->Scale(Fact);
}
/////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
/// This useful routine is written by Phil Stephens
void HistNormalize(TString normname, TFile *file) {
  file->cd();  /// <-- Important!!!
  file->ls();
  TH1D  *hnorm = (TH1D*)file->Get(normname);
  TList *keys  = file->GetListOfKeys();
  TIterator *tit = keys->MakeIterator();
  TKey *key;
  while((key = (TKey*)(*tit)())) {
    if(!strcmp(hnorm->GetName(),  key->GetName())) continue;
    /// This method may not work for older root files
    TObject *obj = key->ReadObj();
    if(obj->IsA()->InheritsFrom("TH1D")) {
      TH1D *h = (TH1D*)obj;
      cout<<"///// Normalizing 1dim //// "<< h->GetName() <<" //// "<< h->GetTitle()<<endl;
      HistNorm(hnorm,h);
   } else if(obj->IsA()->InheritsFrom("TH2D")) {
      TH2D *h2 = (TH2D*)obj;
      cout<<"||||| Normalizing 2dim |||| "<< h2->GetName()<<" |||| "<< h2->GetTitle()<<endl;
      HistNorm2(hnorm,h2);
    }
    /// Another method which does not work any more (SJ)
    /*
    TObject *obj = file->Get(key->GetName());
    if(!obj) { cout << "Could not get object " << key->GetName() << endl; continue; }
    TH1D *h1 = (TH1D*)obj;
    if(h1) {
      cout<<"Normalizing 1dim: "<< h1->GetTitle()<<endl;
      HistNorm(hnorm,h1);
      continue; 
    }
    TH2D *h2 = (TH2D*)obj;
    if(h2) { 
      cout<<"Normalizing 2dim: "<< h2->GetTitle()<<endl;
      HistNorm2(hnorm,h2); 
      continue; 
    }
    *///
  }
}
