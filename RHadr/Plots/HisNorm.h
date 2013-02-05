/////////////////////////////////////////////////////////////////////
void HisNorm(TH1D *NorHst, TH1D *Hst){
  // normalize histogram in nanobarns
  Long_t   Nevt = NorHst->GetEntries();
  Double_t Xsav = NorHst->GetBinContent(0)/Nevt; // NANOBARNS
  cout<<"Xsav = "<<Xsav<<"  Nevt =  "<<Nevt<<endl;
  //
  int      nbt  = Hst->GetNbinsX();
  Double_t tmax = Hst->GetXaxis()->GetXmax();
  Double_t tmin = Hst->GetXaxis()->GetXmin();
  Double_t Fact = nbt*Xsav/(tmax-tmin)/Nevt;
  cout<<"Fact = "<<Fact<<endl;
  Hst->Scale(Fact);
}
/////////////////////////////////////////////////////////////////////
void HisNorm2(TH1D *NorHst, TH1D *Hst, Float_t msize, Int_t mcolor, Int_t mark){
  // normalize histogram in nanobarns and sets plotting params
  Long_t   Nevt = NorHst->GetEntries();
  Double_t Xsav = NorHst->GetBinContent(0)/Nevt/Nevt;
  //
  int      nbt  = Hst->GetNbinsX();
  Double_t tmax = Hst->GetXaxis()->GetXmax();
  Double_t tmin = Hst->GetXaxis()->GetXmin();
  Double_t Fact = nbt*Xsav/(tmax-tmin);
  //
  Hst->Scale(Fact);
  //
  Hst->SetMarkerStyle( mark);
  Hst->SetMarkerColor( mcolor);
  Hst->SetLineColor(   mcolor);
  Hst->SetMarkerSize(  msize);
}
