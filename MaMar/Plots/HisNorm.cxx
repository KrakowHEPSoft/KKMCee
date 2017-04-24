/////////////////////////////////////////////////////////////////////
// Collection of programs for renormalizing histograms
// and projecting 2D scatergrams into 1D histograms
// Also tools for calculating AFB(vmax)
/////////////////////////////////////////////////////////////////////
#include "HisNorm.h"

double sqr( const Double_t x ){ return x*x;};


// This works for 1-dim histograms
void HisNorm0( long   Nevt, double Xsav, TH1 *Hst){
// normalize histogram according to Xsav
  cout<<"HisNorm0: Xsav = "<<Xsav<<"  Nevt =  "<<Nevt<<endl;
  //
  int      nbX  = Hst->GetNbinsX();
  //cout<<"nbt = "<<nbt<<endl;
  Double_t Xmax = Hst->GetXaxis()->GetXmax();
  Double_t Xmin = Hst->GetXaxis()->GetXmin();
  Double_t Fact = nbX*Xsav/(Xmax-Xmin)/Nevt;
  cout<<"HisNorm0: Fact = "<<Fact<<endl;
  Hst->Scale(Fact);
}

/////////////////////////////////////////////////////////////////////
// This works for 1-dim histograms
void HisNorm1(TH1D *NorHst, TH1 *Hst){
  // normalize histogram in nanobarns
  Long_t   Nevt = NorHst->GetEntries();
  Double_t Xsav = NorHst->GetBinContent(0)/Nevt; // NANOBARNS
  cout<<"HisNorm1: Xsav = "<<Xsav<<"  Nevt =  "<<Nevt<<endl;
  //
  int      nbX  = Hst->GetNbinsX();
  //cout<<"nbt = "<<nbt<<endl;
  Double_t Xmax = Hst->GetXaxis()->GetXmax();
  Double_t Xmin = Hst->GetXaxis()->GetXmin();
  Double_t Fact = nbX*Xsav/(Xmax-Xmin)/Nevt;
  cout<<"HisNorm1: Fact = "<<Fact<<endl;
  Hst->Scale(Fact);
}


/////////////////////////////////////////////////////////////////////
void HisNorm1M(TH1D *NorHst, TH1D *Hst, Float_t msize, Int_t mcolor, Int_t mark){
  // normalize histogram in nanobarns and sets plotting params
  Long_t   Nevt = NorHst->GetEntries();
  Double_t Xsav = NorHst->GetBinContent(0)/Nevt/Nevt;
  //
  int      nbX  = Hst->GetNbinsX();
  Double_t Xmax = Hst->GetXaxis()->GetXmax();
  Double_t Xmin = Hst->GetXaxis()->GetXmin();
  Double_t Fact = nbX*Xsav/(Xmax-Xmin);
  //
  Hst->Scale(Fact);
  //
  Hst->SetMarkerStyle( mark);
  Hst->SetMarkerColor( mcolor);
  Hst->SetLineColor(   mcolor);
  Hst->SetMarkerSize(  msize);
}

/////////////////////////////////////////////////////////////////////
// This works for 2-dim histograms
void HisNorm2(TH1D *NorHst, TH2 *Hst){
  // normalize histogram in nanobarns
  Long_t   Nevt = NorHst->GetEntries();
  Double_t Xsav = NorHst->GetBinContent(0)/Nevt; // NANOBARNS
  cout<<"HisNorm2: Xsav = "<<Xsav<<"  Nevt =  "<<Nevt<<endl;
  //
  int      nbX  = Hst->GetNbinsX();
  int      nbY  = Hst->GetNbinsY();
  //cout<<"nbt = "<<nbt<<endl;
  Double_t Xmax = Hst->GetXaxis()->GetXmax();
  Double_t Xmin = Hst->GetXaxis()->GetXmin();
  Double_t Ymax = Hst->GetYaxis()->GetXmax();
  Double_t Ymin = Hst->GetYaxis()->GetXmin();
  Double_t Fact = Xsav/Nevt;
  Fact *= nbX/(Xmax-Xmin);
  Fact *= nbY/(Ymax-Ymin);
  cout<<"HisNorm2: Fact = "<<Fact<<endl;
  Hst->Scale(Fact);
}

///////////////////////////////////////////////////////////////////////////////////
void ProjX1(TH2D *Scat, TH1D *&HstProjX)
{
  // Simple Projection onto X axis taking into account errors!
  int      nbX  = Scat->GetNbinsX();
  int      nbY  = Scat->GetNbinsY();
  Double_t Xmax = Scat->GetXaxis()->GetXmax();
  Double_t Xmin = Scat->GetXaxis()->GetXmin();
  Double_t Ymax = Scat->GetYaxis()->GetXmax();
  Double_t Ymin = Scat->GetYaxis()->GetXmin();
  //
  HstProjX = (TH1D*)Scat->ProjectionX("HstNew",1,nbX,"e"); // option "e" doesnt work
  HstProjX->Reset();
  double sum,sum2, dy;
  dy= (Ymax-Ymin)/nbY; // integration over Y
  for(int ix=0; ix <= nbX+1; ix++){
    sum=0.0; sum2=0.0;
    for(int iy=0; iy <= nbY+1; iy++){
      sum  += Scat->GetBinContent(ix,iy);
      sum2 += sqr(Scat->GetBinError(ix,iy));
    }
    HstProjX->SetBinContent(ix, dy*sum);
    HstProjX->SetBinError(  ix, dy*sqrt(sum2));
  }
}// ProjX1

///////////////////////////////////////////////////////////////////////////////////
void ProjY1(TH2D *Scat, TH1D *&HstProjY)
{
  // Simple Projection onto Y axis taking into account errors!
  int      nbX  = Scat->GetNbinsX();
  int      nbY  = Scat->GetNbinsY();
  Double_t Xmax = Scat->GetXaxis()->GetXmax();
  Double_t Xmin = Scat->GetXaxis()->GetXmin();
  Double_t Ymax = Scat->GetYaxis()->GetXmax();
  Double_t Ymin = Scat->GetYaxis()->GetXmin();
  //
  HstProjY = (TH1D*)Scat->ProjectionY("HstNew",1,nbX,"e"); // option "e" doesnt work
  HstProjY->Reset();
  double sum,sum2, dx;
  dx= (Xmax-Xmin)/nbX; // integration over X
  for(int iy=0; iy <= nbY+1; iy++){
    sum=0.0; sum2=0.0;
    for(int ix=0; ix <= nbX+1; ix++){
      sum  += Scat->GetBinContent(ix,iy);
      sum2 += sqr(Scat->GetBinError(ix,iy));
    }
    HstProjY->SetBinContent(iy, dx*sum);
    HstProjY->SetBinError(  iy, dx*sqrt(sum2));
  }
}// ProjY1

///////////////////////////////////////////////////////////////////////////////////
void ProjV(TH2D *Scat, TH1D *&hxTot, TH1D *&hxAfb, int NbMax)
{
  //  Projection onto v axis, suming over Y=cos(theta) up to a limit
  int      nbX  = Scat->GetNbinsX();
  int      nbY  = Scat->GetNbinsY();
  Double_t Xmax = Scat->GetXaxis()->GetXmax();
  Double_t Xmin = Scat->GetXaxis()->GetXmin();
  Double_t Ymax = Scat->GetYaxis()->GetXmax();
  Double_t Ymin = Scat->GetYaxis()->GetXmin();
  //
  hxTot = (TH1D*)Scat->ProjectionX("HstNew",1,nbX,"e");
  hxTot->Reset();
  hxAfb = (TH1D*)hxTot->Clone("HstNew2");
  //
  double forw,forw2, back, back2, dx,dy;
  double Forw,Forw2, Back, Back2;
  dx= (Xmax-Xmin)/nbX; // integration over X
  dy= (Ymax-Ymin)/nbY; // integration over Y
  Forw=0.0; Forw2=0.0;
  Back=0.0; Back2=0.0;
  int nbYhalf = nbY/2;
  int nbY2;
  if( (NbMax>0) && (NbMax<nbYhalf) )
    nbY2 = NbMax;
  else
    nbY2 = nbYhalf;
  for(int ix=0; ix <= nbX+1; ix++){
    forw=0.0; forw2=0.0;
    back=0.0; back2=0.0;
    // loop over cos(theta) bins
    for(int iy=1; iy <= nbY2; iy++){
      forw  += Scat->GetBinContent(  ix, nbYhalf+iy);
      forw2 += sqr(Scat->GetBinError(ix, nbYhalf+iy));
      back  += Scat->GetBinContent(  ix, nbYhalf-iy+1);
      back2 += sqr(Scat->GetBinError(ix, nbYhalf-iy+1));
    }// iy
    Forw  += forw;  Forw2 += forw2;
    Back  += back;  Back2 += back2;
    //hxTot->SetBinContent(ix, dy*sum);
    //hxTot->SetBinError(  ix, dy*sqrt(sum2));
    hxTot->SetBinContent(ix, dx*dy*    (Forw +Back));
    hxTot->SetBinError(  ix, dx*dy*sqrt(Forw2+Back2));
    hxAfb->SetBinContent(ix, dx*dy*    (Forw -Back));
    hxAfb->SetBinError(  ix, dx*dy*sqrt(Forw2+Back2));
  }
  hxAfb->Divide(hxTot);
}// ProjV

///////////////////////////////////////////////////////////////////////////////////
void ProjC(TH2D *Scat, TH1D *&hTot, TH1D *&hAsy, int NbMax)
{
  // Projection onto c=cos(theta) axis, suming over v up to a limit
  int      nbX  = Scat->GetNbinsX();
  int      nbY  = Scat->GetNbinsY();
  Double_t Xmax = Scat->GetXaxis()->GetXmax();
  Double_t Xmin = Scat->GetXaxis()->GetXmin();
  Double_t Ymax = Scat->GetYaxis()->GetXmax();
  Double_t Ymin = Scat->GetYaxis()->GetXmin();
  //
  hTot = (TH1D*)Scat->ProjectionY("HstNew",1,nbX,"e"); // option "e" doesnt work
  hTot->Reset();
  hAsy = (TH1D*)hTot->Clone("HstNew1");
  TH1D *hSym = (TH1D*)hTot->Clone("HstNew2"); // local histo
  //
  double forw,forw2;
  double back,back2;
  double dx= (Xmax-Xmin)/nbX; // integration over X
  int ivMax;
  //
  if( (NbMax>0) && (NbMax<nbX) )
    ivMax = NbMax;
  else
    ivMax = nbX;
  int nbYhalf = nbY/2;
  for(int ic=1; ic <= nbYhalf; ic++){
    int iforw = nbYhalf+ic;
    int iback = nbYhalf-ic+1;
    forw=0.0;  forw2=0.0;
    back=0.0;  back2=0.0;
    for(int iv=0; iv <= ivMax; iv++){
      forw  += Scat->GetBinContent(  iv,iforw);
      forw2 += sqr(Scat->GetBinError(iv,iforw));
      back  += Scat->GetBinContent(  iv,iback);
      back2 += sqr(Scat->GetBinError(iv,iback));
    }
    // total
    hTot->SetBinContent(iforw,     dx*(forw));
    hTot->SetBinError(  iforw, dx*sqrt(forw2));
    hTot->SetBinContent(iback,     dx*(back));
    hTot->SetBinError(  iback, dx*sqrt(back2));
    // asymetric
    hAsy->SetBinContent(iforw,     dx*(forw-back));
    hAsy->SetBinError(  iforw, dx*sqrt(forw2+back2));
    hAsy->SetBinContent(iback,    -dx*(forw-back));
    hAsy->SetBinError(  iback, dx*sqrt(forw2+back2));
    // symetric
    hSym->SetBinContent(iforw,     dx*(forw+back));
    hSym->SetBinError(  iforw, dx*sqrt(forw2+back2));
    hSym->SetBinContent(iback,     dx*(forw+back));
    hSym->SetBinError(  iback, dx*sqrt(forw2+back2));
  }
  hAsy->Divide(hSym);
  hSym->Delete();
}// ProjC

///////////////////////////////////////////////////////////////////////////////////
void MakeCumul(TH1D *hst1, TH1D *&hcum1)
{
  // makes cumulative distribution
  cout<<"Entering MakeCumul for  ";
  cout<< hst1->GetName() <<endl;
  int      nbX  = hst1->GetNbinsX();
  Double_t Xmax = hst1->GetXaxis()->GetXmax();
  Double_t Xmin = hst1->GetXaxis()->GetXmin();
  //
  hcum1 = (TH1D*)hst1->Clone("hcum1"); // allocate hcum1
  hcum1->Reset();
  double sum=0 ,sum2=0;
  double dx= (Xmax-Xmin)/nbX; // integration over X
  for(int iv=0; iv <= nbX; iv++){
    sum   += hst1->GetBinContent(  iv);
    sum2 += sqr(hst1->GetBinError(iv));
    hcum1->SetBinContent(iv,     dx*(sum));
    hcum1->SetBinError(  iv, dx*sqrt(sum2));
   }
//

}//MakeCumul

  ///////////////////////////////////////////////////////////////////////////////////
void MakeAFB(TH1D *hAll, TH1D *&hAFB)
{
  // makes Afb(c) out of hAll(c)
  int      nbX  = hAll->GetNbinsX();
  Double_t Xmax = hAll->GetXaxis()->GetXmax();
  Double_t Xmin = hAll->GetXaxis()->GetXmin();
  //
  hAFB = (TH1D*)hAll->Clone("HstNew3"); // allocate hAFB
  hAFB->Reset();
  TH1D *hSym = (TH1D*)hAFB->Clone("HstNew4"); // local temporary histo
  hAFB->Reset();
  //
  double forw,forw2;
  double back,back2;
  double dx= (Xmax-Xmin)/nbX; // integration over X
  int nbXhalf = nbX/2;
  for(int ic=1; ic <= nbXhalf; ic++){
    int iforw = nbXhalf+ic;
    int iback = nbXhalf-ic+1;
    //cout<<"************** "<<iforw<<" "<<iback<<endl;
    forw  = hAll->GetBinContent(  iforw);
    forw2 = sqr(hAll->GetBinError(iforw));
    back  = hAll->GetBinContent(  iback);
    back2 = sqr(hAll->GetBinError(iback));
    // asymetric
    hAFB->SetBinContent(iforw,     dx*(forw-back));
    hAFB->SetBinError(  iforw, dx*sqrt(forw2+back2));
    hAFB->SetBinContent(iback,    -dx*(forw-back));
    hAFB->SetBinError(  iback, dx*sqrt(forw2+back2));
    // symetric
    hSym->SetBinContent(iforw,     dx*(forw+back));
    hSym->SetBinError(  iforw, dx*sqrt(forw2+back2));
    hSym->SetBinContent(iback,     dx*(forw+back));
    hSym->SetBinError(  iback, dx*sqrt(forw2+back2));
    cout<<"MakeAFB: (f-b)/2,f,b==="<<ic<<" "<< (forw-back)/2<<" b="<<back<<" f="<<forw<<endl;
  }
  hAFB->Divide(hSym);
  hSym->Delete();
  //for(int ic=1; ic <= nbX; ic++)
  //  cout<<"***"<<ic<<" "<< hAFB->GetBinContent(ic)<<endl;
}// ProjC
