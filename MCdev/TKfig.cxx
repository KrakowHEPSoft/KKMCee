///////////////////////////////////////////////////////////////////////////////

#include"TKfig.h"

#define SP15 setw(15)<<setprecision(9)
#define SP10 setw(10)<<setprecision(7)
#define SW10 setw(10)

///////////////////////////////////////////////////////////////////////////////
TKfig::TKfig(){
  m_PlotSelector = "";
}/// 

/////////////////////////////////////////////////////////////////////
void TKfig::HistNorm(TH1D *NorHst, TH1D *Hst){
  // normalize histogram in nanobarns
  Hst->ls();
  Double_t Nevt = NorHst->GetBinContent(2);
  Double_t Xsav = NorHst->GetBinContent(1)/Nevt; // NANOBARNS
  //
  int      nbt  = Hst->GetNbinsX();
  Double_t tmax = Hst->GetXaxis()->GetXmax();
  Double_t tmin = Hst->GetXaxis()->GetXmin();
  Double_t Fact = nbt*Xsav/(tmax-tmin)/Nevt;
  cout<<"HistNorm: Xsav = "<<Xsav<<"  Nevt =  "<<Nevt<<"  Fact = "<<Fact<<endl;
  Hst->Scale(Fact);
}
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
void TKfig::HistNorm2(TH1D *NorHst, TH2D *Hst){
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
void TKfig::HistNormalize(TString normname, TFile *file) {
  file->cd();  /// <-- Important!!!
  file->ls();
  TH1D  *hnorm = (TH1D*)file->Get(normname);
  cout << "--------HistNormalize:  norma = " << hnorm->GetName()<< endl;
  TList *keys  = file->GetListOfKeys();
  TIterator *tit = keys->MakeIterator();
  TKey *key;
  while((key = (TKey*)(*tit)())) {
    if(!strcmp(hnorm->GetName(),  key->GetName())) continue;
    /// This method may not work for older root files
    TObject *obj = key->ReadObj();
    if(obj->IsA()->InheritsFrom("TH1D")) {
      TH1D *h = (TH1D*)obj;
//if (!h) cout << "--------HistNormalize:  no histo!!!";
      cout<<"///// Normalizing 1dim //// "<< h->GetName() <<" //// "<< h->GetTitle()<<endl;
      HistNorm(hnorm,h);
   } else if(obj->IsA()->InheritsFrom("TH2D")) {
      TH2D *h2 = (TH2D*)obj;
      cout<<"||||| Normalizing 2dim |||| "<< h2->GetName()<<" |||| "<< h2->GetTitle()<<endl;
      HistNorm2(hnorm,h2);
    }

  }
cout << "--------HistNormalize:  the end";
}


//______________________________________________________________________________
TH1D *TKfig::DefHist1D(TString ScaleType, TString FuncType, TH1D *Hst, const char* name, const char* title)
{ //! Define histo and fill in with the analytical function
  //! Binning taken from TH1D *Hst
  //! Mainly for 4-dim MC exercises
 
  TH1D *histo = (TH1D*)Hst->Clone(name);
  histo->SetTitle(title);
  int      nbx  = Hst->GetNbinsX();
  Double_t xmax = Hst->GetXaxis()->GetXmax();
  Double_t xmin = Hst->GetXaxis()->GetXmin();
  Double_t dx   = (xmax-xmin)/nbx;
  //cout<<"nbx,xmin,xmax="<<nbx<<"  "<<xmin<<"  "<<xmax<<endl;
  Double_t x,xl,fun;
  if( ScaleType == "lin_x"){
  for(int i=1; i<=nbx;i++){
    //fun = Dist1Selector6(FuncType, xmin +dx*(i-0.5));
    fun = 0.5*(
          0.88888889*DistSelector1D(xmin +dx*(i-0.5),FuncType)
        + 0.55555555*DistSelector1D( xmin +dx*(i-0.5 +0.5*0.77459667), FuncType)
        + 0.55555555*DistSelector1D( xmin +dx*(i-0.5 -0.5*0.77459667),FuncType));
    histo->SetBinContent(i,fun);
    histo->SetBinError(i,0.0);
    //cout<<"x, fun ="<<x<<"  "<<fun<<endl;
    }
  } else if( ScaleType == "log_x"){
  for(int i=1; i<=nbx;i++){
    //fun = DistSelector1D(FuncType, exp(log(10)* (xmin +dx*(i-0.5)) ));
    fun = 0.5*(
          0.88888889*DistSelector1D( exp(log(10)* (xmin +dx*(i-0.5))),FuncType)
        + 0.55555555*DistSelector1D( exp(log(10)* (xmin +dx*(i-0.5 +0.5*0.77459667))),FuncType)
        + 0.55555555*DistSelector1D( exp(log(10)* (xmin +dx*(i-0.5 -0.5*0.77459667))),FuncType));
    histo->SetBinContent(i,fun);
    histo->SetBinError(i,0.0);
    }
  }

  return histo;
}


//______________________________________________________________________________________
TH2D*  TKfig::DefHist2D(TString ScaleType, TString FuncType, TH2D *histo, const char* name, const char* title)
{
  TH2D *Hst = (TH2D*)histo->Clone(name);
  Hst->SetTitle(title);
/// range taken from *Hst

  cout << "Function type: " << FuncType << endl; 
  
  int      Nbinx  = Hst->GetNbinsX();
  int      Nbiny  = Hst->GetNbinsY();
  Double_t lxmax = Hst->GetXaxis()->GetXmax();
  Double_t lxmin = Hst->GetXaxis()->GetXmin();
  Double_t lymax = Hst->GetYaxis()->GetXmax();
  Double_t lymin = Hst->GetYaxis()->GetXmin();
  
  //Double_t dx   = (lxmax-lxmin)/nbx;
  Double_t x,xl,fun;
if( ScaleType == "lin_x"){
  for (int ix= 1; ix <= Nbinx; ix++)
  {
    for (int iy= 1; iy <= Nbiny; iy++)
    {
      double x = (lxmax-lxmin)*(ix-0.5)/Nbinx + lxmin; //logarithmic coordinates on the plot
      double y = (lymax-lymin)*(iy-0.5)/Nbiny + lymin;
      double alfR = pow(10., y);
      double yy   = pow(10., x);
      double ykT  = yy; 
      double fun=DistSelector2D(x, y, FuncType);

      Hst->SetBinContent(Hst->GetBin(ix, iy), fun);
      Hst->SetBinError(Hst->GetBin(ix, iy),0.0);
    }///for iy
  }///for ix
  cout<< " **************** DefHist2D ends *******************   "<<endl;
  
 }//scaletype
return Hst;
}

double TKfig::DistSelector1D( double x, TString type=" " )
{
  if (type == "") 		return 0;
  else  if(  type=="ONE") 	return 1;
  else{
    cout << endl << "TKfig::DistSelector1D error: wrong function type" << endl;
 				return -1;
  }
}

double TKfig::DistSelector2D( double x, double y, TString type=" ")
{
  if (type == "") 		return 0;
  else  if(  type=="ONE") 	return 1;
  else{
    cout << endl << "TKfig::DistSelector1D error: wrong function type" << endl;
 				return -1;
}
}


//=========================================================================== 
Double_t TKfig::Dilogy(double x)
{
//==========================================================================//
// Translated from FORTRAN to C++ by S.Jadach.                              //
//--------------------------------------------------------------------------//
// Cross-checked against the FORTRAN original version.                      //
// A very good numerical agreement found: relative difference <1e-16.       //
// Wieslaw Placzek,                                         CERN,07.09.2004 //
//==========================================================================//
// DOUBLE PRECISION FUNCTION Mathlib_dilogy(x) 
//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// dilogarithm FUNCTION: dilog(x)=int( -ln(1-z)/z ) , 0 < z < x             //
// this is the cernlib version.                                             //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
   double a, b, s, t, y, z, lgabx;
//------------------------------------------------------------------------------
   z=-1.644934066848226e0;
   if (x < -1.e0)
       goto e1;
   if (x < 0.5e0)
       goto e2;
   if (x == 1.e0)
       goto e3;
   if (x <= 2.e0)
       goto e4;
   z=3.289868133696453e0;
   e1: t=1.e0/x;
   s=-0.5e0;
   lgabx = log(fabs(x));
   z=z-0.5e0*lgabx*lgabx;
   goto e5;
   e2: t=x;
   s=0.5e0;
   z=0.e0;
   goto e5;
   e3: return 1.644934066848226e0;
   e4: t=1.e0-x;
   s=-0.5e0;
   z=1.644934066848226e0-log(x)*log(fabs(t));
   e5: y=2.666666666666667e0*t+0.666666666666667e0;
   b= 0.000000000000001e0;
   a=y*b +0.000000000000004e0;
   b=y*a-b+0.000000000000011e0;
   a=y*b-a+0.000000000000037e0;
   b=y*a-b+0.000000000000121e0;
   a=y*b-a+0.000000000000398e0;
   b=y*a-b+0.000000000001312e0;
   a=y*b-a+0.000000000004342e0;
   b=y*a-b+0.000000000014437e0;
   a=y*b-a+0.000000000048274e0;
   b=y*a-b+0.000000000162421e0;
   a=y*b-a+0.000000000550291e0;
   b=y*a-b+0.000000001879117e0;
   a=y*b-a+0.000000006474338e0;
   b=y*a-b+0.000000022536705e0;
   a=y*b-a+0.000000079387055e0;
   b=y*a-b+0.000000283575385e0;
   a=y*b-a+0.000001029904264e0;
   b=y*a-b+0.000003816329463e0;
   a=y*b-a+0.000014496300557e0;
   b=y*a-b+0.000056817822718e0;
   a=y*b-a+0.000232002196094e0;
   b=y*a-b+0.001001627496164e0;
   a=y*b-a+0.004686361959447e0;
   b=y*a-b+0.024879322924228e0;
   a=y*b-a+0.166073032927855e0;
   a=y*a-b+1.935064300869969e0;
   return s*t*(a-b)+z;
}
//===========================================================================


