///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                Class  TFMAXWT                                             //
//                                                                           //
//      Small auxiliary class for controlling MC weight.                     //
//      It provides certain measure of the "maximum weight"                  //
//      depending on small user-parameter "epsilon"                          //
//      It creates and uses 2 histograms of the TH1D class.                  //
//      User defines no. of bins nBin,  nBin=1000 is  recommended            //
//      wmax defines weight range (1,wmax), it is adjusted "manually"        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include<math.h>
#include<stdlib.h>
#include"TFMAXWT.h"

#ifdef ROOT_DEF
#include"TH1.h"
#else
#include "TFHST.h"
#endif

#ifdef ROOT_DEF
ClassImp(TFMAXWT)
#endif

///////////////////////////////////////////////////////////////////////////////
TFMAXWT::TFMAXWT(){
// Constructor for streamer
  m_Nent = 0;
  m_nBin = 0;
  m_WtHst1 = NULL;
  m_WtHst2 = NULL;
}

///////////////////////////////////////////////////////////////////////////////
TFMAXWT::TFMAXWT(const double wmax, const int nBin){
// practical principal constructor
  m_Nent = 0;
  m_nBin = nBin;
  m_wmax = wmax;
#ifdef ROOT_DEF
  m_WtHst1 = new TH1D("TFMAXWT_hst_Wt1","Histo of weight   ",nBin,0.0,wmax);
  m_WtHst2 = new TH1D("TFMAXWT_hst_Wt2","Histo of weight**2",nBin,0.0,wmax);
  m_WtHst1->SetDirectory(0);// exclude from diskfile
  m_WtHst2->SetDirectory(0);// and enable deleting
#else
  m_WtHst1 = new TFHST(0.0, wmax, m_nBin);      // histogram of wt
  m_WtHst2 = new TFHST(0.0, wmax, m_nBin);      // histogram of wt with wt
#endif
}

///////////////////////////////////////////////////////////////////////////////
TFMAXWT::TFMAXWT(TFMAXWT &From){
// Explicit COPY CONSTRUCTOR (unused, so far)
  m_nBin   = From.m_nBin;
  m_wmax   = From.m_wmax;
  m_WtHst1 = From.m_WtHst1;
  m_WtHst2 = From.m_WtHst2;
  cout<<"+++++ Stop in  COPY CONSTRUCTOR TFMAXWT::TFMAXWT,   NOT TESTED!"<<endl;
  exit(1);  
}

///////////////////////////////////////////////////////////////////////////////
TFMAXWT::~TFMAXWT(){
  delete m_WtHst1; // For this SetDirectory(0) is needed!
  delete m_WtHst2; //
  m_WtHst1=NULL;
  m_WtHst2=NULL;
}
///////////////////////////////////////////////////////////////////////////////
void TFMAXWT::Reset(){
  m_Nent = 0;
  m_WtHst1->Reset();
  m_WtHst2->Reset();
}

///////////////////////////////////////////////////////////////////////////////
TFMAXWT& TFMAXWT::operator =(TFMAXWT &From){
// substitution = 
  if (&From == this) return *this;
  m_nBin = From.m_nBin;
  m_wmax = From.m_wmax;
  m_WtHst1 = From.m_WtHst1;
  m_WtHst2 = From.m_WtHst2;
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
void TFMAXWT::Fill(double wt){
  m_Nent =  m_Nent+1.0; 
  m_WtHst1->Fill(wt,1.0);
  m_WtHst2->Fill(wt,wt);
}

///////////////////////////////////////////////////////////////////////////////
void TFMAXWT::Make(const double eps, double &MCeff){
// Calculates Efficiency= AveWt/WtLim for a given tolerance level epsilon<<1
  double WtLim,AveWt;
  GetMCeff(eps, MCeff, WtLim);
  AveWt = MCeff*WtLim;
  cout<< "00000000000000000000000000000000000000000000000000000000000000000000000"<<endl;
  cout<< "00 -->WtLim: No_evt ="<<m_Nent<<"   <Wt> = "<<AveWt<<"  WtLim=  "<<WtLim<<endl;
  cout<< "00 -->WtLim: For eps = "<<eps  <<"    EFFICIENCY <Wt>/WtLim= "<<MCeff<<endl;
  cout<< "00000000000000000000000000000000000000000000000000000000000000000000000"<<endl;
}

///////////////////////////////////////////////////////////////////////////////
void TFMAXWT::GetMCeff(const double eps, double &MCeff, double &WtLim){
// Calculates Efficiency= AveWt/WtLim for a given tolerance level epsilon<<1
// using information stored in two histograms
  int ib,ibX;
  double LowEdge,Bin,Bin1;
  double AveWt, AveWt1;

  m_WtHst1->Print();
  m_WtHst2->Print();

  // Convention on bin-numbering: nb=1 for 1-st bin, undeflow nb=0, overflow nb=Nb+1
  double sum   = 0.0;
  double sumWt = 0.0;
  for(ib=0;ib<=m_nBin+1;ib++){
    sum   += m_WtHst1->GetBinContent(ib);
    sumWt += m_WtHst2->GetBinContent(ib);
  }
  if( (sum == 0.0) || (sumWt == 0.0) ){
    cout<<"TFMAXWT::Make: zero content of histogram !!!,sum,sumWt ="<<sum<<sumWt<<endl;
  }
  AveWt = sumWt/sum;
  //--------------------------------------
  for( ibX=m_nBin+1; ibX>0; ibX--){
    LowEdge = (ibX-1.0)*m_wmax/m_nBin;
    sum   = 0.0;
    sumWt = 0.0;
    for( ib=0; ib<=m_nBin+1; ib++){
      Bin  = m_WtHst1->GetBinContent(ib);
      Bin1 = m_WtHst2->GetBinContent(ib);
      if(ib >= ibX) Bin1=LowEdge*Bin;
      sum   += Bin;
      sumWt += Bin1;
    }
    AveWt1 = sumWt/sum;
    if( fabs(1.0-AveWt1/AveWt) > eps ) break;
  }
  //---------------------------
  if(ibX == (m_nBin+1) ){
    WtLim = 1.0e200;
    MCeff   = 0.0;
    cout<< "+++++ WtLim undefined. Higher uper limit in histogram"<<endl;
  }else if( ibX == 1){
    WtLim = 0.0;
    MCeff   =-1.0;
    cout<< "+++++ WtLim undefined. Lower uper limit or more bins "<<endl;
  }else{
    WtLim= (ibX)*m_wmax/m_nBin; // We over-estimate WtLim, under-estimate MCeff
    MCeff  = AveWt/WtLim;
  }
}
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//      End of    Class  TFMAXWT                                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
