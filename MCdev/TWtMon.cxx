///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                Class  WtMonit                                             //
//                                                                           //
//      Small auxiliary class for detailed monitoring single MC weight.      //
//      It creates and uses 1 histograms.                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TWtMon.h"

#define SP15 setw(15)<<setprecision(9)
#define SP10 setw(10)<<setprecision(7)
#define SW10 setw(10)

ClassImp(TWtMon);



///---------------------------------------------------------------------------
TWtMon::TWtMon(){
/// explicit default constructor for ROOT streamers
/// should not be used by the user
  cout<< "||||> TWtMon::TWtMon: DEFAULT Constructor (for ROOT only) "<<endl;
  m_HST = NULL;
}///TWtMon


///////////////////////////////////////////////////////////////////////////////
TWtMon::TWtMon(const char* name)
{
  // principal constructor
  cout<< "||||> TWtMon::TWtMon: User Constructor "<<endl;
  m_HST = new TH1D(name, " TWtMon histogram", 10, 0.0, 1.0);
  m_Xup =1.0;
  m_Ntot = 0.0;
  m_Nacc = 0.0;
  m_Nzer = 0.0;
  m_Wtmax= -1e150;
  m_Wtmin= +1e150;  
  m_Wtsum = 0.0;
  m_Wtsu2 = 0.0;
  m_WtsuOv= 0.0;
}
///////////////////////////////////////////////////////////////////////////////
TWtMon::TWtMon(const char* name, const char* title, int nbinsx, const double xup){
  // principal constructor
  cout<< "||||> TWtMon::TWtMon: User principal Constructor "<<endl;
  m_HST = new TH1D(name, title, nbinsx, 0.0, xup);
  m_Xup = xup;
  m_Ntot = 0.0;
  m_Nacc = 0.0;
  m_Nzer = 0.0;
  m_Wtmax= -1e150;
  m_Wtmin= +1e150;  
  m_Wtsum = 0.0;
  m_Wtsu2 = 0.0;
  m_WtsuOv= 0.0;
}
///////////////////////////////////////////////////////////////////////////////
void TWtMon::Reset(){
// reinitialize and reset
  m_HST->Reset();
  m_Ntot = 0.0;
  m_Nacc = 0.0;
  m_Nzer = 0.0;
  m_Wtmax= -1e150;
  m_Wtmin= +1e150;
  m_Wtsum = 0.0;
  m_Wtsu2 = 0.0;
  m_WtsuOv= 0.0;
}
///////////////////////////////////////////////////////////////////////////////
void TWtMon::Fill(const double wt){
// fill in single MC weight
  m_Ntot   += 1.0;
  m_Wtsum  += wt;
  m_Wtsu2  += wt*wt;
  m_Nacc   += 1.0;
  m_WtsuOv += wt-m_Xup;
  if( wt == 0.0 )   m_Nzer += 1.0;
  if(m_Wtmax < wt)  m_Wtmax =wt;
  if(m_Wtmin > wt)  m_Wtmin =wt;
  m_HST->Fill(wt);
}
///////////////////////////////////////////////////////////////////////////////
void TWtMon::Fill(const double wt, const double rand){
// fill in single MC weight, optional arg. rand is rnd.numb. used for rejection
  m_Ntot   += 1.0;
  m_Wtsum  += wt;
  m_Wtsu2  += wt*wt;
  if(wt>m_Xup)
    m_WtsuOv += wt-m_Xup;
  if( wt <= rand )  m_Nacc += 1.0;
  if( wt == 0.0 )   m_Nzer += 1.0;
  if(m_Wtmax < wt)  m_Wtmax =wt;
  if(m_Wtmin > wt)  m_Wtmin =wt;
  m_HST->Fill(wt);
}
///////////////////////////////////////////////////////////////////////////////
void TWtMon::GetAver(double& AveWt, double& ErrAbs){
  // Calculate arerage weight and its error
  AveWt=0.0;
  ErrAbs=0.0;
  double sigma =0.0;
  if(m_Ntot != 0.0) {
    AveWt = m_Wtsum/m_Ntot;
    sigma = sqrt(m_Wtsu2/m_Ntot -AveWt*AveWt);
    ErrAbs = sigma/sqrt(m_Ntot);
  }
}
///////////////////////////////////////////////////////////////////////////////
void TWtMon::GetAll(
      double& AveWt, double& ERela, double& WtMax, double& WtMin, double& AvUnd, double& AvOve,
      double& Ntot,  double& Nacc,  double& Nneg,  double& Nove,  double& Nzer ){
  // Define all parameters of the weight
  // We dont use GetMean() and GetRMS(), they exclude under/overflow!
  AveWt=0.0;
  ERela=0.0;
  double sigma =0.0;
  if(m_Ntot != 0.0) {
    AveWt = m_Wtsum/m_Ntot;
    sigma = sqrt(m_Wtsu2/m_Ntot -AveWt*AveWt);
    ERela = sigma/sqrt(m_Ntot)/AveWt;
  }
  long Nbin = m_HST->GetNbinsX();
  AvUnd=0.0;
  AvOve=0.0;
  if(m_Ntot != 0.0) {
    AvUnd= m_HST->GetBinContent(0)/m_Ntot;
    AvOve=m_WtsuOv/m_Ntot;
  }
  WtMax = m_Wtmax;
  WtMin = m_Wtmin;
  Ntot  = m_Ntot;
  Nacc  = m_Nacc;
  Nzer  = m_Nzer;
  Nneg  = m_HST->GetBinContent(0);
  Nove  = m_HST->GetBinContent(Nbin+1);
}
///////////////////////////////////////////////////////////////////////////////
void TWtMon::PrintAll(){
  // print all parameters of the weight
  double AveWt,ERela,WtMax,WtMin,AvUnd,AvOve;
  double Ntot, Nacc, Nneg, Nove, Nzer;
  GetAll( AveWt,  ERela,  WtMax,  WtMin,  AvUnd,  AvOve,
          Ntot,   Nacc,   Nneg,   Nove,   Nzer );
  const char *title= m_HST->GetTitle();
  cout<<"///////////////////////////// "<<title;
  cout<<" /////////////////////////////"<<endl;
  cout<<"  Ntot  = "<<SP15<< Ntot<<endl;
  cout<<"  Nacc  = "<<SP15<< Nacc<<endl;
  cout<<"  Nneg  = "<<SP15<< Nneg<<endl;
  cout<<"  Nove  = "<<SP15<< Nove<<endl;
  cout<<"  Nzer  = "<<SP15<< Nzer<<endl;
  cout<<"  AveWt = "<<SP15<<AveWt<<endl;
  cout<<"  ERela = "<<SP15<<ERela<<endl;
  cout<<"  WtMax = "<<SP15<<WtMax<<endl;
  cout<<"  WtMin = "<<SP15<<WtMin<<endl;
  cout<<"  AvOve = "<<SP15<<AvOve<<endl;
  cout<<"  AvUnd = "<<SP15<<AvUnd<<endl;
}
///////////////////////////////////////////////////////////////////////////////
void TWtMon::PrintLine(){
  // print all parameters of the weight
  double AveWt,ERela,WtMax,WtMin,AvUnd,AvOve;
  double Ntot, Nacc, Nneg, Nove, Nzer;
  GetAll( AveWt,  ERela,  WtMax,  WtMin,  AvUnd,  AvOve,
          Ntot,   Nacc,   Nneg,   Nove,   Nzer );
  const char *title= m_HST->GetTitle();
  cout<<"///////////////////////////// "<<title;
  cout<<" /////////////////////////////"<<endl;
  cout<<setw(15)<< "Ntot";
  cout<<setw(15)<< "Nacc";
  cout<<setw(15)<< "Nneg";
  cout<<setw(15)<< "Nove";
  cout<<setw(15)<< "Nzer"<<endl;
  cout<<setw(15)<< Ntot;
  cout<<setw(15)<< Nacc;
  cout<<setw(15)<< Nneg;
  cout<<setw(15)<< Nove;
  cout<<setw(15)<< Nzer<<endl;
  cout<<setw(15)<< "AveWt";
  cout<<setw(15)<< "ERela";
  cout<<setw(15)<< "WtMax";
  cout<<setw(15)<< "WtMin";
  cout<<setw(15)<< "AvOve";
  cout<<setw(15)<< "AvUnd"<<endl;
  cout<<SP15<<AveWt;
  cout<<SP15<<ERela;
  cout<<SP15<<WtMax;
  cout<<SP15<<WtMin;
  cout<<SP15<<AvOve;
  cout<<SP15<<AvUnd<<endl;
}
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//         Endo of Class WtMonit                                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

