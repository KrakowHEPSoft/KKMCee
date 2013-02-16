//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//               Class   KKMC                                               //
//                                                                          //
//              Interface (wrapper)  to MC event generator KKMC             //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
#include "KKMC.h"


//////////////////////////////////////////////////////////////////////////////
//========================================================
//   COMMON/LUJETS/N,     K(4000,5),P(4000,5),V(4000,5)
//   COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
//--------------------------------------------------------
typedef struct {
//    long   n; 
//    long   npad;
//    long   k[5][4000];
    int    n; 
    int    npad;
    int    k[5][4000];
    double p[5][4000];
    double v[5][4000];
                } CommonPYJETS;
extern CommonPYJETS &pyjets ;
extern CommonPYJETS pyjets_ ;
#define cb_PYjets pyjets_
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//============================================================
//  COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
//  COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
//------------------------------------------------------------
typedef struct {
//    long   mstu[200]; 
    int    mstu[200]; 
    double paru[200];
//    long   mstj[200];
    int    mstj[200];
    double parj[200];
                } CommonPYDAT1;
extern CommonPYDAT1 &pydat1 ;
extern CommonPYDAT1 pydat1_ ;
#define cb_PYdat1 pydat1_
//////////////////////////////////////////////////////////////////////////////


///      SUBROUTINE KK2f_Initialize(xpar)
extern "C" void kk2f_initialize_(double xpar[]);
extern "C" void kk2f_make_();
extern "C" void kk2f_finalize_();
extern "C" void kk2f_print1_(    const long&);
extern "C" void kk2f_fort_open_( const long&, char*, long);
extern "C" void kk2f_fort_close_(const long&);
///////////////////////////////////////////////////////////////////////////////
extern "C" void kk2f_getwt_(     double&, double&);
extern "C" void kk2f_getxsecmc_( double&, double&);
extern "C" void kk2f_getbeams_(   double [], double []);
extern "C" void kk2f_getfermions_(double [], double []);
extern "C" void kk2f_getnphot_(  long&);
extern "C" void kk2f_getphoton1_( long&, double []);
extern "C" void kk2f_getprimanorma_( double&, const int&);
extern "C" void kk2f_getxsnormpb_( double&, double&);
extern "C" void kk2f_getwtalter_( long&, double&);
///////////////////////////////////////////////////////////////////////////////
extern "C" void pylist_(const long&);
extern "C" void pygive_(char *directive, long s1);
//
///////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
void KKMC::Initialize(double ypar[])
{
  m_NevTot = 0;
  m_EvenCounter = 0;

// Initialisation input data before generation
  cout << "*****************************" << endl;
  cout << "**   KKMC   Initialize     **" << endl;
  cout << "*****************************" << endl;

  //=============================================================
  //   opening disk fime for fortran part of code
  m_out = ypar[3];
  char *output_file = "./pro.output";
  long sl2 = strlen(output_file);
  kk2f_fort_open_(m_out,output_file,sl2);

  //*******************//
  kk2f_initialize_(ypar);
  //*******************//

  //PyGive("MDCY(113,1)=1;");            // allow rho0 decay
  //PyGive("MDCY(333,1)=0;");            // inhibit phi decay
  //PyGive("MDCY(111,1)=0;");            // inhibit pi0 decay
  //PyGive("MSTP(41)=1;");               // allow all resoanace decay (2->1)
  //PyGive("MSTP(41)=0;");               // inhibit all resoanace decay (2->1)

  PyGive("MSTU(21)=1;");               // no stop due to errors !!!!

}
///////////////////////////////////////////////////////////////////////////////
void KKMC::Make()
{
  kk2f_make_();
}
///////////////////////////////////////////////////////////////////////////////
void KKMC::Finalize( double &XsNormPb, double &XsErroPb)
{
// XsNormPb +- XsErroPb is xsection to be readily used to normalize
// histograms both for wt=1 and wt-ed events

  cout << "*****************************" << endl;
  cout << "**   KKMC   Finalize       **" << endl;
  cout << "*****************************" << endl;
  //
  kk2f_finalize_();
  //
  kk2f_getxsnormpb_( XsNormPb, XsErroPb);
  //
  kk2f_fort_close_(m_out);
}
///////////////////////////////////////////////////////////////////////////////
long KKMC::GetPyNpart()
{
// provides no. of entries in Lund/Pythia common block
  long npart = cb_PYjets.n;
  //cout<<"KKMC::Make: npart ="<<npart<<endl;
  return npart;
}
///////////////////////////////////////////////////////////////////////////////
void KKMC::GetPyParticle( const long j, PartLund &Event)
{
// Export one particle from /PYJETS/
  Event= PartLund(j+1,
         cb_PYjets.k[0][j], cb_PYjets.k[1][j],
	 cb_PYjets.k[2][j], cb_PYjets.k[3][j],  cb_PYjets.k[4][j],
	 cb_PYjets.p[0][j], cb_PYjets.p[1][j],
	 cb_PYjets.p[2][j], cb_PYjets.p[3][j],  cb_PYjets.p[4][j],
	 cb_PYjets.v[0][j], cb_PYjets.v[1][j],
	 cb_PYjets.v[2][j], cb_PYjets.v[3][j],  cb_PYjets.v[4][j]);
  //cout<<" here we are!!! j= "<<j<<endl;
  //Event.Print(1);
}
///////////////////////////////////////////////////////////////////////////////
void KKMC::Print1()
{
// print event using KKMC
  kk2f_print1_(m_out);  
}
///////////////////////////////////////////////////////////////////////////////
void KKMC::PyList(long lev)
{
// print event using pythia
  PyGive("MSTU(11)=16");
  pylist_(lev);
  PyGive("MSTU(11)=6");
  pylist_(lev);
}
///////////////////////////////////////////////////////////////////////////////
void KKMC::PyGive(char *directive)
{
// set pythia directive
  long s1;
  s1 = strlen(directive);
  pygive_(directive, s1);
}
///////////////////////////////////////////////////////////////////////////////
void KKMC::GetWt( double &WtMain,  double &WtCrude)
{
// get MC weight of event
  kk2f_getwt_(WtMain, WtCrude);
}
///////////////////////////////////////////////////////////////////////////////
void KKMC::GetBeams( TLorentzVector &B1,  TLorentzVector &B2)
{
// get 4-momenta of beams
  double p1[4];
  double p2[4];
  kk2f_getbeams_(p1,p2);
  B1.SetPxPyPzE(p1[0],p1[1],p1[2],p1[3]);
  B2.SetPxPyPzE(p2[0],p2[1],p2[2],p2[3]);
}
///////////////////////////////////////////////////////////////////////////////
void KKMC::GetFermions( TLorentzVector &F1,  TLorentzVector &F2)
{
// get 4-momenta of final fermions
  double p1[4];
  double p2[4];
  kk2f_getfermions_(p1,p2);
  F1.SetPxPyPzE(p1[0],p1[1],p1[2],p1[3]);
  F2.SetPxPyPzE(p2[0],p2[1],p2[2],p2[3]);
  //for(int i=0;i<4;i++) cout<<B1[i]<<"  "; cout<<endl;
  //for(int i=0;i<4;i++) cout<<B2[i]<<"  "; cout<<endl;
}
///////////////////////////////////////////////////////////////////////////////
void KKMC::GetXsecMC( double &xSecPb,  double &xErrPb)
{
// get MC final xsection PICOBARNS
  kk2f_getxsecmc_(xSecPb, xErrPb);
}
///////////////////////////////////////////////////////////////////////////////
void KKMC::GetPrimaNorma(double &XsPrim, long &NevPrim)
{
// get normalization elements NANOBARNS
  int NevPrim1;
  kk2f_getprimanorma_( XsPrim, NevPrim1);
  NevPrim=NevPrim1;
}
///////////////////////////////////////////////////////////////////////////////
void KKMC::GetPhoton1(const long iphot, TLorentzVector &phot)
{
// get one photon 4-vector from KKMC
  double p1[4];
  long iphot1=iphot;
  kk2f_getphoton1_(iphot1, p1);
  phot.SetPxPyPzE(p1[0],p1[1],p1[2],p1[3]);
}
///////////////////////////////////////////////////////////////////////////////
void KKMC::GetNphot(long &Nphot)
{
// get photon multiplicity from KKMC
  kk2f_getnphot_( Nphot);
}
///////////////////////////////////////////////////////////////////////////////
double KKMC::GetWtAlter(const long id)
{
  double WtAlter;
  long id1 =id;
  kk2f_getwtalter_( id1, WtAlter);
  return WtAlter;
}
//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//              End of  Class   KKMC                                        //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
