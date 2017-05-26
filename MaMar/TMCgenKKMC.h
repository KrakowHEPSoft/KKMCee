#ifndef TMCgenKKMC_H
#define TMCgenKKMC_H

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//               Class   TMCgenKKMC                                         //
//                                                                          //
//              Interface (wrapper)  to MC event generator KKMC             //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

#include "TMCgen.h"
#include "TLorentzVector.h"
#include "TPartLund.h"


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
extern "C" void karlud_getfermions_(double [], double []);
///////////////////////////////////////////////////////////////////////////////
extern "C" void pylist_(const long&);
extern "C" void pygive_(char *directive, long s1);
extern "C" void hepevt_getkffin_(  long&);
//
///////////////////////////////////////////////////////////////////////////////

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


//________________________________________________________________________
class TMCgenKKMC: public TMCgen
{
 public:
  long      m_NevTot;         //   total numer of events to be generated
  long      m_EvenCounter;    //   event serial counter
  int       m_jmax;           //   lenght of ymar
  double    m_ypar[10001];    //   xpar input/output aprameters of KKMC
  double    m_xpar[10001];    //   xpar input/output aprameters of KKMC
  double    m_XsNormPb;       //   normalization
  double    m_XsErroPb;       //   normalization
  long      m_out;            //   output unit number
 public:
//------ header of event_lu constructor -------
 public:
 TMCgenKKMC();                // explicit default constructor for streamer
 TMCgenKKMC(const char*);     // user constructor
 ~TMCgenKKMC();               // explicit destructor
 public:
/////////////////////////////////////////////////////////////////////////////
/// methods obligatory
  void Initialize(TRandom*, ofstream*, TH1D*);
  void Finalize();
  //void Make();
  void Generate();
  double Density(int, double *){;};   /// Dummy method of the abstract class TFOAM_INTEGRAND
/////////////////////////////////////////////////////////////////////////////
  //void Initialize(double []);
  //void Finalize(double&,  double&);
///////////////////////////////////////////////////////////////////////////////
  void Print1();
  long GetPyNpart();
  void GetPyParticle( const long, TPartLund &Event);
///////////////////////////////////////////////////////////////////////////////
  void GetWt( double&,  double&);
  void GetBeams(    TLorentzVector&,  TLorentzVector&);
  void GetFermions( TLorentzVector&,  TLorentzVector&);
  void GetPhoton1(const long, TLorentzVector&);
  void GetFermKarlud( TLorentzVector &,  TLorentzVector &);
//  void GetPhoton1(long & , TLorentzVector&);
  void GetNphot(  long &);
  void GetXsecMC( double &,  double &);
  void GetPrimaNorma( double &,  long &);
  double GetWtAlter(const long );
//  double GetWtAlter(long &);
///////////////////////////////////////////////////////////////////////////////
  void PyList(long);
  void PyGive(char *directive);
  void GetKFfin(  long &);
  void ReaData(char DiskFile[], int imax, double xpar[]);
////////////////////////////////////////////////////////////////////////////
    ClassDef(TMCgenKKMC,2); // Monte Carlo generator
};
////////////////////////////////////////////////////////////////////////////
#endif
