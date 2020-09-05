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
extern "C" void kk2f_print1_(    const int&);
extern "C" void kk2f_fort_open_( const int&, const char*, int);
extern "C" void kk2f_fort_close_(const int&);
///////////////////////////////////////////////////////////////////////////////
extern "C" void kk2f_getwt_(     double&, double&);
extern "C" void kk2f_getxsecmc_( double&, double&);
extern "C" void kk2f_getbeams_(   double [], double []);
extern "C" void kk2f_getfermions_(double [], double []);
//extern "C" void kk2f_getnphot_(  int&);
extern "C" void kk2f_getnphot_(  int&);
extern "C" void kk2f_getphoton1_( int&, double []);
extern "C" void kk2f_getprimanorma_( double&, const int&);
extern "C" void kk2f_getxsnormpb_( double&, double&);
extern "C" void kk2f_getwtalter_( int&, double&);
//      SUBROUTINE KK2f_Make_WT
extern "C" void kk2f_make_wt_();
//         SUBROUTINE BornV_SetQEDmodif(lambda)
extern "C" void bornv_setqedmodif_( double&);
///////////////////////////////////////////////////////////////////////////////
extern "C" void karlud_getfermions_(double [], double []);
///////////////////////////////////////////////////////////////////////////////
//extern "C" void pylist_(const int&);
extern "C" void pygive_(const char *directive, int s1);
extern "C" void hepevt_getkffin_(  int&);
//////////////////////////////////////////////////////////////////////////////

extern "C" void pseumar_initialize_( const int&, const int&, const int&);

///////////////////////////////////////////////////////////////////////////////

/*[[[[
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
*/

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
  // Seeds for Pseumar r.n. in KKMC
  int       m_ijkl_new;       // = 54217137;
  int       m_ntot_new;       // = 0;
  int       m_ntot2_new;      // = 0;

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
  double Density(int, double *){return 0;};   /// Dummy method of the abstract class TFOAM_INTEGRAND
/////////////////////////////////////////////////////////////////////////////
  //void Initialize(double []);
  //void Finalize(double&,  double&);
///////////////////////////////////////////////////////////////////////////////
  void Print1();
//[[[  long GetPyNpart();
//[[[  void GetPyParticle( const long, TPartLund &Event);
///////////////////////////////////////////////////////////////////////////////
  void GetWt( double&,  double&);
  void GetBeams(    TLorentzVector&,  TLorentzVector&);
  void GetFermions( TLorentzVector&,  TLorentzVector&);
  void GetPhoton1(const int, TLorentzVector&);
  void GetFermKarlud( TLorentzVector &,  TLorentzVector &);
//  void GetPhoton1(int & , TLorentzVector&);
//  void GetNphot(  int &);
  void GetNphot(  int &);
  void GetXsecMC( double &,  double &);
  void GetPrimaNorma( double &,  int &);
  double GetWtAlter(const int );
//  double GetWtAlter(int &);
///////////////////////////////////////////////////////////////////////////////
//  void PyList(int);
  void PyGive(const char *directive);
  void GetKFfin(  int &);
  void ReaData(const char*, int, double[]);
////////////////////////////////////////////////////////////////////////////
    ClassDef(TMCgenKKMC,2); // Monte Carlo generator
};
////////////////////////////////////////////////////////////////////////////
#endif
