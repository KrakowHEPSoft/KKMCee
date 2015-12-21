//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//               Class   KKMC                                               //
//                                                                          //
//              Interface (wrapper)  to MC event generator KKMC             //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

#include "TLorentzVector.h"
#include "PartLund.h"


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
extern "C" void hepevt_getkffin_(  long&);
//
///////////////////////////////////////////////////////////////////////////////


class KKMC{
//class KoralwMaker {
 public:
  long      m_NevTot;         //   total numer of events to be generated
  long      m_EvenCounter;    //   event serial counter
  double    m_ypar[10001];    //   xpar input/output aprameters of koralw
  long      m_out;            //   output unit number
 public:
//------ header of event_lu constructor -------
  KKMC(){;};
  ~KKMC(){;};
 public:
/////////////////////////////////////////////////////////////////////////////
//                    generation of basic MC event                         //
//------ generator initialization
  void Initialize(double []);
//------ generator finalization
  void Finalize(double&,  double&);
//------ read input data
  void Make();
//  void ReadData( long &ntot );
  void Print1();
///////////////////////////////////////////////////////////////////////////////
  long GetPyNpart();
  void GetPyParticle( const long, PartLund &Event);
///////////////////////////////////////////////////////////////////////////////
  void GetWt( double&,  double&);
  void GetBeams(    TLorentzVector&,  TLorentzVector&);
  void GetFermions( TLorentzVector&,  TLorentzVector&);
  void GetPhoton1(const long, TLorentzVector&);
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
};
////////////////////////////////////////////////////////////////////////////
