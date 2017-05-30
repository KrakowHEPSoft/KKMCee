#ifndef KKMC_H
#define KKMC_H

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
extern "C" void kk2f_print1_(    const int);
extern "C" void kk2f_fort_open_( const int&, const char*, int);
extern "C" void kk2f_fort_close_(const int&);
///////////////////////////////////////////////////////////////////////////////
extern "C" void kk2f_getwt_(     double&, double&);
extern "C" void kk2f_getxsecmc_( double&, double&);
extern "C" void kk2f_getbeams_(   double [], double []);
extern "C" void kk2f_getfermions_(double [], double []);
extern "C" void kk2f_getnphot_(  int&);
extern "C" void kk2f_getphoton1_( int&, double []);
extern "C" void kk2f_getprimanorma_( double&, const int&);
extern "C" void kk2f_getxsnormpb_( double&, double&);
extern "C" void kk2f_getwtalter_( int&, double&);
///////////////////////////////////////////////////////////////////////////////
extern "C" void karlud_getfermions_(double [], double []);
///////////////////////////////////////////////////////////////////////////////
extern "C" void pylist_(const int&);
extern "C" void pygive_(const char *directive, int s1);
extern "C" void hepevt_getkffin_(  int&);
//
///////////////////////////////////////////////////////////////////////////////


class KKMC{
//class KoralwMaker {
 public:
  int      m_NevTot;         //   total numer of events to be generated
  int      m_EvenCounter;    //   event serial counter
  double    m_ypar[10001];    //   xpar input/output aprameters of koralw
  int      m_out;            //   output unit number
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
//  void ReadData( int &ntot );
  void Print1();
///////////////////////////////////////////////////////////////////////////////
  int GetPyNpart();
  void GetPyParticle( const int, PartLund &Event);
///////////////////////////////////////////////////////////////////////////////
  void GetWt( double&,  double&);
  void GetBeams(    TLorentzVector&,  TLorentzVector&);
  void GetFermions( TLorentzVector&,  TLorentzVector&);
  void GetPhoton1(const int, TLorentzVector&);
  void GetFermKarlud( TLorentzVector &,  TLorentzVector &);
//  void GetPhoton1(int & , TLorentzVector&);
  void GetNphot(  int &);
  void GetXsecMC( double &,  double &);
  void GetPrimaNorma( double &,  int &);
  double GetWtAlter(const int );
//  double GetWtAlter(int &);
///////////////////////////////////////////////////////////////////////////////
  void PyList(int);
  void PyGive(const char *directive);
  void GetKFfin(  int &);
};
////////////////////////////////////////////////////////////////////////////
#endif
