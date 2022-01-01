///////////////////////////////////////////////////////////////////////////////
//         Template of the class with ROOT persistency
///////////////////////////////////////////////////////////////////////////////

#ifndef TauPair_H
#define TauPair_H

#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

#include "BXFORMAT.h"
#include "TObject.h"


extern "C" {
//
// SUBROUTINE Taupair_Initialize(xpar)
   void taupair_initialize_(double[]);
   void taupair_getisinitialized_(const int&);
   void taupair_make1_();        // generates tau decay
   void taupair_imprintspin_();  // introduces spin effects by rejection
   void taupair_make2_();        // book-keeping, Photos, HepEvt
}//extern

//________________________________________________________________________
class TauPair: public TObject{
 public:
 ofstream *m_Out;     //! pointer to external Logfile for messages
// class member data
 public:
 double   CMSene;
//------------------------------------
// Obligatory members
  public:
  TauPair();                    // explicit default constructor for streamer
  TauPair(ofstream *OutFile);   // user constructor
  ~TauPair();                   // explicit destructor
  public:
/////////////////////////////////////////////////////////////////////////////
// class member functions
double sqr( const double x );

void Initialize(double[]);

int IsTauInitialized(){
  int IsTau;
  taupair_getisinitialized_(IsTau);
  return IsTau;
}

void Make(){
  taupair_make1_();        // generates tau decay
  taupair_imprintspin_();  // introduces spin effects by rejection
  taupair_make2_();        // book-keeping, Photos, HepEvt
}

////////////////////////////////////////////////////////////////////////////
       ClassDef(TauPair,1); // Data base
};// TauPair class
////////////////////////////////////////////////////////////////////////////
#endif
