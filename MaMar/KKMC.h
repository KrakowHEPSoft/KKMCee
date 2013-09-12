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
};
////////////////////////////////////////////////////////////////////////////
