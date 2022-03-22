///////////////////////////////////////////////////////////////////////////////
//  Replacement for HepEvt.f
//  Class with ROOT persistency
//
///////////////////////////////////////////////////////////////////////////////

#ifndef HepFace_H
#define HepFace_H

#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

#include "KKevent.h"

#include "HepMC3/GenEvent.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/Print.h"
#include "HepMC3/Selector.h"

using namespace HepMC3;

#include "BXFORMAT.h"
#include "TObject.h"
//________________________________________________________________________
class HepFace: public TObject{
 public:
 ofstream *m_Out;     //! pointer to external Logfile for messages
// class member data
 public:
 KKevent   *m_Event;           //!  MC event ISR+FSR in KKMC format
 GenEvent  *m_Hvent;           //! HEPMC3 event (no persistency)
//------------------------------------
// Obligatory members
  public:
  HepFace();                    // explicit default constructor for streamer
  HepFace(ofstream *OutFile);   // user constructor
  ~HepFace();                   // explicit destructor
  public:
/////////////////////////////////////////////////////////////////////////////
// class member functions
double sqr( const double x );

void SetEvent( KKevent  *Event){ m_Event = Event;};
void SetHvent( GenEvent *Event){ m_Hvent = Event;};
FourVector Vect4( TLorentzVector X){ return FourVector(X.Px(),X.Py(),X.Pz(),X.E()); };

void Initialize();
void make1();
void FillHep3(int N, int IST, int ID, int JMO1, int JMO2, int JDA1, int JDA2, float P4[], float &PINV, bool PHFLAG);
////////////////////////////////////////////////////////////////////////////
       ClassDef(HepFace,1); // Data base
};// HepFace class
////////////////////////////////////////////////////////////////////////////
#endif
