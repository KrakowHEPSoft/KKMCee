///////////////////////////////////////////////////////////////////////////////
//   This is ISR multi-photon generator
///////////////////////////////////////////////////////////////////////////////

#ifndef KKarFin_H
#define KKarFin_H

#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

#include "TString.h"
#include "TRandom.h"

#include "KKdbase.h"
#include "KKevent.h"
#include "KKbvir.h"

#include "BXFORMAT.h"
#include "TObject.h"
//________________________________________________________________________
class KKarFin: public TObject{
 public:
 ofstream *m_Out;     //! pointer to external Logfile for messages
// class objects
 TRandom  *m_RNgen;                // r.n.generator
 KKdbase  *DB;                     // Database
 KKevent  *m_Event;                // MC event record
 KKbvir   *m_BVR;                  // virtual formfactors
 // class member data
  public:
  int    m_KFfin;             // KF of final fermion
  int    m_IsFSR;             // is FSR active?
  double m_vvmin;             // input IR cutoff
  double m_delta;             // reduced IR cutoff
  double m_Emin;              // FSR/FSR absolute common IR cutoff
  double m_MasPhot;           // dummy photon mass
  double m_WtMass;
  double m_VoluMC;
  double m_YFS_IR;
  double m_YFSkon;
  // in order to avoid allocation for each event, move m_phot into m_Event?
  static const int maxPhot =   101;// max. num. KKMC photons +1
  TLorentzVector  m_phot[maxPhot]; // FSR photons, f77 indexing, m_phot[0] is SUM!
  TLorentzVector  m_q1, m_q2;      // final fermions
  TLorentzVector  m_r1, m_r2;      // final fermions
  int    m_nphot;                  // no. of photons
  double m_yfin[maxPhot];          // Sudakov variables
  double m_zfin[maxPhot];          // Sudakov variables
  //
  int    m_NevGen;                 // event counter for debug
  int    m_icont;                  // event counter for debug
  int    m_nfail;                  // control of maxPhot
  int    m_MarTot;                 // counting vetoed photons
//------------------------------------
// Obligatory members
  public:
  KKarFin();                    // explicit default constructor for streamer
  KKarFin(ofstream *OutFile);   // user constructor
  ~KKarFin();                   // explicit destructor
  public:
/////////////////////////////////////////////////////////////////////////////
// class member functions
double sqr( const double x );

void SetDB(    KKdbase *DBase){ DB = DBase;};
void SetRNgen( TRandom *RNgen){ m_RNgen = RNgen;};
void SetEvent( KKevent *Event){ m_Event = Event;};
void SetBVR(    KKbvir *BVR){   m_BVR = BVR;};

void Initialize();
void  Make(TLorentzVector *PX, double *wt);
void  YFSfin( TLorentzVector *PX, double Mas1, double Mas2, double CharSq, double *WtFin);
void  Piatek(double Mas1, double Mas2, double CharSq, double WtMlist[], double *Wt3);
void  PoissGen(double average, int *mult, double rr[]);
void  AngBre(double am2, double *ddel1, double *ddel2, double *costhg, double *sinthg, double *dist0, double *dist1);
void  MomPrint1(TString text, TLorentzVector &Vect);

////////////////////////////////////////////////////////////////////////////
       ClassDef(KKarFin,1); // Data base
};// KKarFin class
////////////////////////////////////////////////////////////////////////////
#endif
