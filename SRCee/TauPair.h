/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
//                         CLASS  TauPair                                          //
//       Purpose:                                                                  //
//       (a) Interface to Simulation of TWO tau decays                             //
//       (b) Calculates spin weight wt using KKceex::MakeRho2 and introduces       //
//           spin effects in tau decays by means of rejection with <wt>=1          //
//       (c) Transforms decay products to CMS frame                                //
//       (d) Interfaces Photos to Tauola                                           //
//                                                                                 //
//   Notes:                                                                        //
//   The class is initialized by KKee2f::Initialize                                //
//   It is called from KKee2f::generate                                            //
//   It needs KKceex to be initialized in order to calculate spin weight (final)   //
//                                                                                 //
//   For the moment this file contains the interface to tauola and Photo           //
//   The rest of code is in tauface.f and tauola.f                                 //
//                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////

#ifndef TauPair_H
#define TauPair_H

#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

#include "BXFORMAT.h"
#include "TObject.h"
#include "TLorentzVector.h"
#include "TRandom.h"

#include "KKceex.h"
#include "KKevent.h"
#include "KKdbase.h"

#include "HepMC3/GenEvent.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/Print.h"
#include "HepMC3/Selector.h"
using namespace HepMC3;

// PHOTOS header files
#include "Photos/Photos.h"
#include "Photos/PhotosHepMC3Event.h"
#include "Photos/Log.h"
using namespace Photospp;

extern "C" {
/////////////////////////////
// HepEvt
  void hepevt_getf_( const int&);                    // fermion is here
  void hepevt_getfbar_(const int&);                  // antifermion is here
// tauface
  void tauface_setfermpos_(const int&, const int&);  // set ffbar positions in Tauola
  void tauface_print_();                    // printing event using pythia
// PYTHIA
  void pyhepc_(const int&);     // HepEvt-->Pythia
  void pylist_(const int&);
  void pygive_(const char*, long int);
/////////////////////////////////////////
//      TAUOLA
// SUBROUTINE INIETC(ITAUXPAR,xpar)
  void inietc_( int*, double[]);
// SUBROUTINE INIMAS(ITAUXPAR,xpar)
  void inimas_( int *ITAUXPAR, double xpar[]);
// SUBROUTINE INITDK(ITAUXPAR,xpar)
  void initdk_( int *ITAUXPAR, double xpar[]);
// SUBROUTINE INIPHY(XK00)
  void iniphy_( double *xk0qed);
//  SUBROUTINE DEKAY(KTO,HX)
  void dekay_( int* , double Hvec[]);
/////////////////////////////////////////
//       PHOTOS
//  SUBROUTINE PHOINI
  void phoini_();
  void photos_(int&);
}//extern


//________________________________________________________________________
class TauPair: public TObject{
// class member data
 public:
 ofstream  *m_Out;             //! pointer to external Logfile for messages
 KKdbase   *DB;                // Database
 KKceex    *m_GPS;             //  CEEX matrix element
 KKevent   *m_Event;           //!  MC event ISR+FSR in KKMC format (no persistency)
 GenEvent  *m_Hvent;           //! HEPMC3 event (no persistency)!
 TRandom   *m_RNgen;           //  central r.n. generator
 //
 double     m_HvecTau1[4];     //! Spin Polarimeter vector first  Tau
 double     m_HvecTau2[4];     //! Spin Polarimeter vector second Tau
 double     m_HvClone1[4];     //! Clone Spin Polarimeter vector first  Tau
 double     m_HvClone2[4];     //! Clone Spin Polarimeter vector second Tau
 TLorentzVector  m_H1;     //!
 TLorentzVector  m_H2;     //!
 TLorentzVector  m_PP;     //!
// TLorentzVector  m_HLveclo1;   //!
// TLorentzVector  m_HLveclo2;   //!
 double     m_beta1;           // Random Euler angle for cloning 1-st tau
 double     m_alfa1;           // Random Euler angle for cloning 1-st tau
 double     m_gamma1;          // Random Euler angle for cloning 1-st tau
 double     m_beta2;           // Random Euler angle for cloning 2-nd tau
 double     m_alfa2;           // Random Euler angle for cloning 2-nd tau
 double     m_gamma2;          // Random Euler angle for cloning 2-nd tau
 double     m_phi1;            // phi   of HvecTau1
 double     m_thet1;           // theta of HvecTau1
 double     m_phi2;            // phi   of HvecTau2
 double     m_thet2;           // theta of HvecTau2
 int        m_IsInitialized;   // super key, for inhibiting all tauola activity
 int        m_itdkRC;          // old key for QED in tauola leptonic decays
 int        m_KeyClone;        // switch for cloning procedure =1,2
//------------------------------------
// Obligatory members
  public:
  TauPair();                    // explicit default constructor for streamer
  TauPair(ofstream *OutFile);   // user constructor
  ~TauPair();                   // explicit destructor
/////////////////////////////////////////////////////////////////////////////
// class member functions
double sqr( const double x );

void SetDB(    KKdbase  *DBase){ DB = DBase;};
void SetEvent( KKevent  *Event){ m_Event = Event;};
void SetHvent( GenEvent *Hvent){ m_Hvent = Hvent;};
void SetGPS(   KKceex   *GPS){   m_GPS   = GPS;};

void Initialize(double[]);

int IsTauInitialized(){
//  int IsTau;
//  taupair_getisinitialized_(IsTau);
  return m_IsInitialized;
}
void SetRNgen(TRandom *RNgen){ m_RNgen= RNgen;};

void Make1();
void Clone();
void ImprintSpin();
void Make2();
void RunPhotosPP();
void Tralo4(int Kto, float P[], float Q[], float &AM);
void Finalize();
////////////////////////////////////////////////////////////////////////////
       ClassDef(TauPair,1); // Data base
};// TauPair class
////////////////////////////////////////////////////////////////////////////
#endif
