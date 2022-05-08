///////////////////////////////////////////////////////////////////////////////
//   This is ISR multi-photon generator
///////////////////////////////////////////////////////////////////////////////

#ifndef KKarLud_H
#define KKarLud_H

#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

#include "TString.h"
#include "TRandom.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "KKdbase.h"
#include "KKevent.h"

#include "BXFORMAT.h"
#include "TObject.h"
//________________________________________________________________________
class KKarLud: public TObject{
 public:
 ofstream *m_Out;     //! pointer to external Logfile for messages
// class objects
 KKdbase  *DB;                     //! Database
 TRandom  *m_RNgen;                //! r.n.generator
 KKevent  *m_Event;                //! MC event record
// class member data
 public:
 static const int maxPhot =   101; // max. num. KKMC photons +1
 double   m_vvmin;                 // IR cutoff
 // variables event dependent
 double   m_CMSene;                // reduced CMS energy
 double   m_XXXene;                // reduced CMS energy
 double   m_amel;                  // mass of beam fermion
 int      m_KFini;                 // incoming fermion ID
 int      m_KFfin;                 // outgoing fermion ID
 double   m_vv;                    // vv from Foam
 double   m_r1, m_r2;              // z variables from Foam
 double   m_AvMult;                // average photon multiplicity provided by Foam
 int      m_nphot;                 // Photon multiplicity
 double   m_yini[maxPhot];         // Sudakov y-variables
 double   m_zini[maxPhot];         // Sudakov z-variables
 int      m_icont;                 // event counter for debug
 int      m_nfail;                 // control of maxPhot
// Obligatory members
  public:
  KKarLud();                    // explicit default constructor for streamer
  KKarLud(ofstream *OutFile);   // user constructor
  ~KKarLud();                   // explicit destructor
  public:
/////////////////////////////////////////////////////////////////////////////
// class member functions
double sqr( const double x );

void SetDB(    KKdbase *DBase){ DB = DBase;};
void SetRNgen( TRandom *RNgen){ m_RNgen = RNgen;};
void SetEvent( KKevent *Event){ m_Event = Event;};

//void SetFoamVars(int KFini, int KFfin, double vv, double z1, double z2){
//	m_KFini=KFini;  m_KFfin=KFfin ;
//	m_vv= vv;  m_r1=z1;  m_r2=z2;
//};

//void SetAvMult( double AvMult){ m_AvMult = AvMult;};

void Make(TLorentzVector *PX, double *wt_ISR);

void Initialize();
void PoissGen(double average, int *mult, double rr[]);

void AngBre(double am2,
		double *del1, double *del2, double *costhg, double *sinthg, double *dist0, double *dist1);

void YFSgen(double XXXene, double vv, TLorentzVector *PX, double *WtIni);


void MomPrint1(TString text, TLorentzVector &Vect);
//void MomPrint1(TLorentzVector &Vect);

////////////////////////////////////////////////////////////////////////////
       ClassDef(KKarLud,1); // Data base
};// KKarLud class
////////////////////////////////////////////////////////////////////////////
#endif
