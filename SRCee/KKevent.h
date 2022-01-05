///////////////////////////////////////////////////////////////////////////////
//         Template of the class with ROOT persistency
///////////////////////////////////////////////////////////////////////////////

#ifndef KKevent_H
#define KKevent_H

#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

#include "BXFORMAT.h"
#include "TObject.h"

#include "TLorentzVector.h"


//________________________________________________________________________
class KKevent: public TObject{
 public:
 ofstream *m_Out;     //! pointer to external Logfile for messages
// class member data
 public:
 static const int maxPhot =   101;   // max. num. KKMC photons +1
 static const int maxWT   =  1001;   // max. num. KKMC wt list +1
 double         m_CMSene;            // hadron-hadron CMS energy
 double         m_XXXene;            // parto-parton  CMS energy
 // internal MC generation variables changing event per event
 long     m_EventCounter;            // event serial counter
 int      m_KFini;                   // incoming fermion ID (negative if 1st is anti)
 int      m_KFfin;                   // outgoing fermion ID (negative if 1st is anti)
 int      m_HasFSR;                  // =0,1 FSR absent,present
 double   m_Xenph;                   // ISR enhancement in case of IFI, otherwise =1
 double   m_Mbeam1, m_Mbeam2;        // masses of parton (quark) beams
 double   m_vv;                      // vv ISR variable of KKMC
 double   m_CosTheta;                // cos(theta) at the crude level
 double   m_r1, m_r2;                // z variables for initial Energy Beam Spread (EBS)
 // event 4-momenta
 TLorentzVector m_Pf1;               // initial quark m_p1
 TLorentzVector m_Pf2;               // initial quark m_p2
 TLorentzVector m_Qf1;               // final fermion m_q1
 TLorentzVector m_Qf2;               // final fermion m_q2
 TLorentzVector m_PX;                // Z boson after ISR before FSR (unphysical)
 TLorentzVector m_Rem1;              // (BES) remnant 1
 TLorentzVector m_Rem2;              // (BES) remnant 2
 int            m_nPhot;             // No of generated photons ISR+FSR
 int      m_nPhotISR, m_nPhotFSR;    // No of generated photons ISR, FSR
 int            m_isr[maxPhot];      // marker of ISR photons, f77 indexing
 TLorentzVector m_PhotAll[maxPhot];  // all ISR+FSR photons,   f77 indexing
 TLorentzVector m_PhotISR[maxPhot];  // ISR photons, f77 indexing
 TLorentzVector m_PhotFSR[maxPhot];  // FSR photons, f77 indexing
 //------------------------------------
 // Formfactors
 double   m_AvMult;                // Average ISR photon multiplicity
 double   m_YFS_IR_ini;            // YFS formfactor part from Foam Density
 double   m_YFSkon_ini;            // YFS formfactor part from Foam Density
 double   m_YFS_IR_fin;            // YFS formfactor part from KKarfin (Piatek)
 double   m_YFSkon_fin;            // YFS formfactor part from KKarfin (Piatek)
 // MC Weights
 double   m_WtFoam;                // Foam weight
 double   m_WT_ISR;                // ISR weight
 double   m_WT_FSR;                // FSR weight
 double   m_WtCrude;               // Crude weight
 double   m_WtMain;                // Main weight
 double   m_WtSet[maxWT];          // list of alternative weights f77 indexing
 double   m_WtSetNew[maxWT];       // list of alternative weights f77 indexing
 double   m_WtAlter[maxWT];        // list of alternative weights f77 indexing
//------------------------------------
// Obligatory members
  public:
  KKevent();                    // explicit default constructor for streamer
  KKevent(ofstream *OutFile);   // user constructor
  ~KKevent();                   // explicit destructor
  public:
/////////////////////////////////////////////////////////////////////////////
// class member functions
double sqr( const double x );

void Initialize(double );
void Merge();

void DefPair(double CMSene, double qm1, double qm2, TLorentzVector *Pf1, TLorentzVector *Pf2);
void GivePair(double CMSene, double qm1, double qm2,
     TLorentzVector *Pf1, TLorentzVector *Pf2, double *beta, double *eta1, double *eta2 );


void RotEuler(double alf, double bet, double gam, TLorentzVector *Pf);
void RotEul(  double the, double phi, TLorentzVector *Pf);
void BoostQ(int Mode, TLorentzVector *Q, TLorentzVector *Pf);

void BoostEul(double the, double phi, TLorentzVector *QQk, TLorentzVector *PX, TLorentzVector *pvec);
void RotEulInv(double the, double phi, TLorentzVector *Pf);

void PhaSpac2(TLorentzVector *PX, double the, double phi, double amfin,
		      TLorentzVector *Qf1, TLorentzVector *Qf2);

void ThetaR(double *cth11, double *cth12, double *cth21, double *cth22);
void ThetaD(TLorentzVector &PX, double &CosTheta);

void ZBoostALL();
void ZReverseALL();

void MomPrint( TLorentzVector &Vect);
void PrintISR();
void PrintISR_FSR();
//
void MomPrint(ofstream *Out, TLorentzVector &Vect);
void PrintISR_FSR(ofstream *Out);
//
void EventPrintAll();
void EventPrintAll(ofstream *Out);
//
////////////////////////////////////////////////////////////////////////////
       ClassDef(KKevent,1); // Data base
};// KKevent class
////////////////////////////////////////////////////////////////////////////
#endif
