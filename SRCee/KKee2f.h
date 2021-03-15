#ifndef KKee2f_H
#define KKee2f_H

#include <stdlib.h>
#include <math.h>
#include <complex>
using namespace std;

#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

#include "TString.h"

#include "TRandom.h"

#include "TLorentzVector.h"
#include "TFile.h"
#include "TH1.h"

#include "TMCgen.h"
#include "KKdbase.h"

#include "KKborn.h"
#include "KKdizet.h"
#include "KKevent.h"
#include "KKarLud.h"
#include "KKqed3.h"

#include "KKlasa.h"


class KKee2f: public TMCgen{
 public:
 KKdbase *DB;                     // Database
 KKborn  *m_BornDist;             // Born differential distribution
 KKdizet *m_EWtabs;               // EW formfactors
 KKevent *m_Event;                // MC event ISR+FSR in KKMC format
 KKarLud *m_GenISR;               // ISR YFS generator
// KKarFin *m_GenFSR;               // FSR YFS generator
 KKqed3  *m_QED3;                 // EEX matrix element
// KKceex  *m_GPS;                  // CEEX matrix element
// KKbvir  *m_BVR;                  // Library of virtual corrections
 KKlasa  *m_KKexamp;              // Template for new class

// Dimensionality
 static const int maxPar  = 10001;    // max. num. KKMC parameters +1
 static const int maxPhot =   101;    // max. num. KKMC photons +1
 static const int maxWT   =  1001;    // max. num. KKMC wt list +1
//----------------------------------
// INPUT data
 long     m_NevTot;           //   total numer of events to be generated, obsolete
 long     m_EventCounter;     //   event serial counter
 long     m_icont;            //   counter for debug printouts
 double   m_ypar[maxPar];     //   xpar input parameters of KKMC, c++ indexing
 double   m_xpar[maxPar];     //   xpar input parameters of KKMC, f77 indexing
 long     m_out;              //   output unit number for f77
// Seeds for Pseumar r.n. in KKMC
 int      m_ijkl_new;         // = 54217137;
 int      m_ntot_new;         // = 0;
 int      m_ntot2_new;        // = 0;
 int      m_iseed;                // KKMC random number generator seed
 //
 int      m_KFlist[6];            // list of final fermion KF indices to generate
 int      m_nKF;                  // number of fermions selected to generate

 double   m_MminCEEX[17];         // minimum fermion pair mass for CEEX

 int      m_nCallsFoam0;          // function calls during initialization
 double   m_XsPrim;               // primary cross section for overall normalization

 double   m_ParCIRCE[5];          // List of CIRCE parameters
 double   m_ParBES[5];            // List of BES  parameters
//-------------------------------------------
// Switches controlling operation
 int      m_FoamMode;             // Foam Density mode, <0 generation, >0 initialization
 int      m_RhoMode;              // Type of Density function =3,4,5
 int      m_Icont;                // density call counter for initialization
//-------------------------------------------
// some variables of the event
 TLorentzVector m_XXf;             // input for karfin
 double   m_CMSene;
 double   m_XXXene;
 double   m_Ebeam1;
 double   m_Ebeam2;
 double   m_CosTheta;
 int      m_KFini;
 int      m_KFfin;
 double   m_vv;
 double   m_r1;
 double   m_r2;
 // MC Weights
  double   m_WtSet[maxWT];           // list of alternative weights f77 indexing
 // double   m_WtSet[1001];             // list of alternative weights f77 indexing
  double   m_WtSetNew[maxWT];         // list of alternative weights f77 indexing
  double   m_WtAlter[maxWT];          // list of alternative weights f77 indexing
  double   m_WtMain;                  // Main weight
  double   m_WtCrude;                 // Crude weight
  double   m_WtFoam;                  // Foam weight
//-------------------------------------------
 public:
 KKee2f();                // explicit default constructor for streamer
 KKee2f(const char*);     // user constructor
 ~KKee2f();               // explicit destructor
 public:

/////////////////////////////////////////////////////////////////////////////
/// methods obligatory

  double sqr( const double x ){ return x*x;};
  //
  void   Initialize(TRandom*, ofstream*, TH1D*);

  void   ReaData(const char *DiskFile, int imax, double xpar[]);

  void   InitParams();
  void   FoamInitA();

  double Density(int, Double_t*);
  double RhoFoam5(double *Xarg);

  void   MakeGami(int KFini, double CMSene, double &gamiCR, double &gami, double &alfi);
  double MakeRhoISR(double gamiCR, double gami, double alfi, double vv, double vvmin, double vvmax);

  void   Generate();

  ////////////////////////////////////////////////////////////////////////////
      ClassDef(KKee2f,2); // Monte Carlo generator
  };
  ////////////////////////////////////////////////////////////////////////////
  #endif

