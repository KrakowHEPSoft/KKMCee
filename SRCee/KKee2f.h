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

#include "KKlasa.h"


class KKee2f: public TMCgen{
//class KKee2f: public TFOAM_INTEGRAND {
 public:
 KKdbase *DB;                     // Database
 KKborn  *m_BornDist;             // Born differential distribution
// KKdizet *m_EWforms;              // EW formfactors
 KKdizet *m_EWtabs;               // EW formfactors
 KKlasa  *m_KKexamp;              // Template for new class

// Dimensionality
 static const int maxPar  = 10001;    // max. num. KKMC parameters +1
 static const int maxPhot =   101;    // max. num. KKMC photons +1
 static const int maxWT   =  1001;    // max. num. KKMC wt list +1
 double   m_pi;               // pi
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
 int      m_iseed;                 // KKMC random number generator seed

//-------------------------------------------
// Switches controlling operation
 int      m_FoamMode;              // Foam Density mode, <0 generation, >0 initialization
 int      m_RhoMode;               // Type of Density function =3,4,5
 int      m_Icont;                 // density call counter for initialization


 public:
 KKee2f();                // explicit default constructor for streamer
 KKee2f(const char*);     // user constructor
 ~KKee2f();               // explicit destructor
 public:

/////////////////////////////////////////////////////////////////////////////
/// methods obligatory


  double sqr( const Double_t x ){ return x*x;};

  double Density(int, Double_t*);
  //
  void Initialize(TRandom*, ofstream*, TH1D*);

  void ReaData(const char *DiskFile, int imax, double xpar[]);
  ////////////////////////////////////////////////////////////////////////////
      ClassDef(KKee2f,2); // Monte Carlo generator
  };
  ////////////////////////////////////////////////////////////////////////////
  #endif

