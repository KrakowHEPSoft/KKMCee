///////////////////////////////////////////////////////////////////////////////
//                NAMING CONVENTION
// The class has no local variables hence prefix m_ for data members is
// not necessary. Moreover in all object refering to data members of THIS class
// prefix like DB->CMSene will be there anyway.
///////////////////////////////////////////////////////////////////////////////

#ifndef KKdbase_H
#define KKdbase_H

#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

#include "TMCgen.h"
//________________________________________________________________________
class KKdbase: public TObject{
 public:
 static const int maxPar  = 10001;    // max. num. KKMC parameters +1
 double   m_xpar[maxPar];      //   xpar input parameters of KKMC, f77 indexing
//----------------------------------
 ofstream *m_Out;              //! External Logfile for messages
// INPUT data
 public:
 double   CMSene;
 int      KeyGPS;
 int      KeyINT;
 int      KeyISR;
 int      KeyFSR;
 int      KeyBES;
 int      KeyWgt;
 int      KeyElw;
 int      KeyZet;
 int      KeyWtm;
 int      KeyPia;
 int      KeyThe;
 int      KeyPhts;
 // FOAM
 int      Foam_nCells;
 int      Foam_Vopt;
 int      Foam_nSampl;
 int      Foam_nBins;
 int      Foam_EvPerBin;
 int      Foam_OptRej;
 int      Foam_OptDrive;
 double   Foam_WtMaxRej;
// Electroweak
 double   MZ;
 double   GamZ;
 double   swsq;
 double   MW;
 double   GamW;
 double   GFermi;
 int      KeyMasQ;
 double   MasPhot;
 double   WTmax;
 double   delfac;
 double   vvmin;
 double   vvmax;
 double   Xenph;
 double   Vcut[3];
 double   VQcut;
 double   Alfinv0;
 double   AlfinvZ;
 int      KeyFixAlf;
 double   XXXmin;
 double   XXXmax;
 int      Nalgo;
 double   gnanob;
 int      IsGenerated[20]; // initial/final state generation flag
 int      Nc[20];
 double   Qf[20];          // electric charge
 double   T3f[20];         // izospin, L-hand component
 double   fmcon[20];       // constituent mass
 double   fmass[20];       // current fermion mass
 // for DEBUG
 int      KeyDebug;       // =0,1,2
 int      LevPri;         // printout level
 int      Ie1Pri;         // printout from generator start
 int      Ie2Pri;         // printout from generator end
//------------------------------------
 public:
 KKdbase();                    // explicit default constructor for streamer
 KKdbase(ofstream *OutFile);   // user constructor
 ~KKdbase();                   // explicit destructor
 public:
/////////////////////////////////////////////////////////////////////////////
//

 double sqr( const Double_t x ){ return x*x;};


 void Initialize(double []);
////////////////////////////////////////////////////////////////////////////
      ClassDef(KKdbase,1); // Data base
};
////////////////////////////////////////////////////////////////////////////
#endif
