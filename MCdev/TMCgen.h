#ifndef MCGEN_H
#define MCGEN_H
///////////////////////////////////////////////////////////////////////////////
///    Base Class   TMCgen for severa types of MC generators

/// C++ headers
using namespace std;
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>

/// ROOT headers
#include "TRandom.h"
#include "TROOT.h"
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TString.h"
///
#include "TFOAM.h"
#include "TFOAM_INTEGRAND.h"
/// My headers
#include"BXFORMAT.h"

class TMCgen   : public TFOAM_INTEGRAND {
 public:
  char  f_Name[64];       // Name of a give instance of the class
  float f_Version;         // Actual VERSION of the program
  char  f_Date[40];        // Release DATE of the program
  //-------------- USER Input parameters--------------------------
 public:
  /// Engines and services
  TRandom  *f_RNgen;            //!  External RN event generator
  TFOAM    *f_FoamI;            //  Foam object for generating Initial density
  TH1D     *f_TMCgen_NORMA;     //! special histo keeping overall normalization
   /// data members
 public:
  int       f_IsInitialized;    //  prevents repeating initialization
  TFile    *f_GenFile;          //! ROOT DiskFile with MC generator object and data
  TFile    *f_HstFile;          //! ROOT DiskFile with all histos
  ofstream *f_Out;              //! External Logfile for messages
  double    f_NevGen;           //  event serial number
  TString m_MEtype;
 public:
  TMCgen();                // explicit default constructor for streamer
  TMCgen(const char*);     // user constructor
  TMCgen(const char*, const char*, double);     // user constructor
  ~TMCgen();               // explicit destructor
 public:
/// Methods
  virtual void Initialize(TRandom*, ofstream*, TH1D*);
  virtual void Redress(   TRandom*, ofstream*, TH1D*);
  virtual void Generate();
  virtual void Finalize();
  virtual double Density(int nDim, double *Xarg){return 0;};
  double sqr( const double x ){ return x*x;};
  int Min( int &i, int &j){ if(i<j) return i; else return j;};
/// for debug
  void StopM(const char* message){
    cout <<"++++ TMCgen: "<< message << endl; exit(5);}    //Error message
///  getters
  int GetIsNewRun() const
    { if(f_IsInitialized == 0 ) return 1; else return 0;}
///
  ClassDef(TMCgen,2); // Monte Carlo generator
};
/////////////////////////////////////////////////////////////////////////////
//                End of the class TMCgen                                  //
/////////////////////////////////////////////////////////////////////////////
#endif
