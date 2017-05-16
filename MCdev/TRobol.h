#ifndef TROBOL_H
#define TROBOL_H
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                CLASS ROBOL                                                //
//                                                                           //
//     Three-stroke engine for analysis: Initialize-Production-Finalize      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
// C++ headers
using namespace std;
#include <iostream>
#include <fstream>
#include <iomanip>

// ROOT headers
#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
//#include "TNtuple.h"
// OUR headers
#include "TRandom.h"
#include "TMCgen.h"

class TRobol : public TObject {
/// Base class part
 public:
  char     f_Name[64];          // Name of a give instance of the class
  TRandom  *f_RNgen;            //! Central RN event generator
  TMCgen   *f_MCgen;            //! MC event generator
  TFile    *f_HstFile;          //! ROOT DiskFile with all histos
  TFile    *f_GenFile;          //! ROOT DiskFile with MC generator
  ofstream *f_Out;              //! Central Logfile for messages
  ofstream  f_TraceFile;        //! Special DiskFile for debug
  double    f_NevGen;           //  event serial number
  double    f_count1;           //  auxiliary event counter (debug)
  long      f_isNewRun;         //  "is this new run?"
//////////////////////////////////////////////////////////////
 public:
  TRobol();
  TRobol(const char* );
  ~TRobol();
private:
  TRobol(const TRobol &org) { /// dummy copy constructor
  }
 public:
/// Methods
  virtual void Initialize(ofstream*, TFile*, TFile*);
  virtual void Production(double&);
  virtual void FileDump();
  virtual void Finalize();
  TH1D *TH1D_UP(const char*, const char*, int, double, double);
  TH2D *TH2D_UP(const char*, const char*, int, double, double, int, double, double);
  double sqr( const double x ){ return x*x;};
/// for debug
  void StopM(char* message){
    cout <<"++++ TRobol: "<< message << endl; exit(5);}    //Error message
////////////////////////////////////////////////////////////////////////////
  ClassDef(TRobol,1); // MC Analysis module
};
/////////////////////////////////////////////////////////////////////////////
//                End of the class TRobol                                 //
/////////////////////////////////////////////////////////////////////////////
#endif
