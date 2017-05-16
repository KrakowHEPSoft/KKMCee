#include "TMCgen.h"

#define SP21 setw(21)<<setprecision(13)
#define SP15 setw(15)<<setprecision(9)
#define SP10 setw(10)<<setprecision(5)

ClassImp(TMCgen);


///______________________________________________________________________________________
TMCgen::TMCgen(){
/// explicit default constructor for ROOT streamers
/// should not be used by the user
  cout<< "====> TMCgen::TMCgen DEFAULT Constructor (for ROOT only) in MCdev"<<endl;
  f_Out     = NULL;
  f_GenFile = NULL;
  f_HstFile = NULL;
  f_RNgen   = NULL;
  f_FoamI   = NULL;
  f_TMCgen_NORMA = NULL;
}

///______________________________________________________________________________________
TMCgen::~TMCgen(){
  //!Explicit destructor
  cout<< "====> TMCgen::TMCgen default DESTRUCTOR !!!! "<<endl;
  //
}///destructor

//______________________________________________________________________________________
//TMCgen::TMCgen(const TMCgen &From){
//// Copy Constructor  NOT IMPLEMENTED (NEVER USED)
//  cout  <<"+++++ Stop in Copy Constructor TMCgen::TMCgen  NOT IMPLEMENTED"<<endl;
//  *f_Out<<"+++++ Stop in Copy Constructor TMCgen::TMCgen  NOT IMPLEMENTED"<<endl;
//  exit(1);
//  //m_nDim = From.m_nDim; // this is to apease insure
//}

///______________________________________________________________________________________
TMCgen::TMCgen(const char* Name){
//! Constructor to be used by the user!!!
//! Its sole important role is to define ALL DEFAULTS.
//! such that these defaults can be changed 
//! by the user before calling TMCgen::Initialize

  cout<< "====> TMCgen::TMCgen USER Constructor in MCdev"<<endl;
  if(strlen(Name)  >65) {
    cout<< " ++++STOP in TMCgen::TMCgen, Name too long "<<strlen(Name)<<endl;
    exit(2);
  }
  sprintf(f_Name,"%s",Name);         // Class name
  sprintf(f_Date,"%s","  Release date:  yyyy.mm.dd   "); // Release date
  f_Version  = 1.00;                                      // Release version
  ///
  f_RNgen   = NULL;
  f_FoamI   = NULL;
  f_Out     = NULL;
  f_GenFile = NULL;
  f_HstFile = NULL;
  f_TMCgen_NORMA = NULL;
  
  ///
  f_IsInitialized = 0; /// prevents Initialization when reset to 1
  f_NevGen  =   0;
  
}//! user constructor


///______________________________________________________________________________________
TMCgen::TMCgen(const char* Name, const char* Date, double Version){
//! Constructor to be used by the user!!!
//! Its sole important role is to define ALL DEFAULTS.
//! such that these defaults can be changed 
//! by the user before calling TMCgen::Initialize

  cout<< "====> TMCgen::TMCgen USER Constructor in MCdev"<<endl;
  if(strlen(Name)  >65 || strlen(Date) > 30 ) {
    cout<< " ++++STOP in TMCgen::TMCgen, Name or date too long "<<strlen(Name)<<endl; exit(2);
  }
  sprintf(f_Name,"%s",Name);         // Class name
  sprintf(f_Date,"%s",Date);         // Release date
  f_Version  = Version;                                      // Release version
  ///
  f_RNgen   = NULL;
  f_FoamI   = NULL;
  f_Out     = NULL;
  f_GenFile = NULL;
  f_HstFile = NULL;
  f_TMCgen_NORMA = NULL;
  ///
  f_IsInitialized = 0; /// prevents Initialization when reset to 1
  f_NevGen  =   0;
  
}//! user constructor



///______________________________________________________________________________________
void TMCgen::Initialize(TRandom *RNgen, ofstream *OutFile, TH1D* h_NORMA)
{
  //////////////////////////////////////////////////////////////////////////////////////
  /// Initializes MC generator
  /// Input parameters of the generator are set to the default values
  /// in the "user constructor" TMCgen(const char*)
  /// and the user has a chance to reset them before calling THIS method
  /// and after calling "user constructor". Later on,
  /// once this method is called, the user should not touch input params!
  ///         *** IMPORTANT ***
  /// Repeated of this Initialize method is allowed
  cout<< "====> TMCgen::Initialize, user initializator "<<endl;
  f_Out     = OutFile;
  f_RNgen   = RNgen;
  f_TMCgen_NORMA = h_NORMA; /// normalization histo taken from outside
  ///----------------------------------------
  /// physics params, redefined defaults
  ///----------------------------------------
  ///=================================================================
  BXOPE(*f_Out);
  BXTXT(*f_Out,"========================================");
  BXTXT(*f_Out,"======     TMCgen::Initialize     ======");
  BXTXT(*f_Out,"========================================");
  BXTXT(*f_Out,f_Name);
  BX1F(*f_Out,"  Version",f_Version,  f_Date);
  BXTXT(*f_Out,"============== INPUT ===================");
  BXCLO(*f_Out);
}//! TMCgen::Initialize


///______________________________________________________________________________________
void TMCgen::Redress(TRandom *RNgen, ofstream *OutFile, TH1D* h_NORMA)
{
  StopM("TMCgen::Redress: +++ stop, not implemented in base class");
} /// TMCgen::Redress


///______________________________________________________________________________________
void TMCgen::Finalize()
{
  ///===================================================
  ///   Finalize MC  run, final printouts, cleaning etc.
  double nevtot = f_TMCgen_NORMA->GetBinContent(2);
  ///
  BXOPE(*f_Out);
  BXTXT(*f_Out,"========================================");
  BXTXT(*f_Out,"======     TMCgen::Finalize       ======");
  BXTXT(*f_Out,"========================================");
  BX1I(*f_Out,"   nevtot", nevtot, " number of generated events                      =");
}//!Finalize


///______________________________________________________________________________________
void TMCgen::Generate()
{
 StopM("TMCgen::Generate: +++ stop, not implemented in base class");
}//! Generate

