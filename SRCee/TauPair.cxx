///////////////////////////////////////////////////////////////////////////////
#include "TauPair.h"

ClassImp(TauPair);


TauPair::TauPair()
{
  // This constructor is for ROOT streamers ONLY
  // all pointers has to be NULLed
  cout<< "----> TauPair Default Constructor (for ROOT only) "<<endl;
  m_Out= NULL;
}

///_____________________________________________________________
TauPair::TauPair(ofstream *OutFile)
{
  cout<< "----> TauPair USER Constructor "<<endl;
  m_Out = OutFile;
}//TauPair

///______________________________________________________________________________________
TauPair::~TauPair()
{
  //Explicit destructor
  cout<< "----> TauPair::TauPair !!!! DESTRUCTOR !!!! "<<endl;
}///destructor

double TauPair::sqr( const Double_t x ){ return x*x;};

///______________________________________________________________________________________
void TauPair::Initialize(double ypar[])
{
  cout  << "----> TauPair::Initialize, Entering "<<endl;
//=================================================================
// BX*** macros are in MCdev/BXFORMAT.h
  BXOPE(*m_Out);
  BXTXT(*m_Out,"========================================");
  BXTXT(*m_Out,"======    TauPair::Initialize      ======");
  BXTXT(*m_Out,"========================================");

  taupair_initialize_(ypar);
  ///////////////////////////////////////////////////
}// Initialize
