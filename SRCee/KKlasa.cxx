///////////////////////////////////////////////////////////////////////////////
#include "KKlasa.h"

ClassImp(KKlasa);


KKlasa::KKlasa()
{
  // This constructor is for ROOT streamers ONLY
  // all pointers has to be NULLed
  cout<< "----> KKlasa Default Constructor (for ROOT only) "<<endl;
  m_Out= NULL;
}

///_____________________________________________________________
KKlasa::KKlasa(ofstream *OutFile)
{
  cout<< "----> KKlasa USER Constructor "<<endl;
  m_Out = OutFile;
}//KKlasa

///______________________________________________________________________________________
KKlasa::~KKlasa()
{
  //Explicit destructor
  cout<< "----> KKlasa::KKlasa !!!! DESTRUCTOR !!!! "<<endl;
}///destructor

double KKlasa::sqr( const Double_t x ){ return x*x;};

///______________________________________________________________________________________
void KKlasa::Initialize()
{
  cout  << "----> KKlasa::Initialize, Entering "<<endl;
//=================================================================
// BX*** macros are in MCdev/BXFORMAT.h
  BXOPE(*m_Out);
  BXTXT(*m_Out,"========================================");
  BXTXT(*m_Out,"======    KKlasa::Initialize      ======");
  BXTXT(*m_Out,"========================================");
  ///////////////////////////////////////////////////
}// Initialize
