///////////////////////////////////////////////////////////////////////////////
//         Template of the class with ROOT persistency
///////////////////////////////////////////////////////////////////////////////

#ifndef KKborn_H
#define KKborn_H

#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <math.h>
using namespace std;

#include "TObject.h"
#include "KKdbase.h"
#include "KKdizet.h"

typedef complex<double> dcmplx;

#include "BXFORMAT.h"

//________________________________________________________________________
class KKborn: public TObject{
 public:
 ofstream *m_Out;     //! pointer to external Logfile for messages
 KKdbase  *DB;
 KKdizet  *DZ;
// class member data
 public:
 double   m_CMSene;
 int      m_icont;
//------------------------------------
// Obligatory members
  public:
  KKborn();                    // explicit default constructor for streamer
  KKborn(ofstream *OutFile);   // user constructor
  ~KKborn();                   // explicit destructor
  public:
/////////////////////////////////////////////////////////////////////////////
// class member functions
double sqr( const double x );

void Initialize();

void SetDB( KKdbase *DBase){ DB = DBase;};
void SetDZ( KKdizet *Dizet){ DZ = Dizet;};

double BornSimple(int KFi, int KFf, double svar, double costhe);

double Born_Dizet(int KFi, int KFf, double svar, double CosThe,double eps1, double eps2, double ta, double tb);

double Born_DizetS(int KFi, int KFf, double svar, double CosThe);

////////////////////////////////////////////////////////////////////////////
       ClassDef(KKborn,1); // Data base
};// KKborn class
////////////////////////////////////////////////////////////////////////////
#endif
