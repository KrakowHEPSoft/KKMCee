///////////////////////////////////////////////////////////////////////////////
//         Template of the class with ROOT persistency
///////////////////////////////////////////////////////////////////////////////

#ifndef KKlasa_H
#define KKlasa_H

#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

#include "BXFORMAT.h"
#include "TObject.h"
//________________________________________________________________________
class KKlasa: public TObject{
 public:
 ofstream *m_Out;     //! pointer to external Logfile for messages
// class member data
 public:
 double   CMSene;
//------------------------------------
// Obligatory members
  public:
  KKlasa();                    // explicit default constructor for streamer
  KKlasa(ofstream *OutFile);   // user constructor
  ~KKlasa();                   // explicit destructor
  public:
/////////////////////////////////////////////////////////////////////////////
// class member functions
double sqr( const double x );

void Initialize();
////////////////////////////////////////////////////////////////////////////
       ClassDef(KKlasa,1); // Data base
};// KKlasa class
////////////////////////////////////////////////////////////////////////////
#endif
