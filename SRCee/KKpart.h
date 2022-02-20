
#ifndef KKpart_H
#define KKpart_H

#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

#include "BXFORMAT.h"

#include "TObject.h"
#include "TLorentzVector.h"

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Auxiliary class KKpart of 4-vector, mass etc. for spinor algebra         //
//  used only in CEEX matrix element                                         //
//  (overloaded operators: +,-,*,[]                                          //
//  Used also for photons in case photon momentum is used in spinorproduct   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

class KKpart : public TObject {
  // constructor
  public:
	double  M;                    // Mass
	int     Hel;                  // Helicity +-1 both for fermion and photon
	int     C;                    // =+1 U spinor or =-1 for V spinor
        double  P[4];                 // four-momentum
       
        int pdf_id;                   // Staszek we need to implement this 
        double prod_vertex[3];        // Staszek we need to implement this

  

public:
    KKpart();                        // Constructor for Foam
    ~KKpart();                       // Destructor  for Foam
    KKpart(const KKpart &Part);      // copy constructor, is it needed ???
    KKpart(const double );           // USER Constructor
    KKpart(const double x0, const double x1, const double x2,const double x3); // USER Constructor
 //////////////////////////////////////////////////////////////////////////////////////////////
//                         Overloading operators                                             //
///////////////////////////////////////////////////////////////////////////////////////////////
    KKpart& operator =(const KKpart&);   // = operator; Substitution (const ?)
    double& operator[](int);             // [] provides POINTER to 4-momentum component
//////////////////////////   OTHER METHODS  //////////////////////////////////
    KKpart& operator+=(const  KKpart&);       // +=; add 4-vector u+=v      (FAST)
    KKpart& operator-=(const  KKpart&);       // +=; add 4-vector u+=v      (FAST)
    KKpart& operator*=(const  double&);       // *=; mult. by scalar v*=x   (FAST)

    double operator*(const KKpart &p2){ // operator * ; dot product of 2 vectors;
           return (*this).P[0]*p2.P[0] -(*this).P[1]*p2.P[1] -(*this).P[2]*p2.P[2]-(*this).P[3]*p2.P[3];}

    KKpart  operator+( const  KKpart &p){     // u=v+s, C imherited from v, might be SLOW!!!
        KKpart res = *this; res +=p; return res;
    }
    KKpart  operator-( const  KKpart &p){     // u=v-s, C imherited from v, might be SLOW!!!
        KKpart res = *this; res -=p; return res;
    }
    double M2(){ KKpart Q =(*this) ; return Q*Q;}; // square of fourvector component

    void SetMom(const TLorentzVector &);  // USER constructor
    void SetAll(const int C0, const int Hel0, const double M0, TLorentzVector &P0);
    void Print(void);                  // Prints vector
/////////////////////////////////////////////////////////////////////////////
    ClassDef(KKpart,1) //n-dimensional vector with dynamical allocation
};
////////////////////////////////////////////////////////////////////////////
#endif
