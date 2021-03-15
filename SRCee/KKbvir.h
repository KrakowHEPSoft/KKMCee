#ifndef KKbvir_H
#define KKbvir_H

#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>

using namespace std;

#include "KKpart.h"

#include "BXFORMAT.h"
#include "TObject.h"

typedef complex<double> dcmplx;

//________________________________________________________________________
class KKbvir: public TObject{
 public:
 ofstream *m_Out;     //! pointer to external Logfile for messages
// class member data
 public:
 double   CMSene;
//------------------------------------
// Obligatory members
  public:
  KKbvir();                    // explicit default constructor for streamer
  KKbvir(ofstream *OutFile);   // user constructor
  ~KKbvir();                   // explicit destructor
  public:
/////////////////////////////////////////////////////////////////////////////
// class member functions
double sqr( const double x );

void Initialize();
dcmplx IntReson(double MasPhot, double MassZ, double GammZ, double  s, double t, double u);
dcmplx CDLN(dcmplx X, dcmplx A);
dcmplx IntIR(double MasPhot, double s, double t, double u);
dcmplx Spence(dcmplx Y, dcmplx E);
dcmplx CBoxGG(double MasPhot, double s, double t, double u);
dcmplx CBoxGZ(double MasPhot, double MassZ, double GammZ, double s, double t, double u);
double Dilog(double x);
double TBvirt(double alfpic, double p1p2, double m1, double m2, double MasPhot);
double Btilda(double alfpi, double p1p2, double E1, double E2, double Mas1, double Mas2, double Kmax, double MasPhot);
double A(double p1p2, double Mas1, double Mas2);
double A4sng(double E1, double Mas1);
double A4(double p1p2, double En1, double En2, double xm1, double xm2);
double Yijeta(double yi, double yj, double eta);
double SBvirt(double alfpic, double p1p2, double m1, double m2, double MasPhot);
double Btildc(double alfpi, double p1p2, double E1, double E2, double Mas1, double Mas2, double Kmax, double MasPhot);

inline double etaln(double x1, double x2, double x3, double x4, double z){
    return log(abs( (z-x1)*(z-x2)/(z-x3)/(z-x4) ));
}//etaln

////////////////////////////////////////////////////////////////////////////
       ClassDef(KKbvir,1); // Data base
};// KKbvir class
////////////////////////////////////////////////////////////////////////////
#endif
