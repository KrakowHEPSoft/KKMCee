#ifndef VLorenz_H
#define VLorenz_H
//////////////////////////////////////////////////////////////////////////////
//                     CLASS   VLorenz                                      //
//                                                                          //
//     Lorenz four-vector algebra with natural usage:                       //
//     (overloaded operators: +,-,*,/,[])                                   //
//////////////////////////////////////////////////////////////////////////////

#include "TObject.h"

class VLorenz : public TObject {
//class VLorenz {
public:
                double m_comp[4];       // 4-momentum, comp[0] is energy
public:
                VLorenz()  {;}
                VLorenz(double en, double px,  double py,  double pz );
///////////////////   virtual      ~VLorenz() { printf("VLorenz%x\n",this);}
   virtual      ~VLorenz() {;}
   void         print();
public:
                double & operator[](int index);
  ClassDef(VLorenz,1)   // VLorenz  class
};

//////////////////////////////////////////////////////////////////////////////
//                     Overloading operators                                //
//////////////////////////////////////////////////////////////////////////////
VLorenz operator+(VLorenz p1, VLorenz p2);
VLorenz operator-(VLorenz p1, VLorenz p2);
double  operator*(VLorenz p1, VLorenz p2);
VLorenz operator*(double x, VLorenz p2);
VLorenz operator*(VLorenz p1, double x);
double  operator/(VLorenz p1, VLorenz p2);
//////////////////////////////////////////////////////////////////////////////
#endif
