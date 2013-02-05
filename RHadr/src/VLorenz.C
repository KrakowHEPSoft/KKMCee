//////////////////////////////////////////////////////////////////////////////
//                     CLASS   VLorenz                                      //
//                                                                          //
//     Lorenz four-vector algebra with natural usage:                       //
//     (overloaded operators: +,-,*,/,[])                                   //
//     v=v1+v2 for adding vectors                                           //
//     v=v1-v2 for subtraction                                              //
//     v2=a*v1 left multiplyication by scalar a                             //
//     v2=v1*a left multiplyication by scalar a                             //
//     a = v1*v2 dotproduct of two vectors                                  //
//     theta = v1/v2 calculate angle (rads) between                         //
//             3-vector components of two 4-vectors                         //
//                                                                          //
//                                                                          //
//                 end of description  VLorenz                              //
//////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <iostream.h>
#include <fstream.h>
#include <iomanip.h>
#include <math.h>

#include "VLorenz.h"

# define sw2 setprecision(10) << setw(18)

ClassImp(VLorenz)
//////////////////////////////////////////////////////////////////////////////
VLorenz::VLorenz(double en, double px,  double py,  double pz )
{
// Constructor with explicit components
  m_comp[0]  = en;
  m_comp[1]  = px;
  m_comp[2]  = py;
  m_comp[3]  = pz;
}
//////////////////////////////////////////////////////////////////////////////
double & VLorenz::operator[](int index)
{
// overloading [] is necessary for access to components
  return m_comp[index];
}
//////////////////////////////////////////////////////////////////////////////
void VLorenz::print()
{
// printing entire four-vector in one line
  for ( int k=0; k < 4 ; k++ )   cout << sw2 << m_comp[k];
}
//////////////////////////////////////////////////////////////////////////////
VLorenz operator+(VLorenz p1, VLorenz p2)
{
// overloading operator +, adding 2 vectors
  VLorenz psum;
  for(int i=0;i<4;i++) psum[i] = p1[i] + p2[i];
  return psum;
}
//////////////////////////////////////////////////////////////////////////////
VLorenz operator-(VLorenz p1, VLorenz p2)
{
// overloading operator -, subtracting 2 vectors
  VLorenz pdif;
  for(int i=0;i<4;i++) pdif[i] = p1[i] - p2[i];
  return pdif;
}
//////////////////////////////////////////////////////////////////////////////
double operator*(VLorenz p1, VLorenz p2)
{
// overloading operator *, dot-product of 2 vectors
  double dotprod;
  dotprod = p1[0]*p2[0] -p1[1]*p2[1] -p1[2]*p2[2] -p1[3]*p2[3];
  return dotprod;
}
//////////////////////////////////////////////////////////////////////////////
VLorenz operator*(double x, VLorenz p2)
{
// overloading operator *, left multiplication of vector by scalar
  VLorenz p;
  for(int i=0;i<4;i++) p[i] = x*p2[i];
  return p;
}
//////////////////////////////////////////////////////////////////////////////
VLorenz operator*(VLorenz p1, double x)
{
// overloading operator *, right multiplication of vector by scalar
  VLorenz p;
  for(int i=0;i<4;i++) p[i] = p1[i]*x;
  return p;
}
//////////////////////////////////////////////////////////////////////////////
double operator/(VLorenz p1, VLorenz p2)
{
// overloading operator /, calculate angle between 3-vector components
  double angle;
  angle = (p1[1]*p2[1] +p1[2]*p2[2] +p1[3]*p2[3])
     /sqrt(p1[1]*p1[1] +p1[2]*p1[2] +p1[3]*p1[3])
     /sqrt(p2[1]*p2[1] +p2[2]*p2[2] +p2[3]*p2[3]);
  if(angle>1 ) angle = 1;
  angle = acos(angle);
  return angle;
}
///////////////////////////////////////////////////////////////////////////////
//                   END OF CLASS Lorenz                                     //
///////////////////////////////////////////////////////////////////////////////
