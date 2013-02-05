///////////////////////////////////////////////////////////////////////////////
//    make mainVLorenz                                                       //
//    make mainVLorenz.exe                                                   //
///////////////////////////////////////////////////////////////////////////////
#include <stdlib.h>
#include <iostream.h>
#include <fstream.h>
#include <math.h>
#include <iomanip.h>

// ROOT headers
#include "TROOT.h"

#include "VLorenz.h"


/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
//class VLorenz : public TObject {
//class VLorenz {
//public:
//                double m_comp[4];       // 4-momentum, comp[0] is energy
//public:
//                VLorenz()  {;}
//                VLorenz(double en, double px,  double py,  double pz );
/////////////////////   virtual      ~VLorenz() { printf("VLorenz%x\n",this);}
//   virtual      ~VLorenz() {;}
//   void         VLorenz::print();
//public:
//                double & operator[](int index);
//  //ClassDef(VLorenz,1)   // VLorenz  class
//};

//////////////////////////////////////////////////////////////////////////////
//                     Overloading operators                                //
//////////////////////////////////////////////////////////////////////////////
//VLorenz operator+(VLorenz p1, VLorenz p2);
//VLorenz operator-(VLorenz p1, VLorenz p2);
//double  operator*(VLorenz p1, VLorenz p2);
//VLorenz operator*(double x, VLorenz p2);
//VLorenz operator*(VLorenz p1, double x);
//double  operator/(VLorenz p1, VLorenz p2);
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
//VLorenz::VLorenz(double en, double px,  double py,  double pz )
//{
// Constructor with explicit components
//  m_comp[0]  = en;
//  m_comp[1]  = px;
//  m_comp[2]  = py;
//  m_comp[3]  = pz;
//}
//////////////////////////////////////////////////////////////////////////////
//double & VLorenz::operator[](int index)
//{
// overloading [] is necessary for access to components
//  return m_comp[index];
//}
//////////////////////////////////////////////////////////////////////////////
//void VLorenz::print()
//{
//# define sw2 setprecision(10) << setw(18)
// printing entire four-vector in one line
//  for ( int k=0; k < 4 ; k++ )   cout << sw2 << m_comp[k];
//}
//////////////////////////////////////////////////////////////////////////////
//VLorenz operator+(VLorenz p1, VLorenz p2)
//{
// overloading operator +, adding 2 vectors
//  VLorenz psum;
//  for(int i=0;i<4;i++) psum[i] = p1[i] + p2[i];
//  return psum;
//}
//////////////////////////////////////////////////////////////////////////////
//VLorenz operator-(VLorenz p1, VLorenz p2)
//{
// overloading operator -, subtracting 2 vectors
//  VLorenz pdif;
//  for(int i=0;i<4;i++) pdif[i] = p1[i] - p2[i];
//  return pdif;
//}
//////////////////////////////////////////////////////////////////////////////
//double operator*(VLorenz p1, VLorenz p2)
//{
// overloading operator *, dot-product of 2 vectors
//  double dotprod;
//  dotprod = p1[0]*p2[0] -p1[1]*p2[1] -p1[2]*p2[2] -p1[3]*p2[3];
//  return dotprod;
//}
//////////////////////////////////////////////////////////////////////////////
//VLorenz operator*(double x, VLorenz p2)
//{
// overloading operator *, left multiplication of vector by scalar
//  VLorenz p;
//  for(int i=0;i<4;i++) p[i] = x*p2[i];
//  return p;
//}
//////////////////////////////////////////////////////////////////////////////
//VLorenz operator*(VLorenz p1, double x)
//{
// overloading operator *, right multiplication of vector by scalar
//  VLorenz p;
//  for(int i=0;i<4;i++) p[i] = p1[i]*x;
//  return p;
//}
//////////////////////////////////////////////////////////////////////////////
//double operator/(VLorenz p1, VLorenz p2)
//{
// overloading operator /, calculate angle between 3-vector components
//  double angle;
//  angle = (p1[1]*p2[1] +p1[2]*p2[2] +p1[3]*p2[3])
//     /sqrt(p1[1]*p1[1] +p1[2]*p1[2] +p1[3]*p1[3])
//     /sqrt(p2[1]*p2[1] +p2[2]*p2[2] +p2[3]*p2[3]);
//  if(angle>1 ) angle = 1;
//  angle = acos(angle);
//  return angle;
//}
///////////////////////////////////////////////////////////////////////////////
//                   END OF CLASS Lorenz                                     //
///////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////


TROOT root("lorenz", "Lorenz Test", 0);

main()
{
  VLorenz  v1(12.5,  -0.12,  0.31, -13.0);
  cout << "v1    = ";  v1.print(); cout << "\n" << flush;

  VLorenz  v2(11.0,   0.34,  0.50,   8.3);
  cout << "v2    = ";  v2.print(); cout << "\n" << flush;

  VLorenz  v3 = v1 +v2; 
  cout << "v1+v2 = "; v3.print();  cout << "\n" << flush;

  VLorenz  v4 = v1 -v2;
  cout << "v1-v2 = ";  v4.print(); cout << "\n" << flush;

  VLorenz  v5 = 2*v1;
  cout << "2*v1  = ";  v5.print(); cout << "\n" << flush;

  VLorenz  v6 = v1*2;
  cout << "v1*2  = ";  v6.print(); cout << "\n" << flush;

  double prod = v1*v2;
  cout << "v1*v2 = ";  cout << prod << "\n" << flush;

  double angl = v1/v2;
  cout << "v1/v2 = ";  cout << angl << "\n" << flush;
}
