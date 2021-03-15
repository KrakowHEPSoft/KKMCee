///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//      Auxiliary class KKpart of 4-vector, mass etc. for Spinor algebra     //
//      used only in CEEX matrix element                                     //
//      (overloaded operators: +,-,*,[]                                      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
#include "KKpart.h"

ClassImp(KKpart);


#define SW20 setw(20)<<setprecision(14)

KKpart::KKpart(){
// Default constructor creates "empty vector of zero dimension"
M    =0;
Hel  =0;
C    =0;
P[0] =0;
P[1] =0;
P[2] =0;
P[3] =0;
}//KKpart

///////////////////////////////////////////////////////////////////////////////
void KKpart::SetMom(const TLorentzVector &P0){
// truncated USER constructor
M    =0;
Hel  =0;
C    =0;
P[0] = P0.E();
P[1] = P0.Px();
P[2] = P0.Py();
P[3] = P0.Pz();
}
//KKpart

///////////////////////////////////////////////////////////////////////////////
void KKpart::SetAll( const int C0, const int Hel0,const double M0, TLorentzVector &P0){
// complete USER constructor
M    =M0;
Hel  =Hel0;
C    =C0;
P[0] = P0.E();
P[1] = P0.Px();
P[2] = P0.Py();
P[3] = P0.Pz();
}
//KKpart

///////////////////////////////////////////////////////////////////////////////
KKpart::KKpart(const double x){
// exmple USER constructor
M    =0; Hel  =0; C    =0;
P[0] =x; P[1] =x; P[2] =x; P[3] =x;
}
//KKpart

///////////////////////////////////////////////////////////////////////////////
KKpart::KKpart(const double x0, const double x1, const double x2,const double x3){
// exmple USER constructor
M    =0; Hel  =0; C    =0;
P[0] =x0; P[1] =x1; P[2] =x2; P[3] =x3;
}
//KKpart

//////////////////////////////////////////////////////////////////////////////
KKpart::KKpart(const KKpart &Part){
// "copy constructor" not to be used explicitly???
M    = Part.M;
Hel  = Part.Hel;
C    = Part.C;
for(int i; i<4; i++) P[i]= Part.P[i];
//cout << " Kpart::KKpart copy constructor ???"<<endl;
}//KKpart

//////////////////////////////////////////////////////////////////////////////
KKpart::~KKpart(){ }
//
//////////////////////////////////////////////////////////////////////////////
//                     Overloading operators                                //
//////////////////////////////////////////////////////////////////////////////
KKpart& KKpart::operator =(const KKpart& Part){
// operator = ;  substitution
int i;
if (&Part == this) return *this;
M    = Part.M;
Hel  = Part.Hel;
C    = Part.C;
for(int i; i<4; i++) P[i]= Part.P[i];
return *this;
}//

//============================================================================
double &KKpart::operator[](int i){
// [] is for acces to elements for substitution
// one should use rather use a=b than explicit loop!
if ((i<0) || (i>=4)){
  cout<<"++++>> KKpart::operator[], out of range"<<endl;
  exit(1);
}
return P[i];
}
//============================================================================
KKpart& KKpart::operator*=(const double &x){
// operator *=; multiply vector part by scalar; c*=x,
for(int i=0;i<4;i++)
  P[i] = P[i]*x;
return *this;
}
//============================================================================
KKpart& KKpart::operator+=(const KKpart& Shift){
// operator +=; add vector; c*=x, tested
for(int i=0;i<4;i++)
  P[i] = P[i]+Shift.P[i];
return *this;
}
//============================================================================
KKpart& KKpart::operator-=(const KKpart& Shift){
// operator -=; subtract vector; c*=x, tested
for(int i=0;i<4;i++)
  P[i] = P[i]-Shift.P[i];
return *this;
}

//============================================================================
void KKpart::Print(){
//////////////////////////////////////////////////////////////
// printing entire four-vector in one line (withot endline)
  for ( int k=1; k < 4 ; k++ ) cout << SW20 << P[k];
  cout << SW20 << P[0];  cout << SW20 << M; cout << SW20 << Hel;
  cout<<endl;
}//MomPrint
//============================================================================
