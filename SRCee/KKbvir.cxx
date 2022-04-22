///////////////////////////////////////////////////////////////////////////////
#include "KKbvir.h"

ClassImp(KKbvir);


KKbvir::KKbvir()
{
  // This constructor is for ROOT streamers ONLY
  // all pointers has to be NULLed
  cout<< "----> KKbvir Default Constructor (for ROOT only) "<<endl;
  m_Out= NULL;
}

///_____________________________________________________________
KKbvir::KKbvir(ofstream *OutFile)
{
  cout<< "----> KKbvir USER Constructor "<<endl;
  m_Out = OutFile;
}//KKbvir

///______________________________________________________________________________________
KKbvir::~KKbvir()
{
  //Explicit destructor
  cout<< "----> KKbvir::KKbvir !!!! DESTRUCTOR !!!! "<<endl;
}///destructor

double KKbvir::sqr( const Double_t x ){ return x*x;};

///______________________________________________________________________________________
void KKbvir::Initialize()
{
  cout  << "----> KKbvir::Initialize, Entering "<<endl;
//=================================================================
// BX*** macros are in MCdev/BXFORMAT.h
  BXOPE(*m_Out);
  BXTXT(*m_Out,"========================================");
  BXTXT(*m_Out,"======    KKbvir::Initialize      ======");
  BXTXT(*m_Out,"========================================");
  ///////////////////////////////////////////////////
}// Initialize



dcmplx KKbvir::IntReson(double MasPhot, double MassZ, double GammZ, double  s, double t, double u){
/////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                         //
//   Resonant part of virtual formfactor                                                   //
//   Needed for Greco-Pancheri-Srivastava exponentiation                                   //
//                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////
dcmplx Eps  = dcmplx(-1.0,0.0);
dcmplx MZ2  = dcmplx(MassZ*MassZ, -MassZ*GammZ);
//dcmplx MG2  = dcmplx(MasPhot*MasPhot, 0.0);
dcmplx SS   = dcmplx(s, 0.0);
dcmplx TT   = dcmplx(t, 0.0);
dcmplx UU   = dcmplx(u, 0.0);
dcmplx IntReson =
      -2.0*log(TT/UU) *log((MZ2-SS)/MZ2 );
//    -2.0*CDLN( (TT/UU) ,Eps) *CDLN( ((MZ2-SS)/MZ2 )   ,Eps);
//     $    +CDLN( (TT/UU) ,Eps) *CDLN( MG2/CDSQRT(TT*UU) ,Eps)
//     $     +Spence( ((MZ2+UU)/MZ2) ,Eps)
//     $     -Spence( ((MZ2+TT)/MZ2) ,Eps)
return IntReson;
}// IntReson

dcmplx KKbvir::IntIR(double MasPhot, double s, double t, double u){
/////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                         //
//   Virtual 2*(B(t)-B(u)) Intereference IR part to be subtracted from boxes               //
//                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////
//dcmplx Eps = DCMPLX(-1.0,0.0);
dcmplx MG2 = dcmplx(MasPhot*MasPhot, 0.0);
dcmplx TT  = dcmplx(t, 0.0);
dcmplx UU  = dcmplx(u, 0.0);
dcmplx IntIR = log(TT/UU) *log(MG2/sqrt(TT*UU))
              +dcmplx(0.5)*log(TT/UU);
return IntIR;

//IntIR =
//$      CDLN( (TT/UU) ,Eps) *CDLN( (MG2/CDSQRT(TT*UU)) ,Eps)
//$     +DCMPLX(0.5d0)*CDLN( (TT/UU) ,Eps)
}//

dcmplx KKbvir::CDLN(dcmplx X, dcmplx A){
// Apparently does not work correctly!!!!
/////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                         //
//  Complex functions that take account of the I*Epsilon prescription                      //
//  Complex logarithm of X+I*REAL(A) where a is an infinitesimal                           //
//  Programmed probably by R. Stuart                                                       //
//                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////
//--------
dcmplx PI= dcmplx(3.141592653589793238462643e0,0.0);
dcmplx CDLN;

dcmplx signA= dcmplx(1.0,0.0);
if( real(A) < 0 ) signA= dcmplx(-1.0,0.0);

if( imag(X) == 0.0 && real(X) <= 0.0) {
   CDLN = log(-X) +dcmplx(0.0,1.0)*PI *signA;
} else {
   CDLN = log(X);
}
//if( imag(CDLN) > real(PI) ) CDLN =CDLN -dcmplx(0.0,1.0)*PI;
//if( imag(CDLN) > real(-PI)) CDLN =CDLN +dcmplx(0.0,1.0)*PI;
return CDLN;
}// CDLN



dcmplx KKbvir::Spence(dcmplx Y, dcmplx E){
// argument E disabled, CDLN() replaced with log()
/////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                         //
//  Spence function of y+i*REAL(E) where E is an infinitesimal                             //
//  Programmed (in f77) probably by R. Stuart                                              //
//                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////
//DOUBLE COMPLEX  Y,E
//DOUBLE PRECISION  B(9),FACT
//DOUBLE COMPLEX  A,CLN,PISQ6,PROD,TERM,X,Z,ZSQ
//DOUBLE COMPLEX  CDLN
//INTEGER    j,i1,i2
//---------------------------------------------------------------------------------------------
double B[10];
B[1]=1.0/6.0;
B[2]=-1.0/30.0;
B[3]=1.0/42.0;
B[4]=B[2];
B[5]=5.0/66.0;
B[6]=-691.0/2730.0;
B[7]=7.0/6.0;
B[8]=-3617.0/510.0;
B[9]=43867.0/798.0;
dcmplx PISQ6= dcmplx(1.64493406684822640,0.0);
int I1=0;
int I2=0;
dcmplx X =Y;
dcmplx A =E;
dcmplx Spence;
if(X == dcmplx(0.0,0.0)) {
  Spence=(0.0,0.0);
  return Spence;
}
if( X == dcmplx(1.0,0.0)) {
  Spence=PISQ6;
  return Spence;
}
//  IF X LIES OUTSIDE THE UNIT CIRCLE THEN EVALUATE Spence(1/X)
if( abs(X) > 1.0){
  X=1.0/X;
  A=-A;
  I1=1;
}
//  IF REAL(X)>1/2 THEN EVALUATE Spence(1-X)
if( real(X) > 0.5) {
  X=1.0-X;
  A=-A;
  I2=1;
}
//  EVALUATE SERIES FOR Spence(X)
// Z=-CDLN(1.D0-X,-A)
dcmplx Z=-log(1.0-X);
dcmplx ZSQ=Z*Z;
Spence=Z-ZSQ/4.0;
dcmplx PROD=Z;
double FACT=1.0;
dcmplx TERM;
for(int J=2; J<= 18; J=J+2) {
  FACT=FACT*double((J+1)*J);
  PROD=PROD*ZSQ;
  TERM=B[J/2]/FACT*PROD;
  Spence=Spence+TERM;
  if( abs(TERM/Spence) <  1.0e-20) goto e20;
}//for
// ADD APPROPRIATE LOGS TO OBTAIN SPENCE FUNCTION OF ORIGINAL ARGUEMENT
e20:
  if( I2 == 1) {
//  Spence=-Spence+PISQ6-CDLN(X,A)*CDLN(1.D0-X,-A);
  Spence=-Spence+PISQ6-log(X)*log(1.0-X);
  X=1.0-X;
  A=-A;
  }
if( I1  == 1) {
//  CLN=CDLN(-X,-A)
  dcmplx CLN=log(-X);
  Spence=-Spence-PISQ6-CLN*CLN/2.0;
}
return Spence;
}


dcmplx KKbvir::CBoxGG(double MasPhot, double s, double t, double u){
/////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                         //
//   Box Gamma-Gamma, taken from   KORAZ/KORALB programs                                   //
//                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////
//---------
dcmplx Eps = dcmplx(-1.0,0.0);
dcmplx MG  = dcmplx(MasPhot, 0.0);
dcmplx SS  = dcmplx(s, 0.0);
dcmplx TT  = dcmplx(t, 0.0);
dcmplx UU  = dcmplx(u, 0.0);
dcmplx CiPi =dcmplx(0.0, M_PI);
//ddcmplx Clnt = CDLN( (-TT/SS) ,Eps);
dcmplx Clnt = log( (-TT/SS) );
dcmplx CBoxGG=
//     CDLN( (TT/UU) ,Eps) *(CDLN( MG**2/SS ,Eps) +CiPi) !!!   <-- Infrared part
//     +DCMPLX(0.5d0)*SS*(UU-TT)/UU**2 *( DCMPLX(0.5d0)*Clnt**2 +CiPi*Clnt)
//     -DCMPLX(0.5d0)*SS/UU*( Clnt +Cipi);
     log( (TT/UU) ) *(log( MG*MG/SS ) +CiPi)                //!!!   <-- Infrared part
     +dcmplx(0.5)*SS*(UU-TT)/UU/UU *( dcmplx(0.5)*Clnt*Clnt +CiPi*Clnt)
     -dcmplx(0.5)*SS/UU*( Clnt +CiPi);
return CBoxGG;
}// CBoxGG

dcmplx KKbvir::CBoxGZ(double MasPhot, double MassZ, double GammZ, double s, double t, double u){
/////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                         //
// Box Gamma-Z, From W. Brown, R. Decker, E. Pashos, Phys. Rev. Lett., 52 (1984), 1192     //
// Programmed similarly as in KORALZ                                                       //
//                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////
//
dcmplx Eps  = dcmplx(-1.0,0.0);
dcmplx MZ2  = dcmplx(MassZ*MassZ, -MassZ*GammZ);
dcmplx MG2  = dcmplx(MasPhot*MasPhot, 0.0);
dcmplx SS   = dcmplx(s, 0.0);
dcmplx TT   = dcmplx(t, 0.0);
dcmplx UU   = dcmplx(u, 0.0);
//CBoxGZ =
//$        CDLN( (TT/UU) ,Eps) *CDLN( MG2/CDSQRT(TT*UU) ,Eps) !!!<-- Infrared part
//$     -2*CDLN( (TT/UU) ,Eps) *CDLN( ((MZ2-SS)/MZ2 )   ,Eps)
//$     +Spence( ((MZ2+UU)/MZ2) ,Eps)
//$     -Spence( ((MZ2+TT)/MZ2) ,Eps)
//$     +(MZ2-SS)*(UU-TT-MZ2)/UU/UU *(
//$              CDLN( (-TT/SS) ,Eps) *CDLN( ((MZ2-SS)/MZ2) ,Eps)
//$             +Spence( ((MZ2+TT)/MZ2) ,Eps)
//$             -Spence( ((MZ2-SS)/MZ2) ,Eps)  )
//$     +(MZ2-SS)*(MZ2-SS)/UU/SS *CDLN( ((MZ2-SS)/MZ2) ,Eps)
//$     +(MZ2-SS)/UU             *CDLN( (-TT/MZ2)      ,Eps)
dcmplx CBoxGZ =
        log( (TT/UU) ) *log( MG2/sqrt(TT*UU) )           //!!!<-- Infrared part
     -2.0*log( (TT/UU) ) *log( ((MZ2-SS)/MZ2 ))
     +Spence( ((MZ2+UU)/MZ2) ,Eps)
     -Spence( ((MZ2+TT)/MZ2) ,Eps)
     +(MZ2-SS)*(UU-TT-MZ2)/UU/UU *(
              log( (-TT/SS) ) *log( ((MZ2-SS)/MZ2) )
             +Spence( ((MZ2+TT)/MZ2) ,Eps)
             -Spence( ((MZ2-SS)/MZ2) ,Eps)  )
     +(MZ2-SS)*(MZ2-SS)/UU/SS *log( ((MZ2-SS)/MZ2) )
     +(MZ2-SS)/UU             *log( (-TT/MZ2) );
return CBoxGZ;
}


double KKbvir::SBvirt(double alfpic, double p1p2, double m1, double m2, double MasPhot){
/////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                         //
//   Real part of B-VIRTUAL S-CHANNEL, Exact.                                              //
//   Present version according to eq.(12) in uthep-95-0801                                 //
//   Notation: mu->Nu, mu*(1+rho)->nu+xlam, rho=xlam/nu                                    //
//                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////
if( p1p2 <= m1*m2 ) {
	cout<<"##### STOP in SBvirt: p1p2,m1,m2 = "<< p1p2<<"  "<<m1<<"  "<<m2<<endl; exit(93);
}
double Nu = p1p2;
double s = 2.0*p1p2 +m1*m1 +m2*m2;
double xlam = sqrt( Nu*Nu -sqr(m1*m2) );
double SBvirt = alfpic*(
     (Nu/xlam *log((Nu+xlam)/m1/m2) -1.0) *log(sqr(MasPhot)/m1/m2)
     +xlam/s *log((Nu+xlam)/m1/m2)
     +(m1*m1-m2*m2)/(2.0*s) *log(m1/m2)
//###################################################[[[[[[[[[[[[[[[[[[[[[[[[[
//c     $     +Nu/xlam* m_pi**2                    !!!<---  Important pi**2/beta of Schwinger
//$     +m_pi**2                    !!!<--- temporary solution for the coulomb problem
     +M_PI*M_PI
//###################################################]]]]]]]]]]]]]]]]]]]]]]]]]
     +Nu/xlam*(
         -0.5 *log((Nu+xlam)/m1/m1) *log((Nu+xlam)/m2/m2)
         -0.5 *sqr(log( (m1*m1 +(Nu+xlam))/(m2*m2 + Nu+xlam) ))
         -Dilog( 2.0*xlam/(m1*m1 +Nu+xlam) )
         -Dilog( 2.0*xlam/(m2*m2 +Nu+xlam) )
     -1.0)
     );
return SBvirt;
}// SBvirt!!!


double KKbvir::TBvirt(double alfpic, double p1p2, double m1, double m2, double MasPhot){
//DOUBLE PRECISION   FUNCTION TBvirt(alfpic,p1p2,m1,m2,MasPhot)
/////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                         //
//   Real part of B-VIRTUAL T-CHANNEL  m1 is assumed to be very small.                     //
//   Present version according to eq.(14) in uthep-95-0801 with COSMETIC changes.          //
//                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////
double t  = m1*m1 +m2*m2 -2.0*p1p2;
double ta = abs(t);
double zeta = 1+ m2*m2/ta;
double TBvirt = alfpic*(
     (log(2.0*p1p2/(m1*m2)) -1.0)*log(sqr(MasPhot)/(m1*m2))
     +0.5*zeta*log(ta*zeta/(m1*m2))
     -0.5*log(ta/m1/m1)*log(ta/m2/m2)
     +Dilog(1/zeta) -1.0
     +0.5*(zeta -1.0)*log(m1/m2)
     -log(zeta)*(log(ta/(m1*m2)) +0.5*log(zeta))
     );
return TBvirt;
}// TBvirt !!!

double KKbvir::Btilda(double alfpi, double p1p2, double E1, double E2,
		              double Mas1, double Mas2, double Kmax, double MasPhot){
/////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                         //
//  Full/Complete/exact 2-particle Btilde function.                                        //
//  Exact treatment of masses, numericaly stable in high energy limit.                     //
//  Equivalent of routine Btilde programmed by W. Placzek (S.J.)                           //
//                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////
double Btilda = alfpi*(
     (p1p2*A( p1p2,Mas1,Mas2) -1 )*2*log(2*Kmax/MasPhot)
     +p1p2*A4(p1p2,E1,E2,Mas1,Mas2)
     -0.5*Mas1*Mas1*A4sng(E1,Mas1)
     -0.5*Mas2*Mas2*A4sng(E2,Mas2)
     );
return Btilda;
}//Btilda

double KKbvir::A(double p1p2, double Mas1, double Mas2){
/////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                         //
//   Function A(p1,p2) real version appropriate for B-tilde calculation                    //
//                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////
double Mas12 = Mas1*Mas2;
if ( (p1p2-Mas12) < 1e-10) {
   cout<<"+++++++ A:: WARNING, p1p2 = "<<p1p2<<endl;
   return 0.0;
}
double xlam = sqrt( (p1p2 - Mas12)*(p1p2 + Mas12) );
double A  = 1/xlam *log( (p1p2 + xlam)/Mas12 );
return A;
}//A


double KKbvir::A4sng(double E1, double Mas1){
/////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                         //
// Function (p1*p1)*A4(p1,p1) equal argument momenta.                                      //
//                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////
double bet1 = sqrt(1-Mas1*Mas1/E1/E1);
double b1ln = 2*log( (1+bet1)*E1/Mas1 );
double A4sng = -1/Mas1/Mas1/bet1 *b1ln;
return A4sng;
}//A4sng

double KKbvir::A4(double p1p2, double En1, double En2, double xm1, double xm2){
/////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                         //
//! This function provides an analytical result for the integral         !                 //
//! A4(p1,p2) being a part of the YFS IR function B-tilde.               !                 //
//! Note: This is a general case without any approximation!              !                 //
//! INPUT: p1p2    - scalar product of the 4-momenta p1 and p2;          !                 //
//!        E1,E2   - particles energies;                                 !                 //
//!        xm1,xm2 - particles masses;                                   !                 //
//!----------------------------------------------------------------------!                 //
//! Written by:  Wieslaw Placzek                Knoxville, January  1996 !                 //
//! Last update: 30.01.1996                by: W.P.                      !                 //
//                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////
// Some auxiliary variables
double E1 = En1;
double E2 = En2;
double Mas1 = xm1;
double Mas2 = xm2;
double p1s = E1*E1 - Mas1*Mas1;
double p2s = E2*E2 - Mas2*Mas2;
if (p1s < p2s) {
  Mas1 = xm2;
  Mas2 = xm1;
  E1 = En2;
  E2 = En1;
}
double Ep  = E1 + E2;
double Em  = E1 - E2;
double sm  = Mas1 + Mas2;
double dm  = Mas1 - Mas2;
double Q2  = 2*p1p2 - Mas1*Mas1 - Mas2*Mas2;
double xl  = sqrt( (Q2 + sm*sm)*(Q2 + dm*dm) );
double xq  = sqrt(Q2 + Em*Em);
double qp = xq + Em;
double qm = xq - Em;
double et0 = sqrt(E2*E2 - Mas2*Mas2);
if (p1p2 > E1*E2) et0 = -et0;
double et1 = sqrt(E1*E1 - Mas1*Mas1) + xq;
double y1  = 0.5*( (xq - Ep) + (sm*dm + xl)/qp );
double y2  = y1 - xl/qp;
double y3  = 0.5*( (xq + Ep) + (sm*dm + xl)/qm );
double y4  = y3 - xl/qm;
// Some auxiliary functions
double Eln;
if (abs(Em) > 1e-10) {
  Eln = log(abs(qm/qp))*( etaln(y1,y4,y2,y3,et1)
                        - etaln(y1,y4,y2,y3,et0) );
} else {
  Eln = 0;
}
double Vet0, Vet1;
Vet0 = Yijeta(y1,y4,et0) + Yijeta(y2,y1,et0)
     + Yijeta(y3,y2,et0) - Yijeta(y3,y4,et0)
     + 0.5*etaln(y1,y2,y3,y4,et0)*etaln(y2,y3,y1,y4,et0);
Vet1 = Yijeta(y1,y4,et1) + Yijeta(y2,y1,et1)
     + Yijeta(y3,y2,et1) - Yijeta(y3,y4,et1)
     + 0.5*etaln(y1,y2,y3,y4,et1)*etaln(y2,y3,y1,y4,et1);
// Function A4(p1,p2)
double A4 = 1/xl*(Eln + Vet1 - Vet0 );
return A4;
}//A4

double KKbvir::Yijeta(double yi, double yj, double eta){
/////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                         //
// !----------------------------------------------------------------------!                //
// ! Some auxiliary function (combination of Logs and Dilogs) used in     !                //
// ! the function A4anal for A4(p1,p2).                                   !                //
// !----------------------------------------------------------------------!                //
// ! Written by:  Wieslaw Placzek                Knoxville, January  1996 !                //
// ! Last update: 30.01.1996                by: W.P.                      !                //
//                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////
//-----
double Yijeta = 2*Dilog( ( yj-yi)/(eta-yi) )
          + 0.5*sqr(log(abs( (eta-yi)/(eta-yj) )));
return Yijeta;
}//Yijeta


double KKbvir::Btildc(double alfpi, double p1p2, double E1, double E2, double Mas1, double Mas2, double Kmax, double MasPhot){
/////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                         //
// Crude/Truncated (crude MC) 2-particle Btilde equivalent S.J.                            //
//                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////
//IMPLICIT NONE
//DOUBLE PRECISION    alfpi,p1p2,E1,E2,Mas1,Mas2,Kmax,MasPhot
//DOUBLE PRECISION    A, A4
//*-------------
double Btildc = alfpi*(
      p1p2*A( p1p2      ,Mas1,Mas2)*2*log(2*Kmax/MasPhot)
     +p1p2*A4(p1p2,E1,E2,Mas1,Mas2)
     );
return Btildc;
}// Btildc


double KKbvir::Dilog(double x){
/////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                         //
// dilogarithm FUNCTION: dilog(x)=int( -ln(1-z)/z ) , 0 < z < x .                          //
// this is the cernlib version.                                                            //
//                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////
double  a,b,y,s,t,z, dilog;
z=-1.644934066848226e0;
if(x  <  -1.0) goto e1;
if(x  <=  0.5) goto e2;
if(x  ==  1.0) goto e3;
if(x  <=  2.0) goto e4;
z=3.289868133696453e0;
e1:
  t=1.e0/x;
  s=-0.5e0;
  z=z-0.5e0*sqr(log(abs(x)));
goto e5;
e2:
  t=x;
  s=0.5e0;
  z=0.e0;
goto e5;
e3:
 dilog = 1.644934066848226e0;
 return dilog;
e4:
  t=1.e0-x;
  s=-0.5e0;
  z=1.644934066848226e0-log(x)*log(abs(t));
e5:
y=2.666666666666667e0*t+0.666666666666667e0;
b=      0.000000000000001e0;
a=y*b  +0.000000000000004e0;
b=y*a-b+0.000000000000011e0;
a=y*b-a+0.000000000000037e0;
b=y*a-b+0.000000000000121e0;
a=y*b-a+0.000000000000398e0;
b=y*a-b+0.000000000001312e0;
a=y*b-a+0.000000000004342e0;
b=y*a-b+0.000000000014437e0;
a=y*b-a+0.000000000048274e0;
b=y*a-b+0.000000000162421e0;
a=y*b-a+0.000000000550291e0;
b=y*a-b+0.000000001879117e0;
a=y*b-a+0.000000006474338e0;
b=y*a-b+0.000000022536705e0;
a=y*b-a+0.000000079387055e0;
b=y*a-b+0.000000283575385e0;
a=y*b-a+0.000001029904264e0;
b=y*a-b+0.000003816329463e0;
a=y*b-a+0.000014496300557e0;
b=y*a-b+0.000056817822718e0;
a=y*b-a+0.000232002196094e0;
b=y*a-b+0.001001627496164e0;
a=y*b-a+0.004686361959447e0;
b=y*a-b+0.024879322924228e0;
a=y*b-a+0.166073032927855e0;
a=y*a-b+1.935064300869969e0;
dilog=s*t*(a-b)+z;
return dilog;
}//Dilog
