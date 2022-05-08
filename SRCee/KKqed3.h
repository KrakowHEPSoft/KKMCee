/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
//                     Pseudo-CLASS  QED3                                          //
//                                                                                 //
//   Calculation of QED matrix element with Yennie-Frautschi-Suura exponentiation  //
//   for s-chanel exchange-exchange fermion-anfifermion production processe.       //
//   Order alpha^1 is complete, beyond O(alf^1) leading-log is mainly exploited.   //
//                                                                                 //
//   e+ e- ---> f + fbar + n gamma                                                 //
//                                                                                 //
//   The following contributions are included:                                     //
//                                                                                 //
//   ISR:  O(L^0*alf^0)                                                            //
//         O(L^1*alf^1)  O(L^0*alf^1)                                              //
//         O(L^2*alf^2)                                                            //
//         O(L^3*alf^3)                                                            //
//   FSR:  O(L^0*alf^0)                                                            //
//         O(L^1*alf^1)  O(L^0*alf^1)                                              //
//         O(L^2*alf^2)                                                            //
//                                                                                 //
//   Neglected:                                                                    //
//      ISR*FSR interferences, spin polarization                                   //
//      t-chanel exchanges                                                         //
//                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////
// HERWIRI2: Added support for different choices of KeyElw to simplify

#ifndef KKqed3_H
#define KKqed3_H

#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

#include "KKdbase.h"
#include "KKevent.h"
#include "KKborn.h"

#include "BXFORMAT.h"
#include "TObject.h"
//________________________________________________________________________
class KKqed3: public TObject{
 public:
 ofstream *m_Out;     //! pointer to external Logfile for messages
// class member data
 public:
 KKdbase  *DB;                     //! Database
 KKevent  *m_Event;                //! MC event record
 KKborn   *m_BornV;                //!
 public:
 static const int maxPhot =   101; // max. no. of KKMC photons +1
 static const int maxWT   =  1001;   // max. num. KKMC wt list +1
 double    m_CMSene;
 int       m_KeyOrd;               //
 int       m_nphox;                // ISR
 int       m_nphoy;                // FSR
 double    m_yini[maxPhot];        // ISR
 double    m_zini[maxPhot];
 double    m_yfin[maxPhot];        // FSR
 double    m_zfin[maxPhot];
 // Elements of beta calculation, Auxiliary/temporary
 double    m_beta20,m_beta21,m_beta30;
 double    m_betx12,m_bety12;
 double    m_beti12,m_beti21;
 double    m_betf01;
 // Contributions from individual photons
 double    m_xSfac[maxPhot];            // ISR soft factors
 double    m_ySfac[maxPhot];            // FSR soft factors
 double m_betx10[maxPhot],              // ISR beta1 tree     times (1+delf1), factorization//
        m_betx11[maxPhot],              // ISR beta1 one-loop times (1+delf1), factorization//
        m_bety10[maxPhot],              // FSR beta1 tree     times (1+deli1), factorization//
        m_bety11[maxPhot],              // FSR beta1 one-loop times (1+deli1), factorization//
        m_betxx20[maxPhot][maxPhot],    // beta2 ISR*ISR tree
        m_betxy20[maxPhot][maxPhot],    // beta2 ISR*FSR tree
        m_betyy20[maxPhot][maxPhot],    // beta2 FSR*FSR tree
        m_beti10[maxPhot],              // beta1 tree   ISR only
        m_beti11[maxPhot],              // beta1 1-loop ISR only
        m_beti20[maxPhot][maxPhot],     // beta2 tree   ISR only
        m_betf10[maxPhot],              // beta1 tree   FSR only, for beta2, beta3
        m_betf11[maxPhot];              // beta1 1-loop FSR only, for beta2, beta3
 double m_Beta00,                               // beta0 O(alf0) ISR+FSR
        m_Beta01,  m_xBet10,  m_yBet10,         // beta1 O(alf1) ISR+FSR
        m_Beta02,                               // beta0 O(alf2) ISR+FSR
        m_xBet11,  m_yBet11,                    // beta1 O(alf2) ISR+FSR
        m_xxBet20, m_xyBet20, m_yyBet20,        // beta2 O(alf2) ISR+FSR
        m_Beta03,                               // beta0 O(alf3) ISR+FSR
        m_xBet12,  m_yBet12,                    // beta1 O(alf3) ISR+FSR
        m_xxBet21, m_xyBet21, m_yyBet21,        // beta2 O(alf3) ISR+FSR
        m_xxxBet30,m_xxyBet30,m_xyyBet30,       // beta3 O(alf3) ISR+FSR
        m_beti00,                               // O(alf0) pure ISR
        m_beti01,m_sbti10,                      // O(alf1) pure ISR
        m_beti02,m_sbti11,m_sbti20,             // O(alf2) pure ISR
        m_beti03,m_sbti12,m_sbti21,m_sbti30;    // O(alf2) pure ISR
 double m_WtSet[maxWT];
//***********************************************************************
 double     m_vlim1;
 double     m_vlim2;
 int        m_icont;              // event counter for debug
// Obligatory members
  public:
  KKqed3();                    // explicit default constructor for streamer
  KKqed3(ofstream *OutFile);   // user constructor
  ~KKqed3();                   // explicit destructor
  public:
/////////////////////////////////////////////////////////////////////////////
// class member functions
//    INLINE functions
//double sqr( const double x );
inline double sqr( double x ){ return x*x;};
inline double chi1( double a  ){ return 0.5*(1-a)*(1-a); };
inline double chi2( double a, double b){ return 0.5*((1-a)*(1-a)+(1-b)*(1-b));};
//
inline double wm0( double del, double a, double b){ return 1 -2*del -del*(a/b+b/a); };
// O(alf1) amplitude, bremsstrahlung factor, del->0 limit
inline double wm1(double del,double a,double b){
	return  1- del*(a/b+b/a)*(1-a)*(1-b) *2/(sqr(1-a)+sqr(1-b)); };
// wmd as in BHLUMI, (what about exact mass terms???)
inline double wmd(double del,double a,double b) {
	return  1 + del*(a/b+b/a)*(a*a+b*b)/(sqr(1-a)+sqr(1-b)); };
// Factorized wm1, as in BHLUMI, to improve numerical stability.
// Identity wm1=wm7=wm0*wmd, is true up to delta**4 terms
inline double wm7(double del,double a,double b) {
	return (1 +2*del -del*(a/b+b/a))
              *(1 + del*(a*a+b*b)*(a*a+b*b)/(sqr(1-a)+sqr(1-b))/(a*b)); }
// end of INLINE functions

void SetDB(    KKdbase *DBase){ DB = DBase;};
void SetEvent( KKevent *Event){ m_Event = Event;};
void SetBornV( KKborn  *BornV){ m_BornV = BornV;};

void Initialize();
void Make();

double Dilog(double x);
void   Disr1a(double gami, int j1, double *g1, double *g2, double *gg1, double *gg2, double *ggg1, double *ggg2);
void   Disr1( double gami, int j1, double *g1, double *g2, double *gg1, double *gg2, double *ggg1, double *ggg2);

void   Disr2a(double gami, int j1,int j2, double *g1, double *g2, double *gg1, double *gg2);
void   Disr2(double gami, int jhard, int j1, int j2, double *g1, double *g2, double *gg1, double *gg2);

void   Disr3(int jhard, int j1, int j2, int j3, double *g1, double *g2);

void   Dfsr1( double gamf, int j1, double *g1, double *g2, double *gg1, double *gg2);
void   Dfsr2( int jhard, int j1, int j2, double *g1, double *g2);

void   bvirt0(double alfinv, double charg2, double svar, double am, double *dels1,  double *dels2,  double *dels3);
////////////////////////////////////////////////////////////////////////////
       ClassDef(KKqed3,1); // Data base
};// KKqed3 class
////////////////////////////////////////////////////////////////////////////
#endif
