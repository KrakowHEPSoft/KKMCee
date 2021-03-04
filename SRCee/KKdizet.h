///////////////////////////////////////////////////////////////////////////////
//         Template of the class with ROOT persistency
///////////////////////////////////////////////////////////////////////////////

#ifndef KKdizet_H
#define KKdizet_H

#include <stdlib.h>
#include <math.h>
#include <complex>
using namespace std;

#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

typedef complex<double> dcmplx;

#include "BXFORMAT.h"
#include "TObject.h"
//________________________________________________________________________
class KKdizet: public TObject{
 public:
 ofstream *m_Out;     //! pointer to external Logfile for messages
// class member data
 public:
 double   CMSene;

 // see KKDdizet.cxx for initialization
 static const double  m_WminLEP1;  // LEP1 basic range
 static const double  m_WmaxLEP1;  // LEP1 basic range
 static const double  m_WdelZ   ;  // Z range (amz +- m_WdelZ)
 static const double  m_WmaxLEP2;  // LEP2 interval (m_WmaxLEP1,m_WmaxLEP2)
 static const double  m_WmaxNLC ;  // LHC/NLC range (m_WmaxLEP2,m_WmaxNLC)
 double m_WminZ, m_WmaxZ;

 static const int m_poinG  =  7;   // No of EW  form-factors
 static const int m_poinQ =   4;   // No of QCD form-factors

 static const int m_poin1  =100;   // No. of sqrt(s) points
 static const int m_poTh1 =   0;   // Just one point for all theta's

 static const int m_poin2 = 120;   // Z range sqrt(s) spacing
 static const int m_poTh2 =  14;   // theta points =14 is it overkill?

 static const int m_poin3 = 145;   // LEP2 interval sqrt(s) spacing
 static const int m_poTh3 =  30;   // Cost(heta) Overkill, but lets keep it

 static const int m_poin4 = 180;   // NLC range sqrt(s)  spacing
 static const int m_poTh4 =  14;   // Cost(heta) spacing

 dcmplx  m_cyys[m_poin1+1]           [m_poinG][5][16];  // EW form-fact. table
 dcmplx  m_czzs[m_poin2+1][m_poTh2+1][m_poinG][5][16];  // EW form-fact. table
 dcmplx  m_ctts[m_poin3+1][m_poTh3+1][m_poinG][5][16];  // EW form-fact. table
 dcmplx  m_clcs[m_poin4+1][m_poTh4+1][m_poinG][5][16];  // EW form-fact. table

 // QCD formaftors are absent in case of leptonic final states and here are not used
 // but let us keep them for the KKMC-ee version
 double  m_syys[m_poin1+1][           m_poinQ][5][16];  // QCD correction, OBSOLETE!!!
 double  m_szzs[m_poin2+1][           m_poinQ][5][16];  // QCD correction, OBSOLETE!!!
 double  m_stts[m_poin3+1][           m_poinQ][5][16];  // QCD correction, OBSOLETE!!!
 double  m_slcs[m_poin4+1][           m_poinQ][5][16];  // QCD correction, OBSOLETE!!!

 dcmplx m_GSW[    m_poinG];    // form-factors,   at the actual energy/angle
 double m_QCDcorR[m_poinQ];    // QCD correction, at the actual energy,    OBSOLETE!!!

 double m_MZ;         // SM input mass of Z
 double m_amh;        // SM input mass of Higgs
 double m_amtop;      // SM input mass of top
 double D_swsq;       // SM calculated EW mixing angle
 double D_GamZ;       // SM calculated Z width
 double D_MW;         // SM calculated mass of W
 double D_GammW;      // SM calculated width of W

 int m_KeyQCD;                 // To be connected to input or removed!!!
//------------------------------------
// Obligatory members
  public:
  KKdizet();                    // explicit default constructor for streamer
  KKdizet(ofstream *OutFile);   // user constructor
  ~KKdizet();                   // explicit destructor
  public:
/////////////////////////////////////////////////////////////////////////////
// class member functions
double sqr( const double x );

void Initialize();
void ReadEWtabs();
void InterpoGSW(int KFi, int KFf, double svar, double CosThe);
void GetGSWxy(double[],double[]);
////////////////////////////////////////
// for tests
//void ImportEWtabs();

////////////////////////////////////////////////////////////////////////////
       ClassDef(KKdizet,1); // Data base
};// KKdizet class
////////////////////////////////////////////////////////////////////////////
#endif
