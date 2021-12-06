///////////////////////////////////////////////////////////////////////////////
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef KKceex_H
#define KKceex_H

#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
using namespace std;

#include "KKdbase.h"
#include "KKevent.h"
#include "KKborn.h"
#include "KKdizet.h"
#include "KKpart.h"
#include "KKbvir.h"

typedef complex<double> dcmplx;

#include "BXFORMAT.h"

#include "TObject.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TRandom.h"


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//   Auxiliary class encapsulating  4-fermion spin anmplitudes               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
class KKcmplx4 : public TObject {
  public:
	dcmplx  m_A[2][2][2][2];    // amplitude
  public:
  KKcmplx4(){                        // Constructor for Foam
    //cout<< "----> KKcmplx4 Default Constructor (for ROOT only) "<<endl;
  }
  ~KKcmplx4(){
    //Explicit destructor for Foam
    //cout<< "----> KKcmplx4::~KKcmplx4 !!!! DESTRUCTOR !!!! "<<endl;
  }
  KKcmplx4(const double );           // example USER Constructor
//////////////////////////////////////////////////////////////////////////////
  #ifdef ROOT_DEF
      ClassDef(KKcmplx4,1) //n-dimensional vector with dynamic allocation
  #endif
};// KKcmplx4

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//   Auxiliary class encapsulating  2-dim. complex array                     //
///////////////////////////////////////////////////////////////////////////////
class KKcmplx2 : public TObject {
  public:
	dcmplx  m_A[2][2];    // amplitude
  public:
  KKcmplx2(){                        // Constructor for Foam
    //cout<< "----> KKcmplx2 Default Constructor (for ROOT only) "<<endl;
  }
  ~KKcmplx2(){
    //Explicit destructor for Foam
    //cout<< "----> KKcmplx2::~KKcmplx2 !!!! DESTRUCTOR !!!! "<<endl;
  }
  KKcmplx2(const double );           // example USER Constructor
//////////////////////////////////////////////////////////////////////////////
  #ifdef ROOT_DEF
      ClassDef(KKcmplx2,1) //n-dimensional vector with dynamic allocation
  #endif
};// KKcmplx2



///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//    Second order CEEX matrix element  as in original KKMC                  //
//    This is version without W t-chanel exchange  for KKMC-hh               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
//________________________________________________________________________
class KKceex: public TObject{
 public:
 ofstream *m_Out;     //! pointer to external Logfile for messages
// class member data
 public:
 KKdbase  *DB;                     // Database
 KKdizet  *m_DZ;                   // Dizet interface
 KKevent  *m_Event;                // MC event record
 KKborn   *m_BornDist;             // Born differential distribution
 KKbvir   *m_BVR;                  // Library of virtual corrections
 TRandom  *m_RNgen;
 public:
 static const int maxPhot =   101;  // max. no. of KKMC photons +1
 static const int maxWT   =  1001;  // max. num. KKMC wt list +1
 double      m_zeta2;
 double      m_zeta3;
 public:
 int         m_icont;              // event counter for debug
 int         m_nPhot;              // No of generated photons ISR+FSR
 int         m_isr[maxPhot];       // marker of ISR photons,  =1,0 for ISR,FSR, randomized!
 //int         m_isr0[maxPhot];      // Copy of partition from crude MC. OBSOLETE!!!
 double      m_Phel[maxPhot];      // photon helicity =1,0 for +,-, randomized!
 dcmplx      m_Sini[2][maxPhot];   // soft factors ISR
 dcmplx      m_Sfin[2][maxPhot];   // soft factors FSR
 dcmplx      m_Pauli[ 4][2][2];    // Pauli matrices
 dcmplx      m_SDMat3[2][2];       // spin density matrix final fermion
 dcmplx      m_SDMat4[2][2];       // spin density matrix final fermion
 KKpart      m_p1,m_p2;            // fermion momenta  initial, variable
 KKpart      m_p3,m_p4;            // fermion momenta  final,   variable
 KKpart      m_r1,m_r2;            // initial quark momenta with zero mass for Born spinors
 KKpart      m_Phot[maxPhot];      // Photons, variable
 double      m_CMSene;
 double      m_e_QED;              // QED coupling
 double      m_Alfpi;              // alpha_QED/pi
 double      m_Emin;               // minimum photon energy in CMS
 double      m_YFS_IR, m_YFSkon;   // IR formactors from Foam Density
 int         m_HasFSR;             // =KeyFSR, but later on varies
 int         m_KeyFSR;             //
 int         m_KeyISR;             //
 int         m_KeyInt;             // IFI in CEEX switch
 int         m_KeyArb;             // Type of spinor scheme
 KKpart      m_Xi;
 KKpart      m_Eta;
 KKpart      m_b;                  // auxiliary vector in spinor construction
 KKpart      m_b1;                 // auxiliary vector in spinor construction
 KKpart      m_b2;                 // auxiliary vector in spinor construction
 KKpart      m_b3;                 // auxiliary vector in spinor construction
 TLorentzRotation m_TraLor1;
 TLorentzRotation m_TraLor2;
 TLorentzRotation m_TraLin1;
 TLorentzRotation m_TraLin2;

 double      m_RhoCrud;
 double      m_ExpoNorm;
 KKcmplx4    m_Boxy;               // Encapsulated amplitude
 KKcmplx4    m_BornB;              // Encapsulated amplitude
 KKcmplx4    m_BornC;              // Encapsulated amplitude
 KKcmplx4    m_AmpExpo0;  // O(alf0)_exp
 KKcmplx4    m_AmpExpo1;  // O(alf1)_exp
 KKcmplx4    m_AmpExpo2;  // O(alf2)_exp
 KKcmplx4    m_AmpBornW;  // W exchange for e-neutrino
 KKcmplx4    m_AmpTemp;   // auxiliary temporary
 KKcmplx4    m_AmpTemp1;  // auxiliary temporary
 KKcmplx4    m_AmpTemp2;  // auxiliary temporary

 dcmplx      m_SpinoTT[2][2][2][2];
 dcmplx      m_SpinoUU[2][2][2][2];

 dcmplx      m_FFacTT[2],   m_FFacUU[2];
 dcmplx      m_FFacTG[2],   m_FFacTZ[2];
 dcmplx      m_FFacUG[2],   m_FFacUZ[2];


 dcmplx      m_F1ini1;
 dcmplx      m_F1fin1;
 dcmplx      m_F1ini2;
 dcmplx      m_F1fin2;

 dcmplx      m_IntReson;
 dcmplx      m_IntIR;
 dcmplx      m_BoxGGtu;
 dcmplx      m_BoxGZtu;
 dcmplx      m_BoxGGut;
 dcmplx      m_BoxGZut;

 double      m_RhoExp0;
 double      m_RhoExp1;
 double      m_RhoExp2;
 double      m_WtSet[maxWT];
 double      m_WtBest;

//------------------------------------
// Obligatory members
  public:
  KKceex();                    // explicit default constructor for streamer
  KKceex(ofstream *OutFile);   // user constructor
  ~KKceex();                   // explicit destructor
  public:
/////////////////////////////////////////////////////////////////////////////
// class member functions
inline double sqr( double x ){ return x*x;};
inline double qub( double x ){ return x*x*x;};

void SetDB(    KKdbase *DBase){ DB = DBase;};
void SetEvent( KKevent *Event){ m_Event = Event;};
void SetDZ(    KKdizet *DZ){    m_DZ = DZ;};
void SetBornV(KKborn *BornDist){ m_BornDist= BornDist;};
void SetBVR(KKbvir *BVR){ m_BVR= BVR;};
void SetRNgen(TRandom *RNgen){ m_RNgen= RNgen;};

void SetEmin(double Emin){      m_Emin   = Emin;};
void SetKeyInt(int KeyInt){     m_KeyInt = KeyInt;};
//void SetIR(double YFS_IR, double YFSkon){ m_YFS_IR= YFS_IR; m_YFSkon= YFSkon;};

inline double XiProd( KKpart &p, KKpart &q){
/////////////////////////////////////////////////////////////////////////////////////
//   auxiliary function called in iProd2                                           //
/////////////////////////////////////////////////////////////////////////////////////
    return sqrt( ( (p.P)[0] -(p.P)[1] ) / ( (q.P)[0] -(q.P)[1] ) );
}//XiProd


void Initialize();
void Make();
void MakeRho();
void MakeRho2(double h1[], double h2[], double &wt0, double &wt1, double &wt2);
double MakeRhoFoam();
void MakeAmpBox();
void Born(int KFini, int KFfin, TLorentzVector &PX, double CosThetD,
         KKpart &p1, KKpart &p2, KKpart &p3, KKpart &p4, KKcmplx4 &AmpBorn);
void BornW(int KFini, int KFfin, TLorentzVector &PX, double s, double t,
         KKpart &p1, KKpart &p2, KKpart &p3, KKpart &p4, KKcmplx4 &AmpBornW);

void BornPlus(int KFi, int KFf, dcmplx Cfac, TLorentzVector &PX);
// NEW!!!
void GPS_EWFFactW(int KFi, int KFf,double s,double t, dcmplx &PropW, dcmplx &WVPi);

void HiniPlus(int KFini, int KFfin, TLorentzVector &PX,
         KKpart &ph, int Hel, dcmplx &Sactu, dcmplx &sProd);
void HiniPlusW(int Ibeta, int KFini, int KFfin, TLorentzVector &PX,
         KKpart &ph, int Hel, dcmplx &Sactu, dcmplx &sProd);

double BornFoam0(int KFini, int KFfin, double SvarX, double CosThetD);
void   BornFoam2(int Mode, int KFini, int KFfin, double SvarX, double CosThetD, double &Yint);

void HfinPlus(int KFini, int KFfin, TLorentzVector &PX,
         KKpart &ph,  int Hel, dcmplx &Sactu, dcmplx &sProd, double &CKine);
void GPS_HiiPlus(dcmplx Cfact2, int KFini, int KFfin, TLorentzVector PX,
         KKpart ph1, int Hel1, KKpart ph2, int Hel2);
void GPS_HiiPlusW(dcmplx Cfact2, int KFini, int KFfin, TLorentzVector PX,
         KKpart ph1, int Hel1, KKpart ph2, int Hel2);
void GPS_HffPlus(dcmplx CNorm, int KFini, int KFfin, TLorentzVector PX,
		 KKpart ph1, int Hel1, KKpart ph2, int Hel2);
void GPS_HifPlus(dcmplx CNorm, int KFini, int KFfin, TLorentzVector PX,
         KKpart ph1, int Hel1, KKpart ph2, int Hel2);

void Amp2Mult(KKcmplx2 &Res, dcmplx Fact, const KKcmplx2 &V, const KKcmplx2 &U);
void Amp2WCplus( KKcmplx2 &V, const double sw2);
void Amp2WCminus(KKcmplx2 &V, const double Sw2);
void Amp2to4(KKcmplx4 &Res, dcmplx Fact, const KKcmplx2 &V, const KKcmplx2 &U);
void Amp2to4add(KKcmplx4 &Res, dcmplx Fact, const KKcmplx2 &V, const KKcmplx2 &U);

void Amp4Zer(KKcmplx4 &Born);
void AmpAdd( KKcmplx4 &Work, dcmplx Fact, const KKcmplx4 &Born);
void AmpAdd(KKcmplx4 &Work, dcmplx Fact,
       const KKcmplx2 &U, const KKcmplx2 &VX, const KKcmplx2 &V, const KKcmplx2 &UX);
void AmpAddI(KKcmplx4 &Work, dcmplx Fact, const KKcmplx4 &Born, const KKcmplx2 &U);
void AmpAddI(KKcmplx4 &Work, dcmplx Fact, const KKcmplx2 &V, const KKcmplx4 &Born);
void AmpAddI(KKcmplx4 &Work, dcmplx Fact, const KKcmplx2 &V, const KKcmplx4 &Born,  const KKcmplx2 &U);
void AmpAddI(KKcmplx4 &Work, dcmplx Fact, const KKcmplx4 &Born, const KKcmplx2 &U1, const KKcmplx2 &U2);
void AmpAddI(KKcmplx4 &Work, dcmplx Fact, const KKcmplx2 &V1, const KKcmplx2 &V2, const KKcmplx4 &Born);

void AmpAddF(KKcmplx4 &Work, dcmplx Fact, const KKcmplx4 &Born, const KKcmplx2 &V);
void AmpAddF(KKcmplx4 &Work, dcmplx Fact, const KKcmplx2 &U, const KKcmplx4 &Born);
void AmpAddF(KKcmplx4 &Work, dcmplx Fact, const KKcmplx2 &U, const KKcmplx4 &Born, const KKcmplx2 &V);
void AmpAddF(KKcmplx4 &Work, dcmplx Fact, const KKcmplx4 &Born, const KKcmplx2 &V1, const KKcmplx2 &V2);
void AmpAddF(KKcmplx4 &Work, dcmplx Fact, const KKcmplx2 &U1, const KKcmplx2 &U2, const KKcmplx4 &Born);

void EWFFact(int KFini, int KFfin, double Svar, dcmplx &Ve, dcmplx &Vf, dcmplx &Ae, dcmplx &Af );
void EWFFact(const int KFini, const int KFfin, const double Svar, const double CosThetD,
		         dcmplx &Ve,    dcmplx &Vf,     dcmplx  &Ae,     dcmplx &Af,
		         dcmplx &VVcor, dcmplx &GamVPi, dcmplx  &ZetVPi, double &RsqV, double & RsqA);

dcmplx iProd1(int L, KKpart &p, KKpart &q);
dcmplx iProd2(KKpart &p, KKpart &q);
dcmplx iProd2(int Lp, KKpart &p, int Lq, KKpart &q);
dcmplx iProd2(int Cp, int Lp, KKpart &p, int Cq, int Lq, KKpart &q);
//another version of iProd2
dcmplx iProd3(int Lamp, KKpart &p, double mp, int Lamq, KKpart &q, double mq);

dcmplx Soft( int sigma, KKpart &ph, KKpart &p1, KKpart &p2);
dcmplx Bfacb(int sigma, KKpart &phot, KKpart &pferm);

void   GPS_MakeU(dcmplx Cfac, KKpart &ph, int sigma, KKpart &p1, KKpart &p2,  KKcmplx2 &U);
void   GPS_MakeV(dcmplx Cfac, KKpart &ph, int sigma, KKpart &p1, KKpart &p2,  KKcmplx2 &V);
dcmplx GPS_Sof1(int sigma, KKpart &ph, KKpart &pf);
dcmplx GPS_Sof1x(int sigma, KKpart &ph, KKpart &pf);
void   GPS_Make_eps(KKpart &ph, int Sig, dcmplx eps[4]);

void GPS_MakeUW(dcmplx Cfac, KKpart &ph, int sigma, KKpart &p1, KKpart &p2,  KKcmplx2 &U);
void GPS_MakeVW(dcmplx Cfac, KKpart &ph, int sigma, KKpart &p1, KKpart &p2,  KKcmplx2 &V);
void GPS_MatrS(KKpart &p1, KKpart &p2, KKcmplx2 &V);
void GPS_MakeUX(dcmplx Cfac, KKpart &ph, KKpart &p1, KKpart &p2,  KKcmplx2 &U);
void GPS_MakeVX(dcmplx Cfac, KKpart &ph, KKpart &p1, KKpart &p2,  KKcmplx2 &V);

// soft formfactors
double SForFac(double alfpic, KKpart &p1, KKpart &p2, double Emin, double MasPhot);
double TForFac(double alfpic, KKpart &p1, KKpart &p2, double Emin, double MasPhot);
// virtual corrections
void   MakeF1ini(double Svar, double Mas1, double Mas2, double Q, dcmplx &F1_1, dcmplx &F1_2);
void   MakeF1fin(double Svar, double Mas1, double Mas2, double Q, dcmplx &F1_1, dcmplx &F1_2);
dcmplx CnuA(double Svar, double Mas1, double Mas2);
void   MakeVini(KKpart p1, KKpart p2, KKpart ph, dcmplx &V1, dcmplx &V2);
void   MakeVfin(KKpart p3, KKpart p4, KKpart ph, dcmplx &V1, dcmplx &V2);

void PhelRandom();
void PartitionStart(int &last);
void PartitionPlus( int &last);

void ZerAmplit();
void TralorPrepare(int KTO);
void TralorDoIt(int KTO, double P[], double Q[]);
////////////////////////////////////////////////////////////////////////////
       ClassDef(KKceex,1); // Data base
};// KKceex class
////////////////////////////////////////////////////////////////////////////
#endif
