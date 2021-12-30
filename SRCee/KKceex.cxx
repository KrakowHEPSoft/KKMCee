///////////////////////////////////////////////////////////////////////////////
#include "KKceex.h"

ClassImp(KKceex);
ClassImp(KKcmplx4);
ClassImp(KKcmplx2);

#define SW20 setw(20)<<setprecision(14)
#define SW208 setw(30)<<setprecision(8)

extern "C" {
//
   void pseumar_initialize_(const int&, const int&, const int&);
   void pseumar_makevec_(float rvec[], const int&);
}//

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//   Auxiliary classes encapsulating  complex arrays                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
KKcmplx4::KKcmplx4(double x)
{
  cout<< "----> KKcmplx4 USER Constructor "<<endl;
  for(int j1 = 0; j1<=1; j1++)
    for(int j2 = 0; j2<=1; j2++)
      for(int j3 = 0; j3<=1; j3++)
        for(int j4 = 0; j4<=1; j4++){
              m_A[j1][j2][j3][j4] =  dcmplx(x,0.0);
        }
//----------------------
}//KKcmplx4

KKcmplx2::KKcmplx2(double x)
{
  cout<< "----> KKcmplx2 USER Constructor "<<endl;
  for(int j1 = 0; j1<=1; j1++)
    for(int j2 = 0; j2<=1; j2++){
              m_A[j1][j2] =  dcmplx(x,0.0);
        }
//----------------------
}//KKcmplx2



//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//                                                                                  //
//                       Main KKceex Class                                          //
//                                                                                  //
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////


KKceex::KKceex()
{
  // This constructor is for ROOT streamers ONLY
  // all pointers has to be NULLed
  cout<< "----> KKceex Default Constructor (for ROOT only) "<<endl;
  m_Out  = NULL;
  DB     = NULL;
  m_Event= NULL;
  m_DZ   = NULL;
  m_BornDist = NULL;
  m_BVR  = NULL;
}

///_____________________________________________________________
KKceex::KKceex(ofstream *OutFile)
{
  cout<< "----> KKceex USER Constructor "<<endl;
  m_Out = OutFile;
  DB     = NULL;
  m_Event= NULL;
  m_DZ   = NULL;
  m_BornDist = NULL;
  m_BVR  = NULL;
}//KKceex

///______________________________________________________________________________________
KKceex::~KKceex()
{
  //Explicit destructor
  cout<< "----> KKceex::KKceex !!!! DESTRUCTOR !!!! "<<endl;
}///destructor




///______________________________________________________________________________________
void KKceex::Initialize()
{
  cout  << "----> KKceex::Initialize, Entering "<<endl;
//=================================================================
// BX*** macros are in MCdev/BXFORMAT.h
  BXOPE(*m_Out);
  BXTXT(*m_Out,"========================================");
  BXTXT(*m_Out,"======    KKceex::Initialize      ======");
  BXTXT(*m_Out,"========================================");
  cout<<       "========================================"<<endl;
  cout<<       "======    KKceex::Initialize      ======"<<endl;
  cout<<       "========================================"<<endl;
  m_icont =0;

  m_zeta2  = M_PI*M_PI/6.0;
  m_zeta3  = 1.2020569031595942854;

  m_KeyInt = DB->KeyINT; // this is overwritten anyway
  m_KeyISR = DB->KeyISR;
  m_KeyFSR = DB->KeyFSR;

  m_e_QED  = sqrt( 4.0*M_PI/DB->Alfinv0);
  m_Alfpi  = 1.0/M_PI/DB->Alfinv0;
/////////////////////////////////////////////////////////////////////////////////////
// Key for switching on/off the use of m_b, can be reset with GPS_SetKeyArb
  m_KeyArb = 0;  // default zero value, m_b=Xi is assumed

  m_BornC = KKcmplx4(0.0);


//     Define Pauli matrices
  for(int k = 0; k<=3; k++)
     for(int j1=0; j1<=1; j1++)
        for(int j2 = 1; j2<=1; j2++)
           m_Pauli[ k][j1][j2] = dcmplx(0.0,0.0);
// Sigma0
  m_Pauli[ 0][0][0] = dcmplx( 1.0, 0.0);
  m_Pauli[ 0][1][1] = dcmplx( 1.0, 0.0);
// SigmaX
  m_Pauli[ 1][0][1] = dcmplx( 1.0, 0.0);
  m_Pauli[ 1][1][0] = dcmplx( 1.0, 0.0);
// SigmaY
  m_Pauli[ 2][0][1] = dcmplx( 0.0,-1.0);
  m_Pauli[ 2][1][0] = dcmplx( 0.0, 1.0);
// SigmaZ
  m_Pauli[ 3][0][0] = dcmplx( 1.0, 0.0);
  m_Pauli[ 3][1][1] = dcmplx(-1.0, 0.0);

//--- basic vectors for spin quantization of spinors
  m_Xi  = KKpart(1.0, 1.0, 0.0, 0.0); // zero component is time
  m_Eta = KKpart(0.0, 0.0, 1.0, 0.0); //
// gauge fixing b-vector for constructing spinors
  m_b1 = KKpart(0.0,  0.8723e0, -0.7683e0, 0.3348e0);
  m_b1[ 0]   =  sqrt( sqr(m_b1[ 1]) + sqr(m_b1[ 2]) +sqr(m_b1[ 3]) );
// another random setting
  m_b2 = KKpart(0.0, -0.78833e0, 34788e0, 33282e0 );
  m_b2[ 0]   =  sqrt( sqr(m_b2[ 1]) + sqr(m_b2[ 2]) +sqr(m_b2[ 3]) );
// this setting is very close to b=Xi, for special tests
  m_b3 = KKpart(0.0,  1.0, 1e-7, 0.0);
  m_b3[ 0]   =  sqrt( sqr(m_b3[ 1]) + sqr(m_b3[ 2]) +sqr(m_b3[ 3]) );
//  This is actual assignment
  m_b = m_b1;

  ///////////////////////////////////////////////////
}// Initialize

///______________________________________________________________________________________
void KKceex::Make(){
  m_icont++;
  int    KFini   = m_Event->m_KFini;
  int    KFfin   = m_Event->m_KFfin;

  m_HasFSR = DB->KeyFSR;          // general FSR switch
  m_HasFSR = m_Event->m_HasFSR;   // overruled by KKarfin (exception for neutrinos)

  double Mbeam = DB->fmass[KFini];  // to be refined, current or constituent?
  double Massf = DB->fmass[KFfin];
  int      NCf = DB->Nc[KFfin];

  // private vectors used as arguments in spinor products
  //SetAll( const int C0, const int Hel0,const double M0, TLorentzVector &P0){
  double fleps = 1e-20;
  int hel1 =1, hel2= -1, hel3 = 1, hel4=1;
  m_p1.SetAll(  1, hel1, Mbeam,  m_Event->m_Pf1);  // e-    U-spinor
  m_p2.SetAll( -1, hel2, Mbeam,  m_Event->m_Pf2);  // e+    V-bar spinor
  m_p3.SetAll(  1, hel3, Massf,  m_Event->m_Qf1);  // f     U-bar-spinor
  m_p4.SetAll( -1, hel4, Massf,  m_Event->m_Qf2);  // f-bar V-spinor
  // Special version of initial momenta with infinitesimal mass for Born spinors
  m_r1.SetAll(  1, hel1, fleps,  m_Event->m_Pf1);  // e-    U-spinor
  m_r2.SetAll( -1, hel2, fleps,  m_Event->m_Pf2);  // e+    V-bar spinor

  double  C = 0; // for photon used as spinors product argument it has to be +-1 !!!
  double MasPhot = DB->MasPhot;
  m_nPhot = m_Event->m_nPhot;
  for(int i=1; i<= m_nPhot; i++){
      m_isr[i]  = m_Event->m_isr[i];    // f77 indexing
      m_Phot[i].SetAll( C, m_Phel[i], MasPhot, m_Event->m_PhotAll[i]);
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /*
  if( m_icont <= 20) {
  // fermion helicities still undefined
  // SetAll(const double M0, const int C0, const int Hel0, TLorentzVector &P0);
  cout<<"##############################################KKceex::Make#############################################"<<endl;
  //  cout<<"m_p1=";  m_p1.Print();
  //  cout<<"m_p2=";  m_p2.Print();
  //  cout<<"m_p3=";  m_p3.Print();
  //  cout<<"m_p4=";  m_p4.Print();
    double  Xp1p3= XiProd(m_p1, m_p3);
    double  Xp2p4= XiProd(m_p2, m_p4);
    double  Xp1p4= XiProd(m_p1, m_p4);
    double  Xp2p3= XiProd(m_p2, m_p3);
    cout<< "XiProd: Xp1p3="<< Xp1p3<< " Xp2p4="<< Xp2p4<<" Xp1p4="<< Xp1p4<<" Xp2p3="<< Xp2p3<<endl;
    //dcmplx KKceex::iProd1(int L, KKpart &p, KKpart &q){
    dcmplx Ip3p1= iProd1(hel3,m_p3, m_p1);
    dcmplx Ip2p4= iProd1(hel2,m_p2, m_p4);
    dcmplx Ip3p2= iProd1(hel3,m_p3, m_p2);
    dcmplx Ip1p4= iProd1(hel1,m_p1, m_p4);
    cout<< "iProd1: Ip3p1="<< Ip3p1<<" Ip2p4="<< Ip2p4<<" Ip3p2="<< Ip3p2<<" Ip1p4="<< Ip1p4<<endl;
    dcmplx Yp3p1= iProd2(m_p3, m_p1);
    dcmplx Yp2p4= iProd2(m_p2, m_p4);
    dcmplx Yp3p2= iProd2(m_p3, m_p2);
    dcmplx Yp1p4= iProd2(m_p1, m_p4);
    cout<< "iProd2: Yp3p1="<< Yp3p1<<" Yp2p4="<<Yp2p4<< " Yp3p2="<< Yp3p2<<" Yp1p4="<< Yp1p4<<endl;
    // helicities as explitit arguments
    // iProd2(int Hp, KKpart &p, int Hq, KKpart &q){
    dcmplx Qp3p1= iProd2(hel3,m_p3, hel1,m_p1);
    dcmplx Qp2p4= iProd2(hel2,m_p2, hel4,m_p4);
    dcmplx Qp3p2= iProd2(hel3,m_p3, hel2,m_p2);
    dcmplx Qp1p4= iProd2(hel1,m_p1, hel4,m_p4);
    cout<< "iProd2: Qp3p1="<< Qp3p1<<" Qp2p4="<<Qp2p4<< " Yp3p2="<< Yp3p2<<" Yp1p4="<< Yp1p4<<endl;
    // All in arguments
    //iProd2(int Cp, int Lp, KKpart &p, int Cq, int Lq, KKpart &q)
    dcmplx Up2p4, Up3p2, Up1p4, Up3p1;
    Up3p1= iProd2( 1, hel3, m_p3,  1, hel1, m_r1);
    Up2p4= iProd2(-1, hel2, m_r2, -1, hel4, m_p4);
    Up3p2= iProd2( 1, hel3, m_p3, -1, hel2, m_r2);
    Up1p4= iProd2( 1, hel1, m_r1, -1, hel4, m_p4);
    cout<< "iProd2: Up3p1="<< Up3p1<<" Up2p4="<<Up2p4<< " Up3p2="<< Up3p2<<" Up1p4="<< Up1p4<<endl;
// initial masses are fleps!!!
    dcmplx Zp2p4, Zp3p2, Zp1p4, Zp3p1;
    double m3 = Massf, m4=Massf;
    Zp3p1 =iProd3( hel3, m_p3, m3,    hel1, m_p1, fleps); //! t
    Zp2p4 =iProd3( hel2, m_p2,-fleps, hel4, m_p4,-m4);    //! t1
    Zp3p2 =iProd3( hel3, m_p3, m3,    hel2, m_p2,-fleps); //! u1
    Zp1p4 =iProd3( hel1, m_p1, fleps, hel4, m_p4,-m4);    //! u
    cout<< "iProd2: Zp3p1="<< Zp3p1<<" Zp2p4="<<Zp2p4<< " Zp3p2="<< Zp3p2<<" Zp1p4="<< Zp1p4<<endl;
    }//if
 */
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////

//-------------------------------------------------------------
  TLorentzVector PP = m_Event->m_Pf1 + m_Event->m_Pf2;
  TLorentzVector QQ = m_Event->m_Qf1 + m_Event->m_Qf2;
  double svar  = PP*PP;
  double svarQ = QQ*QQ;
  double Ene   = sqrt(svar)/2.0;

  m_Emin      = Ene *DB->vvmin; //IR cutoff
//-------------------------------------------------------------
// Overall normalization factors
  double CrudNorm  =  1.0;
  m_ExpoNorm  =  2.0/(4.0*M_PI)*NCf;  //! it is still quasi-empirical...

////////////////////////////////////////////////////////////////////////////////////
//                         YFS  FormFactors                                       //
//  Note that FSR formfactor below cannot be used for KeyPia=0, Emin is in CMS!!! //
////////////////////////////////////////////////////////////////////////////////////
  double ChaIni=DB->Qf[KFini];
  double ChaFin=DB->Qf[KFfin];
  double alfpini  = m_Alfpi* sqr(ChaIni);
  double alfpfin  = m_Alfpi* sqr(ChaFin);
  double alfpmix  = m_Alfpi*ChaFin*ChaIni;
  double YFS_IRfin, YFS_IRini, YFSkonIni, YFSkonFin, YFS_isr, Yisr;
  if( DB->KeyISR != 0) {
// m_YFS_IR and m_YFSkon provided from outside using SetIR()
//    YFS_IRini = m_YFS_IR; // ISR param. imported from Foam Density
//    YFSkonIni = m_YFSkon; // !!!<<**  YFSkon not in WtCrud (historical reasons)
    ///////// Imported for Density function through Event
    YFS_IRini = m_Event->m_YFS_IR_ini;
    YFSkonIni = m_Event->m_YFSkon_ini;   //!!!<<** YFSkon not in WtCrud (historical reasons)
    /////////
    CrudNorm  = CrudNorm *YFS_IRini;
    Yisr= SForFac( alfpini, m_p1, m_p2, m_Emin, MasPhot);
    m_ExpoNorm  = m_ExpoNorm *Yisr;
  }//if KeyISR
//////////////////////////////////////////////////////////////////////
  if( m_HasFSR != 0) {
/*
    double QQk   = QQ*m_Event->m_PhotFSR[0];
    double Delta = DB->vvmin *DB->delfac;
    double Delta1 = Delta*(1+ 2*QQk/svarQ);
    double q1q2  =  m_p3*m_p4;
    YFS_IRfin = -2.0* alfpfin *( q1q2 *m_BVR->A( q1q2, Massf,Massf) -1.0  ) *log(1/Delta1);
    double Eqq   = 0.5*sqrt(svarQ);
    double EminQ  = Eqq*Delta1;
    double DelYFS =
        m_BVR->Btilda(alfpfin, q1q2, m_p3[0],m_p4[0], Massf, Massf,  m_Emin, MasPhot) //!exact
       -m_BVR->Btilda(alfpfin, q1q2,   Eqq,    Eqq,   Massf ,Massf,   EminQ, MasPhot); //!exact
    if(DB->KeyPia == 1 ) {
       YFS_IRfin = exp(YFS_IRfin +DelYFS);
    } else {
        YFS_IRfin = exp(YFS_IRfin);
    }//if
    YFSkonFin =  1/4.0 *2*  alfpfin *(log(svarQ/sqr(Massf))-1)
                          + alfpfin*( -.5  +sqr(M_PI)/3.0); //! Mass<<sqrt(s)
    YFSkonFin =  exp(YFSkonFin);         //!!!<<** YFSkon not in WtCrud (historical reasons)
*/
/////// the above repeats calculations in KKarfin:Piatek correctly,
/////// but it is saver to import them through Event
    YFS_IRfin = m_Event->m_YFS_IR_fin;
    YFSkonFin = m_Event->m_YFSkon_fin;   //!!!<<** YFSkon not in WtCrud (historical reasons)
///////
    CrudNorm  = CrudNorm *YFS_IRfin;
    double Yfsr= SForFac( alfpfin, m_p3, m_p4, m_Emin, MasPhot);
    m_ExpoNorm  = m_ExpoNorm *Yfsr;
  }// if m_HasFSR

// Remember Yint depends on Emin and provides angular asymmetry (MasPhot is dummy)
  if(  m_KeyInt != 0 && DB->KeyISR != 0 &&  m_HasFSR != 0  ) {
    double Yint= TForFac( alfpmix, m_p1, m_p3, m_Emin, MasPhot)
                *TForFac( alfpmix, m_p2, m_p4, m_Emin, MasPhot)
                *TForFac(-alfpmix, m_p1, m_p4, m_Emin, MasPhot)
                *TForFac(-alfpmix, m_p2, m_p3, m_Emin, MasPhot);
    m_ExpoNorm  = m_ExpoNorm *Yint;
  }//if

//////////////////////////////////////////////////////////////////
//                       S-factors                              //
//////////////////////////////////////////////////////////////////
// List of soft-factors, note that we calculate them for helicity=+1
// The other one helicity=-1 is just minus complex conjugate!
/// Sig = 3-2*Hel, Hel=1,2 --> Sig=1,-1
  KKpart ph;
  for(int j=1; j<=m_nPhot;j++){  // f77 indexing
    CrudNorm = CrudNorm *1.0/qub(2.0*M_PI);    //<-- photon phase-space factor
    m_ExpoNorm = m_ExpoNorm *1.0/qub(2.0*M_PI);    //<-- photon phase-space factor
    ph = m_Phot[j];
    m_Sini[0][j]  =  dcmplx(ChaIni*m_e_QED) *Soft(  1,ph,m_p1,m_p2);
    m_Sfin[0][j]  = -dcmplx(ChaFin*m_e_QED) *Soft(  1,ph,m_p3,m_p4);
    m_Sini[1][j]  = -conj(m_Sini[0][j]);
    m_Sfin[1][j]  = -conj(m_Sfin[0][j]);
    //cout<<"******** m_Sini[0][j]="<<m_Sini[0][j]<<"   m_Sini[1][j]="<<m_Sini[1][j]<<endl;
    //cout<<"******** m_Sfin[0][j]="<<m_Sfin[0][j]<<"   m_Sfin[1][j]="<<m_Sfin[1][j]<<endl;
  }//for
//////////////////////////////////////////////////////////////////
//             Define (randomly) photon helicities              //
//////////////////////////////////////////////////////////////////
//[[[[ mooved up
//    PhelRandom();
//]]]]

//////////////////////////////////////////////////////////////////////
//  *************************************************************   //
//              LOOP OVER PARTITIONS STARTS HERE                    //
//  Initialize loop over partitions, m_isr(j)=1,0 denotes isr,fsr   //
//  *************************************************************   //
//////////////////////////////////////////////////////////////////////
  double CrudSum = 0.0;
  double BornCru, fLLux, betaf, DistCru;
  int last, Hel;
  ZerAmplit();
  //
  PartitionStart(last);
//=====================
  TLorentzVector PX;
  dcmplx sProd, Sactu, Cfact0;
  double svarX;
  for(int loop=1; loop<10000000; loop++){
  /////////////////////////////////////////////////////////
  //         ===============================             //
  //         Soft-part, soft-factors, m-zero             //
  //         ===============================             //
  /////////////////////////////////////////////////////////
  PX = m_Event->m_Pf1 + m_Event->m_Pf2;
  sProd=dcmplx(1.0,0.0);
  for(int j=1; j<= m_nPhot; j++){
    Hel    = m_Phel[j]; //    Phel =0,1 !!!
    Sactu  = dcmplx(m_isr[j])*m_Sini[Hel][j] + dcmplx(1-m_isr[j])*m_Sfin[Hel][j];
    sProd  = sProd*Sactu;
//    Calculate reduced 4-momentum to be used in gamma/Z propagator
    PX = PX  -m_isr[j]* m_Event->m_PhotAll[j];
  }// for j
  svarX = PX*PX;
  Cfact0 = sProd  *(svarX/svarQ);

  BornPlus(KFini, KFfin, Cfact0, PX); // this is for KeyELW <=0 !!!
// =======================================================================

  BornCru = 4.0/3.0 *m_BornDist->BornSimple(KFini, KFfin, svarX, 0.0); //  4/3 = <(1+c^2)> where c in (-1,i)
  BornCru = BornCru *(svar/svarX);                                     //<-- Born(svar)*svar
  fLLux   = svarX/svarQ;                                               //<-- extra LL factor
  betaf   = sqrt( 1.0 - 4* sqr(Massf)/svarQ );                         //<-- 2-body phase spase
  DistCru = BornCru/(4.0*M_PI) *fLLux *2.0/betaf * sqr(abs(sProd));    //<-- CRUDE
  CrudSum =       CrudSum  + DistCru;

  if(m_nPhot == 0) goto e300;
  //[[[[[[[[[[[[[[[[[[[[test
  //   Amp4Zer(m_AmpExpo2);
  //]]]]]]]]]]]]]]]]]]]]]
  int Hel1,Hel2;
  KKpart ph1,ph2;
  dcmplx Sactu1,Sactu2,SactuA,SactuB;
  dcmplx Cfact2;
  double svarX1, CKine;
  TLorentzVector QQ1;
  for(int j1=1; j1<=m_nPhot;j1++){
     ph1= m_Phot[j1];
     if ( ph1.P[0]/Ene > (DB->Vcut[0]) ) { // accept 1 hard only
        Hel1    = m_Phel[j1];
        Sactu1  = double(m_isr[j1])*m_Sini[Hel1][j1] + double(1-m_isr[j1])*m_Sfin[Hel1][j1];
        QQ1 = QQ+m_Event->m_PhotAll[j1];
        svarX1 = QQ1*QQ1;
        CKine   = (svarX1/svarQ);
        if( m_isr[j1] == 1) {
          HiniPlus(   KFini, KFfin, PX, ph1, Hel1, Sactu1, sProd);
          HiniPlusW(1,KFini, KFfin, PX, ph1, Hel1, Sactu1, sProd);
        } else {
          HfinPlus(KFini, KFfin, PX, ph1, Hel1, Sactu1, sProd, CKine);
        }
     }// if Vcut
  }//for j1

  for(int j1=1; j1<=m_nPhot;j1++)
     for(int j2=j1+1; j2<=m_nPhot;j2++){
        ph1= m_Phot[j1];
        ph2= m_Phot[j2];
        if( (ph1.P[0]/Ene > (DB->Vcut[1])) && (ph2.P[0]/Ene > (DB->Vcut[1])) ) {
           Hel1  = m_Phel[j1];
           Hel2  = m_Phel[j2];
           Sactu2  = ( double(m_isr[j1])*m_Sini[Hel1][j1] + double(1-m_isr[j1])*m_Sfin[Hel1][j1] )
                    *( double(m_isr[j2])*m_Sini[Hel2][j2] + double(1-m_isr[j2])*m_Sfin[Hel2][j2] );
           Cfact2 = sProd/Sactu2;
           if(       (m_isr[j1] == 1) && (m_isr[j2] == 1) ) {                 // ini-ini
              GPS_HiiPlus( Cfact2, KFini, KFfin, PX, ph1, Hel1, ph2, Hel2);
              SactuA  = double(m_isr[j1])*m_Sini[Hel1][j1];
              SactuB  = double(m_isr[j2])*m_Sini[Hel2][j2];
              HiniPlusW(-1,KFini, KFfin, PX, ph1, Hel1, SactuA, sProd);
              HiniPlusW(-1,KFini, KFfin, PX, ph2, Hel2, SactuB, sProd);
              GPS_HiiPlusW(Cfact2, KFini, KFfin, PX, ph1, Hel1, ph2, Hel2);
           }else if( (m_isr[j1] == 0) && (m_isr[j2] == 0) ) {                 // fin-fin
              GPS_HffPlus(Cfact2, KFini, KFfin, PX, ph1, Hel1, ph2, Hel2);
           }//if
        }//if vcut1
  }// for j1,j2

//   O(alf2) correction to 1 photon ISR and 1 photon FSR, explicit Bose-Einstein symmetrisation
  for(int j1=1; j1<=m_nPhot;j1++)
     for(int j2=1; j2<=m_nPhot;j2++){
        if( j1 != j2) {
          ph1= m_Phot[j1];
          ph2= m_Phot[j2];
          if( (ph1.P[0]/Ene > (DB->Vcut[1])) && (ph2.P[0]/Ene > (DB->Vcut[1])) ) { // accept 2 hard only
            Hel1  = m_Phel[j1];
            Hel2  = m_Phel[j2];
            Sactu2  = ( double(m_isr[j1])*m_Sini[Hel1][j1] + double(1-m_isr[j1])*m_Sfin[Hel1][j1] )
                     *( double(m_isr[j2])*m_Sini[Hel2][j2] + double(1-m_isr[j2])*m_Sfin[Hel2][j2] );
            Cfact2 = sProd/Sactu2;
            if(     (m_isr[j1] == 1) && (m_isr[j2] == 0) ) { // ini-fin
                   //CALL GPS_HifPlus(Cfact2,KFi,KFf,PX,p1,m1,p2,m2,p3,m3,p4,m4,Hel1,ph1,Hel2,ph2,mph,m_AmpExpo2) !
                   GPS_HifPlus(Cfact2, KFini, KFfin, PX, ph1, Hel1, ph2, Hel2);
            }//if
          }//if
        }//if
     }//for j1,j2

/////////////////////////////////////////////////////////
//   Update m_isr, check if it is the last partition   //
/////////////////////////////////////////////////////////
  if(last == 1) goto e300;
  PartitionPlus(last);
//
    if(last == 2) goto e300;
  }//for loop
  cout<<"########### INCORRECT EXIT from loop over partitions"<<endl;
  exit(95);
e300:
//////////////////////////////////////////////////////////////////////
//  *************************************************************   //
//              LOOP OVER PARTITIONS ENDS HERE                      //
//  *************************************************************   //
//////////////////////////////////////////////////////////////////////

  m_RhoCrud = CrudSum *CrudNorm;   //<-- Crude (unpolarized for the time being)

  MakeRho();      //<-- Defines m_RhoExp0, m_RhoExp1, m_RhoExp2

  if( m_KeyInt == 0 ) {
     m_WtSet[ 51] =   m_RhoExp0 /m_RhoCrud;    //!!! Interference OFF
     m_WtSet[ 52] =   m_RhoExp1 /m_RhoCrud;
     m_WtSet[ 53] =   m_RhoExp2 /m_RhoCrud;
     m_WtBest     =   m_WtSet[ 53];
//[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
//     (*m_Out)<< ">>>KKceex::Make:  CrudSum= "<<CrudSum << " CrudNorm= "<<CrudNorm<<endl;
//     (*m_Out)<< ">>>KKceex::Make:  m_RhoCrud= "<<m_RhoCrud<<endl;
//     (*m_Out)<< ">>>KKceex::Make:  m_WtSet[ 51,52,53] ="<< m_WtSet[ 51]<<"  "<< m_WtSet[ 52]<<"  "<< m_WtSet[ 53]<<endl;
//]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
  } else {
     m_WtSet[ 1]  =   m_RhoExp0 /m_RhoCrud;    //!!! Interference ON
     m_WtSet[ 2]  =   m_RhoExp1 /m_RhoCrud;
     m_WtSet[ 3]  =   m_RhoExp2 /m_RhoCrud;
     m_WtBest     =   m_WtSet[ 3];
//[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
//     (*m_Out)<< ">>>KKceex::Make:  m_WtSet[ 1,2,3] ="<< m_WtSet[ 1]<<"  "<< m_WtSet[ 2]<<"  "<< m_WtSet[ 3]<<endl;
//]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
  }
}// Make

void KKceex::Born(int KFini, int KFfin, TLorentzVector &PX, double CosThetD,
                 KKpart &p1, KKpart &p2, KKpart &p3, KKpart &p4, KKcmplx4 &AmpBorn){
/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
//   Used in construction of the hard Non-IR parts: in GPS_HiniPlus, GPS_HfinPlus  //
//   CAREFUL!!! p_i are sometimes be substituted for photons!!!                    //
//   Input:                                                                        //
//   KFi, Kff = beam and final fermion flavour codes (to define charges)           //
//   PX       = s-chanel momentum for gamma and Z propagators (not for spinors)    //
//   pi,    are for spinors,  not for gamma and Z propagators                      //
//   p1,p2    =fermion momentum and mass (beam) <-- m_r1, m_r2                     //
//   p3,p4    =fermion momentum and mass final state                               //
//   Output:                                                                       //
//   AmpBorn   = spin amplitudes                                                   //
//   Notes:                                                                        //
//   Electron mass neglected in spinors, this is why we may use Chisholm!          //
//   Final fermion mass kept exactly.                                              //
//   Gamma and Z in s-channel.                                                      //
//                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////
  dcmplx TT, UU;
  int Hel1,Hel2,Hel3,Hel4;
  dcmplx s31,s24,s32,s14;
  for(int j1 = 0; j1<=1; j1++)
     for(int j2 = 0; j2<=1; j2++){
        for(int j3 = 0; j3<=1; j3++){
           for(int j4 = 0; j4<=1; j4++){
              Hel1 = 1-2*j1;
              Hel2 = 1-2*j2;
              Hel3 = 1-2*j3;
              Hel4 = 1-2*j4;
              TT  = dcmplx(0.0,0.0);
              UU  = dcmplx(0.0,0.0);
              if( Hel2 == -Hel1) { //!!! <--helicity conservation imposed
                 s31 = iProd2( 1,   Hel3, p3,  1,    Hel1, p1);
                 s24 = iProd2(-1,   Hel2, p2, -1,    Hel4, p4);
                 s32 = iProd2( 1,   Hel3, p3,  1,    Hel2, p2);
                 s14 = iProd2(-1,   Hel1, p1, -1,    Hel4, p4);
                 TT  = s31*s24;
                 UU  = s32*s14;
              }//if
              m_SpinoTT[j1][j2][j3][j4] =  TT;
              m_SpinoUU[j1][j2][j3][j4] =  UU;
           }//j4,
        }//j3
  }// j1,j2
//[[[[[[[[[[[[[[[[[[[[[[[[[[[*debug*
//  (*m_Out)<< "----------------------------------------KKceex::Born----------------------------------------------------"<<endl;
// for(int j1=0; j1<=1; j1++) for(int j2=0; j2<=1; j2++){ (*m_Out)<< "m_SpinoTT("<<j1<<","<<j2<<",*,*)=";
//     for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<m_SpinoTT[j1][j2][j][k];(*m_Out)<<endl;}
//  for(int j1=0; j1<=1; j1++) for(int j2=0; j2<=1; j2++){ (*m_Out)<< "m_SpinoUU("<<j1<<","<<j2<<",*,*)=";
//     for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<m_SpinoUU[j1][j2][j][k];(*m_Out)<<endl;}
//]]]]]]]]]]]]]]]]]]]]]]]]]]]

  double SvarX= PX*PX;
  double MZ    = m_DZ->m_MZ;    // the same as DB->_MZ i.e. xpar
  double GammZ = m_DZ->D_GamZ;  // from Dizet !!!
  double Qe  = DB->Qf[ KFini];
  double Qf  = DB->Qf[ KFfin];

// Possibility to switch off Z or gamma, etc.
  dcmplx PropGam, PropZet;
  if(DB->KeyZet == 9) {
     PropGam =  dcmplx(0.0);
  } else {
     PropGam =  dcmplx(  1.0/SvarX,  0.0);
  }//if
  if(DB->KeyZet == 0) {
     PropZet =  dcmplx(0.0);
  }else if((DB->KeyZet == -1) || (DB->KeyZet == -2) ) {
// for KKMC-hh fixed width option:
     PropZet =  1.0/dcmplx(SvarX- MZ*MZ, GammZ*MZ);
  } else {
// normal running width:
     PropZet =  1.0/dcmplx(SvarX-MZ*MZ, GammZ*SvarX/MZ);
  }// if else
  if( DB->KeyZet == -1) {  // using rescaled parameters
     PropZet = PropZet/dcmplx(1.0,GammZ/MZ);
  }//
// Exponentiate Resonance BigLogs according to Greco et al.
  if( m_KeyInt == 2 && DB->KeyISR != 0 &&  m_HasFSR != 0  ) {
       PropZet = PropZet * exp(m_IntReson);
  }

//================================================================
  dcmplx Ve, Vf, Ae, Af ;
  dcmplx GamVPi = dcmplx(1.0);
  dcmplx ZetVPi = dcmplx(1.0);
  dcmplx VVcor  = dcmplx(1.0);
  double RsqV=1.0;
  double RsqA=1.0;
  if( DB->KeyElw <= 0){
     EWFFact(KFini, KFfin, SvarX, Ve, Vf, Ae, Af ); // Born version!!!
  } else {
     EWFFact(KFini, KFfin, SvarX, CosThetD, Ve, Vf, Ae, Af,
                    VVcor, GamVPi, ZetVPi, RsqV, RsqA);
  }//

////////////////////////////////////////////////////////////////////////////////////////////
//     Primitives formfactor-type for construction of spin amplitudes                     //
//     (Ve -Hel1*Ae)*(Vf +Hel1*Af) is expanded because of double-vector f-factor          //
////////////////////////////////////////////////////////////////////////////////////////////
//  dcmplx FFacTT[2],FFacUU[2],FFacTG[2],FFacTZ[2],FFacUG[2],FFacUZ[2];
  if(DB->KeyElw <= 0){
    for(int j1 = 0; j1<=1; j1++){
      Hel1 = 1-2*j1;
      m_FFacTT[j1] = PropGam *dcmplx(Qe*Qf)  +PropZet *(Ve*Vf -dcmplx(Hel1)*Ae*Vf +dcmplx(Hel1)*Ve*Af -Ae*Af);
      m_FFacUU[j1] = PropGam *dcmplx(Qe*Qf)  +PropZet *(Ve*Vf -dcmplx(Hel1)*Ae*Vf -dcmplx(Hel1)*Ve*Af +Ae*Af);
  //      m_FFacTT[j1] = PropGam *Qe*Qf  +PropZet *(Ve- Hel1*Ae)*(Vf+ Hel1*Af);
  //      m_FFacUU[j1] = PropGam *Qe*Qf  +PropZet *(Ve- Hel1*Ae)*(Vf- Hel1*Af);
  //      cout<<"***  m_FFacTT  m_FFacUU ="<< m_FFacTT[j1]<<"   "<<m_FFacUU[j1]<<endl;
    }//for j1
  } else {
    for(int j1 = 0; j1<=1; j1++){
      Hel1 = 1-2*j1;
      m_FFacTG[j1] = PropGam*GamVPi *Qe*Qf *RsqV;
      m_FFacTZ[j1] = PropZet*ZetVPi *(Ve*Vf*VVcor*RsqV -dcmplx(Hel1)*Ae*Vf*RsqV +dcmplx(Hel1)*Ve*Af*RsqA -Ae*Af*RsqA);
      m_FFacUG[j1] = PropGam*GamVPi *Qe*Qf *RsqV;
      m_FFacUZ[j1] = PropZet*ZetVPi *(Ve*Vf*VVcor*RsqV -dcmplx(Hel1)*Ae*Vf*RsqV -dcmplx(Hel1)*Ve*Af*RsqA +Ae*Af*RsqA);
      m_FFacTT[j1] = m_FFacTG[j1]+m_FFacTZ[j1];
      m_FFacUU[j1] = m_FFacUG[j1]+m_FFacUZ[j1];
//      (*m_Out)<<"*** m_FFacTZ[j1]="<< m_FFacTZ[j1]<<"    m_FFacUZ[j1]="<< m_FFacUZ[j1]<<endl;
    }//for
  }//if

//[[[[[[[[[[[[[[[[[[[[[[[[
//  (*m_Out)<< "----------------------------------------KKceex::Born----------------------------------------------------"<<endl;
//      (*m_Out)<<"********************************************************************************"<<endl;
//      (*m_Out)<<"***    SvarX=" << SvarX  <<"  CosThetD="<< CosThetD <<endl;
//      (*m_Out)<<"***    PropGam="<< PropGam<<"    PropZet="<< PropZet<<endl;
//      (*m_Out)<<"***    Ve="<< Ve<<"  Vf="<< Vf<<"  Ae="<< Ae<<" Af="<< Af<<endl;
//      (*m_Out)<<"***    VVcor="<< VVcor<<endl;
//      (*m_Out)<<"***    GamVPi="<< GamVPi<<"  ZetVPi="<< ZetVPi<<"  RsqV="<< RsqV<<" RsqA="<< RsqA<<endl;
//      (*m_Out)<<"********************************************************************************"<<endl;
//]]]]]]]]]]]]]]]]]]]]]]]]

//[[[[[[[[[[[[[[[[[[[[[[[[[[[*debug*
//  (*m_Out)<< "----------------------------------------KKceex::Born----------------------------------------------------"<<endl;
//  (*m_Out)<< "m_FFacTT="; for(int j=0;j<=1;j++)(*m_Out)<<"  "<<SW208<<m_FFacTT[j];(*m_Out)<<endl;
//  (*m_Out)<< "m_FFacUU="; for(int j=0;j<=1;j++)(*m_Out)<<"  "<<SW208<<m_FFacUU[j];(*m_Out)<<endl;
//]]]]]]]]]]]]]]]]]]]]]]]]]]]

/////////////////////////////////////////////////////
// Dresed Born 16 amplitudes
  for(int j1 = 0; j1<=1; j1++)
    for(int j2 = 0; j2<=1; j2++){
      for(int j3 = 0; j3<=1; j3++){
        for(int j4 = 0; j4<=1; j4++){
            AmpBorn.m_A[j1][j2][j3][j4] =
                   m_SpinoTT[j1][j2][j3][j4]* m_FFacTT[j1]
                  +m_SpinoUU[j1][j2][j3][j4]* m_FFacUU[j1];
        }//j4
      }//j3
//[[[[[[[[[[[[[[[[[[[[[[[[[[[
//     (*m_Out)<< "AmpBorn="; for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<AmpBorn.m_A[j1][j2][j][k];(*m_Out)<<endl;
//]]]]]]]]]]]]]]]]]]]]]]]]]]]
    }//j1,j2


//[[[[[[[[[[[[[[[[[[[[[[[[[[[*debug*
//  for(int j1=0; j1<=1; j1++) for(int j2=0; j2<=1; j2++){ (*m_Out)<< "AmpBorn("<<j1<<","<<j2<<",*,*)=";
//     for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<AmpBorn.m_A[j1][j2][j][k];(*m_Out)<<endl;}
//]]]]]]]]]]]]]]]]]]]]]]]]]]]

}//KKceex::Born


void KKceex::BornW(int KFini, int KFfin, TLorentzVector &PX, double s, double t,
                 KKpart &p1, KKpart &p2, KKpart &p3, KKpart &p4, KKcmplx4 &AmpBornW){
/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
//   Used in construction of the hard Non-IR parts: in GPS_HiniPlus, GPS_HfinPlus  //
//   CAREFUL!!! p_i are sometimes be substituted for photons!!!                    //
//   Input:                                                                        //
//   KFini, Kffin = beam and final fermion flavour codes (to define charges)       //
//   PX       = s-chanel momentum for gamma and Z propagators (not for spinors)    //
//   pi,    are for spinors,  not for gamma and Z propagators                      //
//   p1,p2    =fermion momentum and mass (beam) <-- m_r1, m_r2                     //
//   p3,p4    =fermion momentum and mass final state                               //
//   Output:                                                                       //
//   AmpBorn   = spin amplitudes                                                   //
//   Notes:                                                                        //
//   Electron mass neglected in spinors, this is why we may use Chisholm!          //
//   Final fermion mass kept exactly.                                              //
//   Gamma and Z in s-channel.                                                      //
//                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////
Amp4Zer(AmpBornW);
if(abs(KFfin) != 12)  return;
dcmplx TT, UU;
int Hel1,Hel2,Hel3,Hel4;
dcmplx s31,s24,s32,s14;
//[[[[[[[[[[[[[[[[[[[[[[[[[[[
//(*m_Out)<< "----------------------------------------KKceex::BornW----------------------------------------------------"<<endl;
//]]]]]]]]]]]]]]]]]]]]]]]]]]]
for(int j1 = 0; j1<=1; j1++)
   for(int j2 = 0; j2<=1; j2++){
      for(int j3 = 0; j3<=1; j3++){
         for(int j4 = 0; j4<=1; j4++){
            Hel1 = 1-2*j1;
            Hel2 = 1-2*j2;
            Hel3 = 1-2*j3;
            Hel4 = 1-2*j4;
            TT  = dcmplx(0.0,0.0);
            UU  = dcmplx(0.0,0.0);
            if( Hel2 == -Hel1) { //!!! <--helicity conservation imposed
               s31 = iProd2( 1,   Hel3, p3,  1,    Hel1, p1);
               s24 = iProd2(-1,   Hel2, p2, -1,    Hel4, p4);
               s32 = iProd2( 1,   Hel3, p3,  1,    Hel2, p2);
               s14 = iProd2(-1,   Hel1, p1, -1,    Hel4, p4);
               TT  = s31*s24;
               UU  = s32*s14;
            }//if
            m_SpinoTT[j1][j2][j3][j4] =  TT;
            m_SpinoUU[j1][j2][j3][j4] =  UU;
         }//j4,
      }//j3
   }//j1,j2
/////////////////////////
//double  s =  PX*PX;
//double  t = -s*(1.0-CosThetD)/2;
//double  u = -s*(1.0+CosThetD)/2;
////////////////////////////////////////////////////////////////////////////////////////////
double Sw2 = m_DZ->D_swsq;
double Coef  =1.0/2.0/Sw2;
dcmplx PropW,WVPi;
GPS_EWFFactW(KFini,KFfin,s,t,PropW,WVPi);
double Ve= 0.5;
double Vf= 0.5;
double Ae= 0.5;
double Af= 0.5;
dcmplx      FFacTT[2],   FFacUU[2];
for(int j1 =0; j1<=1;j1++){
   Hel1 = 1-2*j1;
   FFacTT[j1] = Coef* PropW*WVPi *(Ve*Vf -Hel1*Ae*Vf +Hel1*Ve*Af -Ae*Af);
   FFacUU[j1] = Coef* PropW*WVPi *(Ve*Vf -Hel1*Ae*Vf -Hel1*Ve*Af +Ae*Af);
}//j1
///////////////////////////////////////////////////////////////////////////////////
//                      Total result = Spinors*Formfactor                        //
///////////////////////////////////////////////////////////////////////////////////
for(int j1 = 0; j1<=1; j1++)
  for(int j2 = 0; j2<=1; j2++){
    for(int j3 = 0; j3<=1; j3++)
      for(int j4 = 0; j4<=1; j4++){
// Born,  Zero order
            AmpBornW.m_A[j1][j2][j3][j4]=
                       m_SpinoTT[j1][j2][j3][j4]* FFacTT[j1]
                      +m_SpinoUU[j1][j2][j3][j4]* FFacUU[j1];
      }//j3,j4
  }//j1,j2
//[[[[[[[[[[[[[[[[[[[[[[[[[[[*debug*
//(*m_Out)<<"*******************************KKceex::BornW****************************************"<<endl;
//(*m_Out)<< "FFacTT="; for(int j=0;j<=1;j++)(*m_Out)<<"  "<<SW208<<FFacTT[j];(*m_Out)<<endl;
//(*m_Out)<< "FFacUU="; for(int j=0;j<=1;j++)(*m_Out)<<"  "<<SW208<<FFacUU[j];(*m_Out)<<endl;
//(*m_Out)<< "s,t,PropW,WVPi="<< s<<"  "<<t<<"  "<<PropW<<"  "<<WVPi <<endl;
//  for(int j1=0; j1<=1; j1++) for(int j2=0; j2<=1; j2++){ (*m_Out)<< "AmpBornW("<<j1<<","<<j2<<",*,*)=";
//     for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<AmpBornW.m_A[j1][j2][j][k];(*m_Out)<<endl;}
//]]]]]]]]]]]]]]]]]]]]]]]]]]]
}//KKceex::BornW

double  KKceex::BornFoam0(int KFini, int KFfin, double SvarX, double CosThetD){
/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////
  TLorentzVector PX;
  PX.SetPxPyPzE( 0, 0, 0, sqrt(SvarX) );
//
  double MZ    = m_DZ->m_MZ;    // the same as DB->_MZ i.e. xpar
  double GammZ = m_DZ->D_GamZ;  // from Dizet !!!
  double Qe  = DB->Qf[ KFini];
  double Qf  = DB->Qf[ KFfin];
  //================================================================
  dcmplx m_BoxGGtu  = dcmplx(0.0,0.0);
  dcmplx m_BoxGGut  = dcmplx(0.0,0.0);
  dcmplx m_BoxGZtu  = dcmplx(0.0,0.0);
  dcmplx m_BoxGZut  = dcmplx(0.0,0.0);

  m_IntReson        = dcmplx(0.0,0.0);
  m_IntIR           = dcmplx(0.0,0.0);
//
  double  s =  SvarX;
  double  t = -s*(1.0-CosThetD)/2;
  double  u = -s*(1.0+CosThetD)/2;
if(  m_KeyInt != 0 && m_KeyISR != 0 &&  m_HasFSR != 0  ) {
  dcmplx  Coef  = dcmplx(m_Alfpi*Qe*Qf);
  m_IntReson = Coef*m_BVR->IntReson(DB->MasPhot, MZ, GammZ, s, t, u); //<- asymetric in (t,u)
  // Virtual 2*(B(t)-B(u)) Intereference IR part to be subtracted from boxes
  m_IntIR    = Coef*m_BVR->IntIR(   DB->MasPhot,s,t,u);               //<- asymetric in (t,u)

// Turn of Boxes in no-EWK xchecks
  if (DB->KeyElw > 0) {
    m_BoxGGtu  = Coef*( m_BVR->CBoxGG(DB->MasPhot,         s,t,u)) -m_IntIR;
    m_BoxGZtu  = Coef*( m_BVR->CBoxGZ(DB->MasPhot,MZ,GammZ,s,t,u)) -m_IntIR;
    m_BoxGGut  = Coef*(-m_BVR->CBoxGG(DB->MasPhot,         s,u,t)) -m_IntIR;
    m_BoxGZut  = Coef*(-m_BVR->CBoxGZ(DB->MasPhot,MZ,GammZ,s,u,t)) -m_IntIR;
   }//
// Exponentiate Resonance BigLogs according to Greco et al.
  if( m_KeyInt == 2) {
    m_BoxGZtu = m_BoxGZtu -m_IntReson;
    m_BoxGZut = m_BoxGZut -m_IntReson;
  }//if m_KeyInt
}// if
//===================================================================
  Born(KFini, KFfin, PX, CosThetD, m_r1, m_r2, m_p3, m_p4, m_BornC);
//===================================================================
  double Sum0 = 0;
  for(int j1 = 0; j1<=1; j1++)
      for(int j2 = 0; j2<=1; j2++)
         for(int j3 = 0; j3<=1; j3++)
            for(int j4 = 0; j4<=1; j4++){
                     Sum0 = Sum0 +real(m_BornC.m_A[j1][j2][j3][j4]*conj(m_BornC.m_A[j1][j2][j3][j4]));
            }
  double Massf = DB->fmass[ KFfin];
  double betaf = sqrt( 1.0 - 4* sqr(Massf)/SvarX );   //<-- 2-body phase spase
  Sum0 *=betaf;
  return Sum0;
}//BornFoam0


void  KKceex::BornFoam2(int Mode, int KFini, int KFfin,
                        double SvarX, double CosThetD, double &Yint){
/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////
  TLorentzVector PX;
  PX.SetPxPyPzE( 0, 0, 0, sqrt(SvarX) );
//
  double MZ    = m_DZ->m_MZ;    // the same as DB->_MZ i.e. xpar
  double GammZ = m_DZ->D_GamZ;  // from Dizet !!!
  double Qe  = DB->Qf[ KFini];
  double Qf  = DB->Qf[ KFfin];
  //================================================================
  dcmplx m_BoxGGtu  = dcmplx(0.0,0.0);
  dcmplx m_BoxGGut  = dcmplx(0.0,0.0);
  dcmplx m_BoxGZtu  = dcmplx(0.0,0.0);
  dcmplx m_BoxGZut  = dcmplx(0.0,0.0);

  m_IntReson        = dcmplx(0.0,0.0);
  m_IntIR           = dcmplx(0.0,0.0);
//
  double  s =  SvarX;
  double  t = -s*(1.0-CosThetD)/2;
  double  u = -s*(1.0+CosThetD)/2;
if(  m_KeyInt != 0 && m_KeyISR != 0 &&  m_HasFSR != 0  ) {
  dcmplx  Coef  = dcmplx(m_Alfpi*Qe*Qf);
  m_IntReson = Coef*m_BVR->IntReson(DB->MasPhot, MZ, GammZ, s, t, u); //<- asymetric in (t,u)
  // Virtual 2*(B(t)-B(u)) Intereference IR part to be subtracted from boxes
  m_IntIR    = Coef*m_BVR->IntIR(   DB->MasPhot,s,t,u);               //<- asymetric in (t,u)

// Gamma-Z boxes (finite part)
  if (DB->KeyElw > 0) {
    m_BoxGGtu  = Coef*( m_BVR->CBoxGG(DB->MasPhot,         s,t,u)) -m_IntIR;
    m_BoxGZtu  = Coef*( m_BVR->CBoxGZ(DB->MasPhot,MZ,GammZ,s,t,u)) -m_IntIR;
    m_BoxGGut  = Coef*(-m_BVR->CBoxGG(DB->MasPhot,         s,u,t)) -m_IntIR;
    m_BoxGZut  = Coef*(-m_BVR->CBoxGZ(DB->MasPhot,MZ,GammZ,s,u,t)) -m_IntIR;
   }//
// Resonance BigLogs exponentiated according to Greco et al.
  if( m_KeyInt == 2) {
    m_BoxGZtu = m_BoxGZtu -m_IntReson; // avoid double counting
    m_BoxGZut = m_BoxGZut -m_IntReson; // avoid double counting
  }//if m_KeyInt
}// if

Amp4Zer(m_Boxy);
if( m_KeyInt == 2) MakeAmpBox(); // making m_Boxy, non-IR part

double ChaIni=DB->Qf[KFini];
double ChaFin=DB->Qf[KFfin];
double alfpmix  = m_Alfpi*ChaFin*ChaIni;
double Ene = sqrt(SvarX)/2.0;
// Remember Yint depends on Emin and provides angular asymmetry (MasPhot is dummy)
Yint= TForFac( alfpmix, m_p1, m_p3, Ene, DB->MasPhot)
     *TForFac( alfpmix, m_p2, m_p4, Ene, DB->MasPhot)
     *TForFac(-alfpmix, m_p1, m_p4, Ene, DB->MasPhot)
     *TForFac(-alfpmix, m_p2, m_p3, Ene, DB->MasPhot);
//[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
//m_icont++;
//if(m_icont<1000) {
//	cout<<"===================================="<<endl;
//	cout<<"m_p1= "; m_p1.Print();cout<<endl;
//	cout<<"m_p2= "; m_p2.Print();cout<<endl;
//	cout<<"m_p3= "; m_p3.Print();cout<<endl;
//	cout<<"m_p4= "; m_p4.Print();cout<<endl;
//	cout<<" BornFoam2: m_Emin ="<<m_Emin<<endl;
//	cout<<" BornFoam2:   Yint ="<<Yint<<endl;
//}
//]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
//===================================================================
dcmplx cOne = dcmplx(1.0, 0.0);
  if(        Mode == 10){
    Born(KFini, KFfin, PX, CosThetD, m_r1, m_r2, m_p3, m_p4, m_BornB);
//    AmpAdd( m_BornB, cOne,m_Boxy); // adding finite part of Boxes
  } else if( Mode == 11){
    Born(KFini, KFfin, PX, CosThetD, m_r1, m_r2, m_p3, m_p4, m_BornC);
//    AmpAdd( m_BornC, cOne,m_Boxy); // adding finite part of Boxes
  } else {
  cout<<"++++ KKceex::BornFoam2: wrong Mode ="<<Mode<<endl;
  }
//===================================================================
}//BornFoam2

//_____________________________________________________________________
void KKceex::MakeAmpBox(){
// Defining Gamma-Z Box amplitudes
dcmplx BoxGG, BoxGZ, AmpBoxy, AmpBorn;
int Hel1,Hel2,Hel3,Hel4;
for(int j1 = 0; j1<=1; j1++)
  for(int j2 = 0; j2<=1; j2++){
    for(int j3 = 0; j3<=1; j3++)
      for(int j4 = 0; j4<=1; j4++){
        Hel1 = 1-2*j1;
        Hel2 = 1-2*j2;
        Hel3 = 1-2*j3;
        Hel4 = 1-2*j4;
//  O(alf1) QED Boxes
        BoxGG = dcmplx(0.0,0.0); BoxGZ= dcmplx(0.0,0.0);
        if((Hel2 == -Hel1) && (Hel4 == -Hel3)) { //<--helicity conserv.
           if( Hel1*Hel3 == 1) {
              BoxGG = m_BoxGGtu;
              BoxGZ = m_BoxGZtu;
           } else {
              BoxGG = m_BoxGGut;
              BoxGZ = m_BoxGZut;
           }//if
        }//if
        AmpBoxy = dcmplx(0.0,0.0);
        if( DB->KeyElw > 0) {
           AmpBoxy =  m_SpinoTT[j1][j2][j3][j4]* m_FFacTG[j1] *BoxGG
                     +m_SpinoTT[j1][j2][j3][j4]* m_FFacTZ[j1] *BoxGZ
                     +m_SpinoUU[j1][j2][j3][j4]* m_FFacUG[j1] *BoxGG
                     +m_SpinoUU[j1][j2][j3][j4]* m_FFacUZ[j1] *BoxGZ;
        }
  //
        m_Boxy.m_A[j1][j2][j3][j4] = AmpBoxy;   // O(alf1) Gamma-Z boxes
     }//for j3,j4
}//for j1,j2
//---------------------
}//MakeBox


void KKceex::BornPlus(int KFini, int KFfin, dcmplx Cfac, TLorentzVector &PX){
/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
//   Two-stroke version of GPS_Born, optimized for summation over partitions       //
//   Virtual corrections (boxes vertices) to ms-one are also here!!!               //
//   Warning! Input masses are TRUE (not signed as in GPS_Born)                    //
//                                                                                 //
//   Born spin amplitudes calculated with spinor methods.                          //
//   Mass of the final fermion kept exactly.                                       //
//                                                                                 //
//   Input:                                                                        //
//   KFi, Kff = beam and final fermion flavour codes (to define charges)           //
//   PX       = s-chanel momentum for gamma and Z propagators (not for spinors)    //
//                                                                                 //
//   Common working space:                                                         //
//   m_AmpExpo*  is working space, used by HiniPlus, HfinPlus, HfinMinus           //
//                                                                                 //
//   Notes:                                                                        //
//   Electron mass neglected in spinors, this is why we may use Chisholm!          //
//   Final fermion mass kept exactly.                                              //
//   Gamma and Z in s-chanel.                                                      //
//                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////
  int Hel1,Hel2,Hel3,Hel4;

  double SvarX= PX*PX;
  double MZ    = m_DZ->m_MZ;    // the same as DB->_MZ i.e. xpar
  double GammZ = m_DZ->D_GamZ;  // from Dizet !!!
  double Qe  = DB->Qf[ KFini];
  double Qf  = DB->Qf[ KFfin];

// IR-subtracted, helicity conservation imposed (small mass approx.)
  if(SvarX <= sqr( m_p3.M + m_p4.M) ) return; //????
  double CosThetD;
  m_Event->ThetaD(PX,CosThetD);

//================================================================
  dcmplx m_BoxGGtu  = dcmplx(0.0,0.0);
  dcmplx m_BoxGGut  = dcmplx(0.0,0.0);
  dcmplx m_BoxGZtu  = dcmplx(0.0,0.0);
  dcmplx m_BoxGZut  = dcmplx(0.0,0.0);

  m_IntReson        = dcmplx(0.0,0.0);
  m_IntIR           = dcmplx(0.0,0.0);

  double  s =  SvarX;
  double  t = -s*(1.0-CosThetD)/2;
  double  u = -s*(1.0+CosThetD)/2;
  if(  m_KeyInt != 0 && m_KeyISR != 0 &&  m_HasFSR != 0  ) {
	dcmplx  Coef  = dcmplx(m_Alfpi*Qe*Qf);
    m_IntReson = Coef*m_BVR->IntReson(DB->MasPhot, MZ, GammZ, s, t, u); //<- asymetric in (t,u)
// Virtual 2*(B(t)-B(u)) Intereference IR part to be subtracted from boxes
    m_IntIR      = Coef*m_BVR->IntIR(   DB->MasPhot,s,t,u);               //<- asymetric in (t,u)

// Turn of Boxes in no-EWK xchecks
    if (DB->KeyElw > 0) {
       m_BoxGGtu  = Coef*( m_BVR->CBoxGG(DB->MasPhot,         s,t,u)) -m_IntIR;
       m_BoxGZtu  = Coef*( m_BVR->CBoxGZ(DB->MasPhot,MZ,GammZ,s,t,u)) -m_IntIR;
       m_BoxGGut  = Coef*(-m_BVR->CBoxGG(DB->MasPhot,         s,u,t)) -m_IntIR;
       m_BoxGZut  = Coef*(-m_BVR->CBoxGZ(DB->MasPhot,MZ,GammZ,s,u,t)) -m_IntIR;
    }//
// Exponentiate Resonance BigLogs according to Greco et al.
    if( m_KeyInt == 2) {
      m_BoxGZtu = m_BoxGZtu -m_IntReson;
      m_BoxGZut = m_BoxGZut -m_IntReson;
    }//if m_KeyInt
  }
//===================================================================
  Born(KFini, KFfin, PX, CosThetD, m_r1, m_r2, m_p3, m_p4, m_BornC);
//===================================================================

///////////////////////////////////////////////////////////////////////////////////
//       QED vertex  FFactor F1 minus B-vrtual (IR removed), exact mass terms    //
///////////////////////////////////////////////////////////////////////////////////
    TLorentzVector PP = m_Event->m_Pf1 + m_Event->m_Pf2;
    TLorentzVector QQ = m_Event->m_Qf1 + m_Event->m_Qf2;
    double SvarP = PP*PP;
    double SvarQ = QQ*QQ;
    m_F1ini1 = dcmplx(0.0);
    m_F1fin1 = dcmplx(0,0);
    m_F1ini2 = dcmplx(0,0);
    m_F1fin2 = dcmplx(0,0);
    if( m_KeyISR != 0) {
      MakeF1ini(SvarP, m_p1.M, m_p2.M, DB->Qf[KFini], m_F1ini1, m_F1ini2);
    }
    if( m_HasFSR != 0) {
      MakeF1fin(SvarQ, m_p3.M, m_p4.M, DB->Qf[KFfin], m_F1fin1, m_F1fin2);
    }
//---------------------------------------------------------------------------
// ZBW: Where is dominant other photon?  It is instead of reduction procedure
// IF ((p3(3)+p4(3))*p1(3).LT.0D0) THEN
double s0, t0;
if( (m_p3.P[3]+m_p4.P[3])*m_p1.P[3] < 0.0){
          s0=  2*(m_p3*m_p4); t0= -2*(m_p4*m_p2);
} else {  s0=  2*(m_p3*m_p4); t0= -2*(m_p3*m_p1);}//if
//===================================================================
BornW(KFini, KFfin, PX, s0, t0, m_r1, m_r2, m_p3, m_p4, m_AmpBornW);
//===================================================================

dcmplx BoxGG, BoxGZ, AmpBoxy, AmpBorn,AmpBornW;
for(int j1 = 0; j1<=1; j1++)
  for(int j2 = 0; j2<=1; j2++){
    for(int j3 = 0; j3<=1; j3++)
      for(int j4 = 0; j4<=1; j4++){
        Hel1 = 1-2*j1;
        Hel2 = 1-2*j2;
        Hel3 = 1-2*j3;
        Hel4 = 1-2*j4;
//  O(alf1) QED Boxes
        BoxGG = dcmplx(0.0,0.0); BoxGZ= dcmplx(0.0,0.0);
        if((Hel2 == -Hel1) && (Hel4 == -Hel3)) { //<--helicity conserv.
           if( Hel1*Hel3 == 1) {
              BoxGG = m_BoxGGtu;
              BoxGZ = m_BoxGZtu;
           } else {
              BoxGG = m_BoxGGut;
              BoxGZ = m_BoxGZut;
           }//if
        }//if
        AmpBoxy = dcmplx(0.0,0.0);
        if( DB->KeyElw > 0) {
           AmpBoxy =  m_SpinoTT[j1][j2][j3][j4]* m_FFacTG[j1] *BoxGG
                     +m_SpinoTT[j1][j2][j3][j4]* m_FFacTZ[j1] *BoxGZ
                     +m_SpinoUU[j1][j2][j3][j4]* m_FFacUG[j1] *BoxGG
                     +m_SpinoUU[j1][j2][j3][j4]* m_FFacUZ[j1] *BoxGZ;
        }
  // Accumulate and assemble results
        AmpBorn  = m_BornC.m_A[j1][j2][j3][j4];
        AmpBornW = m_AmpBornW.m_A[j1][j2][j3][j4];
        AmpBorn += AmpBornW;
        m_BornD.m_A[j1][j2][j3][j4]=AmpBorn;                     // W included
        m_AmpExpo0.m_A[j1][j2][j3][j4] +=  Cfac*AmpBorn;         // O(alf^0)_exp
  //
        m_AmpExpo1.m_A[j1][j2][j3][j4] +=
                 +Cfac*AmpBorn*(1.0 +m_F1ini1)*(1.0 +m_F1fin1 )  // Born, O(alf1) m_FFactors
                 +Cfac*AmpBoxy;                                  // O(alf1) boxes
//[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
//[[[[[test        m_AmpExpo1.m_A[j1][j2][j3][j4] = dcmplx(0.0);
//]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
  //
        m_AmpExpo2.m_A[j1][j2][j3][j4] +=
                 +Cfac*AmpBorn*(1.0 +m_F1ini2)*(1.0 +m_F1fin2 )  // Born, O(alf2) FFactors
                 +Cfac*AmpBoxy;                                  // O(alf1) boxes
     }//for j3,j4
  }//for j1,j2

//[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
//dcmplx PropW, WVPi;
//GPS_EWFFactW(KFini, KFfin, s,t, PropW, WVPi);
//]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]

//[[[[[[[[[[[[[[[[[[[[[[[[[[[*debug*
//(*m_Out)<<"////////////////////////////////////////BornPlus//////////////////////////////////////////////////"<<endl;
//(*m_Out)<< "m_FFacTT="; for(int j=0;j<=1;j++)(*m_Out)<<"  "<<SW208<<m_FFacTT[j];(*m_Out)<<endl;
//(*m_Out)<< "m_FFacUU="; for(int j=0;j<=1;j++)(*m_Out)<<"  "<<SW208<<m_FFacUU[j];(*m_Out)<<endl;
//(*m_Out)<<"-----------------------------------------------------------------------------------------------"<<endl;
//for(int j1=0; j1<=1; j1++) for(int j2=0; j2<=1; j2++){(*m_Out)<< "   m_BornC("<<j1<<","<<j2<<",*,*)=";
//   for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<m_BornC.m_A[j1][j2][j][k];(*m_Out)<<endl;}
//(*m_Out)<<"-----------------------------------------------------------------------------------------------"<<endl;
//for(int j1=0; j1<=1; j1++) for(int j2=0; j2<=1; j2++){(*m_Out)<< "m_AmpBornW("<<j1<<","<<j2<<",*,*)=";
//   for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<m_AmpBornW.m_A[j1][j2][j][k];(*m_Out)<<endl;}
//]]]]]]]]]]]]]]]]]]]]]]]]]]]


}//BornPlus

void KKceex::ZerAmplit(){
//
for(int j1 = 0; j1<=1; j1++)
  for(int j2 = 0; j2<=1; j2++)
    for(int j3 = 0; j3<=1; j3++)
      for(int j4 = 0; j4<=1; j4++){
            m_AmpExpo0.m_A[j1][j2][j3][j4] =  dcmplx(0.0,0.0);
            m_AmpExpo1.m_A[j1][j2][j3][j4] =  dcmplx(0.0,0.0);
            m_AmpExpo2.m_A[j1][j2][j3][j4] =  dcmplx(0.0,0.0);
      }
//----------------------
}//ZerAmplit


void KKceex::MakeRho(){
//////////////////////////////////////////////////////////////////////////////////
//                                                                              //
//   Calculate differential distributions (normalized to LIPS) from spin ampl.  //
//                UNPOLARIZED final fermions                                    //
//                                                                              //
//   To be done:                                                                //
//   One needs to put Wigner rotation for initial polarizations somewhere       //
//   either in setter called in KK2f, or in flight, for every event             //
//                                                                              //
//////////////////////////////////////////////////////////////////////////////////

double  Sum0 = 0;
double  Sum1 = 0;
double  Sum2 = 0;
//    IF( (polar1+polar2) .LT. 1d-6) THEN
//        The case with UNPOLARIZED beams and UNPOLARIZED final fermions
for(int j1 = 0; j1<=1; j1++)
    for(int j2 = 0; j2<=1; j2++)
       for(int j3 = 0; j3<=1; j3++)
          for(int j4 = 0; j4<=1; j4++){
                   Sum0 = Sum0 +real(m_AmpExpo0.m_A[j1][j2][j3][j4]*conj(m_AmpExpo0.m_A[j1][j2][j3][j4]));
                   Sum1 = Sum1 +real(m_AmpExpo1.m_A[j1][j2][j3][j4]*conj(m_AmpExpo1.m_A[j1][j2][j3][j4]));
                   Sum2 = Sum2 +real(m_AmpExpo2.m_A[j1][j2][j3][j4]*conj(m_AmpExpo2.m_A[j1][j2][j3][j4]));
          }

m_RhoExp0 = Sum0 *m_ExpoNorm;
m_RhoExp1 = Sum1 *m_ExpoNorm;
m_RhoExp2 = Sum2 *m_ExpoNorm;
//[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
(*m_Out)<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@KKceex::MakeRho@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endl;
//(*m_Out)<<"m_ExpoNorm= "<<m_ExpoNorm<<endl;
//  for(int j1=0; j1<=1; j1++) for(int j2=0; j2<=1; j2++){ (*m_Out)<< "m_AmpExpo0("<<j1<<","<<j2<<",*,*)=";
//     for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<m_AmpExpo0.m_A[j1][j2][j][k];(*m_Out)<<endl;}
//(*m_Out)<<"-----------------------------------------------------------------------------------------------"<<endl;
//  for(int j1=0; j1<=1; j1++) for(int j2=0; j2<=1; j2++){ (*m_Out)<< "m_AmpExpo1("<<j1<<","<<j2<<",*,*)=";
//     for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<m_AmpExpo1.m_A[j1][j2][j][k];(*m_Out)<<endl;}
//  (*m_Out)<<"-----------------------------------------------------------------------------------------------"<<endl;
//    for(int j1=0; j1<=1; j1++) for(int j2=0; j2<=1; j2++){ (*m_Out)<< "m_AmpExpo2("<<j1<<","<<j2<<",*,*)=";
//       for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<m_AmpExpo2.m_A[j1][j2][j][k];(*m_Out)<<endl;}
(*m_Out)<<"@@@@@@@@ m_RhoExp0= "<<m_RhoExp0<<"  m_RhoExp1= "<<m_RhoExp1<<"  m_RhoExp2= "<<m_RhoExp2<<endl;
//]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
}//MakeRho


double KKceex::MakeRhoFoam(){
//////////////////////////////////////////////////////////////////////////////////
//                                                                              //
//   Calculate differential distributions do be used in KKhhFoam                //
//                                                                              //
//////////////////////////////////////////////////////////////////////////////////
double  Sum0 = 0;
for(int j1 = 0; j1<=1; j1++)
    for(int j2 = 0; j2<=1; j2++)
       for(int j3 = 0; j3<=1; j3++)
          for(int j4 = 0; j4<=1; j4++){
                   Sum0 = Sum0 +  real(m_BornB.m_A[j1][j2][j3][j4]
                                 *conj(m_BornC.m_A[j1][j2][j3][j4]));
          }
return Sum0;
}//MakeRho




void KKceex::MakeRho2(double h1[], double h2[], double &wt0, double &wt1, double &wt2){
//////////////////////////////////////////////////////////////////////////////////
//                                                                              //
//   Used in Taupair_ImprintSpin                                                //
//                                                                              //
//   Calculate differential distributions (normalized to LIPS) from spin ampl.  //
//                  POLARIZED final fermions                                    //
//                                                                              //
//                                                                              //
//////////////////////////////////////////////////////////////////////////////////
// translate indexing in h-vectors, time component has to be 0-th
  double HvecFer1[4], HvecFer2[4];
  HvecFer1[0] = h1[3];
  for(int k=1; k<=3; k++) HvecFer1[k] = h1[k-1];
  HvecFer2[0] = h2[3];
  for(int k=1; k<=3; k++) HvecFer2[k] = h2[k-1];
// Define polarimeter density matrices
  for(int i=0; i<=1; i++)
     for(int j=0; j<=1; j++){
        m_SDMat3[i][j]=dcmplx(0.0);
        m_SDMat4[i][j]=dcmplx(0.0);
        for(int k=0; k<=3; k++){
           m_SDMat3[i][j]=m_SDMat3[i][j]+m_Pauli[ k][i][j] *HvecFer1[k];
           m_SDMat4[i][j]=m_SDMat4[i][j]+m_Pauli[ k][i][j] *HvecFer2[k];
        }// for k
     }// for i,j
//  The case with UNPOLARIZED beams and POLARIZED final fermions
  double Sum0 = 0.0;
  double Sum1 = 0.0;
  double Sum2 = 0.0;
  dcmplx Tensor0, Tensor1, Tensor2, SDMprod;
  for(int j1=0; j1<=1; j1++)
    for(int j2=0; j2<=1; j2++)
        for(int j3=0; j3<=1; j3++)
          for(int i3=0; i3<=1; i3++)
            for(int j4=0; j4<=1; j4++)
              for(int i4=0; i4<=1; i4++){
                SDMprod = m_SDMat3[j3][i3]*m_SDMat4[j4][i4];
                Tensor0 = m_AmpExpo0.m_A[j1][j2][i3][i4]*conj(m_AmpExpo0.m_A[j1][j2][j3][j4]);
                Tensor1 = m_AmpExpo1.m_A[j1][j2][i3][i4]*conj(m_AmpExpo1.m_A[j1][j2][j3][j4]);
                Tensor2 = m_AmpExpo2.m_A[j1][j2][i3][i4]*conj(m_AmpExpo2.m_A[j1][j2][j3][j4]);
                Sum0 = Sum0 + real(Tensor0 *SDMprod);
                Sum1 = Sum1 + real(Tensor1 *SDMprod);
                Sum2 = Sum2 + real(Tensor2 *SDMprod);
              }
// distributions with polarized final fermions
  double Rho0 = Sum0 *m_ExpoNorm;
  double Rho1 = Sum1 *m_ExpoNorm;
  double Rho2 = Sum2 *m_ExpoNorm;
// Spin weight for  polarized final fermions
  wt0 = Rho0/m_RhoExp0;
  wt1 = Rho1/m_RhoExp1;
  wt2 = Rho2/m_RhoExp2;
  //[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
  (*m_Out)<<"----------------------------------KKceex::MakeRho2---------------------------------------------------"<<endl;
  (*m_Out)<<" h1[0-3]="<<h1[0]<<" "<<h1[1]<<" "<<h1[2]<<" "<<h1[3]<<endl;
  (*m_Out)<<" h2[0-3]="<<h2[0]<<" "<<h2[1]<<" "<<h2[2]<<" "<<h2[3]<<endl;
//  (*m_Out)<<"HvecFer1="<<HvecFer1[0]<<" "<<HvecFer1[1]<<" "<<HvecFer1[2]<<" "<<HvecFer1[3]<<endl;
//  (*m_Out)<<"HvecFer2="<<HvecFer2[0]<<" "<<HvecFer2[1]<<" "<<HvecFer2[2]<<" "<<HvecFer2[3]<<endl;
  (*m_Out)<<" Rho0= "<<Rho0<<"  Rho1= "<<Rho1<<"  Rho2= "<<Rho2<<endl;
  //]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
}// GPS_MakeRho2



void KKceex::HiniPlus(int KFini, int KFfin, TLorentzVector &PX,
               KKpart &ph1, int Hel, dcmplx &Sactu, dcmplx &sProd){
////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
//   IR-finite part od 1-photon amplitudes for ISR  (equiv. to GPS_Hini)           //
//   Photon helicity imported from the calling program                             //
//                                                                                 //
//   m_AmpExpo*  is working space                                                  //
//   m_BornC.is hidden INPUT                                                       //
//                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////
double Qe  = DB->Qf[ KFini];
dcmplx Vir1,Vir2;  // Virtual corrections
MakeVini(m_p1, m_p2, ph1, Vir1,Vir2);
// ISR non-infrared two parts: (1) p1 -> photon, contracted with U-matrix
//                             (2) p2 -> photon, contracted with V-matrix
KKcmplx4 AmpBornU, AmpBornV; // dressed Born spin amplitudes
//* Calculate Born spin amplitudes
//      CALL GPS_Born     (    KFi,KFf,PX,       ph,mph,    p2,-Fleps,  p3,m3,   p4,-m4,   AmpBornU)
//      CALL GPS_Born(         KFi,KFf,PX,       p1,Fleps,  ph,-mph,    p3,m3,   p4,-m4,   AmpBornV)
double CosThe=0.0;
KKpart ph1r = ph1; ph1r.C=-1;
Born(KFini, KFfin, PX, CosThe,  ph1, m_r2, m_p3, m_p4, AmpBornU);
Born(KFini, KFfin, PX, CosThe, m_r1, ph1r, m_p3, m_p4, AmpBornV);
//[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
//(*m_Out)<< "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&KKceex::HiniPlus&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<endl;
//for(int j1=0; j1<=1; j1++) for(int j2=0; j2<=1; j2++){ (*m_Out)<< "AmpBornU("<<j1<<","<<j2<<",*,*)=";
//   for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<AmpBornU.m_A[j1][j2][j][k];(*m_Out)<<endl;}
//(*m_Out)<< "----------------------------------------------------------------------------------------------------"<<endl;
//for(int j1=0; j1<=1; j1++) for(int j2=0; j2<=1; j2++){ (*m_Out)<< "AmpBornV("<<j1<<","<<j2<<",*,*)=";
//   for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<AmpBornV.m_A[j1][j2][j][k];(*m_Out)<<endl;}
//]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
// Fermion propagarotors
double pr1 = 1.0/(m_p1*ph1)/2.0;
double pr2 =-1.0/(m_p2*ph1)/2.0;
KKcmplx2 U, V;
int Sig = 1-2*Hel;  // sigma   = photon polarization (+1,-1)
dcmplx gI = dcmplx(Qe *m_e_QED);
//CALL GPS_MakeU(ph,Sig,  ph,mph,  p1,m1,    U)
//CALL GPS_MakeV(ph,Sig,  p2,m2,   ph,mph,   V)
GPS_MakeU(gI, ph1, Sig,  ph1,  m_p1,   U);
GPS_MakeV(gI, ph1, Sig,  m_p2,  ph1,   V);
//[[[[[[[[[[[[[[[[[[[[[
//(*m_Out)<< "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&KKceex::HiniPlus&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<endl;
//(*m_Out)<< "U(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<U.m_A[j][k]/gI;(*m_Out)<<endl;
//(*m_Out)<< "V(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<V.m_A[j][k]/gI;(*m_Out)<<endl;
//]]]]]]]]]]]]]]]]]]]]]
//========================================================================
/*   Version as in f77
// ISR non-infrared two parts: (1) p1 -> photon, contracted with U-matrix
dcmplx Csum1, Csum2;
dcmplx AmpBorn, AmpExpo1, AmpExpo2;
for(int j1 = 0; j1<=1; j1++)
  for(int j2 = 0; j2<=1; j2++){
    for(int j3 = 0; j3<=1; j3++)
      for(int j4 = 0; j4<=1; j4++){
            Csum1=dcmplx(0.0,0.0);
            Csum2=dcmplx(0.0,0.0);
            for(int  j= 0; j<=1; j++){
               Csum1=Csum1 +dcmplx(Qe *m_e_QED) *U.m_A[j][j1]*dcmplx(pr1) *AmpBornU.m_A[ j][j2][j3][j4];
               Csum2=Csum2 +dcmplx(Qe *m_e_QED) *V.m_A[j2][j]*dcmplx(pr2) *AmpBornV.m_A[j1][ j][j3][j4];
            }
            AmpBorn = m_BornC.m_A[j1][j2][j3][j4];
            AmpExpo1 =  sProd/Sactu*(Csum1+Csum2);
// (1+Vir1)*AmpBorn is already included in AmpExpo0 so we drop it to avoid double counting
// Note that remaining Vir2*AmpBorn is IR-finite because Vir2->0 in the IR limit
            AmpExpo2 =
                 sProd/Sactu*(Csum1+Csum2) *(1.0+Vir1+Vir2)*(1.0+m_F1fin1) // non-IR sigle bremss. part
                       +sProd*AmpBorn*Vir2;                             // add virtual_non_IR*Born
            m_AmpExpo1.m_A[j1][j2][j3][j4]  +=  AmpExpo1;
            m_AmpExpo2.m_A[j1][j2][j3][j4]  +=  AmpExpo2;
      }//for  j3,j4
  }//for j1,j2
  */
//========================================================================
//O(alf1)
  dcmplx Fact0= sProd/Sactu;
  AmpAddI(m_AmpExpo1, Fact0*dcmplx(pr1), AmpBornU, U);
  AmpAddI(m_AmpExpo1, Fact0*dcmplx(pr2), V, AmpBornV);
// O(alf2)
  dcmplx Fact2= Fact0*(1.0+Vir1+Vir2)*(1.0+m_F1fin1);
  AmpAddI(m_AmpExpo2, Fact2*dcmplx(pr1), AmpBornU, U);
  AmpAddI(m_AmpExpo2, Fact2*dcmplx(pr2), V, AmpBornV);
//[[[  AmpAdd( m_AmpExpo2, sProd*Vir2,m_BornC);
  AmpAdd( m_AmpExpo2, sProd*Vir2,m_BornD);             // m_BornD inludes W-exchange!
//[[[[[[[[[[[[[[[[[[[[
//(*m_Out)<<"++++++++++++++++++++++++++++++++++++++++++GPS_HiniPlus+++++++++++++++++++++++++++++++++++++++++++++"<<endl;
//(*m_Out)<<" Vir1,Vir2,m_F1fin1"<< Vir1<<"  "<<Vir2<<"  "<<m_F1fin1<<endl;
//(*m_Out)<<"---------------------------------------------------------------------------------------------------"<<endl;
//  for(int j1=0; j1<=1; j1++) for(int j2=0; j2<=1; j2++){ (*m_Out)<< "m_AmpExpo2("<<j1<<","<<j2<<",*,*)=";
//     for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<m_AmpExpo2.m_A[j1][j2][j][k];(*m_Out)<<endl;}
//]]]]]]]]]]]]]]]]]]]]]
}// GPS_HiniPlus

//*NEW*
void KKceex::HiniPlusW(int Ibeta, int KFini, int KFfin, TLorentzVector &PX,
               KKpart &ph1, int Hel, dcmplx &Sactu, dcmplx &sProd){
////////////////////////////////////////////////////////////////////////////////////
//   IR-finite part od 1-photon amplitudes for ISR W-exch only (ext. of GPS_Hini)  //
//   Photon helicity imported from the calling program                             //
//   Ibeta= 1 normal mode of operation                                             //
//   Ibeta=-1 removes action of Ibeta=1 (for m_AmpExpo2, m_AmpExpo2p)  part one    //
//   m_AmpExpo*  is working space                                                  //
/////////////////////////////////////////////////////////////////////////////////////
if(abs(KFfin) != 12)  return;
//[[[[[[[[[[[[[[[[[[[[[[[[[[[
//(*m_Out)<< "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&KKceex::HiniPlusW-&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<endl;
//]]]]]]]]]]]]]]]]]]]]]]]]]]]
double Qe  = DB->Qf[ KFini];
dcmplx Vir1,Vir2;  // Virtual corrections
MakeVini(m_p1, m_p2, ph1, Vir1,Vir2);
double CosThetD;
m_Event->ThetaD(PX,CosThetD);
double s0,t0;
s0 = PX*PX;
t0 = -s0*(1.0-CosThetD)/2.0;
//double u0 = -s0*(1.0+CosThetD)/2.0; //not used
int IFONE;
//Where is dominant other photon?  It is instead of reduction procedure
if( (m_p3.P[3]+m_p4.P[3])*m_p1.P[3] < 0.0) IFONE=1; else IFONE=0;
if(IFONE) { s0=  2*(m_p3*m_p4); t0= -2*(m_p4*m_p2);
} else {    s0=  2*(m_p3*m_p4); t0= -2*(m_p3*m_p1);
}//if
BornW(KFini, KFfin, PX, s0, t0,  m_r1, m_r2, m_p3, m_p4, m_AmpBornW); // for WW graph
//[[[[[[[[[[[[[[[[[[[[[
//(*m_Out)<< "KKceex::HiniPlusW: IFONE= "<< IFONE<<"  s0="<<s0<<"  t0="<<t0<<endl;
//]]]]]]]]]]]]]]]]]]]]]
//--------------------------------------------------------------
KKcmplx4 AmpBornU, AmpBornV; // dressed Born spin amplitudes
KKpart ph1r = ph1; ph1r.C=-1;
double sa,ta;
//Where is dominant other photon?  It is instead of reduction procedure
if( (ph1.P[3]+m_p3.P[3]+m_p4.P[3]) * m_p1.P[3] < 0.0) {
 IFONE= 1;      // We assume all extra photons were emitted from p1
} else {
 IFONE= 0;      // We assume all extra photons were emitted from p2
}
if(IFONE) {sa= (m_p3+m_p4+ph1).M2();     // (p3+p4+ph1)^2
           ta= -2*(m_p4*m_p2);}          // (p4-p2)^2
else      {sa= (m_p3+m_p4+ph1).M2();     // (p3+p4+ph1)^2
           ta= (m_p3-m_p1+ph1).M2();}    // (p3-p1+ph1)^2
BornW(KFini, KFfin, PX, sa, ta,  ph1, m_r2, m_p3, m_p4, AmpBornU); // single W, vector-like
//[[[[[[[[[[[[[[[[[[[[[
//(*m_Out)<< "KKceex::HiniPlusW: IFONE= "<< IFONE<<"  sa="<<sa<<"  ta="<<ta<<endl;
//]]]]]]]]]]]]]]]]]]]]]
//--------------------------------------------------------------------------------
double sb,tb;
if (IFONE) {sb= (m_p3+m_p4+ph1).M2();   // (p3+p4+ph1)**2
            tb= (m_p4-m_p2+ph1).M2();}  // (p4-p2+ph1)**2
else       {sb= (m_p3+m_p4+ph1).M2();   // (p3+p4+ph1)**2
            tb= -2*(m_p3*m_p1);}        // (p3-p1)**2
//[[[[[[[[[[[[[[[[[[[[[
//(*m_Out)<< "KKceex::HiniPlusW: IFONE= "<< IFONE<<"   sb="<<sb<<"  tb="<<tb<<endl;
//]]]]]]]]]]]]]]]]]]]]]
BornW(KFini, KFfin, PX, sb, tb, m_r1, ph1r, m_p3, m_p4, AmpBornV); // single W
//--------------------------------------------------------------------------------
dcmplx PropW0,WVPi0,PropWa,WVPia,PropWb,WVPib;
GPS_EWFFactW(KFini,KFfin,s0,t0,PropW0,WVPi0);
GPS_EWFFactW(KFini,KFfin,sa,ta,PropWa,WVPia);
GPS_EWFFactW(KFini,KFfin,sb,tb,PropWb,WVPib);
WVPia=WVPi0; // to keep gauge invariance we install t-transfer in formfactor at 0 order
WVPib=WVPi0; // to keep gauge invariance we install t-transfer in formfactor at 0 order
// Fermion propagarotors
double pr1  = 1.0/(m_p1*ph1)/2.0;
double pr2  =-1.0/(m_p2*ph1)/2.0;
//[[[[[[[[[[[[[[[[[[[[[
//(*m_Out)<< "KKceex::HiniPlusW:   pr1="<<pr1<<"  pr2="<<pr2<<endl;
//]]]]]]]]]]]]]]]]]]]]]
int Sig = 1-2*Hel;
dcmplx Cone = dcmplx(1.0);
KKcmplx2 U,V, UW,VW, UWX,VWX;
//CALL GPS_MakeU(ph,Sig,  ph,mph,  p1,m1,    U)
//CALL GPS_MakeV(ph,Sig,  p2,m2,   ph,mph,   V)
GPS_MakeU(Cone, ph1, Sig,  ph1,  m_p1,   U);
GPS_MakeV(Cone, ph1, Sig,  m_p2,  ph1,   V);
//CALL GPS_MakeUW(Cnor,ph,Sig,  p3,m3,   p1,m1,    UW)
//CALL GPS_MakeVW(Cnor,ph,Sig,  p2,m2,   p4,m4,    VW)
GPS_MakeUW(Cone, ph1, Sig,  m_p3,  m_p1,   UW);
GPS_MakeVW(Cone, ph1, Sig,  m_p2,  m_p4,   VW);
//CALL GPS_MakeUX(Cnor,ph,Fleps, p3,m3,   p1,m1,    UWX) ! v-a inside
//CALL GPS_MakeVX(Cnor,ph,Fleps, p2,m2,   p4,m4,    VWX) ! v-a inside
GPS_MakeUX(Cone, ph1, m_p3, m_p1, UWX);
GPS_MakeVX(Cone, ph1, m_p2, m_p4, VWX);
//[[[[[[[[[[[[[[[[[[[[[
//(*m_Out)<< "U(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<U.m_A[j][k];(*m_Out)<<endl;
//(*m_Out)<< "V(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<V.m_A[j][k];(*m_Out)<<endl;
//(*m_Out)<< "UW(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<UW.m_A[j][k];(*m_Out)<<endl;
//(*m_Out)<< "VW(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<VW.m_A[j][k];(*m_Out)<<endl;
//(*m_Out)<< "UWX(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<UWX.m_A[j][k];(*m_Out)<<endl;
//(*m_Out)<< "VWX(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<VWX.m_A[j][k];(*m_Out)<<endl;
//]]]]]]]]]]]]]]]]]]]]]
dcmplx s1v[2],s2v[2],EpsDot[2]; // soft factors
//s1v(1)  = -GPS_Sof1( 1,phv,p1v)
//s2v(1)  =  GPS_Sof1( 1,phv,p2v)
s1v[0]  = -GPS_Sof1( 1,ph1,m_p1);
s2v[0]  =  GPS_Sof1( 1,ph1,m_p2);
//IF (IFONE) THEN
//  EpsDot(1)=(+GPS_Sof1x( 1,phv,p2v)-GPS_Sof1x( 1,phv,p4v)) ! minis sign is in pr2/4
//ELSE
//  EpsDot(1)=(-GPS_Sof1x( 1,phv,p1v)+GPS_Sof1x( 1,phv,p3v)) ! minis sign is in pr2/4
//ENDIF
if(IFONE) EpsDot[0]= +GPS_Sof1x( 1,ph1,m_p2)-GPS_Sof1x( 1,ph1,m_p4); // minus is in pr2/4
else      EpsDot[0]= -GPS_Sof1x( 1,ph1,m_p1)+GPS_Sof1x( 1,ph1,m_p3); // minus is in pr2/4
s1v[1]   = -conj(s1v[0]);
s2v[1]   = -conj(s2v[0]);
EpsDot[1]= -conj(EpsDot[0]);
//[[[[[[[[[[[[[[[[[[[[[
//(*m_Out)<< "PropW0="<<PropW0<< "WVPi0="<<WVPi0<<endl;
//(*m_Out)<< "PropWa="<<PropWa<< "WVPia="<<WVPia<<endl;
//(*m_Out)<< "PropWb="<<PropWb<< "WVPib="<<WVPib<<endl;
//(*m_Out)<< "s1v(*)=";for(int j=0;j<=1;j++)  (*m_Out)<<"  "<<SW208<<s1v[j];(*m_Out)<<endl;
//(*m_Out)<< "s2v(*)=";for(int j=0;j<=1;j++)  (*m_Out)<<"  "<<SW208<<s2v[j];(*m_Out)<<endl;
//(*m_Out)<< "EpsDot(*)=";for(int j=0;j<=1;j++)  (*m_Out)<<"  "<<SW208<<EpsDot[j];(*m_Out)<<endl;
//(*m_Out)<< "Hel="<<Hel<< "  Sig="<<Sig<<endl;
//(*m_Out)<<"-----------------------------------------------------------------------------------------------"<<endl;
//  for(int j1=0; j1<=1; j1++) for(int j2=0; j2<=1; j2++){ (*m_Out)<< "m_AmpExpo1("<<j1<<","<<j2<<",*,*)=";
//     for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<m_AmpExpo1.m_A[j1][j2][j][k];(*m_Out)<<endl;}
//  for(int j1=0; j1<=1; j1++) for(int j2=0; j2<=1; j2++){ (*m_Out)<< "m_AmpBornW("<<j1<<","<<j2<<",*,*)=";
//     for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<m_AmpBornW.m_A[j1][j2][j][k];(*m_Out)<<endl;}
//ZerAmplit();
//(*m_Out)<< "m_e_QED="<<m_e_QED<<endl;
//(*m_Out)<< "s1v[Hel]="<<s1v[Hel]<<" s2v[Hel]="<<s2v[Hel]<<" EpsDot[Hel]="<<EpsDot[Hel]<<endl;
//]]]]]]]]]]]]]]]]]]]]]
//-----------------------------------------------------------
//O(alf1)
//    Csum1=Csum1 +Ibeta*DCMPLX(Qe *m_e_QED) *U(j,j1)*pr1 *AmpBornU( j,j2,j3,j4)
//    Csum2=Csum2 +Ibeta*DCMPLX(Qe *m_e_QED) *V(j2,j)*pr2 *AmpBornV(j1, j,j3,j4)
dcmplx gI = dcmplx(Qe *m_e_QED); //??? what about Qe?
dcmplx Fact0= gI*sProd/Sactu;
//
Amp4Zer(m_AmpTemp1);
//
AmpAddI(m_AmpTemp1, Fact0*dcmplx(pr1), AmpBornU, U);
AmpAddI(m_AmpTemp1, Fact0*dcmplx(pr2), V, AmpBornV);
/////////////////
//Csum4=Csum4 +DCMPLX(    m_e_QED)*PropWa*WVPia*PropWb*WVPib/WVpi0 !denominator for gauge invariance see up
//$                    *2            ! from  feynman diagram
//$                    *(-0.5D0)     ! fixup originating from  test in  GPS_BornWPlusT
//$                    *(UW(j3,j1)*VWX(j2,j4)-VW(j2,j4)*UWX(j3,j1)) ! non-infrared part of emission from W
dcmplx Fact4= dcmplx(m_e_QED)*PropWa*WVPia*PropWb*WVPib/WVPi0
                       *2.0           // from  feynman diagram
                       *(-0.5);       // fixup originating from  test in  GPS_BornWPlusT
Fact4 *= sProd/Sactu;
AmpAdd( m_AmpTemp1, Fact4 , UW, VWX, VW, UWX);
/////////////////
//Csum3=Csum3 +DCMPLX(m_e_QED)*AmpBornW(j1,j2,j3,j4)/PropW0/WVPi0*(
//$                      Qe*s1v(Hel) *(PropWa*WVPia-PropW0*WVPi0)         ! t-channel W-prop variation
//$                    + Qe*s2v(Hel) *(PropWb*WVPib-PropW0*WVPi0)         ! t-channel W-prop variation
//$                    + EpsDot(Hel)* PropWa*WVPia*PropWb*WVPib/  WVPib   ! basically IR emis. from W, reduction procedure used only here
dcmplx Fact3= dcmplx(m_e_QED)/PropW0/WVPi0*(
      Qe*s1v[Hel] *(PropWa*WVPia-PropW0*WVPi0)           // t-channel W-prop variation
    + Qe*s2v[Hel] *(PropWb*WVPib-PropW0*WVPi0)           // t-channel W-prop variation
    + EpsDot[Hel]* PropWa*WVPia*PropWb*WVPib/WVPib   );  // basically IR emis. from W, reduction procedure used only here
Fact3 *= sProd/Sactu;
AmpAdd( m_AmpTemp1, Fact3,m_AmpBornW);
//
if(Ibeta ==  1) {
  AmpAdd( m_AmpExpo1, dcmplx(1.0), m_AmpTemp1);
  AmpAdd( m_AmpExpo2, dcmplx(1.0+Vir1+Vir2), m_AmpTemp1);
}
if(Ibeta == -1){
  AmpAdd( m_AmpExpo2, dcmplx(-1.0), m_AmpTemp1);
}

//[[[[[[[[[[[[[[[[[[[[[
//(*m_Out)<< "s1v(*)=";for(int j=0;j<=1;j++)  (*m_Out)<<"  "<<SW208<<s1v[j];(*m_Out)<<endl;
//(*m_Out)<< "s2v(*)=";for(int j=0;j<=1;j++)  (*m_Out)<<"  "<<SW208<<s2v[j];(*m_Out)<<endl;
//(*m_Out)<< "EpsDot(*)=";for(int j=0;j<=1;j++)  (*m_Out)<<"  "<<SW208<<EpsDot[j];(*m_Out)<<endl;
//(*m_Out)<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
//  for(int j1=0; j1<=1; j1++) for(int j2=0; j2<=1; j2++){ (*m_Out)<< "m_AmpExpo1("<<j1<<","<<j2<<",*,*)=";
//     for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<m_AmpExpo1.m_A[j1][j2][j][k];(*m_Out)<<endl;}
//]]]]]]]]]]]]]]]]]]]]]
//
}//HiniPlusW

void KKceex::HfinPlus(int KFini, int KFfin, TLorentzVector &PX,
             KKpart &ph1, int Hel, dcmplx &Sactu, dcmplx &sProd, double &CKine){
/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
//   IR-finite part od 1-photon amplitudes for FSR  (equiv. to GPS_HfinPlus)       //
//   Photon helicity is give by the calling program                                //
//                                                                                 //
//   Missing contribution in FSR non-IR part due to  svarX/svarQ                   //
//   Contribution -svarX/svarQ from HERE cancels exactly with svarX/svarQ in beta0 //
//                                                                                 //
//   m_AmpExpo*  is working space                                                  //
//   m_AmpBorn   is hidden INPUT                                                   //
//                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////
  double Qf  = DB->Qf[ KFfin];
  dcmplx Vir1,Vir2;  // Virtual corrections
  MakeVfin(m_p3, m_p4, ph1, Vir1,Vir2);
// FSR non-infrared two parts: (1) p1 -> photon, contracted with U-matrix
//                             (2) p2 -> photon, contracted with V-matrix
// CALL GPS_Born(KFi,KFf,PX, p1,Fleps, p2,-Fleps,  ph,mph,   p4,-m4,   AmpBornU)
// CALL GPS_Born(KFi,KFf,PX, p1,Fleps, p2,-Fleps,  p3,m3,    ph,-mph,  AmpBornV)
  KKcmplx4 AmpBornU, AmpBornV; // dressed Born spin amplitudes
  double CosThe=0.0;
  KKpart ph1r = ph1; ph1r.C=-1;
  Born(KFini,KFfin, PX, CosThe, m_r1, m_r2,  ph1,  m_p4,  AmpBornU);
  Born(KFini,KFfin, PX, CosThe, m_r1, m_r2, m_p3,  ph1r,  AmpBornV);
// Fermion propagarotors
  double pr1 = 1.0/(m_p3*ph1)/2.0;
  double pr2 =-1.0/(m_p4*ph1)/2.0;
  KKcmplx2 U, V;
  int Sig = 1-2*Hel;  //
  dcmplx gF = dcmplx(Qf *m_e_QED);
  GPS_MakeU(gF, ph1, Sig, m_p3,  ph1, U);
  GPS_MakeV(gF, ph1, Sig,  ph1, m_p4, V);
//=======================================================================================
/*
// version as in f77 code
  dcmplx Csum1, Csum2;
  dcmplx AmpBorn, AmpExpo1, AmpExpo2;
  for(int j1 = 0; j1<=1; j1++)
    for(int j2 = 0; j2<=1; j2++)
      for(int j3 = 0; j3<=1; j3++)
        for(int j4 = 0; j4<=1; j4++){
          Csum1=dcmplx(0.0,0.0);
          Csum2=dcmplx(0.0,0.0);
          for(int j=0; j<=1; j++){
            Csum1=Csum1 +dcmplx(Qf *m_e_QED) *U.m_A[j3][j]*pr1* AmpBornU.m_A[j1][j2][ j][j4];
            Csum2=Csum2 +dcmplx(Qf *m_e_QED) *V.m_A[j][j4]*pr2* AmpBornV.m_A[j1][j2][j3][ j];
          }//for
          AmpBorn = m_BornC.m_A[j1][j2][j3][j4];
///// first order
          AmpExpo1 =
            +sProd/Sactu*(Csum1+Csum2)             // non-IR sigle bremss. part
            -sProd*CKine*AmpBorn +sProd*AmpBorn;   // compensate for (svarX1/svarQ)
///// second order
// (1+Vir1)*AmpBorn is already included in AmpExpo2 so we drop it to avoid double counting
// the remaining Vir2*AmpBorn is IR-finite because Vir2->0 in the IR limit
          AmpExpo2 =
            +sProd/Sactu*(Csum1+Csum2)*(1.0+Vir1+Vir2)*(1.0+m_F1ini1) // non-IR sigle bremss. part
            +sProd*(-CKine*AmpBorn+AmpBorn)*(1.0+Vir1)*(1.0+m_F1ini1) // compensate for (svarX1/svarQ)
            +sProd*AmpBorn*Vir2;                                      // add virtual_non_IR*Born
          m_AmpExpo1.m_A[j1][j2][j3][j4] += AmpExpo1;
          m_AmpExpo2.m_A[j1][j2][j3][j4] += AmpExpo2;
        }// for j1-j4
*/
//=======================================================================================
//O(alf1)
  dcmplx Fact0= sProd/Sactu;
  AmpAddF(m_AmpExpo1, Fact0*dcmplx(pr1), U, AmpBornU);
  AmpAddF(m_AmpExpo1, Fact0*dcmplx(pr2), AmpBornV, V);
  dcmplx Fact1= sProd*(1-CKine);   // (1-CKine)->0 in the IR limit
  AmpAdd( m_AmpExpo1, Fact1 ,m_BornC);
// O(alf2)
  dcmplx Fact2= Fact0*(1.0+Vir1+Vir2)*(1.0+m_F1ini1);
  AmpAddF(m_AmpExpo2, Fact2*dcmplx(pr1), U, AmpBornU );
  AmpAddF(m_AmpExpo2, Fact2*dcmplx(pr2), AmpBornV, V);
  dcmplx Fact3= sProd*(1-CKine)*(1.0+Vir1)*(1.0+m_F1ini1)
		      + sProd*Vir2;      //Vir2->0 in the IR limit
  AmpAdd( m_AmpExpo2, Fact3 ,m_BornC);
}// GPS_HfinPlus !!!

void KKceex::GPS_HffPlus(dcmplx CNorm, int KFini, int KFfin, TLorentzVector PX,
		                 KKpart ph1, int Hel1, KKpart ph2, int Hel2){
/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
//   Genuine IR-finite non 1-photon amplitudes for FSR-FSR are added to AmpWork    //
//   Photon helicity imported from the calling program.                            //
//                                                                                 //
//   For Y_IR=1 IR-finite 1-photon contributions are calculated here as well.      //
//                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
//                                    FSR                                            //
//                                                                                   //
//                            1                  2                                   //
//                            |                  |                                   //
//             c              |                  |          d                        //
//    u  ------<------OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO----<----- v                 //
//                                     |                                             //
//                                     |X                                            //
//                                     |                                             //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////
  int  Y_IR=0;                  // No, IR not included
  int N_IR=1-Y_IR;
  double ChaFin =  DB->Qf[ KFfin];
  dcmplx gF = dcmplx(ChaFin*m_e_QED);
  dcmplx sC[2][2], sD[2][2];
  KKpart pC=m_p3;
  KKpart pD=m_p4;
  sC[0][0] =  gF *GPS_Sof1( 1,ph1,pC);
  sC[1][0] =  gF *GPS_Sof1( 1,ph2,pC);
  sD[0][0] = -gF *GPS_Sof1( 1,ph1,pD);
  sD[1][0] = -gF *GPS_Sof1( 1,ph2,pD);
  sC[0][1] = -conj(sC[0][0]);
  sC[1][1] = -conj(sC[1][0]);
  sD[0][1] = -conj(sD[0][0]);
  sD[1][1] = -conj(sD[1][0]);
// Calculate Born spin amplitudes, also with substitutions
//  CALL GPS_Born(KFi,KFf,PX, pA,mA,  pB,-mB,  pC,   mC,   pD,   -mD, BornABCD) ! Standard
//  CALL GPS_Born(KFi,KFf,PX, pA,mA,  pB,-mB,  ph1, mph,   pD,   -mD, BornAB1D) ! C->1
//  CALL GPS_Born(KFi,KFf,PX, pA,mA,  pB,-mB,  pC,   mC,   ph1, -mph, BornABC1) ! D->1
//  CALL GPS_Born(KFi,KFf,PX, pA,mA,  pB,-mB,  ph2, mph,   pD,   -mD, BornAB2D) ! C->2
//  CALL GPS_Born(KFi,KFf,PX, pA,mA,  pB,-mB,  pC,   mC,   ph2, -mph, BornABC2) ! D->2
//  CALL GPS_Born(KFi,KFf,PX, pA,mA,  pB,-mB,  ph1, mph,   ph2, -mph, BornAB12) ! C->1,D->2
//  CALL GPS_Born(KFi,KFf,PX, pA,mA,  pB,-mB,  ph2, mph,   ph1, -mph, BornAB21) ! C->2,D->1
// Calculate Born spin amplitudes, also with substitutions
  double CosThe=0.0;
  KKpart ph1r = ph1; ph1r.C=-1;
  KKpart ph2r = ph2; ph2r.C=-1;
  KKcmplx4 BornABCD, BornAB1D, BornABC1, BornAB2D,BornABC2, BornAB12, BornAB21;
  Born(KFini,KFfin,PX, CosThe,  m_p1, m_p2,  m_p3, m_p4, BornABCD);  // Standard
  Born(KFini,KFfin,PX, CosThe,  m_p1, m_p2,  ph1,  m_p4, BornAB1D);  // C->1
  Born(KFini,KFfin,PX, CosThe,  m_p1, m_p2,  m_p3,  ph1, BornABC1);  // D->1
  Born(KFini,KFfin,PX, CosThe,  m_p1, m_p2,  ph2,  m_p4, BornAB2D);  // C->2
  Born(KFini,KFfin,PX, CosThe,  m_p1, m_p2,  m_p3,  ph2, BornABC2);  // D->2
  Born(KFini,KFfin,PX, CosThe,  m_p1, m_p2,  ph1,   ph2, BornAB12);  // C->1,D->2
  Born(KFini,KFfin,PX, CosThe,  m_p1, m_p2,  ph2,   ph1, BornAB21);  // C->2,D->1
// propagators
  KKpart QQ   = pC;  QQ +=pD;    // pC+pD
  KKpart PP1  = QQ;  PP1 +=ph1;  // pC+pD+ph1;
  KKpart PP2  = QQ;  PP2 +=ph2;  // pC+pD+ph2;
  KKpart PP12 = PP1; PP12 += ph2; // (pC+pD)+ph1+ph2;
  double svarX12 = PP12*PP12;
  double svarX1  = PP1*PP1;
  double svarX2  = PP2*PP2;
  double svarQ   = QQ*QQ;
// Fermion propagarotors 1
  double prC1=  1.0/(pC*ph1)/2.0;
  double prD1= -1.0/(pD*ph1)/2.0;
// Fermion propagarotors 2
  double prC2=  1.0/(pC*ph2)/2.0;
  double prD2= -1.0/(pD*ph2)/2.0;
// Double propagators
  double prC12= 1.0/( pC*ph1 +pC*ph2 +ph1*ph2)/2.0;
  double prD12=-1.0/( pD*ph1 +pD*ph2 +ph1*ph2)/2.0;
  double Fprop1= (1.0/prC1+1.0/prC2)*prC12 -1.0;
  double Fprop2= (1.0/prD1+1.0/prD2)*prD12 -1.0;
//
  int Sig = 1-2*Hel1;
  KKcmplx2 V11d, Uc11, V21d, Uc12, Vd12, U21c, V112, U211, V212, U212;
  GPS_MakeV( gF, ph1,Sig,  ph1, pD,    V11d);   // <1|{1}|D>
  GPS_MakeU( gF, ph1,Sig,  pC,   ph1,  Uc11);   // <C|[1]|1>
// false second
  GPS_MakeV( gF, ph1,Sig,  ph2, pD,    V21d);   // <2|{1}|D>
  GPS_MakeU( gF, ph1,Sig,  pC,   ph2,  Uc12);   // <C|[1]|2>
// reverse order
  GPS_MakeV( gF, ph1,Sig,  pD,  ph2,   Vd12);   // <D|{1}|2>
  GPS_MakeU( gF, ph1,Sig,  ph2, pC,    U21c);   // <2|[1]|C>
// xk-xk term case, ph2 first
  GPS_MakeV( gF, ph1,Sig,  ph1, ph2,   V112);   // <1|{1}|2>
  GPS_MakeU( gF, ph1,Sig,  ph2, ph1,   U211);   // <2|[1]|1>
// xk-xk term case ph2 first
  GPS_MakeV( gF, ph1,Sig,  ph2, ph2,   V212);   // <2|{1}|2>
  GPS_MakeU( gF, ph1,Sig,  ph2, ph2,   U212);   // <2|[1]|2>
//
  Sig = 1-2*Hel2;
  KKcmplx2 V22d, Uc22, V12d, Uc21, Vd21, U12c, V221, U122, V121, U121;
  GPS_MakeV( gF, ph2,Sig,  ph2, pD,    V22d);   // <2|{2}|D>
  GPS_MakeU( gF, ph2,Sig,  pC,  ph2,   Uc22);   // <C|[2]|2>
// falSe second
  GPS_MakeV( gF, ph2,Sig,  ph1, pD,    V12d);   // <1|{2}|D>
  GPS_MakeU( gF, ph2,Sig,  pC,  ph1,   Uc21);   // <C|[2]|1>
// reverse order
  GPS_MakeV( gF, ph2,Sig,  pD,  ph1,   Vd21);   // <D|{2}|1>
  GPS_MakeU( gF, ph2,Sig,  ph1, pC,    U12c);   // <1|[2]|C>
// xk-xk term, ph1 first
  GPS_MakeV( gF, ph2,Sig,  ph2, ph1,   V221);   // <2|{2}|1>
  GPS_MakeU( gF, ph2,Sig,  ph1, ph2,   U122);   // <1|[2]|2>
// xk-xk term, ph1 first
  GPS_MakeV( gF, ph2,Sig,  ph1, ph1,   V121);   // <1|{2}|1>
  GPS_MakeU( gF, ph2,Sig,  ph1, ph1,   U121);   // <1|[2]|1>

  ///////////////////////////////////////////////////////////////////////////////////////
  //                     1|            2|                                              //
  //               c      |    c+m+1    |    c+m+1+2          -d                       //
  //       u  -----<------S-----<-------U------<-------O-------<------ v               //
  //                                                   |X                              //
  ///////////////////////////////////////////////////////////////////////////////////////
  //    Su1=Su1 +(sC(1,Hel1)*(prC12-prC2*N_IR))*Uc22(j3,j)*BornAB2D(j1,j2,j,j4) !<c|(1)c[2]2|X|d>
  //    Su1=Su1  +sC(1,Hel1)* prC12            *Uc21(j3,j)*BornAB1D(j1,j2,j,j4) !<c|(1)c[2]1|X|d>
  AmpAddF(m_AmpExpo2, CNorm*sC[0][Hel1]*dcmplx(prC12-prC2*N_IR), Uc22, BornAB2D);
  AmpAddF(m_AmpExpo2, CNorm*sC[0][Hel1]*dcmplx(prC12          ), Uc21, BornAB1D);

  ///////////////////////////////////////////////////////////////////////////////////////
  //                     2|            1|                                              //
  //               c      |    c+m+2    |    c+m+1+2          -d                       //
  //       u  -----<------S-----<-------U------<-------O-------<------ v               //
  //                                                   |X                              //
  ///////////////////////////////////////////////////////////////////////////////////////
  //    Su1=Su1+(sC(2,Hel2)*(prC12-prC1*N_IR))*Uc11(j3,j)*BornAB1D(j1,j2,j,j4) !<c|(2)c[1]1|X|d>
  //    Su1=Su1 +sC(2,Hel2)* prC12            *Uc12(j3,j)*BornAB2D(j1,j2,j,j4) !<c|(2)c[1]2|X|d>
  AmpAddF(m_AmpExpo2, CNorm*sC[1][Hel2]*dcmplx(prC12-prC1*N_IR), Uc11, BornAB1D);
  AmpAddF(m_AmpExpo2, CNorm*sC[1][Hel2]*dcmplx(prC12          ), Uc12, BornAB2D);
  ///////////////////////////////////////////////////////////////////////////////////////
  //                                    |1             |2                              //
  //               c          -d+m-1-2  |   -d+m-2     |       -d                      //
  //       u  -----<------O-----<-------V-----<--------S--------<----- v               //
  //                     X|                                                            //
  ///////////////////////////////////////////////////////////////////////////////////////
  //    Su1=Su1 +BornABC1(j1,j2,j3,j)*V11d(j,j4)*( sD(2,Hel2)*(prD12-prD1*N_IR))!<c|X|1{1}d(2)|d>
  //    Su1=Su1 +BornABC2(j1,j2,j3,j)*V21d(j,j4)*  sD(2,Hel2)* prD12            !<c|X|2{1}d(2)|d>
  AmpAddF(m_AmpExpo2, CNorm*sD[1][Hel2]*dcmplx(prD12-prD1*N_IR), BornABC1,V11d);
  AmpAddF(m_AmpExpo2, CNorm*sD[1][Hel2]*dcmplx(prD12          ), BornABC2,V21d);
  ///////////////////////////////////////////////////////////////////////////////////////
  //                                    |2             |1                              //
  //               c          -d+m-1-2  |   -d+m-1     |       -d                      //
  //       u  -----<------O-----<-------V-----<--------S--------<----- v               //
  //                     X|                                                            //
  ///////////////////////////////////////////////////////////////////////////////////////
  //     Su1=Su1 +BornABC2(j1,j2,j3,j)*V22d(j,j4)*( sD(1,Hel1)*(prD12-prD2*N_IR))!<c|X|2{2}d(1)|d>
  //     Su1=Su1 +BornABC1(j1,j2,j3,j)*V12d(j,j4)  *sD(1,Hel1)* prD12            !<c|X|1{2}d(1)|d>
  AmpAddF(m_AmpExpo2, CNorm*sD[0][Hel1]*dcmplx(prD12-prD2*N_IR), BornABC2,V22d);
  AmpAddF(m_AmpExpo2, CNorm*sD[0][Hel1]*dcmplx(prD12          ), BornABC1,V12d);
  //
  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////
  //                     2|                            |1                              //
  //               c      |    c+m+1         c+m+1     |      -d                       //
  //       u  -----<------U-----<-------O------<-------S-------<------ v               //
  //                                    |X                                             //
  ///////////////////////////////////////////////////////////////////////////////////////
  //       Su3=Su3 +Uc22(j3,j)*prC2  *BornAB2D(j1,j2,j,j4) *sD(1,Hel1) *Y_IR !<c|[2]c|X|1(1)|d>
  if( Y_IR )  AmpAddF(m_AmpExpo2, CNorm*dcmplx(prC2)*sD[0][Hel1], Uc22, BornAB2D);
  ///////////////////////////////////////////////////////////////////////////////////////
  //                     1|                            |2                              //
  //               c      |    c+m+1         c+m+2     |      -d                       //
  //       u  -----<------U-----<-------O------<-------S-------<------ v               //
  //                                    |X                                             //
  ///////////////////////////////////////////////////////////////////////////////////////
  //       Su3=Su3 +Uc11(j3,j)*prC1 *BornAB1D(j1,j2,j,j4)  *sD(2,Hel2) *Y_IR !<c|[1]c|X|2(2)|d>
  if( Y_IR )  AmpAddF(m_AmpExpo2, CNorm*dcmplx(prC1)*sD[1][Hel2], Uc11, BornAB1D);
  ///////////////////////////////////////////////////////////////////////////////////////
  //                     2|                            |1                              //
  //               c      |    c+m+2         c+m+1     |      -d                       //
  //       u  -----<------S-----<-------O------<-------U-------<------ v               //
  //                                    |X                                             //
  ///////////////////////////////////////////////////////////////////////////////////////
  //       Su3=Su3 +sC(2,Hel2) *BornABC1(j1,j2,j3,j) *prD1 *V11d(j,j4) *Y_IR !<c|(2)c|X|1{1}|d>
  if( Y_IR )  AmpAddF(m_AmpExpo2, CNorm*sC[1][Hel2]*dcmplx(prD1), BornABC1, V11d);
  ///////////////////////////////////////////////////////////////////////////////////////
  //                     1|                            |2                              //
  //               c      |    c+m+1         c+m+2     |      -d                       //
  //       u  -----<------S-----<-------O------<-------U-------<------ v               //
  //                                    |X                                             //
  ///////////////////////////////////////////////////////////////////////////////////////
  //      Su3=Su3 +sC(1,Hel1) *BornABC2(j1,j2,j3,j) *prD2 *V22d(j,j4) *Y_IR !<c|(1)c|X|2{2}|d>
  if( Y_IR )  AmpAddF(m_AmpExpo2, CNorm*sC[0][Hel1]*dcmplx(prD2), BornABC2, V22d);
//
  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////
  //                      |2                           |1                              //
  //               c      |    c+m+2        -d+m-1     |       -d                      //
  //       u  -----<------U-----<-------O-----<--------V--------<----- v               //
  //                                    |X                                             //
  ///////////////////////////////////////////////////////////////////////////////////////
  //      Su2=Su2 +Uc22( j3,l)*prC2 *BornAB21(j1,j2,l,j ) *V11d( j,j4)*prD1 !<c|[2]2|X|1{1}|d>
  AmpAddF(m_AmpExpo2, CNorm* dcmplx( prC2*prD1 ), Uc22, BornAB21, V11d);
  ///////////////////////////////////////////////////////////////////////////////////////
  //                      |1                           |2                              //
  //               c      |    c+m+1        -d+m-2     |       -d                      //
  //       u  -----<------U-----<-------O-----<--------V--------<----- v               //
  //                                    |X                                             //
  ///////////////////////////////////////////////////////////////////////////////////////
  //     Su2=Su2 +Uc11( j3,l)*prC1 *BornAB12(j1,j2,l,j ) *V22d( j,j4)*prD2 !<c|[1]1|X|2{2}|d>
  AmpAddF(m_AmpExpo2, CNorm* dcmplx( prC1*prD2 ), Uc11, BornAB12, V22d);
  ///////////////////////////////////////////////////////////////////////////////////////
  //                     1|            2|                                              //
  //               c      |    c+m+1    |    c+m+1+2          -d                       //
  //       u  -----<------U-----<-------O------<-------V-------<------ v               //
  //                                                  X|                               //
  ///////////////////////////////////////////////////////////////////////////////////////
  //     Su2=Su2 +Uc11( j3,l)*prC1  *U122(l,j)*prC12  *BornAB2D(j1,j2,j ,j4) ! <c|[1]1[2]2|X|d>
  //     Su2=Su2 +Uc11( j3,l)*prC1  *U121(l,j)*prC12  *BornAB1D(j1,j2,j ,j4) ! <c|[1]1[2]1|X|d>
  //     Su2=Su2 +Uc11( j3,l)*prC1  *U12c(l,j)*prC12  *BornABCD(j1,j2,j ,j4) ! <c|[1]1[2]c|X|d>
  AmpAddF(m_AmpExpo2, CNorm* dcmplx( prC1*prC12 ), Uc11, U122, BornAB2D);
  AmpAddF(m_AmpExpo2, CNorm* dcmplx( prC1*prC12 ), Uc11, U121, BornAB1D);
  AmpAddF(m_AmpExpo2, CNorm* dcmplx( prC1*prC12 ), Uc11, U12c, BornABCD);
  ///////////////////////////////////////////////////////////////////////////////////////
  //                     2|            1|                                              //
  //               c      |    c+m+2    |    c+m+1+2          -d                       //
  //       u  -----<------U-----<-------U------<-------O-------<------ v               //
  //                                                  X|                               //
  ///////////////////////////////////////////////////////////////////////////////////////
  //     Su2=Su2 +Uc22( j3,l)*prC2  *U211(l,j)*prC12  *BornAB1D(j1,j2,j ,j4) ! <c|[2]2[1]1|X|d>
  //     Su2=Su2 +Uc22( j3,l)*prC2  *U212(l,j)*prC12  *BornAB2D(j1,j2,j ,j4) ! <c|[2]2[1]2|X|d>
  //     Su2=Su2 +Uc22( j3,l)*prC2  *U21c(l,j)*prC12  *BornABCD(j1,j2,j ,j4) ! <c|[2]2[1]c|X|d>
  AmpAddF(m_AmpExpo2, CNorm* dcmplx( prC2*prC12 ), Uc22, U211, BornAB1D);
  AmpAddF(m_AmpExpo2, CNorm* dcmplx( prC2*prC12 ), Uc22, U212, BornAB2D);
  AmpAddF(m_AmpExpo2, CNorm* dcmplx( prC2*prC12 ), Uc22, U21c, BornABCD);
  ///////////////////////////////////////////////////////////////////////////////////////
  //                                    |2             |1                              //
  //               c          -d+m-1-2  |   -d+m-1     |       -d                      //
  //       u  -----<------O-----<-------V-----<--------V--------<----- v               //
  //                      |X                                                           //
  ///////////////////////////////////////////////////////////////////////////////////////
  //      Su2=Su2 +BornABC2(j1,j2,j3,j ) *V221(j,l )*prD12 *V11d( l,j4)*prD1  ! <c|X|2{2}1{1}|d>
  //      Su2=Su2 +BornABC1(j1,j2,j3,j ) *V121(j,l )*prD12 *V11d( l,j4)*prD1  ! <c|X|1{2}1{1}|d>
  //      Su2=Su2 +BornABCD(j1,j2,j3,j ) *Vd21(j,l )*prD12 *V11d( l,j4)*prD1  ! <c|X|d{2}1{1}|d>
  AmpAddF(m_AmpExpo2, CNorm* dcmplx( prD12*prD1 ), BornABC2, V221, V11d);
  AmpAddF(m_AmpExpo2, CNorm* dcmplx( prD12*prD1 ), BornABC1, V121, V11d);
  AmpAddF(m_AmpExpo2, CNorm* dcmplx( prD12*prD1 ), BornABCD, Vd21, V11d);
  ///////////////////////////////////////////////////////////////////////////////////////
  //                                    |1             |2                              //
  //               c          -d+m-1-2  |   -d+m-2     |       -d                      //
  //       u  -----<------O-----<-------V-----<--------V--------<----- v               //
  //                      |X                                                           //
  ///////////////////////////////////////////////////////////////////////////////////////
  //      Su2=Su2 +BornABC1(j1,j2,j3,j ) *V112(j,l )*prD12 *V22d( l,j4)*prD2  ! <c|X|1{1}2{2}|d>
  //      Su2=Su2 +BornABC2(j1,j2,j3,j ) *V212(j,l )*prD12 *V22d( l,j4)*prD2  ! <c|X|2{1}2{2}|d>
  //      Su2=Su2 +BornABCD(j1,j2,j3,j ) *Vd12(j,l )*prD12 *V22d( l,j4)*prD2  ! <c|X|d{1}2{2}|d>
  AmpAddF(m_AmpExpo2, CNorm* dcmplx( prD12*prD2 ), BornABC1, V112, V22d);
  AmpAddF(m_AmpExpo2, CNorm* dcmplx( prD12*prD2 ), BornABC2, V212, V22d);
  AmpAddF(m_AmpExpo2, CNorm* dcmplx( prD12*prD2 ), BornABCD, Vd12, V22d);
///////////////////////////////////////////////////////////////////////////////////////
//  sProd = (sC(1,Hel1)+sD(1,Hel1)) *( sC(2,Hel2)+sD(2,Hel2)) !
//  AmpWork(j1,j2,j3,j4) = AmpWork(j1,j2,j3,j4)
//$                 +CNorm*( Su1 +Su2 +Su3)
//$                 +CNorm*BornABCD(j1,j2,j3,j4)*( sC(1,Hel1)*sC(2,Hel2)*Fprop1   !
//$                                               +sD(1,Hel1)*sD(2,Hel2)*Fprop2 ) !
//$                 +CNorm*BornABCD(j1,j2,j3,j4)*sProd *(1d0 -svarX12/svarQ)  *N_IR !
//$                 +CNorm*BornABCD(j1,j2,j3,j4)*sProd *( svarX1/svarQ -1d0)  *N_IR !
//$                 +CNorm*BornABCD(j1,j2,j3,j4)*sProd *( svarX2/svarQ -1d0)  *N_IR !
//$                 +CNorm*BornABCD(j1,j2,j3,j4)*sProd                        *Y_IR !
///////////////////////////////////////////////////////////////////////////////////////
  dcmplx sProd = ( sC[0][Hel1]+sD[0][Hel1])*( sC[1][Hel2]+sD[1][Hel2]);
  AmpAdd(m_AmpExpo2, CNorm*(sC[0][Hel1]*sC[1][Hel2]*Fprop1+sD[0][Hel1]*sD[1][Hel2]*Fprop2), BornABCD);
  AmpAdd(m_AmpExpo2, CNorm*sProd*dcmplx(1.0 -svarX12/svarQ), BornABCD);
  AmpAdd(m_AmpExpo2, CNorm*sProd*dcmplx( svarX1/svarQ -1.0  +svarX2/svarQ -1.0), BornABCD);
  if( Y_IR) AmpAdd(m_AmpExpo2, CNorm*sProd, BornABCD);
//
}//GPS_HffPlus

void KKceex::GPS_HiiPlus(dcmplx CNorm, int KFini, int KFfin, TLorentzVector PX,
                         KKpart ph1, int Hel1, KKpart ph2, int Hel2){
/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
//   Genuine IR-finite non 1-photon amplitudes for ISR-ISR are added to AmpWork    //
//   That is for dip-switch Y_IR=0.                                                //
//   For Y_IR=1 IR-finite 1-photon contributions are calculated here as well       //
//                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////
//                                        |                                        //
//                              1         |          2                             //
//                              |         |X         |                             //
//      _       -b              |         |          |          a                  //
//      v  ------<------OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO----<----- u           //
//                                                                                 //
//                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////
  double ChaIni =  DB->Qf[ KFini];
  dcmplx gI = dcmplx(ChaIni*m_e_QED);
  dcmplx sA[2][2], sB[2][2];
  KKpart pA=m_p1;
  KKpart pB=m_p2;
  KKpart pC=m_p3;
  KKpart pD=m_p4;
//  sA(1,1)  = -gI*GPS_Sof1( 1,ph1,pA)
//  sA(2,1)  = -gI*GPS_Sof1( 1,ph2,pA)
//  sB(1,1)  =  gI*GPS_Sof1( 1,ph1,pB)
//  sB(2,1)  =  gI*GPS_Sof1( 1,ph2,pB)
  sA[0][0]  = -gI*GPS_Sof1( 1,ph1,pA);  // includes propagator
  sA[1][0]  = -gI*GPS_Sof1( 1,ph2,pA);
  sB[0][0]  =  gI*GPS_Sof1( 1,ph1,pB);
  sB[1][0]  =  gI*GPS_Sof1( 1,ph2,pB);
//  sA(1,2) = -DCONJG(sA(1,1))
//  sA(2,2) = -DCONJG(sA(2,1))
//  sB(1,2) = -DCONJG(sB(1,1))
//  sB(2,2) = -DCONJG(sB(2,1))
  sA[0][1] = -conj(sA[0][0]);
  sA[1][1] = -conj(sA[1][0]);
  sB[0][1] = -conj(sB[0][0]);
  sB[1][1] = -conj(sB[1][0]);
// Calculate Born spin amplitudes
  // CALL GPS_Born(KFi,KFf,PX, pA,mA,    pB,-mB,      pC,MC,   pD,-mD,   BornABCD) ! Standard
  // CALL GPS_Born(KFi,KFf,PX, ph1,mph,  pB,-mB,      pC,mC,   pD,-mD,   Born1BCD) ! A->1
  // CALL GPS_Born(KFi,KFf,PX, ph2,mph,  pB,-mB,      pC,mC,   pD,-mD,   Born2BCD) ! A->2
  // CALL GPS_Born(KFi,KFf,PX, pA,mA,    ph1,-mph,    pC,mC,   pD,-mD,   BornA1CD) ! B->1
  // CALL GPS_Born(KFi,KFf,PX, pA,mA,    ph2,-mph,    pC,mC,   pD,-mD,   BornA2CD) ! B->2
  // CALL GPS_Born(KFi,KFf,PX, ph1,mph,  ph2,-mph,    pC,mC,   pD,-mD,   Born12CD) ! A->1,B->2
  // CALL GPS_Born(KFi,KFf,PX, ph2,mph,  ph1,-mph,    pC,mC,   pD,-mD,   Born21CD);  // A->2,B->1
  double CosThe=0.0;
  KKpart ph1r = ph1; ph1r.C=-1;
  KKpart ph2r = ph2; ph2r.C=-1;
// dressed Born spin amplitudes
  KKcmplx4 BornABCD, Born1BCD, BornA1CD, Born2BCD, BornA2CD,Born12CD, Born21CD;
  Born(KFini,KFfin, PX, CosThe, m_p1, m_p2, m_p3, m_p4, BornABCD); //Standard
  Born(KFini,KFfin, PX, CosThe,  ph1, m_p2, m_p3, m_p4, Born1BCD); // A->1
  Born(KFini,KFfin, PX, CosThe,  ph2, m_p2, m_p3, m_p4, Born2BCD); // A->2
  Born(KFini,KFfin, PX, CosThe, m_p1, ph1r, m_p3, m_p4, BornA1CD); // B->1
  Born(KFini,KFfin, PX, CosThe, m_r1, ph2r, m_p3, m_p4, BornA2CD); // B->2
  Born(KFini,KFfin, PX, CosThe,  ph1, ph2r, m_p3, m_p4, Born12CD); // A->1,B->2
  Born(KFini,KFfin, PX, CosThe,  ph2, ph1r, m_p3, m_p4, Born21CD); // A->2,B->1
//Fermion propagarotors ini1
  double prA1= 1.0/(pA*ph1)/2.0;
  double prB1=-1.0/(pB*ph1)/2.0;
//Fermion propagarotors ini2
  double prA2= 1.0/(pA*ph2)/2.0;
  double prB2=-1.0/(pB*ph2)/2.0;
//DOUBLE propagators
  double prA12= 1.0/( pA*ph1 + pA*ph2 - ph1*ph2)/2.0;
  double prB12=-1.0/( pB*ph1 + pB*ph2 - ph1*ph2)/2.0;
  double Fprop1=(1.0/prA1+1.0/prA2)*prA12-1.0;
  double Fprop2=(1.0/prB1+1.0/prB2)*prB12-1.0;
  //
  KKcmplx2 U11a, Vb11, U21a, Vb12, Ua12, V21b, U112, V211, U212, V212;
  int Sig = 1-2*Hel1;
  GPS_MakeU( gI, ph1,Sig,  ph1, pA,   U11a);   // <1|[1]|a>
  GPS_MakeV( gI, ph1,Sig,  pB,  ph1,  Vb11);   // <b|{1}|1>
// falSe second
  GPS_MakeU( gI, ph1,Sig,  ph2, pA,   U21a);   // <2|[1]|a>
  GPS_MakeV( gI, ph1,Sig,  pB,  ph2,  Vb12);   // <b|{1}|2>
// reverse order
  GPS_MakeU( gI, ph1,Sig,  pA,  ph2,  Ua12);   // <a|[1]|2>
  GPS_MakeV( gI, ph1,Sig,  ph2, pB,   V21b);   // <2|{1}|b>
// for the case when there was ph2 first xk-xk term
  GPS_MakeU( gI, ph1,Sig,  ph1, ph2,  U112);   // <1|[1]|2>
  GPS_MakeV( gI, ph1,Sig,  ph2, ph1,  V211);   // <2|{1}|1>
// for the case when there was ph2 first xk-xk term
  GPS_MakeU( gI, ph1,Sig,  ph2, ph2,  U212);   // <2|[1]|2>
  GPS_MakeV( gI, ph1,Sig,  ph2, ph2,  V212);   // <2|{1}|2>
//
  KKcmplx2 U22a,Vb22, U12a, Vb21, Ua21,V12b,U221,V122,U121,V121;
  Sig = 1-2*Hel2;
  GPS_MakeU( gI, ph2,Sig,  ph2, pA,   U22a);   // <2|[2]|a>
  GPS_MakeV( gI, ph2,Sig,  pB,  ph2,  Vb22);   // <b|{2}|2>
// falSe second
  GPS_MakeU( gI, ph2,Sig,  ph1, pA,   U12a);  // <1|[2]|a>
  GPS_MakeV( gI, ph2,Sig,  pB,  ph1,  Vb21);  // <b|{2}|1>
// reverse order
  GPS_MakeU( gI, ph2,Sig,  pA,  ph1,  Ua21);  // <a|[2]|1>
  GPS_MakeV( gI, ph2,Sig,  ph1, pB,   V12b);  // <1|{2}|b>
// for the case when there was ph1 first xk-xk term
  GPS_MakeU( gI, ph2,Sig,  ph2, ph1,  U221);  // <2|[2]|1>
  GPS_MakeV( gI, ph2,Sig,  ph1, ph2,  V122);  // <1|{2}|2>
// for the case when there was ph1 first xk-xk term
  GPS_MakeU( gI, ph2,Sig,  ph1, ph1,  U121);  // <1|[2]|1>
  GPS_MakeV( gI, ph2,Sig,  ph1, ph1,  V121);  // <1|{2}|1>

  int N_IR = 1;  // legacy debugging variable
  int Y_IR = 1-N_IR;
  /////////////////////////////////////////////////////////////////////////////////////
  //                    |                  1               2                         //
  //                    |X                 |               |                         //
  //      _       -b    |      b+m-1-2     |      a+m-2    |      a                  //
  //      v  ------<----O--------<---------U--------<------S------<----- u           //
  /////////////////////////////////////////////////////////////////////////////////////
        //Su1=Su1+Born1BCD(j,j2,j3,j4) *U11a(j,j1)*(prA12-prA1*N_IR)*sA(2,Hel2) !<b|X|1[1]a(2)|a>
        //Su1=Su1+Born2BCD(j,j2,j3,j4) *U21a(j,j1)*   prA12         *sA(2,Hel2) !<b|X|2[1]a(2)|a>
  AmpAddI(m_AmpExpo2, CNorm*dcmplx(prA12-prA1*N_IR)*sA[1][Hel2], Born1BCD, U11a);
  AmpAddI(m_AmpExpo2, CNorm*dcmplx(prA12          )*sA[1][Hel2], Born2BCD, U21a);
  /////////////////////////////////////////////////////////////////////////////////////
  //                    |                  2               1                         //
  //                    |X                 |               |                         //
  //      _       -b    |      b+m-1-2     |      a+m-1    |      a                  //
  //      v  ------<----O--------<---------U-------<-------S------<----- u           //
  /////////////////////////////////////////////////////////////////////////////////////
//   Su1=Su1+Born2BCD(j,j2,j3,j4) *U22a(j,j1)*(prA12-prA2*N_IR)*sA(1,Hel1) !<b|X|2[2]a(1)|a>
//   Su1=Su1+Born1BCD(j,j2,j3,j4) *U12a(j,j1)* prA12           *sA(1,Hel1) !<b|X|1[2]a(1)|a>
  AmpAddI(m_AmpExpo2, CNorm*dcmplx(prA12-prA2*N_IR)*sA[0][Hel1],Born2BCD,U22a);
  AmpAddI(m_AmpExpo2, CNorm*dcmplx(prA12          )*sA[0][Hel1],Born1BCD,U12a);
  /////////////////////////////////////////////////////////////////////////////////////
  //                    2                  1               |                         //
  //                    |                  |               |X                        //
  //      _       -b    |     -b+m+2       |    -b+m+1+2   |      a                  //
  //      v  ------<----S--------<---------V--------<------O------<----- u           //
  /////////////////////////////////////////////////////////////////////////////////////
//   Su1=Su1+( sB(2,Hel2)*(-prB1*N_IR+prB12))*Vb11(j2,j)*BornA1CD(j1,j,j3,j4)!<b|(2)b[1]1|X|a>
//   Su1=Su1+( sB(2,Hel2))           *prB12  *Vb12(j2,j)*BornA2CD(j1,j,j3,j4)!<b|(2)b[1]2|X|a>
  AmpAddI(m_AmpExpo2, CNorm*sB[1][Hel2]*dcmplx(-prB1*N_IR+prB12),Vb11,BornA1CD);
  AmpAddI(m_AmpExpo2, CNorm*sB[1][Hel2]*dcmplx(           prB12),Vb12,BornA2CD);
  /////////////////////////////////////////////////////////////////////////////////////
  //                    1                  2               |                         //
  //                    |                  |               |X                        //
  //      _       -b    |     -b+m+1       |    -b+m+1+2   |      a                  //
  //      v  ------<----S--------<---------V--------<------O------<----- u           //
  /////////////////////////////////////////////////////////////////////////////////////
//   Su1=Su1 +(sB(1,Hel1)*(-prB2*N_IR+prB12))*Vb22(j2,j)*BornA2CD(j1,j,j3,j4)!<b|(1)b[2]2|X|a>
//   Su1=Su1 +(sB(1,Hel1))           *prB12  *Vb21(j2,j)*BornA1CD(j1,j,j3,j4)!<b|(1)b[2]1|X|a>
  AmpAddI(m_AmpExpo2, CNorm*sB[0][Hel1]*dcmplx(-prB2*N_IR+prB12),Vb22,BornA2CD);
  AmpAddI(m_AmpExpo2, CNorm*sB[0][Hel1]*dcmplx(           prB12),Vb21,BornA1CD);
  /////////////////////////////////////////////////////////////////////////////////////
  //                    2                  |               1                         //
  //                    |                  |X              |                         //
  //      _       -b    |     -b+m+2       |     a+m-1     |      a                  //
  //      v  ------<----S--------<---------O--------<------U------<----- u           //
  /////////////////////////////////////////////////////////////////////////////////////
//   Su1=Su1 +sB(2,Hel2)    *Born1BCD(j,j2,j3,j4) *prA1*U11a(j,j1) *Y_IR !<b|(2)2|X|1[1]|a>
  if(Y_IR ) AmpAddI(m_AmpExpo2, CNorm*sB[1][Hel2]*dcmplx( prA1 ),Born1BCD, U11a);
  /////////////////////////////////////////////////////////////////////////////////////
  //                    1                  |               2                         //
  //                    |                  |X              |                         //
  //      _       -b    |     -b+m+1       |     a+m-2     |      a                  //
  //      v  ------<----S--------<---------O--------<------U------<----- u           //
  /////////////////////////////////////////////////////////////////////////////////////
//   Su1=Su1 +sB(1,Hel1)    *Born2BCD(j,j2,j3,j4) *prA2*U22a(j,j1) *Y_IR !<b|(1)1|X|2[2]|a>
  if(Y_IR ) AmpAddI(m_AmpExpo2, CNorm*sB[0][Hel1]*dcmplx( prA2 ), Born2BCD, U22a);
//  /////////////////////////////////////////////////////////////////////////////////////
  //                    1                  |               2                         //
  //                    |                  |X              |                         //
  //      _       -b    |     -b+m+1       |     a+m-2     |      a                  //
  //      v  ------<----V--------<---------O--------<------S------<----- u           //
  /////////////////////////////////////////////////////////////////////////////////////
//   Su1=Su1 +Vb11(j2,j)*prB1 *BornA1CD(j1,j,j3,j4) *sA(2,Hel2)    *Y_IR !<b|[1]1|X|a[2]|a>
  if(Y_IR ) AmpAddI(m_AmpExpo2, CNorm* dcmplx( prB1 )*sA[1][Hel2], Vb11, BornA1CD);
  /////////////////////////////////////////////////////////////////////////////////////
  //                    2                  |               1                         //
  //                    |                  |X              |                         //
  //      _       -b    |     -b+m+2       |     a+m-1     |      a                  //
  //      v  ------<----V--------<---------O--------<------S------<----- u           //
  /////////////////////////////////////////////////////////////////////////////////////
//   Su1=Su1 +Vb22(j2,j)*prB2 *BornA2CD(j1,j,j3,j4) *sA(1,Hel1)    *Y_IR !<b|[2]2|X|a(1)|a>
  if(Y_IR ) AmpAddI(m_AmpExpo2, CNorm* dcmplx( prB2 )*sA[0][Hel1], Vb22, BornA2CD);


  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////
  //                    2                  |               1                         //
  //                    |                  |X              |                         //
  //      _       -b    |     -b+m+2       |     a+m-1     |      a                  //
  //      v  ------<----*--------<---------O--------<------O------<----- u           //
  /////////////////////////////////////////////////////////////////////////////////////
//  Amp4Zer(m_AmpTemp);
//  AmpAddI(m_AmpTemp,         dcmplx(prB2), Vb22, Born12CD);
//  AmpAddI(m_AmpExpo2, CNorm* dcmplx(prA1), m_AmpTemp, U11a);
// Su2=Su2  +Vb22(j2,l)*prB2  *Born12CD(j,l,j3,j4 )  *U11a(j,j1)*prA1 ! <b|[2]2|X|1[1]|a>
  AmpAddI(m_AmpExpo2, CNorm*dcmplx(prB2*prA1), Vb22, Born12CD, U11a);
  /////////////////////////////////////////////////////////////////////////////////////
  //                    1                  |               2                         //
  //                    |                  |X              |                         //
  //      _       -b    |     -b+m+1       |     a+m-2     |      a                  //
  //      v  ------<----*--------<---------O--------<------O------<----- u           //
  /////////////////////////////////////////////////////////////////////////////////////
      // Su2=Su2  +Vb11(j2,l)*prB1  *Born21CD(j,l,j3,j4 )  *U22a(j,j1)*prA2 ! <b|[1]1|X|2[2]|a>
  AmpAddI(m_AmpExpo2, CNorm*dcmplx(prB1*prA2), Vb11, Born21CD, U22a);

  /////////////////////////////////////////////////////////////////////////////////////
  //                    |                  2               1                         //
  //                    |X                 |               |                         //
  //      _       -b    |      b+m-1-2     |      a+m-1    |      a                  //
  //      v  ------<----O--------<---------O--------<------*------<----- u           //
  /////////////////////////////////////////////////////////////////////////////////////
    //   Su2=Su2 +Born1BCD(j,j2,j3,j4) *U121(j,l)*prA12 *U11a(l,j1)*prA1  ! <b|X|1[2]1(1)|a>
    //   Su2=Su2 +Born2BCD(j,j2,j3,j4) *U221(j,l)*prA12 *U11a(l,j1)*prA1  ! <b|X|2[2]1(1)|a>
    //   Su2=Su2 -BornABCD(j,j2,j3,j4) *Ua21(j,l)*prA12 *U11a(l,j1)*prA1  ! <b|X|a[2]1(1)|a>
  AmpAddI(m_AmpExpo2, CNorm*dcmplx( prA12*prA1), Born1BCD, U121, U11a);
  AmpAddI(m_AmpExpo2, CNorm*dcmplx( prA12*prA1), Born2BCD, U221, U11a);
  AmpAddI(m_AmpExpo2, CNorm*dcmplx(-prA12*prA1), BornABCD, Ua21, U11a);
  /////////////////////////////////////////////////////////////////////////////////////
  //                    |                  1               2                         //
  //                    |X                 |               |                         //
  //      _       -b    |      b+m-1-2     |      a+m-2    |      a                  //
  //      v  ------<----O--------<---------O--------<------*------<----- u           //
  /////////////////////////////////////////////////////////////////////////////////////
  //    Su2=Su2 +Born2BCD(j,j2,j3,j4) *U212(j,l)*prA12  *U22a(l,j1)*prA2 ! <b|X|2[1]2(2)|a>
  //    Su2=Su2 +Born1BCD(j,j2,j3,j4) *U112(j,l)*prA12  *U22a(l,j1)*prA2 ! <b|X|1[1]2(2)|a>
  //    Su2=Su2 -BornABCD(j,j2,j3,j4) *Ua12(j,l)*prA12  *U22a(l,j1)*prA2 ! <b|X|a[1]2(2)|a>
  AmpAddI(m_AmpExpo2, CNorm*dcmplx( prA12*prA2), Born2BCD, U212, U22a);
  AmpAddI(m_AmpExpo2, CNorm*dcmplx( prA12*prA2), Born1BCD, U112, U22a);
  AmpAddI(m_AmpExpo2, CNorm*dcmplx(-prA12*prA2), BornABCD, Ua12, U22a);

  /////////////////////////////////////////////////////////////////////////////////////
  //                    1                  2               |                         //
  //                    |                  |               |X                        //
  //      _       -b    |     -b+m+1       |    -b+m+1+2   |      a                  //
  //      v  ------<----*--------<---------O--------<------O------<----- u           //
  /////////////////////////////////////////////////////////////////////////////////////
  //      Su2=Su2 +Vb11(j2,l)*prB1 *V121(l,j)*prB12 *BornA1CD(j1,j,j3,j4) ! <b|(1)1[2]1|X|a>
  //      Su2=Su2 +Vb11(j2,l)*prB1 *V122(l,j)*prB12 *BornA2CD(j1,j,j3,j4) ! <b|(1)1[2]2|X|a>
  //      Su2=Su2 -Vb11(j2,l)*prB1 *V12b(l,j)*prB12 *BornABCD(j1,j,j3,j4) ! <b|(1)1[2]b|X|a>
  AmpAddI(m_AmpExpo2, CNorm*dcmplx( prB1*prB12), Vb11, V121, BornA1CD);
  AmpAddI(m_AmpExpo2, CNorm*dcmplx( prB1*prB12), Vb11, V122, BornA2CD);
  AmpAddI(m_AmpExpo2, CNorm*dcmplx(-prB1*prB12), Vb11, V12b, BornABCD);
  /////////////////////////////////////////////////////////////////////////////////////
  //                    2                  1               |                         //
  //                    |                  |               |X                        //
  //      _       -b    |     -b+m+2       |    -b+m+1+2   |      a                  //
  //      v  ------<----*--------<---------O--------<------O------<----- u           //
  /////////////////////////////////////////////////////////////////////////////////////
  //      Su2=Su2 +Vb22(j2,l)*prB2 *V212(l,j)*prB12  *BornA2CD(j1,j,j3,j4) ! <b|(2)2[1]2|X|a>
  //      Su2=Su2 +Vb22(j2,l)*prB2 *V211(l,j)*prB12  *BornA1CD(j1,j,j3,j4) ! <b|(2)2[1]1|X|a>
  //      Su2=Su2 -Vb22(j2,l)*prB2 *V21b(l,j)*prB12  *BornABCD(j1,j,j3,j4) ! <b|(2)2[2]b|X|a>
  AmpAddI(m_AmpExpo2, CNorm*dcmplx( prB2*prB12), Vb22, V212, BornA2CD);
  AmpAddI(m_AmpExpo2, CNorm*dcmplx( prB2*prB12), Vb22, V211, BornA1CD);
  AmpAddI(m_AmpExpo2, CNorm*dcmplx(-prB2*prB12), Vb22, V21b, BornABCD);
/////////////////////////////////////////////
//  sProd = ( sA(1,Hel1)+sB(1,Hel1))*( sA(2,Hel2)+sB(2,Hel2))
//  AmpWork(j1,j2,j3,j4) =AmpWork(j1,j2,j3,j4)
//$                 +CNorm*( Su1+Su2 )
//$                 +CNorm*BornABCD(j1,j2,j3,j4)*( sA(1,Hel1)*sA(2,Hel2)*Fprop1  !
//$                                               +sB(1,Hel1)*sB(2,Hel2)*Fprop2) !
//$                 +CNorm*BornABCD(j1,j2,j3,j4)* sProd                   *Y_IR  !
//////////////////////////////////////////////////////
  dcmplx sProd = ( sA[0][Hel1]+sB[0][Hel1])*( sA[1][Hel2]+sB[1][Hel2]);
  AmpAdd(m_AmpExpo2, CNorm*(sA[0][Hel1]*sA[1][Hel2]*Fprop1+sB[0][Hel1]*sB[1][Hel2]*Fprop2), BornABCD);
  if( Y_IR) AmpAdd(m_AmpExpo2, CNorm*sProd, BornABCD);
}//GPS_HiiPlus

void KKceex::GPS_HiiPlusW(dcmplx CNorm, int KFini, int KFfin, TLorentzVector PX,
                         KKpart ph1, int Hel1, KKpart ph2, int Hel2){
/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
//                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////
if(abs(KFfin) != 12)  return;
//depending on Level some groups of terms are calculated. only Level=0 is used now.
int Level = 0;
int ICOL,ICOL1,NICOL1,NICOL;
if(        Level == 0){
  ICOL =1;   ICOL1=-1;  NICOL1=-1;   NICOL=-1;
} else if( Level == 1) {
  ICOL =0;   ICOL1=1;   NICOL1=0;    NICOL=0;
} else if( Level == 2) {
  ICOL =0;   ICOL1=0;   NICOL1=1;    NICOL=0;
} else {
  ICOL =0;   ICOL1=-1;  NICOL1=-1 ;  NICOL=-1;
}
// technical switches for some tests. Must be all 1 or, better know what you do.
// numerical stability not yet checked, problems like with single brem. possible.
int  I71=ICOL1;  //   ! single additional loop  1st photon ir-factor     >  ggchkok
int  I72=ICOL1;  //   ! single additional loop  2nd photon ir-factor     >  ggchkok
int  I8= ICOL;   //   ! double additional loop                           >  ggchkok
int  IA= ICOL;   //   ! doble emission from single leg, but selfcanc     >  ggchkok
int  I9= NICOL;  //   ! no additional loop                               >  #######ggchkok
int  I9B=NICOL;  //   ! no additional loop                               >  #######ggchkok
int  I9X=ICOL;   //   !NICOL1    ! first finite upperline U, second for cancel      >  ggchkok
int  I9Y=ICOL;   //   !NICOL1    ! first finite upperline V, second for cancel      >  ggchkok
int  I9Z=ICOL;   //   !NICOL1    ! second finite upperline U, second for cancel     >  ggchkok
int  I9T=ICOL;   //   !NICOL1    ! second finite upperline V, second for cancel     >  ggchkok
int  IVI=ICOL;   //   !NICOL   ! simplified IV i.e. 'true double infrared'        >  ggchkok
int  IV2=ICOL;   //   !ICOL1 ! ICOL1   ! double photon rest (two from 1 leg) k1*k2 ....   >  ggchkok
int  IV1=ICOL;   //   !ICOL1 ! ICOL1   ! double photon rest (two from 1 leg) k1*k2 ....   >  ggchkok
int  I10=NICOL;  //   ! four boson (rest)                                >  #######ggchkok
int  I9s1=ICOL1; //   ! first soft second gaugeinv nontrivial WWgamma   >  ggchkok
int  I9s2=ICOL1; //   ! second soft first gaugeinv nontrivial WWgamma   >  ggchkok
int  I71b=0;     //   !(NICOL)  ! non-ir remnant of I71     >  outed to I9B (ggchkok)
int  I72b=0;     //   !(NICOL)  ! non-ir remnant of I72     >  outed to I9B (ggchkok)
//-----
double ChaIni =  DB->Qf[ KFini];
dcmplx gI = dcmplx(ChaIni*m_e_QED);
//[[[[[[[[[[[[[[[[[[[[[[[[[[[
//(*m_Out)<< "========================================================================================================"<<endl;
//(*m_Out)<< "========================================KKceex::GPS_HiiPlusW============================================="<<endl;
//]]]]]]]]]]]]]]]]]]]]]]]]]]]
//double Qe  = DB->Qf[ KFini];
KKpart pA= m_p1;
KKpart pB= m_p2;
KKpart pC= m_p3;
KKpart pD= m_p4;
//----------------
double s0,t0;
int IFONE;
//Where is dominant other photon?  It is instead of reduction procedure
if( (m_p3.P[3]+m_p4.P[3])*m_p1.P[3] < 0.0) IFONE=1; else IFONE=0;
if(IFONE) { s0=  2*(m_p3*m_p4); t0= -2*(m_p4*m_p2);
} else {    s0=  2*(m_p3*m_p4); t0= -2*(m_p3*m_p1);
}//if
//----------------
if( (ph1.P[3]+ph2.P[3]+m_p3.P[3]+m_p4.P[3]) * m_p1.P[3] < 0.0) {
 IFONE= 1;      // We assume all extra photons were emitted from p1
} else {
 IFONE= 0;      // We assume all extra photons were emitted from p2
}
double s,t,s1,t1,s2,t2,s12,t12;
if(IFONE){
s =(pC+pD+ph1+ph2).M2();  t =(pD-pB+ph1+ph2).M2();
s1=(pC+pD+ph1+ph2).M2();  t1=(pD-pB+ph2).M2();
s2=(pC+pD+ph1+ph2).M2();  t2=(pD-pB+ph1).M2();
s12=(pC+pD+ph1+ph2).M2(); t12=(pD-pB).M2();
} else {
s =(pC+pD+ph1+ph2).M2();  t =(pC-pA).M2();
s1=(pC+pD+ph1+ph2).M2();  t1=(pC-pA+ph1).M2();
s2=(pC+pD+ph1+ph2).M2();  t2=(pC-pA+ph2).M2();
s12=(pC+pD+ph1+ph2).M2(); t12=(pC-pA+ph1+ph2).M2();
}
//[[[[[[[[[[[[[[[[[[[[[
//(*m_Out)<< "KKceex::HiniPlusW: IFONE= "<< IFONE<<"  s0="<<s0<<"  t0="<<t0<<endl;
//(*m_Out)<< "KKceex::HiniPlusW: s="<<s<<"  t="<<t<<" s1="<<s1<<"  t1="  <<t1<<" s2="<<s2<<"  t2="<<t2 <<endl;
//(*m_Out)<< "KKceex::HiniPlusW: s12="<<s12<<"  t12="<<t12<<endl;
//]]]]]]]]]]]]]]]]]]]]]
//CALL GPS_EWFFactW(KFi,KFf,s0,t0,PropW0,WVPi0)
//CALL GPS_EWFFactW(KFi,KFf,s,t,PropW,WVPi)
//CALL GPS_EWFFactW(KFi,KFf,s1,t1,PropW1,WVPi1)
//CALL GPS_EWFFactW(KFi,KFf,s2,t2,PropW2,WVPi2)
//CALL GPS_EWFFactW(KFi,KFf,s12,t12,PropW12,WVPi12)
// all W propagators ... some simplifications on vac pol set
dcmplx PropW0,PropW,PropW1,PropW2,PropW12;
dcmplx WVPi0, WVPi, WVPi1, WVPi2, WVPi12;
GPS_EWFFactW(KFini,KFfin,s0,t0,PropW0,WVPi0);
GPS_EWFFactW(KFini,KFfin,s,t,PropW,WVPi);
GPS_EWFFactW(KFini,KFfin,s1,t1,PropW1,WVPi1);
GPS_EWFFactW(KFini,KFfin,s2,t2,PropW2,WVPi2);
GPS_EWFFactW(KFini,KFfin,s12,t12,PropW12,WVPi12);
WVPi  =1.0;            //!!!  WVPi0
WVPi1 =1.0;            //!!!  WVPi0
WVPi2 =1.0;            //!!!  WVPi0
WVPi12=1.0;            //!!!  WVPi0
dcmplx CPF0 =PropW0 *1.0;     //!!!   *WVPi0
dcmplx CPF  =PropW  *WVPi;
dcmplx CPF1 =PropW1 *WVPi1;
dcmplx CPF2 =PropW2 *WVPi2;
dcmplx CPF12=PropW12*WVPi12;
//[[[[[[[[[[[[[[[[[[[[[
//(*m_Out)<< "PropW0="<<PropW0<<"  PropW="   <<PropW<<endl;
//(*m_Out)<< "PropW1="<<PropW1<<"  PropW2="  <<PropW2<<" PropW12="<<PropW12<<endl;
//(*m_Out)<< "WVPi0="<<WVPi0<<"  WVPi="   <<WVPi<<endl;
//(*m_Out)<< "WVPi1="<<WVPi1<<"  WVPi2="  <<WVPi2<<" WVPi12="<<WVPi12<<endl;
//]]]]]]]]]]]]]]]]]]]]]
//sA(1,1)  = -gI*GPS_Sof1( 1,ph1,pA)
//sA(2,1)  = -gI*GPS_Sof1( 1,ph2,pA)
//sB(1,1)  =  gI*GPS_Sof1( 1,ph1,pB)
//sB(2,1)  =  gI*GPS_Sof1( 1,ph2,pB)
dcmplx sA[2][2],sB[2][2]; // soft factors
sA[0][0]  = -gI*GPS_Sof1( 1,ph1,m_p1);
sA[1][0]  = -gI*GPS_Sof1( 1,ph2,m_p1);
sB[0][0]  =  gI*GPS_Sof1( 1,ph1,m_p2);
sB[1][0]  =  gI*GPS_Sof1( 1,ph2,m_p2);
sA[0][1] = -conj(sA[0][0]);
sA[1][1] = -conj(sA[1][0]);
sB[0][1] = -conj(sB[0][0]);
sB[1][1] = -conj(sB[1][0]);
//IF (IFONE) THEN
//  EpsDot1(1) = -gI*(+GPS_Sof1bx( 1,ph1,pB,mB)-GPS_Sof1bx( 1,ph1,pD,mD)-GPS_Sof1bx( 1,ph1,ph2,mph))
//  EpsDot12(1)= -gI*(+GPS_Sof1bx( 1,ph1,pB,mB)-GPS_Sof1bx( 1,ph1,pD,mD))
//  EpsDot2(1) = -gI*(+GPS_Sof1bx( 1,ph2,pB,mB)-GPS_Sof1bx( 1,ph2,pD,mD)-GPS_Sof1bx( 1,ph2,ph1,mph))
//  EpsDot21(1)= -gI*( GPS_Sof1bx( 1,ph2,pB,mB)-GPS_Sof1bx( 1,ph2,pD,mD))
//ELSE
//  EpsDot1(1) = -gI*(-GPS_Sof1bx( 1,ph1,pA,mA)+GPS_Sof1bx( 1,ph1,pC,mC))
//  EpsDot12(1)= -gI*(-GPS_Sof1bx( 1,ph1,pA,mA)+GPS_Sof1bx( 1,ph1,pC,mC)+GPS_Sof1bx( 1,ph1,ph2,mph))
//  EpsDot2(1) = -gI*(-GPS_Sof1bx( 1,ph2,pA,mA)+GPS_Sof1bx( 1,ph2,pC,mC))
//  EpsDot21(1)= -gI*(-GPS_Sof1bx( 1,ph2,pA,mA)+GPS_Sof1bx( 1,ph2,pC,mC)+GPS_Sof1bx( 1,ph2,ph1,mph))
//ENDIF
dcmplx EpsDot1[2],EpsDot2[2],EpsDot12[2],EpsDot21[2]; // soft factors
if(IFONE){
  EpsDot1[0] = -gI*( GPS_Sof1x( 1,ph1,pB)-GPS_Sof1x( 1,ph1,pD)-GPS_Sof1x( 1,ph1,ph2));
  EpsDot12[0]= -gI*( GPS_Sof1x( 1,ph1,pB)-GPS_Sof1x( 1,ph1,pD));
  EpsDot2[0] = -gI*( GPS_Sof1x( 1,ph2,pB)-GPS_Sof1x( 1,ph2,pD)-GPS_Sof1x( 1,ph2,ph1));
  EpsDot21[0]= -gI*( GPS_Sof1x( 1,ph2,pB)-GPS_Sof1x( 1,ph2,pD));
} else {
  EpsDot1[0] = -gI*(-GPS_Sof1x( 1,ph1,pA)+GPS_Sof1x( 1,ph1,pC));
  EpsDot12[0]= -gI*(-GPS_Sof1x( 1,ph1,pA)+GPS_Sof1x( 1,ph1,pC)+GPS_Sof1x( 1,ph1,ph2));
  EpsDot2[0] = -gI*(-GPS_Sof1x( 1,ph2,pA)+GPS_Sof1x( 1,ph2,pC));
  EpsDot21[0]= -gI*(-GPS_Sof1x( 1,ph2,pA)+GPS_Sof1x( 1,ph2,pC)+GPS_Sof1x( 1,ph2,ph1));
}//
EpsDot1[1]  = -conj(EpsDot1[0]);
EpsDot12[1] = -conj(EpsDot12[0]);
EpsDot2[1]  = -conj(EpsDot2[0]);
EpsDot21[1] = -conj(EpsDot21[0]);
//[[[[[[[[[[[[[[[[[[[[[
//(*m_Out)<< "sA(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<sA[j][k];(*m_Out)<<endl;
//(*m_Out)<< "sB(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<sB[j][k];(*m_Out)<<endl;
//(*m_Out)<< " EpsDot1(*,*)=";for(int j=0;j<=1;j++) (*m_Out)<<"  "<<SW208<<EpsDot1[j];(*m_Out)<<endl;
//(*m_Out)<< " EpsDot2(*,*)=";for(int j=0;j<=1;j++) (*m_Out)<<"  "<<SW208<<EpsDot2[j];(*m_Out)<<endl;
//(*m_Out)<< "EpsDot12(*,*)=";for(int j=0;j<=1;j++) (*m_Out)<<"  "<<SW208<<EpsDot12[j];(*m_Out)<<endl;
//(*m_Out)<< "EpsDot21(*,*)=";for(int j=0;j<=1;j++) (*m_Out)<<"  "<<SW208<<EpsDot21[j];(*m_Out)<<endl;
//]]]]]]]]]]]]]]]]]]]]]
// photon polarization 4-vectors calculated explicitelly ...
//      Sig = 3-2*Hel1
//      CALL GPS_Make_eps(ph1,Sig,eps1)
//      Sig = 3-2*Hel2
//      CALL GPS_Make_eps(ph2,Sig,eps2)
dcmplx eps1[4],eps2[4];
int Sig1 = 1-2*Hel1, Sig2 = 1-2*Hel2;
GPS_Make_eps(ph1, Sig1, eps1);
GPS_Make_eps(ph2, Sig2, eps2);
dcmplx eps1D2=eps1[0]*eps2[0]-eps1[3]*eps2[3]-eps1[2]*eps2[2]-eps1[1]*eps2[1];
dcmplx eps1pA=eps1[0]*pA[0]-eps1[3]*pA[3]-eps1[2]*pA[2]-eps1[1]*pA[1];
dcmplx eps2pA=eps2[0]*pA[0]-eps2[3]*pA[3]-eps2[2]*pA[2]-eps2[1]*pA[1];
dcmplx eps1pB=eps1[0]*pB[0]-eps1[3]*pB[3]-eps1[2]*pB[2]-eps1[1]*pB[1];
dcmplx eps2pB=eps2[0]*pB[0]-eps2[3]*pB[3]-eps2[2]*pB[2]-eps2[1]*pB[1];
dcmplx eps1pC=eps1[0]*pC[0]-eps1[3]*pC[3]-eps1[2]*pC[2]-eps1[1]*pC[1];
dcmplx eps2pC=eps2[0]*pC[0]-eps2[3]*pC[3]-eps2[2]*pC[2]-eps2[1]*pC[1];
dcmplx eps1pD=eps1[0]*pD[0]-eps1[3]*pD[3]-eps1[2]*pD[2]-eps1[1]*pD[1];
dcmplx eps2pD=eps2[0]*pD[0]-eps2[3]*pD[3]-eps2[2]*pD[2]-eps2[1]*pD[1];
dcmplx eps1p2=eps1[0]*ph2[0]-eps1[3]*ph2[3]-eps1[2]*ph2[2]-eps1[1]*ph2[1];
dcmplx eps2p1=eps2[0]*ph1[0]-eps2[3]*ph1[3]-eps2[2]*ph1[2]-eps2[1]*ph1[1];
dcmplx  Cp2p1=ph2[0]*ph1[0]-ph2[3]*ph1[3]-ph2[2]*ph1[2]-ph2[1]*ph1[1];
//[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
//(*m_Out)<<"  eps1pA="  <<SW208<<eps1pA<<" eps2pA="<<SW208<<eps2pA<<endl;
//(*m_Out)<<"  eps1pB="  <<SW208<<eps1pB<<" eps2pB="<<SW208<<eps2pB<<endl;
//(*m_Out)<<"  eps1pC="  <<SW208<<eps1pC<<" eps2pC="<<SW208<<eps2pC<<endl;
//(*m_Out)<<"  eps1pD="  <<SW208<<eps1pD<<" eps2pD="<<SW208<<eps2pD<<endl;
//(*m_Out)<<"  eps1p2="  <<SW208<<eps1p2<<" eps2p1="<<SW208<<eps2p1<<endl;
//(*m_Out)<<"  eps1D2="  <<SW208<<eps1D2<<"  Cp2p1="<<SW208<< Cp2p1<<endl;
//]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]

double CosThe=0.0;
KKpart ph1r = ph1; ph1r.C=-1;
KKpart ph2r = ph2; ph2r.C=-1;
//CALL GPS_BornWPlus(1,0,KFi,KFf,s0,t0,u0, pA,mA,    pB,-mB,      pC,MC,   pD,-mD,   BornABCD) ! Standard
//CALL GPS_BornWPlus(1,0,KFi,KFf,s0,t0,u0, ph1,mph,  pB,-mB,      pC,mC,   pD,-mD,   Born1BCD) ! A->1
//CALL GPS_BornWPlus(1,0,KFi,KFf,s0,t0,u0, pA,mA,    ph1,-mph,    pC,mC,   pD,-mD,   BornA1CD) ! B->1
//CALL GPS_BornWPlus(1,0,KFi,KFf,s0,t0,u0, ph2,mph,  pB,-mB,      pC,mC,   pD,-mD,   Born2BCD) ! A->2
//CALL GPS_BornWPlus(1,0,KFi,KFf,s0,t0,u0, pA,mA,    ph2,-mph,    pC,mC,   pD,-mD,   BornA2CD) ! B->2
//CALL GPS_BornWPlus(1,0,KFi,KFf,s0,t0,u0, ph1,mph,  ph2,-mph,    pC,mC,   pD,-mD,   Born12CD) ! A->1,B->2
//CALL GPS_BornWPlus(1,0,KFi,KFf,s0,t0,u0, ph2,mph,  ph1,-mph,    pC,mC,   pD,-mD,   Born21CD) ! A->2,B->1
// dressed Born spin amplitudes
KKcmplx4 BornABCD, Born1BCD, BornA1CD, Born2BCD, BornA2CD,Born12CD, Born21CD;
BornW(KFini,KFfin, PX, s0,t0, m_p1, m_p2, m_p3, m_p4, BornABCD); //Standard
BornW(KFini,KFfin, PX, s0,t0,  ph1, m_p2, m_p3, m_p4, Born1BCD); // A->1
BornW(KFini,KFfin, PX, s0,t0, m_p1, ph1r, m_p3, m_p4, BornA1CD); // B->1
BornW(KFini,KFfin, PX, s0,t0,  ph2, m_p2, m_p3, m_p4, Born2BCD); // A->2
//[[[BornW(KFini,KFfin, PX, s0,t0, m_r1, ph2r, m_p3, m_p4, BornA2CD); // B->2
BornW(KFini,KFfin, PX, s0,t0, m_p1, ph2r, m_p3, m_p4, BornA2CD); // B->2
BornW(KFini,KFfin, PX, s0,t0,  ph1, ph2r, m_p3, m_p4, Born12CD); // A->1,B->2
BornW(KFini,KFfin, PX, s0,t0,  ph2, ph1r, m_p3, m_p4, Born21CD); // A->2,B->1
//[[[[[[[[[[[[[[[[[[[[[[[[[[[*debug*
/*
(*m_Out)<<"--------------------------------------------------------------------------------------------"<<endl;
(*m_Out)<<" s0="<<s0<<" t0="<<t0<<endl;
for(int j1=0; j1<=1; j1++) for(int j2=0; j2<=1; j2++){ (*m_Out)<< "BornABCD("<<j1<<","<<j2<<",*,*)=  ";
     for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<BornABCD.m_A[j1][j2][j][k];(*m_Out)<<endl;}
for(int j1=0; j1<=1; j1++) for(int j2=0; j2<=1; j2++){ (*m_Out)<< "Born1BCD("<<j1<<","<<j2<<",*,*)=";
     for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<Born1BCD.m_A[j1][j2][j][k];(*m_Out)<<endl;}
for(int j1=0; j1<=1; j1++) for(int j2=0; j2<=1; j2++){ (*m_Out)<< "BornBornA1CD("<<j1<<","<<j2<<",*,*)=  ";
     for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<BornA1CD.m_A[j1][j2][j][k];(*m_Out)<<endl;}
(*m_Out)<<"--------------------------------------------------------------------------------------------"<<endl;
for(int j1=0; j1<=1; j1++) for(int j2=0; j2<=1; j2++){ (*m_Out)<< "Born2BCD("<<j1<<","<<j2<<",*,*)=";
     for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<Born2BCD.m_A[j1][j2][j][k];(*m_Out)<<endl;}
for(int j1=0; j1<=1; j1++) for(int j2=0; j2<=1; j2++){ (*m_Out)<< "BornA2CD("<<j1<<","<<j2<<",*,*)=  ";
     for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<BornA2CD.m_A[j1][j2][j][k];(*m_Out)<<endl;}
(*m_Out)<<"--------------------------------------------------------------------------------------------"<<endl;
for(int j1=0; j1<=1; j1++) for(int j2=0; j2<=1; j2++){ (*m_Out)<< "Born12CD("<<j1<<","<<j2<<",*,*)=";
     for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<Born12CD.m_A[j1][j2][j][k];(*m_Out)<<endl;}
for(int j1=0; j1<=1; j1++) for(int j2=0; j2<=1; j2++){ (*m_Out)<< "Born21CD("<<j1<<","<<j2<<",*,*)=  ";
     for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<Born21CD.m_A[j1][j2][j][k];(*m_Out)<<endl;}
*/
//]]]]]]]]]]]]]]]]]]]]]]]]]]]
// Fermion propagarotors ini1
double prA1= 1.0/(pA*ph1)/2.0;
double prB1=-1.0/(pB*ph1)/2.0;
// Fermion propagarotors ini2
double prA2= 1.0/(pA*ph2)/2.0;
double prB2=-1.0/(pB*ph2)/2.0;
// DOUBLE propagators
double prA12= 1.0/( pA*ph1 +pA*ph2 -ph1*ph2 )/2.0;
double prB12=-1.0/( pB*ph1 +pB*ph2 -ph1*ph2 )/2.0;
dcmplx Fprop1 =(1.0/prA1+1.0/prA2)*prA12*CPF12/CPF0-1.0;
dcmplx Fprop2 =(1.0/prB1+1.0/prB2)*prB12*CPF/CPF0-1.0;
dcmplx Fprop1B=CPF12/CPF0-1.0;
dcmplx Fprop2B=CPF/CPF0-1.0;
//[[[[[[[[[[[[[[[[[[[[[
//(*m_Out)<< "prA1= "<<prA1 <<" prB1 =" <<prB1<<" prA2="<<prA2<<" prB2="<<prB2<<endl;
//(*m_Out)<<" prA12="<<prA12<<" prB12="<<prB12<<endl;
//(*m_Out)<<" Fprop1="<<Fprop1<<" Fprop2="<<Fprop2<<endl;
//(*m_Out)<<" Fprop1B="<<Fprop1B<<" Fprop2B="<<Fprop2B<<endl;
//---------------------------------------------------
//Sig = 3-2*Hel1
//   CALL GPS_MatrU( gI, ph1,Sig,  ph1,mph, pA,mA,     U11a) ! <1|[1]|a>
//   CALL GPS_MatrV( gI, ph1,Sig,  pB,mB,   ph1,mph,   Vb11) ! <b|{1}|1>
//* falSe second
//   CALL GPS_MatrU( gI, ph1,Sig,  ph2,mph, pA,mA,     U21a) ! <2|[1]|a>
//   CALL GPS_MatrV( gI, ph1,Sig,  pB,mB,   ph2,mph,   Vb12) ! <b|{1}|2>
//* reverse order
//   CALL GPS_MatrU( gI, ph1,Sig,  pA,mA,   ph2,mph,   Ua12) ! <a|[1]|2>
//   CALL GPS_MatrV( gI, ph1,Sig,  ph2,mph, pB,mB,     V21b) ! <2|{1}|b>
//* for the case when there was ph2 first xk-xk term
//   CALL GPS_MatrU( gI, ph1,Sig,  ph1,mph, ph2,mph,   U112) ! <1|[1]|2>
//   CALL GPS_MatrV( gI, ph1,Sig,  ph2,mph, ph1,mph,   V211) ! <2|{1}|1>
//* for the case when there was ph2 first xk-xk term
//   CALL GPS_MatrU( gI, ph1,Sig,  ph2,mph, ph2,mph,   U212) ! <2|[1]|2>
//   CALL GPS_MatrV( gI, ph1,Sig,  ph2,mph, ph2,mph,   V212) ! <2|{1}|2>
KKcmplx2 U11a,Vb11,Vb12,U21a,Ua12,V21b,U112,V211,U212,V212;
GPS_MakeU( gI, ph1,Sig1,  ph1, pA,   U11a); // <1|[1]|a>
GPS_MakeV( gI, ph1,Sig1,  pB,  ph1,  Vb11); // <b|{1}|1>
// falSe second
GPS_MakeU( gI, ph1,Sig1,  ph2, pA,   U21a); // <2|[1]|a>
GPS_MakeV( gI, ph1,Sig1,  pB,  ph2,  Vb12); // <b|{1}|2>
// reverse order
GPS_MakeU( gI, ph1,Sig1,  pA,  ph2,  Ua12); // <a|[1]|2>
GPS_MakeV( gI, ph1,Sig1,  ph2, pB,   V21b); // <2|{1}|b>
// for the case when there was ph2 first xk-xk term
GPS_MakeU( gI, ph1,Sig1,  ph1, ph2,  U112); // <1|[1]|2>
GPS_MakeV( gI, ph1,Sig1,  ph2, ph1,  V211); // <2|{1}|1>
// for the case when there was ph2 first xk-xk term
GPS_MakeU( gI, ph1,Sig1,  ph2, ph2,  U212); // <2|[1]|2>
GPS_MakeV( gI, ph1,Sig1,  ph2, ph2,  V212); // <2|{1}|2>
//[[[[[[[[[[[[[[[[[[[[[
//(*m_Out)<< "U11a(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<U11a.m_A[j][k];(*m_Out)<<endl;
//(*m_Out)<< "Vb11(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<Vb11.m_A[j][k];(*m_Out)<<endl;
//
//(*m_Out)<< "U21a(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<U21a.m_A[j][k];(*m_Out)<<endl;
//(*m_Out)<< "Vb12(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<Vb12.m_A[j][k];(*m_Out)<<endl;
//
//(*m_Out)<< "Ua12(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<Ua12.m_A[j][k];(*m_Out)<<endl;
//(*m_Out)<< "V21b(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<V21b.m_A[j][k];(*m_Out)<<endl;
//
//(*m_Out)<< "U112(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<U112.m_A[j][k];(*m_Out)<<endl;
//(*m_Out)<< "V211(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<V211.m_A[j][k];(*m_Out)<<endl;
//
//(*m_Out)<< "U212(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<U212.m_A[j][k];(*m_Out)<<endl;
//(*m_Out)<< "V212(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<V212.m_A[j][k];(*m_Out)<<endl;
//]]]]]]]]]]]]]]]]]]]]]
//   CALL GPS_MatrU( gI, ph2,Sig,  ph2,mph,  pA,mA,    U22a) ! <2|[2]|a>
//   CALL GPS_MatrV( gI, ph2,Sig,  pB,mB,   ph2,mph,   Vb22) ! <b|{2}|2>
//* falSe second
//   CALL GPS_MatrU( gI, ph2,Sig,  ph1,mph,  pA,mA,    U12a) ! <1|[2]|a>
//   CALL GPS_MatrV( gI, ph2,Sig,  pB,mB,   ph1,mph,   Vb21) ! <b|{2}|1>
//* reverse order
//   CALL GPS_MatrU( gI, ph2,Sig,  pA,mA,   ph1,mph,   Ua21) ! <a|[2]|1>
//   CALL GPS_MatrV( gI, ph2,Sig,  ph1,mph, pB,mB,     V12b) ! <1|{2}|b>
//* for the case when there was ph1 first xk-xk term
//   CALL GPS_MatrU( gI, ph2,Sig,  ph2,mph, ph1,mph,   U221) ! <2|[2]|1>
//   CALL GPS_MatrV( gI, ph2,Sig,  ph1,mph, ph2,mph,   V122) ! <1|{2}|2>
//c for the case when there was ph1 first xk-xk term
//   CALL GPS_MatrU( gI, ph2,Sig,  ph1,mph, ph1,mph,   U121) ! <1|[2]|1>
//   CALL GPS_MatrV( gI, ph2,Sig,  ph1,mph, ph1,mph,   V121) ! <1|{2}|1>
KKcmplx2 U22a,Vb22, U12a,Vb21, Ua21,V12b, U221,V122, U121,V121;
GPS_MakeU( gI, ph2,Sig2,  ph2, pA,   U22a); // <1|[1]|a>
GPS_MakeV( gI, ph2,Sig2,  pB,  ph2,  Vb22); // <b|{1}|1>
// falSe second
GPS_MakeU( gI, ph2,Sig2,  ph1, pA,   U12a); // <2|[1]|a>
GPS_MakeV( gI, ph2,Sig2,  pB,  ph1,  Vb21); // <b|{1}|2>
// reverse order
GPS_MakeU( gI, ph2,Sig2,  pA,  ph1,  Ua21); // <a|[1]|2>
GPS_MakeV( gI, ph2,Sig2,  ph1, pB,   V12b); // <2|{1}|b>
// for the case when there was ph2 first xk-xk term
GPS_MakeU( gI, ph2,Sig2,  ph2, ph1,  U221); // <1|[1]|2>
GPS_MakeV( gI, ph2,Sig2,  ph1, ph2,  V122); // <2|{1}|1>
// for the case when there was ph2 first xk-xk term
GPS_MakeU( gI, ph2,Sig2,  ph1, ph1,  U121); // <2|[1]|2>
GPS_MakeV( gI, ph2,Sig2,  ph1, ph1,  V121); // <2|{1}|2>
//[[[[[[[[[[[[[[[[[[[[[
//(*m_Out)<< "U22a(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<U22a.m_A[j][k];(*m_Out)<<endl;
//(*m_Out)<< "Vb22(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<Vb22.m_A[j][k];(*m_Out)<<endl;
//
//(*m_Out)<< "U12a(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<U12a.m_A[j][k];(*m_Out)<<endl;
//(*m_Out)<< "Vb21(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<Vb21.m_A[j][k];(*m_Out)<<endl;
//
//(*m_Out)<< "Ua21(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<Ua21.m_A[j][k];(*m_Out)<<endl;
//(*m_Out)<< "V12b(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<V12b.m_A[j][k];(*m_Out)<<endl;
//
//(*m_Out)<< "U221(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<U221.m_A[j][k];(*m_Out)<<endl;
//(*m_Out)<< "V122(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<V122.m_A[j][k];(*m_Out)<<endl;
//
//(*m_Out)<< "U121(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<U121.m_A[j][k];(*m_Out)<<endl;
//(*m_Out)<< "V121(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<V121.m_A[j][k];(*m_Out)<<endl;
//]]]]]]]]]]]]]]]]]]]]]
//Sig = 3-2*Hel1
//CALL GPS_MakeUW(Cnor,ph1,Sig, pC,  mC,   pA,  mA,     UCAW1) ! v-a inside
//CALL GPS_MakeVW(Cnor,ph1,Sig, pB,  mB,   pD,  mD,     VBDW1) ! v-a inside
//CALL GPS_MakeUW(Cnor,ph1,Sig, pC,  mC,   ph2, mph,    UC2W1) ! v-a inside
//CALL GPS_MakeVW(Cnor,ph1,Sig, ph2,mph,   pD,  mD,     V2DW1) ! v-a inside
//Sig = 3-2*Hel2
//CALL GPS_MakeUW(Cnor,ph2,Sig, pC,  mC,   pA,   mA,    UCAW2) ! v-a inside
//CALL GPS_MakeVW(Cnor,ph2,Sig, pB,  mB,   pD,   mD,    VBDW2) ! v-a inside
//CALL GPS_MakeUW(Cnor,ph2,Sig, pC,  mC,   ph1, mph,    UC1W2) ! v-a inside
//CALL GPS_MakeVW(Cnor,ph2,Sig, ph1,mph,   pD,   mD,    V1DW2) ! v-a inside
KKcmplx2 UCAW1,VBDW1,UC2W1,V2DW1,UCAW2,VBDW2,UC1W2,V1DW2;
dcmplx Cone=1.0;
GPS_MakeUW(Cone, ph1, Sig1,  pC,  pA,   UCAW1);
GPS_MakeVW(Cone, ph1, Sig1,  pB,  pD,   VBDW1);
GPS_MakeUW(Cone, ph1, Sig1,  pC,  ph2,  UC2W1);  // v-a inside
GPS_MakeVW(Cone, ph1, Sig1,  ph2, pD,   V2DW1);  // v-a inside
//
GPS_MakeUW(Cone, ph2, Sig2,  pC,  pA,   UCAW2);  // v-a inside
GPS_MakeVW(Cone, ph2, Sig2,  pB,  pD,   VBDW2);  // v-a inside
GPS_MakeUW(Cone, ph2, Sig2,  pC,  ph1,  UC1W2);  // v-a inside
GPS_MakeVW(Cone, ph2, Sig2,  ph1, pD,   V1DW2);  // v-a inside
//[[[[[[[[[[[[[[[[[[[[[
//(*m_Out)<< "UCAW1(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<UCAW1.m_A[j][k];(*m_Out)<<endl;
//(*m_Out)<< "VBDW1(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<VBDW1.m_A[j][k];(*m_Out)<<endl;
//(*m_Out)<< "UC2W1(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<UC2W1.m_A[j][k];(*m_Out)<<endl;
//(*m_Out)<< "V2DW1(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<V2DW1.m_A[j][k];(*m_Out)<<endl;
//(*m_Out)<< "UCAW2(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<UCAW2.m_A[j][k];(*m_Out)<<endl;
//(*m_Out)<< "VBDW2(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<VBDW2.m_A[j][k];(*m_Out)<<endl;
//(*m_Out)<< "UC1W2(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<UC1W2.m_A[j][k];(*m_Out)<<endl;
//(*m_Out)<< "V1DW2(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<V1DW2.m_A[j][k];(*m_Out)<<endl;
//]]]]]]]]]]]]]]]]]]]]]]]]
//CALL GPS_MakeUX(Cnor,ph1,mph, pC,mC,   pA,mA,    UCAWX1) ! v-a inside
//CALL GPS_MakeVX(Cnor,ph1,mph, pB,mB,   pD,mD,    VBDWX1) ! v-a inside
//CALL GPS_MakeUX(Cnor,ph2,mph, pC,mC,   pA,mA,    UCAWX2) ! v-a inside
//CALL GPS_MakeVX(Cnor,ph2,mph, pB,mB,   pD,mD,    VBDWX2) ! v-a inside
//
//CALL GPS_MakeUX(Cnor,pB,mB, pC,mC,   pA,mA,    UCAWXB) ! v-a inside
//CALL GPS_MakeVX(Cnor,pA,mA, pB,mB,   pD,mD,    VBDWXA) ! v-a inside
//CALL GPS_MakeUX(Cnor,pD,mD, pC,mC,   pA,mA,    UCAWXD) ! v-a inside
//CALL GPS_MakeVX(Cnor,pC,mC, pB,mB,   pD,mD,    VBDWXC) ! v-a inside
KKcmplx2 UCAWX1,VBDWX1,UCAWX2,VBDWX2,UCAWXB,VBDWXA,UCAWXD,VBDWXC;
GPS_MakeUX(Cone, ph1, pC,   pA,  UCAWX1);  // v-a inside
GPS_MakeVX(Cone, ph1, pB,   pD,  VBDWX1);  // v-a inside
GPS_MakeUX(Cone, ph2, pC,   pA,  UCAWX2);  // v-a inside
GPS_MakeVX(Cone, ph2, pB,   pD,  VBDWX2);  // v-a inside
//
GPS_MakeUX(Cone, pB,  pC,   pA,  UCAWXB);   // v-a inside
GPS_MakeVX(Cone, pA,  pB,   pD,  VBDWXA);   // v-a inside
GPS_MakeUX(Cone, pD,  pC,   pA,  UCAWXD);   // v-a inside
GPS_MakeVX(Cone, pC,  pB,   pD,  VBDWXC);   // v-a inside
//[[[[[[[[[[[[[[[[[[[[[
//(*m_Out)<< "UCAWX1(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<UCAWX1.m_A[j][k];(*m_Out)<<endl;
//(*m_Out)<< "VBDWX1(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<VBDWX1.m_A[j][k];(*m_Out)<<endl;
//(*m_Out)<< "UCAWX2(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<UCAWX2.m_A[j][k];(*m_Out)<<endl;
//(*m_Out)<< "VBDWX2(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<VBDWX2.m_A[j][k];(*m_Out)<<endl;
//(*m_Out)<< "UCAWXB(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<UCAWXB.m_A[j][k];(*m_Out)<<endl;
//(*m_Out)<< "VBDWXA(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<VBDWXA.m_A[j][k];(*m_Out)<<endl;
//(*m_Out)<< "UCAWXD(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<UCAWXD.m_A[j][k];(*m_Out)<<endl;
//(*m_Out)<< "VBDWXC(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<VBDWXC.m_A[j][k];(*m_Out)<<endl;
//]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
//CALL GPS_MakeUX(Cnor,ph1,mph, pC, mC,   ph2,mph,    UC2WX1) ! v-a inside
//CALL GPS_MakeUX(Cnor,pA, mA,  pC, mC,   ph2,mph,    UC2WXA) ! v-a inside
//CALL GPS_MakeVX(Cnor,ph1,mph, ph2,mph,  pD, mD,     V2DWX1) ! v-a inside
//CALL GPS_MakeUX(Cnor,ph2,mph, pC, mC,   ph1,mph,    UC1WX2) ! v-a inside
//CALL GPS_MakeVX(Cnor,ph2,mph, ph1,mph,  pD, mD,     V1DWX2) ! v-a inside
KKcmplx2 UC2WX1,UC2WXA,V2DWX1,UC1WX2,V1DWX2;
GPS_MakeUX(Cone,ph1, pC,   ph2,  UC2WX1);   // v-a inside
GPS_MakeUX(Cone,pA,  pC,   ph2,  UC2WXA);   // v-a inside
GPS_MakeVX(Cone,ph1, ph2,  pD,   V2DWX1);   // v-a inside
GPS_MakeUX(Cone,ph2, pC,   ph1,  UC1WX2);   // v-a inside
GPS_MakeVX(Cone,ph2, ph1,  pD,   V1DWX2);   // v-a inside
//[[[[[[[[[[[[[[[[[[[[[
//(*m_Out)<< "UC2WX1(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<UC2WX1.m_A[j][k];(*m_Out)<<endl;
//(*m_Out)<< "UC2WXA(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<UC2WXA.m_A[j][k];(*m_Out)<<endl;
//(*m_Out)<< "V2DWX1(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<V2DWX1.m_A[j][k];(*m_Out)<<endl;
//(*m_Out)<< "V1DWX2(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<V1DWX2.m_A[j][k];(*m_Out)<<endl;
//]]]]]]]]]]]]]]]]]]]]]
//CALL GPS_MakeUX(Cnor,pB,mB, pC, mC,   ph2,mph,    UC2WXB) ! v-a inside
//CALL GPS_MakeVX(Cnor,pA,mA, ph2,mph,  pD, mD,     V2DWXA) ! v-a inside
//CALL GPS_MakeVX(Cnor,pB,mB, ph2,mph,  pD, mD,     V2DWXB) ! v-a inside
//CALL GPS_MakeUX(Cnor,pB,mB, pC, mC,   ph1,mph,    UC1WXB) ! v-a inside
//CALL GPS_MakeUX(Cnor,pA,mA, pC, mC,   ph1,mph,    UC1WXA) ! v-a inside
//CALL GPS_MakeVX(Cnor,pA,mA, ph1,mph,  pD, mD,     V1DWXA) ! v-a inside
KKcmplx2 UC2WXB,V2DWXA,V2DWXB,UC1WXB,UC1WXA,V1DWXA;
GPS_MakeUX(Cone, pB, pC,   ph2,   UC2WXB);   // v-a inside
GPS_MakeVX(Cone, pA, ph2,  pD,    V2DWXA);   // v-a inside
GPS_MakeVX(Cone, pB, ph2,  pD,    V2DWXB);   // v-a inside
GPS_MakeUX(Cone, pB, pC,   ph1,   UC1WXB);   // v-a inside
GPS_MakeUX(Cone, pA, pC,   ph1,   UC1WXA);   // v-a inside
GPS_MakeVX(Cone, pA, ph1,  pD,    V1DWXA);   // v-a inside
//[[[[[[[[[[[[[[[[[[[[[
//(*m_Out)<< "UC2WXB(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<UC2WXB.m_A[j][k];(*m_Out)<<endl;
//(*m_Out)<< "V2DWXA(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<V2DWXA.m_A[j][k];(*m_Out)<<endl;
//(*m_Out)<< "V2DWXB(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<V2DWXB.m_A[j][k];(*m_Out)<<endl;
//(*m_Out)<< "UC1WXB(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<UC1WXB.m_A[j][k];(*m_Out)<<endl;
//(*m_Out)<< "UC1WXA(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<UC1WXA.m_A[j][k];(*m_Out)<<endl;
//(*m_Out)<< "V1DWXA(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<V1DWXA.m_A[j][k];(*m_Out)<<endl;
//]]]]]]]]]]]]]]]]]]]]]
//CALL GPS_MakeUX(Cnor,pD,mD, pC, mC,   ph2,mph,    UC2WXD) ! v-a inside
//CALL GPS_MakeVX(Cnor,pC,mC, ph2,mph,  pD, mD,     V2DWXC) ! v-a inside
//CALL GPS_MakeUX(Cnor,pD,mD, pC, mC,   ph1,mph,    UC1WXD) ! v-a inside
//CALL GPS_MakeVX(Cnor,pC,mC, ph1,mph,  pD, mD,     V1DWXC) ! v-a inside
//CALL GPS_MakeVX(Cnor,pB,mB, ph1,mph,  pD, mD,     V1DWXB) ! v-a inside
//
KKcmplx2 UC2WXD,V2DWXC,UC1WXD,V1DWXC,V1DWXB;
GPS_MakeUX(Cone, pD, pC,   ph2,    UC2WXD);   // v-a inside
GPS_MakeVX(Cone, pC, ph2,  pD,     V2DWXC);   // v-a inside
GPS_MakeUX(Cone, pD, pC,   ph1,    UC1WXD);   // v-a inside
GPS_MakeVX(Cone, pC, ph1,  pD,     V1DWXC);   // v-a inside
GPS_MakeVX(Cone, pB, ph1,  pD,     V1DWXB);   // v-a inside
//[[[[[[[[[[[[[[[[[[[[[
//(*m_Out)<< "UC2WXD(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<UC2WXD.m_A[j][k];(*m_Out)<<endl;
//(*m_Out)<< "V2DWXC(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<V2DWXC.m_A[j][k];(*m_Out)<<endl;
//(*m_Out)<< "UC1WXD(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<UC1WXD.m_A[j][k];(*m_Out)<<endl;
//(*m_Out)<< "V1DWXC(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<V1DWXC.m_A[j][k];(*m_Out)<<endl;
//(*m_Out)<< "V1DWXB(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<V1DWXB.m_A[j][k];(*m_Out)<<endl;
//]]]]]]]]]]]]]]]]]]]]]]]]
//------------------------------------------------------------------------------------------
int Y_IR=1;        //    ! YES, IR included
Y_IR=0;            //    ! No,  IR not included
int Y_IR1=1;       //    ! defunct, must be 1 use I9X-T YES, IR times beta 0 included
int N_IR=1-Y_IR1;  //
N_IR=1;            //    ! #########################################################
//-----------------------------------------------------------------------------------------

// Assembling amplitude
//[[[[[[[[[[[[[[[[[[[
//Amp4Zer(m_AmpExpo2);
//]]]]]]]]]]]]]]]]]]]
/////////////////////////////////////////////////////////////////////////////////////
//                    |                  1               2                         //
//                    |X                 |               |                         //
//      _       -b    |      b+m-1-2     |      a+m-2    |      a                  //
//      v  ------<----O--------<---------U--------<------S------<----- u           //
/////////////////////////////////////////////////////////////////////////////////////
//    Su1=Su1+Born1BCD(j,j2,j3,j4) *U11a(j,j1)*(prA12-prA1)*sA(2,Hel2)*CPF12/CPF0*IV2       !<b|X|1[1]a(2)|a>
//    Su1=Su1+Born1BCD(j,j2,j3,j4) *U11a(j,j1)*(      prA1)*sA(2,Hel2)*CPF12/CPF0*Y_IR1*I9X !<b|X|1[1]a(2)|a>
//    Su1=Su1+Born2BCD(j,j2,j3,j4) *U21a(j,j1)* prA12      *sA(2,Hel2)*CPF12/CPF0*IV2       !<b|X|2[1]a(2)|a>
AmpAddI(m_AmpExpo2, CNorm*(prA12-prA1)*sA[1][Hel2]*CPF12/CPF0*double(IV2),       Born1BCD, U11a); // !<b|X|1[1]a(2)|a>
AmpAddI(m_AmpExpo2, CNorm*(      prA1)*sA[1][Hel2]*CPF12/CPF0*double(Y_IR1*I9X), Born1BCD, U11a); // !<b|X|1[1]a(2)|a>
AmpAddI(m_AmpExpo2, CNorm*(prA12     )*sA[1][Hel2]*CPF12/CPF0*double(IV2),       Born2BCD, U21a); // !<b|X|2[1]a(2)|a>

/////////////////////////////////////////////////////////////////////////////////////
//                    |                  2               1                         //
//                    |X                 |               |                         //
//      _       -b    |      b+m-1-2     |      a+m-1    |      a                  //
//      v  ------<----O--------<---------U-------<-------S------<----- u           //
/////////////////////////////////////////////////////////////////////////////////////
//      Su1=Su1+Born2BCD(j,j2,j3,j4) *U22a(j,j1)*(prA12-prA2)*sA(1,Hel1) *CPF12/CPF0*IV2       !<b|X|2[2]a(1)|a>
//      Su1=Su1+Born2BCD(j,j2,j3,j4) *U22a(j,j1)*(      prA2)*sA(1,Hel1) *CPF12/CPF0*Y_IR1*I9Z !<b|X|2[2]a(1)|a>
//      Su1=Su1+Born1BCD(j,j2,j3,j4) *U12a(j,j1)* prA12      *sA(1,Hel1) *CPF12/CPF0*IV2       !<b|X|1[2]a(1)|a>
AmpAddI(m_AmpExpo2, CNorm*(prA12-prA2)*sA[0][Hel1] *CPF12/CPF0*double(IV2),       Born2BCD, U22a);
AmpAddI(m_AmpExpo2, CNorm*(      prA2)*sA[0][Hel1] *CPF12/CPF0*double(Y_IR1*I9Z), Born2BCD, U22a);
AmpAddI(m_AmpExpo2, CNorm* prA12      *sA[0][Hel1] *CPF12/CPF0*double(IV2),       Born1BCD, U12a);

/////////////////////////////////////////////////////////////////////////////////////
//                    2                  1               |                         //
//                    |                  |               |X                        //
//      _       -b    |     -b+m+2       |    -b+m+1+2   |      a                  //
//      v  ------<----S--------<---------V--------<------O------<----- u           //
/////////////////////////////////////////////////////////////////////////////////////
//      Su1=Su1+( sB(2,Hel2)*(-prB1+prB12))*Vb11(j2,j)*BornA1CD(j1,j,j3,j4)*CPF/CPF0*IV1       !<b|(2)b[1]1|X|a>
//      Su1=Su1+( sB(2,Hel2)*( prB1      ))*Vb11(j2,j)*BornA1CD(j1,j,j3,j4)*CPF/CPF0*Y_IR1*I9Y !<b|(2)b[1]1|X|a>
//      Su1=Su1+( sB(2,Hel2)*(      prB12))*Vb12(j2,j)*BornA2CD(j1,j,j3,j4)*CPF/CPF0*IV1       !<b|(2)b[1]2|X|a>
AmpAddI(m_AmpExpo2, CNorm*(sB[1][Hel2]*(-prB1+prB12))*CPF/CPF0*double(IV1),       Vb11, BornA1CD);
AmpAddI(m_AmpExpo2, CNorm*(sB[1][Hel2]*( prB1      ))*CPF/CPF0*double(Y_IR1*I9Y), Vb11, BornA1CD);
AmpAddI(m_AmpExpo2, CNorm*(sB[1][Hel2]*(      prB12))*CPF/CPF0*double(IV1),       Vb12, BornA2CD);
//(*m_Out)<<"sB[1][Hel2]="<<sB[1][Hel2]<<endl;
//(*m_Out)<<"      prB12="<<prB12<<endl;
//(*m_Out)<<" CPF, CPF0 ="<<CPF<<"  "<<CPF0<<endl;
//(*m_Out)<< "Vb12(*,*)=";for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<Vb12.m_A[j][k];(*m_Out)<<endl;
//for(int j1=0; j1<=1; j1++) for(int j2=0; j2<=1; j2++){ (*m_Out)<< "BornA2CD("<<j1<<","<<j2<<",*,*)=  ";
//     for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<BornA2CD.m_A[j1][j2][j][k];(*m_Out)<<endl;}
/////////////////////////////////////////////////////////////////////////////////////
//                    1                  2               |                         //
//                    |                  |               |X                        //
//      _       -b    |     -b+m+1       |    -b+m+1+2   |      a                  //
//      v  ------<----S--------<---------V--------<------O------<----- u           //
/////////////////////////////////////////////////////////////////////////////////////
//      Su1=Su1 +(sB(1,Hel1)*(-prB2+prB12))*Vb22(j2,j)*BornA2CD(j1,j,j3,j4)*CPF/CPF0*IV1       !<b|(1)b[2]2|X|a>
//      Su1=Su1 +(sB(1,Hel1)*( prB2      ))*Vb22(j2,j)*BornA2CD(j1,j,j3,j4)*CPF/CPF0*Y_IR1*I9T !<b|(1)b[2]2|X|a>
//      Su1=Su1 +(sB(1,Hel1)*(      prB12))*Vb21(j2,j)*BornA1CD(j1,j,j3,j4)*CPF/CPF0*IV1       !<b|(1)b[2]1|X|a>
AmpAddI(m_AmpExpo2, CNorm*(sB[0][Hel1]*(-prB2+prB12))*CPF/CPF0*double(IV1),       Vb22,  BornA2CD);
AmpAddI(m_AmpExpo2, CNorm*(sB[0][Hel1]*( prB2      ))*CPF/CPF0*double(Y_IR1*I9T), Vb22,  BornA2CD);
AmpAddI(m_AmpExpo2, CNorm*(sB[0][Hel1]*(      prB12))*CPF/CPF0*double(IV1),       Vb21,  BornA1CD);

/////////////////////////////////////////////////////////////////////////////////////
//                    2                  |               1                         //
//                    |                  |X              |                         //
//      _       -b    |     -b+m+2       |     a+m-1     |      a                  //
//      v  ------<----S--------<---------O--------<------U------<----- u           //
/////////////////////////////////////////////////////////////////////////////////////
//              Su1=Su1 +sB(2,Hel2)    *Born1BCD(j,j2,j3,j4) *prA1*U11a(j,j1) *CPF1/CPF0*Y_IR1*I9X !<b|(2)2|X|1[1]|a>
AmpAddI(m_AmpExpo2, CNorm*sB[1][Hel2]*prA1*CPF1/CPF0*double(Y_IR1*I9X),   Born1BCD, U11a);

/////////////////////////////////////////////////////////////////////////////////////
//                    1                  |               2                         //
//                    |                  |X              |                         //
//      _       -b    |     -b+m+1       |     a+m-2     |      a                  //
//      v  ------<----S--------<---------O--------<------U------<----- u           //
/////////////////////////////////////////////////////////////////////////////////////
//            Su1=Su1 +sB(1,Hel1)    *Born2BCD(j,j2,j3,j4) *prA2*U22a(j,j1) *CPF2/CPF0*Y_IR1*I9Z !<b|(1)1|X|2[2]|a>
AmpAddI(m_AmpExpo2, CNorm*sB[0][Hel1]*prA2*CPF2/CPF0* double(Y_IR1*I9Z) ,   Born2BCD, U22a);

/////////////////////////////////////////////////////////////////////////////////////
//                    1                  |               2                         //
//                    |                  |X              |                         //
//      _       -b    |     -b+m+1       |     a+m-2     |      a                  //
//      v  ------<----V--------<---------O--------<------S------<----- u           //
/////////////////////////////////////////////////////////////////////////////////////
//             Su1=Su1 +Vb11(j2,j)*prB1 *BornA1CD(j1,j,j3,j4) *sA(2,Hel2) *CPF2/CPF0   *Y_IR1*I9Y !<b|[1]1|X|a[2]|a>
AmpAddI(m_AmpExpo2, CNorm*prB1*sA[1][Hel2]*CPF2/CPF0*double(Y_IR1*I9Y),  Vb11,  BornA1CD);

/////////////////////////////////////////////////////////////////////////////////////
//                    2                  |               1                         //
//                    |                  |X              |                         //
//      _       -b    |     -b+m+2       |     a+m-1     |      a                  //
//      v  ------<----V--------<---------O--------<------S------<----- u           //
/////////////////////////////////////////////////////////////////////////////////////
//             Su1=Su1 +Vb22(j2,j)*prB2 *BornA2CD(j1,j,j3,j4) *sA(1,Hel1)*CPF1/CPF0      *Y_IR1*I9T !<b|[2]2|X|a(1)|a>
AmpAddI(m_AmpExpo2, CNorm*prB2*sA[0][Hel1]*CPF1/CPF0 *double(Y_IR1*I9T), Vb22, BornA2CD);

/////////////////////////////////////////////////////////////////////////////////////
//                                       E--2                                      //
//                                       |               1                         //
//                                       |X              |                         //
//      _       -b                       |     a+m-1     |      a                  //
//      v  ------<-------------<---------O--------<------U------<----- u           //
/////////////////////////////////////////////////////////////////////////////////////
//            Su1=Su1 +EpsDot21(Hel2)  *Born1BCD(j,j2,j3,j4) *prA1*U11a(j,j1) *CPF1*CPF12/CPF0*I9X    !<b|(2)2|X|1[1]|a>
//     $                       -0.5D0*DCMPLX(ChaIni*m_e_QED)*prA1*U11a(j,j1)*CPF1*CPF12*WVpi0             *I71
//     $                       *( UC1W2(j3,j )*(2*VBDWX2(j2,j4) )
//     $                         -VBDW2(j2,j4)*(2*UC1WX2(j3,j ) )                        )
AmpAddI(m_AmpExpo2, CNorm*EpsDot21[Hel2]*prA1*CPF1*CPF12/CPF0*double(I9X),     Born1BCD, U11a);
Amp2to4(m_AmpTemp1, -2*0.5*dcmplx(  ChaIni*m_e_QED)*prA1*CPF1*CPF12*WVPi0, VBDWX2, UC1W2);
Amp2to4(m_AmpTemp2, +2*0.5*dcmplx(  ChaIni*m_e_QED)*prA1*CPF1*CPF12*WVPi0, VBDW2, UC1WX2);
AmpAddI(m_AmpExpo2, CNorm, m_AmpTemp1,U11a);
AmpAddI(m_AmpExpo2, CNorm, m_AmpTemp2,U11a);

/////////////////////////////////////////////////////////////////////////////////////
//                                       E--1                                      //
//                                       |               2                         //
//                                       |X              |                         //
//      _       -b                       |     a+m-2     |      a                  //
//      v  ------<-------------<---------O--------<------U------<----- u           //
/////////////////////////////////////////////////////////////////////////////////////
//                     Su1=Su1 +EpsDot12(Hel1)  *Born2BCD(j,j2,j3,j4) *prA2*U22a(j,j1) *CPF2*CPF12/CPF0*I9Z !<b|(1)1|X|2[2]|a>
//     $                       -0.5D0*DCMPLX(ChaIni*m_e_QED)*prA2*U22a(j,j1)*CPF2*CPF12*WVpi0*I72
//     $                       *( UC2W1(j3,j )*(2*VBDWX1(j2,j4) )
//     $                         -VBDW1(j2,j4)*(2*UC2WX1(j3,j )              ) )
AmpAddI(m_AmpExpo2, CNorm*EpsDot12[Hel1]*prA2*CPF2*CPF12/CPF0*double(I9Z),   Born2BCD, U22a);
Amp2to4(m_AmpTemp1, -2*0.5*dcmplx(ChaIni*m_e_QED)*prA2*CPF2*CPF12*WVPi0*double(I72), VBDWX1, UC2W1);
Amp2to4(m_AmpTemp2, +2*0.5*dcmplx(ChaIni*m_e_QED)*prA2*CPF2*CPF12*WVPi0*double(I72), VBDW1, UC2WX1);
AmpAddI(m_AmpExpo2, CNorm, m_AmpTemp1,U22a);
AmpAddI(m_AmpExpo2, CNorm, m_AmpTemp2,U22a);

/////////////////////////////////////////////////////////////////////////////////////
//                                       E--2                                      //
//                    1                  |                                         //
//                    |                  |X                                        //
//      _       -b    |     -b+m+1       |                      a                  //
//      v  ------<----V--------<---------O--------<-------------<----- u           //
/////////////////////////////////////////////////////////////////////////////////////
//                     Su1=Su1 +Vb11(j2,j)*prB1 *BornA1CD(j1,j,j3,j4) *EpsDot2(Hel2) *CPF *CPF2/CPF0*I9Y  !<b|[1]1|X|a[2]|a>
//     $                       -0.5D0*DCMPLX(ChaIni*m_e_QED)*Vb11(j2,j)*prB1*CPF*CPF2*WVpi0*I71
//     $                       *( UCAW2(j3,j1)*(2*V1DWX2(j ,j4)              )
//     $                         -V1DW2(j ,j4)*(2*UCAWX2(j3,j1)) )
AmpAddI(m_AmpExpo2, CNorm*prB1*EpsDot2[Hel2] *CPF *CPF2/CPF0*double(I9Y) ,  Vb11,  BornA1CD);
Amp2to4(m_AmpTemp1, -2*0.5*dcmplx(ChaIni*m_e_QED)*prB1*CPF*CPF2*WVPi0*double(I71), V1DWX2, UCAW2);
Amp2to4(m_AmpTemp2, +2*0.5*dcmplx(ChaIni*m_e_QED)*prB1*CPF*CPF2*WVPi0*double(I71), V1DW2, UCAWX2);
AmpAddI(m_AmpExpo2, CNorm, Vb11, m_AmpTemp1);
AmpAddI(m_AmpExpo2, CNorm, Vb11, m_AmpTemp2);

/////////////////////////////////////////////////////////////////////////////////////
//                                       E--1                                      //
//                    2                  |                                         //
//                    |                  |X                                        //
//      _       -b    |     -b+m+2       |                      a                  //
//      v  ------<----V--------<---------O--------<-------------<----- u           //
/////////////////////////////////////////////////////////////////////////////////////
//                     Su1=Su1 +Vb22(j2,j)*prB2 *BornA2CD(j1,j,j3,j4) *EpsDot1(Hel1)*CPF *CPF1/CPF0*I9T  !<b|[2]2|X|a(1)|a>
//     $                       -0.5D0*DCMPLX(ChaIni*m_e_QED)*Vb22(j2,j)*prB2*CPF*CPF1*WVpi0*I72
//     $                       *( UCAW1(j3,j1)*(2*V2DWX1(j ,j4)               )
//     $                         -V2DW1(j ,j4)*(2*UCAWX1(j3,j1) ) )
AmpAddI(m_AmpExpo2, CNorm*prB2*EpsDot1[Hel1]*CPF *CPF1/CPF0*double(I9T),   Vb22,  BornA2CD);
Amp2to4(m_AmpTemp1, -2*0.5*dcmplx(ChaIni*m_e_QED) *prB2*CPF*CPF1*WVPi0*double(I72), V2DWX1, UCAW1);
Amp2to4(m_AmpTemp2, +2*0.5*dcmplx(ChaIni*m_e_QED) *prB2*CPF*CPF1*WVPi0*double(I72), V2DW1, UCAWX1);
AmpAddI(m_AmpExpo2, CNorm, Vb22, m_AmpTemp1);
AmpAddI(m_AmpExpo2, CNorm, Vb22, m_AmpTemp2);

//-----------------------------------------------------------------------------------

/////////////////////////////////////////////////////////////////////////////////////
//                    2                  |               1                         //
//                    |                  |X              |                         //
//      _       -b    |     -b+m+2       |     a+m-1     |      a                  //
//      v  ------<----*--------<---------O--------<------O------<----- u           //
/////////////////////////////////////////////////////////////////////////////////////
//           Su2=Su2  +Vb22(j2,l)*prB2  *Born12CD(j,l,j3,j4 )  *U11a(j,j1)*prA1*CPF1/CPF0*I8 ! <b|[2]2|X|1[1]|a>
AmpAddI(m_AmpExpo2, CNorm*prB2*prA1*CPF1/CPF0*double(I8),   Vb22,  Born12CD, U11a);

/////////////////////////////////////////////////////////////////////////////////////
//                    1                  |               2                         //
//                    |                  |X              |                         //
//      _       -b    |     -b+m+1       |     a+m-2     |      a                  //
//      v  ------<----*--------<---------O--------<------O------<----- u           //
/////////////////////////////////////////////////////////////////////////////////////
//           Su2=Su2  +Vb11(j2,l)*prB1  *Born21CD(j,l,j3,j4 )  *U22a(j,j1)*prA2*CPF2/CPF0*I8 ! <b|[1]1|X|2[2]|a>
AmpAddI(m_AmpExpo2, CNorm*prB1*prA2*CPF2/CPF0*double(I8) ,   Vb11,  Born21CD, U22a);

/////////////////////////////////////////////////////////////////////////////////////
//                    |                  2               1                         //
//                    |X                 |               |                         //
//      _       -b    |      b+m-1-2     |      a+m-1    |      a                  //
//      v  ------<----O--------<---------O--------<------*------<----- u           //
/////////////////////////////////////////////////////////////////////////////////////
//           Su2=Su2 +Born1BCD(j,j2,j3,j4) *U121(j,l)*prA12 *U11a(l,j1)*prA1*CPF12/CPF0*IV2  ! <b|X|1[2]1(1)|a>
//           Su2=Su2 +Born2BCD(j,j2,j3,j4) *U221(j,l)*prA12 *U11a(l,j1)*prA1*CPF12/CPF0*IA  ! <b|X|2[2]1(1)|a>
//           Su2=Su2 -BornABCD(j,j2,j3,j4) *Ua21(j,l)*prA12 *U11a(l,j1)*prA1*CPF12/CPF0*IV2  ! <b|X|a[2]1(1)|a>
AmpAddI(m_AmpExpo2, CNorm*prA12*prA1*CPF12/CPF0*double(IV2),   Born1BCD, U121, U11a);
AmpAddI(m_AmpExpo2, CNorm*prA12*prA1*CPF12/CPF0*double(IA),    Born2BCD, U221, U11a);
AmpAddI(m_AmpExpo2,-CNorm*prA12*prA1*CPF12/CPF0*double(IV2),   BornABCD, Ua21, U11a);

/////////////////////////////////////////////////////////////////////////////////////
//                    |                  1               2                         //
//                    |X                 |               |                         //
//      _       -b    |      b+m-1-2     |      a+m-2    |      a                  //
//      v  ------<----O--------<---------O--------<------*------<----- u           //
/////////////////////////////////////////////////////////////////////////////////////
//           Su2=Su2 +Born2BCD(j,j2,j3,j4) *U212(j,l)*prA12  *U22a(l,j1)*prA2*CPF12/CPF0*IV2 ! <b|X|2[1]2(2)|a>
//           Su2=Su2 +Born1BCD(j,j2,j3,j4) *U112(j,l)*prA12  *U22a(l,j1)*prA2*CPF12/CPF0*IA  ! <b|X|1[1]2(2)|a>
//           Su2=Su2 -BornABCD(j,j2,j3,j4) *Ua12(j,l)*prA12  *U22a(l,j1)*prA2*CPF12/CPF0*IV2 ! <b|X|a[1]2(2)|a>
AmpAddI(m_AmpExpo2, CNorm*prA12*prA2*CPF12/CPF0*double(IV2),   Born2BCD, U212, U22a);
AmpAddI(m_AmpExpo2, CNorm*prA12*prA2*CPF12/CPF0*double(IA),    Born1BCD, U112, U22a);
AmpAddI(m_AmpExpo2,-CNorm*prA12*prA2*CPF12/CPF0*double(IV2),   BornABCD, Ua12, U22a);

/////////////////////////////////////////////////////////////////////////////////////
//                    1                  2               |                         //
//                    |                  |               |X                        //
//      _       -b    |     -b+m+1       |    -b+m+1+2   |      a                  //
//      v  ------<----*--------<---------O--------<------O------<----- u           //
/////////////////////////////////////////////////////////////////////////////////////
//            Su2=Su2 +Vb11(j2,l)*prB1 *V121(l,j)*prB12 *BornA1CD(j1,j,j3,j4)*CPF/CPF0*IV1 ! <b|(1)1[2]1|X|a>
//            Su2=Su2 +Vb11(j2,l)*prB1 *V122(l,j)*prB12 *BornA2CD(j1,j,j3,j4)*CPF/CPF0*IA  ! <b|(1)1[2]2|X|a>
//            Su2=Su2 -Vb11(j2,l)*prB1 *V12b(l,j)*prB12 *BornABCD(j1,j,j3,j4)*CPF/CPF0*IV1 ! <b|(1)1[2]b|X|a>
AmpAddI(m_AmpExpo2, CNorm*prB1*prB12*CPF/CPF0*double(IV1),  Vb11, V121, BornA1CD);
AmpAddI(m_AmpExpo2, CNorm*prB1*prB12*CPF/CPF0*double(IA ),  Vb11, V122, BornA2CD);
AmpAddI(m_AmpExpo2,-CNorm*prB1*prB12*CPF/CPF0*double(IV1),  Vb11, V12b, BornABCD);

/////////////////////////////////////////////////////////////////////////////////////
//                    2                  1               |                         //
//                    |                  |               |X                        //
//      _       -b    |     -b+m+2       |    -b+m+1+2   |      a                  //
//      v  ------<----*--------<---------O--------<------O------<----- u           //
/////////////////////////////////////////////////////////////////////////////////////
//            Su2=Su2 +Vb22(j2,l)*prB2 *V212(l,j)*prB12  *BornA2CD(j1,j,j3,j4)*CPF/CPF0*IV1 ! <b|(2)2[1]2|X|a>
//            Su2=Su2 +Vb22(j2,l)*prB2 *V211(l,j)*prB12  *BornA1CD(j1,j,j3,j4)*CPF/CPF0*IA  ! <b|(2)2[1]1|X|a>
//            Su2=Su2 -Vb22(j2,l)*prB2 *V21b(l,j)*prB12  *BornABCD(j1,j,j3,j4)*CPF/CPF0*IV1 ! <b|(2)2[2]b|X|a>
AmpAddI(m_AmpExpo2, CNorm*prB2*prB12*CPF/CPF0*double(IV1),  Vb22, V212, BornA2CD);
AmpAddI(m_AmpExpo2, CNorm*prB2*prB12*CPF/CPF0*double(IA ),  Vb22, V211, BornA1CD);
AmpAddI(m_AmpExpo2,-CNorm*prB2*prB12*CPF/CPF0*double(IV1),  Vb22, V21b, BornABCD);

//===========================
//c      sProd = ( sA(1,Hel1)+sB(1,Hel1))*( sA(2,Hel2)+sB(2,Hel2))
//c      AmpWork(j1,j2,j3,j4) =AmpWork(j1,j2,j3,j4) ! non-infrared part: Fprop1- Fprop1B) is prop to ph1*ph2
//c     $                 +CNorm*BornABCD(j1,j2,j3,j4)*( sA(1,Hel1)*sA(2,Hel2)*(Fprop1- Fprop1B)      *IV2 !
//c     $                                               +sB(1,Hel1)*sB(2,Hel2)*(Fprop2- Fprop2B)      *IV1
//c     $                                                                                      ) !
dcmplx cFactX = CNorm*( sA[0][Hel1]*sA[1][Hel2]*(Fprop1- Fprop1B)*double(IV2)
                       +sB[0][Hel1]*sB[1][Hel2]*(Fprop2- Fprop2B)*double(IV1) );
AmpAdd( m_AmpExpo2,  cFactX ,BornABCD);
//c===========================
//c      AmpWork(j1,j2,j3,j4) =AmpWork(j1,j2,j3,j4) ! infrared-type part
//c     $                 +CNorm*BornABCD(j1,j2,j3,j4)*( sA(1,Hel1)*sA(2,Hel2)*Fprop1B                      *IVI !
//c     $                                               +sA(1,Hel1)*sB(2,Hel2)*(CPF1/CPF0-1d0)              *IVI
//c     $                                               +sB(1,Hel1)*sA(2,Hel2)*(CPF2/CPF0-1d0)              *IVI
//c     $                                               +sB(1,Hel1)*sB(2,Hel2)*Fprop2B                      *IVI
//c     $                                               +sA(1,Hel1)*CPF1*CPF12/CPF0*EpsDot21(Hel2)          *IVI! terms due to diag WWgamma
//c     $                                               +sA(2,Hel2)*CPF2*CPF12/CPF0*EpsDot12(Hel1)          *IVI
//c     $                                               +sB(1,Hel1)*CPF*CPF2/CPF0*EpsDot2(Hel2)             *IVI
//c     $                                               +sB(2,Hel2)*CPF*CPF1/CPF0*EpsDot1(Hel1)             *IVI
//c     $                                               +CPF*CPF1*CPF12/CPF0*EpsDot1(Hel1)*EpsDot21(Hel2)   *IVI
//c     $                                               +CPF*CPF2*CPF12/CPF0*EpsDot2(Hel2)*EpsDot12(Hel1)   *IVI
//c     $                                                                                      ) !
//c     $                 +CNorm*BornABCD(j1,j2,j3,j4)* sProd                   *Y_IR                        !
dcmplx cFactY=
        CNorm*( sA[0][Hel1]*sA[1][Hel2]*Fprop1B                     *double(IVI)
               +sA[0][Hel1]*sB[1][Hel2]*(CPF1/CPF0-1.0)             *double(IVI)
               +sB[0][Hel1]*sA[1][Hel2]*(CPF2/CPF0-1.0)             *double(IVI)
               +sB[0][Hel1]*sB[1][Hel2]*Fprop2B                     *double(IVI)
               +sA[0][Hel1]*CPF1*CPF12/CPF0*EpsDot21[Hel2]          *double(IVI)  // terms due to diag WWgamma
               +sA[1][Hel2]*CPF2*CPF12/CPF0*EpsDot12[Hel1]          *double(IVI)
               +sB[0][Hel1]*CPF*CPF2/CPF0*EpsDot2[Hel2]             *double(IVI)
               +sB[1][Hel2]*CPF*CPF1/CPF0*EpsDot1[Hel1]             *double(IVI)
               +CPF*CPF1*CPF12/CPF0*EpsDot1[Hel1]*EpsDot21[Hel2]    *double(IVI)
               +CPF*CPF2*CPF12/CPF0*EpsDot2[Hel2]*EpsDot12[Hel1]    *double(IVI) );
dcmplx sProd =( sA[0][Hel1]+sB[0][Hel1])*( sA[1][Hel2]+sB[1][Hel2]);
AmpAdd( m_AmpExpo2,  cFactY ,BornABCD);
AmpAdd( m_AmpExpo2,  sProd*double(Y_IR) ,BornABCD);

//===========================
////dcmplx gI = dcmplx(ChaIni*m_e_QED);
//! terms due to ph1 ph2 attached to W
//c      AmpWork(j1,j2,j3,j4) =AmpWork(j1,j2,j3,j4) + CNorm*( 0.5D0) *DCMPLX(-ChaIni*m_e_QED)*CPF*CPF1*CPF12*WVpi0*(
//c     $  ! absorbed to infraed-type      BornABCD(j1,j2,j3,j4)/CPF0/WVpi0*(2*eps1pC-2*eps1pA)*(2*eps2pB-2*eps2pD)* *IV ! 1
//c     $ +EpsDot1 (hel1)*VBDW2(j2,j4)*(               -2*UCAWX2(j3,j1)) *I9s1 ! 2
//c     $ +EpsDot1 (hel1)*VBDW2(j2,j4)*(-UCAWX1(j3,j1)                 ) *I9   ! 2
//c     $  +EpsDot1 (hel1)*2*VBDWX2(j2,j4) * UCAW2 (j3,j1)                *I9s1 ! 3
//c     $  - 2*UCAWX1(j3,j1)*EpsDot21(hel2)* VBDW1 (j2,j4)                *I9s2 ! 4
//c     $ -UCAWX1(j3,j1)* VBDW2 (j2,j4) *(-EpsDot12(hel1)-2*eps1p2*DCMPLX(-ChaIni*m_e_QED))  *I9   ! 5
//c     $ + eps1D2*(-2*UCAWX1(j3,j1))*        2* VBDWX2(j2,j4)    *DCMPLX(-ChaIni*m_e_QED)   *I9   ! 6
//c     $ + UCAW1(j3,j1)*(2*VBDWX1(j2,j4)              )*EpsDot21(hel2)  *I9s2 ! 7
//c     $ + UCAW1(j3,j1)*(                VBDWX2(j2,j4))*EpsDot21(hel2)  *I9   ! 7
//!!!!     $ + UCAW1(j3,j1)*(t1-2*p2p1 )*VBDW2(j2,j4)          *DCMPLX(-ChaIni*m_e_QED)  *I9 ! 8 (out)
//c     $ + UCAW1(j3,j1)*(1D0/CPF2-2*p2p1 )*VBDW2(j2,j4)       *DCMPLX(-ChaIni*m_e_QED)  *I9 ! 8    'Higgs' contrib added
//c     $ + UCAW1(j3,j1)*VBDWX2(j2,j4)*(-EpsDot2(hel2)+2*eps2p1*DCMPLX(-ChaIni*m_e_QED)) *I9    ! 9
//c     $                                                                            )
dcmplx cFactZ= CNorm*0.5 *(-gI)*CPF*CPF1*CPF12*WVPi0;
Amp2to4add( m_AmpExpo2, cFactZ*(-2.0)*EpsDot1[Hel1] *double(I9s1), VBDW2, UCAWX2);   // ! 1 (j2,j4)*(j3,j1)
Amp2to4add( m_AmpExpo2, cFactZ*(-1.0)*EpsDot1[Hel1] *double(I9),   VBDW2, UCAWX1);   // ! 2
Amp2to4add( m_AmpExpo2, cFactZ*( 2.0)*EpsDot1[Hel1] *double(I9s1), VBDWX2, UCAW2);   // ! 3
Amp2to4add( m_AmpExpo2, cFactZ*(-2.0)*EpsDot21[Hel2]*double(I9s2), VBDW1, UCAWX1);   // ! 4
Amp2to4add( m_AmpExpo2, cFactZ*( EpsDot12[Hel1]+2.0*eps1p2*(-gI))*double(I9),   VBDW2,  UCAWX1); //! 5
Amp2to4add( m_AmpExpo2, cFactZ*(-2.0)*eps1D2*2.0*(-gI)*double(I9),    VBDWX2, UCAWX1);           //! 6
Amp2to4add( m_AmpExpo2, cFactZ*( 2.0)*EpsDot21[Hel2]  *double(I9s2),  VBDWX1, UCAW1);            //! 7
Amp2to4add( m_AmpExpo2, cFactZ*( 1.0)*EpsDot21[Hel2]  *double(I9),    VBDWX2, UCAW1);            //! 7
Amp2to4add( m_AmpExpo2, cFactZ*((1.0/CPF2-2.0*Cp2p1 ))*(-gI)*double(I9),      VBDW2, UCAW1);     //! 8
Amp2to4add( m_AmpExpo2, cFactZ*(-EpsDot2[Hel2]+2.0*eps2p1*(-gI))*double(I9), VBDWX2, UCAW1);     //! 9
//c===========================
//! terms due to ph2 ph1 attached to W
//c      AmpWork(j1,j2,j3,j4) =AmpWork(j1,j2,j3,j4) + CNorm*( 0.5D0) *DCMPLX(-ChaIni*m_e_QED)*CPF*CPF2*CPF12*WVpi0*(
//c     $   ! absorbed to infraed-type     BornABCD(j1,j2,j3,j4)/CPF0/WVpi0*(2*eps2pC-2*eps2pA)*(2*eps1pB-2*eps1pD)*IV ! 1
//c     $ +EpsDot2 (hel2)*VBDW1(j2,j4)*(               -2*UCAWX1(j3,j1)) *I9s2 ! 2
//c     $ +EpsDot2 (hel2)*VBDW1(j2,j4)*(-UCAWX2(j3,j1)                 ) *I9   ! 2
//c     $ +EpsDot2 (hel2)*2*VBDWX1(j2,j4) * UCAW1 (j3,j1)                *I9s2 ! 3
//c     $ - 2*UCAWX2(j3,j1)*EpsDot12(hel1)* VBDW2 (j2,j4)                *I9s1 ! 4
//c     $ -UCAWX2(j3,j1)* VBDW1 (j2,j4)*(-EpsDot21(hel2)-2*eps2p1*DCMPLX(-ChaIni*m_e_QED))           *I9  ! 5
//c     $ + eps1D2*(-2*UCAWX2(j3,j1))*        2* VBDWX1(j2,j4)                *I9  *DCMPLX(-ChaIni*m_e_QED) ! 6
//c     $ + UCAW2(j3,j1)*(2*VBDWX2(j2,j4)              )*EpsDot12(hel1)  *I9s1 ! 7
//c     $ + UCAW2(j3,j1)*(                VBDWX1(j2,j4))*EpsDot12(hel1)  *I9   ! 7
//!!!!     $ + UCAW2(j3,j1)*(t2-2*p2p1)*VBDW1(j2,j4)                            *I9  *DCMPLX(-ChaIni*m_e_QED) ! 8 (out)
//c     $ + UCAW2(j3,j1)*(1D0/CPF1-2*p2p1)*VBDW1(j2,j4)                       *I9  *DCMPLX(-ChaIni*m_e_QED) ! 8   'Higgs' contrib added
//c     $ + UCAW2(j3,j1)*VBDWX1(j2,j4)*(-EpsDot1(hel1)+2*eps1p2*DCMPLX(-ChaIni*m_e_QED))            *I9   ! 9
//c     $                                                                            )
cFactZ= CNorm*0.5 *(-gI)*CPF*CPF2*CPF12*WVPi0;
Amp2to4add( m_AmpExpo2, cFactZ*(-2.0)*EpsDot2[Hel2] *double(I9s1), VBDW1, UCAWX1);   // ! 1 (j2,j4)*(j3,j1)
Amp2to4add( m_AmpExpo2, cFactZ*(-1.0)*EpsDot2[Hel2] *double(I9),   VBDW1, UCAWX2);   // ! 2
Amp2to4add( m_AmpExpo2, cFactZ*( 2.0)*EpsDot2[Hel2] *double(I9s1), VBDWX1, UCAW1);   // ! 3
Amp2to4add( m_AmpExpo2, cFactZ*(-2.0)*EpsDot12[Hel1]*double(I9s2), VBDW2, UCAWX2);   // ! 4
Amp2to4add( m_AmpExpo2, cFactZ*( EpsDot21[Hel2]+2.0*eps2p1*(-gI))*double(I9),   VBDW1,  UCAWX2); //! 5
Amp2to4add( m_AmpExpo2, cFactZ*(-2.0)*eps1D2*2.0*(-gI)*double(I9),    VBDWX1, UCAWX2);           //! 6
Amp2to4add( m_AmpExpo2, cFactZ*( 2.0)*EpsDot12[Hel1]  *double(I9s2),  VBDWX2, UCAW2);            //! 7
Amp2to4add( m_AmpExpo2, cFactZ*( 1.0)*EpsDot12[Hel1]  *double(I9),    VBDWX1, UCAW2);            //! 7
Amp2to4add( m_AmpExpo2, cFactZ*((1.0/CPF1-2.0*Cp2p1 ))*(-gI)*double(I9),      VBDW1, UCAW2);     //! 8
Amp2to4add( m_AmpExpo2, cFactZ*(-EpsDot1[Hel1]+2.0*eps1p2*(-gI))*double(I9), VBDWX1, UCAW2);     //! 9
//c===========================
//c      AmpWork(j1,j2,j3,j4) =AmpWork(j1,j2,j3,j4)   ! terms due to 4 boson coupling
//c     $                 +CNorm*(1D0) *(m_e_QED**2)*CPF*CPF12*WVpi0        !!!+CNorm*(-1D0) *DCMPLX(0,-m_e_QED**2)*CPF*CPF12*WVpi0
//c     $                 *(-BornABCD(j1,j2,j3,j4)/CPF0/WVpi0*2*eps1D2         *IVI
//c     $                      +0.5D0 *UCAW1(j3,j1)*VBDW2(j2,j4)               *I10
//c     $                      +0.5D0 *UCAW2(j3,j1)*VBDW1(j2,j4)               *I10
//c     $                                                                )
cFactZ= CNorm *m_e_QED*m_e_QED *CPF*CPF12*WVPi0;
AmpAdd(     m_AmpExpo2,  -cFactZ /CPF0/WVPi0*2.0*eps1D2*dcmplx(IVI)  ,BornABCD);   // IVI
Amp2to4add( m_AmpExpo2, cFactZ *(0.5)*double(I10), VBDW2, UCAW1);                  // I10
Amp2to4add( m_AmpExpo2, cFactZ *(0.5)*double(I10), VBDW1, UCAW2);                  // I10
//c===========================
//c      AmpWork(j1,j2,j3,j4) =AmpWork(j1,j2,j3,j4)+ CNorm*(-0.5D0) *DCMPLX(ChaIni*m_e_QED)*WVpi0 *(  ! terms abcd due rem. single WWgam coupl (redistr)
//c     $  +( UCAW1(j3,j1)*( 2*VBDWX1(j2,j4) ) -VBDW1(j2,j4)*( 2*UCAWX1(j3,j1)             ) )*sA(2,Hel2)      *CPF2*CPF12 *I9s2   ! (a) plus below
//c     $  +( UCAW1(j3,j1)*( 2*VBDWX1(j2,j4)               )-VBDW1(j2,j4)*( 2*UCAWX1(j3,j1)) )*sB(2,Hel2)      *CPF *CPF1  *I9s2   ! (b) plus below
//c     $  +( UCAW2(j3,j1)*( 2*VBDWX2(j2,j4) ) -VBDW2(j2,j4)*( 2*UCAWX2(j3,j1)             ) )*sA(1,Hel1)      *CPF1*CPF12 *I9s1   ! (c) plus below
//c     $  +( UCAW2(j3,j1)*( 2*VBDWX2(j2,j4)               )-VBDW2(j2,j4)*( 2*UCAWX2(j3,j1)) )*sB(1,Hel1)      *CPF*CPF2   *I9s1   ! (d) plus below
//c     $                                      -VBDW1(j2,j4)*(                  UCAW2(j3,j1) )*DCMPLX( m_e_QED)*CPF2*CPF12 *I9B    ! (a) plus  I72b
//c     $  +( UCAW1(j3,j1)*(                 VBDW2(j2,j4)  )                                 )*DCMPLX(-m_e_QED)*CPF *CPF1  *I9B    ! (b) plus  I72b
//c     $                                      -VBDW2(j2,j4)*(                  UCAW1(j3,j1) )*DCMPLX( m_e_QED)*CPF1*CPF12 *I9B    ! (c) plus  I71b
//c     $  +( UCAW2(j3,j1)*(                 VBDW1(j2,j4) )                                  )*DCMPLX(-m_e_QED)*CPF*CPF2   *I9B    ! (d) plus  I71b
//c     $                                                                         )
cFactZ= CNorm*(-0.5) *gI*WVPi0;
Amp2to4add( m_AmpExpo2, cFactZ *(+2.0)*sA[1][Hel2]*CPF2*CPF12*double(I9s2), VBDWX1, UCAW1);           // (a) plus below
Amp2to4add( m_AmpExpo2, cFactZ *(-2.0)*sA[1][Hel2]*CPF2*CPF12*double(I9s2), VBDW1, UCAWX1);           // (a) plus below
Amp2to4add( m_AmpExpo2, cFactZ *(+2.0)*sB[1][Hel2]*CPF *CPF1 *double(I9s1), VBDWX1, UCAW1);           // (b) plus below
Amp2to4add( m_AmpExpo2, cFactZ *(-2.0)*sB[1][Hel2]*CPF *CPF1 *double(I9s2), VBDW1, UCAWX1);           // (b) plus below
Amp2to4add( m_AmpExpo2, cFactZ *(+2.0)*sA[0][Hel1]*CPF1*CPF12*double(I9s2), VBDWX2, UCAW2);           // (c) plus below
Amp2to4add( m_AmpExpo2, cFactZ *(-2.0)*sA[0][Hel1]*CPF1*CPF12*double(I9s2), VBDW2, UCAWX2);           // (c) plus below
Amp2to4add( m_AmpExpo2, cFactZ *(+2.0)*sB[0][Hel1]*CPF *CPF2 *double(I9s1), VBDWX2, UCAW2);           // (d) plus below
Amp2to4add( m_AmpExpo2, cFactZ *(-2.0)*sB[0][Hel1]*CPF *CPF2 *double(I9s2), VBDW2, UCAWX2);           // (d) plus below
Amp2to4add( m_AmpExpo2, cFactZ *(-1.0)*dcmplx( m_e_QED)*CPF2*CPF12 *double(I9B), VBDW1, UCAW2);
Amp2to4add( m_AmpExpo2, cFactZ *(+1.0)*dcmplx(-m_e_QED)*CPF *CPF1  *double(I9B), VBDW2, UCAW1);
Amp2to4add( m_AmpExpo2, cFactZ *(-1.0)*dcmplx( m_e_QED)*CPF1*CPF12 *double(I9B), VBDW2, UCAW1);
Amp2to4add( m_AmpExpo2, cFactZ *(+1.0)*dcmplx(-m_e_QED)*CPF *CPF2  *double(I9B), VBDW1, UCAW2);

//[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
//(*m_Out)<<"////////////////////////////////////////////////////////////////////////////////////////////////////"<<endl;
//(*m_Out)<<"/////////////////////////////////////////GPS_HiiPlusW///////////////////////////////////////////////"<<endl;
//for(int j1=0; j1<=1; j1++) for(int j2=0; j2<=1; j2++){ (*m_Out)<< "m_AmpExpo2("<<j1<<","<<j2<<",*,*)=  ";
//     for(int j=0;j<=1;j++) for(int k=0;k<=1;k++) (*m_Out)<<"  "<<SW208<<m_AmpExpo2.m_A[j1][j2][j][k];(*m_Out)<<endl;}
//]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]

}//GPS_HiiPlusW

void KKceex::GPS_HifPlus(dcmplx CNorm, int KFini, int KFfin, TLorentzVector PX,
                         KKpart ph1, int Hel1, KKpart ph2, int Hel2){
//SUBROUTINE GPS_HifPlus(CNorm,KFi,KFf,PX,pA,mA,pB,mB,pC,mC,pD,mD,Hel1,ph1,Hel2,ph2,mph,AmpWork) !
/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
//   Genuine IR-finite non 1-photon amplitudes for ISR-FSR are added to AmpWork    //
//   1-st photon in Initial  state,  symmetrisation 1<-->2 is required!            //
//                                                                                 //
//   For Y_IR=1 IR-finite 1-photon contributions are calculated here as well       //
//                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
//                              1                2                                 //
//                              |                |                                 //
//                              |                |                                 //
//               c              |    OOOOOOOOOOOOOOOOOO          d                 //
//     u  -------<------------- | ---OOOOOOOOOOOOOOOOOO----------<----- v          //
//                              |         |                                        //
//                              |         |X                                       //
//                              |         |                                        //
//      _       -b          OOOOOOOOOOOOOOOOOOOOO                a                 //
//      v  ------<----------OOOOOOOOOOOOOOOOOOOOO----------------<----- u          //
//                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////
  int  Y_IR=0;                  // No, IR not included
  int N_IR=1-Y_IR;
  double ChaFin =  DB->Qf[ KFfin];
  double ChaIni =  DB->Qf[ KFini];
  dcmplx gI = dcmplx(ChaIni*m_e_QED);
  dcmplx gF = dcmplx(ChaFin*m_e_QED);
  KKpart pA=m_p1;
  KKpart pB=m_p2;
  KKpart pC=m_p3;
  KKpart pD=m_p4;
  dcmplx Sini[2],Sfin[2];
  Sini[0]  =  gI *Soft(  1,ph1,pA,pB);
  Sfin[0]  = -gF *Soft(  1,ph2,pC,pD);
  Sini[1] = -conj(Sini[0]);
  Sfin[1] = -conj(Sfin[0]);
// Calculate Born spin amplitudes
  double CosThe=0.0;
  KKcmplx4 BornABCD, Born1BCD, BornA1CD, BornAB2D;
  KKcmplx4 BornABC2, Born1B2D, Born1BC2, BornA12D, BornA1C2;
  Born(KFini,KFfin,PX, CosThe, pA,  pB,  pC,  pD,  BornABCD);   // standard
  Born(KFini,KFfin,PX, CosThe, ph1, pB,  pC,  pD,  Born1BCD);   // A->1
  Born(KFini,KFfin,PX, CosThe, pA,  ph1, pC,  pD,  BornA1CD);   // B->1
  Born(KFini,KFfin,PX, CosThe, pA,  pB,  ph2, pD,  BornAB2D);   // C->2
  Born(KFini,KFfin,PX, CosThe, pA,  pB,  pC,  ph2, BornABC2);   // D->2
  Born(KFini,KFfin,PX, CosThe, ph1, pB,  ph2, pD,  Born1B2D);   // A->1,C->2
  Born(KFini,KFfin,PX, CosThe, ph1, pB,  pC,  ph2, Born1BC2);   // A->1,D->2
  Born(KFini,KFfin,PX, CosThe, pA,  ph1, ph2, pD,  BornA12D);   // B->1,C->2
  Born(KFini,KFfin,PX, CosThe, pA,  ph1, pC,  ph2, BornA1C2);   // B->1,D->2
// Propagators
  KKpart QQ=pC; QQ+=pD;     // pC+pD
  KKpart PP2=QQ; PP2+=ph2;  // pC+pD+ph2
  double svarX2  =  PP2*PP2;
  double svarQ   =  QQ*QQ;
// Fermion propagarotors ini
  double prA1=  1.0/(pA*ph1)/2.0;
  double prB1= -1.0/(pB*ph1)/2.0;
// Fermion propagarotors fin
  double prC2=  1.0/(pC*ph2)/2.0;
  double prD2= -1.0/(pD*ph2)/2.0;
//
  KKcmplx2 U11a, Vb11, Uc22, V22d;
  int Sig = 1-2*Hel1;
  GPS_MakeU( gI, ph1,Sig, ph1, pA,  U11a); // <1|[1]|a>
  GPS_MakeV( gI, ph1,Sig, pB,  ph1, Vb11); // <b|{1}|1>
  Sig = 1-2*Hel2;
  GPS_MakeU( gF, ph2,Sig, pC,  ph2, Uc22); // <c|[2]|2>
  GPS_MakeV( gF, ph2,Sig, ph2, pD,  V22d); // <2|{2}|d>

/////////////////////////////////////////////////////////////////////////////////////
//                    IR divergent part for tests only
/////////////////////////////////////////////////////////////////////////////////////
//               c                |2                             d                 //
//      u  ------<-------SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS-----         //
//                                        |                                        //
//                                        |X             |1                        //
//      _       -b                        |     a+m-1    |       a                 //
//      v  ------<------------------------O--------------O-------<----- u          //
/////////////////////////////////////////////////////////////////////////////////////
//      Su1=Su1 +Born1BCD( j,j2,j3,j4) *U11a(j,j1)*prA1 *Sfin(Hel2)*Y_IR !<b|X|1[1]|a><c|(+2)|X|(+2)|d>
/////////////////////////////////////////////////////////////////////////////////////
//               c                |2                            -d                 //
//      u  ------<-------SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS-----<----- v          //
//                                        |                                        //
//                          |1            |X                                       //
//      _       -b          |   -b+m+1    |                      a                 //
//      v  ------<----------O-------------O----------------------<----- u          //
/////////////////////////////////////////////////////////////////////////////////////
//      Su1=Su1 +Vb11(j2,j)*prB1 *BornA1CD(j1,j,j3,j4) *Sfin(Hel2)*Y_IR !<b|[1]1|X|a><c|(+2)|X|(+2)|d>
/////////////////////////////////////////////////////////////////////////////////////
//                         |2                                                      //
//               c         |   c+m+2                            -d                 //
//      u  ------<---------O--------------O----------------------<----- v          //
//                                        |X                                       //
//      _       -b                        |       |1             a                 //
//      v  ------<------SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS-----<----- u          //
/////////////////////////////////////////////////////////////////////////////////////
//      Su1=Su1 +Uc22(j3,j)*prC2 *BornAB2D(j1,j2,j,j4) *Sini(Hel1)*Y_IR !<b|(+1)|X|(+1)|a><c|[2]2|X|d>
/////////////////////////////////////////////////////////////////////////////////////
//                                                       |2                        //
//               c                            -d+m+2     |      -d                 //
//      u  ------<------------------------O----------------------<----- v          //
//                                        |X                                       //
//      _       -b                        |       |1             a                 //
//      v  ------<------SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS-----<----- u          //
/////////////////////////////////////////////////////////////////////////////////////
//      Su1=Su1 +BornABC2(j1,j2,j3,j) *V22d(j,j4)*prD2  *Sini(Hel1)*Y_IR !<b|(+1)|X|(+1)|a><c|X|2[2]d>


/////////////////////////////////////////////////////////////////////////////////////
//                  IR finite parts
/////////////////////////////////////////////////////////////////////////////////////
//                         |2                                                      //
//               c         |   c+m+2                             d                 //
//      u  ------<---------O--------------O----------------------<----- v          //
//                                        |                                        //
//                                        |X             |1                        //
//      _       -b                        |     a+m-1    |       a                 //
//      v  ------<------------------------O--------------O-------<----- u          //
/////////////////////////////////////////////////////////////////////////////////////
//      Su2=Su2 +Uc22(j3,l )*prC2 *Born1B2D(j,j2,l,j4) *U11a(j,j1)*prA1 !<b|X|1[1]|a><c|[2]2|X|d>
  Amp4Zer(m_AmpTemp);
  AmpAddF(m_AmpTemp,         dcmplx(prC2), Uc22, Born1B2D);
  AmpAddI(m_AmpExpo2, CNorm* dcmplx(prA1), m_AmpTemp, U11a);
/////////////////////////////////////////////////////////////////////////////////////
//                                                       |2                        //
//               c                             -d+m-2    |       d                 //
//      u  ------<------------------------O--------------O-------<----- v          //
//                                        |                                        //
//                                        |X             |1                        //
//      _       -b                        |     a+m-1    |       a                 //
//      v  ------<------------------------O--------------O-------<----- u          //
/////////////////////////////////////////////////////////////////////////////////////
//      Su2=Su2 +Born1BC2(j,j2,j3,l) *U11a(j,j1)*prA1  *V22d(l,j4)*prD2 !<b|X|1[1]|a><c|X|2{2}|d>
  Amp4Zer(m_AmpTemp);
  AmpAddI(m_AmpTemp,         dcmplx(prA1), Born1BC2,  U11a);
  AmpAddF(m_AmpExpo2, CNorm* dcmplx(prD2), m_AmpTemp, V22d);
/////////////////////////////////////////////////////////////////////////////////////
//                       |2                                                        //
//               c       |    c+m+2                              d                 //
//      u  ------<-------O----------------O----------------------<----- v          //
//                                        |                                        //
//                       |1               |X                                       //
//      _       -b       |   -b+m+1       |                      a                 //
//      v  ------<-------O----------------O----------------------<----- u          //
/////////////////////////////////////////////////////////////////////////////////////
//      Su2=Su2 +Vb11(j2,j)*prB1  *Uc22(j3,l )*prC2 *BornA12D(j1,j,l,j4)!<b|{1}1|X|a><c|[2]2|X|d>
  Amp4Zer(m_AmpTemp);
  AmpAddF(m_AmpTemp,         dcmplx(prC2), Uc22, BornA12D);
  AmpAddI(m_AmpExpo2, CNorm* dcmplx(prB1), Vb11, m_AmpTemp);
/////////////////////////////////////////////////////////////////////////////////////
//                                                       |2                        //
//               c                             -d+m-2    |       d                 //
//      u  ------<------------------------O--------------O-------<----- v          //
//                                        |                                        //
//                       |1               |X                                       //
//      _       -b       |   -b+m+1       |                      a                 //
//      v  ------<-------O----------------O----------------------<----- u          //
/////////////////////////////////////////////////////////////////////////////////////
//      Su2=Su2 +Vb11(j2,j)*prB1  *BornA1C2(j1,j,j3,l) *V22d(l,j4)*prD2 !<b|{1}1|X|a><c|X|2[2]|d>
  Amp4Zer(m_AmpTemp);
  AmpAddI(m_AmpTemp,         dcmplx(prB1), Vb11, BornA1C2);
  AmpAddF(m_AmpExpo2, CNorm* dcmplx(prD2), m_AmpTemp, V22d);

//////////////////////////////////////////////////////////////////////////////
//                    sProd = Sini(Hel1)* Sfin(Hel2)
//                    AmpWork(j1,j2,j3,j4) =AmpWork(j1,j2,j3,j4)
//       $                 +CNorm*( Su1+Su2 )
//       $                 +CNorm*BornABCD(j1,j2,j3,j4)*sProd  *Y_IR !
/////////////////////////////////////////////////////////////////////////////////////
  dcmplx sProd = Sini[Hel1]* Sfin[Hel2];
  if(Y_IR ) AmpAdd(m_AmpExpo2, CNorm*sProd, BornABCD);

}//GPS_HifPlus


void KKceex::Amp2WCplus(KKcmplx2 &V, const double Sw2){
//    Cnor   = sqrt(1.d0/2.d0/m_Sw2)
//    V[0,1] = Cnor*0D0 !(+-)
//    V(2,1) = Cnor*0D0 !(-+)
//    V(2,2) = Cnor*0D0 !(--)
//    V(1,1) = Cnor*1D0 !(++)
    V.m_A[0][1] = dcmplx(0.0); // (+-)
    V.m_A[1][0] = dcmplx(0.0); // (-+)
    V.m_A[1][1] = dcmplx(0.0); // (--)
    V.m_A[0][0] = dcmplx( sqrt(1.0/2.0/Sw2)); //(++)
}//KKceex::Amp2WCplus

void KKceex::Amp2WCminus(KKcmplx2 &V, const double Sw2){
//    Cnor   = sqrt(1.d0/2.d0/m_Sw2)
//    V(1,2) = Cnor*0D0 !(+-)
//    V(2,1) = Cnor*0D0 !(-+)
//    V(2,2) = Cnor*1D0 !(--)
//    V(1,1) = Cnor*0D0 !(++)
    V.m_A[0][1] = dcmplx(0.0); // (+-)
    V.m_A[1][0] = dcmplx(0.0); // (-+)
    V.m_A[1][1] = dcmplx( sqrt(1.0/2.0/Sw2)); // (--)
    V.m_A[0][0] = dcmplx(0.0); //(++)
}//KKceex::Amp2WCminus


void KKceex::Amp2to4(KKcmplx4 &Res, dcmplx Fact, const KKcmplx2 &V, const KKcmplx2 &U){
// Born-like 4-dim. tensor product of two 2-dim. transition matrices
for(int j1 = 0; j1<=1; j1++)
  for(int j2 = 0; j2<=1; j2++)
      for(int j3 = 0; j3<=1; j3++)
         for(int j4 = 0; j4<=1; j4++){
           Res.m_A[j1][j2][j3][j4] = Fact* V.m_A[j2][j4] *U.m_A[j3][j1];
//           VBDWX2(j2,j4) *UC1W2( j3,j ) *U11a(j,j1)
  }// j1,j2,j3,j4
}//KKceex::Amp2to4

void KKceex::Amp2to4add(KKcmplx4 &Res, dcmplx Fact, const KKcmplx2 &V, const KKcmplx2 &U){
// Born-like 4-dim. tensor product of two 2-dim. transition matrices
// Additive version
for(int j1 = 0; j1<=1; j1++)
  for(int j2 = 0; j2<=1; j2++)
      for(int j3 = 0; j3<=1; j3++)
         for(int j4 = 0; j4<=1; j4++){
           Res.m_A[j1][j2][j3][j4] += Fact* V.m_A[j2][j4] *U.m_A[j3][j1];
  }// j1,j2,j3,j4
}//KKceex::Amp2to4

void KKceex::Amp2Mult(KKcmplx2 &Res, dcmplx Fact, const KKcmplx2 &V, const KKcmplx2 &U){
// Res = Fact* V * U
dcmplx Csum;
for(int j1 = 0; j1<=1; j1++)
  for(int j2 = 0; j2<=1; j2++){
       Csum=dcmplx(0.0,0.0);
       for(int  j= 0; j<=1; j++)
          Csum += (V.m_A[j1])[j] *(U.m_A[j])[j2];
       (Res.m_A[j1])[j2] = Fact*Csum;
  }// j1,j2
}//KKceex::Amp2Mult

void KKceex::Amp4Zer(KKcmplx4 &Born){
for(int j1 = 0; j1<=1; j1++)
  for(int j2 = 0; j2<=1; j2++)
    for(int j3 = 0; j3<=1; j3++)
      for(int j4 = 0; j4<=1; j4++)
         (Born.m_A)[j1][j2][j3][j4]= dcmplx(0.0,0.0);
}//AmpZero


void KKceex::AmpAdd(KKcmplx4 &Work, dcmplx Fact, const KKcmplx4 &Born){
dcmplx Csum;
for(int j1 = 0; j1<=1; j1++)
  for(int j2 = 0; j2<=1; j2++)
    for(int j3 = 0; j3<=1; j3++)
      for(int j4 = 0; j4<=1; j4++){
            Work.m_A[j1][j2][j3][j4] +=  Fact *(Born.m_A)[j1][j2][j3][j4];
      }//for
}// KKceex::AmpAddI

void KKceex::AmpAdd(KKcmplx4 &Work, dcmplx Fact,
       const KKcmplx2 &U, const KKcmplx2 &VX, const KKcmplx2 &V, const KKcmplx2 &UX){
//$  *(U(j3,j1)*VX(j2,j4)-V(j2,j4)*UX(j3,j1)) ! non-infrared part of emission from W
for(int j1 = 0; j1<=1; j1++)
  for(int j2 = 0; j2<=1; j2++)
    for(int j3 = 0; j3<=1; j3++)
      for(int j4 = 0; j4<=1; j4++){
         Work.m_A[j1][j2][j3][j4] += Fact*((U.m_A)[j3][j1]*(VX.m_A)[j2][j4]-(V.m_A)[j2][j4]*(UX.m_A)[j3][j1]);
      }
}//KKceex::AmpAdd

void KKceex::AmpAddI(KKcmplx4 &Work, dcmplx Fact, const KKcmplx4 &Born, const KKcmplx2 &U){
dcmplx Csum;
for(int j1 = 0; j1<=1; j1++)
  for(int j2 = 0; j2<=1; j2++)
    for(int j3 = 0; j3<=1; j3++)
      for(int j4 = 0; j4<=1; j4++){
            Csum=dcmplx(0.0,0.0);
            for(int  j= 0; j<=1; j++){
               Csum=Csum +(Born.m_A)[ j][j2][j3][j4] * (U.m_A)[j][j1];
            }//for
            Work.m_A[j1][j2][j3][j4] += Fact*Csum;
      }//for
}// KKceex::AmpAddI

void KKceex::AmpAddI(KKcmplx4 &Work, dcmplx Fact, const KKcmplx2 &V, const KKcmplx4 &Born){
dcmplx Csum;
for(int j1 = 0; j1<=1; j1++)
  for(int j2 = 0; j2<=1; j2++)
    for(int j3 = 0; j3<=1; j3++)
      for(int j4 = 0; j4<=1; j4++){
            Csum=dcmplx(0.0,0.0);
            for(int  j= 0; j<=1; j++){
               Csum=Csum +(V.m_A[j2])[j] * (Born.m_A)[j1][ j][j3][j4];
            }//for
            Work.m_A[j1][j2][j3][j4] += Fact*Csum;
      }//for
}// KKceex::AmpAddI

void KKceex::AmpAddI(KKcmplx4 &Work, dcmplx Fact, const KKcmplx2 &V, const KKcmplx4 &Born, const KKcmplx2 &U){
// V * Born * U
dcmplx Csum;
for(int j1 = 0; j1<=1; j1++)
  for(int j2 = 0; j2<=1; j2++)
    for(int j3 = 0; j3<=1; j3++)
      for(int j4 = 0; j4<=1; j4++){
            Csum=dcmplx(0.0,0.0);
            for(int  j= 0; j<=1; j++)
                for(int  l= 0; l<=1; l++){
                Csum=Csum +(V.m_A[j2])[j] * (Born.m_A)[l][j][j3][j4] *(U.m_A[l])[j1];
            }//for j,l
            Work.m_A[j1][j2][j3][j4] += Fact*Csum;
      }//for j1...j4
}// KKceex::AmpAddI

void KKceex::AmpAddI(KKcmplx4 &Work, dcmplx Fact, const KKcmplx4 &Born, const KKcmplx2 &U1, const KKcmplx2 &U2){
// Born * U1 * U2
dcmplx Csum;
for(int j1 = 0; j1<=1; j1++)
  for(int j2 = 0; j2<=1; j2++)
    for(int j3 = 0; j3<=1; j3++)
      for(int j4 = 0; j4<=1; j4++){
            Csum=dcmplx(0.0,0.0);
            for(int  j= 0; j<=1; j++)
                for(int  l= 0; l<=1; l++){
                Csum=Csum + (Born.m_A)[l][j2][j3][j4] * (U1.m_A[l])[j] *(U2.m_A[j])[j1];
            }//for j,l
            Work.m_A[j1][j2][j3][j4] += Fact*Csum;
      }//for j1...j4
}// KKceex::AmpAddI


void KKceex::AmpAddI(KKcmplx4 &Work, dcmplx Fact, const KKcmplx2 &V1, const KKcmplx2 &V2, const KKcmplx4 &Born){
// Born * U1 * U2
dcmplx Csum;
for(int j1 = 0; j1<=1; j1++)
  for(int j2 = 0; j2<=1; j2++)
    for(int j3 = 0; j3<=1; j3++)
      for(int j4 = 0; j4<=1; j4++){
            Csum=dcmplx(0.0,0.0);
            for(int  j= 0; j<=1; j++)
                for(int  l= 0; l<=1; l++){
                Csum=Csum +  (V1.m_A[j2])[j] *(V2.m_A[j])[l] *(Born.m_A)[j1][l][j3][j4];
            }//for j,l
            Work.m_A[j1][j2][j3][j4] += Fact*Csum;
      }//for j1...j4
}// KKceex::AmpAddI



void KKceex::AmpAddF(KKcmplx4 &Work, dcmplx Fact, const KKcmplx4 &Born, const KKcmplx2 &V){
dcmplx Csum;
for(int j1 = 0; j1<=1; j1++)
  for(int j2 = 0; j2<=1; j2++)
    for(int j3 = 0; j3<=1; j3++)
      for(int j4 = 0; j4<=1; j4++){
            Csum=dcmplx(0.0,0.0);
            for(int  j= 0; j<=1; j++){
                Csum=Csum + Born.m_A[j1][j2][j3][ j] *V.m_A[j][j4];
        //      Csum=Csum + Born.m_A[j1][j2][j3][ j]* V.m_A[j][j4];
            }//for
            Work.m_A[j1][j2][j3][j4] += Fact*Csum;
      }//for
}// KKceex::AmpAddI

void KKceex::AmpAddF(KKcmplx4 &Work, dcmplx Fact, const KKcmplx2 &U, const KKcmplx4 &Born){
dcmplx Csum;
for(int j1 = 0; j1<=1; j1++)
  for(int j2 = 0; j2<=1; j2++)
    for(int j3 = 0; j3<=1; j3++)
      for(int j4 = 0; j4<=1; j4++){
            Csum=dcmplx(0.0,0.0);
            for(int  j= 0; j<=1; j++){
                Csum=Csum +U.m_A[j3][j] * Born.m_A[j1][j2][ j][j4];
            }//for
            Work.m_A[j1][j2][j3][j4] += Fact*Csum;
      }//for
}// KKceex::AmpAddI


void KKceex::AmpAddF(KKcmplx4 &Work, dcmplx Fact, const KKcmplx2 &U, const KKcmplx4 &Born, const KKcmplx2 &V){
// V * Born * U
dcmplx Csum;
for(int j1 = 0; j1<=1; j1++)
  for(int j2 = 0; j2<=1; j2++)
    for(int j3 = 0; j3<=1; j3++)
      for(int j4 = 0; j4<=1; j4++){
            Csum=dcmplx(0.0,0.0);
            for(int  j= 0; j<=1; j++)
                for(int  l= 0; l<=1; l++){
                Csum=Csum +(U.m_A[j3])[j] * (Born.m_A)[j1][j2][j][l] *(V.m_A[l])[j4];
            }//for j,l
            Work.m_A[j1][j2][j3][j4] += Fact*Csum;
      }//for j1...j4
}// KKceex::AmpAddF

void KKceex::AmpAddF(KKcmplx4 &Work, dcmplx Fact, const KKcmplx4 &Born, const KKcmplx2 &V1, const KKcmplx2 &V2){
// Born * U1 * U2
dcmplx Csum;
for(int j1 = 0; j1<=1; j1++)
  for(int j2 = 0; j2<=1; j2++)
    for(int j3 = 0; j3<=1; j3++)
      for(int j4 = 0; j4<=1; j4++){
            Csum=dcmplx(0.0,0.0);
            for(int  j= 0; j<=1; j++)
                for(int  l= 0; l<=1; l++){
                Csum=Csum + (Born.m_A)[j1][j2][j3][l] * (V1.m_A[l])[j] *(V2.m_A[j])[j4];
            }//for j,l
            Work.m_A[j1][j2][j3][j4] += Fact*Csum;
      }//for j1...j4
}// KKceex::AmpAddF


void KKceex::AmpAddF(KKcmplx4 &Work, dcmplx Fact, const KKcmplx2 &U1, const KKcmplx2 &U2, const KKcmplx4 &Born){
// U1 * U2 *Born
dcmplx Csum;
for(int j1 = 0; j1<=1; j1++)
  for(int j2 = 0; j2<=1; j2++)
    for(int j3 = 0; j3<=1; j3++)
      for(int j4 = 0; j4<=1; j4++){
            Csum=dcmplx(0.0,0.0);
            for(int  j= 0; j<=1; j++)
                for(int  l= 0; l<=1; l++){
                Csum=Csum +  (U1.m_A[j3])[l] *(U2.m_A[l])[j] *(Born.m_A)[j1][j2][j][j4];
            }//for j,l
            Work.m_A[j1][j2][j3][j4] += Fact*Csum;
      }//for j1...j4
}// KKceex::AmpAddF




void KKceex::EWFFact(int KFini, int KFfin, double Svar,
                     dcmplx &Ve, dcmplx &Vf, dcmplx &Ae, dcmplx &Af ){
////////////////////////////////////////////////////////////////////////////////////////////
//                   ElectroWeak Couplings for effective Born                             //
////////////////////////////////////////////////////////////////////////////////////////////
//   GamVPi = DCMPLX(1d0) ZetVPi = DCMPLX(1d0) VVCor  = DCMPLX(1d0)
//   RSQV=1d0   RSQA=1d0
//===============================================================================
// Get charges, izospin, color
double Qe  = DB->Qf[ KFini];
double T3e = DB->T3f[KFini];
double Qf  = DB->Qf[ KFfin];
double T3f = DB->T3f[KFfin];
//
double Sw2 = m_DZ->D_swsq;       // from Dizet
Sw2 = DB->swsq;                  //from xpar ????
// Couplings costants
double Deno   = sqrt(16.0*Sw2*(1.0-Sw2));
Ve     = (2*T3e -4*Qe*Sw2)/Deno;
Vf     = (2*T3f -4*Qf*Sw2)/Deno;
Ae     =  2*T3e           /Deno;
Af     =  2*T3f           /Deno;
}//EWFFact LO


void KKceex::EWFFact(const int KFini,const  int KFfin,const  double Svar,const  double CosThetD,
             dcmplx &Ve,    dcmplx &Vf,     dcmplx  &Ae,     dcmplx &Af,
             dcmplx &VVcor, dcmplx &GamVPi, dcmplx  &ZetVPi, double &RsqV, double & RsqA){
////////////////////////////////////////////////////////////////////////////////////////////
//                        ElectroWeak Corrections                                         //
//  They are in Vector Couplings (multiplied by correcting f-factors)                     //
//  Because of cost(theta) depenedence of WW boxes we need to define CosThetD variable    //
////////////////////////////////////////////////////////////////////////////////////////////
//===============================================================================
// Get charges, izospin, color
double Qe  = DB->Qf[ KFini];
double T3e = DB->T3f[KFini];
double Qf  = DB->Qf[ KFfin];
double T3f = DB->T3f[KFfin];
double Sw2 = m_DZ->D_swsq;       // from Dizet
// Get EW form-factors
// CALL GPS_GSWimport(KFini, KFfin, Svar, CosThetD, GSW);
//g_KKhhGen->m_BornDist->DZ->InterpoGSW(*KFi, *KFf, *svar, *CosThe);
//g_KKhhGen->m_BornDist->DZ->GetGSWxy(GSW_re,GSW_im);
  m_DZ->InterpoGSW(KFini, KFfin, Svar, CosThetD);
  dcmplx RhoEW     = m_DZ->m_GSW[1-1];
  dcmplx CorEle    = m_DZ->m_GSW[2-1];
  dcmplx CorFin    = m_DZ->m_GSW[3-1];
  dcmplx CorEleFin = m_DZ->m_GSW[4-1];
  dcmplx VPgamma   = m_DZ->m_GSW[6-1];
// Vacuum polarization factors
  GamVPi = 1.0   /(2.0-VPgamma);
  ZetVPi = DB->GFermi *sqr(DB->MZ) *DB->Alfinv0 /(sqrt(2.0)*8.0*M_PI)
            *(Sw2*(1.0-Sw2)) *16.0
            *RhoEW;
//[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
//  (*m_Out)<<" KKceex::EWFFact: Svar="<<Svar<<"  CosThetD="<<CosThetD<<endl;
//  (*m_Out)<<" KKceex::EWFFact: RhoEW="<<RhoEW<<"  CorEle="<< CorEle<<endl;
//  (*m_Out)<<" KKceex::EWFFact: CorFin="<<CorFin<<"  CorEleFin="<< CorEleFin<<" VPgamma="<<VPgamma<<endl;
//]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
// Coupling costants times EW form-factors
  double Deno   = sqrt(16.0*Sw2*(1.0-Sw2));
  Ve     = (2*T3e -4*Qe*Sw2*CorEle)/Deno;
  Vf     = (2*T3f -4*Qf*Sw2*CorFin)/Deno;
  Ae     =  2*T3e             /Deno;
  Af     =  2*T3f             /Deno;
// Angle dependent double-vector extra-correction
   dcmplx VVCef  = ( (2*T3e)      *(2*T3f)
            -(4*Qe*Sw2) *(2*T3f)      *CorEle
            -(4*Qf*Sw2) *(2*T3e)      *CorFin
            +(4*Qe*Sw2) *(4*Qf*Sw2) *CorEleFin )/sqr(Deno);
  VVcor  = VVCef/(Ve*Vf);
//--- SY: The following looks like a misprint - corrected
// CosThetD = 1d0 is special
//  if( CosThetD == 0.0) VVcor  = dcmplx(1.0); //????
}//


void KKceex::GPS_EWFFactW(int KFi, int KFf,double s,double t, dcmplx &PropW, dcmplx &WVPi)
{ //*new*
////////////////////////////////////////////////////////////////////////////////////////////
//                        ElectroWeak Corrections                                         //
//                        empty prepared for nunu chanel                                  //
////////////////////////////////////////////////////////////////////////////////////////////
// we do it like for photon. I was unable to properly understand normalization
// this is half-guess half observation of the KK-KORALZ differences in code.
// it can be easily wrong. nothing like this tem is present in KORALZ,
// nonetheless KORALZ results are 137/127 higher than KK for pure W exchange.
// in koralz this vacpol factor is installed at born step mode 1/0 ratio, here not
  int KeyElw = DB->KeyElw;
  double MW  = DB->MZ *sqrt(1.0-DB->swsq);
  if( DB->KeyElw > 0) {
//         next will be mass fro common /cdzwsm/ from  DZface.f it is actually raw mass as above ...
    MW = m_DZ->D_MW;
  }
  PropW   =    1.0/dcmplx(t-MW*MW,0.0 );
//  WVPi=DCMPLX(1.06794821,-0.0185004916)
  WVPi=dcmplx(1.0,0.0);
  double Sw2 = m_DZ->D_swsq;
  double Coef  =1.0/2.0/Sw2;
  double Coef1= DB->GFermi *MW*MW* DB->Alfinv0 /(sqrt(2.0)*M_PI);
  double CosThetD;
  if( KeyElw >0) {
    WVPi= WVPi*Coef1/Coef ;  // effective coupling
    CosThetD=1+2*t/s;
// DIRTY trick to get around the problem. Anyway problem appears in non used
//           WVPi's where it is later overwritten.
//           Physically it is ambiguity of order alpha**2 ?????
    if( abs(1+2*t/s) > 1) {
       CosThetD=1.0/(1+2*t/s);
    }
    m_DZ->InterpoGSW(KFi, KFf, s, CosThetD);
    dcmplx GSW5 = m_DZ->m_GSW[5-1];
//     use of rhocc forbiden, disk-mode does not work
//           CALL rhocc(s,-t,-s-t,-1D0,1D0,0D0,0D0,ROW)
//           IF(icont.LE.100) WRITE(*,*) '%%%%%%   ',GSW(5),ROW,'   ',s,t
//           IF (ABS(row-gsw(5)).GT.0.0002.AND.s.GT.10000) WRITE(*,*) GSW(5),ROW,'   ',s,t
//           IF (ABS(row-1).gt.0.016) STOP
// this was added to GSW(5) in BornV_Dizet
    double DelW= 1.0/DB->Alfinv0/M_PI/2*(-3.0/2*log(s/MW/MW)+1.0/2*sqr(log(-t/s))-4.0/6.0*sqr(M_PI)+2.0); // (b);
//ccc       DelW= 1D0/m_AlfInv/m_pi/2*(-3D0/2*LOG(s/m_MW**2)+1D0/2*(LOG(-t/s))**2-1d0/6d0*m_pi**2+2D0) ! (a)
     WVPi= WVPi*(GSW5+dcmplx(DelW,0.0));
  }
//[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
//(*m_Out)<< "///////////////////////////////GPS_EWFFactW//////////////////////////////////////////////"<<endl;
//(*m_Out)<< "KeyElw="<< KeyElw <<endl;
//(*m_Out)<< "MW"<< MW<<" s="<< s<<" t="<<t<<endl;
//(*m_Out)<< "PropW"<< PropW <<" WVPi="<<WVPi<<endl;
//]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
//
}//GPS_EWFFactW


void KKceex::PartitionStart(int &last){
///////////////////////////////////////////////////////////////////////////////
//   Initialize first partition                                              //
///////////////////////////////////////////////////////////////////////////////
//----------
  if(     (m_KeyISR == 1) && (m_HasFSR == 1) ) {
     if(      m_KeyInt != 0) {
// Normal case, ISR+FSR
// cout<<"  KKceex::PartitionStart: Normal case, ISR+FSR"<<endl;
        for(int i = 1; i<=m_nPhot; i++){
           m_isr[i] = 0;     // Start with all FSR
        }
        last=0;              // Run through all partitions !!!
     } else {
// INTerference OFF, Copy partition from crude MC
// cout<<"  KKceex::PartitionStart: INTerference OFF, Copy partition from crude MC"<<endl;
        for(int i = 1; i<= m_nPhot; i++) {
//           m_isr[i] = m_isr0[i];     // Start with all FSR ???
           m_isr[i] = m_Event->m_isr[i];     // Start with copy of partition from crude MC
        }
        last=1;              // Exit next time
     }//if else
  } else if( (m_KeyISR == 1) && (m_HasFSR == 0) ) {
// Special test, ISR ONLY
     for(int i = 1; i<= m_nPhot; i++) {
        m_isr[i] = 1;        // Start with all ISR
     }// for
     last=1;                // Exit next time
// changed for KKMC-hh to allow ISR to be turned off
  } else if(m_KeyISR == 0) {
// Special test, FSR ONLY
     for(int i = 1; i<= m_nPhot; i++) {
        m_isr[i] = 0;        // Start with all FSR
     }//for
     last=1;                 // Exit next time
  } else {
     cout<< "#### GPS_PartitionStart: Wrong KeyISR,HasFSR = "<<m_KeyISR<<"  "<< m_HasFSR<<endl;
     exit(80);
  }//if
}//PartitionStart

void KKceex::PartitionPlus(int &last){
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//   update m_isr, check if it is last partition                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
  if(m_nPhot == 1) last=1;   //!!! Exit next time
  m_isr[1]=m_isr[1]+1;
  for(int i=1; i<= m_nPhot; i++){
     if( m_isr[i] == 2 ) {
        m_isr[i]=0;
        m_isr[i+1] = m_isr[i+1]+1;
        if( m_isr[m_nPhot] == 2 ) last=2; //!!! Immediate exit
     }//if
  }//for
}//PartitionPlus


void KKceex::GPS_MakeU(dcmplx Cfac, KKpart &ph, int sigma, KKpart &p1, KKpart &p2,  KKcmplx2 &U){
/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
//   Transition matrix U, (epsilon-slash sandwiched between ubar-u spinors)        //
//   ph      = photon  4-momentum                                                  //
//   sigma   = photon polarization (+1,-1)                                         //
//   p1,p2   = fermion 4-momenta                                                   //
//                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////
if( abs(sigma) != 1){
	cout<<"++++++++ GPS_MakeU: WRONG sigma= "<<sigma<<endl; exit(97);}
double m1 = p1.M;
double m2 = p2.M;
if( m_KeyArb == 0 ) {
  dcmplx Sqr2 = sqrt(2.0)*Cfac;
  if(     sigma == 1 ) {
     U.m_A[0][0] =  Sqr2*( XiProd(p2,ph) *iProd1(1,ph,p1));           // (++)
     U.m_A[1][1] =  Sqr2*( XiProd(p1,ph) *iProd1(1,ph,p2));           // (--)
     U.m_A[0][1] =  0.0;                                              // (+-)
     U.m_A[1][0] =  Sqr2*( -m1*XiProd(p2,p1) +m2*XiProd(p1,p2));      // (-+)
  }else  if(sigma == -1 ) {
     U.m_A[0][0] =  Sqr2*( XiProd(p1,ph) *iProd1(-1,ph,p2));          // (++)
     U.m_A[1][1] =  Sqr2*( XiProd(p2,ph) *iProd1(-1,ph,p1));          // (--)
     U.m_A[1][0] =  0.0;                                              // (-+)
     U.m_A[0][1] =  Sqr2*( -m1*XiProd(p2,p1) +m2*XiProd(p1,p2));      // (+-)
  }//if sigma
}else{
  dcmplx Cnor;
  if(     sigma == 1 ) {
     Cnor   = dcmplx(sqrt(2.0),0.0)/iProd1(-1,ph,m_b)*Cfac;
     U.m_A[0][0] = Cnor*(     iProd1(  1,p1, ph)*iProd1(-1,m_b,p2)    //(++)
                       +m1*m2*XiProd(    m_b,p1)*XiProd(   ph, p2) );
     U.m_A[1][1] = Cnor*(     iProd1( -1,p1,m_b)*iProd1( 1,ph, p2)    //(--)
                       +m1*m2*XiProd(    m_b,p2)*XiProd(   ph, p1) );
     U.m_A[0][1] = Cnor*(    +m1*XiProd( m_b,p1)*iProd1( 1,ph, p2)    //(+-)
                             +m2*XiProd( m_b,p2)*iProd1( 1,p1, ph) );
     U.m_A[1][0] = Cnor*(    +m1*XiProd(  ph,p1)*iProd1(-1,m_b,p2)    //(-+)
                             +m2*XiProd(  ph,p2)*iProd1(-1,p1,m_b) );
  }else  if(sigma  == -1 ) {
     Cnor   = dcmplx(sqrt(2.0),0.0)/iProd1( 1,ph,m_b);
     U.m_A[0][0] = Cnor*(     iProd1(  1,p1,m_b)*iProd1(-1,ph, p2)    //(++)
                       +m1*m2*XiProd(    m_b,p2)*XiProd(   ph, p1) );
     U.m_A[1][1] = Cnor*(     iProd1( -1,p1, ph)*iProd1( 1,m_b,p2)    //(--)
                       +m1*m2*XiProd(    m_b,p1)*XiProd(   ph, p2) );
     U.m_A[1][0] = Cnor*(    +m1*XiProd( m_b,p1)*iProd1(-1,ph, p2)    //(-+)
                             +m2*XiProd( m_b,p2)*iProd1(-1,p1, ph) );
     U.m_A[0][1] = Cnor*(    +m1*XiProd(  ph,p1)*iProd1( 1,m_b,p2)    //(+-)
                             +m2*XiProd(  ph,p2)*iProd1( 1,p1,m_b) );
  }// if sigma
}// if m_KeyArb
}// GPS_MakeU

void KKceex::GPS_MakeV(dcmplx Cfac, KKpart &ph, int sigma, KKpart &p1, KKpart &p2,  KKcmplx2 &V){
/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
//   Transition matrix V, (epsilon-slash sandwiched between v vbar spinors)        //
//   ph      = photon  4-momentum                                                  //
//   sigma   = photon polarization (+1,-1)                                         //
//   p1,p2   = fermion 4-momenta                                                   //
//                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////
if( abs(sigma) != 1){
     cout<<"++++++++ GPS_MakeV: WRONG sigma= "<<sigma<<endl; exit(97);}
double m1 = p1.M;
double m2 = p2.M;
if( m_KeyArb == 0 ) {
  dcmplx Sqr2 = sqrt(2.0) *Cfac;
  if(     sigma  == 1 ) {
     V.m_A[1][1] =  Sqr2*( XiProd(p2,ph) *iProd1(1,ph,p1));           //(--)
     V.m_A[0][0] =  Sqr2*( XiProd(p1,ph) *iProd1(1,ph,p2));           //(++)
     V.m_A[1][0] =  0.0;                                              //(-+)
     V.m_A[0][1] =  Sqr2*( m1*XiProd(p2,p1) -m2*XiProd(p1,p2));       //(+-)
  }else if(sigma == -1 ) {
     V.m_A[1][1] =  Sqr2*( XiProd(p1,ph) *iProd1(-1,ph,p2));          //(--)
     V.m_A[0][0] =  Sqr2*( XiProd(p2,ph) *iProd1(-1,ph,p1));          //(++)
     V.m_A[0][1] =  0.0;                                              //(+-)
     V.m_A[1][0] =  Sqr2*( m1*XiProd(p2,p1) -m2*XiProd(p1,p2));       //(-+)
  }//if sigma
} else {
  dcmplx Cnor;
  if(   sigma == 1 ) {
    Cnor   = dcmplx(sqrt(2.0),0.0)/iProd1(-1,ph,m_b) *Cfac;
    V.m_A[1][1] = Cnor*(     iProd1(  1,p1, ph)*iProd1(-1,m_b,p2)    //--)
                  +m1*m2*XiProd(    m_b,p1)*XiProd(   ph, p2) );
    V.m_A[0][0] = Cnor*(     iProd1( -1,p1,m_b)*iProd1( 1,ph, p2)    //++)
                   +m1*m2*XiProd(    m_b,p2)*XiProd(   ph, p1) );
    V.m_A[1][0] = Cnor*(    -m1*XiProd( m_b,p1)*iProd1( 1,ph, p2)    //-+)
                         -m2*XiProd( m_b,p2)*iProd1( 1,p1, ph) );
    V.m_A[0][1] = Cnor*(    -m1*XiProd(  ph,p1)*iProd1(-1,m_b,p2)    //+-)
                         -m2*XiProd(  ph,p2)*iProd1(-1,p1,m_b) );
  } else if(sigma == -1 ) {
    Cnor   = dcmplx(sqrt(2.0),0.0)/iProd1( 1,ph,m_b);
    V.m_A[1][1] = Cnor*(     iProd1(  1,p1,m_b)*iProd1(-1,ph, p2)    //--)
                 +m1*m2*XiProd(    m_b,p2)*XiProd(   ph, p1) );
    V.m_A[0][0] = Cnor*(     iProd1( -1,p1, ph)*iProd1( 1,m_b,p2)    //++)
                 +m1*m2*XiProd(    m_b,p1)*XiProd(   ph, p2) );
    V.m_A[0][1] = Cnor*(    -m1*XiProd( m_b,p1)*iProd1(-1,ph, p2)    //+-)
                       -m2*XiProd( m_b,p2)*iProd1(-1,p1, ph) );
    V.m_A[1][0] = Cnor*(    -m1*XiProd(  ph,p1)*iProd1( 1,m_b,p2)    //-+)
                       -m2*XiProd(  ph,p2)*iProd1( 1,p1,m_b) );
  }// if sigma
}//if keyarb
}// GPS_MakeV

void KKceex::GPS_MakeUW(dcmplx Cfac, KKpart &ph, int sigma, KKpart &p1, KKpart &p2,  KKcmplx2 &U){
//SUBROUTINE GPS_MakeUW(Norm,ph,sigma,p1,m1,p2,m2,U)
/////////////////////////////////////////////////////////////////////////////////////
// like GPS_MakeU  but with W v-a coupling included                                //
/////////////////////////////////////////////////////////////////////////////////////
//   CALL GPS_MakeU(ph,Sigma,  p1,m1,   p2,m2,    UW0)
//   CALL  GPS_MatrWm(WC)             ! W-e-nu couplingss
//   CALL  GPS_times(Cnor,WC,UW0,U) ! we add v-a coupl
KKcmplx2 UW0,WC;
double Sw2 = m_DZ->D_swsq;
dcmplx Cone = dcmplx(1.0);
Amp2WCminus(WC,Sw2);           // W-e-nu couplingss
GPS_MakeU(Cfac, ph,sigma, p1, p2, UW0);
Amp2Mult( U, Cone, WC, UW0);  // we add v-a coupl
}//KKceex::GPS_MakeUW

void KKceex::GPS_MakeVW(dcmplx Cfac, KKpart &ph, int sigma, KKpart &p1, KKpart &p2,  KKcmplx2 &V){
//SUBROUTINE GPS_MakeVW(Norm,ph,sigma,p1,m1,p2,m2,V)
/////////////////////////////////////////////////////////////////////////////////////
// like GPS_MakeV  but with W v-a coupling included                                //
/////////////////////////////////////////////////////////////////////////////////////
//   CALL GPS_MakeV(ph,Sigma,  p1,m1,   p2,m2,    VW0)
//   CALL  GPS_MatrW(WC)             ! W-e-nu couplingss
//   CALL  GPS_times(Cnor,VW0,WC,V) ! we add v-a coupl
KKcmplx2 VW0,WC;
double Sw2 = m_DZ->D_swsq;
dcmplx Cone = dcmplx(1.0);
Amp2WCplus(WC,Sw2);           // W-e-nu couplingss
GPS_MakeV(Cfac, ph,sigma, p1, p2, VW0);
Amp2Mult( V, Cone, VW0,WC);  // we add v-a coupl
}//KKceex::GPS_MakeVW

void KKceex::GPS_MatrS(KKpart &p1, KKpart &p2, KKcmplx2 &V){
//SUBROUTINE GPS_MatrS(p1,m1,p2,m2,V)
/////////////////////////////////////////////////////////////////////////////////////
//   S matrix for spinor products (matrix version of  GPS_iProd1                   //
//   p1,p2   = fermion 4-momenta                                                   //
/////////////////////////////////////////////////////////////////////////////////////
//   V(1,2) =      GPS_iProd1(  1,p1, p2)   !(+-)
//   V(2,1) =      GPS_iProd1( -1,p1,p2)    !(-+)
//   V(2,2) = 0D0                            !(--)
//   V(1,1) = 0D0                            !(++)
   V.m_A[0][1] = iProd1(  1, p1, p2); // (+-)
   V.m_A[1][0] = iProd1( -1, p1, p2); // (-+)
   V.m_A[1][1] = dcmplx(0.0); // (--)
   V.m_A[0][0] = dcmplx(0.0); // (++)
}//KKceex::GPS_MatrS

void KKceex::GPS_MakeUX(dcmplx Cfac, KKpart &ph, KKpart &p1, KKpart &p2,  KKcmplx2 &U){
//SUBROUTINE GPS_MakeUX(Cn,ph,mh,p1,m1,p2,m2,U)
///////////////////////////////////////////////////////////////////////////////////////
//   Transition matrix U, (ph-slash sandwiched between ubar-u spinors)             //
//   mass terms not tested   !!!                                                   //
//   ph      = photon  4-momentum                                                  //
//   p1,p2   = fermion 4-momenta                                                   //
/////////////////////////////////////////////////////////////////////////////////////
//DOUBLE COMPLEX   U(2,2),A(2,2),B(2,2),AB(2,2),C(2,2)
//Cnor=1D0
//CALL  GPS_MatrS(p1,m1,ph,mh,A)
//CALL  GPS_MatrWm(C)
//CALL  GPS_MatrS(ph,mh,p2,m2,B)
//CALL  GPS_times(Cnor,A,B,AB)
//CALL  GPS_times(Cn,C,AB,U)
KKcmplx2 A,B,AB,C;
double Sw2 = m_DZ->D_swsq;
dcmplx Cone = dcmplx(1.0);
Amp2WCminus(C,Sw2);           // W-e-nu couplingss
GPS_MatrS(p1,ph,A);
GPS_MatrS(ph,p2,B);
Amp2Mult( AB, Cone, A, B);
Amp2Mult(  U, Cfac, C,AB);
}//KKceex::GPS_MakeUX

void KKceex::GPS_MakeVX(dcmplx Cfac, KKpart &ph, KKpart &p1, KKpart &p2,  KKcmplx2 &V){
//SUBROUTINE GPS_MakeVX(Cn,ph,mh,p1,m1,p2,m2,V)
/////////////////////////////////////////////////////////////////////////////////////
//   Transition matrix V, (ph-slash sandwiched between vbar-v spinors)             //
//   mass terms not studied !!!                                                    //
//   ph      = photon  4-momentum                                                  //
//   p1,p2   = fermion 4-momenta                                                   //
//   massles limit                                                                 //
/////////////////////////////////////////////////////////////////////////////////////
//    DOUBLE COMPLEX   V(2,2),A(2,2),B(2,2),C(2,2),AB(2,2)
//    DOUBLE COMPLEX   Cn,Cnor,XX
//    Cnor=1D0
//    CALL  GPS_MatrS(p1,m1,ph,mh,A)
//    CALL  GPS_MatrWm(C)
//    CALL  GPS_MatrS(ph,mh,p2,m2,B)
//    CALL  GPS_times(Cnor,A,B,AB)
//    CALL  GPS_times(Cn,AB,C,V)
//    XX=V(1,1)
//    V(1,1)=V(2,2)
//    V(2,2)=XX
KKcmplx2 A,B,AB,C;
double Sw2 = m_DZ->D_swsq;
dcmplx Cone = dcmplx(1.0);
Amp2WCminus(C,Sw2);           // W-e-nu couplings ???? Amp2WCplus???
GPS_MatrS(p1,ph,A);
GPS_MatrS(ph,p2,B);
Amp2Mult( AB, Cone,  A, B);
Amp2Mult(  V, Cfac, AB, C);
dcmplx XX = V.m_A[0][0];
V.m_A[0][0]=V.m_A[1][1];
V.m_A[1][1]=XX;
}//KKceex::GPS_MakeUX


dcmplx KKceex::iProd1(int L, KKpart &p, KKpart &q){
/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
// This is basic inner s-product of spinors s_{L}(p,q)=ubar_{L}(p)*u_{-L}(q)       //
// We exploit identity s_{-}(p,q) = -[s_{+}(p,q)]^*                                //
// Four-vectors p,q are in principle massless, however, massive p,q                //
// are also meaningfull. Implicit projection on xi is then assumed.                //
//                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////
dcmplx  Prod;
if(     L ==  1 ) {
   Prod=
        -sqrt( ((p.P)[0]-(p.P)[1]) / ((q.P)[0]-(q.P)[1]) ) *dcmplx((q.P)[2],(q.P)[3])
        +sqrt( ((q.P)[0]-(q.P)[1]) / ((p.P)[0]-(p.P)[1]) ) *dcmplx((p.P)[2],(p.P)[3]);
} else  if( L == -1 ) {
   Prod=
        -sqrt( ((q.P)[0]-(q.P)[1]) / ((p.P)[0]-(p.P)[1]) ) *dcmplx((p.P)[2],(p.P)[3])
        +sqrt( ((p.P)[0]-(p.P)[1]) / ((q.P)[0]-(q.P)[1]) ) *dcmplx((q.P)[2],(q.P)[3]);
   Prod= conj(Prod);
}   else {
   cout<<"##### KKceex::iProd1: Wrong L= "<< L << endl; exit(90);
}
return Prod;
}//iProd1

dcmplx KKceex::iProd2(KKpart &p, KKpart &q){
/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
// General spinor product s_{l1,l2}(p,q) for massive spinors U and/or V            //
// p.M and q.M are spinor masses of p and q.                                       //
//                                                                                 //
// Antiparticle, v-spinor, is recognized according to sign of p.C, q.C             //
// Spin sign p.Hel is flipped for V-spinor.                                        //
//                                                                                 //
// Note also that p.M=0 or q.M=0 for massive p, q, mean that p, q                  //
// are projected on xi.                                                            //
//                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////
 dcmplx  Prod;
 if( abs(p.C) == 1 && abs(q.C) == 1 && abs(p.Hel) == 1 && abs(q.Hel) == 1 ){
   // Helicity conservation and non-conservation separately
   if(        p.C*p.Hel ==  q.C *q.Hel  ) {
      Prod = dcmplx( p.C *p.M *XiProd(q,p) + q.C *q.M *XiProd(p,q), 0.0);
   } else {
      Prod =  iProd1( p.C*p.Hel, p, q); // Helicity for v-spinor is fliped
   }// if p.Hel
 } else {
    cout<<"##### KKceex::iProd2: Wrong p.C="<< p.C<<" q.C= "<<q.C;
    cout<<"#####   or wrong helicities  p.Hel="<<  p.Hel<<"  q.Hel= "<<q.Hel<< endl;exit(92);
 }//if p.C
 return Prod;
}//iProd2

dcmplx KKceex::iProd2(int Hp, KKpart &p, int Hq, KKpart &q){
/////////////////////////////////////////////////////////////////////////////////////
//        The same but helicities  as arguments                                    //
/////////////////////////////////////////////////////////////////////////////////////
 dcmplx  Prod;
 if( abs(p.C) == 1 && abs(q.C) == 1 && abs(Hp) == 1 && abs(Hq) == 1 ){
   // Helicity conservation and non-conservation separately
   if(    p.C *Hp ==  q.C *Hq  ) {
      Prod = dcmplx( p.C *p.M *XiProd(q,p) + q.C *q.M *XiProd(p,q), 0.0);
   } else {
      Prod =  iProd1( p.C*Hp, p, q); // Helicity for v-spinor is fliped
   }// if Hp
 } else {
    cout<<"##### KKceex::iProd2: Wrong p.C="<< p.C<<" q.C= "<<q.C;
    cout<<"#####   or wrong helicities  Hp="<<  Hp<<"  Hq= "<<Hq<< endl;exit(92);
 }//if p.C
 return Prod;
}//iProd2

/*
dcmplx KKceex::iProd2(int Cp, int Hp, KKpart &p, int Cq, int Hq, KKpart &q){
/////////////////////////////////////////////////////////////////////////////////////
//        The same but C and helicities  as arguments                              //
/////////////////////////////////////////////////////////////////////////////////////
 dcmplx  Prod;
 if( abs(Cp) == 1 && abs(Cq) == 1 && abs(Hp) == 1 && abs(Hq) == 1 ){
   // Helicity conservation and non-conservation separately
   if(        Cp*Hp ==  Cq*Hq  ) {
      Prod = dcmplx( Cp *p.M *XiProd(q,p) + Cq *q.M *XiProd(p,q), 0.0);
   } else {
      Prod =  iProd1( Cp*Hp, p, q); // Helicity for v-spinor is fliped
   }// if Hp
 } else {
    cout<<"##### KKceex::iProd2: Wrong p.C="<< p.C<<" q.C= "<<q.C;
    cout<<"#####   or wrong helicities  Hp="<<  Hp<<"  Hq= "<<Hq<< endl;exit(92);
 }//if p.C
 return Prod;
}//iProd2
*/


dcmplx KKceex::iProd2( int Cp, int Lamp, KKpart &p,   int Cq, int Lamq, KKpart &q){

double mp= abs(p.M) * Cp;
double mq= abs(q.M) * Cq;
int Lp = Lamp;
int Lq = Lamq;
dcmplx Prod;
// Helicity for v-spinor is fliped
if( Cp < 0) Lp = -Lp;
if( Cq < 0) Lq = -Lq;
// Helicity conservation and non-conservation separately
if(     Lp == -Lq ) {
   Prod = iProd1(Lp,p,q);
//   cout<<"[1 iProd1]"<<endl;
}else if( Lp ==  Lq ) {
   Prod = dcmplx( mp*XiProd(q,p) +mq*XiProd(p,q), 0.0);
//   cout<<"[2 XiProd]"<<endl;
} else {
   cout<< "##### GPS_iProd2: Wrong Lp,Lq="<< Lp<<"  "<<Lq<<endl;
}
return Prod;
}


dcmplx KKceex::iProd3(int Lamp, KKpart &p, double mp, int Lamq, KKpart &q, double mq){
/////////////////////////////////////////////////////////////////////////////////////
dcmplx  Prod;
int  Lp,Lq;
Lp = Lamp;
Lq = Lamq;
// Helicity for v-spinor is fliped
if( mp < 0.0) Lp = -Lp;
if( mq < 0.0) Lq = -Lq;
// Helicity conservation and non-conservation separately
if(     Lp == -Lq ) {
   Prod = iProd1(Lp,p,q);
}else if( Lp ==  Lq ) {
   Prod = dcmplx( mp*XiProd(q,p) +mq*XiProd(p,q), 0.0);
} else {
   cout<<"##### GPS_iProd2: Wrong Lp,Lq= "<< Lp<<"  "<<Lq<<endl;
}
return Prod;
}//iProd3(


dcmplx KKceex::Soft(int sigma, KKpart &ph, KKpart &p1, KKpart &p2){
/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
//   2-fermion Soft factor at amplitude level                                      //
//   sigma   = photon polarization (+1,-1)                                         //
//   ph      = photon  4-momentum                                                  //
//   p1,p2   = fermion 4-momenta                                                   //
//   for KeyArb =1 dependend of vector auxial-gauge vector b=beta=m_b !!!          //
//                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////
//
dcmplx Soft;
if(m_KeyArb == 0 ){
 dcmplx Sqr2 = dcmplx(sqrt(2.0), 0.0);
 double pk1 = p1*ph;
 double pk2 = p2*ph;
 dcmplx bf1 =  Sqr2* iProd1(sigma,ph,p1)*XiProd(p1,ph); //!!! =GPS_bfact(sigma,ph,p1);
 dcmplx bf2 =  Sqr2* iProd1(sigma,ph,p2)*XiProd(p2,ph); //!!! =GPS_bfact(sigma,ph,p2);
 Soft = -bf1/(2*pk1) +bf2/(2*pk2);
} else {
 double pk1 = p1*ph;
 double pk2 = p2*ph;
 dcmplx bf1 = Bfacb(sigma,ph,p1);
 dcmplx bf2 = Bfacb(sigma,ph,p2);
 Soft = -bf1/(2*pk1) +bf2/(2*pk2);
}//
return Soft;
}//KKceex::Soft


dcmplx KKceex::GPS_Sof1(int sigma, KKpart &ph, KKpart &pf){
/////////////////////////////////////////////////////////////////////////////////////
//   Single soft photon contribution to Soft factor at amplitude level             //
//   sigma   = photon polarization (+1,-1)                                         //
//   ph      = photon  4-momentum                                                  //
//   pf      = fermion 4-momenta                                                   //
/////////////////////////////////////////////////////////////////////////////////////
dcmplx Sof1;
if(m_KeyArb == 0) {
  Sof1 = sqrt(2.0)*iProd1(sigma,ph,pf)*XiProd(pf,ph) /(2.0*(pf*ph));
} else {
  Sof1 = Bfacb(sigma,ph,pf) /(2.0*(pf*ph));
}
return Sof1;
}//GPS_Sof1

dcmplx KKceex::GPS_Sof1x(int sigma, KKpart &ph, KKpart &pf){
/////////////////////////////////////////////////////////////////////////////////////
//   Single soft photon contribution to Soft factor at amplitude level             //
//   sigma   = photon polarization (+1,-1)                                         //
//   ph      = photon  4-momentum                                                  //
//   pf      = fermion 4-momenta                                                   //
/////////////////////////////////////////////////////////////////////////////////////
dcmplx Sof1x;
if(m_KeyArb == 0) {
  Sof1x = sqrt(2.0)*iProd1(sigma,ph,pf)*XiProd(pf,ph);
} else {
  Sof1x = Bfacb(sigma,ph,pf);
}
return Sof1x;
}//GPS_Sof1x

void KKceex::GPS_Make_eps(KKpart &ph, int Sig, dcmplx eps[4]){
//SUBROUTINE GPS_Make_eps(ph,Sig,eps)
/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
//  photon polarization 4-vector explicitelly calculated                           //
//                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////
/*
IMPLICIT NONE
INCLUDE 'GPS.h'
SAVE                      !!! <-- necessary !!!
*
DOUBLE PRECISION      k1(4),k2(4),k3(4),k4(4),ph(4)
DOUBLE PRECISION      k1m(4),k2m(4),k3m(4),k4m(4)
*
DOUBLE PRECISION      Fleps
*-----------------------------------------------------------------------------
INTEGER    j,k
INTEGER    Sig
DOUBLE COMPLEX  GPS_Sof1x,GPS_Sof1bx
DOUBLE COMPLEX  eps(4)
DOUBLE COMPLEX  x1,x2,x3,x1m,x2m,x3m
*/
//   k=4
//   k1(k)=sqrt(2D0)
//   k2(k)=1D0
//   k3(k)=1D0
//   k4(k)=1D0
//   k1m(k)=sqrt(2D0)
//   k2m(k)=1D0
//   k3m(k)=1D0
//   k4m(k)=1D0
//
//   k1(1)=1D0
//   k1(2)=1D0
//   k2(2)=1D0
//   k3(3)=1D0
//   k4(4)=1D0
//   k1m(1)=-1D0
//   k1m(2)= 1D0
//   k2m(2)=-1D0
//   k3m(3)=-1D0
//   k4m(4)=-1D0
KKpart k1  = KKpart(sqrt(2.0),  1.0,  0.0,  0.0);
KKpart k2  = KKpart(      1.0,  0.0,  1.0,  0.0);
KKpart k3  = KKpart(      1.0,  0.0,  0.0,  1.0);
KKpart k4  = KKpart(      1.0,  0.0,  0.0,  0.0);
KKpart k1m = KKpart(sqrt(2.0), -1.0,  0.0,  0.0);
KKpart k2m = KKpart(      1.0,  0.0, -1.0,  0.0);
KKpart k3m = KKpart(      1.0,  0.0,  0.0, -1.0);
KKpart k4m = KKpart(     -1.0,  0.0,  0.0,  0.0);
//   x1 =   GPS_Sof1x( 1,ph,k1 )
//   x1m=   GPS_Sof1x( 1,ph,k1m)
//   x2 =   GPS_Sof1x( 1,ph,k2 )
//   x2m=   GPS_Sof1x( 1,ph,k2m)
//   x3 =   GPS_Sof1x( 1,ph,k3 )
//   x3m=   GPS_Sof1x( 1,ph,k3m)
dcmplx x1 =   GPS_Sof1x( 1,ph,k1 );
dcmplx x1m=   GPS_Sof1x( 1,ph,k1m);
dcmplx x2 =   GPS_Sof1x( 1,ph,k2 );
dcmplx x2m=   GPS_Sof1x( 1,ph,k2m);
dcmplx x3 =   GPS_Sof1x( 1,ph,k3 );
dcmplx x3m=   GPS_Sof1x( 1,ph,k3m);
//   eps(4)=(x3+x3m)/2D0
//   eps(3)=-(x3-x3m)/2D0
//   eps(2)=-(x2-x2m)/2D0
//   eps(1)=-(x1-x1m)/2D0
//DO j=1,4
// eps(j)=eps(j)/2.0
//ENDDO
eps[0]= (x3+x3m)/2.0/2.0;
eps[3]=-(x3-x3m)/2.0/2.0;
eps[2]=-(x2-x2m)/2.0/2.0;
eps[1]=-(x1-x1m)/2.0/2.0;
//   IF (sig.LT.0d0) THEN
//     DO j=1,4
//       eps(j)=-DCONJG(eps(j))
//     ENDDO
//   ENDIF
if(Sig<0) for(int j=0;j<4;j++) eps[j]=-conj(eps[j]);
//   write(*,*) 'vec  ph=',ph
//   write(*,*) 'sig=',sig
//   write(*,*) 'eps(i)=',eps
//   write(*,*) 'consistency check; product=',
//$               ((eps(4)*ph(4)-eps(3)*ph(3)-eps(2)*ph(2)-eps(1)*ph(1)))
//   write(*,*) 'consistency check; square=',
//$               eps(4)*DCONJG(eps(4))-eps(3)*DCONJG(eps(3))-eps(2)*DCONJG(eps(2))-eps(1)*DCONJG(eps(1))
}//GPS_make_eps!!!



dcmplx KKceex::Bfacb(int sigma, KKpart &phot, KKpart &pferm){
/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
//   diagonal element of U-matrix for massive fermion = denominator in s-factor    //
//   dependend of vector auxial-gauge vector b=beta=m_b !!!                        //
//   sigma  = photon polarization (+1,-1)                                          //
//   phot   = photon 4-momentum                                                    //
//   pferm  = fermion 4-momentum                                                   //
//                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////
double mass = pferm.M;
dcmplx bfac = dcmplx(sqrt(2.0),0.0)/iProd1(-sigma, phot, m_b)
        *(     iProd1( -sigma, m_b, pferm) *iProd1( sigma, pferm, phot)
           +sqr(mass) *XiProd( m_b, pferm )  *XiProd( phot, pferm) );
return bfac;
}//KKceex::bfacb



void KKceex::PhelRandom(){
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//   Generate photon helicities randomly                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

m_nPhot = m_Event->m_nPhot;
//[[[[[[[[[[[[[
  float              rvec[101];
  pseumar_makevec_(rvec,m_nPhot);
//----------------
//double rvec[101];
//m_RNgen->RndmArray(m_nPhot, rvec);
//]]]]]]]]]]]]]
for(int i=1; i<=m_nPhot; i++){      //f77 indexing
   if( rvec[i-1] > 0.5 ) {
      m_Phel[i] =  0;
   } else {
      m_Phel[i] = +1;
   }
}//for
//[[[[[[[[[[[[[[
(*m_Out) << "%%%%%%%%% m_Phel= ";
for(int i=1; i<=m_nPhot;i++ ) (*m_Out)<<" "<<m_Phel[i]; (*m_Out)<<endl;
//]]]]]]]]]]]]]]
//if(m_icont <200){for(int i=1; i<=m_nPhot;i++ ) cout<<" "<<m_Phel[i]; cout<<endl;}
}//PhelRandom

void KKceex::TralorPrepare(int KTO){
// Prepare Lorentz transformations from tau rest frame with
// GPS quantization axes to LAB/CMS frame
  if( KTO<1 || KTO>2) { cout<<"Kceex::TralorPrepare: Wrong KTO"<<endl; exit(90);}
  //
  TVector3 b1;
  if( KTO == 1) b1 = (m_Event->m_Qf1).BoostVector();
  if( KTO == 2) b1 = (m_Event->m_Qf2).BoostVector();
//
  TLorentzRotation TraLor= TLorentzRotation(); // set to unity
  TLorentzRotation TraLin= TLorentzRotation(); // set to unity
//
  TraLor.Boost(b1);
  TraLin.Boost(-b1);
//
  TLorentzVector Xi(  m_Xi.P[1], m_Xi.P[2], m_Xi.P[3], m_Xi.P[0]);
  TLorentzVector Eta(m_Eta.P[1],m_Eta.P[2],m_Eta.P[3],m_Eta.P[0]);
//
  Xi.Transform( TraLin); // from LAB to tau rest frame
  Eta.Transform(TraLin); // from LAB to tau rest frame
//
  TVector3 xi  = Xi.Vect();  // extract 3-vector
  TVector3 eta = Eta.Vect(); // extract 3-vector
//
// GPS rule 1: z-axis antiparalel to xi
// GPS rule 2: x-axis in plane (+eta, -xi),
  TRotation Rot3;
  Rot3.SetZAxis(-xi,eta);

  TLorentzRotation Rot = TLorentzRotation(Rot3);
  TLorentzRotation RotInv = Rot.Inverse();

  TraLor = TraLor*Rot;
  TraLin = RotInv*TraLin;

  if( KTO == 1) {m_TraLor1 = TraLor; m_TraLin1 = TraLin;}
  if( KTO == 2) {m_TraLor2 = TraLor; m_TraLin2 = TraLin;}

}//TralorPrepare

void KKceex::TralorDoIt(int KTO, double P[], double Q[]){
// Execute Lorentz transformations from tau rest frame with
// GPS quantization axes down to LAB/CMS frame
TLorentzVector p,q;
//
for(int i=0; i<=3;i++) p[i]=P[i];
//
//cout<<"KKceex::TralorDoIt: KTO="<<KTO<<endl;
if( KTO ==1){
   q = m_TraLor1 * p; // transform
//   cout<<"p[xyzt]= "; cout<<" "<<p.X()<<" "<<p.Y()<<" "<<p.Z()<<" "<<p.T()<<endl;
//   q=p; q.Boost(b1);
 //  cout<<"q[o123]= "; cout<<" "<<q.X()<<" "<<q.Y()<<" "<<q.Z()<<" "<<q.T()<<endl;
} else {
   q = m_TraLor2 * p; // transform
//   q=p; q.Boost(b2);
}
//
for(int i=0; i<=3;i++) Q[i]=q[i];
//

}// TralorDoit



double  KKceex::SForFac(double alfpic, KKpart &p1, KKpart &p2, double Emin, double MasPhot){
/////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                         //
//   s-channel YFS formfactor for  acollinear fermion pair.                                //
//   Mass effects are eaxct.                                                               //
//                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////
double m1 = p1.M;
double m2 = p2.M;
double p1p2  =  p1*p2;
double Breal = m_BVR->Btilda( alfpic, p1p2, p1[0],p2[0], m1, m2,  Emin, MasPhot); // !! Exact
double Bvirt = m_BVR->SBvirt( alfpic, p1p2, m1, m2, MasPhot);                     // !! Exact
double SForFac = exp( Breal + Bvirt);
return SForFac;
}// SForFac !!!


double KKceex::TForFac(double alfpic, KKpart &p1, KKpart &p2, double Emin, double MasPhot){
/////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                         //
//   T-channel YFS formfactor for acollinear fermion pair.                                 //
//   m1 is assumed to be small, m2 can be finite.                                          //
//                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////
double m1   =  p1.M;
double m2   =  p2.M;
double p1p2 =  p1*p2;
double Breal = m_BVR->Btilda( alfpic, p1p2, p1[0],p2[0], m1, m2,  Emin, MasPhot);
double Bvirt = m_BVR->TBvirt( alfpic, p1p2, m1, m2, MasPhot);
double TForFac = exp( Breal + Bvirt);
//[[[[[[[[[[[[[[[[[
//if(m_icont<1000) cout<<"KKceex::TForFac:  Breal="<<Breal<<"  Bvirt="<<Bvirt<<endl;
//if(m_icont<1000) cout<<"alfpic="<<alfpic<<"  p1p2="<<p1p2<<" p1[0]="<<p1[0]<<" p2[0]="<<p2[0]<<endl;
//if(m_icont<1000) cout<<" m1="<<m1<<" m2="<<m2<<" Emin="<<Emin<<" MasPhot="<<MasPhot<<endl;
//]]]]]]]]]]]]]]]]]
return TForFac;
}// TForFac !!!


void  KKceex::MakeF1ini(double Svar, double Mas1, double Mas2, double Q, dcmplx &F1_1, dcmplx &F1_2){
/////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                         //
//   Initial state vertex correction, IR subtracted, in small mass approximation.          //
//   Second order LL+NLL is included.                                                      //
//                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////
//double Alfpi = 1.0/(DB->Alfinv0)/M_PI;
dcmplx cL = CnuA(Svar,Mas1,Mas2) - dcmplx(1.0);    // <-- this is just ln(s/m**2)-i*pi -1
F1_1 = (m_Alfpi*Q*Q)   *0.5*cL;
F1_2 = F1_1
     +sqr(m_Alfpi*Q*Q) *(
              +cL*cL/8.0
              +cL*( 3.0/32 -3.0/4*m_zeta2 +3.0/2*m_zeta3 ));
}//MakeF1ini


void  KKceex::MakeF1fin(double Svar, double Mas1, double Mas2, double Q, dcmplx &F1_1, dcmplx &F1_2){
/////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                         //
//   Final state vertex correction, second order                                           //
//   finite mass  only for first order part                                                //
//                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////
//double Alfpi = 1.0/(DB->Alfinv0)/M_PI;
dcmplx cL = CnuA(Svar,Mas1,Mas2) - dcmplx(1.0);    // <-- this is just ln(s/m**2)-i*pi -1
F1_1 = (m_Alfpi*Q*Q)   *0.5*cL;
F1_2 = F1_1
     +sqr(m_Alfpi*Q*Q) *(
              +cL*cL/8.0
              +cL*( 3.0/32 -3.0/4*m_zeta2 +3.0/2*m_zeta3 ) );
}//MakeF1fin

dcmplx KKceex::CnuA(double Svar, double Mas1, double Mas2){
/////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                         //
//   Function nu*A(nu) Complex version appropriate for s and t-chanels.                    //
//   No small mass approximation.                                                          //
//                                                                                         //
//       s-chanel:  Nu = (-s+m1**2+m2**2)/2 = -p1p2 < 0   p1 and p2 incoming or outgoing   //
//       t-chanel:  Nu = (-t+m1**2+m2**2)/2 =  p1p2 > 0   p1 incoming p2 outgoing          //
//                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////

if( abs(Mas1*Mas2) < 1e-10 ) {
             cout<<" ++++++ CnuA: STOP due to zero mass"<<endl; exit(90);}
dcmplx Mas12 = dcmplx(Mas1*Mas2);
dcmplx Nu    = dcmplx( (-Svar +sqr(Mas1)+ sqr(Mas2))/2.0 );
dcmplx xlam  = sqrt( ( Nu - Mas12)*( Nu + Mas12) );
// take care of numerical stability
dcmplx z;
if( Svar > 0.0) {
   z = Mas12/(Nu-xlam);
} else {
   z = (Nu + xlam)/Mas12;
}
//dcmplx CnuA  = Nu/xlam *log( z ); // wrong sign of i*PI

//Eps = DCMPLX( +1d0,0.d0)   ! plus was adjusted empiricaly
//CnuA  = Nu/xlam *BVR_CDLN( z ,Eps)

dcmplx Eps= dcmplx( +1.0,0.0);   // plus was adjusted empiricaly
dcmplx CnuA  = Nu/xlam *m_BVR->CDLN( z ,Eps);

return CnuA;
}//


void KKceex::MakeVini(KKpart p1, KKpart p2, KKpart ph, dcmplx &V1, dcmplx &V2){
/////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                         //
//   Initial state virtual correction to single bremsstrahlung                             //
//   IR subtracted, in small mass approximation.                                           //
//   Second order LL  +NLL is included                                                     //
//   See Phys.Rev. D65 (2002) 073030  eq. 2.18                                             //
//                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////
//
  int  KFini = m_Event->m_KFini;
  double Qe  = DB->Qf[ KFini];
  double m1= p1.M;
  double m2= p2.M;
  double p1p2  =  p1*p2;
  double Svar  = 2*(p1*p2) +m1*m1 +m2*m2;
  double r1 = (p1*ph)/p1p2;
  double r2 = (p2*ph)/p1p2;
  dcmplx cL = CnuA(Svar,m1,m2)-1.0;           // <-- this is just ln(s/m**2)-i*pi -1
  V1 = (m_Alfpi*Qe*Qe) *0.5*cL;                   // constant LL part
  V2 = (m_Alfpi*Qe*Qe) *0.5*(
          +log(r1)*log(1-r2)  +log(r2)*log(1-r1)           // LL part
          +m_BVR->Dilog(r1)          +m_BVR->Dilog(r2)             // NLL this and all the rest
          -1.0/2*sqr(log(1-r1))-1.0/2*sqr(log(1-r2))
          +3.0/2*log(1-r1)     +3.0/2*log(1-r2)
          +1.0/2*r1*(1-r1)/(1+sqr(1-r1))
          +1.0/2*r2*(1-r2)/(1+sqr(1-r2)));
}//MakeVini

void KKceex::MakeVfin(KKpart p3, KKpart p4, KKpart ph, dcmplx &V1, dcmplx &V2){
/////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                         //
//   Final state virtual correction to single bremsstrahlung                               //
//   IR subtracted, in small mass approximation.                                           //
//   Second order LL  (+NLL???) is included                                                //
//                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////
//
  int  KFfin = m_Event->m_KFfin;
  double Qf  = DB->Qf[ KFfin];
  double m3 = p3.M;
  double m4 = p4.M;
  double p3p4  =  p3*p4;
  double Svar  = 2*p3p4 +m3*m3+m4*m4;
  double r3 = (p3*ph)/p3p4;
  double r4 = (p4*ph)/p3p4;
// normal definition as in O(alf1) single-photon case
  double s3 = r3/( 1.0 +r3 +r4 );
  double s4 = r4/( 1.0 +r3 +r4 );
  dcmplx cL = CnuA(Svar,m3,m4)-1.0;           // <-- this is just ln(s/m**2)-i*pi -1
  V1 = (m_Alfpi*Qf*Qf) *0.5*cL;                   // constant LL part
// FSR formula of Eq (3.31) of Phys. Rev. D65, p.073030
  V2 = (m_Alfpi*Qf*Qf) *0.5*(
          -log(s3)*log(1-s4)      -log(s4)*log(1-s3)   // LL part
          -3.0/2*log(1-s3)        -3.0/2*log(1-s4)     // NLL this and all the rest
          -m_BVR->Dilog(s3)              -m_BVR->Dilog(s4)
          -1.0/2*s3/(1+sqr(1-s3)) -1.0/2*s4/(1+sqr(1-s4)) );
// ***  V2= (Alfpi*Q**2) *( -0.25d0*cL*DLOG((1d0-s3)*(1d0-s4)) )  !! LL formula averaged over s3,s4
// *** $   +(Alfpi*Q**2) *( +0.50d0*cL*DLOG((1d0-s3)*(1d0-s4)) )  !! corr. due YFS formfactor
}// MakeVfin

