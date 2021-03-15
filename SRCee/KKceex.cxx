///////////////////////////////////////////////////////////////////////////////
#include "KKceex.h"

ClassImp(KKceex);
ClassImp(KKcmplx4);
ClassImp(KKcmplx2);

#define SW20 setw(20)<<setprecision(14)


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

  // axial vectors (arbitrary lightlike ) for photon polarization

  m_HasFSR = DB->KeyFSR;     // temporary arangement
  int    KFini   = m_Event->m_KFini;
  int    KFfin   = m_Event->m_KFfin;
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
    PhelRandom();

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
  int Hel1,Hel2;
  KKpart ph1,ph2;
  dcmplx Sactu1,Sactu2;
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
          HiniPlus(KFini, KFfin, PX, ph1, Hel1, Sactu1, sProd);
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
              GPS_HiiPlus(Cfact2, KFini, KFfin, PX, ph1, Hel1, ph2, Hel2);
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
  } else {
     m_WtSet[ 1]  =   m_RhoExp0 /m_RhoCrud;    //!!! Interference ON
     m_WtSet[ 2]  =   m_RhoExp1 /m_RhoCrud;
     m_WtSet[ 3]  =   m_RhoExp2 /m_RhoCrud;
     m_WtBest     =   m_WtSet[ 3];
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
        for(int j3 = 0; j3<=1; j3++)
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
           }//j3,j4
  }// j1,j2

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
    }//for
  }//if

/////////////////////////////////////////////////////
// Dresed Born 16 amplitudes
  for(int j1 = 0; j1<=1; j1++)
    for(int j2 = 0; j2<=1; j2++)
      for(int j3 = 0; j3<=1; j3++)
        for(int j4 = 0; j4<=1; j4++){
            AmpBorn.m_A[j1][j2][j3][j4] =
                   m_SpinoTT[j1][j2][j3][j4]* m_FFacTT[j1]
                  +m_SpinoUU[j1][j2][j3][j4]* m_FFacUU[j1];
        }// for j1...j4

}//KKceex::Born

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

dcmplx BoxGG, BoxGZ, AmpBoxy, AmpBorn;
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
        AmpBorn = m_BornC.m_A[j1][j2][j3][j4];
        m_AmpExpo0.m_A[j1][j2][j3][j4] +=  Cfac*AmpBorn;         // O(alf^0)_exp
  //
        m_AmpExpo1.m_A[j1][j2][j3][j4] +=
                 +Cfac*AmpBorn*(1.0 +m_F1ini1)*(1.0 +m_F1fin1 )  // Born, O(alf1) m_FFactors
                 +Cfac*AmpBoxy;                                                    // O(alf1) boxes
  //
        m_AmpExpo2.m_A[j1][j2][j3][j4] +=
                 +Cfac*AmpBorn*(1.0 +m_F1ini2)*(1.0 +m_F1fin2 )  // Born, O(alf2) FFactors
                 +Cfac*AmpBoxy ;                                                   // O(alf1) boxes
     }//for j3,j4

}//for j1,j2
//---------------------
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
  AmpAdd( m_AmpExpo2, sProd*Vir2,m_BornC);
}// GPS_HiniPlus



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




void KKceex::EWFFact(int KFini, int KFfin, double Svar, dcmplx &Ve, dcmplx &Vf, dcmplx &Ae, dcmplx &Af ){
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


void KKceex::EWFFact(int KFini, int KFfin, double Svar, double CosThetD,
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
  if( CosThetD == 0.0) VVcor  = dcmplx(1.0); //????
}//


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
     U.m_A[0][1] =  0.0;                                                      // (+-)
     U.m_A[1][0] =  Sqr2*( -m1*XiProd(p2,p1) +m2*XiProd(p1,p2));      // (-+)
  }else  if(sigma == -1 ) {
     U.m_A[0][0] =  Sqr2*( XiProd(p1,ph) *iProd1(-1,ph,p2));          // (++)
     U.m_A[1][1] =  Sqr2*( XiProd(p2,ph) *iProd1(-1,ph,p1));          // (--)
     U.m_A[1][0] =  0.0;                                                      // (-+)
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
//                                                                                 //
//   Single soft photon contribution to Soft factor at amplitude level             //
//   sigma   = photon polarization (+1,-1)                                         //
//   ph      = photon  4-momentum                                                  //
//   pf      = fermion 4-momenta                                                   //
//                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////
dcmplx Sof1;
if(m_KeyArb == 0) {
  Sof1 = sqrt(2.0)*iProd1(sigma,ph,pf)*XiProd(pf,ph) /(2.0*(pf*ph));
} else {
  Sof1 = Bfacb(sigma,ph,pf) /(2.0*(pf*ph));
}
return Sof1;
}//GPS_Sof1



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

double rvec[101];
m_RNgen->RndmArray(m_nPhot, rvec);

for(int i=1; i<=m_nPhot; i++){      //f77 indexing
   if( rvec[i-1] > 0.5 ) {
      m_Phel[i] =  0;
   } else {
      m_Phel[i] = +1;
   }
}//for
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

