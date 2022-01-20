///////////////////////////////////////////////////////////////////////////////
#include "KKborn.h"

ClassImp(KKborn);


KKborn::KKborn()
{
  // This constructor is for ROOT streamers ONLY
  // all pointers has to be NULLed
  cout<< "----> KKborn Default Constructor (for ROOT only) "<<endl;
  m_Out= NULL;
  DB = NULL;
  DZ = NULL;
}

///_____________________________________________________________
KKborn::KKborn(ofstream *OutFile)
{
  cout<< "----> KKborn USER Constructor "<<endl;
  m_Out = OutFile;
  DB = NULL;
  DZ = NULL;
  m_icont = 0;
}//KKborn

///______________________________________________________________________________________
KKborn::~KKborn()
{
  //Explicit destructor
  cout<< "----> KKborn::KKborn !!!! DESTRUCTOR !!!! "<<endl;
}///destructor

double KKborn::sqr( const Double_t x ){ return x*x;};

///______________________________________________________________________________________
void KKborn::Initialize()
{
  cout  << "----> KKborn::Initialize, Entering "<<endl;
//=================================================================
// BX*** macros are in MCdev/BXFORMAT.h
  BXOPE(*m_Out);
  BXTXT(*m_Out,"========================================");
  BXTXT(*m_Out,"======    KKborn::Initialize      ======");
  BXTXT(*m_Out,"========================================");
  m_CMSene = DB->CMSene;
  m_KeyDBG = 0;
  *m_Out<< " Testing acces to DBase "<<endl;
  *m_Out<< " KKborn::Initialize: m_CMSene = "<<m_CMSene<<endl;

  ///////////////////////////////////////////////////
}// Initialize



double KKborn::BornSimple(int KFi, int KFf, double svar, double costhe){
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// This routine provides unsophisticated Born differential cross section     //
// at the crude x-section level, with Z and gamma s-chanel exchange.         //
// KKMC-hh: Changed to use swsq from input parameter, not Dizet, if KeyElw<1.//
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
double ss = svar;
// Z and gamma couplings to beams (electrons)
double MZ    = DB->MZ;
double MW    = DB->MW;
//////////
double swsq  = DZ->D_swsq; // from Dizet!!!
double gammz = DZ->D_GamZ; // from Dizet!!!
//////////
double deno= 4*sqrt(swsq*(1-swsq));
// initial state
double T3e = DB->T3f[KFi];  // isospin, L-hand component
double Qe  = DB->Qf[ KFi];  // electric charge
double Ve  = (2*T3e -4*Qe*swsq)/deno;
double Ae  =  2*T3e            /deno;
// final state
int    NCf = DB->Nc[ KFf];   // number of colours
double T3f = DB->T3f[KFf];   // isospin, L-hand component
double Qf  = DB->Qf[ KFf];   // electric charge
double Vf    = (2*T3f -4*Qf*swsq)/deno;
double Af    =  2*T3f            /deno;
// Switch off Z or gamma
if(DB->KeyZet == 0) { Ve=0; Ae=0; }
if(DB->KeyZet == 9) { Qe=0; Qf=0; }
double BWD;
if ( (DB->KeyZet == -2) || (DB->KeyZet == -1) )
   BWD = sqr(ss-sqr(MZ)) + sqr(gammz*MZ);    // <--! fixed width
else
   BWD = sqr(ss-sqr(MZ)) + sqr(gammz*ss/MZ); // <--! running width
//
double chi2,rechi;
if (DB->KeyZet == -1) {   // <--! fixed width with redefined parameters
   BWD = BWD * (1 + sqr(gammz/MZ));
   chi2 = sqr(ss)/BWD;
   rechi=(ss-sqr(MZ)- sqr(gammz))*ss/BWD;
} else {                  // <--! running width or fixed without redefinition
   chi2 = sqr(ss)/BWD;
   rechi=(ss-sqr(MZ))*ss/BWD;
}
double xe= Ve*Ve +Ae*Ae;
double xf= Vf*Vf +Af*Af;
double ye= 2*Ve*Ae;
double yf= 2*Vf*Af;
double ff0= Qe*Qe *Qf*Qf +2*rechi*Qe*Qf*Ve*Vf +chi2*xe*xf;
double ff1=              +2*rechi*Qe*Qf*Ae*Af +chi2*ye*yf;
double Born    = (1+ sqr(costhe) )*ff0 +2*costhe*ff1;

////////////////////////////////////////////////////////
// electron neutrino is special,
// eq. (9)http://arxiv.org/abs/hep-ph/0110371v1
double ve, tt, ChiW, sig_s, sig_st, sig_t, Born2;
dcmplx ChiZ;
if( KFf == 12){
  ve = 1.0 -4.0*abs(Qe)*swsq;
  tt = -ss*(1.0 - costhe )/2.0;
  ChiZ= sqr(MZ)/dcmplx(-ss+sqr(MZ), -ss/MZ*gammz);
  ChiW =       MW*MW/(-tt+MW*MW);
  sig_s  = real(ChiZ*conj(ChiZ)) *( (1.0+sqr(costhe))*(1.0+sqr(ve)) +4.0*costhe* ve  );
  sig_st = -4.0/3.0*real(ChiZ)*ChiW*sqr(1.0+costhe)*(1+ve);
  sig_t  = 8.0/3.0 *sqr(ChiW)*sqr(1.0+costhe);
  Born    = 2.0*sqr(ss/sqr(MZ))/sqr(deno)/sqr(deno) *( sig_s +sig_st +sig_t);
}//Kff
//     Colour factor
Born = NCf*Born;
if( abs(costhe) > 1) cout<< "------------> BornV: costhe="<<costhe<<endl;
// This is a bit crude method of introducing threshold behaviour
// cos(theta) depencence incorrect!!!
double amfin = DB->fmass[KFf];     // mass
double thresh;
if(    svar <=  4*sqr(amfin) )
   thresh=0;
else if(svar <= 160*sqr(amfin) ) {
   double amx2=4*sqr(amfin)/svar;
   thresh=sqrt(1 -amx2 )*(1+amx2/2);
} else {
   thresh=1;
}
Born= Born*thresh;
return Born;
}//KKbornSimple

//_________________________________________________________________________
double KKborn::Born_DizetS(int KFi, int KFf, double svar, double CosThe){

  DZ->InterpoGSW(KFi, KFf, svar, CosThe);
  double Born = Born_Dizet(KFi, KFf, svar, CosThe, 0.0, 0.0, 0.0,0.0);
  return Born;

}//Born_DizetS

//_________________________________________________________________________
double KKborn::Born_Dizet(int KFi, int KFf, double svar, double CosThe,
                          double eps1, double eps2, double ta, double tb)
{
//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//   Calculates differential born cross section.                            //
//   For Mode=0 pure Born and for Mode=1 electroweak corrs. are added.      //
//   KFi,KFf can be also negative for antiparticle, in this case it is      //
//   important to produce tables with correct input KFini, KFfin !!!        //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
  m_icont++;
  double xupgi[3],    // Left/Right coupling gamma initial
         xupzi[3],    // Left/Right coupling Z     initial
         xupgf[3],    // Left/Right coupling gamma final
         xupzf[3];    // Left/Right coupling Z     final
//-------------------------------------------------------------------------------------
  double  t3e,     // Left izospin initial ???
          t3f;     // Left izospin final   ???
  int kolor=1;     // Color final fermion, here =1 for leptons
//-------------------------------------------------------------------------------------
// All dimensions increased by one in order to keep f77 indexing
  dcmplx  aborn[3][3], aphot[3][3], azett[3][3];
  dcmplx  xupzfp[3],xupzip[3];
  dcmplx  abornm[3][3],aphotm[3][3],azettm[3][3];
  dcmplx  propa,propz;
  dcmplx  xupf,xupi,xff[5],xfem,xfota,xrho,xke,xkf,xkef;
  dcmplx  xthing,xve,xvf,xvef;
////////////////////////////////////////////////////////////////////////
//               Coupling constants                                   //
////////////////////////////////////////////////////////////////////////
// Z and gamma couplings
  double MZ    = DB->MZ;   // from input (xpar)
  double swsq  = DZ->D_swsq; // from Dizet!!!
  double gammz = DZ->D_GamZ; // from Dizet!!!
// initial state
  double amin = DB->fmass[KFi]; // not used
  double T3e  = DB->T3f[KFi];   // isospin, L-hand component
  double qe   = DB->Qf[ KFi];   // electric charge
  double aizor =0;   // ok
  double aizol =T3e; // ok
  xupgi[1]=qe;
  xupgi[2]=qe;
////////////////////
//   t3e    = aizol+aizor // not available
  xupzi[1]=(aizor-qe*swsq)/sqrt(swsq*(1-swsq));
  xupzi[2]=(aizol-qe*swsq)/sqrt(swsq*(1-swsq));
// final state
  double amfin = DB->fmass[KFf];
  double T3f  = DB->T3f[KFf];  // isospin, L-hand component
  double qf   = DB->Qf[ KFf];  // electric charge
  aizor=0;
  aizol=T3f;
  xupgf[1]=qf;
  xupgf[2]=qf;
//  t3f    =  aizol+aizor // not available
  xupzf[1]=(aizor -qf*swsq)/sqrt(swsq*(1-swsq));
  xupzf[2]=(aizol -qf*swsq)/sqrt(swsq*(1-swsq));
//
  double sinthe = sqrt(1-sqr(CosThe));
  double beta   = sqrt( std::max(0.0 , 1-4*amfin*amfin/svar) );
// Multiply axial coupling by beta factor.
  xupzfp[1]= 0.5*(xupzf[1]+xupzf[2])+0.5*beta*(xupzf[1]-xupzf[2]);
  xupzfp[2]= 0.5*(xupzf[1]+xupzf[2])-0.5*beta*(xupzf[1]-xupzf[2]);
  xupzip[1]= 0.5*(xupzi[1]+xupzi[2])     +0.5*(xupzi[1]-xupzi[2]);
  xupzip[2]= 0.5*(xupzi[1]+xupzi[2])     -0.5*(xupzi[1]-xupzi[2]);
// Final state vector coupling
  xupf     = 0.5*(xupzf[1]+xupzf[2]);
  xupi     = 0.5*(xupzi[1]+xupzi[2]);
  xthing   = 0;
////////////////////////////////////////////////////////////////////////
//                          Propagators                               //
////////////////////////////////////////////////////////////////////////
// Multiply axial coupling by beta factor.
// Add formfactors initialisation of s-dependent electro-weak form factors and
// photonic vacuum polarisation
// (electro-weak box contributions left out here, they depend on acos)
//  CALL BornV_GetQCDcor2(KFf,RSQV,RSQA)
  double RSQV=1, RSQA=1;
//  bornv_getqcdcor2_(KFf,RSQV,RSQA);  //?????
  xff[1]=DZ->m_GSW[1-1];
  xff[2]=DZ->m_GSW[2-1];
  xff[3]=DZ->m_GSW[3-1];
  xff[4]=DZ->m_GSW[4-1];
//xffa  =UNDEFINED !!!!
  xfem  =DZ->m_GSW[6-1];
  xfota =DZ->m_GSW[7-1];
//-------------------------------------------------
  xrho =xff[1];
  xke  =xff[2];
  xkf  =xff[3];
  xkef =xff[4];
  double qfm  = abs(qf);
  double qem  = abs(qe);
  double xe   =  1 -4*swsq*qem;
  double xf   =  1 -4*swsq*qfm;
  double xef  = -1 +xe +xf +16*qem*qfm*swsq*swsq; // xef=xe*xf !!!
  xve  =  1.0 -4*swsq*qem*xke;
  xvf  =  1.0 -4*swsq*qfm*xkf;
  xvef = -1.0 +xve +xvf +16*qem*qfm*swsq*swsq*xkef;
// Multiply axial  coupling by beta factor.
// Multiply vector coupling by form-factor.
// Multiply final vector by RSQV and final axial by RSQA (QCD corrections)
  xupgf[1]=xupgf[1]*RSQV;
  xupgf[2]=xupgf[2]*RSQV;
  xupzfp[1]=0.5*(xupzf[1]+xupzf[2])*xvf/xf*RSQV  +0.5*(xupzf[1]-xupzf[2])*beta*RSQA;
  xupzfp[2]=0.5*(xupzf[1]+xupzf[2])*xvf/xf*RSQV  -0.5*(xupzf[1]-xupzf[2])*beta*RSQA;
  xupzip[1]=0.5*(xupzi[1]+xupzi[2])*xve/xe  +0.5*(xupzi[1]-xupzi[2]);
  xupzip[2]=0.5*(xupzi[1]+xupzi[2])*xve/xe  -0.5*(xupzi[1]-xupzi[2]);
// Final state vector coupling
  xupf     =0.5*(xupzf[1]+xupzf[2])*xvf/xf*RSQV;
// Double vector formfactor thing
  xthing=0.25*(xupzf[1]+xupzf[2])*(xupzi[1]+xupzi[2])*(xvef/xef-xvf*xve/xe/xf)*RSQV;
  propa =1.0/svar/(2.0-xfem);
  if ((DB->KeyZet == -1) ||(DB->KeyZet == -2)) {
      // KKMC-hh: implement fixed width option...
      // but note: fixed width isn't really compatible with Dizet.
      propz =1.0/dcmplx(svar-MZ*MZ, MZ*gammz);
  } else {
      propz =1.0/dcmplx(svar-MZ*MZ, svar/MZ*gammz);
         }
      if (DB->KeyZet == -1) { // fixed width parameter redefinition
      propz = propz/dcmplx(1.0,gammz/MZ);
  }
// Replace Born normalization of Z propagator by the better one
  double del1 = DB->GFermi *MZ*MZ *DB->Alfinv0 /(sqrt(2)*8 *M_PI); // AlfInvZ ???
  double del0 =1/(swsq*(1-swsq))/16;
  propz = propz*del1/del0*xrho;
////////////////////////////////////////////////////////////////////////
//             Spin amplitudes   Z+gamma case                         //
////////////////////////////////////////////////////////////////////////
  for(int  i=1; i<=2;i++){
    for(int j=1; j<=2; j++){
      double regula= (3-2*i)*(3-2*j) + CosThe;
      double regulm=-(3-2*i)*(3-2*j) * sinthe *2*amfin/sqrt(svar);
      aphot[i][j]=propa*(xupgi[i] *xupgf[j]*regula);
      azett[i][j]=propz*(xupzip[i]*xupzfp[j]+xthing)*regula;
      aborn[i][j]=aphot[i][j]+azett[i][j];
      aphotm[i][j]= propa*dcmplx(0,1)  *xupgi[i]*xupgf[j]    *regulm;
      azettm[i][j]= propz*dcmplx(0,1)*(xupzip[i]*xupf+xthing)*regulm;
      abornm[i][j]=aphotm[i][j]+azettm[i][j];
    }// for j
  }// for i
////////////////////////////////////////////////////////////////////////
//           Differential X-section out of spin amplituds             //
//  Helicity conservation explicitly obeyed:                          //
//  Only diagonal elements of the spin density matrices.              //
//  (Only longitudinal polarizations)                                 //
////////////////////////////////////////////////////////////////////////
  double polar1 =  (eps1);
  double polar2 = (-eps2);
  double Born   =  0;
  for(int i=1; i<=2; i++){
    int helic= 3-2*i;
    for(int j=1; j<=2; j++){
      int helit=3-2*j;
      double factor = kolor*(1+helic*polar1)*(1-helic*polar2)/4;
      double factom = factor*(1+helit*ta)*(1-helit*tb);
      factor = factor*(1+helit*ta)*(1+helit*tb);
// Normal case (mass terms included in Born. is it better ??????)
      Born=Born +sqr( abs(aborn[i][j])  ) *factor;
      Born=Born +sqr( abs(abornm[i][j]) ) *factom;
    }// for j
  }// for i
// phase space threshold factor, and multiply by svar**2 to get R-units!
  if (svar > 4*amfin*amfin) {
    double thresh=sqrt(1-4*amfin*amfin/svar);
    Born = Born*svar*svar*thresh;
  } else {
    Born=0;
  }
//
   return Born;
}//Born_Dizet
