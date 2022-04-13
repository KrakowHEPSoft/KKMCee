///////////////////////////////////////////////////////////////////////////////
// Replaces f77 module TauPair
///////////////////////////////////////////////////////////////////////////////
#include "TauPair.h"

ClassImp(TauPair);

TauPair::TauPair()
{
  // This constructor is for ROOT streamers ONLY
  // all pointers has to be NULLed
  cout<< "----> TauPair Default Constructor (for ROOT only) "<<endl;
  m_Out  = NULL;
  DB     = NULL;
  m_Event    = NULL;
  m_GPS      = NULL;
  m_RNgen    = NULL;

}

///_____________________________________________________________
TauPair::TauPair(ofstream *OutFile)
{
  cout<< "----> TauPair USER Constructor "<<endl;
  m_Out  = OutFile;
  DB     = NULL;
  m_Event    = NULL;
  m_GPS      = NULL;
  m_RNgen    = NULL;
}//TauPair

///______________________________________________________________________________________
TauPair::~TauPair()
{
  //Explicit destructor
  cout<< "----> TauPair::TauPair !!!! DESTRUCTOR !!!! "<<endl;
}///destructor

double TauPair::sqr( const Double_t x ){ return x*x;};

///______________________________________________________________________________________
void TauPair::Initialize(double xpar[])
{
  cout  << "----> TauPair::Initialize, Entering "<<endl;
//=================================================================
// BX*** macros are in MCdev/BXFORMAT.h
  BXOPE(*m_Out);
  BXTXT(*m_Out,"========================================");
  BXTXT(*m_Out,"======    TauPair::Initialize     ======");
  BXTXT(*m_Out,"========================================");


  int ITAUXPAR=2000;

  m_IsInitialized = xpar[415-1];  // General mask for tau channel

// switches of tau+ tau- decay modes !!
  m_itdkRC        = xpar[ITAUXPAR+4-1];   // QED internal rad. in leptonic decays
  int Jak1        = xpar[ITAUXPAR+1-1];   // Decay Mask for first tau
  int Jak2        = xpar[ITAUXPAR+2-1];   // Decay Mask for second tau
  if( (Jak1 == -1) && (Jak2 == -1) ) m_IsInitialized = 0;

  m_KeyClone      = 1;       // dip-switch for cloning procedure, =1,2
  m_KeyClone      = 2;       // dip-switch for cloning procedure, =1,2
  BXTXT(*m_Out, " KK interface of Tauola                 ");
  BX1I( *m_Out, "  IsInit", m_IsInitialized, "xpar[415]       =");
  BX1I( *m_Out, "    Jak1",            Jak1, "xpar[2001]      =");
  BX1I( *m_Out, "    Jak2",            Jak2, "xpar[2002]      =");
  BX1I( *m_Out, "  itdkRC",        m_itdkRC, "xpar[2004]      =");
  BX1I( *m_Out, "KeyClone",      m_KeyClone, "Cloning proc.   =");
  BXCLO(*m_Out);

// Initialisation of tau decay package TAUOLA; ITAUXPAR is for indirect adressing.
  inietc_(&ITAUXPAR,xpar);
  if( m_IsInitialized == 0) {
     BXOPE(*m_Out);
     BXTXT(*m_Out, " !!!!! Tauola inhibited !!!!    ");
     BXCLO(*m_Out);
  } else {
// Initialisation of TAUOLA
    inimas_(&ITAUXPAR,xpar);
    initdk_(&ITAUXPAR,xpar);
    double xk0qed = 0.1;            // <=== It seems to be never used
    iniphy_(&xk0qed);
    int JAK =-1;
    dekay_(&JAK, m_HvecTau1);

//[[[[[[[[[[[[[[[[[[[[[[[[
// old fortran PHOTOS to be removed
// Initialization of old fortran PHOTOS
//    if(DB->KeyPhts >1) phoini_();
//]]]]]]]]]]]]]]]]]]]]]]]]

// Initialization of PHOTOS++
// KeyPhts =0 for off; =1 in non-leptonic; =2 in all decays
    if(DB->KeyPhts ==1 ){
      Photos::initialize();
    } else if( DB->KeyPhts ==2){
// Suppressing Photos for leptonic decays
      Photos::initialize();
      Photos::suppressAll();
      Photos::forceBremForBranch(0, 15);
      Photos::forceBremForBranch(0, -15);
      Photos::suppressBremForDecay(3, 15, 16, 11, -12);
      Photos::suppressBremForDecay(3, -15, -16, -11, 12);
      Photos::suppressBremForDecay(3, 15, 16, 13, -14);
      Photos::suppressBremForDecay(3, -15, -16, -13, 14);
    }//KeyPhts
  }//IsInitialized
///////////////////////////////////////////////////
}// end if Initialize

///______________________________________________________________________________________
void TauPair::Make1(){
  int J;
  if( m_IsInitialized != 0) {
    J=1; dekay_(&J,m_HvecTau1); // TAUOLA
    J=2; dekay_(&J,m_HvecTau2); // TAUOLA
  }
}//Make1

///______________________________________________________________________________________
void TauPair::Clone(){
/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
//   This routine is strongly interrelated with Tralor  !!!                        //
//                                                                                 //
//   Cloning tau decays by additional rotation tau decay products with respect     //
//   to frames  initially used in the decay simulation.                            //
/////////////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////////////////////////
//   Generation of random two independent Euler rotations                          //
/////////////////////////////////////////////////////////////////////////////////////
  double rrr[10];
  m_RNgen->RndmArray(3, rrr);
  m_alfa1  = 2.0*M_PI*rrr[2];        // azimuthal angle in (0,2*pi)
  m_beta1  = acos(2.0*rrr[0]-1.0);   // polar angle     in (0,  pi)
  m_gamma1 = 2.0*M_PI*rrr[1];        // azimuthal angle in (0,2*pi)
//------------------------------------------------
  m_RNgen->RndmArray(3, rrr);
  m_alfa2  = 2.0*M_PI*rrr[2];        // azimuthal angle in (0,2*pi)
  m_beta2  = acos(2.0*rrr[0]-1.0);   // polar angle     in (0,  pi)
  m_gamma2 = 2.0*M_PI*rrr[1];        // azimuthal angle in (0,2*pi)
//------------------------------------------------
  m_H1.SetPxPyPzE(m_HvecTau1[0],m_HvecTau1[1],m_HvecTau1[2],m_HvecTau1[3]);
  m_H2.SetPxPyPzE(m_HvecTau2[0],m_HvecTau2[1],m_HvecTau2[2],m_HvecTau2[3]);
  if(m_KeyClone == 1) {
/////////////////////////////////////////////////////////////////////////////////////
//   Cloning tau decay with help of  Euler rotations FIRST method                  //
    double Habs1 = (m_H1.Vect()).Mag();
    double Habs2 = (m_H2.Vect()).Mag();
    // Standart phi, theta for polarimeter fectors, phi in (0,2*pi), theta in (0,pi)
    m_phi1  =0.0; m_thet1 =0.0;
    m_phi2  =0.0; m_thet2 =0.0;
    if(Habs1 > 1e-5) { m_phi1  = m_H1.Phi(); m_thet1 = m_H1.Theta(); }
    if(Habs2 > 1e-5) { m_phi2  = m_H2.Phi(); m_thet2 = m_H2.Theta(); }
    m_H1.SetPxPyPzE(0.0, 0.0, Habs1, 1.0);
    m_H2.SetPxPyPzE(0.0, 0.0, Habs2, 1.0);
    m_Event->RotEul(m_beta1, m_gamma1, &m_H1);
    m_Event->RotEul(m_beta2, m_gamma2, &m_H2);
    for(int i=0; i<4;i++) m_HvClone1[i]=m_H1[i];
    for(int i=0; i<4;i++) m_HvClone2[i]=m_H2[i];
    //
  } else if(m_KeyClone == 2) {
/////////////////////////////////////////////////////////////////////////////////////
//   Cloning tau decay with help of  Euler rotations, SECOND method                //
     m_Event->RotEuler(m_alfa1, m_beta1, m_gamma1, &m_H1);
     m_Event->RotEuler(m_alfa2, m_beta2, m_gamma2, &m_H2);
     for(int i=0; i<4;i++) m_HvClone1[i]=m_H1[i];
     for(int i=0; i<4;i++) m_HvClone2[i]=m_H2[i];
 } else {
     (*m_Out)<< " ##### STOP in Taupair_Clone: wrong KeyClone= "<< m_KeyClone<< endl;
     cout <<    " ##### STOP in Taupair_Clone: wrong KeyClone= "<< m_KeyClone<< endl;
     exit(9);
  }
}//Clone


///______________________________________________________________________________________
void TauPair::ImprintSpin(){
//////////////////////////////////////////////////
//     introduces spin effects by rejection     //
//////////////////////////////////////////////////
  int loop=0;
  double rn,wt,wt0,wt1,wt2, wtmax=4.0;
e1099:
  loop=loop+1;
  Clone();   // Cloning tau decay by Euler rotation
  m_GPS->MakeRho2(m_HvClone1,m_HvClone2,wt0,wt1,wt2);
  wt = wt1;                         // why not wt2???
  rn = m_RNgen->Rndm();
  if (wt < wtmax*rn  && loop<100) goto e1099;
}//ImprintSpin

///______________________________________________________________________________________
void TauPair::Make2(){
//  taupair_make2_();        // book-keeping, Photos, HepEvt
  int ih1,ih2;
  hepevt_getf_(   ih1);         // fermion is here
  hepevt_getfbar_(ih2);         // antifermion is here
  tauface_setfermpos_(ih1,ih2);  // set ffbar positions in Tauola

  int J;
  J=11;  dekay_(&J,m_HvClone1);
  J=12;  dekay_(&J,m_HvClone1);

/////////////////////////////////////////////////////////////////////////////////////
//                    Photos comes last                                            //
/////////////////////////////////////////////////////////////////////////////////////
//[[[[[[[[[[[[[[[[[[[[[[[[[[
// temporary printout before photos
  if( m_Event->m_EventCounter <=3){
    J=2; pyhepc_(J);       // HepEvt-->Pythia
    tauface_print_();
  }
//]]]]]]]]]]]]]]]]]]]]]]]]]]

//[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
// fortran photos to be removed
//  if(DB->KeyPhts >0 ) {
//    photos_(ih1);  // Photos works on HepEvt
//    photos_(ih2);
//  }//IFPHOT
//]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]

//[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
// PRINTOUT of /hepevt/ USIGNG Pythia
/////////////////////////////////////////////
  J=2; pyhepc_(J);       // HepEvt-->Pythia
  if( m_Event->m_EventCounter <=30){
  tauface_print_();
  }//EventCounter
//]]]]]]]]]]]]]]]]]]]]
}//Make2

/////////////////////////////////////////////////////////////////////////////////////
/// Run Photos
void TauPair::RunPhotosPP(){
// Flag for Photos c++ 
  if(DB->KeyPhts > 0) {
// Process HEPMC3 event by PHOTOS++
  int buf = -m_Hvent->particles().size();
// test print before photos
  cout <<   "!!!!!!!!!!!!!!!!!!!!!!!!!! Before Photos !!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
//  (*m_Out) <<   "!!!!!!!!!!!!!!!!!!!!!!!!!! Before Photos !!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;

  Print::listing(*m_Hvent);
  PhotosHepMC3Event photosEvent(m_Hvent);
  photosEvent.process();
// test print after photos
  cout<<        "!!!!!!!!!!!!!!!!!!!!!!!!!! After Photos !!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
//  (*m_Out) <<   "!!!!!!!!!!!!!!!!!!!!!!!!!! After Photos !!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;

  Print::listing(*m_Hvent);

  buf += m_Hvent->particles().size();
  if(buf>0 && m_Event->m_EventCounter <=200){
     cout<<     ">>>>>>> TauPair::Make2: ["<<m_Event->m_EventCounter<< "] PHOTOS++ added "<<buf<<" new photons !!!!!!"<<endl;
    (*m_Out) << ">>>>>>> TauPair::Make2: ["<<m_Event->m_EventCounter<< "] PHOTOS++ added "<<buf<<" new photons !!!!!!"<<endl;
    }
  }
} //end Run Photos

void TauPair::Finalize(){
////////////////////////////////
// Final printout from Tauola
////////////////////////////////
//  CALL DEKAY(100,HvecDummy)
//
    int JAK =100;
    dekay_(&JAK, m_HvecTau1);
//
}//Finalize

void TauPair::Tralo4(int Kto, float P[], float Q[], float &AM){
/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
//   This routine is strongly interrelated with Taupair::Clone!                    //
//                                                                                 //
//  SUBSITUTE OF TRALO4                                                            //
//  TRALO4 is called in TAUOLA --> /hepevt/ interface to boost from tau+-          //
//  restframe to lab. It includes rotations in tau rest frame due to spin effect   //
//  implementation                                                                 //
//                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////
// locals
double  Pd[4];
//* ------------------------------------------------------------
AM = sqrt(abs( P[3]*P[3] -P[2]*P[2] -P[1]*P[1] -P[0]*P[0] ));  // Mass
//
for(int k=0; k<4; k++) Pd[k]=P[k]; // from REAL to DOUBLE PRECISION
//
if(m_KeyClone == 1) {
   m_PP.SetPxPyPzE(Pd[0],Pd[1],Pd[2],Pd[3]);
   if(   Kto == 1) {
      m_Event->RotEulInv( m_thet1, m_phi1,   &m_PP);
      m_Event->RotEul(    m_beta1, m_gamma1, &m_PP);
   }else if(Kto == 2) {
      m_Event->RotEulInv( m_thet2, m_phi2,   &m_PP);
      m_Event->RotEul(    m_beta2, m_gamma2, &m_PP);
   } else {
     (*m_Out)<<"###### STOP in TRALO4: Wrong Kto = "<<Kto<<endl;
     cout    <<"###### STOP in TRALO4: Wrong Kto = "<<Kto<<endl; exit(9);
   }
   for(int i=0; i<4;i++) Pd[i]=m_PP[i];
} else if(m_KeyClone == 2) {
   m_PP.SetPxPyPzE(Pd[0],Pd[1],Pd[2],Pd[3]);
   if(     Kto == 1) {
     m_Event->RotEuler(m_alfa1, m_beta1, m_gamma1, &m_PP);
   } else if( Kto == 2) {
     m_Event->RotEuler(m_alfa2, m_beta2, m_gamma2, &m_PP);
   } else {
   (*m_Out)<<"###### STOP in TauPair::Tralo4: Wrong Kto = "<<Kto<<endl;
   cout    <<"###### STOP in TauPair::Tralo4: Wrong Kto = "<<Kto<<endl; exit(9);
   }
   for(int i=0; i<4;i++) Pd[i]=m_PP[i];
} else {
   (*m_Out)<<"##### STOP in Taupair_Tralo4: wrong KeyClone="<<m_KeyClone<<endl;
   cout    <<"##### STOP in Taupair_Tralo4: wrong KeyClone="<<m_KeyClone<<endl; exit(9);
}
m_GPS->TralorDoIt(Kto,Pd,Pd);
// Translation from DOUBLE PRECISION  to REAL
for(int k=0; k<4; k++) P[k]=Pd[k];
//----------------------------------------------
}//Tralo4
