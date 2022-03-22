#ifndef Globux_H
#define Globux_H
////////////////////////////////////////////////////////////////////////
//
//  Temporary global interface between C++ and F77 parts of the code
//
////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <unistd.h>
#include <iostream>
#include <string>
using namespace std;

#include "KKee2f.h"
#include "KKdizet.h"

KKee2f   *g_KKeeGen;
KKdizet  *g_EWtabs;

//-----------------------------------------------------------
//  Emulating f77 externals using c++ objects/ members
//-----------------------------------------------------------
extern "C" {
    // getting accsess to MC generator object in f77 code
    void globux_setup_( KKee2f*);
    // used in HepEvt to import c++ event object
    void globux_getevent_( int *KFini, int *KFfin, double *MFini, double *MFfin,
                           double pf1[], double pf2[],double qf1[], double qf2[],
                           int *nphx, int *nphy,
                           double xphot[][100], double yphot[][100],
                           double rem1[], double rem2[]);
    // wrapper to c++ code used in Tauola
    void tralo4_(int *KTO, float P[], float Q[], float *AM);
    void fillhep3_(int *N, int *IST, int *ID, int *JMO1, int *JMO2, int *JDA1, int *JDA2, float P4[], float *PINV, bool *PHFLAG);
    void ranmar_(float rvec[], int *lenv);
}

//////////////////////////////////////////////
// Genuine new global functions (temporary)
// getting accsess to MC gener. object
void globux_setup_(     KKee2f *KKhhObj){   g_KKeeGen = KKhhObj;}

////////////////////////////////////////////////////////////////////////
/// Importing c++ event object into /hepevt/ of f77 HepEvt
void globux_getevent_( int *KFini, int *KFfin, double *MFini, double *MFfin,
                       double pf1[], double pf2[],double qf1[], double qf2[],
                       int *nphx, int *nphy,
                       double xphot[][100], double yphot[][100],
                       double rem1[], double rem2[]){
   *KFini = g_KKeeGen->m_Event->m_KFini;
   *KFfin = g_KKeeGen->m_Event->m_KFfin;
   *MFini = g_KKeeGen->DB->fmass[*KFini];
   *MFfin = g_KKeeGen->DB->fmass[*KFfin];
// DOUBLE PRECISION  pf1(4),pf2(4),qf1(4),qf2(4),xphot(100,4),yphot(100,4)
   for(int i=0; i<=3;i++) pf1[i] = g_KKeeGen->m_Event->m_Pf1[i];
   for(int i=0; i<=3;i++) pf2[i] = g_KKeeGen->m_Event->m_Pf2[i];
   for(int i=0; i<=3;i++) qf1[i] = g_KKeeGen->m_Event->m_Qf1[i];
   for(int i=0; i<=3;i++) qf2[i] = g_KKeeGen->m_Event->m_Qf2[i];
   for(int i=0; i<=3;i++) rem1[i] = g_KKeeGen->m_Event->m_Rem1[i];
   for(int i=0; i<=3;i++) rem2[i] = g_KKeeGen->m_Event->m_Rem2[i];
   *nphx = g_KKeeGen->m_Event->m_nPhotISR;
   *nphy = g_KKeeGen->m_Event->m_nPhotFSR;
   for(int j=1;j<=*nphx;j++)
      for(int i=0; i<=3;i++) xphot[i][j-1] = g_KKeeGen->m_Event->m_PhotISR[j][i];
   for(int j=1;j<=*nphy;j++)
      for(int i=0; i<=3;i++) yphot[i][j-1] = g_KKeeGen->m_Event->m_PhotFSR[j][i];
}//globux_getevent_


////////////////////////////////////////////
//SUBROUTINE Tralo4(Kto,P,Q,AM)
void tralo4_(int *KTO, float P[], float Q[], float *AM){
    g_KKeeGen->m_TauGen->Tralo4(*KTO, P, Q, *AM);
}//tralo4

void fillhep3_(int *N, int *IST, int *ID, int *JMO1, int *JMO2, int *JDA1, int *JDA2, float P4[], float *PINV, bool *PHFLAG){
    g_KKeeGen->m_HEPMC->FillHep3(*N, *IST, *ID, *JMO1, *JMO2, *JDA1, *JDA2, P4, *PINV, *PHFLAG);
}//fillhep3_


// SUBROUTINE RANMAR(RVEC,LENV)
void ranmar_(float rvec[], int *lenv){
   if(*lenv >1000) {
       cout<<"+++ Stop in RANMAR wrapper in globux.h: lenv>1000"<<endl;exit(9);};
   double drvec[1000];
   g_KKeeGen->f_RNgen->RndmArray(*lenv, drvec);
   for(int i=0;i<*lenv;i++) rvec[i]=drvec[i];
}
#endif
