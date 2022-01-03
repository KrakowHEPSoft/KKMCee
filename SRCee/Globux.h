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
//  Emulating f77 externals using c++ objects
//-----------------------------------------------------------
extern "C" {
    // getting accsess to MC generator object in f77 code
    void globux_setup_( KKee2f*);
    // used in HepEvt
    void globux_getevent_( int *KFini, int *KFfin, double *MFini, double *MFfin,
                           double pf1[], double pf2[],double qf1[], double qf2[],
                           int *nphx, int *nphy,
                           double xphot[][100], double yphot[][100],
                           double rem1[], double rem2[]);
    //
    // wrappers to c++ code used in Tauola
    // already obsolete
    //void gps_makerho2_(double HvecFer1[], double HvecFer2[], double *wt0, double *wt1, double *wt2);
    //void gps_tralordoit_(int *KTO, double P[], double Q[]);
    //void gps_tralorprepare_(int *);
    // remains for some time
    void tralo4_(int *KTO, float P[], float Q[], float *AM);
}

//////////////////////////////////////////////
// Genuine new global functions (temporary)
// getting accsess to MC gener. object
void globux_setup_(     KKee2f *KKhhObj){   g_KKeeGen = KKhhObj;}

/////////////////////////////////////////////
// CALL GlobuxGetEvent(KFini,KFfin)
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

/*
void gps_makerho2_(double HvecFer1[], double HvecFer2[], double *wt0, double *wt1, double *wt2){
	g_KKeeGen->m_GPS->MakeRho2( HvecFer1, HvecFer2, *wt0, *wt1, *wt2);
}

void gps_tralorprepare_(int *KTO){
	g_KKeeGen->m_GPS->TralorPrepare(*KTO);
}

void gps_tralordoit_(int *KTO, double P[], double Q[]){
	g_KKeeGen->m_GPS->TralorDoIt(*KTO, P, Q);
}
*/

//SUBROUTINE Tralo4(Kto,P,Q,AM)
void tralo4_(int *KTO, float P[], float Q[], float *AM){
    g_KKeeGen->m_TauGen->Tralo4(*KTO, P, Q, *AM);
}

#endif
