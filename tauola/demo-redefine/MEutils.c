#include <math.h>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <complex>
#include "MEutils.h"
using namespace std;



//User may access parameters taken from TAUOLA fortan code declaring  commons/structs for C++ use
/*
extern "C" struct {
    float GFERMI, GV, GA, CCABIB, SCABIB, GAMEL;
} decpar_;
extern "C" struct {
    float AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1,AMK,AMKZ,AMKST,GAMKST;
} parmas_;
*/

//*******************************************************************************
//Functions below are copies of routines from TAUOLA-FORTRAN/tauola-bbb/tauola.f
//*******************************************************************************

void clvec(complex<float> *HJ, const float *PN, float *PIV)    {
    // ----------------------------------------------------------------------
    // CALCULATES THE "VECTOR TYPE"  PI-VECTOR  PIV
    // NOTE THAT THE NEUTRINO MOM. PN IS ASSUMED TO BE ALONG Z-AXIS
    // ----------------------------------------------------------------------
    complex<float> HN;
    float HH;
    HN= HJ[3]*PN[3]-HJ[2]*PN[2];
    HH= real(HJ[3]*conj(HJ[3])-HJ[2]*conj(HJ[2])-HJ[1]*conj(HJ[1])-HJ[0]*conj(HJ[0]));
    int i;
    for (i=0;i<4;i++) PIV[i]=4.*real(HN*conj(HJ[i]))-2.*HH*PN[i];
    return;
}



void claxi(complex<float> *HJ, const float *PN, float *PIA)    {
    // ----------------------------------------------------------------------
    // CALCULATES THE "AXIAL TYPE"  PI-VECTOR  PIA
    // NOTE THAT THE NEUTRINO MOM. PN IS ASSUMED TO BE ALONG Z-AXIS
    // SIGN is chosen +/- for decay of TAU +/- respectively
    // ----------------------------------------------------------------------
    complex<float> HJC[4];
    float DET2[4][4];
    int i, j;
    for (i=0;i<4;i++) HJC[i]=conj(HJ[i]);    
    for (i=0;i<4;i++){
        for (j=0;j<4;j++) DET2[i][j]=imag(HJC[i]*HJ[j]-HJC[j]*HJ[i]);
    }

     float SIGN = 0;
    // -- SIGN  was affecting sign of A_LR asymmetry in a1 decay.
    // -- note also collision of notation of gamma_va as defined in
    // -- TAUOLA paper and J.H. Kuhn and Santamaria Z. Phys C 48 (1990) 445
    // -----------------------------------

    taugetsign_(SIGN);
                    // SIGN depend on whether we generate decay of tau+ or of tau-
                    // the information has to be obtained form Tauola common blocks
                    // special function has to be provided for that purpose.
                    // This function is the only one which requires information from Tauola

    PIA[0]= -2*PN[2]*DET2[1][3]+2*PN[3]*DET2[1][2];
    PIA[1]= -2*PN[3]*DET2[0][2]+2*PN[2]*DET2[0][3];
    PIA[2]=  2*PN[3]*DET2[0][1];
    PIA[3]=  2*PN[2]*DET2[0][1];
    // ALL FOUR INDICES ARE UP SO  PIA(3) AND PIA(4) HAVE SAME SIGN
    for (i=0;i<4;i++) { 
    PIA[i]=PIA[i]*SIGN;
    }
    return;
}

void clnut(complex<float> *HJ, float &B, float *HV)    {
    // ----------------------------------------------------------------------
    // CALCULATES THE CONTRIBUTION BY NEUTRINO MASS
    // NOTE THE TAU IS ASSUMED TO BE AT REST
    // ----------------------------------------------------------------------
    float P[4];
    claxi(HJ,P,HV);
    
    B = real(HJ[3]*imag(HJ[3])-HJ[2]*imag(HJ[2])-HJ[1]*imag(HJ[1])-HJ[0]*imag(HJ[0]));
    return;
}


// this ME is C++ copy of FORTRAN ME existing in tauola it is universal, can be used for any of hadronic tau decay channels 
void tauDecay_ME(const int &ifcosCabibo, const float *pt, const float *pn, complex<float> *hadcur,  float &amplit, float *hv) {
    float GFERMI = 0.000016637; // float &GFERMI = decpar_.GFERMI;
    float gv     = 1;           // float &gv     = decpar_.GV;
    float ga     = -1;          // float &ga     = decpar_.GA;
    float AMNUTA = 0;           // float &AMNUTA = parmas_.AMNUTA;
    float AMTAU  = 1.777;       // float &AMTAU  = parmas_.AMTAU;
    
    float pivec[4], piaks[4], hvm[4];
    float brak = 0, brakm = 0,  cabi = 0;
    int i;


    if(ifcosCabibo==1) cabi = 0.975;    // cabi= decpar_.CCABIB;
    else               cabi = 0.222205; // cabi= decpar_.SCABIB;
    clvec(hadcur, pn, pivec);
    claxi(hadcur, pn, piaks);
    clnut(hadcur, brakm, hvm);
    brak=(gv*gv+ga*ga)*pt[3]*pivec[3]+2*gv*ga*pt[3]*piaks[3]+2*(gv*gv-ga*ga)*AMNUTA*AMTAU*brakm;
    amplit=(cabi*GFERMI)*(cabi*GFERMI)*brak;
    for (i=0;i<3;i++)
    hv[i]=-(-(AMTAU*((gv*gv+ga*ga)*piaks[i]+2*gv*ga*pivec[i]))+(gv*gv-ga*ga)*AMNUTA*AMTAU*hvm[i])/brak;
    hv[3]=1.0;
    return;
}


complex<float> bwig(float s, float m, float g)    { 
    // **********************************************************
    //     P-WAVE BREIT-WIGNER  FOR RHO
    // **********************************************************
    float PI, PIM, QS, QM, W, GS;
    // ------------ PARAMETERS --------------------
    PI   =3.141592654;
    PIM  =.139;
    // -------  BREIT-WIGNER -----------------------
    if (s>(4*PIM*PIM))    {
        QS=sqrt(abs(abs(s/4-PIM*PIM)+(s/4-PIM*PIM))/2.0);
        QM=sqrt(m*m/4-PIM*PIM);
        W =sqrt(s);
        GS=g*(m/W)*(QS/QM)*(QS/QM)*(QS/QM);
    }
    else GS=0.0;
    complex<float> BWIG;
    complex<float> cmplx=complex<float>(m*m-s, -m*GS);
    return BWIG=m*m/cmplx;
}
      
complex<float> fpik(float W)    {
    // **********************************************************
    //     PION FORM FACTOR
    // **********************************************************
    complex<float> BWIG;
    float ROM, ROG, ROM1, ROG1, BETA1, PI, PIM, S;
    // ------------ PARAMETERS --------------------
    PI   =3.141592654;
    PIM  =.140;
    // next constants are from BaBar
    ROM  =0.773;
    ROG  =0.145;
    ROM1 =1.370;
    ROG1 =0.510;
    BETA1=-0.145;
    // -----------------------------------------------
    S=W*W;
    complex<float> fpik;
    fpik = (bwig(S,ROM,ROG)+BETA1*bwig(S,ROM1,ROG1))/(1+BETA1);
    return fpik;
}
