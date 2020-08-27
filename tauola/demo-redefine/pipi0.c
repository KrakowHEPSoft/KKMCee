#include <math.h>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <complex>
#include "MEutils.h"
#include "pipi0.h"
using namespace std;

    // this current is C++ copy of FORTRAN current existing in tauola
    // It is described in: S. Jadach, Z. Was, J. H. Kuhn, 
    // TAUOLA - A library of MC programs to simulate decays of polarized tau leptons
    // Comp. Phys. Commun. 64 (1991) 275, (CERN-TH-5856 preprint)
    // references to exact equation are given in code bellow
void pipi0_curr(const float *pc, const float *pn, complex<float> *hadcur)    {
    float pks[4] ={0,0,0,0};
    float qq[4]  ={0,0,0,0};
    float pksd   =0;
    float qqpks  =0;
    int i;
    for (i=0;i<4;i++){
        pks[i]=pc[i]+pn[i];  // as defined on page 13, 
        qq[i] =pc[i]-pn[i];  // between eq (3.30) and (3.31)
    }
    pksd  = pks[3]*pks[3]-pks[2]*pks[2]-pks[1]*pks[1]-pks[0]*pks[0]; 
    qqpks = pks[3]*qq[3]-pks[2]*qq[2]-pks[1]*qq[1]-pks[0]*qq[0];
    complex<float> fpirho;
    fpirho=abs(fpik(sqrt(pksd)));          // fpik = Fpi from eq (3.34) 
                                           // cos Cabibo and sqrt(2) presten in paper with above eq. 
                                           // are moved to definition of ME (in def. of "amplit" and "brak", see MEutils.c )
    for (i=0;i<4;i++){                      
        qq[i]    =qq[i]-pks[i]*qqpks/pksd; // ortogonalization ( QQni*(Gmini-QQniQQmi/QQ^2) )
        hadcur[i]=qq[i]*fpirho;            // mentioned on page 15
    }
    return;
}

// this ME is C++ copy of FORTRAN ME existing in tauola 
void pipi0_ME(const float *pt, const float *pn, const float *pim1, const float *pim2, float &amplit, float *hv) {

    complex<float> hadcur[4]  = { 0.0 };

    pipi0_curr(pim1, pim2, hadcur);
    int ifcosCabibo=1;  // we declare whether cosine or sine of Cabibo angle should be used. In this case it is Cosine
    tauDecay_ME(ifcosCabibo,pt,pn,hadcur,amplit,hv);
}
