#include <cstdlib>
#include <cstdio>
#include <complex>
#include "ChannelForTauola.h"
#include "ChannelForTauolaInterface.h"
typedef std::complex <float> complex;
using namespace Tauolapp;

// prints for technical tests
void print_4vector(const char *name, const float *v) {
    printf("%8s %+11.8f %+11.8f %+11.8f %+11.8f\n",name,v[0],v[1],v[2],v[3]);
}

void print_4vector_double(const char *name, const double *v) {
    printf("%8s %+11.8f %+11.8f %+11.8f %+11.8f\n",name,v[0],v[1],v[2],v[3]);
}

void print_4vectorC(const char *name, const complex *v) {
    printf("%8s ( %+11.8f %+11.8f ) ( %+11.8f %+11.8f )  ( %+11.8f %+11.8f ) ( %+11.8f %+11.8f ) \n",name,
	   v[0].real(),v[0].imag(),
	   v[1].real(),v[1].imag(),
	   v[2].real(),v[2].imag(),
	   v[3].real(),v[3].imag());
}

/*  
   ========================================================================
   The complete list  of all methods which are called from Fortran follow. 
   ========================================================================
*/

// Hidden flag that blocks redefinition outside of 'void iniofc_()' below
namespace Tauolapp { extern bool isInChannelRedefinitionFunction; }

/** Called by user program during initialization.
    See e.g. demo-redefine, patch-KK-face, patch-tauolapp */
extern "C" void iniofc_() {
    if(Tauolapp::channelRedefinitionFunction) {
        Tauolapp::isInChannelRedefinitionFunction = true;
        Tauolapp::channelRedefinitionFunction();
        Tauolapp::isInChannelRedefinitionFunction = false;
    }
}

extern "C" void dampry_wrap_(int *ITDKRC, double *XK0DEC, double *XK, double *XA, double *QP, double *XN, double *AMPLIT, double *HV) {
    // input parameters
    const int    &itdkrc = *ITDKRC;// radiative correction switch
    const double &xk0dec = *XK0DEC;
    const double *xk     = XK;     // photon       4-momentum
    const double *xa     = XA;     // antineutrino 4-momentum
    const double *qp     = QP;     // lepton       4-momentum
    const double *xn     = XN;     // neutrino     4-momentum

    // output parameters
    double &amplit = *AMPLIT;   // amplitude
    double *hv     = HV;        // polarimetric 4-vector

    int    mnum = 1;
    double mass = sqrt( qp[3]*qp[3] - qp[2]*qp[2] - qp[1]*qp[1] - qp[0]*qp[0] );

    if( mass > 0.1 ) mnum = 2;

    if( Tauolapp::leptonChannelsMEpointers[mnum-1] ) {
        // WARNING! Right now there is no type safety. We ASSUME pointer is always of valid type
        DAMPRY_POINTER_TYPE xsec = reinterpret_cast<DAMPRY_POINTER_TYPE>(Tauolapp::leptonChannelsMEpointers[mnum-1]);

        if(!xsec) {
            printf("ERROR: dampry_wrap: invalid pointer to xsec for mnum %i\n",mnum);
            exit(-1);
        }

        xsec(itdkrc,xk0dec,xk,xa,qp,xn,amplit,hv);
    }
    else {
        printf("ERROR: dampry_wrap: pointer to xsec for mnum %i not set!\n",mnum);
        exit(-1);
    }

    // check for NaN-s
    if( amplit!=amplit || hv[0]!=hv[0] || hv[1]!=hv[1] || hv[2]!=hv[2] || hv[3]!=hv[3] ) {
        printf("(C++) dampry wrapper. MNUM: %i.\n",mnum);
        printf("(C++) Calculation failed. Input:\n");
        printf(" ITDKRC: %i\n",itdkrc);
        printf(" XK0DEC: %+11.8f\n",xk0dec);
        print_4vector_double("XK:",xk);
        print_4vector_double("XA:",xa);
        print_4vector_double("QP:",qp);
        print_4vector_double("XN:",xn);
        printf("\n");
        printf("(C++) Output:\n");
        printf(" AMPLIT: %+11.8f\n",amplit);
        print_4vector_double("HV:",hv);
    }
}

extern "C" void dam1pi_wrap_(int *MNUM, float *PNU, float *AMF0, float *PKK, float *AMF1, float *GAMM, float *HV) {
    // input parameters
    const int   &mnum = *MNUM; // decay mode identifier
    const float *pnu   = PNU;  // 'neutrino' 4-momentum
    const float &amf0 = *AMF0; // its mass
    const float *pkk   = PKK;  // second product 4-momentum
    const float &amf1 = *AMF1; // its mass

    // output parameters
    float &gamm = *GAMM;   // amplitude
    float *hv   = HV;      // polarimetric 4-vector

    int channel = NM4+NM5+NM6+NM3+NM2+mnum;

    ChannelForTauola *userChannel = Tauolapp::taubra_userChannels[channel];
    if( userChannel ) {
        DAM1PI_POINTER_TYPE xsec = reinterpret_cast<DAM1PI_POINTER_TYPE>(userChannel->getFunctionPointer());

        // This should never happen but we check anyway
        if(!xsec) {
            printf("ERROR: dam1pi_wrap: invalid pointer to xsec for channel %i (mnum %i)\n",channel, mnum);
            exit(-1);
        }

        xsec(pnu,amf0,pkk,amf1,gamm,hv);
    }
    else {
        printf("ERROR: dam1pi_wrap: pointer to xsec for channel %i (mnum %i) not set!\n",channel,mnum);
        exit(-1);
    }

    // check for NaN-s
    if( gamm!=gamm || hv[0]!=hv[0] || hv[1]!=hv[1] || hv[2]!=hv[2] || hv[3]!=hv[3] ) {
        printf("(C++) DAM1PI wrapper. MNUM: %i.\n",mnum);
        printf("(C++) Calculation failed. Input:\n");
        print_4vector("PNU:",pnu);
        printf("AMF0: %f\n",amf0);
        print_4vector("PKK:",pkk);
        printf("AMF1: %f\n",amf1);
        printf("\n");
        printf("(C++) Output:\n");
        printf(" GAMM: %+11.8f\n",gamm);
        print_4vector("HV:",hv);
    }
}

extern "C" void dam2pi_wrap_(int *MNUM, float *PT, float *PN, float *PIM1, float *PIM2, float *AMPLIT, float *HV) {
    // input parameters
    const int   &mnum = *MNUM; // decay mode identifier
    const float *pt   = PT;    // tau 4-momentum
    const float *pn   = PN;    // tau neutrino 4-momentum
    const float *pim1 = PIM1;  // 1st particle 4-momentum
    const float *pim2 = PIM2;  // 2nd particle 4-momentum

    // output parameters
    float &amplit = *AMPLIT;   // amplitude
    float *hv     = HV;        // polarimetric 4-vector

    int channel = NM4+NM5+NM6+NM3+mnum;

    ChannelForTauola *userChannel = Tauolapp::taubra_userChannels[channel];
    if( userChannel ) {
        DAM2PI_POINTER_TYPE xsec = reinterpret_cast<DAM2PI_POINTER_TYPE>(userChannel->getFunctionPointer());

        // This should never happen but we check anyway
        if(!xsec) {
            printf("ERROR: dam2pi_wrap: invalid pointer to xsec for channel %i (mnum %i)\n",channel, mnum);
            exit(-1);
        }

        xsec(pt,pn,pim1,pim2,amplit,hv);
    }
    else {
        printf("ERROR: dam2pi_wrap: pointer to xsec for channel %i (mnum %i) not set!\n",channel,mnum);
        exit(-1);
    }

    // check for NaN-s
    if( amplit!=amplit || hv[0]!=hv[0] || hv[1]!=hv[1] || hv[2]!=hv[2] || hv[3]!=hv[3] ) {
        printf("(C++) DAM2PI wrapper. MNUM: %i.\n",mnum);
        printf("(C++) Calculation failed. Input:\n");
        print_4vector("PT:",pt);
        print_4vector("PN:",pn);
        print_4vector("PIM1:",pim1);
        print_4vector("PIM2:",pim2);
        printf("\n");
        printf("(C++) Output:\n");
        printf(" AMPLIT: %+11.8f\n",amplit);
        print_4vector("HV:",hv);
    }
}

extern "C" void dam3pi_wrap_(int *MNUM, float *PT, float *PN, float *PIM1, float *PIM2, float *PIM3, float *AMPLIT, float *HV) {
    // input parameters
    const int   &mnum = *MNUM; // decay mode identifier
    const float *pt   = PT;    // tau 4-momentum
    const float *pn   = PN;    // tau neutrino 4-momentum
    const float *pim1 = PIM1;  // 1st particle 4-momentum
    const float *pim2 = PIM2;  // 2nd particle 4-momentum
    const float *pim3 = PIM3;  // 3rd particle 4-momentum

    // output parameters
    float &amplit = *AMPLIT;   // amplitude
    float *hv     = HV;        // polarimetric 4-vector

    int channel = NM4+NM5+NM6+mnum;

    ChannelForTauola *userChannel = Tauolapp::taubra_userChannels[channel];
    if( userChannel ) {
        DAM3PI_POINTER_TYPE xsec = reinterpret_cast<DAM3PI_POINTER_TYPE>(userChannel->getFunctionPointer());

        // This should never happen but we check anyway
        if(!xsec) {
            printf("ERROR: dam3pi_wrap: invalid pointer to xsec for channel %i (mnum %i)\n",channel, mnum);
            exit(-1);
        }

        xsec(pt,pn,pim1,pim2,pim3,amplit,hv);
    }
    else {
        printf("ERROR: dam3pi_wrap: pointer to xsec for channel %i (mnum %i) not set!\n",channel,mnum);
        exit(-1);
    }

    // check for NaN-s
    if( amplit!=amplit || hv[0]!=hv[0] || hv[1]!=hv[1] || hv[2]!=hv[2] || hv[3]!=hv[3] ) {
        printf("(C++) DAM3PI wrapper. MNUM: %i.\n",mnum);
        printf("(C++) Calculation failed. Input:\n");
        print_4vector("PT:",pt);
        print_4vector("PN:",pn);
        print_4vector("PIM1:",pim1);
        print_4vector("PIM2:",pim2);
        print_4vector("PIM3:",pim3);
        printf("\n");
        printf("(C++) Output:\n");
        printf(" AMPLIT: %+11.8f\n",amplit);
        print_4vector("HV:",hv);
    }
}

extern "C" void dam4pi_wrap_(int *MNUM, float *PT, float *PN, float *PIM1, float *PIM2, float *PIM3, float *PIM4, float *AMPLIT, float *HV) {
    // input parameters
    const int   &mnum = *MNUM; // decay mode identifier
    const float *pt   = PT;    // tau 4-momentum
    const float *pn   = PN;    // tau neutrino 4-momentum
    const float *pim1 = PIM1;  // 1st particle 4-momentum
    const float *pim2 = PIM2;  // 2nd particle 4-momentum
    const float *pim3 = PIM3;  // 3rd particle 4-momentum
    const float *pim4 = PIM4;  // 4th particle 4-momentum

    // output parameters
    float &amplit = *AMPLIT;   // amplitude
    float *hv     = HV;        // polarimetric 4-vector

    int channel = mnum;

    ChannelForTauola *userChannel = Tauolapp::taubra_userChannels[channel];
    if( userChannel ) {
        DAM4PI_POINTER_TYPE xsec = reinterpret_cast<DAM4PI_POINTER_TYPE>(userChannel->getFunctionPointer());

        // This should never happen but we check anyway
        if(!xsec) {
            printf("ERROR: dam4pi_wrap: invalid pointer to xsec for channel %i (mnum %i)\n",channel, mnum);
            exit(-1);
        }

        xsec(pt,pn,pim1,pim2,pim3,pim4,amplit,hv);
    }
    else {
        printf("ERROR: dam4pi_wrap: pointer to xsec for channel %i (mnum %i) not set!\n",channel,mnum);
        exit(-1);
    }

    // check for NaN-s
    if( amplit!=amplit || hv[0]!=hv[0] || hv[1]!=hv[1] || hv[2]!=hv[2] || hv[3]!=hv[3] ) {
        printf("(C++) DAM4PI wrapper. MNUM: %i.\n",mnum);
        printf("(C++) Calculation failed. Input:\n");
        print_4vector("PT:",pt);
        print_4vector("PN:",pn);
        print_4vector("PIM1:",pim1);
        print_4vector("PIM2:",pim2);
        print_4vector("PIM3:",pim3);
        print_4vector("PIM4:",pim4);
        printf("\n");
        printf("(C++) Output:\n");
        printf(" AMPLIT: %+11.8f\n",amplit);
        print_4vector("HV:",hv);
    }
}

extern "C" void dam5pi_wrap_(int *MNUM, float *PT, float *PN, float *PIM1, float *PIM2, float *PIM3, float *PIM4, float *PIM5, float *AMPLIT, float *HV) {
    // input parameters
    const int   &mnum = *MNUM; // decay mode identifier
    const float *pt   = PT;    // tau 4-momentum
    const float *pn   = PN;    // tau neutrino 4-momentum
    const float *pim1 = PIM1;  // 1st particle 4-momentum
    const float *pim2 = PIM2;  // 2nd particle 4-momentum
    const float *pim3 = PIM3;  // 3rd particle 4-momentum
    const float *pim4 = PIM4;  // 4th particle 4-momentum
    const float *pim5 = PIM5;  // 5th particle 4-momentum

    // output parameters
    float &amplit = *AMPLIT;   // amplitude
    float *hv     = HV;        // polarimetric 4-vector

    int channel = mnum+NM4; 

    ChannelForTauola *userChannel = Tauolapp::taubra_userChannels[channel];
    if( userChannel ) {
        DAM5PI_POINTER_TYPE xsec = reinterpret_cast<DAM5PI_POINTER_TYPE>(userChannel->getFunctionPointer());

        // This should never happen but we check anyway
        if(!xsec) {
            printf("ERROR: dam5pi_wrap: invalid pointer to xsec for channel %i (mnum %i)\n",channel, mnum);
            exit(-1);
        }

        xsec(pt,pn,pim1,pim2,pim3,pim4,pim5,amplit,hv);
    }
    else {
        printf("ERROR: dam5pi_wrap: pointer to xsec for channel %i (mnum %i) not set!\n",channel,mnum);
        exit(-1);
    }

    // check for NaN-s
    if( amplit!=amplit || hv[0]!=hv[0] || hv[1]!=hv[1] || hv[2]!=hv[2] || hv[3]!=hv[3] ) {
        printf("(C++) DAM5PI wrapper. MNUM: %i.\n",mnum);
        printf("(C++) Calculation failed. Input:\n");
        print_4vector("PT:",pt);
        print_4vector("PN:",pn);
        print_4vector("PIM1:",pim1);
        print_4vector("PIM2:",pim2);
        print_4vector("PIM3:",pim3);
        print_4vector("PIM4:",pim4);
        print_4vector("PIM5:",pim5);
        printf("\n");
        printf("(C++) Output:\n");
        printf(" AMPLIT: %+11.8f\n",amplit);
        print_4vector("HV:",hv);
    }
}

extern "C" void curr2_wrap_(int *MNUM, float *PIM1, float *PIM2, complex *HADCUR) {
    // input parameters
    const int   &mnum = *MNUM; // decay mode identifier
    const float *pim1 = PIM1;  // 1st particle 4-momentum
    const float *pim2 = PIM2;  // 2nd particle 4-momentum

    // output parameters
    complex *hadcur = HADCUR;   // hadronic current

    int channel = NM4+NM5+NM6+NM3+mnum;

    ChannelForTauola *userChannel = Tauolapp::taubra_userChannels[channel];
    if( userChannel ) {
        CURR2_POINTER_TYPE current = reinterpret_cast<CURR2_POINTER_TYPE>(userChannel->getFunctionPointer());

        // This should never happen but we check anyway
        if(!current) {
            printf("ERROR: curr2_wrap: invalid pointer to current for channel %i (mnum %i)\n",channel,mnum);
            exit(-1);
        }

        current(pim1,pim2,hadcur);
    }
    else {
        printf("ERROR: curr2_wrap: pointer to current for channel %i (mnum %i) not set!\n",channel,mnum);
        exit(-1);
    }

    // check for NaN-s
    if( hadcur[0]!=hadcur[0] || hadcur[1]!=hadcur[1] || hadcur[2]!=hadcur[2] || hadcur[3]!=hadcur[3] ) {
        printf("(C++) curr2 wrapper. MNUM: %i.\n",mnum);
        printf("(C++) Calculation failed. Input:\n");
        print_4vector("PIM1:",pim1);
        print_4vector("PIM2:",pim2);
        printf("\n");
        printf("(C++) Output:\n");
        print_4vectorC("complex HADCUR:",hadcur);
    }
}

extern "C" void curr3pi_wrap_(int *MNUM, float *PIM1, float *PIM2, float *PIM3, complex *HADCUR) {
    // input parameters
    const int   &mnum = *MNUM; // decay mode identifier
    const float *pim1 = PIM1;  // 1st particle 4-momentum
    const float *pim2 = PIM2;  // 2nd particle 4-momentum
    const float *pim3 = PIM3;  // 3rd particle 4-momentum

    // output parameters
    complex *hadcur = HADCUR;   // hadronic current

    int channel = NM4+NM5+NM6+mnum;

    ChannelForTauola *userChannel = Tauolapp::taubra_userChannels[channel];
    if( userChannel ) {
        CURR3_POINTER_TYPE current = reinterpret_cast<CURR3_POINTER_TYPE>(userChannel->getFunctionPointer());

        // This should never happen but we check anyway
        if(!current) {
            printf("ERROR: curr3pi_wrap: invalid pointer to current for channel %i (mnum %i)\n",channel,mnum);
            exit(-1);
        }

        current(pim1,pim2,pim3,hadcur);
    }
    else {
        printf("ERROR: curr3pi_wrap: pointer to current for channel %i (mnum %i) not set!\n",channel,mnum);
        exit(-1);
    }

    // check for NaN-s
    if( hadcur[0]!=hadcur[0] || hadcur[1]!=hadcur[1] || hadcur[2]!=hadcur[2] || hadcur[3]!=hadcur[3] ) {
        printf("(C++) curr3 wrapper. MNUM: %i.\n",mnum);
        printf("(C++) Calculation failed. Input:\n");
        print_4vector("PIM1:",pim1);
        print_4vector("PIM2:",pim2);
        print_4vector("PIM3:",pim3);
        printf("\n");
        printf("(C++) Output:\n");
        print_4vectorC("complex HADCUR:",hadcur);
    }
}

extern "C" void curr4_wrap_(int *MNUM, float *PIM1, float *PIM2, float *PIM3, float *PIM4, complex *HADCUR) {
    // input parameters
    const int   &mnum = *MNUM; // decay mode identifier
    const float *pim1 = PIM1;  // 1st particle 4-momentum
    const float *pim2 = PIM2;  // 2nd particle 4-momentum
    const float *pim3 = PIM3;  // 3rd particle 4-momentum
    const float *pim4 = PIM4;  // 4th particle 4-momentum

    // output parameters
    complex *hadcur = HADCUR;   // hadronic current

    int channel = mnum;

    ChannelForTauola *userChannel = Tauolapp::taubra_userChannels[channel];
    if( userChannel ) {
        CURR4_POINTER_TYPE current = reinterpret_cast<CURR4_POINTER_TYPE>(userChannel->getFunctionPointer());

        // This should never happen but we check anyway
        if(!current) {
            printf("ERROR: curr4_wrap: invalid pointer to current for channel %i (mnum %i)\n",channel,mnum);
            exit(-1);
        }

        current(pim1,pim2,pim3,pim4,hadcur);
    }
    else {
        printf("ERROR: curr4_wrap: pointer to current for channel %i (mnum %i) not set!\n",channel,mnum);
        exit(-1);
    }

    // check for NaN-s
    if( hadcur[0]!=hadcur[0] || hadcur[1]!=hadcur[1] || hadcur[2]!=hadcur[2] || hadcur[3]!=hadcur[3] ) {
        printf("(C++) curr4 wrapper. MNUM: %i.\n",mnum);
        printf("(C++) Calculation failed. Input:\n");
        print_4vector("PIM1:",pim1);
        print_4vector("PIM2:",pim2);
        print_4vector("PIM3:",pim3);
        print_4vector("PIM4:",pim4);
        printf("\n");
        printf("(C++) Output:\n");
        print_4vectorC("complex HADCUR:",hadcur);
    }
}

extern "C" void curr5_wrap_(int *MNUM, float *PIM1, float *PIM2, float *PIM3, float *PIM4, float *PIM5, complex *HADCUR) {
    // input parameters
    const int   &mnum = *MNUM; // decay mode identifier
    const float *pim1 = PIM1;  // 1st particle 4-momentum
    const float *pim2 = PIM2;  // 2nd particle 4-momentum
    const float *pim3 = PIM3;  // 3rd particle 4-momentum
    const float *pim4 = PIM4;  // 4th particle 4-momentum
    const float *pim5 = PIM5;  // 5th particle 4-momentum

    // output parameters
    complex *hadcur = HADCUR;   // hadronic current

    int channel = mnum+NM4; 

    ChannelForTauola *userChannel = Tauolapp::taubra_userChannels[channel];
    if( userChannel ) {
        CURR5_POINTER_TYPE current = reinterpret_cast<CURR5_POINTER_TYPE>(userChannel->getFunctionPointer());

        // This should never happen but we check anyway
        if(!current) {
            printf("ERROR: curr5_wrap: invalid pointer to current for channel %i (mnum %i)\n",channel,mnum);
            exit(-1);
        }

        current(pim1,pim2,pim3,pim4,pim5,hadcur);
    }
    else {
        printf("ERROR: curr5_wrap: pointer to current for channel %i (mnum %i) not set!\n",channel,mnum);
        exit(-1);
    }

    // check for NaN-s
    if( hadcur[0]!=hadcur[0] || hadcur[1]!=hadcur[1] || hadcur[2]!=hadcur[2] || hadcur[3]!=hadcur[3] ) {
        printf("(C++) curr5 wrapper. MNUM: %i.\n",mnum);
        printf("(C++) Calculation failed. Input:\n");
        print_4vector("PIM1:",pim1);
        print_4vector("PIM2:",pim2);
        print_4vector("PIM3:",pim3);
        print_4vector("PIM4:",pim4);
        print_4vector("PIM5:",pim5);
        printf("\n");
        printf("(C++) Output:\n");
        print_4vectorC("complex HADCUR:",hadcur);
    }
}

extern "C" float sigee_wrap_(float *AMX2, int *MNUM) {
    // input parameters
    const float &amx2 = *AMX2; // virtuality
    const int   &mnum = *MNUM; // decay mode identifier
    
    float result = 0.0;

    int channel = NM4+NM5+mnum-3; // in tauola.f: jn=JNPI-nm4-nm5+3
                                  // this 3 is a consequence of historical
                                  // choices and interfaces to LabAccelLinear
                                  // of Orsay and of 70's sotware, 
                                  // first 3 channels are of ME treatment 
                                  // since decades. Here is a shadow.
    ChannelForTauola *userChannel = Tauolapp::taubra_userChannels[channel];
    if( userChannel ) {
        SIGEE_POINTER_TYPE sigee = reinterpret_cast<SIGEE_POINTER_TYPE>(userChannel->getFunctionPointer());

        // This should never happen but we check anyway
        if(!sigee) {
            printf("ERROR: sigee_wrap: invalid pointer to Q2-spectrum for channel %i (mnum %i)\n",channel,mnum);
            exit(-1);
        }

        result = sigee(amx2);
    }
    else {
        printf("ERROR: sigee_wrap: pointer to Q2-spectrum for channel %i (mnum %i) not set!\n",channel,mnum);
        exit(-1);
    }

    // check for NaN-s
    if( result!=result ) {
        printf("(C++) sigee wrapper. MNUM: %i.\n",mnum);
        printf("(C++) Calculation failed.\n");
        printf(" Input AMX2: %+11.8f\n",amx2);
        printf(" Output:     %+11.8f\n",result);
    }
    
    return result;
}

extern "C" float fconst_wrap_(int *MNUM) {
    // input parameters
    const int   &mnum = *MNUM; // decay mode identifier

    int channel = NM4+NM5+NM6+NM3+NM2+mnum;

    ChannelForTauola *userChannel = Tauolapp::taubra_userChannels[channel];
    if( userChannel ) {
        FCONST_POINTER_TYPE coupl = reinterpret_cast<FCONST_POINTER_TYPE>(userChannel->getFunctionPointer());

        // This should never happen but we check anyway
        if(!coupl) {
            printf("ERROR: fconst_wrap: invalid pointer to norm. const. for channel %i (mnum %i)\n",channel,mnum);
            exit(-1);
        }

        return coupl();
    }
    else {
        printf("ERROR: fconst_wrap: pointer to current for norm. const. %i (mnum %i) not set!\n",channel,mnum);
        exit(-1);
    }
}
