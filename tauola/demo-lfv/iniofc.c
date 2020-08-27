#include <cstdlib>
#include <cstdio>
#include <complex>
#include "../tauola-c/ChannelForTauolaInterface.h"
#include "../tauola-c/ChannelForTauola.h"
using namespace std;
using namespace Tauolapp;

void lfv_dam2pi(const float *pt, const float *pn, const float *pim1, const float *pim2, float &amplit, float *hv);

// Channel redefinition function. This function adds the lfv_dam2pi channel
// defined in lfv.c file.
void AddLFV() {

    vector<int> products(3);

    //****************************************************************************************************************
    // Change initialization for the already existing channel, activate  user provided ME                            *
    //****************************************************************************************************************
    products[0] = -13;
    products[1] = -13;
    products[2] =  13; // note that this is third decay products id for the 
                       // essentially 2-scalar decay channel. This id will 
                       // override id of tau-neutrino (first tau decay product)
  
    ChannelForTauola *lfv_test1 = new ChannelForTauola(0.0, products, "  Dalitz for Tau->3mu LFV", lfv_dam2pi );

    Tauolapp::RegisterChannel( 90, lfv_test1 ); // NOTE: in 'prod/dane.dat' this number will be
                                                //       90 + 2 leptonic channels = 92

    if( lfv_test1->getChannelNo() == 0 ) printf("Couldn't register new channel!\n");
    else{       printf("\nRegistered new channel. It is now: \n");
                Tauolapp::PrintChannelInfo(lfv_test1->getChannelNo());
    }
}

// Set pointer for user channel redefinition
// This function is called from taumain.f
extern "C" void tauolaredef_() {
    Tauolapp::SetUserRedefinitions(AddLFV);
}

