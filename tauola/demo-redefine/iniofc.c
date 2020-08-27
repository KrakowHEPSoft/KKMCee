#include <cstdlib>
#include <cstdio>
#include <complex>
#include "../tauola-c/ChannelForTauolaInterface.h"
#include "../tauola-c/ChannelForTauola.h"
#include "TestCommunication.c"
#include "pipi0.h" //demo 2pi current+ME
using namespace std;
using Tauolapp::ChannelForTauola;

void lfv_dam2pi(const float *pt, const float *pn, const float *pim1, const float *pim2, float &amplit, float *hv);

// Examples for channel (re-)initialization
void RedefExample() {
    vector<int> products(9);

    
    //**************************************************************************
    // demo pi- pi0 current, it's a copy fortran current using pipi0_curr      *
    //**************************************************************************
/*    
    products[0] = -1; //pi-
    products[1] =  2; //pi0 
    ChannelForTauola *pipi0_test1 = new ChannelForTauola(9.0, products, "  demo pi- pi0  " , pipi0_curr );
    
    Tauolapp::RegisterChannel( 86, pipi0_test1 );
    if( pipi0_test1->getChannelNo() == 0 ) printf("Couldn't register demo pi- pi0 channel!\n");
    else{       printf("\nRegistered new channel. It is now: \n");
                Tauolapp::PrintChannelInfo( pipi0_test1->getChannelNo());
    }               
 */
    //************************************************************************************
    // demo pi- pi0 ME, it's a copy fortran ME implemented in tauola using pipi0_ME      *
    //************************************************************************************
     
    products[0] = -1; //pi-
    products[1] =  2; //pi0 
    ChannelForTauola *pipi0_test2 = new ChannelForTauola(9.0, products, "  demo pi- pi0  " , pipi0_ME );

    Tauolapp::RegisterChannel( 86, pipi0_test2 );
    if( pipi0_test2->getChannelNo() == 0 ) printf("Couldn't register demo pi- pi0 ME channel!\n");
    else{       printf("\nRegistered new channel. It is now: \n");
                Tauolapp::PrintChannelInfo(pipi0_test2->getChannelNo());
    }               


    //************************************************************************************
    // Register new channel no matrix element flat phase space                           *
    //************************************************************************************
    products[0] = -13;
    products[1] = -13;
    products[2] =  13; // note that this is third decay products id for the 
                       // essentially 2-scalar decay channel. This id will 
                       // override id of tau-neutrino (first tau decay product)
    ChannelForTauola *mumumu_test1 = new ChannelForTauola(2, 1, 0.1, products, "TEST_CHANNEL" );
    
    Tauolapp::RegisterChannel( -1, mumumu_test1 );
    if( mumumu_test1->getChannelNo() == 0 ) printf("Couldn't register new channel!\n");
    else{       printf("\nRegistered new channel. It is now: \n");
                Tauolapp::PrintChannelInfo(mumumu_test1->getChannelNo());
    }

    //****************************************************************************************************************
    // Change initialization for the already existing channel, activate SM ME using user provided hadronic current   *
    // WARNING: it is unphysical choice to use SM ME for 3 muon tau_neutrino-less channel                            *
    //****************************************************************************************************************
    products[0] = -13;
    products[1] = -13;
    products[2] =  13; // note that this is third decay products id for the 
                       // essentially 2-scalar decay channel. This id will 
                       // override id of tau-neutrino (first tau decay product)


    ChannelForTauola *mumumu_test2 = new ChannelForTauola(9.999, products, "TEST_CHANNEL_2", test_curr2 );
    
    Tauolapp::RegisterChannel( 154, mumumu_test2 );
    if( mumumu_test2->getChannelNo() == 0 ) printf("Couldn't register new channel!\n");
    else{       printf("\nRegistered new channel. It is now: \n");
                Tauolapp::PrintChannelInfo(mumumu_test2->getChannelNo());
    }



    //****************************************************************************************************************
    // Change initialization for the already existing channel, activate  user provided ME                            *
    //****************************************************************************************************************
    products[0] = -13;
    products[1] = -13;
    products[2] =  13; // note that this is third decay products id for the 
                       // essentially 2-scalar decay channel. This id will 
                       // override id of tau-neutrino (first tau decay product)
  

    ChannelForTauola *lfv_test1 = new ChannelForTauola(9.9991, products, "Dalitz for Tau->3mu LFV",    lfv_dam2pi );
    Tauolapp::RegisterChannel( 90, lfv_test1 );

    // alternative example: use position `92'
      //ChannelForTauola *lfv_test1 = new ChannelForTauola(9.9991, products, "Dalitz for Tau->3mu OS LFV", lfv_dam2pi );  
      //Tauolapp::RegisterChannel( 92, lfv_test1 );

    // alternative example: use position `91'
      //ChannelForTauola *lfv_test1 = new ChannelForTauola(9.9991, products, "Dalitz for Tau->3mu SS LFV", lfv_dam2pi );  
      //Tauolapp::RegisterChannel( 91, lfv_test1 );
    
    if( lfv_test1->getChannelNo() == 0 ) printf("Couldn't register new channel!\n");
    else{       printf("\nRegistered new channel. It is now: \n");
                Tauolapp::PrintChannelInfo(lfv_test1->getChannelNo());
    }


    //****************************************************************************************
    // Example of use basic functionalities allowing modification of existing channels
    // each of the individual operations can be used separtately
    //****************************************************************************************
    ChannelForTauola *demo_modify = GetChannel(87); //Here we create ChannelForTauola object from already existing FORTRAN channel
                                                    //Note: changes won't be included until modified channel registration
    // 1)
    string name = demo_modify->getName();           
    name ="renamed "+name;                          //Example of renaming channel
    demo_modify->setName(name);                         

    // 2)
    double br = demo_modify->getBr();               //Importing and modification of Branching Ratio
    demo_modify->setBr(br*1234);

    // 3)
    products = demo_modify->getProducts();          //Importing and modification of decay products
    products[0] = -3; //K-
    products[1] =  4; //K0
    demo_modify->setProducts(products);  

    // 4)           
    demo_modify->setName("TAU->K-K0 using demo pi-pi0 ME");   //let's re-name it again to be more accurate

    // 5)  THIS IS NOT RECOMMENDED/BLOCKED:    if pointer is used,  use previous group of examples 
    //demo_modify->setFunctionPointer( (void*)pipi0_ME );      //seting pointer to our demo pi- pi0 ME
    //demo_modify->setMeType(5);                               //ME type has to be set accordingly to used ME (5 since we use user-defined ME)

    Tauolapp::RegisterChannel( 87, demo_modify );    //registration back the modified demo_modify into its original postion 
    demo_modify->print();                            //prints channel info
    
    // 6)
    demo_modify->setMeType(1);                       //seting flat phase space, required for registration of FORTRAN channel on empty slot

    Tauolapp::RegisterChannel( -1, demo_modify );    //registration of new channel using demo_modify
    demo_modify->print();                            //prints channel info


    //*********************************************************************************************
    //Example of presampler modification, case of 3 scalars
    //*********************************************************************************************
    ChannelForTauola *demo_presampler = GetChannel(74); // TAU -> pi- pi0 gamma channel
    // Probabilities of two generation channels: prob1, prob2. Third one uses 1-prob1-prob2.
    // Channel presampler for 3-scalar invariant mass uses amrx, gamrx.
    // For 2-scalar invariant masses respectively amra, gamra and amrb, gamrb are used.
    // Name of parameters as used in F77: prob1, prob2, amrx, gamrx,  amra, gamra,  amrb, gamrb
    SetPresampler3( demo_presampler,        0.0,   1.0, 1.12,  0.37, 0.777,  0.15, 0.789,  0.01);



    //*********************************************************************************************
    // Example of modification, case of multiple scalars.
    // This case differs from all others. Sigee is a one-parameter function regardless
    // of multiplicity of the decay channel. NOTE: flat phase space in multiscalar system.
    // Important: if number of products is reduced below 6 those channels  
    // should be registered stating exact channel number. Available numbers have to be 
    // in range NM4+NM5+1 to NM4+NM5+NM6.  See ../TAUDCDsize.inc for numerical value of NM4,NM5,NM6
    //*********************************************************************************************

    ChannelForTauola *test1 = GetChannel(54);
    test1->setBr(9.0);
    test1->setName("TAU -> K+ K+ pi- pi- pi-");
 
    products = test1->getProducts();
    products[3]= 3;  // K+   we replace two of the pions with K's 
    products[4]= 3;  // K+
    int mult=5; // multiplicity of channel
    test1->setMultiplicity(mult);

 
    /*
    // altenatively we can set channel of 7 products
    test1->setName("TAU -> pi+ pi+ pi- pi- pi- pi0 pi0");
    products[3]= 1;  // pi+
    products[4]= 1;  // pi+
    products[5]= 2;  // pi0
    products[6]= 2;  // pi0
    test1->setProducts(products);
    mult=7;
    test1->setMultiplicity(mult);
    */

    Tauolapp::RegisterChannel( 54, test1 );

    ChannelForTauola *sigee_test_channel = new ChannelForTauola(9.99, products,"TAU sigee test", test_sigee, mult);

    Tauolapp::RegisterChannel( 55, sigee_test_channel);
    if( sigee_test_channel->getChannelNo() == 0 ) printf("Couldn't register new channel!\n");
    else{       printf("\nRegistered new channel. It is now: \n");
                Tauolapp::PrintChannelInfo(sigee_test_channel->getChannelNo());
    }

    printf("\n         =============================================================== \n");
    printf("         =  INFO FOR ALL NON-LEPTONIC  CHANNELS (initialized or not):  = \n");
    printf("         =============================================================== \n\n");
    for(int i=1; i<=NMODE; ++i) {
        Tauolapp::PrintChannelInfo(i);
    }
}


// Set pointer for user channel redefinition
// This function is called from taumain.f
extern "C" void tauolaredef_() {

    // uncomment to run examples for re-initialization from C++:
    Tauolapp::SetUserRedefinitions(RedefExample);      // basic examples 
    //Tauolapp::SetUserRedefinitions(TestCommunication); //more technical examples
}
