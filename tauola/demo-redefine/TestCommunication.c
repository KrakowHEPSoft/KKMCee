#include <cstdlib>
#include <cstdio>
#include <complex>
#include "../tauola-c/ChannelForTauolaInterface.h"
#include "../tauola-c/ChannelForTauola.h"
using std::complex;
using namespace Tauolapp;

void test_dam1pi(const float *pnu, const float amf0, const float *pkk, const float amf1, float &gamm, float *hv) {
    static bool printed = false;
    if(!printed) { printf("test_dam1pi turned on\n"); printed = true; }
    hv[0] = 1.1;
    hv[1] = 2.2;
    hv[2] = 3.3;
    hv[3] = 4.4;
    gamm  = 5.5;
}

void  test_dam2pi(const float *pt, const float *pn, const float *pim1, const float *pim2, float &amplit, float *hv) {
    static bool printed = false;
    if(!printed) { printf("test_dam2pi turned on\n"); printed = true; }
    hv[0] = 1.1;
    hv[1] = 2.2;
    hv[2] = 3.3;
    hv[3] = 4.4;
    amplit= 5.5;
}

void  test_dam3pi(const float *pt, const float *pn, const float *pim1, const float *pim2, const float *pim3, float &amplit, float *hv) {
    static bool printed = false;
    if(!printed) { printf("test_dam3pi turned on\n"); printed = true; }
    hv[0] = 1.1;
    hv[1] = 2.2;
    hv[2] = 3.3;
    hv[3] = 4.4;
    amplit= 5.5;
}

void  test_dam4pi(const float *pt, const float *pn, const float *pim1, const float *pim2, const float *pim3, const float *pim4, float &amplit, float *hv) {
    static bool printed = false;
    if(!printed) { printf("test_dam4pi turned on\n"); printed = true; }
    hv[0] = 1.1;
    hv[1] = 2.2;
    hv[2] = 3.3;
    hv[3] = 4.4;
    amplit= 5.5;
}

void  test_dam5pi(const float *pt, const float *pn, const float *pim1, const float *pim2, const float *pim3, const float *pim4, const float *pim5, float &amplit, float *hv) {
    static bool printed = false;
    if(!printed) { printf("test_dam5pi turned on\n"); printed = true; }
    hv[0] = 1.1;
    hv[1] = 2.2;
    hv[2] = 3.3;
    hv[3] = 4.4;
    amplit= 5.5;
}

void test_dampry(const int itdkrc, const double xk0dec, const double *xk, const double *xa, const double *qp, const double *xn, double &amplit, double *hv) {
    static bool printed = false;
    if(!printed) { printf("test_dampry turned on\n"); printed = true; }
    hv[0]  = 1.1;
    hv[1]  = 2.2;
    hv[2]  = 3.3;
    hv[3]  = 4.4;
    amplit = 5.5;
}

void test_curr2(const float *pim1, const float *pim2, complex<float> *hadcur) {
    static bool printed = false;
    if(!printed) { printf("test_curr2 turned on\n"); printed = true; }
    hadcur[0]  = complex<float>(1.,5.);
    hadcur[1]  = complex<float>(2.,6.);
    hadcur[2]  = complex<float>(3.,7.);
    hadcur[3]  = complex<float>(4.,8.);
}

void  test_curr3pi(const float *pim1, const float *pim2, const float *pim3, complex<float> *hadcur) {
    static bool printed = false;
    if(!printed) { printf("test_curr3pi turned on\n"); printed = true; }
    hadcur[0]  = complex<float>(1.,5.);
    hadcur[1]  = complex<float>(2.,6.);
    hadcur[2]  = complex<float>(3.,7.);
    hadcur[3]  = complex<float>(4.,8.);
}

void  test_curr4(const float *pim1, const float *pim2, const float *pim3, const float *pim4, complex<float> *hadcur) {
    static bool printed = false;
    if(!printed) { printf("test_curr4 turned on\n"); printed = true; }
    hadcur[0]  = complex<float>(1.,5.);
    hadcur[1]  = complex<float>(2.,6.);
    hadcur[2]  = complex<float>(3.,7.);
    hadcur[3]  = complex<float>(4.,8.);
}

void  test_curr5(const float *pim1, const float *pim2, const float *pim3, const float *pim4, const float *pim5, complex<float> *hadcur) {
    static bool printed = false;
    if(!printed) { printf("test_curr5 turned on\n"); printed = true; }
    hadcur[0]  = complex<float>(1.,5.);
    hadcur[1]  = complex<float>(2.,6.);
    hadcur[2]  = complex<float>(3.,7.);
    hadcur[3]  = complex<float>(4.,8.);
}

float test_const() {
    static bool printed = false;
    if(!printed) { printf("test_const turned on\n"); printed = true; }
    return 1.0;
}

float test_sigee(const float q2) {
    static bool printed = false;
    if(!printed) { printf("test_sigee turned on\n"); printed = true; }
    return 1.35E-03*(1+ 0.01*q2);
}

void TestCommunication() {
    vector<int> products(9);



    // =====================
    // two body channels,     brief parameter explanation at  ``three body channels``
    // =====================

    products[0] = -1;
    RegisterChannel(193, new ChannelForTauola(9.9991, products, "test_dam1pi 193", test_dam1pi) );
    RegisterChannel( -1, new ChannelForTauola(9.9991, products, "test_dam1pi  -1", test_dam1pi) );


    // =====================
    // three body channels
    // =====================

    products[0] = -1;
    products[1] =  2;
    // to replace initialization for channel 148  with the new one  of new hadronic current pointer test_curr2
    RegisterChannel(148, new ChannelForTauola(9.9991, products, "test_curr2  148", test_curr2)  );

    // to add new channel of matrix element: pointer   test_dam2pi
    RegisterChannel( -1, new ChannelForTauola(9.9991, products, "test_dam2pi  -1", test_dam2pi) );


    // MIND DIFFERENT PARAMETER LIST, in this, and next,  case we use flat phase space, 
    // no pointer to current/matrix element previously used to decipher the channel multiplicity !

    // Register new channel no. matrix element, flat phase space
    // NOTE: In this example, tau neutrino is replaced with muon.
    //       its PDGID identifier is at products[2] which is normally not used by 3-body channel
    //       When replacing tau neutrino, use products[2] for charged particle PDGID
    //       and products[3] for neutral particles (see next example)
    //       Also for 2 body decay modes products[2] can be used to replace tau-neutrino
    products[2] =  13; 
    products[3] =  0;
    RegisterChannel( -1, new ChannelForTauola(2, 1, 9.7, products, "FLAT PHASE-SPACE  -1" ) );

    // Modify matrix element for channel 152 to flat phase space
    // NOTE: Tau neutrino is replaced with photon. As mentioned above,
    //       its PDGID identifier is at products[3] which is normally not used by 3-body channel
    //       Also for 2 body decay modes products[2] can be used to replace tau-neutrino
    products[2] =  0;
    products[3] =  22;
    RegisterChannel(152, new ChannelForTauola(2, 1, 9.6, products, "FLAT PHASE-SPACE 152" ) );

 



    // =====================
    // four body channels,   brief parameter explanation at  ``three body channels``
    // =====================

    products[0] = -1;
    products[1] =  2;
    products[2] =  2;
    products[3] =  0;

    RegisterChannel( 80, new ChannelForTauola(9.9991, products, "test_curr3pi 80", test_curr3pi));
    RegisterChannel( -1, new ChannelForTauola(9.9991, products, "test_curr3pi -1", test_curr3pi));
    RegisterChannel( 81, new ChannelForTauola(9.9991, products, "test_dam3pi  81", test_dam3pi) );

    // NOTE: there is only 1 free space for 3-pi channels,
    //       so 2nd use of '-1' will fail. We check for that
    ChannelForTauola *test_dam3pi_minus1 = new ChannelForTauola(9.9991, products, "test_dam3pi  -1", test_dam3pi);
    RegisterChannel( -1, test_dam3pi_minus1 );

    // =====================
    // special test
    // =====================
    if( test_dam3pi_minus1->getChannelNo() == 0 ) {
        printf("ChannelForTauola test: test_dam3pi -1 not registered as expected (no free space on a list of subchannels for multp=3)\n");
    }
    else {
        printf("ChannelForTauola test: surprisingly, test_dam3pi -1 registered. Check for possible errors.\n");
        exit(-1);
    }

    // =====================
    // five body channels,   brief parameter explanation at  ``three body channels``
    // =====================
    products[0] = -1;
    products[1] =  2;
    products[2] =  2;
    products[3] =  2;
    RegisterChannel( 20, new ChannelForTauola(9.9991, products, "test_dam4pi  20", test_dam4pi) );
    RegisterChannel( -1, new ChannelForTauola(9.9991, products, "test_dam4pi  -1", test_dam4pi) );
    RegisterChannel( 21, new ChannelForTauola(9.9991, products, "test_curr4   21", test_curr4)  );
    RegisterChannel( -1, new ChannelForTauola(9.9991, products, "test_curr4   -1", test_curr4)  );

    // =====================
    // six body channels,    brief parameter explanation at  ``three body channels``
    // =====================
    products[0] = -1;
    products[1] =  2;
    products[2] =  2;
    products[3] =  2;
    products[4] =  2;
    RegisterChannel( 35, new ChannelForTauola(9.9991, products, "test_dam5pi  35", test_dam5pi) );
    RegisterChannel( -1, new ChannelForTauola(9.9991, products, "test_dam5pi  -1", test_dam5pi) );
    RegisterChannel( 36, new ChannelForTauola(9.9991, products, "test_curr5   36", test_curr5)  );
    RegisterChannel( -1, new ChannelForTauola(9.9991, products, "test_curr5   -1", test_curr5)  );

    // =====================
    // seven- to nine- body channels, brief parameter explanation at  ``three body channels``
    // =====================
    products[0] = -1;
    products[1] =  2;
    products[2] =  2;
    products[3] =  2;
    products[4] =  2;
    products[5] =  2;
    RegisterChannel( 55, new ChannelForTauola(9.9991, products, "test_sigee   55   7-body", test_sigee, 6)  );

    products[6] =  2;
    RegisterChannel( 57, new ChannelForTauola(9.9991, products, "test_sigee   57   8-body", test_sigee, 7)  );

    products[7] =  2;
    RegisterChannel( 59, new ChannelForTauola(9.9991, products, "test_sigee   59   9-body", test_sigee, 8)  );

    // =====================
    // leptonic  channel
    // =====================

    ModifyLeptonic( 2, 5, 9.9999, (void*)test_dampry );  // do not forget consequences of itdkrc which must be on argument list 
                                                         // for (void*)test_dampry
}
