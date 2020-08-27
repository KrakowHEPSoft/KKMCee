#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "ChannelForTauolaInterface.h"

namespace Tauolapp {

ChannelForTauola* taubra_userChannels[500] = { NULL };
void*             leptonChannelsMEpointers[2] = { NULL };
void            (*channelRedefinitionFunction)() = NULL;

// Hidden flag that blocks redefinition outside of 'void iniofc_()' (file: channel_wrappers.c)
bool isInChannelRedefinitionFunction = false;

void SetUserRedefinitions(void (*function)()) {
    channelRedefinitionFunction = function;
}

int GetChannelNumber(const int Mgroup, const int subchannel) {
    if( subchannel<=0 || Mgroup<=0 ) return 0;

    int channel = subchannel;

    if( Mgroup != 4)              channel += NM4;
    if( Mgroup != 4 && Mgroup!=5) channel += NM5;
    if( Mgroup <  4)              channel += NM6;
    if( Mgroup <  3)              channel += NM3;
    if( Mgroup <  2)              channel += NM2;

    return channel;
}

void GetMultipAndSubchannelNumbers(const int channel, int &mulpik, int &subchannel, int &Mgroup) {
    if( channel<=0 || channel > NMODE ) {
        printf("Tauolapp::GetMultipAndSubchannelNumbers: channel out of range [1,%i]: %i\n",NMODE,channel);

        mulpik = subchannel = 0;
    }

    mulpik = taudcd_.MULPIK[channel-1];  // we get multip multi-p(pions).
                                         // it is multiplicity minus one:  
                                         // for SM there is always tau neutrino
                     // At present, only when Mgroup= 6, Mgroup may not be  equal mulpik. 
                     // taudcd_.MULPIK Originally it means multiplicity of pions and kaons for the channel [channel-1].
                     // At present of all products which are not tau neutrino. Complication arise 
                     // because of new neutrinoless channels. We copy MULPIK to mulpik for the further use.

    if     (channel<=NM4)                     { Mgroup=4; subchannel = channel; }
    else if(channel<=NM4+NM5)                 { Mgroup=5; subchannel = channel - (NM4); }
    else if(channel<=NM4+NM5+NM6)             { Mgroup=6; subchannel = channel - (NM4+NM5); }
    else if(channel<=NM4+NM5+NM6+NM3)         { Mgroup=3; subchannel = channel - (NM4+NM5+NM6); }
    else if(channel<=NM4+NM5+NM6+NM3+NM2)     { Mgroup=2; subchannel = channel - (NM4+NM5+NM6+NM3); }
    else if(channel<=NM4+NM5+NM6+NM3+NM2+NM1) { Mgroup=1; subchannel = channel - (NM4+NM5+NM6+NM3+NM2); }
    else { printf("Tauolapp::GetMultipAndSubchannelNumbers: channel out of range [1,%i]: %i\n",NMODE,channel); exit(-1);}

}

void PrintChannelInfo(int channel) {
    if( channel<=0 || channel > NMODE ) {
        printf("Tauolapp::PrintChannelInfo: channel out of range [1,%i]: %i\n",NMODE,channel);
        return;
    }

    ChannelForTauola *pointer = GetChannel(channel);
    
    pointer->print();

    delete pointer;
}

int OverwriteChannel(int channel, ChannelForTauola *pointer) {
    if( taubra_userChannels[channel] == NULL ) {
        printf("Tauolapp::OverwriteChannel: INFO: channel %i was not registered.\n",channel);
    }
    else {
        delete taubra_userChannels[channel];
        taubra_userChannels[channel] = NULL;
    }

    return RegisterChannel(channel,pointer);
}

int RegisterChannel(int channel, ChannelForTauola *pointer) {
    if(!isInChannelRedefinitionFunction) {
        printf("Tauolapp::ERROR: channel redefinition/overwriting not permitted\n");
        printf("          outside of user redefinition function. Use:\n");
        printf("          Tauolapp::SetUserRedefinitions( void (*pointer)() )\n");
        printf("          to register function that redefines Tauola channels.\n");
        exit(-1);
    }

    if( !pointer ){ 
        printf("Tauolapp::RegisterChannel: ERROR: NULL pointer provided");
        exit(-1);
    }

    int mulpik = pointer->getMultiplicity();

    if( channel == -1 ) {
        // Find appropriate, empty channel
        
        int *key  = NULL;
        int range = 0;
        int subchannel = 0;
        switch(mulpik) {  // Here Mgroup is be defined. At present Mgroup=6 should be set for mulpik>5 as default.
                           // If we want to have Mgroup=6 for smaller mulpik we have to overwrite the
                          // existing Mgroup=6 channel. This is the case which is not expected to be often in use.
            case 1:  key = metyp_.KEY1; range = NM1; break;
            case 2:  key = metyp_.KEY2; range = NM2; break;
            case 3:  key = metyp_.KEY3; range = NM3; break; 
            case 4:  key = metyp_.KEY4; range = NM4; break; 
            case 5:  key = metyp_.KEY5; range = NM5; break; 
            default: key = metyp_.KEY6; range = NM6; break; // 6 or more
        }
        
        for(int i=0; i<range; ++i) {
            if( key[i] == 0 ) {
                subchannel = i+1;
                break;
            }
        }

        // subchannel out of range
        if( subchannel == 0 ) {
            printf("Tauolapp::RegisterChannel: no free space for channel in Mgroup for multip-1= %i\n",mulpik);
            return 0;
        }

        channel = GetChannelNumber(mulpik,subchannel);
        
        if( channel <= 0 ) return 0;
    }
    
    if( channel<=0 || channel > NMODE ) {
        printf("Tauolapp::RegisterChannel: channel out of range [1,%i]: %i\n",NMODE,channel);
        return 0;
    }
    
    int old_mulpik = 0; 
    int subchannel = 0;
    int Mgroup     = 0;

    // old_mulpik needed for the cross-check of multiplicity 
    GetMultipAndSubchannelNumbers(channel,old_mulpik,subchannel, Mgroup);


    if( Mgroup != 6 && old_mulpik != mulpik ) {
        printf("Tauolapp::RegisterChannel: ERROR: channel multiplicity-1 mismatch old %i supposed new %i\n",old_mulpik,mulpik);
        exit(-1);
    }

    if( taubra_userChannels[channel] != NULL ) {
        printf("Tauolapp::RegisterChannel: the %i channel was already (re)defined in this run.\n",channel);
        printf("                           If it is really necessary, use OverwriteChannel.\n");
        exit(-1);
    }

    if( channel == 1 || channel == 2 ) {
        printf("Tauolapp::RegisterChannel WARNING: order of particles in channel %2i\n",channel);
        printf("          does not match order of arguments for currents and/or matrix\n");
        printf("          element calculation. See README in tauola-c directory for details.\n");
    }
    
    //
    // Set channel info
    //

    string      name = pointer->getName();
    int         me   = pointer->getMeType();
    vector<int> ki   = pointer->getProducts();

    // pointer
    taubra_userChannels[channel] = pointer;

    // name
    memset (taudcd_.NAMES[channel-1],' ',31);
    strncpy(taudcd_.NAMES[channel-1],name.c_str(),name.length());

    // branching ratio
    taubra_.GAMPRT[channel+NLT-1] = pointer->getBr();

    // multiplicity-1 (at present only for group with multiplicity 6+)
    if( Mgroup == 6 ) taudcd_.MULPIK[channel-1] = mulpik;

    // ME type
    switch(Mgroup) {
        case 1:  metyp_.KEY1[subchannel-1] = me; break;
        case 2:  metyp_.KEY2[subchannel-1] = me; break;
        case 3:  metyp_.KEY3[subchannel-1] = me; break; 
        case 4:  metyp_.KEY4[subchannel-1] = me; break;
        case 5:  metyp_.KEY5[subchannel-1] = me; break;
        default: metyp_.KEY6[subchannel-1] = me; break;
    }

    // decay products
    for(unsigned int i=0; i<9; ++i) {
        if( i<ki.size() ) taudcd_.IDFFIN[channel-1][i] = ki[i];
        else              taudcd_.IDFFIN[channel-1][i] = 0;
    }
    
    pointer->setChannelNo(channel);
    pointer->setSubchannelNo(subchannel);

    return channel;
}

int ModifyLeptonic(int channel, int me, float br, void *pointer ) {
    if( channel<=0 || channel > 2 ) {
        printf("Tauolapp::ModifyLeptonic: channel out of range: %i\n",channel);
        return 0;
    }

    if( me > 3 ) {
        leptonChannelsMEpointers[channel-1] = pointer;
        metyp_.KEY0[channel-1]    = me;
        taubra_.GAMPRT[channel-1] = br;
        return channel;
    }
    else {
        printf("Tauolapp::ModifyLeptonic failed. me=  %i\n",me);
        return 0;
    }
}

ChannelForTauola* GetChannel(int channel) {
    if( channel <= 0 || channel > NMODE ) {
        printf("Tauolapp::GetChannel: channel out of range [1,%i]: %i\n",NMODE,channel);
        return NULL;
    }
    
    int    mulpik     = 0;
    int    Mgroup     = 0;
    int    subchannel = 0;
    int    me_type    = 0;
    double br         = taubra_.GAMPRT[channel+NLT-1];
    
    /** @note FORTRAN strings do not end with '\0' delimiter.
              We have to explicitly write its range */
    string name(taudcd_.NAMES[channel-1],taudcd_.NAMES[channel-1]+31);
    vector<int> ki(9);

    // get multip, subchannel
    GetMultipAndSubchannelNumbers(channel,mulpik,subchannel, Mgroup);

    // get me type
    switch(Mgroup) { 
        case 1:  me_type = metyp_.KEY1[subchannel-1]; break;
        case 2:  me_type = metyp_.KEY2[subchannel-1]; break;
        case 3:  me_type = metyp_.KEY3[subchannel-1]; break; 
        case 4:  me_type = metyp_.KEY4[subchannel-1]; break; 
        case 5:  me_type = metyp_.KEY5[subchannel-1]; break; 
        default: me_type = metyp_.KEY6[subchannel-1]; break; // 6 or more
    }

    // get decay products
    for(int i=0;i<9;++i) ki[i] = taudcd_.IDFFIN[channel-1][i];

    // create ChannelForTauola
    ChannelForTauola *pointer = new ChannelForTauola(mulpik,me_type,br,ki,name);
    
    pointer->setChannelNo(channel);
    pointer->setSubchannelNo(subchannel);

    if( me_type > 3 ) {
        pointer->setFunctionPointer(taubra_userChannels[channel]);
    }

    return pointer;
}

int SetPresampler2(ChannelForTauola *pointer, float prob1, float prob2, float am2, float gam2, float am3, float gam3) {
    if( !pointer ) {
        printf("Tauolapp::SetPresampler2: NULL pointer provided.");
        return -1;
    }
    
    int subchannel = pointer->getSubchannelNo();

    if( pointer->getMultiplicity() != 2 || subchannel<1 || subchannel>NM2) {  
        printf("Tauolapp::SetPresampler2: wrong SetPresampler method used for channel %i.",pointer->getChannelNo());
        return -1;
    }

    printf("SetPresampler2: values for channel %i before:  prob1=%9.6f prob2=%9.6f am2=%9.6f gam2=%9.6f am3=%9.6f gam3=%9.6f\n",
            pointer->getChannelNo(),
            sampl2_.PROB1[subchannel-1],sampl2_.PROB2[subchannel-1],
            sampl2_.AM2[subchannel-1],  sampl2_.GAM2[subchannel-1],
            sampl2_.AM3[subchannel-1],  sampl2_.GAM3[subchannel-1]);

    sampl2_.PROB1[subchannel-1] = prob1;
    sampl2_.PROB2[subchannel-1] = prob2;
    sampl2_.AM2  [subchannel-1] = am2;
    sampl2_.GAM2 [subchannel-1] = gam2;
    sampl2_.AM3  [subchannel-1] = am3;
    sampl2_.GAM3 [subchannel-1] = gam3;

    printf("SetPresampler2: values for channel %i changed: prob1=%9.6f prob2=%9.6f am2=%9.6f gam2=%9.6f am3=%9.6f gam3=%9.6f\n",
            pointer->getChannelNo(),prob1,prob2,am2,gam2,am3,gam3);
    return 0;
}

int SetPresampler3(ChannelForTauola *pointer, float prob1, float prob2, float amrx, float gamrx, float amra, float gamra, float amrb, float gamrb) {
    if( !pointer ) {
        printf("Tauolapp::SetPresampler3: NULL pointer provided.");
        return -1;
    }
    
    int subchannel = pointer->getSubchannelNo();

    if( pointer->getMultiplicity() != 3 || subchannel<1 || subchannel>NM3) {  
        printf("Tauolapp::SetPresampler3: wrong SetPresampler method used for channel %i.",pointer->getChannelNo());
        return -1;
    }

    printf("SetPresampler3: values for channel %i before:  prob1=%9.6f prob2=%9.6f amrx=%9.6f gamrx=%9.6f amra=%9.6f gamra=%9.6f amrb=%9.6f gamrb=%9.6f\n",
            pointer->getChannelNo(),
            sampl3_.PROB1[subchannel-1],sampl3_.PROB2[subchannel-1],
            sampl3_.AMRX[subchannel-1], sampl3_.GAMRX[subchannel-1],
            sampl3_.AMRA[subchannel-1], sampl3_.GAMRA[subchannel-1],
            sampl3_.AMRB[subchannel-1], sampl3_.GAMRB[subchannel-1]);

    sampl3_.PROB1[subchannel-1] = prob1;
    sampl3_.PROB2[subchannel-1] = prob2;
    sampl3_.AMRX [subchannel-1] = amrx;
    sampl3_.GAMRX[subchannel-1] = gamrx;
    sampl3_.AMRA [subchannel-1] = amra;
    sampl3_.GAMRA[subchannel-1] = gamra;
    sampl3_.AMRB [subchannel-1] = amrb;
    sampl3_.GAMRB[subchannel-1] = gamrb;

    printf("SetPresampler3: values for channel %i changed. prob1=%9.6f prob2=%9.6f amrx=%9.6f gamrx=%9.6f amra=%9.6f gamra=%9.6f amrb=%9.6f gamrb=%9.6f\n",
            pointer->getChannelNo(),prob1,prob2,amrx,gamrx,amra,gamra,amrb,gamrb);

    return 0;
}

int SetPresampler4(ChannelForTauola *pointer, float prob1, float prob2, float amrx, float gamrx, float amra, float gamra) {
    if( !pointer ) {
        printf("Tauolapp::SetPresampler4: NULL pointer provided.");
        return -1;
    }
    
    int subchannel = pointer->getSubchannelNo();

    if( pointer->getMultiplicity() != 4 || subchannel<1 || subchannel>NM4) {  
        printf("Tauolapp::SetPresampler4: wrong SetPresampler method used for channel %i.",pointer->getChannelNo());
        return -1;
    }

    printf("SetPresampler4: values for channel %i before:  prob1=%9.6f prob2=%9.6f amrx=%9.6f gamrx=%9.6f amra=%9.6f gamra=%9.6f\n",
            pointer->getChannelNo(),
            sampl4_.PROB1[subchannel-1],sampl4_.PROB2[subchannel-1],
            sampl4_.AMRX[subchannel-1], sampl4_.GAMRX[subchannel-1],
            sampl4_.AMRA[subchannel-1], sampl4_.GAMRA[subchannel-1]);

    sampl4_.PROB1[subchannel-1] = prob1;
    sampl4_.PROB2[subchannel-1] = prob2;
    sampl4_.AMRX [subchannel-1] = amrx;
    sampl4_.GAMRX[subchannel-1] = gamrx;
    sampl4_.AMRA [subchannel-1] = amra;
    sampl4_.GAMRA[subchannel-1] = gamra;

    printf("SetPresampler4: values for channel %i changed. prob1=%9.6f prob2=%9.6f amrx=%9.6f gamrx=%9.6f amra=%9.6f gamra=%9.6f\n",
            pointer->getChannelNo(),prob1,prob2,amrx,gamrx,amra,gamra);

    return 0;
}

int SetPresampler5(ChannelForTauola *pointer, float proba2, float probom, float ama2, float gama2, float amom, float gamom) {
    if( !pointer ) {
        printf("Tauolapp::SetPresampler5: NULL pointer provided.");
        return -1;
    }
    
    int subchannel = pointer->getSubchannelNo();

    if( pointer->getMultiplicity() != 5 || subchannel<1 || subchannel>NM5) {  
        printf("Tauolapp::SetPresampler5: wrong SetPresampler method used for channel %i.",pointer->getChannelNo());
        return -1;
    }

    printf("SetPresampler5: values for channel %i before:  proba2=%9.6f probom=%9.6f ama2=%9.6f gama2=%9.6f amom=%9.6f gamom=%9.6f\n",
            pointer->getChannelNo(),
            sampl5_.PROBa2[subchannel-1],sampl5_.PROBOM[subchannel-1],
            sampl5_.ama2[subchannel-1],  sampl5_.gama2[subchannel-1],
            sampl5_.AMOM[subchannel-1],  sampl5_.GAMOM[subchannel-1]);

    sampl5_.PROBa2[subchannel-1] = proba2;
    sampl5_.PROBOM[subchannel-1] = probom;
    sampl5_.ama2  [subchannel-1] = ama2;
    sampl5_.gama2 [subchannel-1] = gama2;
    sampl5_.AMOM  [subchannel-1] = amom;
    sampl5_.GAMOM [subchannel-1] = gamom;

    printf("SetPresampler5: values for channel %i changed. proba2=%9.6f probom=%9.6f ama2=%9.6f gama2=%9.6f amom=%9.6f gamom=%9.6f\n",
            pointer->getChannelNo(),proba2,probom,ama2,gama2,amom,gamom);
    return 0;
}

} // namespace Tauolapp
