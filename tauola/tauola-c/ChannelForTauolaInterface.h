#ifndef TAUOLAPP_CHANNELFORTAUOLAINTERFACE_H
#define TAUOLAPP_CHANNELFORTAUOLAINTERFACE_H

#include "ChannelForTauola.h"
#include "TauolaStructs.h"
#include <vector>
#include <string>
using std::vector;
using std::string;

namespace Tauolapp {

/** @brief Supplement to struct taubra_ defined in TauolaStructs.h
           Holds pointers to user defined ChannelForTauola objects */
extern ChannelForTauola* taubra_userChannels[500];

/** @brief Pointers to ME for lepton channels */
extern void*             leptonChannelsMEpointers[2];

/** @brief Pointer to tauola channel redefinition function */
extern void            (*channelRedefinitionFunction)();

/** @brief Print current definition of a channel */
void PrintChannelInfo(int channel);

/** @brief Set tauola channel redefinition function pointer
 *
 *  Sets pointer to function reinitializing Tauola channels.
 *  @note This is the only place where Tauola channels can be reinitialized.
 */
void SetUserRedefinitions(void (*function)());

/** @brief Overwrite channel
 *
 *  Removes previous information about the channel and overwrites
 *  with new one.
 */
int OverwriteChannel(int channel, ChannelForTauola *pointer);


/** @brief Register channel
 *
 *  Changes information about selected channel or adds new channel
 *
 *  @param  channel Channel to be changed. -1 means new channel has to be added
 *  @param  pointer Pointer to class defining new or modified channel
 *  @return added or modified channel number or 0 in case of failure
 *
 *  @note Same channel number cannot be registered twice. Use OverwriteChannel
 *        instead.
 */
int RegisterChannel(int channel, ChannelForTauola *pointer);


/** @brief Modify leptonic channel
 *
 *  Changes information about selected leptonic channel
 *  This functionality is by far more delicate than RegisterChannel
 *  User is requested to perform detailed tests
 *
 *  @param  channel Channel number (1 or 2) to be changed
 *  @param  me      Matrix element type. See README for the meaning
 *  @param  br      Branching ratio of the selected channel
 *  @param  current Pointer to function calculating the xsection or current
            @warning pointer has to be of type DAMPRY_POINTER_TYPE
                     defined in ChannelForTauola.h
 *  @return modified channel number or 0 in case of failure
 */
int ModifyLeptonic(int channel, int me, float br, void *xsec = NULL);


/** @brief Get channel information
 *
 *  Get channel information from FORTRAN common block
 *  Returns ChannelForTauola object with channel information filled
 *  ready to register without any changes or to update and register.
 *  Note: Such ChannelForTauola object can be registered only back on its original place,
 *        unless matrix element calculation is reduced to constant (flat phase-space)
 *        or ponters to the user methods are introduced, me_type != 2, is required.
 */
ChannelForTauola* GetChannel(int channel);


/** @brief Set parameters of presampler for 2 scalar mode
 *
 *  Set parameters used to optimize efficiency of 2-scalar mode phase space generator.
 *  First argument is a pointer for ChannelForTauola for which presampler parameters will be changed
 *  All following arguments are new numerical values of parameters.
 *
 *  @note 0<=prob1; 0<=prob2; prob1+prob2<=1
 */
int SetPresampler2(ChannelForTauola *pointer, float prob1, float prob2, float am2, float gam2, float am3, float gam3);


/** @brief Set parameters of presampler for 3 scalar mode
 *
 *  Set parameters used to optimize efficiency of 3-scalar mode phase space generator.
 *  First argument is a pointer for ChannelForTauola for which presampler parameters will be changed
 *  All following arguments are new numerical values of parameters.
 *
 *  @note 0<=prob1; 0<=prob2; prob1+prob2<=1
 */
int SetPresampler3(ChannelForTauola *pointer, float prob1, float prob2, float amrx, float gamrx, float amra, float gamra, float amrb, float gamrb);


/** @brief Set parameters of presampler for 4 scalar mode
 *
 *  Set parameters used to optimize efficiency of 4-scalar mode phase space generator.
 *  First argument is a pointer for ChannelForTauola for which presampler parameters will be changed
 *  All following arguments are new numerical values of parameters.
 *
 *  @note 0<=prob1; 0<=prob2; prob1+prob2<=1; only prob1+prob2 is used 
 */
int SetPresampler4(ChannelForTauola *pointer, float prob1, float prob2, float amrx, float gamrx, float amra, float gamra);


/** @brief Set parameters of presampler for 5 scalar mode
 *
 *  Set parameters used to optimize efficiency of 5-scalar mode phase space generator.
 *  First argument is a pointer for ChannelForTauola for which presampler parameters will be changed
 *  All following arguments are new numerical values of parameters.
 *
 *  @note 0<=proba2<=1; 0<=probom<=1 
 */
int SetPresampler5(ChannelForTauola *pointer, float proba2, float probom, float ama2, float gama2, float amom, float gamom);

} // namespace Tauolapp

#endif
