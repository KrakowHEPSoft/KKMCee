#ifndef  TAUOLAPP_CHANNELFORTAUOLA_H
#define  TAUOLAPP_CHANNELFORTAUOLA_H
/**
 *  @file ChannelForTauola.h
 *  @brief Definition and implementation of \b class ChannelForTauola
 *
 *  @class Tauolapp::ChannelForTauola
 *  @brief Interface for adding new decay channels and modifying existing ones
 *
 *  This class ensures that SIMPLE ASPECTS for user channel definition is correct and complete.
 *  Pointer to the function provided by the user must match one of the
 *  types defined at the top of this header file. This type is checked
 *  by the constructors.
 *
 *  This by no means help to verify if user channel definition is correct and complete.
 *  For that task one has to inspect physics sense of each parameter of the  ChannelForTauola method
 *  and of relation between function provided by the user, clalculating hadronic current  or matrix element calculation.
 *  To understand conventions see our pipi0.c  example, in case of doubts Fortran code for the routines of names used 
 *  in pointer type names can be helpful.
 */
#include <vector>
#include <string>
#include <cstdio>
#include <cstring>
#include <complex>
using std::string;
using std::vector;

namespace Tauolapp {

//
// Function pointer types
//
typedef void  (*DAM1PI_POINTER_TYPE)(const float*, const float,  const float*,  const float,  float&, float*);                                           // PNU AMF0 PKK AMF1 GAMM HV
typedef void  (*DAM2PI_POINTER_TYPE)(const float*, const float*, const float*,  const float*, float&, float*);                                           // PT PN PIM1 PIM2                AMPLIT HV
typedef void  (*DAM3PI_POINTER_TYPE)(const float*, const float*, const float*,  const float*, const float*, float&, float*);                             // PT PN PIM1 PIM2 PIM3           AMPLIT HV
typedef void  (*DAM4PI_POINTER_TYPE)(const float*, const float*, const float*,  const float*, const float*, const float*, float&, float*);               // PT PN PIM1 PIM2 PIM3 PIM4      AMPLIT HV
typedef void  (*DAM5PI_POINTER_TYPE)(const float*, const float*, const float*,  const float*, const float*, const float*, const float*, float&, float*); // PT PN PIM1 PIM2 PIM3 PIM4 PIM5 AMPLIT HV
typedef void  (*DAMPRY_POINTER_TYPE)(const int,    const double, const double*, const double*, const double*, const double*, double&, double*);          // ITDKRC XK0DEC XK XA QP XN      AMPLIT HV
typedef void  (*CURR2_POINTER_TYPE) (const float*, const float*, std::complex<float>*);                                                                  // PIM1 PIM2                HADCUR
typedef void  (*CURR3_POINTER_TYPE) (const float*, const float*, const float*,std:: complex<float>*);                                                    // PIM1 PIM2 PIM3           HADCUR
typedef void  (*CURR4_POINTER_TYPE) (const float*, const float*, const float*, const float*, std::complex<float>*);                                      // PIM1 PIM2 PIM3 PIM4      HADCUR
typedef void  (*CURR5_POINTER_TYPE) (const float*, const float*, const float*, const float*, const float*, std::complex<float>*);                        // PIM1 PIM2 PIM3 PIM4 PIM5 HADCUR
typedef float (*SIGEE_POINTER_TYPE) (const float); // AMX2
typedef float (*FCONST_POINTER_TYPE)();

class ChannelForTauola {
//
// Constructors
//
public:
    /** @brief Default constructor
     *
     *  Provides all necessary information about the channel.
     *  Channel and sub-channel numbers are used to define position of the
     *  channel on Tauola channel list. It will be attributed when channel
     *  is registered in Tauola, unless this info is needed
     *  to put the modified channel back in place.
     *
     *  @note This constructor is also used by constructors
     *        providing pointer to user function
     */
    ChannelForTauola(int mulpik, int me_type, float br, const vector<int> &ki, string name = string("UNNAMED")) { init(mulpik,me_type,br,ki,name); }

    ChannelForTauola(float br, const vector<int> &ki, string name, SIGEE_POINTER_TYPE  pointer, int mulpik) { init(mulpik,4,br,ki,name); m_pointer = (void*)pointer; }
    ChannelForTauola(float br, const vector<int> &ki, string name, FCONST_POINTER_TYPE pointer) { init(1,5,br,ki,name); m_pointer = (void*)pointer; }
    ChannelForTauola(float br, const vector<int> &ki, string name, CURR2_POINTER_TYPE  pointer) { init(2,4,br,ki,name); m_pointer = (void*)pointer; }
    ChannelForTauola(float br, const vector<int> &ki, string name, CURR3_POINTER_TYPE  pointer) { init(3,4,br,ki,name); m_pointer = (void*)pointer; }
    ChannelForTauola(float br, const vector<int> &ki, string name, CURR4_POINTER_TYPE  pointer) { init(4,4,br,ki,name); m_pointer = (void*)pointer; }
    ChannelForTauola(float br, const vector<int> &ki, string name, CURR5_POINTER_TYPE  pointer) { init(5,4,br,ki,name); m_pointer = (void*)pointer; }
    ChannelForTauola(float br, const vector<int> &ki, string name, DAM1PI_POINTER_TYPE pointer) { init(1,5,br,ki,name); m_pointer = (void*)pointer; }
    ChannelForTauola(float br, const vector<int> &ki, string name, DAM2PI_POINTER_TYPE pointer) { init(2,5,br,ki,name); m_pointer = (void*)pointer; }
    ChannelForTauola(float br, const vector<int> &ki, string name, DAM3PI_POINTER_TYPE pointer) { init(3,5,br,ki,name); m_pointer = (void*)pointer; }
    ChannelForTauola(float br, const vector<int> &ki, string name, DAM4PI_POINTER_TYPE pointer) { init(4,5,br,ki,name); m_pointer = (void*)pointer; }
    ChannelForTauola(float br, const vector<int> &ki, string name, DAM5PI_POINTER_TYPE pointer) { init(5,5,br,ki,name); m_pointer = (void*)pointer; }

//
// Functions
//
public:
    /** @brief Print information about the channel */
    void print() const;

private:
    /** @brief initializes all fields */
    void init(int mulpik, int me_type, float br, const vector<int> &ki, string name) {
        m_channel_no    = 0;
        m_subchannel_no = 0;
        m_pointer       = NULL;
        m_mulpik        = mulpik;
        m_me_type       = me_type;

        setBr(br);
        setProducts(ki);
        setName(name);
    }
//
// Accessors
//
public:
    int  getChannelNo()       const { return m_channel_no;     } //!< Get     channel position on Tauola channel list
    void setChannelNo(int channel)  { m_channel_no = channel;  } //!< Set     channel position on Tauola channel list
    int  getSubchannelNo()    const { return m_subchannel_no;  } //!< Get sub-channel position on Tauola sub-channel list
    void setSubchannelNo(int subch) { m_subchannel_no = subch; } //!< Set sub-channel position on Tauola sub-channel list

    int  getMultiplicity()     const { return m_mulpik;         } //!< Get multiplicity
    void setMultiplicity(int new_m)  { m_mulpik = new_m;        } //!< Set multiplicity

    int  getMeType()           const { return m_me_type;        } //!< Get ME type
    void setMeType(int metype);                                   //!< Set ME type

    double getBr()            const { return m_br;             } //!< Get branching ratio
    void   setBr(const double br);                               //!< Set branching ratio

    const string& getName()   const { return m_name;            } //!< Get channel name
    void          setName(const string &name);                    //!< Set channel name

    const vector<int>& getProducts() const { return m_products; } //!< Get list of decay products
    void               setProducts(const vector<int>& products ); //!< Set list of decay products

    void* getFunctionPointer() const { return m_pointer;        } //!< Get pointer for function of matrix element/current calculation
    void  setFunctionPointer(void* pointer);                      //!< Set pointer for function of matrix element/current calculation

//
// Fields
//
private:
    /** @brief   m_channel_no m_subchannel_no store position on the lists of decay channels of Tauola and are not used in the class
     *  These numbers are used only if channel of  Tauola defaults  is to be   modified by the user.
     *  Verification can be then performed whether it set back to its original plase in list of Tauola decay channels
     */

    int m_channel_no;       //!< variable in which number on Tauola list of decay channels for selected multiplicity can be stored
    int m_subchannel_no;    //!< variable in which  number on Tauola  sub-list of decay channels for selected multiplicity

    //

    int m_mulpik;           //!< Multiplicity
    int m_me_type;          //!< Matrix element type

    double m_br;            //!< Branching ratio
    string m_name;          //!< Name of the decay channel

    vector<int> m_products; //!< List of decay products
    void *m_pointer;        //!< Function pointer
};

//
// Implementation of print function and setters
//
inline void ChannelForTauola::print() const {

    printf("Channel %3i: '%-31s' b.ratio: %6.5f, multipl-1: %i, subch: %2i ME-type: %i Final state:",
            m_channel_no,
            m_name.c_str(),
            m_br,
            m_mulpik,
            m_subchannel_no,
            m_me_type);

    for(int i=0; i<9; ++i) {
        int precision = 6;
        if( i > 3 ) precision = 3;

        if( i == m_mulpik ) printf("|");
        else                printf(" ");

        if( i < (int)m_products.size() ) printf("%*i",precision,m_products[i]);
        else                             printf("%*i",precision,0);
    }
    
    if( m_me_type > 3 ) {
        if( m_pointer == NULL ) {
            printf("\n             ERROR: me type is %i and pointer to user function not set!",m_me_type);
            exit(-1);
        }
        else printf("\n            Pointer to user function set to: %p",m_pointer);
    }

    printf("\n");
}

inline void ChannelForTauola::setMeType(int metype) {
    if( metype < 0 || metype > 3 ) {
      printf("ChannelForTauola::setMeType: For this method only ME types 0,1,2,3 are allowed.");
        return;
    }

    m_me_type = metype;
}

inline void ChannelForTauola::setBr(double br) {
    if( br < 0.0 ) {
        printf("ChannelForTauola::setBr: branching ratio must be >=0. Set to 0.");
        br = 0.0;
    }

    m_br = br;
}

inline void ChannelForTauola::setName(const string &name) {
    if( name.size() > 31 ) {
        printf("ChannelForTauola::setName: channel name will be truncated to 31 characters: '%s'.\n",name.c_str());
    }

    m_name = name.substr(0,31);
}

inline void ChannelForTauola::setProducts(const vector<int> &products) {

    m_products = products;

    if( m_products.size() > 9 ) {
        printf("ChannelForTauola::setProducts: only first 9 decay products will be used.\n");
        m_products.resize(9);
    }
}

inline void ChannelForTauola::setFunctionPointer(void *pointer) {
    if( pointer == NULL ) {
        printf("ChannelForTauola::setFunctionPointer: NULL pointer provided. Ignored.\n");
        return;
    }

    if( m_pointer != NULL ) {
        printf("ChannelForTauola::setPointer: pointer already set. New pointer ignored.\n");
        return;
    }

    m_pointer = pointer;
}

} // namespace Tauolapp

#endif
