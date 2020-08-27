//################################################################################
//
// Prepared by: 
// M.Chrzaszcz(UZH, IFJ), Z.Was (IFJ)
// with theory from S.Turczyk ( arXiv:0812.3830 )
//
// Names: LLLL: 4 leptons left-handed 
// LLRR: 2 left-handed, 2 right-handed leptons
// Rad: radiative operator
//
// The purpose of these examples is to demonstrate how the program can be used
// Same models were used in LHCb for model dependece studies for tau->3mu
// All models can be used also in Belle(2)
//#################################################################################

#include "iostream"

using namespace std;

static const double mumass = 0.105752;
static const double taumass = 1.77641;
static const double PI = 3.14159265359;
static const double emass = 0.000511;

// Dalitz density function for tau->3mu

double gammaVLLLL(const double m12sq, const double m23sq, const double m13sq) ;
double gammaVLLRR(const double m12sq, const double m23sq, const double m13sq) ;
double gammaRAD(const double m12sq, const double m23sq, const double m13sq);
double interferenceVLLLL_RAD(const double m12sq, const double m23sq, const double m13sq);
double interferenceVLLRR_RAD(const double m12sq, const double m23sq, const double m13sq);
///###################################################

// Dalitz density function for tau- -> mu- mu+ e-

double gammaVLLLL_mumueOS(const double m12sq, const double m23sq, const double m13sq) ;
double gammaVLLRR_mumueOS(const double m12sq, const double m23sq, const double m13sq) ;
double gammaRAD_mumueOS(const double m12sq, const double m23sq, const double m13sq) ;
double interferenceVLLLL_RAD_mumueOS(const double m12sq, const double m23sq, const double m13sq);
double interferenceVLLRR_RAD_mumueOS(const double m12sq, const double m23sq, const double m13sq);
//#####################################################

// Dalitz density function for tau- -> mu- mu- e+  

double gammaVLLLL_mumueSS(const double m12sq, const double m23sq, const double m13sq) ;
double gammaVLLRR_mumueSS(const double m12sq, const double m23sq, const double m13sq) ;
//#####################################################

// Dalitz density function for tau- -> e- e+ mu-

double gammaVLLLL_mueeOS(const double m12sq, const double m23sq, const double m13sq) ;
double gammaVLLRR_mueeOS(const double m12sq, const double m23sq, const double m13sq) ;
double gammaRAD_mueeOS(const double m12sq, const double m23sq, const double m13sq) ;
double interferenceVLLLL_RAD_mueOS(const double m12sq, const double m23sq, const double m13sq);
double interferenceVLLRR_RAD_mueOS(const double m12sq, const double m23sq, const double m13sq);

//#####################################################  

// Dalitz density function for tau- -> e- e- mu+ 

double gammaVLLLL_mueeSS(const double m12sq, const double m23sq, const double m13sq) ;
double gammaVLLRR_mueeSS(const double m12sq, const double m23sq, const double m13sq) ;
//#####################################################

// Dalitz density function for tau- -> 3e

double gammaVLLLL_3e(const double m12sq, const double m23sq, const double m13sq) ;
double gammaVLLRR_3e(const double m12sq, const double m23sq, const double m13sq) ;
double gammaRAD_3e(const double m12sq, const double m23sq, const double m13sq);
double interferenceVLLLL_RAD_3e(const double m12sq, const double m23sq, const double m13sq);
double interferenceVLLRR_RAD_3e(const double m12sq, const double m23sq, const double m13sq);




double mass_sq(const float *pim1, const float *pim2);


void lfv_dam2pi(const float *pt, const float *pn, const float *pim1, const float *pim2, float &amplit, float *hv) {

  //########################################################
  //  For tau-  -> mu- mu- mu+   pm1 pm2  pn                 
  // declaration -13 -13  13
  //######################################################## 

  
  
   double m12sq=mass_sq(pim1,pim2);
   double m13sq=mass_sq(pim1,pn);       
   double m23sq=mass_sq(pim2,pn);     
   
   //   amplit =  gammaVLLLL(m12sq,m23sq,m13sq);        
   //   amplit = gammaVLLRR(m12sq,m23sq,m13sq);    
      amplit = gammaRAD(m12sq,m23sq,m13sq);
   //   amplit = interferenceVLLLL_RAD(m12sq,m23sq,m13sq); 
   //   amplit = interferenceVLLRR_RAD(m12sq,m23sq,m13sq);
  
   //########################################################
   //  For tau-  -> e- mu - mu+   pm1 pm2 pn
   // in tauola the ordering is:  mu- e- mu+ (WRONG)
   //########################################################
   /*
   double m12sq=mass_sq(pim2,pim1);
   double m13sq=mass_sq(pim2,pn);
   double m23sq=mass_sq(pim1,pn);
   
   //  amplit =  gammaVLLLL_mumueOS(m12sq,m23sq,m13sq); 
   //  amplit =  gammaVLLRR_mumueOS(m12sq,m23sq,m13sq);
   //  amplit =  gammaRAD_mumueOS(m12sq,m23sq,m13sq);
   //  amplit = interferenceVLLLL_RAD_mumueOS(m12sq,m23sq,m13sq); 
   amplit = interferenceVLLRR_RAD_mumueOS(m12sq,m23sq,m13sq);
   */
   
   /*
     cout<<"Electron: "<<(double)(pim1[3]*pim1[3]-pim1[2]*pim1[2]-pim1[1]*pim1[1]-pim1[0]*pim1[0])<<endl;
     cout<<"Electron: "<<(double)(pim2[3]*pim2[3]-pim2[2]*pim2[2]-pim2[1]*pim2[1]-pim2[0]*pim2[0])<<endl;
     cout<<"Electron: "<<(double)(pn[3]*pn[3]-pn[2]*pn[2]-pn[1]*pn[1]-pn[0]*pn[0])<<endl<<endl;
   */
   
   //########################################################      
   //  For tau-  -> e+ mu- mu-   pm1 pm2  pn                        
   // in tauola the ordering is:  mu- mu- e+ (WRONG)               
   //########################################################      
  // email from sascha (16.07.2014):  
  //Here m12sq = m12^2 is the invariant mass squared of the two muons, and
  //m23sq = m23^2 is the invariant mass of one of the muons and the electron.
  
  //##############################################################
  /*
  double m12sq=mass_sq(pim1,pim2); 
  double m13sq=mass_sq(pim1,pn);   
  double m23sq=mass_sq(pim2,pn);   

  //  amplit = gammaVLLLL_mumueSS(m12sq,m23sq,m13sq);
  amplit = gammaVLLRR_mumueSS(m12sq,m23sq,m13sq); 
  
  */



  hv[0]  = 1.;
  hv[1]  = 0.;
  hv[2]  = 0.;
  hv[3]  = 0.;
}

// REALISTIC EXAMPLE 
// implem,entation of  EFT aproach for tau->3mu
// see 0707.0988 for theory details
//##################################################################################
//   Tau- >3 mu modelts:
//##################################################################################
// Purely left-handed operator:
double gammaVLLLL(const double m12sq, const double m23sq, const double m13sq) {
  /// notation of 0707.0988                                                                                       
                                                                                                                                  
  double denom = taumass;
  denom = denom*denom*denom;
  double num = m12sq*(taumass*taumass + (mumass*mumass - m12sq  )) -  2* (mumass*mumass * (taumass*taumass + mumass*mumass) );
  return num;

}

// Mixed operator:
double gammaVLLRR(const double m12sq, const double m23sq, const double m13sq) {
  ///  notation of 0707.0988  
  
  double denom = taumass;
  denom = denom*denom*denom;
  double num = -24.*mumass*mumass*(taumass*taumass+mumass*mumass - m12sq/3.) - 4.*(m13sq*m13sq+m23sq*m23sq) + 4.*(taumass*taumass+3.*mumass*mumass) * (m13sq+m23sq) ;
  return num;
}

// Radiative operator:
double gammaRAD(const double m12sq,const  double m23sq, const double m13sq) { 

  double m12sq_MeV=m12sq*1.e6;
  double m23sq_MeV=m23sq*1.e6;
  double m13sq_MeV=m13sq*1.e6;
  
  double denom = taumass*taumass*taumass*(1.e9);
  double a = taumass*taumass*(1e6);
  double b = mumass*mumass*(1e6);
  double c = m13sq_MeV;
  double d = m23sq_MeV;
  double e = m12sq_MeV;
  denom = denom*c*d*c*d;

  double num =
    -6*b*c*c*d*d - 6*b*b*c*d*(c + d) + 2*b*b*b*(c + d)*(c + d) +
    2*a*a*b*(c*c + c*d + d*d) -
    2*a*b*(3*c*d*(c + d) + b*(2*c*c + 3*c*d + 2*d*d)) +
    c*d*(c*c*c + c*c*d + c*d*d + d*d*d + 4*c*d*e + (c + d)*e*e);

  return num/denom;

}

//interference between LLLL and RAD:
double interferenceVLLLL_RAD(const double m12sq, const double m23sq, const double m13sq) {
  /// notation of 0707.0988                                                   
  double m12sq_MeV=m12sq*1.e6;
  double m23sq_MeV=m23sq*1.e6;
  double m13sq_MeV=m13sq*1.e6;

  double denom = taumass*taumass*taumass*(1.e9);
  double a = taumass*taumass*(1e6);
  double b = mumass*mumass*(1e6);
  double c = m13sq_MeV;
  double d = m23sq_MeV;
  double e = m12sq_MeV;
  denom = denom*c*d;
  double num = a*b*(c+d) - b*(6*c*d+b*(c+d)) + 2*c*d*e;
  return num/denom;

}
//interference between LLRR and RAD:              
double interferenceVLLRR_RAD(const double m12sq, const double m23sq, const double m13sq) {
  /// notation of 0707.0988 
  double m12sq_MeV=m12sq*1.e6;
  double m23sq_MeV=m23sq*1.e6;
  double m13sq_MeV=m13sq*1.e6;

  double denom = taumass*taumass*taumass*(1.e9);
  double a = taumass*taumass*(1e6);
  double b = mumass*mumass*(1e6);
  double c = m13sq_MeV;
  double d = m23sq_MeV;
  double e = m12sq_MeV;
  denom = denom*c*d;
  double num = a*b*c + a*(b+c)*d - b*(3*c*d+b*(c+d)) - c*d*e;
  return num/denom;
}

//##########################################################################################    
//##########################################################################################
// TAU--> e- mu- mu+ 
// function for mass square of two particles
//##########################################################################################

// Purely left-handed operator: 
double gammaVLLLL_mumueOS(const double m12sq, const double m23sq, const double m13sq) {
  /// notation of 0707.0988                                                                                                            

  double num = ( (taumass*taumass*taumass*taumass) - (2.*m12sq-taumass*taumass-2.*mumass*mumass)*(2.*m12sq-taumass*taumass-2.*mumass*mumass))/(512.*3.1415*taumass*taumass*taumass);
  return num;

}
// Mixed operator:
double gammaVLLRR_mumueOS(const double m12sq, const double m23sq, const double m13sq) {
  ///  notation of 0707.0988     
  //ATTENTION!!!!!!!
  // In this framework couplig to e an mu are different, here we assume:
  double g_mu=0.;
  double g_e=1.;
  
  double num = g_e* ( (taumass*taumass*taumass*taumass) - (2.*m13sq-taumass*taumass-2.*mumass*mumass)*(2.*m13sq-taumass*taumass-2.*mumass*mumass))/(512.*3.1415*3.1415*3.1415*taumass*taumass*taumass);

  num+= g_mu* ( (taumass*taumass - 2.*mumass*mumass)*(taumass*taumass - 2.*mumass*mumass) - (2.*m23sq-taumass*taumass-2.*mumass*mumass)*(2.*m23sq-taumass*taumass-2.*mumass*mumass))/(512.*3.1415*3.1415*3.1415*taumass*taumass*taumass);
  
  return num;

}
// Radiative operator:
double gammaRAD_mumueOS(const double m12sq, const double m23sq, const double m13sq) {

  double num = ( mumass*mumass*(m23sq-taumass*taumass)*(m23sq-taumass*taumass)  )/(64.*PI*PI*PI*taumass*taumass*taumass*m23sq*m23sq) + ( m12sq*m12sq + m13sq*m13sq - 2.*mumass*mumass*mumass*mumass )/(128.*PI*PI*PI*taumass*taumass*taumass*m23sq) + (taumass*taumass-m23sq)/(128.*PI*PI*PI*taumass*taumass*taumass);
  return num;

}
//interference between LLLL and RAD: 
double interferenceVLLLL_RAD_mumueOS(const double m12sq, const double m23sq, const double m13sq){

  double num = (m12sq-2.*mumass*mumass)/(128.*PI*PI*PI*taumass*taumass) + (mumass*mumass)/(128.*PI*PI*PI*m23sq);  
  return num;

}
//interference between LLRR and RAD:
double interferenceVLLRR_RAD_mumueOS(const double m12sq, const double m23sq, const double m13sq){

  double num = (m13sq-2.*mumass*mumass)/(128.*PI*PI*PI*taumass*taumass) + (mumass*mumass)/(128.*PI*PI*PI*m23sq);
  return num;
  
}
//##########################################################################################
//##########################################################################################
// TAU--> e+ mu- mu-                         
// function for mass square of two particles 
//##########################################################################################

// Purely left-handed operator:  
double gammaVLLLL_mumueSS(const double m12sq, const double m23sq, const double m13sq){
  double num = ((m12sq - 2*mumass*mumass)*(-m12sq + emass*emass + taumass*taumass))/(64.*taumass*taumass*taumass*PI*PI*PI);
  return num;
}

// Mixed operator:
double gammaVLLRR_mumueSS(const double m12sq, const double m23sq, const double m13sq){
  double num = -(m12sq*m12sq - m12sq*(-2*m23sq + emass*emass + 4.*mumass*mumass + taumass*taumass) +
		 2.*(m23sq*m23sq + mumass*mumass*mumass*mumass + 2.*mumass*mumass*taumass*taumass + emass*emass*(2*mumass*mumass + taumass*taumass) -
		     m23sq*(emass*emass + 2*mumass*mumass + taumass*taumass)))/(256.*taumass*taumass*taumass*PI*PI*PI);
    return num;
}


//##########################################################################################
//##########################################################################################
// TAU--> e- e+  mu-                          
// function for mass square of two particles  
//########################################################################################## 

// Purely left-handed operator:
double gammaVLLLL_mueeOS(const double m12sq, const double m23sq, const double m13sq) {
  /// notation of 0707.0988                                                                                                                                   

  double num = ( (taumass*taumass*taumass*taumass) - (2.*m12sq-taumass*taumass-2.*emass*emass)*(2.*m12sq-taumass*taumass-2.*emass*emass))/(512.*3.1415*taumass*taumass*taumass);
  return num;

}
// Mixed operator:
double gammaVLLRR_mueeOS(const double m12sq, const double m23sq, const double m13sq) {
  ///  notation of 0707.0988                                                                                                                                  
  //ATTENTION!!!!!!!                                                                                                                                          
  // In this framework couplig to e an mu are different, here we assume:                                                                                      
  double g_mu=0.;
  double g_e=1.;

  double num = g_e* ( (taumass*taumass*taumass*taumass) - (2.*m13sq-taumass*taumass-2.*emass*emass)*(2.*m13sq-taumass*taumass-2.*emass*emass))/(512.*3.1415*3.1415*3.1415*taumass*taumass*taumass);
  num+= g_mu* ( (taumass*taumass - 2.*emass*emass)*(taumass*taumass - 2.*emass*emass) - (2.*m23sq-taumass*taumass-2.*emass*emass)*(2.*m23sq-taumass*taumass-2.*emass*emass))/(512.*3.1415*3.1415*3.1415*taumass*taumass*taumass);

  return num;

}
// Radiative operator: 
double gammaRAD_mueeOS(const double m12sq, const double m23sq, const double m13sq) {

  double num = ( emass*emass*(m23sq-taumass*taumass)*(m23sq-taumass*taumass)  )/(64.*PI*PI*PI*taumass*taumass*taumass*m23sq*m23sq) + ( m12sq*m12sq + m13sq*m13sq - 2.*emass*emass*emass*emass )/(128.*PI*PI*PI*taumass*taumass*taumass*m23sq) + (taumass*taumass-m23sq)/(128.*PI*PI*PI*taumass*taumass*taumass);
  return num;
}
//interference between LLLL and RAD: 
double interferenceVLLLL_RAD_mueeOS(const double m12sq, const double m23sq, const double m13sq){

  double num = (m12sq-2.*emass*emass)/(128.*PI*PI*PI*taumass*taumass) + (emass*emass)/(128.*PI*PI*PI*m23sq);
  return num;
}
//interference between LLRR and RAD: 
double interferenceVLLRR_RAD_mueeOS(const double m12sq, const double m23sq, const double m13sq){

  double num = (m13sq-2.*emass*emass)/(128.*PI*PI*PI*taumass*taumass) + (emass*emass)/(128.*PI*PI*PI*m23sq);
  return num;

}

//##########################################################################################                                                                
//##########################################################################################                                                                
// TAU--> e+ e- mu-                                                                                                                                        
// function for mass square of two particles                                                                                                                
//##########################################################################################

// Purely left-handed operator: 
double gammaVLLLL_mueeSS(const double m12sq, const double m23sq, const double m13sq){
  double num = ((m12sq - 2*emass*emass)*(-m12sq + mumass*mumass + taumass*taumass))/(64.*taumass*taumass*taumass*PI*PI*PI);
  return num;
}

double gammaVLLRR_mueeSS(const double m12sq, const double m23sq, const double m13sq){
  double num = -(m12sq*m12sq - m12sq*(-2*m23sq + mumass*mumass + 4.*emass*emass + taumass*taumass) +
                 2.*(m23sq*m23sq + emass*emass*emass*emass + 2.*emass*emass*taumass*taumass + mumass*mumass*(2*emass*emass + taumass*taumass) -
                     m23sq*(mumass*mumass + 2*emass*emass + taumass*taumass)))/(256.*taumass*taumass*taumass*PI*PI*PI);
  return num;
}


//##########################################################################################    
//##########################################################################################    
// TAU--> e+ e- e-                                                                             
// function for mass square of two particles                                                    
//########################################################################################## 

// Purely left-handed operator: 
double gammaVLLLL_3e(const double m12sq, const double m23sq, const double m13sq) {
  /// notation of 0707.0988                                                                                                                                          

  double denom = taumass;
  denom = denom*denom*denom;
  double num = m12sq*(taumass*taumass + (emass*emass - m12sq  )) -  2* (emass*emass * (taumass*taumass + emass*emass) );
  return num;

}

// Mixed operator:     
double gammaVLLRR_3e(const double m12sq, const double m23sq, const double m13sq) {
  ///  notation of 0707.0988                                                                                                                                         

  double denom = taumass;
  denom = denom*denom*denom;
  double num = -24.*emass*emass*(taumass*taumass+emass*emass - m12sq/3.) - 4.*(m13sq*m13sq+m23sq*m23sq) + 4.*(taumass*taumass+3.*emass*emass) * (m13sq+m23sq) ;
  return num;
}

// Radiative operator:                                                                                                                                               
double gammaRAD_3e(const double m12sq,const  double m23sq, const double m13sq) {

  double m12sq_MeV=m12sq*1.e6;
  double m23sq_MeV=m23sq*1.e6;
  double m13sq_MeV=m13sq*1.e6;

  double denom = taumass*taumass*taumass*(1.e9);
  double a = taumass*taumass*(1e6);
  double b = emass*emass*(1e6);
  double c = m13sq_MeV;
  double d = m23sq_MeV;
  double e = m12sq_MeV;
  denom = denom*c*d*c*d;

  double num =
    -6*b*c*c*d*d - 6*b*b*c*d*(c + d) + 2*b*b*b*(c + d)*(c + d) +
    2*a*a*b*(c*c + c*d + d*d) -
    2*a*b*(3*c*d*(c + d) + b*(2*c*c + 3*c*d + 2*d*d)) +
    c*d*(c*c*c + c*c*d + c*d*d + d*d*d + 4*c*d*e + (c + d)*e*e);

  return num/denom;

}

//interference between LLLL and RAD:
double interferenceVLLLL_RAD_3e(const double m12sq, const double m23sq, const double m13sq) {
  /// notation of 0707.0988                                                                    
  double m12sq_MeV=m12sq*1.e6;
  double m23sq_MeV=m23sq*1.e6;
  double m13sq_MeV=m13sq*1.e6;

  double denom = taumass*taumass*taumass*(1.e9);
  double a = taumass*taumass*(1e6);
  double b = emass*emass*(1e6);
  double c = m13sq_MeV;
  double d = m23sq_MeV;
  double e = m12sq_MeV;
  denom = denom*c*d;
  double num = a*b*(c+d) - b*(6*c*d+b*(c+d)) + 2*c*d*e;
  return num/denom;

}

//interference between LLRR and RAD:                                                           
double interferenceVLLRR_RAD_3e(const double m12sq, const double m23sq, const double m13sq) {
  /// notation of 0707.0988                                                                    
  double m12sq_MeV=m12sq*1.e6;
  double m23sq_MeV=m23sq*1.e6;
  double m13sq_MeV=m13sq*1.e6;

  double denom = taumass*taumass*taumass*(1.e9);
  double a = taumass*taumass*(1e6);
  double b = emass*emass*(1e6);
  double c = m13sq_MeV;
  double d = m23sq_MeV;
  double e = m12sq_MeV;
  denom = denom*c*d;
  double num = a*b*c + a*(b+c)*d - b*(3*c*d+b*(c+d)) - c*d*e;
  return num/denom;
}

double mass_sq(const float *pim1, const float *pim2)
{
  double tmp[4];
  tmp[0]=(double)((double)pim1[0]+(double)pim2[0]);
  tmp[1]=(double)((double)pim1[1]+(double)pim2[1]);
  tmp[2]=(double)((double)pim1[2]+(double)pim2[2]);
  tmp[3]=(double)((double)pim1[3]+(double)pim2[3]);

  return (double)(tmp[3]*tmp[3]-tmp[2]*tmp[2]-tmp[1]*tmp[1]-tmp[0]*tmp[0]);
}



