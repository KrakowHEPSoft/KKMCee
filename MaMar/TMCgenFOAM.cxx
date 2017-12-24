#include "TMCgenFOAM.h"
#include <stdio.h>

////////////////////////////////////////////////////////////////////////////////
///     TMCgenFOAM class
/// This is class for axiliary exercises,  mainly integration with Monte Carlo

ClassImp(TMCgenFOAM);

TMCgenFOAM::TMCgenFOAM():
  TMCgen()
{
  /// This constructor is for ROOT streamers ONLY
  cout<< "----> TMCgenFOAM Default Constructor (for ROOT only) "<<endl;
  m_Foam3 = NULL;
  m_Foam2 = NULL;
  m_Foam1 = NULL;
}

///______________________________________________________________________________________
TMCgenFOAM::~TMCgenFOAM()
{
  //!Explicit destructor
  cout<< "----> TMCgenFOAM::TMCgenFOAM !!!! DESTRUCTOR !!!! "<<endl;
}///destructor

///_____________________________________________________________
TMCgenFOAM::TMCgenFOAM(const char* Name):
  TMCgen(Name)
{
//! all defaults defined here can be changed by the user
//! before calling TMCgen::Initialize
  m_Foam3 = NULL;
  m_Foam2 = NULL;
  m_Foam1 = NULL;
  m_IsFoam5 = 1;   // Foam5 ON
  m_IsFoam3 = 1;   // Foam3 ON
  m_IsFoam2 = 0;   // Foam2 OFF
  m_IsFoam1 = 0;   // Foam1 OFF
///////////////////////////////////////////////////
// Physics
  m_gnanob  = 389.37966e3;
  m_pi      = 3.1415926535897932;
  m_ceuler  = 0.57721566;
  m_alfinv  = 137.035989;
  m_alfpi   = 1/m_alfinv/m_pi;
  m_amel    = 0.510999e-3;
//
  m_beam    = 0.510999e-3;  // electron
  m_chini   =-1.0;          // electron
//
//  m_fin     = 0.105;      // WRONG by 0.7 MeV !!!! costed 1 week to find it out !!!
  m_fin     = 0.1056583;    // final ferm. muon mass
  m_chfin   =-1.0;          // final ferm. muon charge
//
  m_KFini   = 11;           // electron
  m_KFf     = 13;           // muon
//
  m_KeyISR  = 2;            // Type of ISR/QED switch, 0,1,2
  m_jmax    = 10000;        // length of xpar
///////////////////////////////////////////////////
/// Foam setup
  m_kDim    =    5;         // No. of dim. for Foam, =2,3 Machine energy spread OFF/ON
  m_nCells  = 2000;         // No. of cells, optional, default=2000
  m_nSampl  =  200;         // No. of MC evts/cell in exploration, default=200
//  m_del     = 1e-4;         // limit for |gamma*ln(eps)| in IFI mapping
  m_del     = 1e-6;         // limit for |gamma*ln(eps)| in IFI mapping $$$
//  m_eps = 1e-6;             // IR regulator
  m_eps = 1e-8;             // IR regulator, test $$$
  m_vvcut = 0.020;          // auxiliary vmax for soft limit test
  m_Mode    = 5;
///////////////////////////////////////////////////
// debug
  m_count   =0;

cout<< "----> TMCgenFOAM::TMCgenFOAM USER Constructor "<<endl;
}///

///______________________________________________________________________________________
void TMCgenFOAM::Initialize(TRandom *RNgen, ofstream *OutFile, TH1D* h_NORMA)
{
  cout<< "----> TMCgenFOAM::Initialize, Entering "<<endl;
  ///	      SETTING UP RANDOM NUMBER GENERATOR
  TMCgen::Initialize(  RNgen, OutFile, h_NORMA);

/////////////////////////////////////////////////////////
//  m_NevGen=0;
  const int jmax =m_jmax;
  ReaData("../../.KK2f_defaults", jmax, m_xpar);  // numbering as in input!!!
  ReaData("./pro.input",         -jmax, m_xpar);  // jmax<0 means no-zeroing
  double ypar[jmax];
  for(int j=0;j<jmax;j++) ypar[j]=m_xpar[j+1];    // ypar has c++ numbering
  //
  //NevTot = (long)m_xpar[0];                       // NevTot hidden in xpar[0] !!!
  m_CMSene  = m_xpar[ 1];
  m_vvmax   = m_xpar[17];

  cout<<" TMCgen::Initialize: m_CMSene="<<m_CMSene<<endl;
  cout<<" TMCgen::Initialize: m_vvmax="<<m_vvmax<<endl;

  const char *output_file = "./kkmc.output";
  long stl2 = strlen(output_file);
  int mout = 16;
  kk2f_fort_open_(mout,output_file,stl2);
  kk2f_initialize_(ypar);
  kksem_initialize_(ypar);

  double errel;
  /////////////////////////////////////////////////////////
  if(f_IsInitialized == 0)
  {
  /// ******  SETTING UP FOAM of base class  *****
  if( m_IsFoam5 == 1 ){
    f_FoamI   = new TFOAM("FoamI");   // new instance of MC generator FOAM
    m_kDim    = 5;
    m_nCells  =  10000;
    m_nSampl  = 100000;
    f_FoamI->SetkDim(m_kDim);         // No. of dims. Obligatory!
    f_FoamI->SetnCells(m_nCells);     // No. of cells, optional, default=2000
    f_FoamI->SetnSampl(m_nSampl);     // No. of MC evts/cell in exploration, default=200
    f_FoamI->SetnBin(        16);     // No. of bins default 8
    f_FoamI->SetOptRej(0);            // wted events (=0), default wt=1 events (=1)
    m_Mode    = 5;
    f_FoamI->Initialize( f_RNgen, this);     // Initialize FOAM
    f_FoamI->GetIntNorm(m_Xnorm,errel);   // universal normalization
  }// m_IsFoam5
  //////////////////////////////////////////////////////////////
  /// ******  SETTING UP additional FOAM of the user class *****
  if( m_IsFoam3 == 1 ){
    m_Foam3   = new TFOAM("Foam3");   // new instance of MC generator FOAM
    m_kDim    = 3;
    m_Foam3->SetkDim(m_kDim);         // No. of dims. Obligatory!
    m_Foam3->SetnCells(m_nCells);     // No. of cells, optional, default=2000
    m_Foam3->SetnSampl(m_nSampl);     // No. of MC evts/cell in exploration, default=200
    m_Foam3->SetnBin(        16);     // No. of bins default 8
    m_Foam3->SetOptRej(0);            // wted events (=0), default wt=1 events (=1)
    m_Mode    = 3;
    m_count =0;
    m_Foam3->Initialize( f_RNgen, this);     // Initialize FOAM
    m_Foam3->GetIntNorm(m_Xsav3,errel);
  }// m_IsFoam3
  //////////////////////////////////////////////////////////////
  /// ******  SETTING UP additional FOAM of the user class *****
  if( m_IsFoam1 == 1 ){
  //[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
  // checking EW implementation in Born xsection
  // SUBROUTINE KKsem_Afb_Calc(KeyDist,KFi,KFf,CMSene,vv,Result)
    double xBorn,xBorn1;
    kksem_ord1v_(  1,m_KFini, m_KFf, m_CMSene, 0e0, xBorn);  // Born [nb]
    kksem_ord1v_(501,m_KFini, m_KFf, m_CMSene, 0e0, xBorn1);  // Born [nb]
    cout<< "|||| xBorn   Gmu scheme = "<< xBorn  << endl;
    cout<< "|||| xBorn alpha scheme = "<< xBorn1 << "   "<< xBorn1/xBorn <<endl;
    double AfbBorn;
    kksem_afb_calc_(  1,m_KFini, m_KFf, m_CMSene, 0e0, AfbBorn);  // AFB
    cout<< "**** KKsem_Afb_Calc AFB = "<< AfbBorn <<endl;
    double svar = sqr(m_CMSene);
    bornv_interpogsw_(m_KFf,svar, 0.0);
    double dSig_EEX0 = bornv_dizet_( 1, m_KFini, m_KFf, svar,  0.0, 0.0, 0.0, 0.0, 0.0);
    double sig0nb = 4*m_pi* sqr(1/m_alfinv)/(3.0*svar )*m_gnanob;
    //double sig_EEX =  dSig_EEX0   *3.0/8.0 *sig0nb;  // Born of EEX
    double sig_EEX =  dSig_EEX0   *sig0nb;  // Born of EEX
    cout<< "|||| Born Dizet   = "<< sig_EEX << "   "<< sig_EEX/xBorn <<endl;
    //exit(-5);
    //]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
    m_Foam1   = new TFOAM("Foam1");   // new instance of MC generator FOAM
    m_kDim    = 2;
    m_Foam1->SetkDim(m_kDim);         // No. of dims. Obligatory!
    m_Foam1->SetnCells(m_nCells);     // No. of cells, optional, default=2000
    m_Foam1->SetnSampl(m_nSampl);     // No. of MC evts/cell in exploration, default=200
    m_Foam1->SetnBin(        16);     // No. of bins default 8
    m_Foam1->SetOptRej(0);            // wted events (=0), default wt=1 events (=1)
    m_Mode    = 1;
    m_count =0;
    m_Foam1->Initialize( f_RNgen, this);     // Initialize FOAM
    m_Foam1->GetIntNorm(m_Xsav1,errel);
  }// m_IsFoam1
  ////////////////////////////////
  if( m_IsFoam2 == 1 ){
    m_Foam2   = new TFOAM("Foam2");   // new instance of MC generator FOAM
    m_kDim    = 1;
    m_Foam2->SetkDim(m_kDim);         // No. of dims. Obligatory!
    m_Foam2->SetnCells(m_nCells);     // No. of cells, optional, default=2000
    m_Foam2->SetnSampl(m_nSampl);     // No. of MC evts/cell in exploration, default=200
    m_Foam2->SetnBin(        16);     // No. of bins default 8
    m_Foam2->SetOptRej(0);            // wted events (=0), default wt=1 events (=1)
    m_Mode    = 2;
    m_count =0;
    m_Foam2->Initialize( f_RNgen, this);     // Initialize FOAM
    m_Foam2->GetIntNorm(m_Xsav2,errel);
  }// m_IsFoam2
 //////////////////////////////////////////////////////////////
  //screen output
  BXOPE(*f_Out);
  BXTXT(*f_Out,"========================================");
  BXTXT(*f_Out,"======   TMCgenFOAM::Initialize   ======");
  BXTXT(*f_Out,"========================================");
  BXTXT(*f_Out,"========================================");
  f_IsInitialized = 1;  /// re-initialization inhibited
  } else {
	  cout<< "----> TMCgenFOAM::Initialize, already initialized "<<endl;
  }
}/// Initialize

///______________________________________________________________________________________
void TMCgenFOAM::Generate()
{
  f_NevGen++;

  if( m_IsFoam5 == 1){
    m_Mode = -5;
    f_FoamI->MakeEvent();         // Foam of base class
    m_WT   = f_FoamI->GetMCwt();  // get weight
    f_TMCgen_NORMA->Fill(-1, m_Xnorm);    // New style
  }//

}//! Generate

///______________________________________________________________________________________
void TMCgenFOAM::Finalize()
{
  TMCgen::Finalize();
  ///   Finalize MC  run, final printouts, cleaning etc.
  BXOPE(*f_Out);
  BXTXT(*f_Out,"****************************************");
  BXTXT(*f_Out,"******     TMCgenFOAM::Finalize   ******");
  BXTXT(*f_Out,"****************************************");
  ///------------------------
  if( m_IsFoam5 == 1){
    Double_t MCresult, MCerror, MCnorm, Errel;
    f_FoamI->Finalize( MCnorm, Errel);  //!
    f_FoamI->GetIntegMC( MCresult, MCerror);  //! get MC integral, should be one
    cout << "**************************************************************"<<endl;
    cout << "**************** TMCgenFOAM::Finalize  ***********************"<<endl;
    cout << "Directly from FOAM: MCresult= " << MCresult << " +- "<<MCerror <<endl;
    cout << "**************************************************************"<<endl;
  }// m_IsFoam5
  ///------------------------
}//!Finalize



///------------------------------------------------------------------------
double TMCgenFOAM::Fyfs(double gam){
	  return exp(-m_ceuler*gam)/TMath::Gamma(1+gam);
}

///------------------------------------------------------------------------
double TMCgenFOAM::gamISR( double svar){
	  return  sqr(m_chini)*2*m_alfpi*( log(svar/sqr(m_beam)) -1);
}

///------------------------------------------------------------------------
double TMCgenFOAM::gamFSR( double svar){
	  return  sqr(m_chfin)*2*m_alfpi*( log(svar/sqr(m_fin)) -1);
}

///------------------------------------------------------------------------
double TMCgenFOAM::gamIFI( double costhe){
	  return  m_chini*m_chfin*2*m_alfpi* log( (1-costhe)/(1+costhe));
}

///------------------------------------------------------------------------
void TMCgenFOAM::MapPlus( double r, double gam, double &v, double &dJac){
// Maping for POSITIVE gam
// Input r in (0,1) is uniform random number or FOAM vriable
// Returned v is distributed according to gam*v^{gam-1}
  double eps = m_eps;
  double Reps,Rat,RV;
  if( fabs(gam*log(eps)) > m_del){
	  Reps = exp(gam*log(eps));
	  RV   = exp(gam*log(m_vvmax));
	  Rat  = Reps/RV;
	  if( r< Rat ){
		  v = 0;  dJac=1/Rat;
	  } else {
		  v = exp((1/gam)*log(r)); // mapping
		  v *= m_vvmax;
		  dJac = RV/(gam/v*exp(gam*log(v)));      // jacobian
	  }
  } else {
	  Reps = 1+gam*log(eps);
	  RV   = 1+gam*log(m_vvmax);
	  Rat  = Reps/RV;
	  if( r< Rat ){
		  v = 0; dJac=1/Rat;
	  } else {
		  v = exp(-(1/gam)*(1-r)); // mapping
		  v *= m_vvmax;
		  dJac = RV/(gam/v);        // jacobian
	  }
  }
  if( v<0 || v>1) {
	  cout<<"STOP in TMCgenFOAM::MapPlus: +++ v = "<<v<<endl;
	  exit(11);
  }
}// MapPlus

///------------------------------------------------------------------------
void TMCgenFOAM::MapMinus( double r, double gam, double &v, double &dJac){
// Maping for NEGATIVE gam
// Input r in (0,1) is uniform random number of FOAM variable
// Returned v is distributed according to gam*v^{gam-1}
// dJac is normalization (part of Jacobian) factor
  double eps = m_eps;
  double Rat, Reps, RV, R1;
  if( fabs(gam*log(eps)) > m_del){
	  Reps = exp(gam*log(eps)); // R(eps)
//	  R1   = 2*Reps-1;          // R(1)
	  RV   = 2*Reps-exp(gam*log(m_vvmax)); // R(vmax)
	  Rat  = Reps/RV;
	  if( r< Rat ){
		  v = 0; dJac= 1/Rat;
	  } else {
		  v    = exp( (1/gam)*log( 2*Reps -RV*r ) ); // mapping
		  dJac = RV/(-gam/v*exp(gam*log(v)));        // jacobian
	  }
  } else {
	  Reps = 1+gam*log(eps);        // R(eps)
	  R1   = 1+2*gam*log(eps);      // R(1)
	  RV   = R1 +gam*log(m_vvmax);  // R(V)
	  Rat  = Reps/RV;
	  if( r< Rat){
		  v = 0; dJac= 1/Rat;
	  } else {
		  v    = exp( (1/gam)*(R1 -r*RV) ); // mapping
		  dJac = RV/(-gam/v);               // jacobian
	  }
  }
  if( v<0 || v>1) {
	  cout<<"STOP in TMCgenFOAM::MapMinus: +++ v = "<<v<<endl;
	  exit(11);
  }
}// MapMinus

///------------------------------------------------------------------------
void TMCgenFOAM::MapIFI( double r, double gam, double &v, double &R){
//// mapping for IFI
if(gam > 0)
	MapPlus(  r, gam, v, R);
else
    MapMinus( r, gam, v, R);
}// MapIFI


///------------------------------------------------------------------------
void TMCgenFOAM::GetRhoISR1(double svar, double vv, double &rho, double &rho2){
/// First order rho-function dsigma/dv for ISR
  double gami   = gamISR(svar);
  double dels   = 3.0/4.0*gami +m_alfpi*(-0.5  +sqr(m_pi)/3.0) ;
  if( vv > m_eps){
     rho  = gami *(1 +sqr(1-vv))/(2*vv);
     rho2 = rho -sqr(m_chini)*m_alfpi *vv;
  }else{
	 rho = 1.0 + log(m_eps)*(gami) +dels;
	 rho2 = rho;
  }
}//GetRhoISR1


///------------------------------------------------------------------------
void TMCgenFOAM::GetRhoFSR1(double svar, double vv, double &rho, double &rho2){
/// First order rho-function dsigma/dv for FSR
  double gamf   = gamFSR(svar*(1.0-vv));
  double dels   = 3.0/4.0*gamf +m_alfpi*(-0.5  +sqr(m_pi)/3.0) ;
  if( vv > m_eps){
     rho  = gamf *(1 +sqr(1-vv))/(2*vv);
     rho2 = rho -sqr(m_chfin)*m_alfpi *vv;
 }else{
	 rho = 1.0 + log(m_eps)*(gamf) +dels;
	 rho2 = rho;
  }
}//GetRhoFSR1


///------------------------------------------------------------------------
void TMCgenFOAM::GetRhoIFI1(double svar, double vv, double &rhov, double &rhoc){
/// First order functions dsigma/dv and <2c>dsigma/dv for IFI
  double alfpi   = m_alfpi*m_chfin*m_chini;
  if( vv > m_eps){
     rhov  = 2*alfpi *(-3)*(2-vv)/(2*vv);
     rhoc  = 2*alfpi *(-1)/(2-vv)/vv*(10*(1-vv)+3*vv*vv);
 }else{
	 rhov  = 2*alfpi*3.0*log(1/m_eps);
	 rhoc  = 2*alfpi*5.0*log(1/m_eps);
  }
}//GetRhoIFI1


///------------------------------------------------------------------------
void TMCgenFOAM::GetRhoIFI1c(int KeyDist, double cc, double &rho){
// First order dsigma/dc at v=0 for IFI [nb]
//
  double svar = sqr(m_CMSene);
  double sig0nb  = 4*m_pi* sqr(1/m_alfinv)/(3.0*svar )*m_gnanob; // Born Pointlike [nb]
  double alfpi   = m_alfpi*m_chfin*m_chini;

//SUBROUTINE kksem_ord1c(KeyDist,KFi,KFf,CMSene,cc,reps,Result)
  double RhoSoft;
  kksem_ord1c_(KeyDist,m_KFini, m_KFf, m_CMSene, cc, m_eps, RhoSoft);
//

//*******************************************
// gamma-gamma Box plus soft real interference
  double cp = (1+cc)/2e0;
  double cm = (1-cc)/2e0;
  double RhoGG  =
		  4*(1+sqr(cc))*log(cm/cp)*log(m_eps)
		   +(1+sqr(cc))*( sqr(log(cm)) -sqr(log(cp))
                       -2*bvr_dilog_(cm) +2*bvr_dilog_(cp) )
		+2*cp*log(cm)     -2*cm*log(cp)
        -cc*sqr(log(cm)) -cc*sqr(log(cp)) ;
// rho = RhoGG *alfpi* 3e0/8e0*sig0nb;  // gamma-gamm box only
// gamma-gamma nad gamma-Z boxes plus soft real interference
  rho = RhoSoft *alfpi* 3e0/8e0*sig0nb;
//********************************************

//  if(m_count <200 ){
//	  cout<<"GetRhoIFI1c: cc="<< cc<< " RhoG="<<RhoGG <<"  RhoSoft="<< RhoSoft<< "  rat="<< RhoSoft/RhoGG<<endl;
//  }

}//Rho_FSR1


///------------------------------------------------------------------------
double TMCgenFOAM::RhoISR(int KeyISR, double svar, double vv, double eps){
/// ISR rho-function for ISR
  double alf1   = m_alfpi;
  double gami   = gamISR(svar);
  double gamfac = Fyfs(gami);
  double delb   = gami/4 +alf1*(-0.5  +sqr(m_pi)/3.0);
  double ffact  = gamfac*exp(delb);
//
  double rho,dels,delh;
  if(       KeyISR == 0){
/// zero   order exponentiated
	dels = 0;
	delh = -gami/4.0 *log(1.0-vv);
  }else if( KeyISR == 1){
/// first  order
	dels = gami/2;   /// NLO part =0 as for vector boson???
    delh = vv*(-1 +vv/2);
  }else if( KeyISR == 2){
/// second order without NLO part
    dels = gami/2 +sqr(gami)/8;
    delh = vv*(-1+vv/2.0)
          +gami*0.5*(-0.25*(4.0-6.0*vv+3.0*vv*vv)*log(1-vv)-vv);
  }else{
	 cout<<"+++++TMCgenFOAM::KKdistr: Wrong KeyISR = " << m_KeyISR<<endl; exit(5);
  }
  if( vv > eps){
     rho = ffact*gami* exp( log(vv)*(gami-1) ) *(1 +dels +delh);
  }else{
	 rho = ffact* exp( log(eps)*(gami) ) *(1 +dels);
  }
  return rho;
}//Rho_ISR



///------------------------------------------------------------------------
double TMCgenFOAM::RhoFSR(int KeyFSR, double svar, double uu, double eps){
/// ISR+FSR rho-function

//////// from KKsem_uurho(), keyd=302
//  sprim  = svar*(1-uu)
//  bilg   = DLOG(sprim/amfin**2)
//  betf   = Chf2* 2d0*alf1*(bilg-1d0)
//  delb   = Chf2* alf1*(0.5d0*bilg -1d0  +pi**2/3d0)
//  delb   = delb -betf/2 *dlog(1-uu)
//  gamfac = exp(-ceuler*betf)/KKsem_Gamma(1d0+betf)
//  ffact  = gamfac*exp(delb)
////////////////
  double alf1   = m_alfpi;
  double gamf   = gamFSR(svar*(1-uu));
  double delb   = gamf/4 +alf1*(-0.5  +sqr(m_pi)/3.0)
		         -gamf/2 *log(1-uu);
  double ffact  = Fyfs(gamf)*exp(delb);

  double rho,dels,delh;
  if(       KeyFSR == 0){
/// zero   order exponentiated
	dels = 0;
    delh = -gamf/4.0 *log(1.0-uu);
//	rho  = ffact*gamf* exp( log(uu)*(gamf-1) ) *(1 +dels +delh);
  }else if( KeyFSR == 1){
/// first  order
	dels = gamf/2;   /// NLO part =0 as for vector boson???
    delh = uu*(-1 +uu/2);
//    rho = ffact*gamf* exp( log(uu)*(gamf-1) ) *(1 +dels +delh);
  }else if( KeyFSR == 2){
/// from KKsem_uurho(), keyd=302, 2nd ord. without NLO part
//      dels  = betf/2d0 +betf**2/8d0
//      delh  = uu*(-1d0+uu/2d0)
//     $        +betf*(-0.5d0*uu-0.25d0*uu*(-1d0+0.5d0*uu)*log(1d0-uu))
	dels = gamf/2 +sqr(gamf)/8;
    delh = uu*(-1+uu/2.0)
          +gamf*0.5*( -0.5*uu -0.25*uu*(-1.0 +0.5*uu)*log(1-uu));
//    rho = ffact*gamf* exp( log(uu)*(gamf-1) ) *(1 +dels +delh);
  }else{
	  cout<<"+++++TMCgenFOAM::KKdistr: Wrong KeyISR = " << KeyFSR<<endl; exit(5);
  }
  if( uu > eps) {
	  rho = ffact*gamf* exp( log(uu)*(gamf-1) ) *(1 +dels +delh);
  }else{
	  rho = ffact*exp( log(eps)*gamf ) *(1 +dels);
  }
  return rho;
}//RhoFSR


///------------------------------------------------------------------------
double TMCgenFOAM::Soft_yfs(double gam){
/// YFS soft limit coefficient for testing
 double ffact  = Fyfs(gam)*exp(  gam/4 +m_alfpi*(-0.5  +sqr(m_pi)/3.0) );
 double rho = ffact * (1+ gam/2 +sqr(gam)/8);
 return rho;
}


///--------------------------------------------------------------
double TMCgenFOAM::RhoIFI(double costhe, double uu, double eps){
/// ISR+FSR rho-function
  double rho, gami;
  gami   = gamIFI( costhe );
  if( uu > eps){
	rho = gami*exp( log(uu)*(gami-1) );
  } else {
	rho = exp( log(eps)*gami );
  }
  rho *= Fyfs(gami);
  return rho;
}//RhoIFI


///________________________________________________________________________
double TMCgenFOAM::Density(int nDim, double *Xarg){
	//
	if( abs(m_Mode) == 5 ){
	    return Density5(nDim, Xarg);
	} else if( abs(m_Mode) == 3 ){
		return Density3(nDim, Xarg);
	} else if( abs(m_Mode) == 1 ){
		return Density1(nDim, Xarg);
	} else if( abs(m_Mode) == 2 ){
		return Density2(nDim, Xarg);
	} else {
		cout<<" TMCgenFOAM::Density: wrong Mode ="<<m_Mode<<endl;
		exit(-9);
	}
}// Density

///////////////////////////////////////////////////////////////
double TMCgenFOAM::Density5(int nDim, double *Xarg)
{ // density distribution for Foam
	m_count++;  // counter for debug
	//
	Double_t Dist=1;
	double svar = sqr(m_CMSene);
	double svarCum = svar;

// ******** mapping for ISR *******
	double R= Xarg[0];
	double gami = gamISR(svar);
	double dJacISR;
	MapPlus(  R, gami, m_vv, dJacISR);
	double RhoIsr  = RhoISR(2,svar,m_vv,m_eps);
	double DistISR = dJacISR * RhoIsr;
	double RhoIsr0 = RhoISR(0,svar,m_vv,m_eps);
	svarCum *= (1-m_vv);
	double svar2 = svar*(1-m_vv);
// ******** mapping for FSR *******
    double rr= Xarg[1];
    double gamf   = gamFSR(svar2);
    double dJacFSR;
	MapPlus(  rr, gamf, m_uu, dJacFSR);
 	double RhoFsr  = RhoFSR(2, svar2,m_uu,m_eps);
 	double DistFSR = dJacFSR *RhoFsr;
 	double RhoFsr0 = RhoFSR(0, svar2,m_uu,m_eps);
    svarCum *= (1-m_uu);
    // ******** mapping for polar angle *******
    double cmax = 0.99999;
    m_CosTheta = cmax*( -1.0 + 2.0* Xarg[2] );
    Dist *= 2.0*cmax;
    // ******** mapping for IFI variable *******
    double gamint = gamIFI(m_CosTheta);
    double dJacInt1, dJacInt2;
    MapIFI( Xarg[3], gamint, m_r1, dJacInt1);           // mapping eps-dependent !!!
    MapIFI( Xarg[4], gamint, m_r2, dJacInt2);           // mapping eps-dependent !!!
    double RhoInt1 = RhoIFI( m_CosTheta, m_r1,m_eps);  // implicitly eps-dependent !!!
    double RhoInt2 = RhoIFI( m_CosTheta, m_r2,m_eps);  // implicitly eps-dependent !!!
    double DistIFI1 = dJacInt1 *RhoInt1;
    double DistIFI2 = dJacInt2 *RhoInt2;
    Dist *= DistISR* DistFSR* DistIFI1 *DistIFI2;
// ******* MC event *******
    double zz = (1-m_vv)*(1-m_uu)*(1-m_r1)*(1-m_r2);
    m_xx = 1-zz;
//    m_xx = m_vv + m_uu - m_vv*m_uu;  // numerically more stable
// effective masses
	m_Mka = sqrt(svar2);   // after ISR obsolete !!!
// =============== Sigm/dOmega from spin amplitudes ===============
// Effective 4-momenta, KKMC convention: p={px,py,pz,E)
	double Ene = sqrt(svar2)/2;
	double Pmb  = sqrt( (Ene-m_beam)*(Ene+m_beam) ); // modulus
	Vdef(m_p1, 0, 0 , Pmb, Ene);  // beam
	Vdef(m_p2, 0, 0 ,-Pmb, Ene);  // beam
	double Pmf  =sqrt( (Ene-m_fin)*(Ene+m_fin) ); // modulus
	Vdef(m_p3, Pmf*sqrt(1-sqr(m_CosTheta)), 0 , Pmf*m_CosTheta,  Ene); // final
	Vdef(m_p4,-Pmf*sqrt(1-sqr(m_CosTheta)), 0 ,-Pmf*m_CosTheta,  Ene); // final
	double PX[4] = {0, 0, 0, 2*Ene};
	double dSig_GPSF1,dSig_GPSF2, Misr1,Misr2;
	Misr1 = sqrt((1-m_vv)*(1-m_r1)*svar);
	Misr2 = sqrt((1-m_vv)*(1-m_r2)*svar);
//
// Three-stroke calculation of Re(M M^*) including boxes
	gps_bornfoam_( 20,m_KFini,m_KFf,Misr1,m_CosTheta,dSig_GPSF1);
	gps_bornfoam_( 21,m_KFini,m_KFf,Misr2,m_CosTheta,dSig_GPSF2);
    double dBorn_GPS = gps_makerhofoam_(1.0);
//
// Re(M M^*) including only leading part on gamma-Z box
    gps_bornfoam_(  0,m_KFini,m_KFf,Misr1,m_CosTheta,dSig_GPSF1);
    gps_bornfoam_(  1,m_KFini,m_KFf,Misr2,m_CosTheta,dSig_GPSF2);
    double dBorn_GPS0 = gps_makerhofoam_(1.0);
//************ Debug*** Debug*** Debug*** Debug*** Debug ***********
//    if( m_count <10 && fabs(svar/svar2-1)>0.20 ){  // debug
    if( m_count <10 ){  // debug
    	double Rat;
    	Rat = dSig_GPSF1/( dSig_GPSF2 );
    	cout<<" =============================================== "<< m_count<< endl;
    	cout<<" Density5 debug m_count= "<< m_count<< endl;
    	cout<<" dSig_GPSF1    = "<< dSig_GPSF1;
    	cout<<" dSig_GPSF2    = "<< dSig_GPSF2;
    	cout<<" svar/svar2 = "<< svar/svar2;
    	cout<<" Rat = "<<Rat<<endl;
    } //
    if( m_count <10 ){  // debug
//    if( m_count <1 && m_r1 > 0 && gamint <0 ){  // debug
    	cout<<" Density5 debug m_count= "<< m_count<< endl;
    	cout<<" m_r1= "<< m_r1 <<"  m_r2="<< m_r2<<"  m_xx="<< m_xx <<endl;
    	cout<<" m_CosTheta ="<< m_CosTheta <<" gamint= "<<gamint<<endl;
    	cout<<" WT1 ="<< DistIFI1 <<"  WT2="<< DistIFI2<<endl;
   }
//   **********  end debug **********
//	Dist *=  dBorn_GPS;  // This basic distr. includes boxes!!!
    Dist = dJacISR*RhoIsr *dJacFSR*RhoFsr *DistIFI1*DistIFI2 *2.0*cmax *dBorn_GPS;
//****************************************
//          Model weights
//****************************************
// Pure CEEX alf^0 distribution without finite parts of boxes
    double Dist0 = dJacISR*RhoIsr0 *dJacFSR*RhoFsr0 *DistIFI1*DistIFI2 *2.0*cmax *dBorn_GPS0;
    m_WTmodel[50] = 0.0;
    if( Dist != 0.0){
      m_WTmodel[50] = Dist0/Dist; // Auxiliary model weight, pure CEEX0
    }//
//****************************************
//*****  principal distribution for FOAM
    double sig0nb = 4*m_pi* sqr(1/m_alfinv)/(3.0*svar2 )*m_gnanob;
	Dist *=  3.0/8.0 *sig0nb;
	if( svarCum < sqr(2*m_fin)) Dist = 1e-100;
	if(m_Mode > 0 ) Dist = fabs(Dist); // For initialization mode
	return Dist;
}// Density5



///////////////////////////////////////////////////////////////
Double_t TMCgenFOAM::Density2(int nDim, Double_t *Xarg)
{ // density distribution for Foam
	m_count++;  // counter for debug
	//
	Double_t Dist=1;
	double svar = sqr(m_CMSene);
/////////////////////////////////////////////////////////
// ******** ISR *******
	double gami   = gamISR(svar);
    double gamf   = gamFSR(svar);
    // ******** mapping for polar angle *******
    m_CosTheta = -1.0 + 2.0* Xarg[0];
    Dist *= 2.0;
    // bremsstrahlung part
    double RhoIsr2,RhoFsr2,vvcut;
    vvcut = 0.02;
	RhoIsr2 = RhoISR(2, svar,vvcut*0.99999,vvcut);
 	RhoFsr2 = RhoFSR(2, svar,vvcut*0.99999,vvcut);
 	double Rho_cut02= RhoIsr2*RhoFsr2;
    vvcut = 0.002;
	RhoIsr2 = RhoISR(2, svar,vvcut*0.99999,vvcut);
 	RhoFsr2 = RhoFSR(2, svar,vvcut*0.99999,vvcut);
 	double Rho_cut002= RhoIsr2*RhoFsr2;
    vvcut = 0.0002;
	RhoIsr2 = RhoISR(2, svar,vvcut*0.99999,vvcut);
 	RhoFsr2 = RhoFSR(2, svar,vvcut*0.99999,vvcut);
 	double Rho_cut0002= RhoIsr2*RhoFsr2;
 	Dist *= Fyfs(gami+gamf)/Fyfs(gami)/Fyfs(gamf);
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    double Dist_EEX, Dist_GPS;
    double sig0nb = 4*m_pi* sqr(1/m_alfinv)/(3.0*svar )*m_gnanob;
    double BetaFin = sqrt(1-4*sqr(m_fin)/svar ); // phase space factor
	bornv_interpogsw_(m_KFf,svar, m_CosTheta);
	double dSig_EEX = bornv_dizet_( 1, m_KFini, m_KFf, svar, m_CosTheta, 0.0, 0.0, 0.0, 0.0);
	Dist_EEX =  dSig_EEX   *3.0/8.0 *sig0nb;  // Born of EEX
//Dist_EEX = BetaFin*(1+ sqr(m_CosTheta)) *3.0/8.0 *sig0nb;
/////////////////////////////////////////////////////////////////
// =============== Sigm/dOmega from spin amplitudes ===============
// Effective 4-momenta, KKMC convention: p={px,py,pz,E)
	double Ene = sqrt(svar)/2;
	double Pmb  = sqrt( (Ene-m_beam)*(Ene+m_beam) ); // modulus
	Vdef(m_p1, 0, 0 , Pmb, Ene);  // beam
	Vdef(m_p2, 0, 0 ,-Pmb, Ene);  // beam
	double Pmf  =sqrt( (Ene-m_fin)*(Ene+m_fin) ); // modulus
	Vdef(m_p3, Pmf*sqrt(1-sqr(m_CosTheta)), 0 , Pmf*m_CosTheta,  Ene); // final
	Vdef(m_p4,-Pmf*sqrt(1-sqr(m_CosTheta)), 0 ,-Pmf*m_CosTheta,  Ene); // final
	double PX[4] = {0, 0, 0, 2*Ene};
//***** pure Born of CEEX, boxes included
	double dSig_GPS;
    gps_bornf_(m_KFini, m_KFf ,PX, m_CosTheta, m_p1,m_beam, m_p2, -m_beam,
                                               m_p3,m_fin,  m_p4, -m_fin,   dSig_GPS);
	Dist_GPS =  Dist* dSig_GPS   *3.0/8.0 *sig0nb *BetaFin;  // Born of CEEX2
////////////////////////////////////////////////////////////////
	Dist *= Dist_EEX *Rho_cut02;
    m_WTmodel[72] = 0.0;
    m_WTmodel[73] = 0.0;
    if( Dist != 0.0){
      m_WTmodel[72] = Rho_cut002 /Rho_cut02;
      m_WTmodel[73] = Rho_cut0002/Rho_cut02;
    }//

    if(m_Mode > 0 ) Dist = fabs(Dist); // For initialization mode
    return Dist; // principal distribution for FOAM
}//Density2

///////////////////////////////////////////////////////////////
Double_t TMCgenFOAM::Density3(int nDim, Double_t *Xarg)
{ // density distribution for Foam
	m_count++;  // counter for debug
	//
	Double_t Dist=1;
	double svar = sqr(m_CMSene);
	double svarCum = svar;
/////////////////////////////////////////////////////////
// ******** mapping for ISR *******
	double gami = gamISR(svar);
	double R= Xarg[0];
	double dJac,RhoIsr2,RhoFsr2,RhoIsr0,RhoFsr0;
	MapPlus(  R, gami, m_vv, dJac);
	RhoIsr2 = RhoISR(2, svar,m_vv,m_eps);
	RhoIsr0 = RhoISR(0, svar,m_vv,m_eps);
	Dist *= dJac *RhoIsr2;
	svarCum *= (1-m_vv);
	double svar2 = svar*(1-m_vv);
///////////////////////////////////////////////////////
// ******** mapping for FSR *******
    double rr= Xarg[1];
    double gamf   = gamFSR(svar2);
	MapPlus(  rr, gamf, m_uu, dJac);
 	RhoFsr2 = RhoFSR(2, svar2,m_uu,m_eps);
 	RhoFsr0 = RhoFSR(0, svar2,m_uu,m_eps);
 	Dist *= dJac* RhoFsr2;
    svarCum *= (1-m_uu);
////////////////////////////////////////////////////////////
// ******** mapping for polar angle *******
    m_CosTheta = -1.0 + 2.0* Xarg[2];
    Dist *= 2.0;
//
    double zz = (1-m_vv)*(1-m_uu);
    m_xx = 1-zz;
    m_xx = m_vv + m_uu - m_vv*m_uu;  // numerically more stable
    //!!!
    if( m_xx > m_vvmax) return 1e-100;
// ******** Born in early version *******
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    long KeyFob =0;
    KeyFob=   10; // BornV_Dizet, with EW and without integration ???
    KeyFob=  -11; // BornV_Simple, for KeyLib=0, NO EW, NO integration OK
    KeyFob=  -10; // KKsem_BornV, NO EW, NO integration OK!
    KeyFob= -100; // KKsem_BornV, NO EW, WITH integration, OK
    KeyFob=    0; // With EW (BornV_Dizet) With integration OK!
//  -----------------
	kksem_setkeyfob_( KeyFob );
	double sigBornEEX, sigBornEEX0;
//***** Integrated EEX Born from KKMC
	kksem_makeborn_( svar2, sigBornEEX);
	kksem_makeborn_( svar,  sigBornEEX0); // svar2<-svar
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// In BornV_Differential we see:
//    CALL BornV_InterpoGSW( ABS(KFf),  svar, CosThe)
//    Born= BornV_Dizet( 1,m_KFini,KFf, svar, CosThe, eps1,eps2,ta,tb)
//    Born = 4*pi*alfa**2/(3d0*svar )*BornY  *m_gnanob
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	bornv_interpogsw_(m_KFf,svar2, m_CosTheta);
	double dSig_EEX = bornv_dizet_( 1, m_KFini, m_KFf, svar2, m_CosTheta, 0.0, 0.0, 0.0, 0.0);
//[[[[
//	bornv_interpogsw_(m_KFf,svar, m_CosTheta);
//	double dSig_EEX  = bornv_dizet_( 1, m_KFini, m_KFf, svar, m_CosTheta, 0.0, 0.0, 0.0, 0.0); // svar2->svar
	double dSig_EEX0 = bornv_dizet_( 1, m_KFini, m_KFf, svar2,       0.0, 0.0, 0.0, 0.0, 0.0);
//]]]
// ******* effective masses *********
	m_Mka = sqrt(svar2);   // after ISR obsolete
// =============== Sigm/dOmega from spin amplitudes ===============
// Effective 4-momenta, KKMC convention: p={px,py,pz,E)
	double Ene = sqrt(svar2)/2;
	double Pmb  = sqrt( (Ene-m_beam)*(Ene+m_beam) ); // modulus
	Vdef(m_p1, 0, 0 , Pmb, Ene);  // beam
	Vdef(m_p2, 0, 0 ,-Pmb, Ene);  // beam
	double Pmf  =sqrt( (Ene-m_fin)*(Ene+m_fin) ); // modulus
	Vdef(m_p3, Pmf*sqrt(1-sqr(m_CosTheta)), 0 , Pmf*m_CosTheta,  Ene); // final
	Vdef(m_p4,-Pmf*sqrt(1-sqr(m_CosTheta)), 0 ,-Pmf*m_CosTheta,  Ene); // final
	double PX[4] = {0, 0, 0, 2*Ene};
//***** pure Born of CEEX
	double BetaFin = sqrt(1-4*sqr(m_fin)/svar2 ); // missing phase space factor
	double dSig_GPS;
    gps_bornf_(m_KFini, m_KFf ,PX, m_CosTheta, m_p1,m_beam, m_p2, -m_beam,
                                               m_p3,m_fin,  m_p4, -m_fin,   dSig_GPS);
    dSig_GPS *= BetaFin;  // missing phase space factor
/////////////////////////////////////////////////////////////////
    double Dist_EEX, Dist_GPS;
    double sig0nb = 4*m_pi* sqr(1/m_alfinv)/(3.0*svar2 )*m_gnanob;
 	Dist_EEX =  Dist* dSig_EEX   *3.0/8.0 *sig0nb;  // Born of EEX
	Dist_GPS =  Dist* dSig_GPS   *3.0/8.0 *sig0nb;  // Born of CEEX
	if( svarCum < sqr(2*m_fin)) Dist = 1e-100;
/////////////////////////////////////////////////////////////////
//    double dSigRef = bornv_dizet_( 1, m_KFini, m_KFf, svar2, 0.0 , 0.0, 0.0, 0.0, 0.0); // at cos(theta)=0
//************ Debug*** Debug*** Debug*** Debug*** Debug ***********
      if( m_count < 1 ){  // debug
//    if( m_count <10000 && m_xx<1e-14 ){  // debug
//    if( m_count <100000 && fabs(dSig_GPS/dSig_EEX -1) >0.10 ){  // debug
//    if( m_count <10000 && fabs(dSig_GPS-dSig_EEX)/dSigRef >0.002 ){  // debug
    	cout<<" ******************** Density3 debug m_count= "<< m_count<< "  Dist= " << Dist<< endl;
        double SoftIni = Soft_yfs(gami);
//    	cout<<" (dSig_GPS-dSig_EEX)/ref  = "<< (dSig_GPS-dSig_EEX)/dSigRef ;
//   	  // Born+boxes, WARNING Z-box may be modified for KeyZet=2
//      double dSig_GPSF0,dSig_GPSF1;
//    	gps_bornfoam_( 20,m_KFini,m_KFf,m_Mka,m_CosTheta,dSig_GPSF0);
//    	cout<<" dSig_GPSF0/dSig_EEX = "<< dSig_GPSF/dSig_EEX;
//    	gps_bornfoam_( 21,m_KFini,m_KFf,m_Mka,m_CosTheta,dSig_GPSF1);
//      double dSig_GPSFR = gps_makerhofoam_(1.0);
//      cout<<" // dSig_GPSFR: "<< (dSig_GPSFR-dSig_EEX)/dSigRef;
//      cout<<"    dSig_GPSF0: "<< (dSig_GPSF0-dSig_EEX)/dSigRef;
//
    	double CosTheta = 0.0;
    	bornv_interpogsw_(m_KFf,svar2, CosTheta);
    	double dSig_EEXc = bornv_dizet_( 1, m_KFini, m_KFf, svar2, CosTheta, 0.0, 0.0, 0.0, 0.0);
    	cout<< "dSig_EEXc = "<< dSig_EEXc <<endl;
    	cout<< "dSig_EEX0 = "<< dSig_EEX0 <<endl;
    	cout<< "             Sig_EEX0   = "<< dSig_EEX0 *sig0nb <<endl;
  	    cout<< " From KKsem: sigBornEEX = "<< sigBornEEX <<endl;
 	    cout<< " From KKsem: sigBornEEX0= "<< sigBornEEX0 <<endl;
    	//
    	cout<<" m_CosTheta= "<< m_CosTheta;
    	cout<<" m_xx= "<< m_xx<<" m_vv= "<< m_vv<<" m_uu= "<< m_uu << "  m_vvmax="<<m_vvmax;
        cout<<endl;
        ///////////// testing soft limit /////////////
     	double SoftFin = Soft_yfs(gamf);
        gamf =  gamFSR(svar);
    	double    softISR = RhoISR(2, svar,m_vv,m_eps)  / exp( (gami-1)* log(m_vv))/gami;
    	cout<<"  gami="<< gami<<"  softISR = "<< softISR<< " SoftIni= "<<SoftIni<<endl;
    	double    softFSR = RhoFSR(2, svar2,m_uu,m_eps) / exp( (gamf-1)* log(m_uu))/gamf;
    	cout<<"  gamf="<< gamf<<"  softFSR = "<< softFSR<< " SoftFin= "<<SoftFin;
    	cout<<"  SoftIni* SoftFin = "                <<  SoftIni* SoftFin  <<endl;
    	cout<<"  exp( (gami+gamf) * log(m_vvmax)) = "<<  exp( (gami+gamf) * log(m_vvmax)) <<endl;
    	double fudge = TMath::Gamma(1+gami)*TMath::Gamma(1+gamf)/TMath::Gamma(1+gami+gamf);
    	cout<<"  Gamma(1+gami)*Gamma(1+gamf)/Gamma(1+gami+gamf)= " << fudge <<endl;
    	double DistEstim  = sigBornEEX0 *SoftIni*SoftFin *exp( (gami+gamf) * log(m_vvmax));
    	double DistEstim2 = DistEstim *fudge;
    	double SoftKKsem  = sigBornEEX0 *SoftIni*SoftFin * fudge;
    	cout<<"          SoftKKsem  = "<< SoftKKsem <<endl;
    	cout<<" !!! Estimated Dist  = "<< DistEstim << "   "<< DistEstim2 <<endl;
    	cout<<" !!! gami+gamf  = "<< gami+gamf<< " " <<endl;
    } // end debug **********
//****************************************
//          Model weights
//****************************************
    m_WTmodel[ 2] = 0.0;
    m_WTmodel[ 3] = 0.0;
    m_WTmodel[52] = 0.0;
    if( Dist_EEX != 0.0){
    	m_WTmodel[ 2] = Dist_GPS/Dist_EEX; // Auxiliary model weight
    	m_WTmodel[ 3] = RhoIsr0/RhoIsr2 *RhoFsr0/RhoFsr2; // Auxiliary model weight
    	m_WTmodel[52] = Dist_GPS/Dist_EEX *RhoIsr0/RhoIsr2 *RhoFsr0/RhoFsr2;
    }//
// principal distribution for FOAM, always positive
//[[[	return Dist_EEX; // principal distribution for FOAM
    return Dist;
//
}// Density3

///////////////////////////////////////////////////////////////
Double_t TMCgenFOAM::Density1(int nDim, Double_t *Xarg)
{ // density distribution for Foam
	m_count++;  // counter for debug
	//
	Double_t Dist, xDist, yDist, xDist1, yDist1;
	double svar = sqr(m_CMSene);
/////////////////////////////////////////////////////////
// ******** mapping for ISR *******
	double gami   = gamISR(svar);
    double gamf0  = gamFSR(svar);
	double R= Xarg[0];
	double dJac,gamm;
	gamm = gami+gamf0;
	MapPlus(  R, gamm, m_vv, dJac); // with m_eps IR cut-off
//[[[
//	m_vv = R *m_vvmax; dJac = m_vvmax;
//]]]
	double cc= (2*Xarg[1]-1)*0.99999999;
//*****************************************
	double xRhoI, xRhoF, xRhoIFI;
	double yRhoI, yRhoF, yRhoIFI;
	GetRhoISR1(svar, m_vv, xRhoI,   yRhoI);     // ISR complete 1st ord
	GetRhoFSR1(svar, m_vv, xRhoF,   yRhoF);     // FSR complete 1st ord
	GetRhoIFI1(svar, m_vv, xRhoIFI, yRhoIFI);   // IFI 1st order without virt. cors.
//*****************************************
//  Born convolution with ISR and FSR
	double xBorn00,xBornv, xBorn0, xBornv2, yBornv, yBorn0, yBornv2;
//  x_Born((1-v)*s) and x_Born(s) provided by KKsem
//  SUBROUTINE KKsem_Ord1v(KeyDist,KFi,KFf,CMSene,vv,Result)
	kksem_ord1v_(  1,m_KFini, m_KFf, m_CMSene, m_vv, xBornv);   // Born [nb]
	kksem_ord1v_(101,m_KFini, m_KFf, m_CMSene, m_vv, xBornv2);  // Born [nb]
	kksem_ord1v_(  1,m_KFini, m_KFf, m_CMSene,  0e0, xBorn0);   // Born [nb]
//  sigmaBorn<2c> provided by KKsem
	kksem_ord1v_(  2,m_KFini, m_KFf, m_CMSene, m_vv, yBornv);   // sigma<2c> [nb]
	kksem_ord1v_(102,m_KFini, m_KFf, m_CMSene, m_vv, yBornv2);  // sigma<2c> [nb]
	kksem_ord1v_(  2,m_KFini, m_KFf, m_CMSene,  0e0, yBorn0);   // sigma<2c> [nb]
//*******************************************************************
//  Virtual+soft IFI contributions (at v=0) integrated analyticaly over cos(theta)
	double xVirt, yVirt;
	kksem_ord1v_( 10,m_KFini, m_KFf, m_CMSene,  0e0, xVirt);  // virt+soft for sigma
	kksem_ord1v_( 20,m_KFini, m_KFf, m_CMSene,  0e0, yVirt);  // virt+soft for sigma<2c>
//*******************************************************************
//  Virtual+soft contributions at vv=0 with live cos(theta) dependence
    double xRhoIFIc,xRhoIFIcIR;
    GetRhoIFI1c(  1, cc, xRhoIFIc);
    GetRhoIFI1c(111, cc, xRhoIFIcIR);  // IR subtraction term
//*******************************************************************
    double xIRv2, yIRv2;   // IR subtraction term
    kksem_ord1v_(111,m_KFini, m_KFf, m_CMSene, m_vv, xIRv2);    // Born [nb]
    kksem_ord1v_(112,m_KFini, m_KFf, m_CMSene, m_vv, yIRv2);    // sigma<2c> [nb]
//*******************************************************************
    double yDist2=0, yDist3=0, yDist4=0;
// Additive combination of RhoI and RhoF (ISR+FSR)
	if( m_vv > m_eps){     // HARD part, integrated over c
		 //[[[[[[[[[[[[
		 //xRhoF=0;yRhoF=0;    // FSR off
		 //xRhoI=0;yRhoI=0;    // ISR off
		 //]]]]]]]]]]]
	     xDist  = xRhoI*xBornv + xRhoF*xBorn0;    // ISR+FSR dsig/dv
	     yDist  = yRhoI*yBornv + yRhoF*yBorn0;    // ISR+FSR <2c>dsig/dv
	     xDist1 = xDist + xRhoIFI*yBornv2;        // ISR+FSR+IFI dsig/dv
	     yDist1 = yDist + yRhoIFI*xBornv2;        // ISR+FSR+IFI <2c>dsig/dv
	     yDist2 =         yRhoIFI*xBornv2 -xIRv2; // IFI hard part with IR subtracted
	     yDist4 =         yRhoIFI*xBornv2;        // IFI alone
	}else{                //  SOFT+VIRT
		 //[[[[[[[[[[[[
		 //xRhoF=1;yRhoF=1;   // FSR off
		 //xRhoI=1;yRhoI=1;   // ISR off
		 //]]]]]]]]]]]
		 xDist  = (1 +(xRhoI-1) +(xRhoF-1) )*xBorn0;  // ISR+FSR sig0
		 yDist  = (1 +(yRhoI-1) +(yRhoF-1) )*yBorn0;  // ISR+FSR <2c>sig0
// version with fully analytical c-integration (PLB219), does not work
//	     xDist1 = xDist + xRhoIFI*yBorn +xVirt;       // ISR+FSR+IFI sig0
//	     yDist1 = yDist + yRhoIFI*xBorn +yVirt;       // ISR+FSR+IFI <2c>sig0
// version with partial analytical integration over c of IR part
//	     xDist1 = xDist + xRhoIFI*yBorn         +2*xRhoIFIc;   // ISR+FSR+IFI sig0
//	     yDist1 = yDist + yRhoIFI*xBorn   +2*(2*c)*xRhoIFIc;   // ISR+FSR+IFI <2c>sig0
//	     xDist1 = xDist + xRhoIFI*yBorn;   // ISR+FSR+IFI sig0
//	     yDist1 = yDist + yRhoIFI*xBorn;   // IxRhoIFIcIRSR+FSR+IFI <2c>sig0
// version with MC integration over c of IR part and the rest
	     xDist1 = xDist        +2*xRhoIFIc;    // ISR+FSR+IFI 0sig: 2=jacob=d(c)/d(r)
	     yDist1 = yDist +(2*cc)*2*xRhoIFIc;    // ISR+FSR+IFI <2c>0sig
	     yDist3 =        (2*cc)*2*(xRhoIFIc -xRhoIFIcIR); // with IR subtraction
	     yDist4 =        (2*cc)*2*xRhoIFIc;    // IFI alone
		 //dJac *= 1/m_eps;
	}//
	// model WT for AFB without IFI
	for(int i=0; i<100; i++) m_WTmodel[i] = 0.0;
// Final weight is typicaly yDist1/xDist*xDist/fabs(xDist)= yDist/fabs(xDist) !
    if( xDist != 0.0){
    	m_WTmodel[10] = yDist /xDist;  // WT for AFB without IFI
    	m_WTmodel[11] = xDist1/xDist;  // WT for sig with IFI on
    	m_WTmodel[12] = yDist1/xDist;  // WT for AFB with IFI on
    	m_WTmodel[22] = yDist2/xDist;  // WT for AFB with IFI, non-IR hard part
    	m_WTmodel[23] = yDist3/xDist;  // WT for AFB with IFI, non-IR soft part
    	m_WTmodel[24] = yDist4/xDist;  // WT for AFB from IFI alone
    }//
//************  principal distribution for FOAM
	Dist = xDist* dJac;
	if(m_Mode > 0 ) Dist = fabs(Dist); // For initialization mode
	return Dist;
	////////////////////////////////////////////////////////////
}// Density1

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                             UTILITIES                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

void TMCgenFOAM::ReaData(const char *DiskFile, int imax, double xpar[])
//////////////////////////////////////////////////////////////
//    subprogram reading input data file and packing        //
//    entries into matrix xpar                              //
//    WARNING: input file cannot include empty lines        //
//    it cannot handle entries like 1d-30, has to be 1e-30! //
//////////////////////////////////////////////////////////////
{
  char trail[200];
  char ch1;
  int  foundB=0, foundE=0, line, indx;
  int  line_max =2000;
  double value;
  cout<<"============================ReaData=============================="<<endl;
  cout<<"===                     "<< DiskFile <<"               =========="<<endl;
  cout<<"================================================================="<<endl;
  ifstream InputFile;
  InputFile.open(DiskFile);
  for(indx=0;indx<imax; indx++) xpar[indx]=0.0;
  for(line=0;line<line_max; line++){
    InputFile.get(ch1);
    if( ch1 == 'B') foundB=1;
    InputFile.getline(trail,200);
    if(foundB) break;
  }
  for(line=0;line<line_max; line++){
    InputFile.get(ch1);
    if( ch1 == 'E'){
      foundE=1;
      InputFile.getline(trail,200);
      cout<<ch1<<trail<<"["<<line<<"]"<<endl;
      break;
    }
    if( ch1 == '*'){
      InputFile.getline(trail,200);
      cout<<ch1<<trail<<endl;
    }else{
      InputFile>>indx>>value;
      if(indx<0 || indx>abs(imax) ){
	cout<<" ++++++++ReaData: wrong indx = "<<indx<<endl;
	exit(0);
      }
      xpar[indx] = value;
      //xpar[indx-1] = value; // correction for fortran indexing in input file
      InputFile.getline(trail,200);
      cout<<ch1;
      cout<<setw(4)<<indx<<setw(15)<<value<<" ";
      cout<<trail<<endl;
    }
  }
  cout<<"================================================================="<<endl;
  InputFile.close();
}

void TMCgenFOAM::Vdef(double v[4], const double v1, const double v2, const double v3, const double v4)
  { // define a 4-vector (avoids initialization warnings)
       v[0] = v1; v[1] = v2; v[2] = v3; v[3] = v4;
  } 
