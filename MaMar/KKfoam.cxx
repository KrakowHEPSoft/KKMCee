//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//               Class   KKfoam                                             //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
// KKfoam is multipurpose toolbox for KKMC testing.
//  1. Interfaces (erappers) to KKMC and KKsem F77 subrograms
//  2. Integrand for Foam
//  3. A few routines for producing tatex table out of histograms
//////////////////////////////////////////////////////////////////////////////

#include "KKfoam.h"

double sqr( const Double_t x ){ return x*x;};

KKfoam::KKfoam(const char* Name)
{
		cout<< "----> KKfoam USER Constructor "<<endl;
}

///////////////////////////////////////////////////////////////////////////////////
void KKfoam::Initialize( double ypar[10000] ){
  //------------------------------------------------------------------------
  cout<<"==================================================================="<<endl;
  cout<<"================ KKfoam initialization begin ======================"<<endl;
  m_jmax =10000;
  for(int j=1; j<=m_jmax; j++)
	  m_ypar[j-1]= ypar[j-1];   // ypar has c++ numbering

  m_CMSene  = m_ypar[ 1 -1];
  m_vvmax   = m_ypar[17 -1];

  m_alfinv   = m_ypar[30 -1];     // 1/alphaQED at Q^2=0
  m_gnanob   = m_ypar[31 -1];     // GeV^2 -> nanobarns

  m_pi      = 3.1415926535897932;
  m_ceuler  = 0.57721566;
//
  m_alfpi   = 1/m_alfinv/m_pi;
//
  m_beam    = 0.510999e-3;  // electron
  m_chini   = 1.0;          // electron

  m_fin     = 0.105;        // final ferm. muon mass
  m_chfin   = 1.0;          // final ferm. muon charge

  m_KFini   = 11; // electron
  m_KFf     = 13; // muon

  m_kDim    =    3;         // No. of dim. for Foam, =2,3 Machine energy spread OFF/ON
  m_nCells  = 2000;         // No. of cells, optional, default=2000
  m_nSampl  =  200;         // No. of MC evts/cell in exploration, default=200

  m_KeyISR  = 2;            // Type of ISR/QED switch, 0,1,2

  m_Mode    = 3;            // Operation mode for Foam
  m_del     = 0.0001;       // limit for gamma*ln(eps) in IFI mapping

  m_count   = 0;            // counter for debug

  cout<<"   KKfoam::Initialize:  m_CMSene= "<< m_CMSene<<endl;
  cout<<"================ KKfoam initialization end  ======================"<<endl;
}// KKfoam::Initialize


///------------------------------------------------------------------------
double KKfoam::Fyfs(double gam){
	  return exp(-m_ceuler*gam)/TMath::Gamma(1+gam);
}

///------------------------------------------------------------------------
double KKfoam::gamISR( double svar){
	  return  sqr(m_chini)*2*m_alfpi*( log(svar/sqr(m_beam)) -1);
}

///------------------------------------------------------------------------
double KKfoam::gamFSR( double svar){
	  return  sqr(m_chfin)*2*m_alfpi*( log(svar/sqr(m_fin)) -1);
}

///------------------------------------------------------------------------
double KKfoam::gamIFI( double costhe){
	  return  m_chini*m_chfin*2*m_alfpi* log( (1-costhe)/(1+costhe));
}

///------------------------------------------------------------------------
void KKfoam::MapIFI1( double r, double gam, double eps, double &v, double &dJac){
// Maping for POSITIVE gam
// Input r in (0,1) is random number
// Returned v is distributed according to gam*v^{gam-1}
  if( fabs(gam*log(eps)) > m_del){
	  double eg = exp(gam*log(eps));
	  if( r< eg ){
		  v = 0;  dJac=1/eg;
	  } else {
		  v = exp((1/gam)*log(r)); // mapping
		  dJac = 1/(r*gam/v);      // jacobian
	  }
  } else {
	  double eg = 1+gam*log(eps);
	  if( r< eg ){
		  v = 0; dJac=1/eg;
	  } else {
		  v = exp(-(1/gam)*(1-r)); // mapping
		  dJac = 1/(gam/v);        // jacobian
	  }
  }
  if( v<0 || v>1) {
	  cout<<"STOP in KKfoam::MapIFI: +++ v = "<<v<<endl;
	  exit(11);
  }
}// MapIFI

///------------------------------------------------------------------------
void KKfoam::MapIFI2( double r, double gam, double eps, double &v, double &dJac){
// Maping for NEGATIVE gam
// Input r in (0,1) is random number
// Returned v is distributed according to gam*v^{gam-1}
// dJac is normalization (part of Jacobian) factor
  double r0, eg;
  if( fabs(gam*log(eps)) > m_del){
	  eg = exp(gam*log(eps));
	  r0 =eg/(2*eg-1);
	  if( r< r0 ){
		  v = 0; dJac= 1/r0;
	  } else {
		  v = exp( (1/gam)*log( 2*eg -(2*eg-1)*r ) ); // mapping
		  dJac = (2*eg-1)/(-gam/v*exp(gam*log(v)));   // jacobian
	  }
  } else {
	  eg = 1+gam*log(eps);
	  r0 = eg/(2*eg-1);
	  if( r< r0){
		  v = 0; dJac= 1/r0;
	  } else {
		  v = exp( (1/gam)*(1-r)*(2*eg-1) ); // mapping
		  dJac = (2*eg-1)/(-gam/v);          // jacobian
	  }
  }
  if( v<0 || v>1) {
	  cout<<"STOP in KKfoam::MapIFI2: +++ v = "<<v<<endl;
	  exit(11);
  }
}// MapIFI2

///------------------------------------------------------------------------
void KKfoam::MapIFI( double r, double gam, double eps, double &v, double &R){
//// mapping for IFI
if(gam > 0)
	MapIFI1( r, gam, eps, v, R);
else
    MapIFI2( r, gam, eps, v, R);
}// MapIFI

///------------------------------------------------------------------------
double KKfoam::Rho_isr(double svar, double vv){
/// ISR rho-function for ISR

  double alf1   = m_alfpi;
  double gami   = gamISR(svar);
  //gami = sqr(m_chini)*2*m_alfpi*( log(svar/sqr(m_beam)) -1);
///
  double gamfac = Fyfs(gami);
  double delb   = gami/4 +alf1*(-0.5  +sqr(m_pi)/3.0);
  double ffact  = gamfac*exp(delb);

  double rho,dels,delh;
  if(       m_KeyISR == 0){
/// zero   order exponentiated
	dels = 0;
	delh = 0;
	rho  = ffact*gami* exp( log(vv)*(gami-1) ) *(1 +dels +delh);
  }else if( m_KeyISR == 1){
/// first  order
	dels = gami/2;   /// NLO part =0 as for vector boson???
    delh = vv*(-1 +vv/2);
    rho = ffact*gami* exp( log(vv)*(gami-1) ) *(1 +dels +delh);
  }else if( m_KeyISR == 2){
/// second order without NLO part
    dels = gami/2 +sqr(gami)/8;
    delh = vv*(-1+vv/2.0)
          +gami*0.5*(-0.25*(4.0-6.0*vv+3.0*vv*vv)*log(1-vv)-vv);
    rho = ffact*gami* exp( log(vv)*(gami-1) ) *(1 +dels +delh);
  }else{
	  cout<<"+++++TMCgenFOAM::KKdistr: Wrong KeyISR = " << m_KeyISR<<endl;
	  exit(5);
  }
///
  return rho;
}//Rho_isr


///------------------------------------------------------------------------
double KKfoam::Rho_fsr(double svar, double uu){
/// ISR+FSR rho-function

  double alf1   = m_alfpi;
  double gamf   = gamFSR(svar*(1-uu));
///
  /////delb   = Chf2* alf1*(0.5d0*bilg -1d0  +pi**2/3d0)
  /////delb   = delb -betf/2 *dlog(1-uu)
  /////double gamfac = exp(-m_ceuler*gamf)/TMath::Gamma(1+gamf);
  double delb   = gamf/4 +alf1*(-0.5  +sqr(m_pi)/3.0)
		         -gamf/2 *log(1-uu);
  double ffact  = Fyfs(gamf)*exp(delb);

  double rho,dels,delh;
  if(       m_KeyISR == 0){
/// zero   order exponentiated
	dels = 0;
	delh = 0;
	rho  = ffact*gamf* exp( log(uu)*(gamf-1) ) *(1 +dels +delh);
  }else if( m_KeyISR == 1){
/// first  order
	dels = gamf/2;   /// NLO part =0 as for vector boson???
    delh = uu*(-1 +uu/2);
    rho = ffact*gamf* exp( log(uu)*(gamf-1) ) *(1 +dels +delh);
  }else if( m_KeyISR == 2){
/// second order without NLO part
    /////dels  = betf/2d0 +betf**2/8d0
    /////delh  = uu*(-1d0+uu/2d0)
    /////$        +betf*(-0.5d0*uu-0.25d0*uu*(-1d0+0.5d0*uu)*log(1d0-uu))
    dels = gamf/2 +sqr(gamf)/8;
    delh = uu*(-1+uu/2.0)
          +gamf*0.5*( -0.5*uu -0.25*uu*(-1.0 +0.5*uu)*log(1-uu));
    rho = ffact*gamf* exp( log(uu)*(gamf-1) ) *(1 +dels +delh);
  }else{
	  cout<<"+++++TMCgenFOAM::KKdistr: Wrong KeyISR = " << m_KeyISR<<endl;
	  exit(5);
  }
//
  //if(rho<0) rho=0;
  return rho;
}//Rho_fsr

///--------------------------------------------------------------
double KKfoam::Rho_ifi(double costhe, double uu, double eps){
/// ISR+FSR rho-function
  double rho, gami;
  gami   = gamIFI( costhe );
  if( fabs(gami*log(eps)) > m_del){
	  if( uu < eps){
		  rho = exp( log(eps)*gami );
	  } else {
		  rho = gami*exp( log(uu)*(gami-1) );
	  }
  } else {
	  if( uu < eps){
		  rho = 1+ gami*log(eps);
	  }  else {
		  rho = gami/uu;
	  }
  }
  rho *= Fyfs(gami);
  return rho;
}//Rho_ifi

///////////////////////////////////////////////////////////////
Double_t KKfoam::Density(int nDim, Double_t *Xarg)
{ // density distribution for Foam
    if (     abs(m_Mode) == 3 )
    	Density3(nDim, Xarg);
    else if( abs(m_Mode) == 5 )
    	Density5(nDim, Xarg);
    else {
    	cout<<"+++ KKfoam::Density: Wrong m_Mode = "<< m_Mode<<endl;
    	exit(20);
    }
}// Density


///////////////////////////////////////////////////////////////
Double_t KKfoam::Density3(int nDim, Double_t *Xarg)
{ // density distribution for Foam
	m_count++;  // counter for debug
	//
	Double_t Dist=1;

	double svar = sqr(m_CMSene);
	double svarCum = svar;

// ******** mapping for ISR *******
	double gamiCR,gami,alfi;
	double CMSene1= sqrt(svar);
	bornv_makegami_( CMSene1, gamiCR,gami,alfi);   // from KKMC
	//[[[[ debug
	//gami = gamISR(CMSene1);
	//]]]]
	// cout<<" CMSene1,gami= "<< CMSene1 <<"  "<< gami <<endl;
	double R= Xarg[0];
	m_vv = exp(1.0/gami *log(R)) *m_vvmax; // mapping
	Dist *= m_vv/R/gami ;                  // Jacobian
	if( gami < 0 )      return 0.0;    // temporary fix
    if( m_vv < 1e-200 ) return 0.0;    // temporary fix
    // ISR photonic distribution
	double Rho2 = Rho_isr(svar,m_vv);  // remember take care of m_mbeam!!!
	Dist *= Rho2;
	svarCum *= (1-m_vv);
	double svar2 = svar*(1-m_vv);

// ******** mapping for FSR *******
    double rr= Xarg[1];
    double gamf   = gamFSR(svar2);
    if( gamf <0 )       return 0.0;      // just in case
    m_uu = exp(1.0/gamf *log(rr));       // mapping
    if( m_uu < 1e-200 ) return 0.0;      // temporary fix
    Dist *= m_uu/rr/gamf;                // Jacobian
// FSR photonic distribution
  	double Rho3 = Rho_fsr(svar2,m_uu);   // remember take care of m_mbeam!!!
  	if( Rho3 <0 ) return 1e-100;
 	Dist *= Rho3;
    svarCum *= (1-m_uu);

    // ******** mapping for polar angle *******
    m_CosTheta = -1.0 + 2.0* Xarg[2];
    Dist *= 2.0;

    double zz = (1-m_vv)*(1-m_uu);
    m_xx = 1-zz;
    m_xx = m_vv + m_uu - m_vv*m_uu;  // numerically more stable

// ******** finally Born factor *******
    long KeyFob;
    KeyFob=   10; // BornV_Dizet, with EW and without integration ???
    KeyFob=  -11; // BornV_Simple, for KeyLib=0, NO EW, NO integration OK
    KeyFob=  -10; // KKsem_BornV, NO EW, NO integration OK!
    KeyFob= -100; // KKsem_BornV, NO EW, WITH integration, OK
    KeyFob=    0; // With EW (BornV_Dizet) With integration OK!
//  -----------------
//	kksem_setkeyfob_( KeyFob );

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//double xBorn;
	//Integrated Born from KKMC
	//kksem_makeborn_( svar2, xBorn);
	//Dist *= xBorn/2.0;
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//    In BornV_Differential:
//    CALL BornV_InterpoGSW( ABS(KFf),  svar, CosThe)
//    Born= BornV_Dizet( 1,m_KFini,KFf, svar, CosThe, eps1,eps2,ta,tb)
//    Born = 4*pi*alfa**2/(3d0*svar )*BornY  *m_gnanob

//
	bornv_interpogsw_(m_KFf,svar2, m_CosTheta);
	double dSig_dCos = bornv_dizet_( 1, m_KFini, m_KFf, svar2, m_CosTheta, 0.0, 0.0, 0.0, 0.0);

// ******* effective masses *********
	m_Mka = sqrt(svar2);   // after ISR

// =============== Sigm/dOmega from spin amplitudes ===============
// Effective 4-momenta, KKMC convention: p={px,py,pz,E)
	double Ene = m_Mka/2;
	double Pmb  = sqrt( (Ene-m_beam)*(Ene+m_beam) ); // modulus
	m_p1 = { 0, 0 , Pmb, Ene};  // beam
	m_p2 = { 0, 0 ,-Pmb, Ene};  // beam
	double Pmf  =sqrt( (Ene-m_fin)*(Ene+m_fin) ); // modulus
	m_p3 = { Pmf*sqrt(1-sqr(m_CosTheta)), 0 , Pmf*m_CosTheta,  Ene}; // final
	m_p4 = {-Pmf*sqrt(1-sqr(m_CosTheta)), 0 ,-Pmf*m_CosTheta,  Ene}; // final
	double PX[4] = {0, 0, 0, 2*Ene};
//***** pure Born of CEEX
	double dSigAng;
    gps_bornf_(m_KFini, m_KFf ,PX, m_CosTheta, m_p1,m_beam, m_p2, -m_beam,
                                               m_p3,m_fin,  m_p4, -m_fin,   dSigAng);
//[[[[[[[
    double dSigRef = bornv_dizet_( 1, m_KFini, m_KFf, svar2, 0.0 , 0.0, 0.0, 0.0, 0.0); // at cos(theta)=0
//************ Debug*** Debug*** Debug*** Debug*** Debug ***********
//    if( m_count <10 && m_xx>0.6 ){  // debug
//    if( m_count <100000 && fabs(dSigAng/dSig_dCos -1) >0.10 ){  // debug
    if( m_count <10000 && fabs(dSigAng-dSig_dCos)/dSigRef >0.002 ){  // debug
    	cout<<" ******************** Density3 debug m_count= "<< m_count<< endl;
    	cout<<" (dSigAng-dSig_dCos)/ref  = "<< (dSigAng-dSig_dCos)/dSigRef ;
   	  // Born+boxes, WARNING Z-box may be modified for KeyZet=2
      double dSigAngF0,dSigAngF1;
    	gps_bornfoam_( 0,m_KFini,m_KFf,m_Mka,m_CosTheta,dSigAngF0);
//    	cout<<" dSigAngF/dSig_dCos = "<< dSigAngF/dSig_dCos;
    	gps_bornfoam_( 1,m_KFini,m_KFf,m_Mka,m_CosTheta,dSigAngF1);
      double dSigAngFF = gps_makerhofoam_(1.0);
      cout<<" // dSigAngFF: "<< (dSigAngFF-dSig_dCos)/dSigRef;
      cout<<"    dSigAngF0: "<< (dSigAngF0-dSig_dCos)/dSigRef;
  	  cout<<" m_CosTheta= "<< m_CosTheta;
      cout<<" m_Mka= "<< m_Mka;
////    	cout<<" m_vv= "<< m_vv;
////    	cout<<" m_uu= "<< m_uu;
      cout<<endl;
    } // end debug **********
//
	double sig0nb = 4*m_pi* sqr(1/m_alfinv)/(3.0*svar2 )*m_gnanob;
//	Dist *=  dSig_dCos *3.0/8.0 *sig0nb;  // Born of EEX
	Dist *=  dSigAng   *3.0/8.0 *sig0nb;  // Born of CEEX

	if( svarCum < sqr(2*m_fin)) Dist = 1e-100;

	return Dist;
}// Density3


///////////////////////////////////////////////////////////////
Double_t KKfoam::Density5(int nDim, Double_t *Xarg)
{ // density distribution for Foam
	m_count++;  // counter for debug
	//
	Double_t Dist=1;

	double svar = sqr(m_CMSene);
	double svarCum = svar;

// ******** mapping for ISR *******
	double gamiCR,gami,alfi;
	double CMSene1= sqrt(svar);
	bornv_makegami_( CMSene1, gamiCR,gami,alfi);   // from KKMC

	double R= Xarg[0];
	m_vv = exp(1.0/gami *log(R)) *m_vvmax; // mapping
	Dist *= m_vv/R/gami ;                  // Jacobian
	if( gami < 0 )      return 0.0;    // temporary fix
    if( m_vv < 1e-200 ) return 0.0;    // temporary fix
    // ISR photonic distribution
	double Rho2 = Rho_isr(svar,m_vv);  // remember take care of m_mbeam!!!
	Dist *= Rho2;
	svarCum *= (1-m_vv);
	double svar2 = svar*(1-m_vv);

// ******** mapping for FSR *******
    double rr= Xarg[1];
    double gamf   = gamFSR(svar2);
    m_uu = exp(1.0/gamf *log(rr));     // mapping
    Dist *= m_uu/rr/gamf;              // Jacobian
	if( gamf < 0 )      return 0.0;    // temporary fix
    if( m_uu < 1e-200 ) return 0.0;    // temporary fix
    // FSR photonic distribution
  	double Rho3 = Rho_fsr(svar2,m_uu);           // remember take care of m_mbeam!!!
  	Dist *= Rho3;
    svarCum *= (1-m_uu);

    // ******** mapping for polar angle *******
    double cmax = 0.999;
    m_CosTheta = cmax*( -1.0 + 2.0* Xarg[2] );
    Dist *= 2.0*cmax;

    // ******** mapping for IFI variaable *******
    double eps =1e-6;
    double gamint = gamIFI(m_CosTheta);
//
    double R1, R2;
    MapIFI( Xarg[3], gamint, eps, m_r1, R1);     // mapping
    MapIFI( Xarg[4], gamint, eps, m_r2, R2);     // mapping
    double RhoIFI1 = Rho_ifi( m_CosTheta, m_r1, eps);
    double RhoIFI2 = Rho_ifi( m_CosTheta, m_r2, eps);
//
//    m_r1 =0;
//    m_r2 =0;
    double WT1 = R1 *RhoIFI1;
    double WT2 = R2 *RhoIFI2;
    Dist *= WT1*WT2;
//
// ******* MC event *******
    double zz = (1-m_vv)*(1-m_uu)*(1-m_r1)*(1-m_r2);
    m_xx = 1-zz;
//    m_xx = m_vv + m_uu - m_vv*m_uu;  // numerically more stable
// effective masses
	m_Mka = sqrt(svar2);   // after ISR

// =============== Sigm/dOmega from spin amplitudes ===============
// Effective 4-momenta, KKMC convention: p={px,py,pz,E)
	double Ene = m_Mka/2;
	double Pmb  = sqrt( (Ene-m_beam)*(Ene+m_beam) ); // modulus
	m_p1 = { 0, 0 , Pmb, Ene};  // beam
	m_p2 = { 0, 0 ,-Pmb, Ene};  // beam
	double Pmf  =sqrt( (Ene-m_fin)*(Ene+m_fin) ); // modulus
	m_p3 = { Pmf*sqrt(1-sqr(m_CosTheta)), 0 , Pmf*m_CosTheta,  Ene}; // final
	m_p4 = {-Pmf*sqrt(1-sqr(m_CosTheta)), 0 ,-Pmf*m_CosTheta,  Ene}; // final
	double PX[4] = {0, 0, 0, 2*Ene};
	double dSigAngF,dSigAngF1,dSigAngF2, Misr1,Misr2;
	Misr1 = sqrt((1-m_vv)*(1-m_r1)*svar);
	Misr2 = sqrt((1-m_vv)*(1-m_r2)*svar);
	gps_bornfoam_( 0,m_KFini,m_KFf,Misr1,m_CosTheta,dSigAngF1);
	gps_bornfoam_( 1,m_KFini,m_KFf,Misr2,m_CosTheta,dSigAngF2);
    dSigAngF = gps_makerhofoam_(1.0);
//************ Debug*** Debug*** Debug*** Debug*** Debug ***********
    if( m_count <1 && fabs(svar/svar2-1)>0.20 ){  // debug
//    if( m_count <1000 ){  // debug
    	double Rat;
    	Rat = dSigAngF1/( dSigAngF2 );
    	cout<<" Density5 debug m_count= "<< m_count<< endl;
    	cout<<" dSigAngF1    = "<< dSigAngF1;
    	cout<<" dSigAngF2    = "<< dSigAngF2;
    	cout<<" svar/svar2 = "<< svar/svar2;
    	cout<<" Rat = "<<Rat<<endl;
    } //
//    if( m_count <10000 && m_r1 > 0 && m_r2 >0 ){  // debug
    if( m_count <1 && m_r1 > 0 && gamint <0 ){  // debug
    	cout<<" Density5 debug m_count= "<< m_count<< endl;
    	cout<<" m_r1= "<< m_r1 <<"  m_r2="<< m_r2<<"  m_xx="<< m_xx <<endl;
    	cout<<" m_CosTheta ="<< m_CosTheta <<" gamint= "<<gamint<<endl;
    	cout<<" WT1 ="<< WT1 <<"  WT2="<< WT2<<endl;
   }
//   **********  end debug **********
	double sig0nb = 4*m_pi* sqr(1/m_alfinv)/(3.0*svar2 )*m_gnanob;
	Dist *=  dSigAngF *3.0/8.0 *sig0nb;

//	if( Dist < 0){
//		cout<< " Density5: Dist= "<< Dist <<endl;
//	}

	if( svarCum < sqr(2*m_fin)) Dist = 1e-100;

	if(m_Mode > 0 ) Dist = fabs(Dist); // For initialization mode

	return Dist;
}// Density5


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                End of Class KKfoam                                        //
///////////////////////////////////////////////////////////////////////////////

