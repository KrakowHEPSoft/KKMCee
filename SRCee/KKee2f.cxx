//               Class   KKee2f                                               //
#include "KKee2f.h"

#include "Globux.h"

ClassImp(KKee2f);


extern "C" {
//
   void fort_open_( const int&, const char*, int);
   void fort_close_(const int&);
//   void pseumar_initialize_(const int&, const int&, const int&);
//   void pseumar_makevec_(float rvec[], const int&);
// SUBROUTINE HepEvt_Fill
   void hepevt_fill_();
}//

#define SW20 setw(20)<<setprecision(14)
#define SP15 setw(15)<<setprecision(9)
#define SP10 setw(10)<<setprecision(7)
#define SW10 setw(10)


//FFFFFF  BoX-FORMATs for nice and flexible outputs
#define BOXOPE cout<<\
"KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK"<<endl<<\
"K                                                                              K"<<endl
#define BOXTXT(text) cout<<\
"K                   "<<setw(40)<<         text           <<"                   K"<<endl
#define BOX1I(name,numb,text) cout<<\
"K "<<setw(10)<<name<<" = "<<setw(10)<<numb<<" = "          <<setw(50)<<text<<" K"<<endl
#define BOX1F(name,numb,text)     cout<<"K "<<setw(10)<<name<<\
        " = "<<setw(15)<<setprecision(8)<<numb<<"   =    "<<setw(40)<<text<<" K"<<endl
#define BOX2F(name,numb,err,text) cout<<"K "<<setw(10)<<name<<\
" = "<<setw(15)<<setprecision(8)<<numb<<" +- "<<setw(15)<<setprecision(8)<<err<<\
                                                    "  = "<<setw(25)<<text<<" K"<<endl
#define BOXCLO cout<<\
"K                                                                              K"<<endl<<\
"KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK"<<endl
//FFFFFF  BoX-FORMATs ends here
KKee2f::KKee2f()
{
  // This constructor is for ROOT streamers ONLY
  // All pointers should be NULLed
  cout<< "----> KKee2f Default Constructor (for ROOT only) "<<endl;
  DB         = NULL;
  m_WtMainMonit = NULL;
  m_BornDist = NULL;
  m_EWtabs   = NULL;
  m_Event    = NULL;
  m_GenISR   = NULL;
  m_GenFSR   = NULL;
  m_QED3     = NULL;
  m_GPS      = NULL;
  m_BVR      = NULL;
  m_TauGen   = NULL;
  m_KKexamp  = NULL;
}

///_____________________________________________________________
KKee2f::KKee2f(const char* Name): TMCgen(Name)
//KKee2f::KKee2f(const char* Name)
{
//! all defaults defined here can be changed by the user
//! before calling TMCgen::Initialize
  cout<< "----> KKee2f USER Constructor "<<endl;
///////////////////////////////////////////////////
  DB         = NULL;
  m_WtMainMonit = NULL;
  m_BornDist = NULL;
  m_EWtabs   = NULL;
  m_Event    = NULL;
  m_GenISR   = NULL;
  m_GenFSR   = NULL;
  m_QED3     = NULL;
  m_GPS      = NULL;
  m_BVR      = NULL;
  m_TauGen   = NULL;
  m_KKexamp  = NULL;
  //
  m_EventCounter = 0;
  m_icont = 0;
// Seeds for Pseumar
  m_ijkl_new  = 54217137;
  m_ntot_new  = 0;
  m_ntot2_new = 0;
  m_FoamMode = 0;
  m_RhoMode  = 0;

}///KKee2f

///______________________________________________________________________________________
KKee2f::~KKee2f()
{
  //!Explicit destructor
  cout<< "----> KKee2f::KKee2f !!!! DESTRUCTOR !!!! "<<endl;
}///destructor



///______________________________________________________________________________________
void KKee2f::Initialize(TRandom *RNgen, ofstream *OutFile, TH1D* h_NORMA)
{
  cout<< "----> KKee2f::Initialize, Entering "<<endl;
// Base class assigns  f_Out,  f_RNgen,  f_TMCgen_NORMA
// Foam Object pointer f_FoamI is also defined, not assigned
  TMCgen::Initialize(  RNgen, OutFile, h_NORMA);

//  m_NevTot = 0;
  m_EventCounter = 0;
  m_WtMainMonit = new TWtMon("WtMainMonit","WtMain",100,5.0);

// Initialisation input data before generation

  *f_Out<< "   *******************************" << endl;
  *f_Out<< "   ****   KKee2f   Initialize ****" << endl;
  *f_Out<< "   *******************************" << endl;
  cout  << "   *******************************" << endl;
  cout  << "   ****   KKee2f   Initialize ****" << endl;
  cout  << "   *******************************" << endl;

//////////////////////////////////////////////////////////////
// It is exporting MC generator pointer to global visibility:
  globux_setup_( this );
//////////////////////////////////////////////////////////////


//=============================================================
//   opening disk fime for fortran part of code
  m_out = 16;
  const char *output_file = "./pro77.output";
  int sl2 = strlen(output_file);
  fort_open_(m_out,output_file,sl2);
//  int seed = 54217317; pseumar_initialize_(seed, 0, 0);
//=============================================================

//////////////////////////////////////////////////////////////
// Presently globux is eliminated, may come back temporarily
// It is exporting MC generator pointer to global visibility:
//  globux_setup_( this );
//////////////////////////////////////////////////////////////
  const int jmax = maxPar;
  ReaData("./KKMCee_defaults", jmax, m_xpar);       // f77 indexing in xpar
  ReaData("./pro.input",      -jmax, m_xpar);       // jmax<0 means no-zeroing
  for(int j=0;j<jmax;j++) m_ypar[j]=m_xpar[j+1];    // c++ indexing in ypar

//////////////////////////////////////////////////////////////
//  Data base of static input data
  DB = new KKdbase(OutFile);
  DB->Initialize(  m_xpar );
  m_CMSene =  DB->CMSene;
  m_KFini = 11;

  InitParams();
////////////////////////////////////////////
  cout<<"***** Reading EW tables from DIZET-table1-KK and DIZET-table2-KK *****"<<endl;
  m_EWtabs   = new KKdizet(OutFile); // EW tables from the disk file
  m_EWtabs->Initialize();
  m_EWtabs->ReadEWtabs();    // reads EW tables from the disk file

//  globux_setewtabs_( m_EWtabs);
//////////////////////////////////////////////////////////////
// This replaces BornV class of original KKMC
  m_BornDist = new KKborn(OutFile);
  m_BornDist->SetDB(DB);
  m_BornDist->SetDZ(m_EWtabs);  // EW tables from the disk
  m_BornDist->Initialize();
//////////////////////////////////////////////////////////////
// MC event ISR+FSR record in KKMC format
  m_Event= new KKevent(OutFile);
  m_Event->Initialize(DB->CMSene);

//++++++++++++++++++++++++++++++++++++++++++++++++++
  m_BVR = new KKbvir(OutFile); // Lib of virtual functions
  m_BVR->Initialize();
//++++++++++++++++++++++++++++++++++++++++++++++++++
  m_GenISR= new KKarLud(OutFile);
  m_GenISR->SetDB(DB);
  m_GenISR->SetRNgen(f_RNgen);
  m_GenISR->SetEvent(m_Event);
  m_GenISR->Initialize();

//++++++++++++++++++++++++++++++++++++++++++++++++++
  m_GenFSR= new KKarFin(OutFile);
  m_GenFSR->SetDB(DB);
  m_GenFSR->SetRNgen(f_RNgen);
  m_GenFSR->SetEvent(m_Event);
  m_GenFSR->SetBVR(m_BVR);
  m_GenFSR->Initialize();
//============================
  m_QED3= new KKqed3(OutFile);
  m_QED3->SetDB(DB);
  m_QED3->SetEvent(m_Event);
  m_QED3->SetBornV(m_BornDist);
  m_QED3->Initialize();
//============================
  m_GPS = new KKceex(OutFile);
  m_GPS->SetDB(DB);         // input database
  m_GPS->SetDZ(m_EWtabs);   // EW tables from the disk
  m_GPS->SetEvent(m_Event); // MC event record
  m_GPS->SetBornV(m_BornDist);
  m_GPS->SetBVR(m_BVR);
  m_GPS->SetRNgen(f_RNgen);
  m_GPS->Initialize();


  cout<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
  cout<<"%%%%%  f_RNgen->Rndm() = "<< f_RNgen->Rndm()  <<endl;
  cout<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;

  cout<<"%%%%%  Beam Eenergy Spread (BES) parameters %%%%% "<<endl;
  cout<<" KeyBES = "<< DB->KeyBES<<endl;
  for(int i=0;i<5;i++) m_ParBES[i] = m_xpar[80+i]; // CIRCE params
  cout<<" ParBES[*]=";
  for(int i=0;i<5;i++) cout<<"  "<<m_ParBES[i]; cout<<endl;
  cout<<"%%%%%  Beamstrahlung parameters %%%%% "<<endl;
  for(int i=0;i<4;i++) m_ParCIRCE[i] = m_xpar[76+i]; // CIRCE params
  cout<<" ParCIRCE[*]=";
  for(int i=0;i<4;i++) cout<<"  "<<m_ParCIRCE[i]; cout<<endl;

  f_FoamI = new TFOAM("FoamI"); // Foam object of base class

  FoamInitA();
//===============

//  double Err, dummy;
//  f_FoamI->GetIntNorm(m_XCrude,Err);

//////////////////////////////////////////////////////////////
// testing object of a new KK template class
  m_KKexamp= new KKlasa(OutFile);
  m_KKexamp->Initialize();

//////////////////////////////////////////////////////////////
// TAUOLA+PHOTOS corner
  m_TauGen= new TauPair(OutFile);
  m_TauGen->Initialize(m_ypar);
  m_TauGen->SetEvent(m_Event);
  m_TauGen->SetGPS(  m_GPS);
  m_TauGen->SetRNgen(f_RNgen);
//////////////////////////////////////////////////////////////

  for(int j=1; j<=jmax; j++)  h_NORMA->SetBinContent(j,  m_xpar[j]  );    // xpar encoded
  cout<<" TMCgen::Initialize:  xpar filled into h_NORMA  "<<endl;

  cout  << "   *******************************" << endl;
  cout  << "   **** KKee2f:Initialize END ****" << endl;
  cout  << "   *******************************" << endl;

}// Initialize


//_______________________________________________________________________________
void KKee2f::InitParams(){

for(int j=1; j<=16; j++ ) {
	m_MminCEEX[j] = m_xpar[500+10*j +8];
//	cout<<"KKhh2f::InitParams: j="<<j<<"  m_MminCEEX[j]="<<m_MminCEEX[j]<<endl;
}
//----------------------------------------------------------------------------
//-- The final particle flavor switches in KKMCee are parameters
//--   401 - 406 = d,u,s,c,b,t. Unsupported
//--   411 - 416 = e, nu_e, mu, nu_mu, tau, nu_tau Supported.
//----------------------------------------------------------------------------
cout<< "%%%%%%%%%%%%%% SetFinalStates start %%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
// Only leptonic final states are implemented:
// Two kinds of lepton numbering (temporary)
m_nKF = 0;
for(int KF = 1; KF<= 16; KF++){ // loop over final fermions, 1,2,..16
   if ( m_xpar[400+KF] == 1 ) {
      m_nKF = m_nKF + 1;
      m_KFlist[m_nKF-1] = KF; // c++ indexing
      cout<<"KKee2f will make final state,  KF= " << KF << endl;
   }//if
}//for i
if( m_nKF == 0 ) {
    cout<<"KKee2f::InitParams: No final state fermions were selected."<<endl;
    exit(4);
}//if
cout<< "%%%%%%%%%%%%%% SetFinalStates end %%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
} //InitParams

void KKee2f::FoamInitA(){
	// Foam object initialization
    cout<<"-----------------------FoamInitA start-------------------------------"<<endl;
    int nDim, kDim;
    nDim = 0;                   // simplices
	m_RhoMode = 5;              // Foam integrand type
    kDim = 5;                   // hyperrectangles:  KF+BES1+BES2+ISR+Theta=5
    if(DB->KeyBES == 0) kDim=3;  // BES off
    f_FoamI->SetnDim(nDim);     // simplicial dimensions
    f_FoamI->SetkDim(kDim);     // hypercubic dimensions
//---------------------------------------------
// set fixed cells for final fermion flavor KFfin
    int nDivL = m_nKF-1;
    double xDivL[nDivL];
    for(int k = 0; k< nDivL; k++) xDivL[k] = ((k+1)*1.0)/nDivL;
    int ibL=0;    // <-- first variable !!!
    f_FoamI->SetInhiDiv(ibL, 1);
    if( nDivL>0 ) f_FoamI->SetXdivPRD(ibL, nDivL, xDivL);
    // cos(theta) generation is frozen
    f_FoamI->SetInhiDiv(1, 1);
//----------------------------------------------
    int nCells   = m_xpar[3021]; // move to database?
    int Vopt     = m_xpar[3022];
    int nSampl   = m_xpar[3023];
    int nBins    = m_xpar[3024];
    int EvPerBin = m_xpar[3025];
    int KeyWgt   = m_xpar[3026]; // move to database?
    double WtMaxRej = m_xpar[3027];
    f_FoamI->SetnCells(nCells);     // Maximum number of cells
    f_FoamI->SetOptRej(1-KeyWgt);   // KeyWgt = 0 constant wt, 1 = variable wt
    f_FoamI->SetMaxWtRej(WtMaxRej); // Maximum weight for rejection.
    f_FoamI->SetOptVert(Vopt);      // 0 = Vertices stored, 1 = not stored.
    f_FoamI->SetnSampl(nSampl);     // Number of MC sampling inside single cell
    f_FoamI->SetnBin(nBins);        // Number of bins for edge explorations
    f_FoamI->SetEvPerBin(EvPerBin); // Events per bin during buildup.
    f_FoamI->SetOptEdge(    0);     // OptEdge excludes vertices
    f_FoamI->SetOptDrive(   2);     // Drive = 0, 1, 2 (TrueVol, Sigma, WtMax)
    f_FoamI->SetOptOrd(     0);     // 0: nDim// simplices, 1: single simplex
    f_FoamI->SetOptPeek(    0);     // Choose max. cell (0) or random (1) in build
    f_FoamI->SetOptMCell(   1);     // 1: Megacell = slim memory
    f_FoamI->SetChat(       1);     // printout level =0,1,2

    cout<<"Creating FOAM grid in "<<nDim + kDim<< " dimensions"<<endl;
	m_FoamMode = 1;  // initialization mode
	m_Icont =0;
    f_FoamI->Initialize(f_RNgen,this);

    // Needed to determine nCalls from generation only
    m_nCallsFoam0   = f_FoamI->GetnCalls();
    cout<<"||||| No of density calls in Foam initialization="<< m_nCallsFoam0<<endl;
    //
    m_XsPrim = f_FoamI->GetPrimary();  // Primary (crude) Integral from initialization

    if (KeyWgt == 0) {
        cout<<"FOAM Initialization: fixed weight    "<<endl;
    } else {
        cout<<"FOAM Initialization: variable weight "<<endl;
    }
    cout<<"f_FoamI: simplicial  nDim= "<<nDim<<endl;
    cout<<"f_FoamI: hypercubic  kDim= "<<kDim<<endl;
    cout<<"f_FoamI: m_XsPrim,       = "<<m_XsPrim<<endl;
    //
    cout<<"-----------------------FoamInitA end-------------------------------"<<endl;
}//FoamInitA


//-----------------------------------------------------------------------------
double KKee2f::Density(int nDim, double *Xarg){
  double rho= 1e-100;
  //
  rho = RhoFoam5(Xarg);      // New c++
  // Foam requires density to be positive.
  // This is compensated later on by the weight,
  // but this trick works only for weighted events
  if(m_FoamMode<0 ) rho = abs(rho);
  //
  return rho;
}

double KKee2f::RhoFoam5(double *Xarg){
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  5-dimensional version with generation of KFfin generation                //
//                                                                           //
//  BornSimple is in R-units (pointlike xsection at  sqrt(s)=m_XXXene!       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
int iarg=0;
double Rho = 1.0;
//==========================================================
// Final fermion type
int iKF = 1 + m_nKF*Xarg[iarg]; iarg++;
m_KFfin = m_KFlist[iKF-1];
Rho = Rho *m_nKF; // to get sum over final fermions, not average
// FSR on/off switch for KKarfin and KKceex matrix element
m_Event->m_HasFSR = DB->KeyFSR;        // general FSR switch
if( m_KFfin == 12 || m_KFfin == 12 || m_KFfin == 16)
               m_Event->m_HasFSR  = 0; // exception for neutrinos

//==========================================================
///////////////  Polar angle Theta ///////////////
// it is dummy variable, cos(theta) is generated by KKarlud or KKarfin
double cmax = 0.99999999;
m_CosTheta = cmax*( -1.0 + 2.0* Xarg[iarg] ); iarg++;
//Rho *= 2.0*cmax;  // jacobian
//==========================================================
m_r1=0.0; m_r2=0.0;
int optGauss = 2;   // Gaussian BES with/without mapping
//////////////////  Beam spread ///////////////////////
if(        DB->KeyBES == 0){
  m_Ebeam1  = 0.5*m_CMSene;
  m_Ebeam2  = 0.5*m_CMSene;
}
else if(   DB->KeyBES == 1){
  double E1,E2,sigma1,sigma2,sigma,corho,delE1,delE2,dGauss;
  E1    = m_ParBES[0];
  E2    = m_ParBES[1];
  sigma1= m_ParBES[2]*E1;
  sigma2= m_ParBES[3]*E2;
  corho = m_ParBES[4];
  if( optGauss==1){
   sigma = sqrt(sigma1*sigma2);
   delE1 = 10*sigma*(2*Xarg[iarg]-1.0); iarg++; // range is +-10sigma
   delE2 = 10*sigma*(2*Xarg[iarg]-1.0); iarg++; // range is +-10sigma
   Rho *= sqr(20*sigma);  // Jacobian
   m_r1 = delE1/E1; // can be negative
   m_r2 = delE2/E2; // can be negative
   m_Ebeam1  = E1+delE1; // =E1*(1+m_r1)
   m_Ebeam2  = E2+delE2; // =E2*(1+m_r2)
   dGauss = sqr(delE1/sigma1)+ sqr(delE2/sigma2)-2*corho*(delE1/sigma1)*(delE2/sigma2);
   dGauss = exp(-0.5/(1-sqr(corho))*dGauss);
   dGauss *= 1/(2.0*M_PI)/(sigma1*sigma2)/sqrt(1-sqr(corho)); // Normalization factor
   Rho *= dGauss;
  }else{
   // mapping of Patrick Janot for 2-dim Gaussian BES with optional correlation
   // in this case Jacobian*distribution=1 is omitted.
   double r1 = Xarg[iarg]; iarg++;
   double r2 = Xarg[iarg]; iarg++;
   double x1 = sqrt(-2.*log(r1)) * cos(2.*M_PI*r2);
   double x2 = sqrt(-2.*log(r1)) * sin(2.*M_PI*r2);
   double y1 = x1;
   double y2 = corho * x1 + sqrt(1.-corho*corho) * x2;
   m_r1= y1 * m_ParBES[2]; // can be negative
   m_r2= y2 * m_ParBES[3]; // can be negative
   m_Ebeam1 = E1 * (1.0 + y1 * m_r1);
   m_Ebeam2 = E2 * (1.0 + y2 * m_r2);
  }//if optGauss
} else if( DB->KeyBES == 2){
// Linear Collider case. Beamstrahlung following CIRCE parametrization.
// delta is the "infinitesimal" width of Gaussian representation
// of the Dirac delta in the circe spectrum. It is adjusted empirically.
  double delta = 0.0005;
  m_r1  = Xarg[iarg]; iarg++; // always positive
  m_r2  = Xarg[iarg]; iarg++; // always positive
  m_Ebeam1 = 0.5*m_CMSene*(1-m_r1);
  m_Ebeam2 = 0.5*m_CMSene*(1-m_r2);
  double GS1 = m_ParCIRCE[0]*sqrt( 2/M_PI)/delta *exp (-sqr(m_r1/delta)/2 ); //half of Gauss
  double GS2 = m_ParCIRCE[0]*sqrt( 2/M_PI)/delta *exp (-sqr(m_r2/delta)/2 ); //half of Gauss
  double SF1 = m_ParCIRCE[1]*exp(log(m_r1)*m_ParCIRCE[3]) *exp(log(1-m_r1)*m_ParCIRCE[2]);   // the same as circee
  double SF2 = m_ParCIRCE[1]*exp(log(m_r2)*m_ParCIRCE[3]) *exp(log(1-m_r2)*m_ParCIRCE[2]);   // the same as circee
  double SF12 = (GS1+SF1)*(GS2+SF2);
  Rho *=SF12;
}// if KeyBES
m_XXXene =  2*sqrt(m_Ebeam1*m_Ebeam2);
//--------------------------------------------

///////////////////////  ISR ////////////////////////////////////
// ISR enhancement factor >1 for IFI on in RhoISRold and KKarlud
m_Event->m_Xenph = DB->Xenph;      // default value, >1 for IFI on
if( m_Event->m_HasFSR == 0) m_Event->m_Xenph = 1.0; // no enhancement for FSR off
if( DB->KeyINT == 0)        m_Event->m_Xenph = 1.0; // no enhancement for IFI off
double R,gamiCR,gami,alfi;
double amfin = DB->fmass[m_KFfin];
double vvmax = min(DB->vvmax, 1.0 - sqr(amfin/m_XXXene));
double RhoISR=1, dJac;
if ( DB->KeyISR == 1) {
  MakeGami(m_KFini,m_XXXene,gamiCR,gami,alfi);
  R = Xarg[iarg]; iarg++; // last dimension
  double GamI = gami;
  //GamI = gamiCR; // about the same efficiency
  //Straightforward mapping R->vv
  //m_vv = vvmax* exp((1.0/GamI)*log(R));
  //dJac   = exp(GamI*log(vvmax ))/(GamI/m_vv*exp(GamI*log(m_vv) ));
  //RhoISR *= dJac;
  ///////////////
  // alternative mapping from KKeeFoam
  // probably safer, more stable numerically
  MapPlus( R, GamI, vvmax, m_vv, dJac);
  //[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
  //if( m_EventCounter == 1 ) m_vv=0.4;
  //if( m_EventCounter == 2 ) m_vv=0.2;
  //if( m_EventCounter == 3 ) m_vv=0.6;
  //if( m_EventCounter == 4 ) m_vv=0.8;
  //if( m_EventCounter == 5 ) m_vv=0.1;
  //if( m_EventCounter == 6 ) m_vv=0.3;
  //if( m_EventCounter == 7 ) m_vv=0.5;
  //if( m_EventCounter == 8 ) m_vv=0.7;
  //if( m_EventCounter == 9 ) m_vv=0.0001;
  //if( m_EventCounter >0 && m_EventCounter<= 9 ){
//  int seed = 54217317+10000000*m_vv;
//  int seed = 54217317;
//  pseumar_initialize_(seed, 0, 0);
//  for(int i=0; i<5;i++) *f_Out<< "============================================================================================"<< endl;
//  *f_Out<< "$$$ KKee2f::RhoFoam5: m_vv="<< m_vv<< endl;}
  //]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
  RhoISR *= dJac;
  /////////////// No mapping at all
  //m_vv = R*vvmax;
  //RhoISR *=vvmax;
  ///////////////
  RhoISR *= RhoISRold(m_vv, m_XXXene);
}else {   // ISR off
  m_vv = 0;
  m_Event->m_AvMult     = 0;
  m_Event->m_YFSkon_ini = 1;
  m_Event->m_YFS_IR_ini = 1;
  RhoISR = 1;
}//KeyISR
Rho *= RhoISR;

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Basic variables to be used later on in generation
if( m_FoamMode<0) {  // Only in generation mode !!!!!
  m_Event->m_KFini   = m_KFini;
  m_Event->m_KFfin   = m_KFfin;
  m_Event->m_vv      = m_vv;
  m_Event->m_CosTheta= m_CosTheta;
  m_Event->m_r1      = m_r1;
  m_Event->m_r2      = m_r2;
  m_Event->m_XXXene  = m_XXXene;
}// Event
//------------------
double svar1  = (1-m_vv)*sqr(m_XXXene);      // = (1-m_vv)*x1*x2*sqr(DB->CMSene)
double BornCR  = m_BornDist->BornSimple(m_KFini, m_KFfin, svar1, 0.0);  // costTheta=0
double sig0nb = 4.0*M_PI/( 3.0 *sqr(DB->Alfinv0)) *1.0/(sqr(m_XXXene )) *DB->gnanob;
BornCR  *= sig0nb/(1.0-m_vv);
Rho  *= BornCR;
//
if( std::isnan(Rho)) cout<<"  m_EventCounter="<<m_EventCounter<<" Rho="<<Rho<<endl;
if( abs(Rho) < 1e-150) Rho=1e-150;  // Foam does not tolerate zero integrand
//
return Rho;
}//RhoFoam5


double KKee2f::RhoISRold(double vv, double CMSene){
//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//   This procedure is tightly related to ISR photon generation in Karlud   //
//   It calculates QED crude radiator Function                              //
//                                                                          //
//   m_AvMult is later used in KarLud_YFSini                                //
//   m_YFSkon m_YFS_IR are later used in GPS_Make  and QED3_Make           //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
double Rho,gamiCR,gami,alfi;
MakeGami(m_KFini,CMSene,gamiCR,gami,alfi);
double DilJac0, AvMult,YFS_IR,VoluMC;
double vvmin=DB->vvmin;
if(vv > vvmin) {
   DilJac0   = (1.0+1.0/sqrt(1.0-vv))/2.0;
   AvMult    = gamiCR*log(vv/vvmin);
   VoluMC    = gamiCR/vv *exp( AvMult );   //!!! Phase space Volume CRUDE
   YFS_IR    = -gami*log(1.0/vvmin);       //!!! IR part of YFS formfactor
   Rho       = VoluMC *exp(YFS_IR);
} else {
   DilJac0   = 1.0;
   AvMult    = 0.0;
   VoluMC    = 1.0;
// IMPORTANT:     The integral over Rho(v<vvmin) = YFS_IR = EXP(-gami*LOG(1/vvmin))
   YFS_IR    = -gami*log(1.0/vvmin);        //!!! IR part of YFS formfactor
//   Rho       = 1.0/vv *gami*exp(log(vv)*gami);
   Rho       =  exp(log(vvmin)*gami); // this is for MapPlus
}
Rho =  Rho * DilJac0;
//* YFS formfactor, finite part, YFS_form_Factor = EXP(YFS_IR + YFSkon)
//* YFSkon is delegated/exported to QED3 and GPS (not used here).
m_Event->m_AvMult     = AvMult;
m_Event->m_YFSkon_ini = exp(1/4.0 *gami + alfi*( -0.5  +M_PI*M_PI/3.0) );
m_Event->m_YFS_IR_ini = exp(YFS_IR);
return Rho;
}// RhoISR


void KKee2f::MakeGami(int KFini, double CMSene, double &gamiCR, double &gami, double &alfi){
//////////////////////////////////////////////////////////////////////////////
//   Crude Gami as a function of CMSene                                     //
//////////////////////////////////////////////////////////////////////////////
double amel = DB->fmass[KFini];
double am2  = sqr(2*amel/CMSene);
double chini2 = sqr(DB->Qf[KFini]);  // electric charge of beam fermion
alfi   = chini2 /DB->Alfinv0 /M_PI;
if( am2 < 1 ) {
  double beta = sqrt(1-am2);
  gami    = 2*alfi *( log( sqr(1+beta)/am2) -1);
  gamiCR  = 2*alfi *  log( sqr(1+beta)/am2);
  gamiCR  = gamiCR * m_Event->m_Xenph;      // enhancement of crude photon multiplicity for IFI
  if(DB->KeyWtm == 1) gamiCR=gami;   // new, for very special tests
} else {
  gamiCR = 0;
  gami   = 0;
}//if
}//MakeGami

///------------------------------------------------------------------------
void KKee2f::MapPlus( double r, double gam, double vvmax, double &v, double &dJac){
// Maping for POSITIVE gam
// Input r in (0,1) is uniform random number or FOAM variable
// Returned v is distributed according to gam*v^{gam-1}
  double del = 1e-4;
  double eps =DB->vvmin;
  double Reps,Rat,RV;
  if( fabs(gam*log(eps)) > del){
      Reps = exp(gam*log(eps));
      RV   = exp(gam*log(vvmax));
      Rat  = Reps/RV;
      if( r< Rat ){
          v = 0;  dJac=1/Rat;
      } else {
          v = exp((1/gam)*log(r)); // mapping
          v *= vvmax;
          dJac = RV/(gam/v*exp(gam*log(v)));      // jacobian
      }
  } else {
      Reps = 1+gam*log(eps);
      RV   = 1+gam*log(vvmax);
      Rat  = Reps/RV;
      if( r< Rat ){
          v = 0; dJac=1/Rat;
      } else {
          v = exp(-(1/gam)*(1-r)); // mapping
          v *= vvmax;
          dJac = RV/(gam/v);        // jacobian
      }
  }
  if( v<0 || v>1) {
      cout<<"STOP in TMCgenFOAM::MapPlus: +++ v = "<<v<<endl;
      exit(11);
  }
}// MapPlus

//_______________________________________________________________________________
void KKee2f::Generate()
{
m_EventCounter++;
f_NevGen++;
m_Event->m_EventCounter = m_EventCounter;
////////////////////////////////////////////////
// Foam here starts MC event generation
////////////////////////////////////////////////
m_FoamMode = -1;   // generation mode
////////////////////////////////////////////////
//[[[[[[[[[[[[[[[[[[[[[
//  int seed = 54217317;
//  pseumar_initialize_(seed, 0, 0);
//]]]]]]]]]]]]]]]]]]]]]
f_FoamI->MakeEvent();
f_FoamI->GetMCwt(m_WtFoam);

m_Event->m_nPhotISR=0;
m_Event->m_nPhotFSR=0;
m_Event->m_nPhot=0;

double Mbeam= DB->fmass[m_KFini];  // electron or muon beam mass
m_Event->DefPair(m_XXXene,Mbeam,Mbeam, &(m_Event->m_Pf1), &(m_Event->m_Pf2));

m_Event->m_Rem1.SetPxPyPzE( 0, 0, 0.5*DB->CMSene*m_r1, 0.5*DB->CMSene*m_r1); //obsolete
m_Event->m_Rem2.SetPxPyPzE( 0, 0,-0.5*DB->CMSene*m_r2, 0.5*DB->CMSene*m_r2); //obsolete

m_XXf = m_Event->m_Pf1 + m_Event->m_Pf2;

// Final fermion moment defined, to be possibly overwritten
double the= acos(m_CosTheta);
double phi= 2*M_PI*f_RNgen->Rndm();
double amfi  =DB->fmass[m_KFfin];
m_Event->PhaSpac2(&m_XXf,the,phi,amfi, &(m_Event->m_Qf1), &(m_Event->m_Qf2) );
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
double WT_ISR=1;      // m_GenISR
double WT_FSR=1;      // m_GenISR
m_WtCrude = m_WtFoam;

//////////////////////////////////////
//          ISR   KKarlud           //
//   KeyISR inside m_GenISR????
m_Event->m_Mbeam1 = Mbeam; // input for m_GenISR
m_Event->m_Mbeam2 = Mbeam; // input for m_GenISR
m_GenISR->Make(&m_XXf,&WT_ISR);
m_WtCrude *=  WT_ISR;
//////////////////////////////////////
if(m_WtCrude != 0){
     DB->KeyPia = 0;  // for DEBUG !!!!
     DB->KeyPia = 1;  // default!!!!
     m_GenFSR->Make(&m_XXf,&WT_FSR);
     m_WtCrude *= WT_FSR;
}// if WtCrude
//////////////////////////////////////
if(m_WtCrude!= 0){
    m_Event->Merge(); // merging ISR and FSR photon momenta
}//if wtcrude

//=============================================================
//                    Model weight
// =============================================================
// WtSet reseting to zero
double WtBest, WtBest2;

for(int j=0; j< maxWT;j++){
   m_WtSet[j]   =0;
   m_WtAlter[j] =0;  // formely WtList
}//for

//m_WtMain  =m_WtCrude;
if(m_WtCrude != 0 ) {
// 4-momenta are transfered to QED directly through m_Event
      m_QED3->Make(); // f77 indexing!!!
      for(int j=0; j< maxWT;j++) m_WtSet[j]=m_QED3->m_WtSet[j];
      WtBest = m_WtSet[74];   // wtset[74]
//[[[[[[[[[[[[[[[
      if(m_EventCounter <10 ){
        (*f_Out)<<"KKee2f::Generate: m_WtSet[71-74]="<<m_WtSet[71]<<"  "<<m_WtSet[72]<<"  "<<m_WtSet[73]<<"  "<<m_WtSet[74]<<endl;
        cout    <<"KKee2f::Generate: m_WtSet[71-74]="<<m_WtSet[71]<<"  "<<m_WtSet[72]<<"  "<<m_WtSet[73]<<"  "<<m_WtSet[74]<<endl;
      }
//]]]]]]]]]]]]]]]
// New CEEX matrix element is now default for leptons and for quarks.
// Its use is controled by auxiliary parameter MinMassCEEX variable [GeV]
// CEEX is calculated twice, with ISR*FSR interference OFF and ON
   double SvarQ;
   SvarQ = (m_Event->m_Qf1 + m_Event->m_Qf2)*(m_Event->m_Qf1 + m_Event->m_Qf2);
   if( (DB->KeyGPS != 0) && (SvarQ > sqr(m_MminCEEX[m_KFfin])) ) {
      int KeyInt0=0;
      m_GPS->SetKeyInt(KeyInt0);
      //[[[[[[[[[[ moved out from Make()
      m_GPS->PhelRandom();
      //]]]]]]]]]]
      m_GPS->Make();   // IFI OFF
      m_WtSetNew[51]=m_GPS->m_WtSet[51];
      m_WtSetNew[52]=m_GPS->m_WtSet[52];
      m_WtSetNew[53]=m_GPS->m_WtSet[53];
      WtBest = m_GPS->m_WtSet[53];
      if( DB->KeyINT != 0 && m_Event->m_HasFSR != 0 ) {
//      if( DB->KeyINT != 0 ) {
        m_GPS->SetKeyInt(DB->KeyINT);
        //[[[[[[[[[[
        //m_GPS->PhelRandom();
        //]]]]]]]]]]
        m_GPS->Make(); // IFI ON
        m_WtSetNew[1]=m_GPS->m_WtSet[1];
        m_WtSetNew[2]=m_GPS->m_WtSet[2];
        m_WtSetNew[3]=m_GPS->m_WtSet[3];
        WtBest = m_GPS->m_WtSet[3];
      }// if KeyInt
// +200 shift in the index should be better done inside gps_make() !!?
   for(int j=1; j<200;j++) m_WtSet[j+200] = m_WtSetNew[j];
   }//if
}//if wtcrude
//[[[[[[[[[[[[[[[
if(m_EventCounter <20 ){
  (*f_Out)<<"KKee2f::Generate: m_WtSet[201-203]="<<m_WtSet[201]<<"  "<<m_WtSet[202]<<"  "<<m_WtSet[203]<<endl;
  cout    <<"KKee2f::Generate: m_WtSet[201-203]="<<m_WtSet[201]<<"  "<<m_WtSet[202]<<"  "<<m_WtSet[203]<<endl;
  (*f_Out)<<"KKee2f::Generate: m_WtSet[251-253]="<<m_WtSet[251]<<"  "<<m_WtSet[252]<<"  "<<m_WtSet[253]<<endl;
  cout    <<"KKee2f::Generate: m_WtSet[251-253]="<<m_WtSet[251]<<"  "<<m_WtSet[252]<<"  "<<m_WtSet[253]<<endl;
 }
//]]]]]]]]]]]]]]]

/////////// VARIABLE-WEIGHT events ///////////////
// collection of the weights for the advanced user
// This is version for WTed events only!!!!
   for(int j=1; j< maxWT; j++)   m_WtAlter[j] = m_WtSet[j]*m_WtCrude;  // indexing f77!!!!

/////////////////////////////////////////////////////////////////////////

m_WtMain=WtBest*m_WtCrude;

m_Event->m_WtCrude= m_WtCrude;
m_Event->m_WtMain = m_WtMain;
m_Event->m_WT_FSR = WT_FSR;
m_Event->m_WT_ISR = WT_ISR;
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
m_Event->ZBoostALL();

/////////////////////////////////////////////////////////////
// hadronization interface
if ( m_WtCrude  != 0) hepevt_fill_();

/////////////////////////////////////////////////////////////
// Tau pair generation using TAUOLA+PHOTOS
if ( (m_WtCrude  != 0) && ( m_KFfin == 15) ) {
   int IsTauInitialized = m_TauGen->IsTauInitialized();
   if( IsTauInitialized != 0) {
     if( DB->KeyGPS == 0 ) {
         cout<< " #### STOP in KKhh2f_Generate: for tau decays GPS not activated !!!"<<endl;
         exit(10); }//if
     m_GPS->TralorPrepare(1);
     m_GPS->TralorPrepare(2);
//
     m_TauGen->Make1();
     m_TauGen->ImprintSpin();
     m_TauGen->Make2();
   }//TauIsInitialized
}//if


//////////////////////////////////////////////////////////////////////
// Special histogram used for control of overall normalization
// in case of both weighted and unweighted events
// NevPrim = no of events at the very bottom of all rejection loops
// Note that crude integral XsPrimPb has zero statistical error!
// what about WTmax ???
  double XsPrimPb = f_FoamI->GetPrimary();
  int NevPrim     = f_FoamI->GetnCalls() -m_nCallsFoam0; // Generation only
  f_TMCgen_NORMA->SetBinContent(0,XsPrimPb*NevPrim);  // Picobarns
  f_TMCgen_NORMA->SetEntries(NevPrim);

/////////////////////////////////////////////////////////////
//[[[[[ DEBUG
  if(m_EventCounter <= 20) {
     m_Event->PrintISR();
     m_Event->PrintISR_FSR();
     //m_Event->PrintISR_FSR(f_Out);
     m_Event->EventPrintAll();
     m_Event->EventPrintAll(f_Out);
  }//m_EventCounter
//////////////////////////
// monitoring main WT
  m_WtMainMonit->Fill(m_WtMain);
  //
  if( std::isnan(m_WtCrude) || std::isnan(m_WtMain) ) {
	  cout<<"+++ STOP in KKee2f::Generate:  m_EventCounter="<<m_EventCounter<<endl;
	  cout<<" m_WtCrude="<<m_WtCrude<<"   m_WtMain="<<m_WtMain<<endl;
	  exit(99);
  }// if nan

}// Generate

void KKee2f::ReaData(const char *DiskFile, int imax, double xpar[])
//////////////////////////////////////////////////////////////
//    subprogram reading input data file and packing        //
//    entries into matrix xpar                              //
//    WARNING: input file cannot include empty lines        //
//    it cannot handle entries like 1d-30, has to be 1e-30! //
//////////////////////////////////////////////////////////////
{
  char trail[200];
  char ch1;
  int  foundB=0, line, indx;
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
}// ReaData

///////////////////////////////////////////////////////////////////////////////
void KKee2f::Finalize()
{
*f_Out<< "   *****************************" << endl;
*f_Out<< "   ****   KKMCee   Finalize ****" << endl;
*f_Out<< "   *****************************" << endl;

// True Crude/Primary of Foam inside all rejection loops
double XsPrimPb = f_FoamI->GetPrimary();
////////////////////////////////////////////////////
double IntNorm, Errel;
f_FoamI->Finalize(IntNorm, Errel );      // with printouts
// f_FoamI->GetIntNorm(IntNorm, Errel ); // without printouts
// Warning: IntNorm not always = Primary
cout<<"||||||||| KKMCee::Finalize: IntNorm, Errel="<< IntNorm<<"  "<< Errel<<endl;
double MCresult,MCerror;
f_FoamI->GetIntegMC(MCresult,MCerror);   // true Foam integral
cout<<"||||||||| KKMCee::Finalize: FOAM MCresult, MCerror ="<< MCresult<<" +- "<< MCerror<<endl;
///////////////////////////////////////////////////
fort_close_(m_out);
///////////////////////////////////////////////////
double AveWt, ErrAbs;
m_WtMainMonit->GetAver(AveWt, ErrAbs);
m_XsMainPb   = XsPrimPb*AveWt;
m_XEMainPb   = XsPrimPb*ErrAbs;
BOXOPE;
BOXTXT("****************************************");
BOXTXT("******      KKMCee::Finalize      ******");
BOXTXT("****************************************");
BOX1F(" XsPrimPb",XsPrimPb,     " Crude from Foam [pb]     ");
BOX1F(" <WtMain>",   AveWt,     " average WtMain           ");
BOX1F("       +-",  ErrAbs,     " error abs.               ");
BOX1F("   XsMain",m_XsMainPb,   " Xsection main [pb]       ");
BOX1F("       +-",m_XEMainPb,   " error abs.               ");
BOXCLO;
///////////////////////////////////
BXOPE(*f_Out);
BXTXT(*f_Out,"****************************************");
BXTXT(*f_Out,"******      KKMCee::Finalize      ******");
BXTXT(*f_Out,"****************************************");
BX1F(*f_Out," XsPrimPb",XsPrimPb,     " Crude from Foam [pb]     ");
BX1F(*f_Out," <WtMain>",   AveWt,     " average WtMain           ");
BX1F(*f_Out,"       +-",  ErrAbs,     " error abs.               ");
BX1F(*f_Out,"   XsMain",m_XsMainPb,   " Xsection main [pb]       ");
BX1F(*f_Out,"       +-",m_XEMainPb,   " error abs.               ");
BXCLO(*f_Out);
///////////////////////////////////
// Finilize the LHE output
//  if (DB->KeyLHE == 1 ) {
//      m_Event->LH_Fin();
//  }
}//Finalize

