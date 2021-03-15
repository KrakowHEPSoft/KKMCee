#include "KKeeFoam.h"
#include <stdio.h>

////////////////////////////////////////////////////////////////////////////////
///     KKeeFoam class
/// This is class for axiliary exercises,  mainly integration with Monte Carlo

ClassImp(KKeeFoam);

KKeeFoam::KKeeFoam():
  TMCgen()
{
  /// This constructor is for ROOT streamers ONLY
  cout<< "----> KKeeFoam Default Constructor (for ROOT only) "<<endl;;
  DB         = NULL;             // Database
  m_DZ       = NULL;             // Dizet interface
  m_BornDist = NULL;             // Born differential distribution
  m_Event    = NULL;             // MC event ISR+FSR in KKMC format
  m_BVR      = NULL;             // Library of virtual corrections
  m_GPS      = NULL;             // CEEX matrix element
  m_Foam9    = NULL;             // Foam object
}

///______________________________________________________________________________________
KKeeFoam::~KKeeFoam()
{
  //!Explicit destructor
  cout<< "----> KKeeFoam::KKeeFoam !!!! DESTRUCTOR !!!! "<<endl;
}///destructor

///_____________________________________________________________
KKeeFoam::KKeeFoam(const char* Name):
  TMCgen(Name)
{
//! all defaults defined here can be changed by the user
//! before calling TMCgen::Initialize
  DB         = NULL;             // Database
  m_DZ       = NULL;             // Dizet interface
  m_BornDist = NULL;             // Born differential distribution
  m_Event    = NULL;             // MC event ISR+FSR in KKMC format
  m_BVR      = NULL;             // Library of virtual corrections
  m_GPS      = NULL;             // CEEX matrix element
  m_Foam9    = NULL;             // Foam object
///////////////////////////////////////////////////
///////////////////////////////////////////////////
/// Foam setup
  m_kDim    =    5;         // No. of dim. for Foam, =2,3 Machine energy spread OFF/ON
  m_nCells  = 2000;         // No. of cells, optional, default=2000
  m_nSampl  =  200;         // No. of MC evts/cell in exploration, default=200

  m_del     = 1e-4;         // limit for |gamma*ln(eps)| in IFI mapping
//  m_del     = 1e-6;         // limit for |gamma*ln(eps)| in IFI mapping $$$
  m_eps = 1e-6;             // IR regulator
//  m_eps = 1e-8;             // IR regulator, test $$$
  m_FoamMode    = 7;
///////////////////////////////////////////////////
// debug
  m_count7   =0;
  m_count9   =0;

  m_ceuler  = 0.57721566;

cout<< "----> KKeeFoam::KKeeFoam USER Constructor "<<endl;
}///

///______________________________________________________________________________________
void KKeeFoam::Initialize(TRandom *RNgen, ofstream *OutFile, TH1D* h_NORMA)
{
  cout<< "----> KKeeFoam::Initialize, Entering "<<endl;
  ///	      SETTING UP RANDOM NUMBER GENERATOR
  TMCgen::Initialize(  RNgen, OutFile, h_NORMA);

//////////////////////////////////////////////////////////////
  const int jmax = maxPar;
  ReaData("./KKMChh_defaults", jmax, m_xpar);       // f77 indexing in xpar
  ReaData("./pro.input",      -jmax, m_xpar);       // jmax<0 means no-zeroing
  for(int j=0;j<jmax;j++) m_ypar[j]=m_xpar[j+1];    // c++ indexing in ypar

//////////////////////////////////////////////////////////////
//  Data base of static input data
  DB = new KKdbase(OutFile);
  DB->Initialize(  m_xpar );
//
  m_CMSene = DB->CMSene;
  m_XXXmin = DB->XXXmin; // initial value from input
  m_XXXmax = DB->XXXmax; // initial value from input
  if ((m_XXXmax > DB->CMSene) || (m_XXXmax < m_XXXmin)) m_XXXmax = DB->CMSene;
  m_alfpi  = 1.0/DB->Alfinv0/M_PI;

  double minMfl = m_xpar[656];  // tau mass
  m_Mffmin = 2*minMfl;

  //=============================================================
//   opening disk fime for fortran part of code
  m_out = 16;
  const char *output_file = "./pro77.output";
  int sl2 = strlen(output_file);
  fort_open_(m_out,output_file,sl2);

// initialization of LHAPDF library
  hhpdf_initialize_(m_ypar);
//===========================

////////////////////////////////////////////
  cout<<"***** Reading EW tables from DIZET-table1-KK and DIZET-table2-KK *****"<<endl;
  m_DZ   = new KKdizet(OutFile); // EW tables from the disk file
  m_DZ->Initialize();
  m_DZ->ReadEWtabs();    // reads EW tables from the disk file

//////////////////////////////////////////////////////////////
// This replaces BornV class of original KKMC
  m_BornDist = new KKborn(OutFile);
  m_BornDist->SetDB(DB);
  m_BornDist->SetDZ(m_DZ);  // EW tables from the disk
  m_BornDist->Initialize();

//////////////////////////////////////////////////////////////
// MC event ISR+FSR record in KKMC format
  m_Event= new KKevent(OutFile);
  m_Event->Initialize(DB->CMSene, DB->PDG_H1, DB->PDG_H2);

// Lib of virtual functions
  m_BVR = new KKbvir(OutFile);
  m_BVR->Initialize();

  m_GPS = new KKceex(OutFile);
  m_GPS->SetDB(DB);         // input database
  m_GPS->SetDZ(m_DZ);   // EW tables from the disk
  m_GPS->SetEvent(m_Event); // MC event record
  m_GPS->SetBornV(m_BornDist);
  m_GPS->SetBVR(m_BVR);
  m_GPS->SetRNgen(f_RNgen);
  m_GPS->Initialize();

//////////////////////////////////////////////////////////////
// list of initial state quarks
  m_nQuarks = 0;
  for(int i = 1; i<= 6; i++){
       if(m_xpar[3400+i] == 1) {
          m_nQuarks = m_nQuarks + 1;
          m_QuarkList[m_nQuarks-1] = i;  // KF of the quark
       }//if
 }// for
  // list of final state leptons
  m_nLeptons = 0;
  for(int i = 1; i<= 6; i++){ // loop over leptons, 1,2,..6
     int KF = 10+i;
     if ( m_xpar[400+KF] == 1 ) {
        m_nLeptons = m_nLeptons + 1;
        m_LeptonList[m_nLeptons-1] = KF; // c++ indexing
      }//if
  }//for i


  /////////////////////////////////////////////////////////
  if(f_IsInitialized == 0)
  {
////////////////////////////////////////////////////////////////
/// ********  SETTING UP FOAM object of base class  *******   //
////////////////////////////////////////////////////////////////
  f_FoamI   = new TFOAM("FoamI");   // new instance of MC generator FOAM
  m_kDim    = 7;  // KFi+KFf flavous + z1,z2 of PDFs + ISR + FSR + cosTheta=7
  m_FoamMode    = 7;
  if(DB->KeyISR == 0) m_kDim=m_kDim -1;
  if(DB->KeyFSR == 0) m_kDim=m_kDim -1;
  int nDim=0;
  f_FoamI->SetnDim(nDim);     // simplicial dimensions not used
  f_FoamI->SetkDim(m_kDim);   // hypercubic dimensions
//---------------------------------------------
// set fixed cells for final lepton flavors
  int nDivL = m_nLeptons-1;
  double xDivL[nDivL];
  for(int k = 0; k< nDivL; k++) xDivL[k] = ((k+1)*1.0)/nDivL;
  int ibL=0;    // <-- first variable !!!
  f_FoamI->SetInhiDiv(ibL, 1);
  if( nDivL>0 ) f_FoamI->SetXdivPRD(ibL, nDivL, xDivL);
  // set fixed cells for initial (five) quark flavors
  int nDivQ = m_nQuarks-1;  // q+qbar common 5 fixed bins
  double xDivQ[nDivQ];
  for(int k = 0; k< nDivQ; k++) xDivQ[k] = ((k+1)*1.0)/nDivQ;
  int ibQ=1;     // <--second variable !!!
  f_FoamI->SetInhiDiv(ibQ, 1);
  if( nDivQ>0 ) f_FoamI->SetXdivPRD(ibQ, nDivQ, xDivQ);
//----------------------------------------------
  m_nCells     = m_xpar[3021];
  m_nSampl     = m_xpar[3023];
  int Vopt     = m_xpar[3022];
  int nBins    = m_xpar[3024];
//  int EvPerBin = m_xpar[3025];
  int KeyWgt   = m_xpar[3026];
  double WtMaxRej = m_xpar[3027];
  m_nCells  =   2000; // for test
  m_nSampl  =  20000; // for test
////////////////////////////////////////////////////////////////////
  f_FoamI->SetnCells(m_nCells);     // Maximum number of cells
  f_FoamI->SetnSampl(m_nSampl);     // Number of MC sampling inside single cell
  f_FoamI->SetnBin(        16);     // Number of bins for edge explorations
  f_FoamI->SetOptRej(       0);     // wted events (=0), default
//  f_FoamI->SetMaxWtRej(WtMaxRej); // Maximum weight for rejection.
//  f_FoamI->SetnBin(nBins);        // Number of bins for edge explorations
//  f_FoamI->SetEvPerBin(EvPerBin); // Events per bin during buildup.
//  f_FoamI->SetOptVert(Vopt);      // 0 = Vertices stored, 1 = not stored.
//  f_FoamI->SetOptEdge(    0);     // OptEdge excludes vertices
//  f_FoamI->SetOptDrive(   2);     // Drive = 0, 1, 2 (TrueVol, Sigma, WtMax)
//  f_FoamI->SetOptOrd(     0);     // 0: nDim// simplices, 1: single simplex
//  f_FoamI->SetOptPeek(    0);     // Choose max. cell (0) or random (1) in build
//  f_FoamI->SetOptMCell(   1);     // 1: Megacell = slim memory
//  f_FoamI->SetChat(       1);     // printout level =0,1,2
  f_FoamI->Initialize( f_RNgen, this);  // Initialize FOAM
  double errel;
  f_FoamI->GetIntNorm(m_Xnorm,errel);     // universal normalization
  m_nCallsFoam0   = f_FoamI->GetnCalls(); // Needed for nCalls from generation ONLY
  cout<<"||||| No of density calls in FoamI initialization="<< m_nCallsFoam0<<endl;
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// Additional Foam instance with IFI, two more dimensions
//////////////////////////////////////////////////////////////
  m_Foam9   = new TFOAM("Foam9");   // new instance of MC generator FOAM
  m_kDim    = 9; // 7 + 2 of IFI variables = 9
  if(DB->KeyISR == 0) m_kDim=m_kDim -1;
//  if(DB->KeyFSR == 0) m_kDim=m_kDim -1;
//  if(DB->KeyINT == 0) m_kDim=m_kDim -1;
  m_FoamMode    = 9; // Density function switch
  m_Foam9->SetnDim(nDim);     // simplicial dimensions not used
  m_Foam9->SetkDim(m_kDim);   // hypercubic dimensions
  //---------------------------------------------
  // set fixed cells for final lepton flavors
  nDivL = m_nLeptons-1;
  for(int k = 0; k< nDivL; k++) xDivL[k] = ((k+1)*1.0)/nDivL;
  ibL=0;    // <-- first variable !!!
  m_Foam9->SetInhiDiv(ibL, 1);
  if( nDivL>0 ) m_Foam9->SetXdivPRD(ibL, nDivL, xDivL);
  // set fixed cells for initial (five) quark flavors
  nDivQ = m_nQuarks-1;  // q+qbar common 5 fixed bins
  for(int k = 0; k< nDivQ; k++) xDivQ[k] = ((k+1)*1.0)/nDivQ;
  ibQ=1;     // <--second variable !!!
  m_Foam9->SetInhiDiv(ibQ, 1);
  if( nDivQ>0 ) m_Foam9->SetXdivPRD(ibQ, nDivQ, xDivQ);
//----------------------------------------------
//  m_nCells  =  10000; // for production
//  m_nSampl  = 100000; // for production
  m_Foam9->SetnCells(m_nCells);     // No. of cells, optional, default=2000
  m_Foam9->SetnSampl(m_nSampl);     // No. of MC evts/cell in exploration, default=200
  m_Foam9->SetnBin(        16);     // No. of bins default 8
  m_Foam9->SetOptRej(       0);            // wted events (=0)is default, (=1) for wt=1 events
  m_Foam9->Initialize( f_RNgen, this);     // Initialize FOAM
  m_Foam9->GetIntNorm(m_Xnorm9,errel);     // universal normalization
  m_nCallsFoam9   = m_Foam9->GetnCalls();  // Needed for nCalls from generation ONLY
  //************* special normalization histos  *************
//  m_TMCgen_NORMA9 = TH1D_UP("m_TMCgen_NORMA9","Normalization and xpar",jmax,0.0,10000.0);
  f_HstFile->cd();
  h_TMCgen_NORMA9 = new TH1D("h_TMCgen_NORMA9","Normalization and xpar",jmax,0.0,10000.0);
  for(int j=1; j<=jmax; j++) h_TMCgen_NORMA9->SetBinContent(j,  m_xpar[j]  );    // xpar encoded
  h_TMCgen_NORMA9->SetEntries(0); // Important!!!
  cout<<"||||| No of density calls in Foam9 initialization="<< m_nCallsFoam0<<endl;
//////////////////////////////////////////////////////////////
//screen output
  BXOPE(*f_Out);
  BXTXT(*f_Out,"========================================");
  BXTXT(*f_Out,"======   KKeeFoam::Initialize   ======");
  BXTXT(*f_Out,"========================================");
  BXTXT(*f_Out,"========================================");
  f_IsInitialized = 1;  /// re-initialization inhibited
  } else {
    cout<< "----> KKeeFoam::Initialize, already initialized "<<endl;
  }
}//Initialize


///______________________________________________________________________________________
void KKeeFoam::Generate()
//////////////////////////////////////////////////////////////////////
// Special histogram is used for control of overall normalization
// in case of both weighted and unweighted events
// NevPrim = no of events INSIDE all rejection loops
// NB. Crude integral XsPrimPb has zero statistical error!
//////////////////////////////////////////////////////////////////////
{
if(m_FoamMode == -7) {
  f_NevGen++;
  f_FoamI->MakeEvent();         // Foam of base class
  m_WTfoam   = f_FoamI->GetMCwt();  // get weight
  double XsPrimPb = f_FoamI->GetPrimary();
  int NevPrim     = f_FoamI->GetnCalls() -m_nCallsFoam0; // Generation only
  f_TMCgen_NORMA->SetBinContent(0,XsPrimPb*NevPrim);     // Picobarns
  f_TMCgen_NORMA->SetEntries(NevPrim);
}else if( m_FoamMode == -9){
  f_NevGen++;
  m_Foam9->MakeEvent();         // Foam of base class
  m_WTfoam   = m_Foam9->GetMCwt();  // get weight
  double XsPrimPb = m_Foam9->GetPrimary();
  int NevPrim     = m_Foam9->GetnCalls() -m_nCallsFoam9; // Generation only
  h_TMCgen_NORMA9->SetBinContent(0,XsPrimPb*NevPrim);     // Picobarns
  h_TMCgen_NORMA9->SetEntries(NevPrim);
}else{
  cout<<"+++++ KKeeFoam::Generate: Wrong m_FoamMode = "<<m_FoamMode<<endl; exit(33);
}
  ///////////////////
}//! Generate

///______________________________________________________________________________________
void KKeeFoam::Finalize()
{
  TMCgen::Finalize();
  ///   Finalize MC  run, final printouts, cleaning etc.
  BXOPE(*f_Out);
  BXTXT(*f_Out,"****************************************");
  BXTXT(*f_Out,"******     KKeeFoam::Finalize   ******");
  BXTXT(*f_Out,"****************************************");
  ///------------------------
  Double_t MCresult, MCerror, MCnorm, Errel;
  f_FoamI->Finalize( MCnorm, Errel);  //!
  f_FoamI->GetIntegMC( MCresult, MCerror);  //! get MC integral, should be one
  cout << "**************************************************************"<<endl;
  cout << "**************** KKeeFoam::Finalize  ***********************"<<endl;
  cout << "Directly from FOAM: MCresult= " << MCresult << " +- "<<MCerror <<endl;
  cout << "**************************************************************"<<endl;
  m_Foam9->Finalize( MCnorm, Errel);  //!
  m_Foam9->GetIntegMC( MCresult, MCerror);  //! get MC integral, should be one
  ///------------------------
}//!Finalize


///________________________________________________________________________
double KKeeFoam::Density(int nDim, double *Xarg){
//
  if(        abs(m_FoamMode) == 7 ){
    return Density4(nDim, Xarg);
  } else if( abs(m_FoamMode) == 9 ){
    return Density9(nDim, Xarg);
  } else {
   cout<<" KKeeFoam::Density: wrong Mode ="<<m_FoamMode<<endl;
    exit(-9);
  }
}// Density


////////////////////////////////////////////////////////////////
double KKeeFoam::Density9(int nDim, double *Xarg)
{ // density distribution for Foam
///////////////////////////////////////////////////////////////
m_count9++;
double Rho = 1.0;
int iarg=0;
//--------------------------------------------
// Initial quark type
int Iq = 1 + m_nQuarks*Xarg[iarg]; iarg++;
// factor 2 because q+qbar bins
Rho  = Rho*2*m_nQuarks; // to get sum over quarks, not average
m_KFini = m_QuarkList[Iq-1];
m_chini = DB->Qf[ m_KFini];
m_Mbeam = DB->fmass[m_KFini];  // to be refined, current or constituent?
//==========================================================
// Final lepton type
int Ilep = 1 + m_nLeptons*Xarg[iarg]; iarg++;
m_KFfin = m_LeptonList[Ilep-1];        // !!!!!!
Rho = Rho *m_nLeptons; // to get sum over leptons, not average
m_Mfin   = DB->fmass[m_KFfin];
m_chfin  = DB->Qf[ m_KFfin];
//==========================================================
// PDF variables
double x1, x2, RhoDY;
MaperDY2(Xarg[iarg], Xarg[iarg+1], x1, x2, RhoDY);  // Nalgo=1 onl
Rho  *=RhoDY;
iarg = iarg+2;
m_XXXene = DB->CMSene * sqrt(x1* x2);    // qqbqr energy
double svar = sqr(m_XXXene);             // qqbar system before ISR
// Structure functions
double SF1 = hhpdf_strucfunc_(1,  m_KFini, m_XXXene, x1);
double SF2 = hhpdf_strucfunc_(2, -m_KFini, m_XXXene, x2);
Rho = Rho*SF1*SF2;
//==========================================================
// Polar angle Theta
double cmax = 0.99999;
m_CosTheta = cmax*( -1.0 + 2.0* Xarg[iarg] ); iarg++;
Rho *= 2.0*cmax;  // jacobian
//==========================================================
// ISR
m_vv=0.0;
m_vvmax = std::min(DB->vvmax, 1.0 - sqr(m_Mffmin/m_XXXene));
double svar1 = sqr(m_XXXene);
double RhoISR0=1;
if( DB->KeyISR != 0){
// ******** mapping for ISR *******
double gami = gamISR(svar1);
double dJacISR;
MapPlus(  Xarg[iarg], gami, m_vv, dJacISR); iarg++;
Rho *= dJacISR;  // jacobian
RhoISR0  = RhoISR(0,svar1,m_vv,m_eps); // radiator function
Rho *= RhoISR0;
m_XXXene *= sqrt(1-m_vv);     // Energy after ISR
}
double svar2 = sqr(m_XXXene); // Energy after ISR
//===========================================================
// IFI
//if( DB->KeyINT != 0){ not tested
// ******** mapping for IFI variable *******
double gamint = gamIFI(m_CosTheta);
double dJacInt1, dJacInt2;
MapIFI( Xarg[iarg], gamint, m_r1, dJacInt1); iarg++;  // mapping eps-dependent !!!
MapIFI( Xarg[iarg], gamint, m_r2, dJacInt2); iarg++;  // mapping eps-dependent !!!
double RhoInt1 = RhoIFI( m_CosTheta, m_r1,m_eps);  // implicitly eps-dependent !!!
double RhoInt2 = RhoIFI( m_CosTheta, m_r2,m_eps);  // implicitly eps-dependent !!!
double DistIFI1 = dJacInt1 *RhoInt1;
double DistIFI2 = dJacInt2 *RhoInt2;
Rho *= DistIFI1*DistIFI2;
//}
//===========================================================
//FSR
m_uu=0.0;
double RhoFSR0=1;
//if( DB->KeyFSR != 0){ not tested
// ******** mapping for FSR *******
//double gamf   = gamFSR(svar2);
double gamf   = gamFSR(svar2*(1-m_r1)*(1-m_r2));
double dJacFSR;
MapPlus(  Xarg[iarg], gamf, m_uu, dJacFSR); iarg++;
Rho *= dJacFSR;  // jacobian
RhoFSR0 = RhoFSR(0,svar2,m_uu,m_eps);
Rho *= RhoFSR0;
//}//KeyFSR
m_xx= x1* x2 *(1-m_vv)*(1-m_uu)*(1-m_r1)*(1-m_r2);
//[[[[[[[[[[[[[[[[[[[
// Soft cut-off on QED ISR+FSR+IFI for tests
//double yy =  1.0 - (1-m_vv)*(1-m_uu)*(1-m_r1)*(1-m_r2);
//if(  yy>0.2 ) return 0;
//]]]]]]]]]]]]]]]]]]]
//===========================================================
// Crude distribution
int Nc =3;  // colour
double svarZ = sqr(DB->CMSene) * x1* x2 *(1-m_vv);
double dBornCrude = m_BornDist->BornSimple(m_KFini, m_KFfin, svarZ, 0.0);
dBornCrude *= 1/(2.0*cmax);  // to be compatible with KKMChh (KORALZ) convention
Rho *= dBornCrude/Nc;
double sig0nb = 4*M_PI* 1.0/(3.0*sqr(DB->Alfinv0)) *1.0/(svarZ )*DB->gnanob;
Rho *=  sig0nb;
Rho *=  1/(x1*x2);    // PDFs normalized as momentum distribution instead luminosities
if( svarZ*(1-m_uu) < sqr(2*m_Mfin)) Rho = 1e-100;
///////////////////////////////////////////////////////////////////////////
//// The basic/crude integrand Rho for Foam is complete at this point !! //
//// Weight for the alternative distribution are next            //////////
///////////////////////////////////////////////////////////////////////////
SetEvent( svarZ, m_CosTheta); // set input for BornFoam0
m_GPS->m_KeyInt=0;
double dSigFoam0 = m_GPS->BornFoam0(   m_KFini, m_KFfin, svarZ, m_CosTheta); // CEEX Born
//======================================================================
//!!!!!!!!!!!!!!! NEW !!!!!!!!!!!!!!!!!!
double Yint, sisr1,sisr2;
sisr1 = (1-m_vv)*(1-m_r1)*svar; // ISR + one IFI variable
sisr2 = (1-m_vv)*(1-m_r2)*svar; // ISR + one IFI variable
// Re(M M^*) including only leading part of gamma-Z box
m_GPS->m_KeyInt=2;              // important!!!
SetEvent( sisr1, m_CosTheta); // set input for BornFoam2
double Emin = 0.5 * m_XXXene * DB->vvmin;
m_GPS->SetEmin(Emin);                      // not needed!!
m_GPS->BornFoam2( 10,m_KFini,m_KFfin,sisr1,m_CosTheta,Yint);
SetEvent( sisr2, m_CosTheta); // set input for BornFoam2
m_GPS->BornFoam2( 11,m_KFini,m_KFfin,sisr2,m_CosTheta,Yint);
double dSigFoam2 = m_GPS->MakeRhoFoam() *Yint;
//[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
//if( f_NevGen>0 && f_NevGen<500) {
//cout<<"KKeeFoam::Density9: gamint="<< gamint<<" vv="<<m_vv<<" uu="<<m_uu<<" r1="<<m_r1<<" r2="<<m_r2<<endl;
//cout<<"KKeeFoam::Density9:  sisr1="<<sisr1<<" sisr2="<<sisr2<<" Yint="<<Yint<<sisr2<<" dSigFoam2="<<dSigFoam2<<endl;
//}
//]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//======================================================================
// alternative weights for CEEX
double wt0  = 3.0/8.0* dSigFoam2/dBornCrude; // for CEEX0 with IFI
//double wt0  = 3.0/8.0* dSigFoam0/dBornCrude; // for CEEX0 with IFI
m_WTset[201]   = wt0; // CEEX0, includes RhoISR(0,...) and/or RhoFSR(0,...)
m_WTset[202]   = wt0; // CEEX1
m_WTset[203]   = wt0; // CEEX1
if(DB->KeyISR != 0) m_WTset[202] *= RhoISR(1,svar1,m_vv,m_eps)/RhoISR0; // ISR CEEX1
if(DB->KeyISR != 0) m_WTset[203] *= RhoISR(2,svar1,m_vv,m_eps)/RhoISR0; // ISR CEEX2
if(DB->KeyFSR != 0) m_WTset[202] *= RhoFSR(1,svar2,m_uu,m_eps)/RhoFSR0; // FSR CEEX1
if(DB->KeyFSR != 0) m_WTset[203] *= RhoFSR(2,svar2,m_uu,m_eps)/RhoFSR0; // FSR CEEX2
//////////////////////////////////////////////////////////////
// random swap of q and qbar
if ( f_RNgen->Rndm() < 0.5) m_AntiQ = 0; else m_AntiQ = 1;
if( m_AntiQ) {double x0; x0=x1; x1=x2; x2=x0;} // swap x1,x2 of PDFs
if( m_AntiQ) m_CosTheta *= -1.0; // and reverse z=axis for final fermions
/////////////////////////////////////////////////////////////
//  Finishing kinematics
m_y1 = 1.0-x1;
m_y2 = 1.0-x2;
if( (m_y1 == 1) || (m_y2 == 1) ) return 0; // implausible but safe
/////////////////////////////////////////////////////////////
if(m_FoamMode > 0 ) Rho = fabs(Rho); // For initialization mode
return Rho;
}//Density9

//_____________________________________________________________________
double KKeeFoam::Density4(int nDim, double *Xarg)
{ // density distribution for Foam
///////////////////////////////////////////////////////////////
m_count7++;
double Rho = 1.0;
int iarg=0;
//--------------------------------------------
// Initial quark type
int Iq = 1 + m_nQuarks*Xarg[iarg]; iarg++;
// factor 2 because q+qbar bins
Rho  = Rho*2*m_nQuarks; // to get sum over quarks, not average
m_KFini = m_QuarkList[Iq-1];
m_chini = DB->Qf[ m_KFini];
m_Mbeam = DB->fmass[m_KFini];  // to be refined, current or constituent?
//--------------------------------------------
// Final lepton type
int Ilep = 1 + m_nLeptons*Xarg[iarg]; iarg++;
m_KFfin = m_LeptonList[Ilep-1];        // !!!!!!
Rho = Rho *m_nLeptons; // to get sum over leptons, not average
m_Mfin   = DB->fmass[m_KFfin];
m_chfin  = DB->Qf[ m_KFfin];
//----------------------------------------------
double x1, x2, RhoDY;
MaperDY2(Xarg[iarg], Xarg[iarg+1], x1, x2, RhoDY);  // Nalgo=1 onl
Rho  *=RhoDY;
//[[[[[[[[[[[[ debug
//x1 = Xarg[iarg]; x2 = Xarg[iarg+1];
//double QQ= DB->CMSene*sqrt(x1*x2);
//if(QQ< m_XXXmin || QQ>m_XXXmax) Rho=1e-100;
//]]]]]]]]]]]]

iarg = iarg+2;
m_XXXene = DB->CMSene * sqrt(x1* x2);

// Structure functions
double SF1 = hhpdf_strucfunc_(1,  m_KFini, m_XXXene, x1);
double SF2 = hhpdf_strucfunc_(2, -m_KFini, m_XXXene, x2);
//[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
//      cout<< "m_KFini, m_XXXene, x1, x2 = "<< m_KFini<<"  "<< m_XXXene<<"  "<< x1<<"  "<< x2<< endl;
//      cout<< "SF1, SF2 = "<< SF1<<"  "<< SF2<< endl;
//]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
Rho = Rho*SF1*SF2;

// ******** mapping for polar angle *******
double cmax = 0.99999;
m_CosTheta = cmax*( -1.0 + 2.0* Xarg[iarg] ); iarg++;
Rho *= 2.0*cmax;  // jacobian

m_vv=0.0;
m_vvmax = std::min(DB->vvmax, 1.0 - sqr(m_Mffmin/m_XXXene));
double svar1 = sqr(m_XXXene);
double RhoISR0;
if( DB->KeyISR != 0){
// ******** mapping for ISR *******
  double gami = gamISR(svar1);
  double dJacISR;
  MapPlus(  Xarg[iarg], gami, m_vv, dJacISR); iarg++;
  Rho *= dJacISR;  // jacobian
  RhoISR0  = RhoISR(0,svar1,m_vv,m_eps); // radiator function
//  double RhoISR0  = RhoISR(2,svar1,m_vv,m_eps); // radiator function
  Rho *= RhoISR0;
  m_XXXene *= sqrt(1-m_vv);
}// if KeyISR

m_uu=0.0;
double svar2 = sqr(m_XXXene);
double RhoFSR0;
if( DB->KeyFSR != 0){
// ******** mapping for FSR *******
  double gamf   = gamFSR(svar2);
  double dJacFSR;
  MapPlus(  Xarg[iarg], gamf, m_uu, dJacFSR); iarg++;
  Rho *= dJacFSR;  // jacobian
  RhoFSR0 = RhoFSR(0,svar2,m_uu,m_eps);
//  double RhoFsr0 = RhoFSR(2, svar2,m_uu,m_eps);
  Rho *= RhoFSR0;
//[[[[[[[[[[[[[[[[[
//  Simplified version, works approximately the same
//  gamf = sqr(m_chfin)*2*m_alfpi*( log(svar2/sqr(m_Mfin)) -1);
//  double rr = Xarg[iarg]; iarg++;
//  m_uu = exp(1.0/gamf *log(rr));     // mapping
//  Rho *= m_uu/rr/gamf;              // Jacobian
//  if( gamf < 0 )      return 0.0;    // temporary fix
//  if( m_uu < 1e-200 ) return 0.0;    // temporary fix
//  Rho *= gamf *exp(gamf*log(m_uu))/m_uu; // FSR distribution
//]]]]]]]]]]]]]]]]]
}// if KeyISR
m_xx= x1* x2 *(1-m_vv)*(1-m_uu);

int Nc =3;  // colour
double svarZ = sqr(DB->CMSene) * x1* x2 *(1-m_vv);
//double svarZ = sqr(m_XXXene); // the same
double dBornCrude = m_BornDist->BornSimple(m_KFini, m_KFfin, svarZ, 0.0);
dBornCrude *= 1/(2.0*cmax);  // to be compatible with KKMChh (KORALZ) convention
//[[[[[[[
//double dBornCrude = m_BornDist->BornSimple(m_KFini, m_KFfin, svarZ, m_CosTheta);
//dBornCrude *= 3.0/8.0;      // 3/8 corrects for KORALZ convention in BornSimple
Rho *= dBornCrude/Nc;
double sig0nb = 4*M_PI* 1.0/(3.0*sqr(DB->Alfinv0)) *1.0/(svarZ )*DB->gnanob;
Rho *=  sig0nb;
Rho *=  1/(x1*x2);    // PDFs normalized as momentum distribution instead luminosities
if( svarZ*(1-m_uu) < sqr(2*m_Mfin)) Rho = 1e-100;
//[[[[[[[[[[[[[[[[[[[
// Temporary cut-off on QED ISR+FSR+IFI for tests
//double yy =  1.0 - (1-m_vv)*(1-m_uu);
//if(  yy>0.2 ) return 0;
//]]]]]]]]]]]]]]]]]]]
////=====================================================================//
//// The basic/crude integrand Rho for Foam is complete at this point !! //
//// Weight for the alternative distribution are next            //////////
///////////////////////////////////////////////////////////////////////////
SetEvent( svarZ, m_CosTheta);
double dSigSimple  = m_BornDist->BornSimple(  m_KFini, m_KFfin, svarZ, m_CosTheta); // Crude(theta)
double dSigDizetS  = m_BornDist->Born_DizetS( m_KFini, m_KFfin, svarZ, m_CosTheta); // EEX
m_GPS->m_KeyInt=0;
double dSigFoam0   =      m_GPS->BornFoam0(   m_KFini, m_KFfin, svarZ, m_CosTheta); // CEEX
//[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[******************
//if( f_NevGen>0 && f_NevGen<500 && m_KFini>2 && m_KFfin==15) {
if( f_NevGen>0 && f_NevGen<30) {
    double ratio =  dSigFoam0/dSigSimple;
    double ratio2 = dSigFoam0/dSigDizetS;
//	cout<<"KKeeFoam::Density4: *** dSigFoam0="<<dSigFoam0<<" dSigSimple"<< dSigSimple    <<" CEEX/EEX="<<ratio<< endl;
	cout<<"KKeeFoam::Density4: KFini="<< m_KFini<<" KFfin="<<m_KFfin;
	cout<<"KKeeFoam::Density4:  Mll= "<<sqrt(svarZ)<<" CosTheta= "<< m_CosTheta<<" Foam0/Simple= "<<ratio<<" Foam0/DizetS= "<<ratio2<< endl;
}
//]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]******************
/////////////////////////////////////////////////////////////
// alternative weights for EEX
m_wt0= 3.0/8.0* dSigDizetS/dBornCrude; // mainly for EEX0, also for keyISR=0, keyFSR=0
m_WTset[ 71]  = m_wt0; // EEX0, includes RhoISR(0,...) and/or RhoFSR(0,...)
m_WTset[ 72]  = m_wt0; // EEX1
m_WTset[ 73]  = m_wt0; // EEX2
if(DB->KeyISR != 0) m_WTset[ 72] *= RhoISR(1,svar1,m_vv,m_eps)/RhoISR0; // ISR EEX1
if(DB->KeyISR != 0) m_WTset[ 73] *= RhoISR(2,svar1,m_vv,m_eps)/RhoISR0; // ISR EEX2
if(DB->KeyFSR != 0) m_WTset[ 72] *= RhoFSR(1,svar1,m_uu,m_eps)/RhoFSR0; // FSR EEX1
if(DB->KeyFSR != 0) m_WTset[ 73] *= RhoFSR(2,svar2,m_uu,m_eps)/RhoFSR0; // FSR EEX2
// alternative weights for CEEX
double wt0  = 3.0/8.0* dSigFoam0/dBornCrude; // for CEEX0
m_WTset[251]   = wt0; // CEEX0, includes RhoISR(0,...) and/or RhoFSR(0,...)
m_WTset[252]   = wt0; // CEEX1
m_WTset[253]   = wt0; // CEEX1
if(DB->KeyISR != 0) m_WTset[252] *= RhoISR(1,svar1,m_vv,m_eps)/RhoISR0; // ISR CEEX1
if(DB->KeyISR != 0) m_WTset[253] *= RhoISR(2,svar1,m_vv,m_eps)/RhoISR0; // ISR CEEX2
if(DB->KeyFSR != 0) m_WTset[252] *= RhoFSR(1,svar2,m_uu,m_eps)/RhoFSR0; // FSR CEEX1
if(DB->KeyFSR != 0) m_WTset[253] *= RhoFSR(2,svar2,m_uu,m_eps)/RhoFSR0; // FSR CEEX2
//////////////////////////////////////////////////////////////
// random swap of q and qbar
if ( f_RNgen->Rndm() < 0.5) m_AntiQ = 0; else m_AntiQ = 1;
if( m_AntiQ) {double x0; x0=x1; x1=x2; x2=x0;} // swap x1,x2 of PDFs
if( m_AntiQ) m_CosTheta *= -1.0; // and reverse z=axis for final fermions
/////////////////////////////////////////////////////////////
//  Finishing kinematics
m_y1 = 1.0-x1;
m_y2 = 1.0-x2;
if( (m_y1 == 1) || (m_y2 == 1) ) return 0; // implausible but safe
/////////////////////////////////////////////////////////////
if(m_FoamMode > 0 ) Rho = fabs(Rho); // For initialization mode
return Rho;
}// Density4

void KKeeFoam::SetEvent( double svarZ, double CosTheta){
/////////////////////////////////////////////////////////////
// Set fermion momenta in Z frame to be used in CEEX Born
double fleps = 1e-20;
double Eini = sqrt(svarZ)/2.0;
double Pini = sqrt((Eini-m_Mbeam)*(Eini+m_Mbeam));
double Pfin = sqrt((Eini-m_Mfin) *(Eini+m_Mfin));
double sinTheta = sqrt(1-sqr(CosTheta));
m_Event->m_Pf1.SetPxPyPzE( 0, 0, Pini, Eini);
m_Event->m_Pf2.SetPxPyPzE( 0, 0,-Pini, Eini);
m_Event->m_Qf1.SetPxPyPzE( 0, sinTheta*Pfin, CosTheta*Pfin, Eini);
m_Event->m_Qf2.SetPxPyPzE( 0,-sinTheta*Pfin,-CosTheta*Pfin, Eini);
int hel1 =1, hel2= -1, hel3 = 1, hel4=1;
m_GPS->m_p1.SetAll(  1, hel1, m_Mbeam, m_Event->m_Pf1);  // e-    U-spinor
m_GPS->m_p2.SetAll( -1, hel2, m_Mbeam, m_Event->m_Pf2);  // e+    V-bar spinor
m_GPS->m_p3.SetAll(  1, hel3, m_Mfin,  m_Event->m_Qf1);  // f     U-bar-spinor
m_GPS->m_p4.SetAll( -1, hel4, m_Mfin,  m_Event->m_Qf2);  // f-bar V-spinor
// Special version of initial momenta with infinitesimal mass for Born spinors
m_GPS->m_r1.SetAll(  1, hel1, fleps,   m_Event->m_Pf1);  // e-    U-spinor
m_GPS->m_r2.SetAll( -1, hel2, fleps,   m_Event->m_Pf2);  // e+    V-bar spinor
}//SetEvent

///------------------------------------------------------------------------
double KKeeFoam::RhoISR(int KeyISR, double svar, double vv, double eps){
/// ISR rho-function for ISR
  double alf1   = m_alfpi;
  double gami   = gamISR(svar);
  double gamfac = Fyfs(gami);
  double delb   = gami/4 +alf1*(-0.5  +sqr(M_PI)/3.0);
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
	 cout<<"+++++TMCgenFOAM::KKdistr: Wrong KeyISR = "<<KeyISR<<endl; exit(5);
  }
  if( vv > eps){
     rho = ffact*gami* exp( log(vv)*(gami-1) ) *(1 +dels +delh);
  }else{
	 rho = ffact* exp( log(eps)*(gami) ) *(1 +dels);
  }
  return rho;
}//Rho_ISR


///------------------------------------------------------------------------
double KKeeFoam::RhoFSR(int KeyFSR, double svar, double uu, double eps){
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
  double delb   = gamf/4 +alf1*(-0.5  +sqr(M_PI)/3.0)
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


///--------------------------------------------------------------
double KKeeFoam::RhoIFI(double costhe, double uu, double eps){
/// ISR+FSR rho-function amlitute level
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


///------------------------------------------------------------------------
double KKeeFoam::Fyfs(double gam){
	  return exp(-m_ceuler*gam)/TMath::Gamma(1+gam);
}

///------------------------------------------------------------------------
double KKeeFoam::gamISR( double svar){
	  return  sqr(m_chini)*2*m_alfpi*( log(svar/sqr(m_Mbeam)) -1);
}

///------------------------------------------------------------------------
double KKeeFoam::gamFSR( double svar){
	  return  sqr(m_chfin)*2*m_alfpi*( log(svar/sqr(m_Mfin)) -1);
}

///------------------------------------------------------------------------
double KKeeFoam::gamIFI( double costhe){
	  return  m_chini*m_chfin*2*m_alfpi* log( (1-costhe)/(1+costhe));
}


void KKeeFoam::MaperDY2(double r1, double r2, double &x1, double &x2, double &Rho){
//----------------------------------------------------------------------------
// Maps r1, r2 into parton momentum fractions x1, x2,
// power map for Q^2, Nalgo = 1 only.
//----------------------------------------------------------------------------
  double svar = sqr(DB->CMSene);
  double Qsqmin = sqr(m_XXXmin);
  double Qsqmax = sqr(m_XXXmax);
  double Qsq = 1/((1-r1)/Qsqmin + r1/Qsqmax);
  Rho = sqr(Qsq) * (1/Qsqmin - 1/Qsqmax);
  double z = Qsq/svar;
  x1 = exp(r2*log(z));
  x2 = z/x1;
  Rho = -log(z) * Rho/svar;
  }// MaperDY2


///------------------------------------------------------------------------
void KKeeFoam::MapPlus( double r, double gam, double &v, double &dJac){
// Maping for POSITIVE gam
// Input r in (0,1) is uniform random number or FOAM variable
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
void KKeeFoam::MapMinus( double r, double gam, double &v, double &dJac){
// Maping for NEGATIVE gam
// Input r in (0,1) is uniform random number of FOAM variable
// Returned v is distributed according to gam*v^{gam-1}
// dJac is normalization (part of Jacobian) factor
  double eps = m_eps;
  double Rat, Reps, RV, R1;
  if( fabs(gam*log(eps)) > m_del){
      Reps = exp(gam*log(eps)); // R(eps)
//         R1   = 2*Reps-1;          // R(1)
      RV   = 2*Reps-exp(gam*log(m_vvmax)); // R(vmax)
      Rat  = Reps/RV;
      if( r< Rat ){
          v = 0; dJac= 1/Rat;
      } else {
          v    = exp( (1/gam)*log( 2*Reps -RV*r ) ); // mapping
          dJac = RV/(-gam/v*exp(gam*log(v)));        // jacobian
//        cout<<"MapMinus://////////////////v="<<v<<endl;
      }
  } else {
      Reps = 1+gam*log(eps);        // R(eps)
      R1   = 1+2*gam*log(eps);      // R(1)
////      RV   = R1 +gam*log(m_vvmax);  // R(V) // Error!!! for r=1 implies v=1/vvmax !!!
      RV   = R1 -gam*log(m_vvmax);  // R(V)
      Rat  = Reps/RV;
      if( r< Rat){
          v = 0; dJac= 1/Rat;
      } else {
          v    = exp( (1/gam)*(R1 -r*RV) ); // mapping
          dJac = RV/(-gam/v);               // jacobian
//        cout<<"MapMinus:!!!!!!!!!!!!!!!!!v="<<v<<endl;
      }
  }
  //[[[[[[[[[[[[[[[[[[[
  //if(m_count9<1000) cout<<"MapMinus:gam="<<gam<<" eps="<<eps<<" del="<<m_del<<" vvmax="<<m_vvmax<<" v="<<v<<endl;

  if( v<0 || v>1) {
      cout<<"STOP in KKeeFoam::MapMinus: +++ v = "<<v<<endl;
      cout<<" m_vvmax = "<<m_vvmax<<" gam = "<<gam<<" r= "<<r<<endl;
      cout<<" Rat = "<<Rat<<" Reps = "<<Reps<<" RV= "<<RV<<" R1= "<<R1<<endl;
      cout<<" m_del = "<<m_del<<" m_eps = "<<m_eps<<endl;
      exit(11);
  }
}// MapMinus

///------------------------------------------------------------------------
void KKeeFoam::MapIFI( double r, double gam, double &v, double &R){
//// mapping for IFI
if(gam > 0)
    MapPlus(  r, gam, v, R);
else
    MapMinus( r, gam, v, R);
}// MapIFI


double KKeeFoam::Sig0nb(double &CMSene){
//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  provides pointlike muon x-section in nanobarns                          //
//  for normalization purpose                                               //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
//Sig0nb =  4d0*pi/(m_AlfInv**2*3d0*CMSene**2)*m_gnanob
return 4*M_PI/( sqr(DB->Alfinv0) *3 *sqr(CMSene ) ) *DB->gnanob;
}//Sig0nb


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                             UTILITIES                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

void KKeeFoam::ReaData(const char *DiskFile, int imax, double xpar[])
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

