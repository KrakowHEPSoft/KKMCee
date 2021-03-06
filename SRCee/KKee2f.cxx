//               Class   KKee2f                                               //
#include "KKee2f.h"

//#include "Globux.h"

ClassImp(KKee2f);

extern "C" {
// SRChh/ffff_aux.f
   void fort_open_( const int&, const char*, int);
   void fort_close_(const int&);
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
#define BX2F(name,numb,err,text) cout<<"K "<<setw(10)<<name<<\
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
//  m_WtMainMonit = NULL;
  DB         = NULL;
  m_BornDist = NULL;
  m_EWtabs   = NULL;
  m_Event    = NULL;
//  m_GenISR   = NULL;
//  m_GenFSR   = NULL;
//  m_QED3     = NULL;
//  m_GPS      = NULL;
//  m_BVR      = NULL;
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
//  m_WtMainMonit = NULL;
  DB         = NULL;
  m_BornDist = NULL;
  m_EWtabs   = NULL;
  m_Event    = NULL;
//  m_GenISR   = NULL;
//  m_GenFSR   = NULL;
//  m_QED3     = NULL;
//  m_GPS      = NULL;
//  m_BVR      = NULL;
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
//  m_WtMainMonit = new TWtMon("WtMainMonit","WtMain",100,5.0);

// Initialisation input data before generation

  *f_Out<< "   *******************************" << endl;
  *f_Out<< "   ****   KKee2f   Initialize ****" << endl;
  *f_Out<< "   *******************************" << endl;
  cout  << "   *******************************" << endl;
  cout  << "   ****   KKee2f   Initialize ****" << endl;
  cout  << "   *******************************" << endl;

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

  cout<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
  cout<<"%%%%%  f_RNgen->Rndm() = "<< f_RNgen->Rndm()  <<endl;
  cout<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;

  f_FoamI = new TFOAM("FoamI"); // Foam object of base class

  FoamInitA();
//===============

//  double Err, dummy;
//  f_FoamI->GetIntNorm(m_XCrude,Err);


//////////////////////////////////////////////////////////////
// testing object of a new KK template class
  m_KKexamp= new KKlasa(OutFile);
  m_KKexamp->Initialize();

  cout  << "   *******************************" << endl;
  cout  << "   ****   KKee2f   END        ****" << endl;
  cout  << "   *******************************" << endl;

}// Initialize


//-----------------------------------------------------------------------------
double KKee2f::Density(int nDim, double *Xarg){
  double rho= 1e-100;
  //
  //rho = RhoFoam5(Xarg);      // New c++
  // Foam requires density to be positive.
  // This is compensated later on by the weight,
  // but this trick works only for weighted events
  if(m_FoamMode<0 ) rho = abs(rho);
  //
  return rho;
}


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
    nDim = 0;                 // simplices
	m_RhoMode = 5;            // Foam integrand type
    kDim = 5;                 // hyperrectangles:  KF+BES1+BES2+ISR+Theta=5
    f_FoamI->SetnDim(nDim);   // simplicial dimensions
    f_FoamI->SetkDim(kDim);   // hypercubic dimensions
//---------------------------------------------
// set fixed cells for final fermion flavor KFfin
    int nDivL = m_nKF-1;
    double xDivL[nDivL];
    for(int k = 0; k< nDivL; k++) xDivL[k] = ((k+1)*1.0)/nDivL;
    int ibL=0;    // <-- first variable !!!
    f_FoamI->SetInhiDiv(ibL, 1);
    if( nDivL>0 ) f_FoamI->SetXdivPRD(ibL, nDivL, xDivL);
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

void KKee2f::MakeGami(int KFini, double CMSene, double &gamiCR, double &gami, double &alfi){
//////////////////////////////////////////////////////////////////////////////
//   Crude Gami as a function of CMSene                                     //
//////////////////////////////////////////////////////////////////////////////

double amel = DB->fmass[KFini];
double am2  = sqr(2*amel/CMSene);
double chini2 = sqr(DB->Qf[KFini]);  // electric charge of beam fermion
alfi   = chini2 /DB->Alfinv0 /m_pi;
if( am2 < 1 ) {
  double beta = sqrt(1-am2);
  gami    = 2*alfi *( log( sqr(1+beta)/am2) -1);
  gamiCR  = 2*alfi *  log( sqr(1+beta)/am2);
  gamiCR  = gamiCR * DB->Xenph;         // enhancement of crude photon multiplicity
  if(DB->KeyWtm == 1) gamiCR=gami;   // new, for very special tests
} else {
  gamiCR = 0;
  gami   = 0;
}//if
}//MakeGami

double KKee2f::MakeRhoISR(double gamiCR, double gami, double alfi, double vv, double vvmin, double vvmax){
//////////////////////////////////////////////////////////////////////////
// This function calculates the entire RhoISR, including the MakeISR    //
// factor and the Jacobian from RhoFoam together, permitting similar    //
// numerical factors to be canceled analytically. Should be more stable.//
// The factor m_vv must have been generated before calling this.        //
// This function is side-effect free and has no hidden inputs.          //
//////////////////////////////////////////////////////////////////////////
double Rho;
if (vv > vvmin) {
    Rho = (gamiCR/gami)* exp( (gamiCR-gami)*log(vv/vvmin) );
    double DilJac0 = 0.5*(1 + 1/sqrt(1-vv));
    Rho = Rho*DilJac0;
} else {
    Rho = 1;
}
return Rho * exp(gami*log(vvmax ));
}//RhoISR

