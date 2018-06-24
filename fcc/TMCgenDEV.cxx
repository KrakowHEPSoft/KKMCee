#include "TMCgenDEV.h"
#include <stdio.h>

////////////////////////////////////////////////////////////////////////////////
///     TMCgenDEV class
/// This is class for axiliary exercises,  mainly integration with Monte Carlo

ClassImp(TMCgenDEV);

TMCgenDEV::TMCgenDEV():
  TMCgen()
{
  /// This constructor is for ROOT streamers ONLY
  cout<< "----> TMCgenDEV Default Constructor (for ROOT only) "<<endl;
}

///______________________________________________________________________________________
TMCgenDEV::~TMCgenDEV()
{
  //!Explicit destructor
  cout<< "----> TMCgenDEV::TMCgenDEV !!!! DESTRUCTOR !!!! "<<endl;
}///destructor

///_____________________________________________________________
TMCgenDEV::TMCgenDEV(const char* Name):
  TMCgen(Name)
{
//! all defaults defined here can be changed by the user
//! before calling TMCgen::Initialize
///////////////////////////////////////////////////
/// Foam setup
  m_kDim    =    5;         // No. of dim. for Foam, =2,3 Machine energy spread OFF/ON
  m_nCells  = 2000;         // No. of cells, optional, default=2000
  m_nSampl  =  200;         // No. of MC evts/cell in exploration, default=200
///////////////////////////////////////////////////
// debug
  m_count   =0;

cout<< "----> TMCgenDEV::TMCgenDEV USER Constructor "<<endl;
}///

///______________________________________________________________________________________
void TMCgenDEV::Initialize(TRandom *RNgen, ofstream *OutFile, TH1D* h_NORMA)
{
  cout<< "----> TMCgenDEV::Initialize, Entering "<<endl;
  ///	      SETTING UP RANDOM NUMBER GENERATOR
  TMCgen::Initialize(  RNgen, OutFile, h_NORMA);

/////////////////////////////////////////////////////////
/*
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
  */

  double errel;
  /////////////////////////////////////////////////////////
  f_FoamI   = new TFOAM("FoamI");   // new instance of MC generator FOAM
  m_kDim    = 2;
  m_nCells  =  10000;
  m_nSampl  = 100000;
  f_FoamI->SetkDim(m_kDim);         // No. of dims. Obligatory!
  f_FoamI->SetnCells(m_nCells);     // No. of cells, optional, default=2000
  f_FoamI->SetnSampl(m_nSampl);     // No. of MC evts/cell in exploration, default=200
  f_FoamI->SetnBin(        16);     // No. of bins default 8
  f_FoamI->SetOptRej(0);            // wted events (=0), default wt=1 events (=1)

  f_FoamI->Initialize( f_RNgen, this);     // Initialize FOAM
  f_FoamI->GetIntNorm(m_Xnorm,errel);   // universal normalization
 }/// Initialize

///______________________________________________________________________________________
void TMCgenDEV::Generate()
{
  f_NevGen++;
  f_FoamI->MakeEvent();         // Foam of base class
  m_WT   = f_FoamI->GetMCwt();  // get weight
  f_TMCgen_NORMA->Fill(-1, m_Xnorm);    // New style

}//! Generate

///______________________________________________________________________________________
void TMCgenDEV::Finalize()
{
  TMCgen::Finalize();
  ///   Finalize MC  run, final printouts, cleaning etc.
  BXOPE(*f_Out);
  BXTXT(*f_Out,"****************************************");
  BXTXT(*f_Out,"******     TMCgenDEV::Finalize   ******");
  BXTXT(*f_Out,"****************************************");
  ///-----------------------
    Double_t MCresult, MCerror, MCnorm, Errel;
    f_FoamI->Finalize( MCnorm, Errel);  //!
    f_FoamI->GetIntegMC( MCresult, MCerror);  //! get MC integral, should be one
    cout << "**************************************************************"<<endl;
    cout << "**************** TMCgenDEV::Finalize  ***********************"<<endl;
    cout << "Directly from FOAM: MCresult= " << MCresult << " +- "<<MCerror <<endl;
    cout << "**************************************************************"<<endl;
  ///------------------------
}//!Finalize



///________________________________________________________________________
double TMCgenDEV::Density(int nDim, double *Xarg){
	//
	m_count++;  // counter for debug
	double R= Xarg[0];
	double S= Xarg[1];

	double Dist = R*S;

	return Dist;

}// Density
