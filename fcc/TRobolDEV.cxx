#include "TRobolDEV.h"


////////////////////////////////////////////////////////////////////////////////
/// MC generator class for Two Parton exercises

ClassImp(TRobolDEV);

TRobolDEV::TRobolDEV():
  TRobol()
{
  /// This constructor is for ROOT streamers ONLY
  cout<< "@@@@> TRobolDEV DEFAULT Constructor (for ROOT only) "<<endl;
}

///_____________________________________________________________
TRobolDEV::TRobolDEV(const char* Name):
  TRobol(Name)
{
//! Constructor to be used by the user!!!
//! Its important role is to define ALL DEFAULTS.
//! to changed by the user before calling TMCgen::Initialize
  cout<< "@@@@> TRobolDEV::TRobolDEV USER Constructor "<<endl;
  m_xmin    = 0.01;        //! x range
  m_xmax    = 0.99;        //! x range

}///

///______________________________________________________________________________________
TRobolDEV::~TRobolDEV()
{
  //!Explicit destructor
  cout<< "@@@@> TRobolDEV::TRobolDEV !!!! DESTRUCTOR !!!! "<<endl;
}///destructor

//______________________________________________________________________________
void TRobolDEV::Initialize(
        ofstream *OutFile, /// Central log-file for messages
        TFile *GenFile,    /// ROOT disk file for CRNG and MC gen.
        TFile *HstFile)    /// ROOT disk file for histograms
{
  cout<< "****> TRobolDEV::Initialize starts"<<endl;
  TRobol::Initialize(OutFile,GenFile,HstFile);
  ///
  //[[[[[[[[[[[[[[[[
    const char *output_file = "./kkmc.output";
    long stl2 = strlen(output_file);
    kk2f_fort_open_(16,output_file,stl2);
  //]]]]]]]]]]]]]]]]

  /// Book histograms or read them from the disk
  Hbooker();

  cout<< "****> TRobolDEV::Initialize: finished"<<endl;
}///Initialize

///______________________________________________________________________________
void TRobolDEV::Hbooker()
{
  ///
  cout<< "****> TRobolDEV::Hbooker: histogram booking STARTS"<<endl;
  BXOPE(*f_Out);
  BXTXT(*f_Out,"========================================");
  BXTXT(*f_Out,"======    TRobol::Hbooker    ===========");
  BXTXT(*f_Out,"========================================");
  BXCLO(*f_Out);
  f_HstFile->cd();

/////////////////////////////////////////////////////////////////////////////////////////
//  ************* user histograms  *************
  hst_weight1 = TH1D_UP("hst_weight1" ,  "MC weight",      100, -1.5, 2.0);


///////////////////////////////////////////////////////////////////////////////////////////
  BXOPE(*f_Out);
  BXTXT(*f_Out,"========================================");
  BXTXT(*f_Out,"====== END of TRobol::Hbooker ==========");
  BXTXT(*f_Out,"========================================");
  BXCLO(*f_Out);
  cout<< "****> TRobolDEV::Hbooker: histogram booking FINISHED"<<endl;
}///Hbooker

///______________________________________________________________________________
void TRobolDEV::Production(double &iEvent)
{
/////////////////////////////////////////////////////////////
  double wt5;
  TMCgenDEV *MCgen = (TMCgenDEV*)f_MCgen;

/////////////////////////////////////////////////////////////
/// MC generation in base class, ISR+FSR+IFI event
  TRobol::Production(iEvent);  // It invokes MCgen->Generate
  /// filling in histos
  MCgen->f_FoamI->GetMCwt(wt5);


///
}///Production

///______________________________________________________________________________
void TRobolDEV::Finalize()
{
/////////////////////////////////////////////////////////////
///------------------------

  Double_t MCresult, MCerror, MCnorm, Errel;
  TMCgenDEV *MCgen = (TMCgenDEV*)f_MCgen;
  cout << "**************************************************************"<<endl;
  cout << "**************** TRobolDEV::Finalize  ***********************"<<endl;

  MCgen->Finalize();

 cout << "**************************************************************"<<endl;

}
