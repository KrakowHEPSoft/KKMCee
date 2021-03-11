//                CLASS ROBOL                                                //

#include "TRobolKKMC.h"

# define sw2 setprecision(10) << setw(18)

ClassImp(TRobolKKMC);

TRobolKKMC::TRobolKKMC():
  TRobol()
{
  /// This constructor is for ROOT streamers ONLY
  cout<< "@@@@> TRobolKKMC DEFAULT Constructor (for ROOT only) "<<endl;
}


///_____________________________________________________________
TRobolKKMC::TRobolKKMC(const char* Name):
  TRobol(Name)
{
//! Constructor to be used by the user!!!
//! Its important role is to define ALL DEFAULTS.
//! to changed by the user before calling TMCgen::Initialize
  cout<< "@@@@> TRobolKKMC::TRobolFOAM USER Constructor "<<endl;
  m_NevGen=0;
  m_count1=0;

}//TRobolKKMC


///______________________________________________________________________________________
TRobolKKMC::~TRobolKKMC()
{
  //!Explicit destructor
  cout<< "@@@@> TRobolKKMC::TRobolFOAM !!!! DESTRUCTOR !!!! "<<endl;
}///destructor


//______________________________________________________________________________
//////////////////////////////////////////////////////////////
//   Initialize MC generator and analysis programs          //
//////////////////////////////////////////////////////////////
void TRobolKKMC::Initialize(
        ofstream *OutFile, /// Central log-file for messages
        TFile *GenFile,    /// ROOT disk file for CRNG and MC gen.
        TFile *HstFile)    /// ROOT disk file for histograms
{
  cout<< "****> TRobolKKMC::Initialize starts"<<endl;
  TRobol::Initialize(OutFile,GenFile,HstFile);
  ///
  /// Book histograms or read them from the disk
  Hbooker();

  cout<< "****> TRobolKKMC::Initialize: finished"<<endl;
}///Initialize

///////////////////////////////////////////////////////////////////////////////
void TRobolKKMC::Hbooker()
{
  ///
  cout<< "****> TRobolFOAM::Hbooker: histogram booking STARTS"<<endl;
  BXOPE(*f_Out);
  BXTXT(*f_Out,"========================================");
  BXTXT(*f_Out,"======    TRobol::Hbooker    ===========");
  BXTXT(*f_Out,"========================================");
  BXCLO(*f_Out);
  f_HstFile->cd();
  //  ************* user histograms  *************
  int nbin=1000;
  hst_weight  = TH1D_UP("hst_weight" ,  "MC weight",      100, 0.000 , 2.0);

}//Hbooker
///////////////////////////////////////////////////////////////////////////////
void TRobolKKMC::Production(double &iEvent)
{
/////////////////////////////////////////////////////////////////////////
//
//   GENERATE AND ANALYZE SINGLE EVENT
//
/////////////////////////////////////////////////////////////////////////
  // ****************************************************************
  // ************ Generate event and import it here  ****************
  m_NevGen++;
  KKee2f *KKMC_generator = (KKee2f*)f_MCgen;
  KKMC_generator->Generate(); // done in user class
 // IMPORT KKMC event and weights
  double WtMain,WtCrude;
 // KKMC_generator->GetWt(WtMain,WtCrude);
  WtMain = KKMC_generator->m_WtFoam; // temporary

  hst_weight->Fill(WtMain);              // histogramming

}//Production


///////////////////////////////////////////////////////////////////////////////
void TRobolKKMC::Finalize()
{
//   Finalize MC  run, final printouts, cleaning etc., xcheck of normalization
//   Plotting histograms is done independently using root file
  KKee2f *KKMC_generator = (KKee2f*)f_MCgen;
//
  double XsNormPb, XsErroPb;
  KKMC_generator->Finalize();
  /*
  XsNormPb =KKMC_generator->m_XsNormPb;
  XsErroPb =KKMC_generator->m_XsErroPb;
  cout << " KKMC: XsNormPb [pb] = "<<  XsNormPb << "  +-  "<< XsErroPb <<endl;
  double xSecPb,xErrPb,xSecNb;
  KKMC_generator->GetXsecMC( xSecPb, xErrPb);
  xSecNb=xSecPb/1000;
  cout << " KKMC: xSecPb   [pb] = "<<  xSecPb << "  +-  "<< xErrPb <<endl;
  cout << " KKMC: xSecNb   [nb] = "<<  xSecNb << "  +-  "<< xErrPb/1000 <<endl;
  */
}//Finalize



