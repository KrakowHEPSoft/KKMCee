///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                Class TRobol                                              //
//                                                                           //
//    This class generates MC events and collects material for analysis      //
//    Not that this class is not really ment for persistency, so far...      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TRobol.h"

#include"BXFORMAT.h"
#define SP12 setprecision(7) << setw(12)
#define SP18 setprecision(9) << setw(18)
#define SW5  setw(5)
#define SW9  setw(9)
#define SW12 setw(12)
#define SW18 setw(18)
#define SP21 setw(21)<<setprecision(13)
#define SP15 setw(15)<<setprecision(9)
#define SP10 setw(10)<<setprecision(5)
 
ClassImp(TRobol);

//______________________________________________________________________________
TRobol::TRobol()
{
/// explicit default constructor for ROOT streamers
/// should not be used by the user
cout<< "****> TRobol DEFAULT constructor (only for ROOT)"<<endl;
  f_RNgen    = NULL;
  f_MCgen    = NULL;
  f_HstFile  = NULL;
  f_GenFile  = NULL;
  f_Out      = NULL;
}/// DEFAULT constructor

///______________________________________________________________________________________
TRobol::TRobol(const char* Name){
//! Constructor to be used by the user!!!
  cout<< "****> TRobol::TRobol USER Constructor "<<endl;
  ///
  if(strlen(Name)  >65) StopM("TRobol::TRobol: +++ stop, Name too long");
  sprintf(f_Name,"%s",Name);         // Class name
  ///
  f_HistNormName = "h_TMCgen_NORMA";
  f_RNgen    = NULL;
  f_MCgen    = NULL;
  f_HstFile  = NULL;
  f_GenFile  = NULL;
  f_Out      = NULL;

}/// USER constructor

///______________________________________________________________________________
TRobol::~TRobol(){
/// Explicit destructor
  cout<< "****> TRobol:: default DESTRUCTOR !!!! "<<endl;
}


//______________________________________________________________________________
TH1D *TRobol::TH1D_UP(const char* name, const char* title, 
                       int nbins, double xmin, double xmax)
{
  TH1D *h;
  if (f_isNewRun) {
    h = new TH1D(name,title,nbins,xmin,xmax);
    h->Sumw2();
  } else {
    h = (TH1D*)f_HstFile->Get(name);
  }
  return h;
}

//______________________________________________________________________________
TH2D *TRobol::TH2D_UP(const char* name, const char* title, 
                       int nbinx, double xmin, double xmax,
                       int nbiny, double ymin, double ymax)
{
  TH2D *h;
  if (f_isNewRun) {
    h = new TH2D(name,title,nbinx,xmin,xmax,nbiny,ymin,ymax);
    h->Sumw2();
  } else {
    h = (TH2D*)f_HstFile->Get(name);
  }
  return h;
}


//______________________________________________________________________________
void TRobol::Initialize(
        ofstream *OutFile, /// Central log-file for messages
        TFile *GenFile,    /// ROOT disk file for CRNG and MC gen.
        TFile *HstFile)    /// ROOT disk file for histograms
{
  //////////////////////////////////////////////////////////////
  //   Initialize MC generator and analysis programs          //
  //////////////////////////////////////////////////////////////
  //
  f_Out        = OutFile;
  f_GenFile    = GenFile;
  f_HstFile    = HstFile;
  //
  f_NevGen=0;
  f_count1=0;

  BXOPE(*f_Out);
  BXTXT(*f_Out,"========================================");
  BXTXT(*f_Out,"======    TRobol::Initialize     ======");
  BXTXT(*f_Out,"========================================");
  BXCLO(*f_Out);
  ///
  /// book histogram keeping track of overall normalization
  f_HstFile->cd(); /// to keep h_TMCgen_NORMA in  f_HstFile!!!
  TH1D *h_TMCgen_NORMA = TH1D_UP(f_HistNormName,"Normalization histo",100,0,100);
  //TH1D *h_TMCgen_NORMA = new TH1D("h_TMCgen_NORMA","Normalization histo",100,0,100);
  ///
  /// Read prepared object MC of generator from the disk
  f_GenFile->cd();
  f_MCgen = (TMCgen*)f_GenFile->Get("MCgen");
  f_RNgen =(TRandom*)f_GenFile->Get("RN_gen");  // read r.n. generator
  ///
  /// This is virgin run if MC generator NOT initialized
  f_isNewRun = f_MCgen->GetIsNewRun();
  ///
  /// Initialize/Redress MC generator
  if( f_isNewRun == 1){
    f_MCgen->f_GenFile = f_GenFile; /// Adding access to disk files,
    f_MCgen->f_HstFile = f_HstFile; /// just in case it is needed.
    f_MCgen->Initialize(f_RNgen,f_Out,h_TMCgen_NORMA);
  } else {
    f_MCgen->f_GenFile = f_GenFile;
    f_MCgen->Redress(   f_RNgen,f_Out,h_TMCgen_NORMA);
    f_MCgen->f_HstFile = f_HstFile;
  }
  /// and write MC generator into disk file.
  f_MCgen->Write("MCgen",TObject::kOverwrite);   // MC generator status
  ///
  BX1I(*f_Out," isNewRun", f_isNewRun," is it new MC run?       =");
  BXTXT(*f_Out,"========================================");
  BXTXT(*f_Out,"====== END of TRobol::Initialize ======");
  BXTXT(*f_Out,"========================================");
  BXCLO(*f_Out);
  cout<< "####> TRobol::Initialize FINISHED"<<endl;

}///Initialize


///______________________________________________________________________________
void TRobol::Production(double &iEvent)
{
  //////////////////////////////////////////////////////////////
  ///      GENERATE SINGLE EVENT
  f_NevGen++;
  f_MCgen->Generate();
  ///
}///Production

//______________________________________________________________________________
void TRobol::FileDump()
{
  f_HstFile->cd();
  f_HstFile->Write("",TObject::kOverwrite);       // All histograms
  f_HstFile->Save();                              // may be helps??
  f_HstFile->Flush();                             // may be helps??
  f_GenFile->cd();
  f_RNgen->Write("RN_gen",TObject::kOverwrite);  // status of r.n. generator
  f_GenFile->Flush();
}///FileDump

//______________________________________________________________________________
void TRobol::Finalize()
{
  //////////////////////////////////////////////////////////////
  //   Finalize MC  run, final printouts, cleaning etc.,      //
  //   Plotting histograms done independently using root file //
  //////////////////////////////////////////////////////////////
  //f_TraceFile.close();
  //
  f_MCgen->Finalize();
  //
  BXOPE(*f_Out);
  BXTXT(*f_Out,"========================================");
  BXTXT(*f_Out,"======     TRobol::Finalize      ======");
  BXTXT(*f_Out,"========================================");
  BXCLO(*f_Out);
  //============================================================================
}///Finalize
