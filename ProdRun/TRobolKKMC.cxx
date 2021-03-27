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
  double delv = 0.012;
  int nbin=100;
  hst_WtMain   = TH1D_UP("hst_WtMain" ,  "MC weight Main",   nbin, 0.00 , 2.0);
  hst_WtFoam   = TH1D_UP("hst_WtFoam" ,  "MC weight Foam",   nbin, 0.00 , 2.0);

  hst_vvBES    = TH1D_UP("hst_vvBES" ,   "BES distr",   nbin, -delv , delv);
  hst_vvTrue   = TH1D_UP("hst_vvTrue" ,  "vv distr",    nbin,  0.0 , 1.0);
  hst_nPhot    = TH1D_UP("hst_nPhot" ,   "nPhot",         20,  0.0 ,  20);
  hst_CosTheta = TH1D_UP("hst_CosTheta", "CosTheta",    nbin,  -1.0 , 1.0);

  int nbv =100;
  int nbc = 50;

  sca_vTcPR_Ceex2  = TH2D_UP("sca_vTcPR_Ceex2",   "dSig/dc/dv ", nbv, 0.0 ,1.0, nbc, -1.0 ,1.0);
  sca_vTcPR_Ceex2n = TH2D_UP("sca_vTcPR_Ceex2n",  "dSig/dc/dv ", nbv, 0.0 ,1.0, nbc, -1.0 ,1.0);
  sca_vTcPR_Eex0   = TH2D_UP("sca_vTcPR_Eex0",    "dSig/dc/dv ", nbv, 0.0 ,1.0, nbc, -1.0 ,1.0);
  sca_vTcPR_Eex2   = TH2D_UP("sca_vTcPR_Eex2",    "dSig/dc/dv ", nbv, 0.0 ,1.0, nbc, -1.0 ,1.0);

  sca_r1r2     = TH2D_UP("sca_r1r2" ,    "BES spectrum", 100, -delv , delv, 100, -delv , delv);

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
  KKevent *Event = KKMC_generator->m_Event;

  double WtMain,WtCrude,WtFoam;
 // KKMC_generator->GetWt(WtMain,WtCrude);
 // WtMain = KKMC_generator->m_WtFoam; // temporary
  WtMain = KKMC_generator->m_WtMain;
  WtFoam = KKMC_generator->m_WtFoam;

  TLorentzVector Pf1 = Event->m_Pf1;      //! initial beams
  TLorentzVector Pf2 = Event->m_Pf2;      //! initial beams
  TLorentzVector Qf1 = Event->m_Qf1;
  TLorentzVector Qf2 = Event->m_Qf2;

  double CosTheta = KKMC_generator->m_CosTheta;
  double r1 = Event->m_r1;
  double r2 = Event->m_r2;
//  double vv = Event->m_vv;
  int nPhot = Event->m_nPhot;

// ****************************************************************
// Kinematic variables to be monitored
// ****************************************************************
  double s  =(Pf1+Pf2)*(Pf1+Pf2);
  double s1 =(Qf1+Qf2)*(Qf1+Qf2);
  double CMSene = sqrt(s);
  double Mff    = sqrt(s1);
  double vv     = 1-s1/s;
  double CosThe1 = Qf1.CosTheta();
  double Theta1  = Qf1.Theta();
  double Phi1    = Qf1.Phi();
  double E1      = Qf1.Energy();

  double CosThe2 = Qf2.CosTheta();
  double Phi2    = Qf2.Phi();
  double Theta2  = Qf2.Theta();
  double E2      = Qf2.Energy();

  double Acol    = fabs(Phi1-(M_PI));

  double SinThe1,SinThe2,yy1,yy2,CosThePL,CosPRD,zAleph,s1Aleph,xe1,xe2;
//--------------------------------------------------------------------
//* Various definitions of Theta and s-propagator
//*--------------------------------------------------------------------
//** definition of P.L. B219, 103 (1989)
  CosThePL = ( E1*CosThe1 -E2*CosThe2)/(E1+E2);
//* definition of P.R. D41, 1425 (1990)
  SinThe1 = sqrt(fabs((1-CosThe1)*(1+CosThe1)));
  SinThe2 = sqrt(fabs((1-CosThe2)*(1+CosThe2)));
  yy1 = SinThe2/(SinThe1+SinThe2);
  yy2 = SinThe1/(SinThe1+SinThe2);
  CosPRD = yy1*CosThe1 - yy2*CosThe2;
  xe1 = 2*E1/CMSene;
  xe2 = 2*E2/CMSene;
//*-------------------------------
//* LL formula for s'/s from angles according to ALEPH note 1996
  zAleph =  (sin(Theta1)+sin(Theta2) -fabs(sin(Theta1+Theta2)))
           /(sin(Theta1)+sin(Theta2) +fabs(sin(Theta1+Theta2)));
  s1Aleph  = s*zAleph;
//*--------------------------------------------------------------------

//  if( m_NevGen<100 ) cout<<"*vv="<<vv<<endl;

// ****************************************************************
// HISTOGRAMMING
// ****************************************************************
  double WtEEX2  = KKMC_generator->m_WtAlter[73];
  double WtEEX0  = KKMC_generator->m_WtAlter[71];
  double WtCEEX1 = KKMC_generator->m_WtAlter[202];    //  CEEX Weight O(alf1)
  double WtCEEX1n= KKMC_generator->m_WtAlter[252];    //  CEEX Weight O(alf1)
  double WtCEEX2 = KKMC_generator->m_WtAlter[203];    //  CEEX Weight O(alf2)
  double WtCEEX2n= KKMC_generator->m_WtAlter[253];    //  CEEX Weight O(alf2) IFI off
  double WtCEEX0 = KKMC_generator->m_WtAlter[201];    //  CEEX Weight O(alf0)
  double WtCEEX0n= KKMC_generator->m_WtAlter[251];    //  CEEX Weight O(alf0) IFI off
//  WtEEX2=WtEEX0;    //!!!!!!!!!!!! DEBUG
//  WtEEX2=WtCEEX2n;  //!!!!!!!!!!!! DEBUG

  hst_WtMain->Fill(WtMain);
  hst_WtFoam->Fill(WtFoam);

  hst_nPhot->Fill(nPhot,WtMain);

  hst_vvTrue->Fill(vv,WtEEX2);

// big scatergrams, range vv< 1.0
  sca_vTcPR_Eex0->Fill(    vv, CosPRD, WtEEX0);
  sca_vTcPR_Eex2->Fill(    vv, CosPRD, WtEEX2);
  sca_vTcPR_Ceex2->Fill(   vv, CosPRD, WtCEEX2);
  sca_vTcPR_Ceex2n->Fill(  vv, CosPRD, WtCEEX2n); // IFI off

//----------------- BES corner ------------------------
  double vvBES = r1 + r2;
  hst_vvBES->Fill(vvBES,WtMain);

  hst_CosTheta->Fill(CosTheta,WtMain);

  sca_r1r2->Fill(r1, r2, WtMain);

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



