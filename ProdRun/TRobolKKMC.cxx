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
  //   for muon pairs
  double delv = 0.007;
  delv = 0.0025; // for beamstrahlung
  int nbin=200;
  hst_WtMain   = TH1D_UP("hst_WtMain" ,  "MC weight Main",   nbin, 0.00 , 8.0);
  hst_WtMain200  = TH1D_UP("hst_WtMain200" , "MC weight Main",10*nbin, 0.00 , 100.0);
  hst_WtFoam   = TH1D_UP("hst_WtFoam" ,  "MC weight Foam",   nbin, 0.00 , 8.0);
  hst_WtCeex2n = TH1D_UP("hst_WtCeex2n" ,"WTmain IFI off",   nbin, 0.00 , 8.0);

//
  hst_vvBES    = TH1D_UP("hst_vvBES" ,   "BES distr",   nbin, -delv , delv);
  hst_vvTrue   = TH1D_UP("hst_vvTrue" ,  "vv distr",    nbin,  0.0 , 1.0);
  hst_nPhot    = TH1D_UP("hst_nPhot" ,   "nPhot",         20,  0.0 ,  20);
  //
  hst_CosTheta = TH1D_UP("hst_CosTheta", "CosTheta",    nbin,  -1.0 , 1.0);
  hst_CosThOve = TH1D_UP("hst_CosThOve", "CosTheta",    nbin,  -1.0 , 1.0);
//
  int nbv =100;
  int nbc = 50;
//
  sca_vTcPR_Ceex2  = TH2D_UP("sca_vTcPR_Ceex2",   "dSig/dc/dv ", nbv, 0.0 ,1.0, nbc, -1.0 ,1.0);
  sca_vTcPR_Ceex2n = TH2D_UP("sca_vTcPR_Ceex2n",  "dSig/dc/dv ", nbv, 0.0 ,1.0, nbc, -1.0 ,1.0);
  sca_vTcPR_Eex0   = TH2D_UP("sca_vTcPR_Eex0",    "dSig/dc/dv ", nbv, 0.0 ,1.0, nbc, -1.0 ,1.0);
  sca_vTcPR_Eex2   = TH2D_UP("sca_vTcPR_Eex2",    "dSig/dc/dv ", nbv, 0.0 ,1.0, nbc, -1.0 ,1.0);

  sca_r1r2     = TH2D_UP("sca_r1r2" ,    "BES spectrum", 100, -delv , delv, 100, -delv , delv);
////////////////////////////////////////////////////////
// --------------------  for neutrinos ----------------
  KKee2f *KKMC_generator = (KKee2f*)f_MCgen;
  double CMSene = KKMC_generator->m_xpar[1];
  double MZ     = KKMC_generator->m_xpar[502];
  double vvZ    = 1-(MZ*MZ)/(CMSene*CMSene);
  //double rat =1.0;
  //rat = 0.715/0.679215;      // later version for Roy
  double vv2    = vvZ+0.020;
  double vv1    = vvZ-0.020;
  //vv2 = vvZ*rat;
  //vv1 = vvZ/rat;
//
  hst_LnThPhAll = TH1D_UP("hst_LnThPhAll", "ln10(sin(theta)) all phot.",  60, -6.0 ,0.0);
  hst_LnThPhVis = TH1D_UP("hst_LnThPhVis", "ln10(sin(theta)) vis. phot.", 60, -6.0 ,0.0);
//
  nbv =100;   // too small
  hst_vtNuCeex2 = TH1D_UP("hst_vtNuCeex2",  "dSig/dvTrue ",    nbv, 0.000 ,1.000);
  hst_vaNuCeex2 = TH1D_UP("hst_vaNuCeex2",  "dSig/dv WTmain ", nbv, 0.000 ,1.000);
//
  int nbv2 =40;   // small range?
  hst_vPhotNuel = TH1D_UP("hst_vPhotNuel", "dSig/dv CEEX1",       nbv2, vv1 , vv2);
  hst_vPhotNumu = TH1D_UP("hst_vPhotNumu", "dSig/dv CEEX2",       nbv2, vv1 , vv2);
  //  ================= all neutrino
  hst_vvNuCeex1 = TH1D_UP("hst_vvNuCeex1", "dSig/dv CEEX1",       nbv2, vv1 , vv2);
  hst_vvNuCeex2 = TH1D_UP("hst_vvNuCeex2", "dSig/dv CEEX2",       nbv2, vv1 , vv2);
  hst_vvNuCeex12= TH1D_UP("hst_vvNuCeex12","dSig/dv CEEX1-CEEX2", nbv2, vv1 , vv2);

  hst_nPhAll  = TH1D_UP("hst_nPhAll" , "No. of photons, all",   8, -0.5 ,7.5);
  hst_nPhVis  = TH1D_UP("hst_nPhVis" , "No. photons, E>10MeV",  8, -0.5 ,7.5);

  nbv =100;
  hst_vaNuElCeex2 = TH1D_UP("hst_vaNuElCeex2",  "dSig/dv WTmain ", nbv, 0.000 ,1.000);
  hst_vaNuMuCeex2 = TH1D_UP("hst_vaNuMuCeex2",  "dSig/dv WTmain ", nbv, 0.000 ,1.000);
  hst_vaNuTaCeex2 = TH1D_UP("hst_vaNuTaCeex2",  "dSig/dv WTmain ", nbv, 0.000 ,1.000);

  hst_vxNuElCeex2 = TH1D_UP("hst_vxNuElCeex2",  "dSig/dv WTmain ", nbv, 0.000 ,1.000);
  hst_vxNuMuCeex2 = TH1D_UP("hst_vxNuMuCeex2",  "dSig/dv WTmain ", nbv, 0.000 ,1.000);
  hst_vxNuTaCeex2 = TH1D_UP("hst_vxNuTaCeex2",  "dSig/dv WTmain ", nbv, 0.000 ,1.000);

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
  int KFfin = Event->m_KFfin;          // final fermion ident
  int NuYes = 0;
  if( KFfin==12 || KFfin==14 || KFfin==16 ) NuYes=1;

  double CosTheta = KKMC_generator->m_CosTheta; // dummy variable, not used
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
// Various definitions of Theta and s-propagator
//--------------------------------------------------------------------
// definition of P.L. B219, 103 (1989)
  CosThePL = ( E1*CosThe1 -E2*CosThe2)/(E1+E2);
// definition of P.R. D41, 1425 (1990)
  SinThe1 = sqrt(fabs((1-CosThe1)*(1+CosThe1)));
  SinThe2 = sqrt(fabs((1-CosThe2)*(1+CosThe2)));
  yy1 = SinThe2/(SinThe1+SinThe2);
  yy2 = SinThe1/(SinThe1+SinThe2);
  CosPRD = yy1*CosThe1 - yy2*CosThe2;
  xe1 = 2*E1/CMSene;
  xe2 = 2*E2/CMSene;
//-------------------------------
// LL formula for s'/s from angles according to ALEPH note 1996
  zAleph =  (sin(Theta1)+sin(Theta2) -fabs(sin(Theta1+Theta2)))
           /(sin(Theta1)+sin(Theta2) +fabs(sin(Theta1+Theta2)));
  s1Aleph  = s*zAleph;
//--------------------------------------------------------------------
//[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
//if( Event->m_KeyDBG == 1){
  if(m_NevGen<5){
  (*f_Out)<< ">>>RobolKKMC::Production: CosPRD= "<<CosPRD<<" CosThePL= "<<CosThePL<<"  vv="<<vv<<endl;
  cout    << ">>>RobolKKMC::Production: CosPRD= "<<CosPRD<<" CosThePL= "<<CosThePL<<"  vv="<<vv<<endl;
  }//NevGen
//}//m_KeyDBG
//]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]

// ****************************************************************
//          various weights
// ****************************************************************
  double WtEEX2  = KKMC_generator->m_WtAlter[73];
  double WtEEX0  = KKMC_generator->m_WtAlter[71];
  double WtCEEX1 = KKMC_generator->m_WtAlter[202];    //  CEEX Weight O(alf1)
  double WtCEEX1n= KKMC_generator->m_WtAlter[252];    //  CEEX Weight O(alf1)
  double WtCEEX2 = KKMC_generator->m_WtAlter[203];    //  CEEX Weight O(alf2)
  double WtCEEX2n= KKMC_generator->m_WtAlter[253];    //  CEEX Weight O(alf2) IFI off
  double WtCEEX0 = KKMC_generator->m_WtAlter[201];    //  CEEX Weight O(alf0)
  double WtCEEX0n= KKMC_generator->m_WtAlter[251];    //  CEEX Weight O(alf0) IFI off

// ****************************************************************
//          HISTOGRAMMING
// ****************************************************************
  hst_WtMain->Fill(WtMain);
  hst_WtMain200->Fill(WtMain);
//
  hst_WtFoam->Fill(WtFoam);
  hst_WtCeex2n->Fill(WtCEEX2n);
  hst_nPhot->Fill(nPhot,WtMain);
  hst_vvTrue->Fill(vv,WtEEX2);
//
  hst_CosTheta->Fill(CosPRD,WtMain);
  if(WtMain>1.0) hst_CosThOve->Fill(CosPRD,WtMain-1.0);
//======================================================
//                Muon pairs
  //======================================================
if( KFfin == 13 ) {
// big scatergrams, range vv< 1.0
  sca_vTcPR_Eex0->Fill(    vv, CosPRD, WtEEX0);
  sca_vTcPR_Eex2->Fill(    vv, CosPRD, WtEEX2);
  sca_vTcPR_Ceex2->Fill(   vv, CosPRD, WtCEEX2);
  sca_vTcPR_Ceex2n->Fill(  vv, CosPRD, WtCEEX2n); // IFI off
}// Muon pairs
//======================================================
//       all 3 nu+nubar + visible photons
//======================================================
if( NuYes==1 ){
/// photon acceptance data
double XEneMin = 0.10;  /// Emin/Ebeam  for visible photon
double XTraMin = 0.02;  /// kTmin/Ebeam for visible photon
double ThetaMin =  15;  /// theta minimum  for visible photon
/// photon acceptance params
int    nph_vis=0;       /// No. of visible (triggered) photons
double phEneVis = 0;    /// Energy of triggered photon
int    TrigPho=1;       /// =1 for accepted, =0 for rejected
//----------------------------------------------------
// Loop over photons. defining visible photons
double phEne,phCosth,phTheta,phPT;
for(int iphot=0;iphot<nPhot;iphot++){
  phEne   = Event->m_PhotAll[iphot+1].Energy();
  phCosth = Event->m_PhotAll[iphot+1].CosTheta();
  phTheta = Event->m_PhotAll[iphot+1].Theta()*180/M_PI;
  phPT    = Event->m_PhotAll[iphot+1].Pt();
  if( phEne < XEneMin*CMSene/2 ) TrigPho=0;
  if( phPT  < XTraMin*CMSene/2 ) TrigPho=0;
  if( phTheta < ThetaMin       ) TrigPho=0;
  if( phTheta > (180-ThetaMin) ) TrigPho=0;
  if( TrigPho)  nph_vis++;
  if( TrigPho)  phEneVis +=phEne;
  /// histogramming photon angle, inclusive, neutrino channels only
  hst_LnThPhAll->Fill(               log10(phPT/phEne), WtMain); // All
  if( TrigPho)  hst_LnThPhVis->Fill( log10(phPT/phEne), WtMain); // triggered
}// for iphot
//--------------------------
/// true nu-pair mass, full range
hst_vtNuCeex2->Fill( vv, WtCEEX2);
/// photon multiplicity, triggered and not
hst_nPhAll->Fill(  nPhot,WtMain);
if(  nph_vis > 0 )hst_nPhVis->Fill(  nph_vis,WtMain);
//-----------------
/// Energy distribution for triggered photon
/// Warning: vPhot not necessarily for the hardest triggered photon!!!
double vPhot = 2*phEneVis/CMSene;
//-----------------------------------------
if( KFfin==12 ) hst_vxNuElCeex2->Fill( vv, WtCEEX2); // wide v-range
if( KFfin==14 ) hst_vxNuMuCeex2->Fill( vv, WtCEEX2); // wide v-range
if( KFfin==16 ) hst_vxNuTaCeex2->Fill( vv, WtCEEX2); // wide v-range
//-----------------------------------------
// photon tagging
if( nph_vis >= 1){
  hst_vaNuCeex2->Fill( vPhot, WtCEEX2); /// tagged, full v-range
/// nu_el and nu_mu separately
  if( KFfin==12 ) hst_vaNuElCeex2->Fill( vPhot, WtCEEX2); // wide v-range
  if( KFfin==14 ) hst_vaNuMuCeex2->Fill( vPhot, WtCEEX2); // wide v-range
  if( KFfin==16 ) hst_vaNuTaCeex2->Fill( vPhot, WtCEEX2); // wide v-range
/// comparing Nuel with Numu
  if( KFfin==12 ) hst_vPhotNuel->Fill( vPhot, WtCEEX2); // v-narrow
  if( KFfin==14 ) hst_vPhotNumu->Fill( vPhot, WtCEEX2); // v-narrow
/// histos with restricted v-range
  hst_vvNuCeex1->Fill( vPhot, WtCEEX1);
  hst_vvNuCeex2->Fill( vPhot, WtCEEX2);
  hst_vvNuCeex12->Fill(vPhot, WtCEEX1-WtCEEX2);
/// absolute energy
//  hst_evNuCeex1->Fill( phEneVis, WtCEEX1);
//  hst_evNuCeex2->Fill( phEneVis, WtCEEX2);
//  hst_evNuCeex12->Fill(phEneVis, WtCEEX1-WtCEEX2);
  }// nph_vis
}// NuYes
//======================================================

//-----------------------------------------------------
//----------------- BES corner ------------------------
  double vvBES = r1 + r2;
  hst_vvBES->Fill(vvBES,WtMain);
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



