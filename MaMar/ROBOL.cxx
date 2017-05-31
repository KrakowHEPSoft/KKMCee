///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                Class ROBOL                                                //
//                                                                           //
//    It contains Makers for MC generation and Analysis event per event      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
#include "ROBOL.h"

# define sw2 setprecision(10) << setw(18)

//ClassImp(ROBOL)

///////////////////////////////////////////////////////////////////////////////
//      *************** temporary entries from KKMC ****************
//      SUBROUTINE KarLud_GetVVxx(vv,x1,x2)
extern "C" void  karlud_getvvxx_(double&, double&, double&);
//extern "C" void  pyhepc_(int&);
//extern "C" void  photos_(int&);
//extern "C" void  phoini_();
//extern "C" void  hepevt_setphotosflagtrue_(int&);
//extern "C" void  hepevt_getnhep_(int&);
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
void ROBOL::Initialize(long &NevTot)
{
  //////////////////////////////////////////////////////////////
  //   Initialize MC generator and analysis programs          //
  //////////////////////////////////////////////////////////////
  m_NevGen=0;
  m_count1=0;
  const int jmax =10000;
  ReaData("../../.KK2f_defaults", jmax, m_xpar);  // numbering as in input!!!
  ReaData("./pro.input",         -jmax, m_xpar);  // jmax<0 means no-zeroing
  double ypar[jmax];
  for(int j=0;j<jmax;j++) ypar[j]=m_xpar[j+1];    // ypar has c++ numbering
  //
  NevTot = (long)m_xpar[0];                       // NevTot hidden in xpar[0] !!!
  KKMC_generator = new KKMC();
  KKMC_generator->Initialize(ypar);
  cout<<"ROBOL::Initialize:  NevTot = "<<NevTot<<endl;
  //  ************* user histograms  *************
  double CMSene = m_xpar[1];
  int nbin =NevTot/100;
  if(nbin>1000) nbin=1000;
  hst_weight  = new TH1D("hst_weight" ,  "MC weight",      100, 0.000 , 2.0);
  hst_Mff     = new TH1D("hst_Mff"    ,  "Mass(f-fbar)",  nbin, 0.000 ,CMSene);  
  hst_weight->Sumw2();
  hst_Mff->Sumw2();
  hst_nPhAll  = new TH1D("hst_nPhAll" , "No. of photons, all",   8, -0.5 ,7.5);
  hst_nPhVis  = new TH1D("hst_nPhVis" , "No. photons, E>10MeV",  8, -0.5 ,7.5);
  hst_nPhAll->Sumw2();
  hst_nPhVis->Sumw2();
  int nbv = 50;
  hst_vTrueMain    = new TH1D("hst_vTrueMain",   "dSig/dvTrue ", nbv, 0.000 ,1.000);
  hst_vTrueCeex2   = new TH1D("hst_vTrueCeex2",  "dSig/dvTrue ", nbv, 0.000 ,1.000);
  hst_vXGenCeex2   = new TH1D("hst_vXGenCeex2",  "dSig/dvTrue ", nbv, 0.000 ,1.000);
  hst_vTrueMain->Sumw2();
  hst_vTrueCeex2->Sumw2();
  hst_vXGenCeex2->Sumw2();
  hst_vAlepCeex2   = new TH1D("hst_vAlepCeex2",  "dSig/dvTrue ", nbv, 0.000 ,1.000);
  hst_vAlepCeex2->Sumw2();
  int nbc =50;
  hst_Cost1Ceex2= new TH1D("hst_Cost1Ceex2",  "dSig/cThet1   ", nbc, -1.000 ,1.000);
  hst_CosPLCeex2= new TH1D("hst_CosPLCeex2", "dSig/cThetPL  ", nbc, -1.000 ,1.000);
  hst_CosPRCeex2= new TH1D("hst_CosPRCeex2", "dSig/cThetPRD ", nbc, -1.000 ,1.000);
  hst_Cost1Ceex2->Sumw2();
  hst_CosPLCeex2->Sumw2();
  hst_CosPRCeex2->Sumw2();
  hst_CosPREex2= new TH1D("hst_CosPREex2", "dSig/cThetPRD ", nbc, -1.000 ,1.000);
  hst_CosPREex2->Sumw2();
  // scatergrams
  sca_vTcPR_Ceex2= new TH2D("sca_vTcPR_Ceex2",  "dSig/dc/dv ", nbv, 0.0 ,1.0, nbc, -1.0 ,1.0);
  sca_vTcPR_Eex2 = new TH2D("sca_vTcPR_Eex2",   "dSig/dc/dv ", nbv, 0.0 ,1.0, nbc, -1.0 ,1.0);
  sca_vTcPR_Ceex2->Sumw2();
  sca_vTcPR_Eex2->Sumw2();
  sca_vTcPR_Ceex2n= new TH2D("sca_vTcPR_Ceex2n",  "dSig/dc/dv ", nbv, 0.0 ,1.0, nbc, -1.0 ,1.0);
  sca_vTcPR_Ceex2n->Sumw2();
  sca_vXcPR_Ceex2= new TH2D("sca_vXcPR_Ceex2",  "dSig/dc/dv ", nbv, 0.0 ,1.0, nbc, -1.0 ,1.0);
  sca_vXcPR_Eex2 = new TH2D("sca_vXcPR_Eex2",   "dSig/dc/dv ", nbv, 0.0 ,1.0, nbc, -1.0 ,1.0);
  sca_vXcPR_Ceex2->Sumw2();
  sca_vXcPR_Eex2->Sumw2();
  //  New bigger scatergrams, restricted vmax
  int NBv =100; int NBc = 100;
  double vmx2= 0.20;
  sct_vTcPR_Ceex2 = new TH2D("sct_vTcPR_Ceex2",  "dSig/dc/dv ", NBv, 0.0 ,vmx2, NBc, -1.0 ,1.0);
  sct_vTcPR_Ceex2n= new TH2D("sct_vTcPR_Ceex2n", "dSig/dc/dv ", NBv, 0.0 ,vmx2, NBc, -1.0 ,1.0);
  sct_vTcPR_Ceex2->Sumw2();
  sct_vTcPR_Ceex2n->Sumw2();
  //
  sct_vAcPR_Ceex2= new TH2D("sct_vAcPR_Ceex2",  "dSig/dc/dv ", NBv, 0.0 ,vmx2, NBc, -1.0 ,1.0);
  sct_vAcPR_Ceex2->Sumw2();
  sct_vAcPR_Ceex2n= new TH2D("sct_vAcPR_Ceex2n","dSig/dc/dv ", NBv, 0.0 ,vmx2, NBc, -1.0 ,1.0);
  sct_vAcPR_Ceex2n->Sumw2();
  sct_vTcPL_Ceex2= new TH2D("sct_vTcPL_Ceex2",  "dSig/dc/dv ", NBv, 0.0 ,vmx2, NBc, -1.0 ,1.0);
  sct_vTcPL_Ceex2->Sumw2();
  sct_vKcPL_Ceex2= new TH2D("sct_vKcPL_Ceex2",  "dSig/dc/dv ", NBv, 0.0 ,vmx2, NBc, -1.0 ,1.0);
  sct_vKcPL_Ceex2->Sumw2();
  sct_vAcPL_Ceex2= new TH2D("sct_vAcPL_Ceex2",  "dSig/dc/dv ", NBv, 0.0 ,vmx2, NBc, -1.0 ,1.0);
  sct_vAcPL_Ceex2->Sumw2();
  // for xcheck and h.o ISR
  hst_vACeex2   = new TH1D("hst_vACeex2",  "dSig/dvTrue ", NBv, 0.000 ,vmx2);
  hst_vACeex21F = new TH1D("hst_vACeex21F","dSig/dvTrue ", NBv, 0.000 ,vmx2);
  hst_vACeex21B = new TH1D("hst_vACeex21B","dSig/dvTrue ", NBv, 0.000 ,vmx2);
  hst_vACeex2->Sumw2();
  hst_vACeex21F->Sumw2();
  hst_vACeex21B->Sumw2();

  //  ************* special histo  *************
  HST_KKMC_NORMA = new TH1D("HST_KKMC_NORMA","KKMC normalization &xpar",jmax,0.0,10000.0);
  for(int j=1; j<=jmax; j++)
    HST_KKMC_NORMA->SetBinContent(j,m_xpar[j]);    // xpar encoded
  m_YSum= 0.0;
  m_YSum2=0.0;
}
///////////////////////////////////////////////////////////////////////////////
void ROBOL::KKMC_NORMA()
{
  // Transfer normalization Record of KKMC to local histogram.
  // For later use in re-normalizing histostograms
  //
  double XsPrim; int NevPrim;
  KKMC_generator->GetPrimaNorma(XsPrim, NevPrim);
  HST_KKMC_NORMA->SetBinContent(0,XsPrim*NevPrim);
  HST_KKMC_NORMA->SetEntries(NevPrim);
  cout<<" ROBOL::KKMC_NORMA: XsPrim, NevPrim ="<< XsPrim <<"  "<<NevPrim << endl;
}
///////////////////////////////////////////////////////////////////////////////
void ROBOL::Production(long &iEvent)
{
/////////////////////////////////////////////////////////////////////////
//
//   GENERATE AND ANALYZE SINGLE EVENT
//
//
/////////////////////////////////////////////////////////////////////////
  // ****************************************************************
  // ************ Generate event and import it here  ****************
  m_NevGen++;
  KKMC_generator->Make();
  //
  PartImport(); // Import Pythia common block content into local LuPart matrix
  // IMPORT KKMC event and weights
  double WtMain,WtCrude;
  KKMC_generator->GetWt(WtMain,WtCrude);
  KKMC_generator->GetBeams(m_pbea1,m_pbea2);
  KKMC_generator->GetFermions(m_pfer1,m_pfer2);
  KKMC_generator->GetNphot(m_Nphot);                  // photon multiplicity
  TLorentzVector VSumPhot;    // By default all components are initialized by zero. 
  int iphot,iphot1;
  for(iphot=0;iphot<m_Nphot;iphot++){
    KKMC_generator->GetPhoton1(iphot+1,m_phot[iphot]);  // photon 4-momenta
    VSumPhot+= m_phot[iphot];
  }
  if(iEvent<10){
    cout<<"-----------------------------------------------------------  "<<iEvent;
    cout<<"  -----------------------------------------------------------"<<endl;
    cout<<" VSumPhot= "; MomPrint( VSumPhot );
    //KKMC_generator->Print1();
    KKMC_generator->PyList(2);
    PyPrint(1);
  }
  // ****************************************************************
  double s  =(m_pbea1+m_pbea2)*(m_pbea1+m_pbea2);
  double s1 =(m_pfer1+m_pfer2)*(m_pfer1+m_pfer2);
  double CMSene = sqrt(s);
  double Mff    = sqrt(s1);
  double vv     = 1-s1/s;
  // ********************************************************************
  // ***   Photon trigger TrigPho is for everybory, all pions, muons etc
  double Pi=4*atan(1.0);
  double phEne,phTheta,phCosth;
  double XEnePho  = 0.010;              // Emin for visible photom
  //**************************************************
  // Loop over photons, just in case
  //**************************************************
  int nph_ene=0;
  for(iphot=0;iphot<m_Nphot;iphot++){
    phEne   = m_phot[iphot].Energy();
    phCosth = m_phot[iphot].CosTheta();
    phTheta = m_phot[iphot].Theta()*180/Pi;
    if(phEne>XEnePho){
      nph_ene++;
    }
  }
  //********************************************************************
  // Muon trigger, it is not realy necessary if MC ir run for mu only
  //********************************************************************
  int TrigMu  = 0;
  // find muons, excluding muons from phi decays!!!
  int jMu1 =PartFindStable( 13);    // fortran numbering!!!
  int jMu2 =PartFindStable(-13);    // fortran numbering!!!
  m_pMu1  = m_Event[jMu1-1].fMom;    // fortran numbering!!!
  m_pMu2  = m_Event[jMu2-1].fMom;    // fortran numbering!!!
  int par1=m_Event[jMu1-1].fParent; // fortran numbering!!!
  int par2=m_Event[jMu2-1].fParent; // fortran numbering!!!
  if( (jMu1*jMu1)  && (par1 == par2) && (par1 == 3) ) TrigMu  = 1; // exclude backgr.
  //**************************************************************
  if( TrigMu && (m_count1<17) ){
    m_count1++;
    cout<<"**************************>>> two muons <<<****************************"<<endl;
    KKMC_generator->PyList(2);      
  }
  // muons,  vv, Q^2 costheta, etc
  double CosThe1 = m_pMu1.CosTheta();
  double Theta1  = m_pMu1.Theta();
  double E1      = m_pMu1.Energy();
  double CosThe2 = m_pMu2.CosTheta();
  double Theta2  = m_pMu2.Theta();
  double E2      = m_pMu2.Energy();
  double SinThe1,SinThe2,yy1,yy2,CosThePL,CosPRD,zAleph,s1Aleph;
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
//*-------------------------------
//* LL formula for s'/s from angles according to ALEPH note 1996
  zAleph =  (sin(Theta1)+sin(Theta2) -fabs(sin(Theta1+Theta2)))
           /(sin(Theta1)+sin(Theta2) +fabs(sin(Theta1+Theta2)));
  s1Aleph  = s*zAleph;
  //
  double vvK,x1,x2;
  karlud_getvvxx_(vvK,x1,x2);
  // *********************************************************************
  //          Most of histogramming starts here
  // *********************************************************************
  double WtEEX2  = KKMC_generator->GetWtAlter( 73);    //  Second ord. EEX2 O(alf2)
  double WtEEX3  = KKMC_generator->GetWtAlter( 74);    //  Third order EEX3 O(alf3)
  //cout<< "&&&&&& WtEEX2,3= "<<WtEEX2<<"  "<<WtEEX3<<endl;
  double WtCEEX1 = KKMC_generator->GetWtAlter(202);    //  CEEX Weight O(alf1)
  double WtCEEX1n= KKMC_generator->GetWtAlter(252);    //  CEEX Weight O(alf1)
  double WtCEEX2 = KKMC_generator->GetWtAlter(203);    //  CEEX Weight O(alf2)
  double WtCEEX2n= KKMC_generator->GetWtAlter(253);    //  CEEX Weight O(alf2) IFI off
  //
  double vvA = 1-zAleph;
  hst_nPhAll->Fill(  m_Nphot,WtMain);
  hst_nPhVis->Fill(  nph_ene,WtMain);
  hst_vTrueMain->Fill(   vv, WtMain);
  hst_vTrueCeex2->Fill(  vv, WtCEEX2);          // M(2f) of mun pair
  hst_vXGenCeex2->Fill(  vvK, WtCEEX2);         // M^star from MC (illegal)
  // semirealistic
  hst_vAlepCeex2->Fill(  vvA, WtCEEX2);         // M^star guessed
  //****[[[
  hst_vACeex2->Fill(       vvA, WtCEEX2);         // M^star guessed
  //hst_vACeex2->Fill(       vv, WtEEX3);         // M^star guessed
  if( CosPRD>0)
      hst_vACeex21F->Fill( vvA, WtCEEX2n-WtCEEX1); // M^star guessed
      //hst_vACeex21F->Fill( vv, WtEEX3-WtEEX2); // M^star guessed
  else
      hst_vACeex21B->Fill( vvA, WtCEEX2n-WtCEEX1); // M^star guessed
      //hst_vACeex21B->Fill( vv, WtEEX3-WtEEX2); // M^star guessed
  //****]]]
  if(vv<0.9){
    hst_Cost1Ceex2->Fill( CosThe1, WtCEEX2);
    hst_CosPLCeex2->Fill(CosThePL, WtCEEX2);
    hst_CosPRCeex2->Fill(  CosPRD, WtCEEX2);
  }
  hst_CosPREex2->Fill(  CosPRD, WtEEX2);
  // big scatergrams
  sca_vTcPR_Ceex2->Fill(   vv, CosPRD, WtCEEX2); 
  sca_vTcPR_Eex2->Fill(    vv, CosPRD, WtEEX2);
  sca_vXcPR_Ceex2->Fill(  vvK, CosPRD, WtCEEX2);
  sca_vXcPR_Eex2->Fill(   vvK, CosPRD, WtEEX2);
  // IFI off
  sca_vTcPR_Ceex2n->Fill(  vv, CosPRD, WtCEEX2n); // true v, IFI off
  //-------------------------------------------
  // New very BIG scaterplots
  sct_vTcPR_Ceex2->Fill(   vv, CosPRD, WtCEEX2);  // true v, IFI on
  sct_vTcPR_Ceex2n->Fill(  vv, CosPRD, WtCEEX2n); // true v, IFI off
  sct_vAcPR_Ceex2->Fill(   vvA, CosPRD,   WtCEEX2);  // Main CEEX2 KKMC , ISR+FSR
  sct_vAcPR_Ceex2n->Fill(  vvA, CosPRD,   WtCEEX2n); // IFI  off
  sct_vKcPL_Ceex2->Fill(   vvK, CosThePL, WtCEEX2);  // vv of Karlud (unphysical, pure ISR) thetaPL
  sct_vTcPL_Ceex2->Fill(    vv, CosThePL, WtCEEX2);  // vv bare muons
  sct_vAcPL_Ceex2->Fill(   vvA, CosThePL, WtCEEX2);  // Main CEEX2 KKMC , ISR+FSR

  // Miscelaneous
  m_YSum  += WtMain;
  m_YSum2 += WtMain*WtMain;
  hst_weight->Fill(WtMain);              // histogramming
  hst_Mff->Fill(Mff,WtMain);             // histogramming
  // debug debug debug debug debug debug debug
  if(iEvent<15){
    cout<<"-----------------------------------------------------------  "<<iEvent;
    cout<<"  -----------------------------------------------------------"<<endl;
    cout<< "vv, 1-zAleph    = "<<     vv<<"  "<<1-zAleph<<endl;
    cout<< "CosThe1,CosThe2 = "<<CosThe1<<"  "<<CosThe2<<endl;
    cout<< "CosThePL,CosPRD = "<<CosThePL<<"  "<<CosPRD<<endl;
  }
} //


///////////////////////////////////////////////////////////////////////////////
void ROBOL::Finalize()
{
//   Finalize MC  run, final printouts, cleaning etc., xcheck of normalization
//   Plotting histograms is done independently using root file
  double XsNormPb, XsErroPb;
  KKMC_generator->Finalize(XsNormPb, XsErroPb);
  cout << " KKMC: XsNormPb [pb] = "<<  XsNormPb << "  +-  "<< XsErroPb <<endl;
  double xSecPb,xErrPb,xSecNb;
  KKMC_generator->GetXsecMC( xSecPb, xErrPb);
  xSecNb=xSecPb/1000;
  cout << " KKMC: xSecPb   [pb] = "<<  xSecPb << "  +-  "<< xErrPb <<endl;
  cout << " KKMC: xSecNb   [nb] = "<<  xSecNb << "  +-  "<< xErrPb/1000 <<endl;
  // *********************************************************************
  // **** examples of normalizing histogram (will not work on farm)  *****
  int      nbt  = hst_Mff->GetNbinsX();
  Double_t tmax = hst_Mff->GetXaxis()->GetXmax();
  Double_t tmin = hst_Mff->GetXaxis()->GetXmin();
  Double_t Fact = nbt*XsNormPb/1000/(tmax-tmin)/m_NevGen; // now [nb]
  // **** re-normalized histogram as a clone
  TH1D *hstC_Mff =(TH1D*)hst_Mff->Clone();
  hstC_Mff->Sumw2();               // is it necessary???
  hstC_Mff->SetName("hstC_Mff"); // otherwise you get 2 histograms with same name
  hstC_Mff->SetTitle("dsigma/dQ2 [nb/GeV^2] ");
  hstC_Mff->Scale(Fact);
  // **** alternatively, re-normalized histo can be defined as a new one ****
  TH1D *hstN_Mff =new TH1D("hstN_Mff","dsigma/dQ2 [nb/GeV^2]",nbt,tmin,tmax);
  hstN_Mff->Sumw2();
  hstN_Mff->Add(hstN_Mff, Fact);
  // *********************************************************************
  // **** alternatively HST_KKMC_NORMA is used at later stage (plotting)
  int   NevPrim = HST_KKMC_NORMA->GetEntries();
  double XsPrima = HST_KKMC_NORMA->GetBinContent(0)/NevPrim;
  cout << "HST_KKMC_NORMA: XsPrima [nb] = "<< XsPrima << " NevPrim= "<< NevPrim <<endl;
  // *********************************************************************
  // Integrated Xsection canculated on-line
  cout << "///////////////////////////////////////////////////////////////"<<endl;
  double  XsecY = m_YSum/m_NevGen;                                // average weight
  double dXsecY = sqrt((m_YSum2/m_NevGen-XsecY*XsecY)/m_NevGen) ; // dispersion of wt.
  XsecY  *= XsNormPb/1000.0;  //  nanob.
  dXsecY *= XsNormPb/1000.0;  //  nanob.
  cout << " XsecY [nb] = "<< XsecY <<" +- "<< dXsecY <<endl;
  cout << "///////////////////////////////////////////////////////////////"<<endl;
}
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                             UTILITIES                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void ROBOL::PartImport(){
/////////////////////////////////////////////////////////////////////////
// Import Pythia common block content into local LuPart matrix
/////////////////////////////////////////////////////////////////////////
  m_Npart= KKMC_generator->GetPyNpart();
  if( (m_Npart<0) || (m_Npart>=4000) ){
    cout<<"++++ ROBOL::Production: STOP m_Npart= "<<m_Npart<<endl;
    exit(5);
  }
  for(int j=0; j<m_Npart;j++){
    KKMC_generator->GetPyParticle( j, m_Event[j]);  // import one particle
    //m_Event[j].Print(1);
  }
}
///////////////////////////////////////////////////////////////////////////////
int ROBOL::PartCount(const int flav){
  int jCount=0;
  for(int j=0; j<m_Npart;j++)
    if(m_Event[j].fFlafor == flav){
      jCount++;
    }
  return jCount;
}
///////////////////////////////////////////////////////////////////////////////
int ROBOL::PartFindAny(const int flav){
// fortran numbering!!!
  int jPosition=0;
  for(int j=0; j<m_Npart;j++)
    if(m_Event[j].fFlafor == flav){
      jPosition=j+1; break;
    }
  return jPosition;
}
///////////////////////////////////////////////////////////////////////////////
int ROBOL::PartFindStable(const int flav){
// fortran numbering!!!
  int jPosition=0;
  for(int j=0; j<m_Npart;j++)
    if((m_Event[j].fStatus  == 1)&&(m_Event[j].fFlafor == flav)){
      jPosition=j+1; break;
    }
  return jPosition;
}
////////////////////////////////////////////////////////////////////////////////
void ROBOL::PyPrint(const int mode){
//
// PRINT entire Pythia EVENT
//
  TLorentzVector Sum;
  cout<<"  "<<endl;
  cout<<"lser status flavor parent child1 child2";
  cout<<"                Px                Py                Pz               Ene";
  cout<<"              Mass"<<endl;
  for(int j=0; j<m_Npart;j++){
    m_Event[j].Print(mode);
    if(m_Event[j].fStatus == 1) Sum += m_Event[j].fMom;
  }
  cout<<"                  Total 4-momentum --> ";
  MomPrint(Sum);
}
void ROBOL::MomPrint( TLorentzVector &Vect){
//////////////////////////////////////////////////////////////
// printing entire four-vector in one line (with endline)
  for ( int k=0; k < 4 ; k++ )   cout << sw2 << Vect[k];
  cout<<endl;
}
void ROBOL::ReaData(const char *DiskFile, int imax, double xpar[])
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
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//           End of Class ROBOL                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
