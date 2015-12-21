///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                Class RoboFSR                                              //
//                                                                           //
//    It contains Makers for MC generation and Analysis event per event      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
#include "RoboFSR.h"

# define sw2 setprecision(10) << setw(18)

//ClassImp(RoboFSR)

///////////////////////////////////////////////////////////////////////////////
//      *************** temporary entries from KKMC ****************
//      SUBROUTINE KarLud_GetVVxx(vv,x1,x2)
extern "C" void  karlud_getvvxx_(double&, double&, double&);
extern "C" void  pyhepc_(long&);
extern "C" void  photos_(long&);
extern "C" void  phoini_();
extern "C" void  hepevt_setphotosflagtrue_(long&);
extern "C" void  hepevt_getnhep_(long&);
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
void RoboFSR::Initialize(long &NevTot)
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
  cout<<"RoboFSR::Initialize:  NevTot = "<<NevTot<<endl;
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
  int nbv =50;
  hst_vTrueMain = new TH1D("hst_vTrueMain",  "dSig/dvTrue ", nbv, 0.000 ,1.000);
  hst_vTrueCeex2= new TH1D("hst_vTrueCeex2", "dSig/dvTrue ", nbv, 0.000 ,1.000);
  hst_vAlepCeex2= new TH1D("hst_vAlepCeex2", "dSig/dvTrue ", nbv, 0.000 ,1.000);
  hst_vXGenCeex2= new TH1D("hst_vXGenCeex2", "dSig/dvTrue ", nbv, 0.000 ,1.000);
  hst_vTrueMain->Sumw2();
  hst_vTrueCeex2->Sumw2();
  hst_vAlepCeex2->Sumw2();
  hst_vXGenCeex2->Sumw2();
  int nbs =50;
  hst_s1Ceex2  = new TH1D("hst_s1Ceex2",  "dSig/ds1 ", nbs, 86.0 , 96.0);
  hst_svk      = new TH1D("hst_svk",      "dSig/ds1 ", nbs, 86.0 , 96.0);
  hst_s1Ceex2->Sumw2();
  hst_svk->Sumw2();
  //  ******************** NEW for LHC ***********************
  int nbm =100;
  hst_M100mu  = new TH1D("hst_M100mu",  "dSig/dMll ", nbm, 60 , 160);
  hst_M100mu->Sumw2();

  //  ************* special histo  *************
  HST_KKMC_NORMA = new TH1D("HST_KKMC_NORMA","KKMC normalization &xpar",jmax,0.0,10000.0);
  for(int j=1; j<=jmax; j++)
    HST_KKMC_NORMA->SetBinContent(j,m_xpar[j]);    // xpar encoded
  m_YSum= 0.0;
  m_YSum2=0.0;
}
///////////////////////////////////////////////////////////////////////////////
void RoboFSR::KKMC_NORMA()
{
  // Transfer normalization Record of KKMC to local histogram.
  // For later use in re-normalizing histostograms
  //
  double XsPrim; long NevPrim;
  KKMC_generator->GetPrimaNorma(XsPrim, NevPrim);
  HST_KKMC_NORMA->SetBinContent(0,XsPrim*NevPrim);
  HST_KKMC_NORMA->SetEntries(NevPrim);
  cout<<" RoboFSR::KKMC_NORMA: XsPrim, NevPrim ="<< XsPrim <<"  "<<NevPrim << endl;
}
///////////////////////////////////////////////////////////////////////////////
void RoboFSR::Production(long &iEvent)
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
  long iphot,iphot1;
  for(iphot=0;iphot<m_Nphot;iphot++){
    KKMC_generator->GetPhoton1(iphot+1,m_phot[iphot]);  // photon 4-momenta
    VSumPhot+= m_phot[iphot];
  }
  if(iEvent<10){
    cout<<"-----------------------------------------------------------  "<<iEvent;
    cout<<"  -----------------------------------------------------------"<<endl;
    cout<<" VSumPhot= "; MomPrint( VSumPhot );
    KKMC_generator->Print1();
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
  // ***   Photon trigger TrigPho is for everybody, all pions, muons etc
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
  long jMu1 =PartFindStable( 13);    // fortran numbering!!!
  long jMu2 =PartFindStable(-13);    // fortran numbering!!!
  m_pMu1  = m_Event[jMu1-1].fMom;    // fortran numbering!!!
  m_pMu2  = m_Event[jMu2-1].fMom;    // fortran numbering!!!
  long par1=m_Event[jMu1-1].fParent; // fortran numbering!!!
  long par2=m_Event[jMu2-1].fParent; // fortran numbering!!!
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
  double vvk,x1,x2;
  karlud_getvvxx_(vvk,x1,x2);
  // *********************************************************************
  //          Most of histogramming starts here
  // *********************************************************************
  double WtEEX2  = KKMC_generator->GetWtAlter( 73);    //  Second ord. EEX2 O(alf2)
  double WtEEX3  = KKMC_generator->GetWtAlter( 74);    //  Third order EEX3 O(alf3)
  //cout<< "&&&&&& WtEEX2,3= "<<WtEEX2<<"  "<<WtEEX3<<endl;
  double WtCEEX1 = KKMC_generator->GetWtAlter(252);    //  CEEX O(alf1) Interf. off
  double WtCEEX2 = KKMC_generator->GetWtAlter(253);    //  CEEX O(alf2) Interf. off
  //
  hst_nPhAll->Fill(  m_Nphot,WtMain);
  hst_nPhVis->Fill(  nph_ene,WtMain);
  hst_vTrueMain->Fill(       vv, WtMain);
  hst_vAlepCeex2->Fill(1-zAleph, WtCEEX2); // M^star guessed
  //
  //hst_vTrueCeex2->Fill(      vv, WtCEEX2); // M(2f) of mun pair
  //hst_vXGenCeex2->Fill(     vvk, WtCEEX2); // M^star from MC (illegal)
  hst_vTrueCeex2->Fill(      vv, WtEEX3); // M(2f) of mun pair
  hst_vXGenCeex2->Fill(     vvk, WtEEX3); // M^star from MC (illegal)
  // NEW!!!
  double svk= s*(1-vvk);
  //hst_s1Ceex2->Fill(   sqrt(s1), WtCEEX2); // M^2(2f) of mun pair
  //hst_svk->Fill(      sqrt(svk), WtCEEX2);
  hst_s1Ceex2->Fill(   sqrt(s1), WtEEX3); // M^2(2f) of mun pair
  hst_svk->Fill(      sqrt(svk), WtEEX3);
  /// ********************** NEW **********************
  hst_M100mu->Fill(   sqrt(s1), WtEEX3); // M^2(2f) of muon pair
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
void RoboFSR::Finalize()
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
  long   NevPrim = HST_KKMC_NORMA->GetEntries();
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
void RoboFSR::PartImport(){
/////////////////////////////////////////////////////////////////////////
// Import Pythia common block content into local LuPart matrix
/////////////////////////////////////////////////////////////////////////
  m_Npart= KKMC_generator->GetPyNpart();
  if( (m_Npart<0) || (m_Npart>=4000) ){
    cout<<"++++ RoboFSR::Production: STOP m_Npart= "<<m_Npart<<endl;
    exit(5);
  }
  for(long j=0; j<m_Npart;j++){
    KKMC_generator->GetPyParticle( j, m_Event[j]);  // import one particle
    //m_Event[j].Print(1);
  }
}
///////////////////////////////////////////////////////////////////////////////
long RoboFSR::PartCount(const long flav){
  long jCount=0;
  for(long j=0; j<m_Npart;j++)
    if(m_Event[j].fFlafor == flav){
      jCount++;
    }
  return jCount;
}
///////////////////////////////////////////////////////////////////////////////
long RoboFSR::PartFindAny(const long flav){
// fortran numbering!!!
  long jPosition=0;
  for(long j=0; j<m_Npart;j++)
    if(m_Event[j].fFlafor == flav){
      jPosition=j+1; break;
    }
  return jPosition;
}
///////////////////////////////////////////////////////////////////////////////
long RoboFSR::PartFindStable(const long flav){
// fortran numbering!!!
  long jPosition=0;
  for(long j=0; j<m_Npart;j++)
    if((m_Event[j].fStatus  == 1)&&(m_Event[j].fFlafor == flav)){
      jPosition=j+1; break;
    }
  return jPosition;
}
////////////////////////////////////////////////////////////////////////////////
void RoboFSR::PyPrint(const int mode){
//
// PRINT entire Pythia EVENT
//
  TLorentzVector Sum;
  cout<<"  "<<endl;
  cout<<"lser status flavor parent child1 child2";
  cout<<"                Px                Py                Pz               Ene";
  cout<<"              Mass"<<endl;
  for(long j=0; j<m_Npart;j++){
    m_Event[j].Print(mode);
    if(m_Event[j].fStatus == 1) Sum += m_Event[j].fMom;
  }
  cout<<"                  Total 4-momentum --> ";
  MomPrint(Sum);
}
void RoboFSR::MomPrint( TLorentzVector &Vect){
//////////////////////////////////////////////////////////////
// printing entire four-vector in one line (with endline)
  for ( int k=0; k < 4 ; k++ )   cout << sw2 << Vect[k];
  cout<<endl;
}
void RoboFSR::ReaData(char DiskFile[], int imax, double xpar[])
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
//           End of Class RoboFSR                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
