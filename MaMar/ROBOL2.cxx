///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                Class ROBOL2                                               //
//                                                                           //
//    It contains Makers for MC generation and Analysis event per event      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
#include "ROBOL2.h"

# define sw2 setprecision(10) << setw(18)

//ClassImp(ROBOL)

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
void ROBOL2::Initialize(long &NevTot)
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
  cout<<"ROBOL2::Initialize:  NevTot = "<<NevTot<<endl;
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
  ///
  int nbv =50;
  hst_vTrueMain = new TH1D("hst_vTrueMain",  "dSig/dvTrue ", nbv, 0.000 ,1.000);
  hst_vTrueCeex2= new TH1D("hst_vTrueCeex2", "dSig/dvTrue ", nbv, 0.000 ,1.000);
  hst_vXGenCeex2= new TH1D("hst_vXGenCeex2", "dSig/dvTrue ", nbv, 0.000 ,1.000);
  hst_vTrueMain->Sumw2();
  hst_vTrueCeex2->Sumw2();
  hst_vXGenCeex2->Sumw2();

  hst_vPhotMain= new TH1D("hst_vPhotMain",  "dSig/dv WTmain ", nbv, 0.000 ,1.000);
  hst_vPhotMain->Sumw2();
  ///
  int nbc =50;
  hst_CosPLCeex2= new TH1D("hst_CosPLCeex2", "dSig/cThetPL  ", nbc, -1.000 ,1.000);
  hst_CosPLCeex2->Sumw2();
  // ln10(sin(theta))
  hst_LnThPhAll = new TH1D("hst_LnThPhAll", "ln10(sin(theta)) all phot.",  60, -6.0 ,0.0);
  hst_LnThPhVis = new TH1D("hst_LnThPhVis", "ln10(sin(theta)) vis. phot.", 60, -6.0 ,0.0);
  hst_LnThPhAll->Sumw2();
  hst_LnThPhVis->Sumw2();
  //  ************* special histo  *************
  HST_KKMC_NORMA = new TH1D("HST_KKMC_NORMA","KKMC normalization &xpar",jmax,0.0,10000.0);
  for(int j=1; j<=jmax; j++)
    HST_KKMC_NORMA->SetBinContent(j,m_xpar[j]);    // xpar encoded
  m_YSum= 0.0;
  m_YSum2=0.0;
}
///////////////////////////////////////////////////////////////////////////////
void ROBOL2::KKMC_NORMA()
{
  // Transfer normalization Record of KKMC to local histogram.
  // For later use in re-normalizing histostograms
  //
  double XsPrim; long NevPrim;
  KKMC_generator->GetPrimaNorma(XsPrim, NevPrim);
  HST_KKMC_NORMA->SetBinContent(0,XsPrim*NevPrim);
  HST_KKMC_NORMA->SetEntries(NevPrim);
  cout<<" ROBOL2::KKMC_NORMA: XsPrim, NevPrim ="<< XsPrim <<"  "<<NevPrim << endl;
}
///////////////////////////////////////////////////////////////////////////////
void ROBOL2::Production(long &iEvent)
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
  double Ebeam  = CMSene/2;
  double Mff    = sqrt(s1);
  double vv     = 1-s1/s;
  // ********************************************************************
  // ***   Photon trigger TrigPho is for everybory, all pions, muons etc
  double Pi=4*atan(1.0);
  double phEne,phTheta,phCosth,phPT,yy;
  double XEneMin = 0.20;  /// Emin for visible photon
  double XTraMin = 0.05;  /// kTmin for visible photon
  double ThetaMin = 15;                   /// theta minimum  for visible photon
  int nph_vis=0;  /// No. of visible (triggered) photons
  double phEneVis = 0;  /// Energy of triggered photon
  int TrigPho=1;  /// =1 for accepted, =0 for rejected
  //**************************************************
  /// Loop over photons
  //**************************************************
  for(iphot=0;iphot<m_Nphot;iphot++){
    phEne   = m_phot[iphot].Energy();
    phCosth = m_phot[iphot].CosTheta();
    phTheta = m_phot[iphot].Theta()*180/Pi;
    phPT    = m_phot[iphot].Pt();
    if( phEne < XEneMin*CMSene/2 ) TrigPho=0;
    if( phPT  < XTraMin*CMSene/2 ) TrigPho=0;
    if( phTheta < ThetaMin       ) TrigPho=0;
    if( phTheta > (180-ThetaMin) ) TrigPho=0;
    if( TrigPho)  nph_vis++;
    if( TrigPho)  phEneVis=phEne;
    /// histogramming, inclusive
    hst_LnThPhAll->Fill( log10(phPT/phEne), WtMain);
    if( TrigPho)  hst_LnThPhVis->Fill( log10(phPT/phEne), WtMain);
  }
  /// triggered photon, only single photon counts
  if(  nph_vis > 0 )hst_nPhVis->Fill(  nph_vis,WtMain);
  double vPhot = phEneVis/Ebeam;
  if( nph_vis == 1 ) hst_vPhotMain->Fill( vPhot, WtMain);

  /// final fermions,  vv, Q^2 costheta, etc
  double CosThe1 = m_pfer1.CosTheta();
  double Theta1  = m_pfer1.Theta();
  double E1      = m_pfer1.Energy();
  double CosThe2 = m_pfer2.CosTheta();
  double Theta2  = m_pfer2.Theta();
  double E2      = m_pfer2.Energy();
//--------------------------------------------------------------------
//** definition cos(theta) of P.L. B219, 103 (1989)
  double CosThePL = ( E1*CosThe1 -E2*CosThe2)/(E1+E2);
  //
  double vvk,x1,x2;
  karlud_getvvxx_(vvk,x1,x2);
  // *********************************************************************
  //          Most of histogramming starts here
  // *********************************************************************
  double WtEEX2  = KKMC_generator->GetWtAlter( 73);    //  Second ord. EEX2 O(alf2)
  double WtEEX3  = KKMC_generator->GetWtAlter( 74);    //  Third order EEX3 O(alf3)
  //cout<< "&&&&&& WtEEX2,3= "<<WtEEX2<<"  "<<WtEEX3<<endl;
  double WtCEEX1 = KKMC_generator->GetWtAlter(202);    //  CEEX Weight O(alf1)
  double WtCEEX2 = KKMC_generator->GetWtAlter(203);    //  CEEX Weight O(alf2)
  //
  hst_nPhAll->Fill(  m_Nphot,WtMain);
  hst_vTrueMain->Fill(       vv, WtMain);
  hst_vTrueCeex2->Fill(      vv, WtCEEX2); // M(2f) of mun pair
  hst_vXGenCeex2->Fill(     vvk, WtCEEX2); // M^star from MC (illegal)
  //
  if(vv<0.9){ hst_CosPLCeex2->Fill(CosThePL, WtCEEX2);}
  // Miscelaneous
  m_YSum  += WtMain;
  m_YSum2 += WtMain*WtMain;
  hst_weight->Fill(WtMain);              // histogramming
  hst_Mff->Fill(Mff,WtMain);             // histogramming
  // debug debug debug debug debug debug debug
  if(iEvent<15){
    cout<<"-----------------------------------------------------------  "<<iEvent;
    cout<<"  -----------------------------------------------------------"<<endl;
    cout<< "vv, vvk         = "<<     vv<<"  "<<vvk<<endl;
    cout<< "CosThe1,CosThe2 = "<<CosThe1<<"  "<<CosThe2<<endl;
    cout<< "CosThePL= "<<CosThePL<<endl;
  }
} //


///////////////////////////////////////////////////////////////////////////////
void ROBOL2::Finalize()
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
void ROBOL2::PartImport(){
/////////////////////////////////////////////////////////////////////////
// Import Pythia common block content into local LuPart matrix
/////////////////////////////////////////////////////////////////////////
  m_Npart= KKMC_generator->GetPyNpart();
  if( (m_Npart<0) || (m_Npart>=4000) ){
    cout<<"++++ ROBOL2::Production: STOP m_Npart= "<<m_Npart<<endl;
    exit(5);
  }
  for(long j=0; j<m_Npart;j++){
    KKMC_generator->GetPyParticle( j, m_Event[j]);  // import one particle
    //m_Event[j].Print(1);
  }
}
///////////////////////////////////////////////////////////////////////////////
long ROBOL2::PartCount(const long flav){
  long jCount=0;
  for(long j=0; j<m_Npart;j++)
    if(m_Event[j].fFlafor == flav){
      jCount++;
    }
  return jCount;
}
///////////////////////////////////////////////////////////////////////////////
long ROBOL2::PartFindAny(const long flav){
// fortran numbering!!!
  long jPosition=0;
  for(long j=0; j<m_Npart;j++)
    if(m_Event[j].fFlafor == flav){
      jPosition=j+1; break;
    }
  return jPosition;
}
///////////////////////////////////////////////////////////////////////////////
long ROBOL2::PartFindStable(const long flav){
// fortran numbering!!!
  long jPosition=0;
  for(long j=0; j<m_Npart;j++)
    if((m_Event[j].fStatus  == 1)&&(m_Event[j].fFlafor == flav)){
      jPosition=j+1; break;
    }
  return jPosition;
}
////////////////////////////////////////////////////////////////////////////////
void ROBOL2::PyPrint(const int mode){
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
void ROBOL2::MomPrint( TLorentzVector &Vect){
//////////////////////////////////////////////////////////////
// printing entire four-vector in one line (with endline)
  for ( int k=0; k < 4 ; k++ )   cout << sw2 << Vect[k];
  cout<<endl;
}
void ROBOL2::ReaData(char DiskFile[], int imax, double xpar[])
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
//           End of Class ROBOL2                                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
