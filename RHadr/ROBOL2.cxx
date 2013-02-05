///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                Class ROBOL2                                                //
//                                                                           //
//    It contains Makers for MC generation and Analysis event per event      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "ROBOL2.h"

# define sw2 setprecision(10) << setw(18)

//ClassImp(ROBOL2)

TH1D *ROBOL2::HistoUP(const char* name, const char* title, int nbins, 
	      double xmin, double xmax) {
  TH1D *h;
  if (m_iBook) {
    h = new TH1D(name,title,nbins,xmin,xmax);
    h->Sumw2();
  } else {
    h = (TH1D*)m_DiskFile->Get(name);
  }
  return h;
}

///////////////////////////////////////////////////////////////////////////////
void ROBOL2::Initialize(KKMC *KKgen, TFile *DiskFile, long &NevTot, int iBook)
{
  //////////////////////////////////////////////////////////////
  //   read MC input data for local use.... to be kicked???
  //////////////////////////////////////////////////////////////
  KKMC_generator = KKgen;
  m_DiskFile     = DiskFile;
  m_iBook        = iBook;
  m_NevGen=0;
  m_count1=0;
  const int jmax =10000;
  ReaData("../../.KK2f_defaults", jmax, m_xpar);  // numbering as in input!!!
  ReaData("./pro.input",         -jmax, m_xpar);  // jmax<0 means no-zeroing
  //
  NevTot = (long)m_xpar[0];                       // NevTot hidden in xpar[0] !!!
  double CMSene = m_xpar[1];
  cout<<"ROBOL2::Initialize:  NevTot = "<<NevTot<<endl;
  double smax=1.0;
  if( smax > CMSene*CMSene) smax = 0.99*CMSene*CMSene;
  //  ************* user histograms  *************
  Hst_Q2hadA = HistoUP("Hst_Q2hadA" , "Q2 npi inclusive" ,    100, 0.00 , smax);  
  Hst_Q2piA  = HistoUP("Hst_Q2piA" ,  "Q2 2pi inclusive" ,    100, 0.00 , smax);  
  Hst_Q2piB  = HistoUP("Hst_Q2piB" ,  "Q2 2pi ph. tagged",    100, 0.00 , smax);  
  Hst_Q2muA  = HistoUP("Hst_Q2muA" ,  "Q2 mu inclusive"  ,    100, 0.00 , smax);  
  Hst_piCosA = HistoUP("Hst_piCosA" , "pi costhe inclusive",  100,-1.00 , smax);  
  Hst_piCosB = HistoUP("Hst_piCosB" , "pi costhe tagged   ",  100,-1.00 , smax);  
  Hst_phCosA = HistoUP("Hst_phCosA" , "phot costhe inclusive",100,-1.00 , smax);  
  Hst_phCosB = HistoUP("Hst_phCosB" , "phot costhe tagged   ",100,-1.00 , smax);  
  //
  Hst_Vmu01  = HistoUP("Hst_Vmu01" ,  "v mu inclusive"  ,    100, 0.00 , 1.0);  
  Hst_Vmu71  = HistoUP("Hst_Vmu71" ,  "v mu inclusive"  ,    100, 0.00 , 1.0);  
  Hst_Vmu72  = HistoUP("Hst_Vmu72" ,  "v mu inclusive"  ,    100, 0.00 , 1.0);  
  Hst_Vmu73  = HistoUP("Hst_Vmu73" ,  "v mu inclusive"  ,    100, 0.00 , 1.0);  
  Hst_Vmu74  = HistoUP("Hst_Vmu74" ,  "v mu inclusive"  ,    100, 0.00 , 1.0);  
  //
  Hst_Vmu10  = HistoUP("Hst_Vmu10" ,  "m_beti01"  ,    100, 0.00 , 1.0);  
  Hst_Vmu11  = HistoUP("Hst_Vmu11" ,  "m_beti10"  ,    100, 0.00 , 1.0);  
  Hst_Vmu20  = HistoUP("Hst_Vmu20" ,  "m_beti02"  ,    100, 0.00 , 1.0);  
  Hst_Vmu21  = HistoUP("Hst_Vmu21" ,  "m_beti11"  ,    100, 0.00 , 1.0);  
  Hst_Vmu22  = HistoUP("Hst_Vmu22" ,  "m_beti20"  ,    100, 0.00 , 1.0);
  //
  Hst_Vmu202  = HistoUP("Hst_Vmu202" ,  "CEEX O1"  ,    100, 0.00 , 1.0);
  Hst_Vmu203  = HistoUP("Hst_Vmu203" ,  "CEEX O1"  ,    100, 0.00 , 1.0);
  Hst_Vmu252  = HistoUP("Hst_Vmu252" ,  "CEEX O1"  ,    100, 0.00 , 1.0);
  Hst_Vmu253  = HistoUP("Hst_Vmu253" ,  "CEEX O1"  ,    100, 0.00 , 1.0);
}//ROBOL2::Initialize


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
  // *******************  import event form KKMC  *******************
  m_NevGen++;
  //
  PartImport(); // Import Pythia common block content into local LuPart matrix
  // IMPORT KKMC event and weights
  double WtMain,WtCrude;
  KKMC_generator->GetWt(WtMain,WtCrude);
  KKMC_generator->GetBeams(m_pbea1,m_pbea2);
  KKMC_generator->GetFermions(m_pfer1,m_pfer2);
  KKMC_generator->GetNphot(m_Nphot);                  // photon multiplicity
  TLorentzVector VSumPhot;    // By default all components are initialized by zero. 
  for(long iph=0;iph<m_Nphot;iph++){
    KKMC_generator->GetPhoton1(iph+1,m_phot[iph]);  // photon 4-momenta
    VSumPhot+= m_phot[iph];
  }
  if(iEvent<10){
    cout<<"-----------------------------------------------------------  "<<iEvent;
    cout<<"  -----------------------------------------------------------"<<endl;
    cout<<" VSumPhot= "; MomPrint( VSumPhot );
    KKMC_generator->Print1();
    KKMC_generator->PyList(2);
    PyPrint(1);
  }
  // ************* Generator level variables ****************************
  double s  =(m_pbea1+m_pbea2)*(m_pbea1+m_pbea2);
  double s1 =(m_pfer1+m_pfer2)*(m_pfer1+m_pfer2);
  double CMSene = sqrt(s);
  double Mff    = sqrt(s1);
  double vv     = 1-s1/s;
  // ********************************************************************
  // ***                  Trigger parameters                          ***
  // ********************************************************************
  double Pi=4*atan(1.0);
  double XEnePho  = 0.010;              // 2003 KLOE
  double XCosPho1 = cos(Pi* 0.0/180.0); // 2003 KLOE
  double XCosPho2 = cos(Pi*15.0/180.0); // 2003 KLOE
  double XCosPion = cos(Pi*40.0/180.0); // 2003 KLOE
  double XPtPion  = 0.200;              // 2003 KLOE
  // ********************************************************************
  // ***                  Trigger logics                              ***
  // ********************************************************************
  // ************ Pion triggers ************ 
  int NumbPiRho = 0; // No of pions from rho
  double PiCosth1=0;
  double PiCosth2=0;
  double PiPT1=0;
  double PiPT2=0;
  long   parent, grandparent;
  long   RhoParent=1;
  int nPiTrg=0;
  // find pi0, works if pi0 is not decayed
  long jPize = PartFindStable(111);    // NOT USED
  // ************ find pi+ ************ 
  long jPlus = PartFindStable(211);    // fortran numbering!!!
  if(jPlus){
    parent  = m_Event[jPlus-1].fParent;
    grandparent = m_Event[parent-1].fParent;
    m_PiPl  = m_Event[jPlus-1].fMom;
    PiPT1    = m_PiPl.Pt();
    PiCosth1 = m_PiPl.CosTheta();
    if( m_Event[parent-1].fFlafor      != 113) RhoParent=0; // parent has to be rho
    if( m_Event[grandparent-1].fFlafor == 333) RhoParent=0; // grandpar. is not phi
    nPiTrg++;
    //if( (PiPT>XPtPion) && (fabs(PiCosth)<XCosPion) ) nPiTrg++;
    if( m_Event[parent-1].fFlafor == 113 ) NumbPiRho++;
  }
  // ************ find pi- ************ 
  long jMinu =PartFindStable(-211);   // fortran numbering!!!
  if(jMinu){
    parent  = m_Event[jMinu-1].fParent;
    grandparent = m_Event[parent-1].fParent;
    m_PiMn  = m_Event[jMinu-1].fMom;
    PiPT2    = m_PiMn.Pt();
    PiCosth2 = m_PiMn.CosTheta();
    if( m_Event[parent-1].fFlafor      != 113) RhoParent=0; // parent has to be rho
    if( m_Event[grandparent-1].fFlafor == 333) RhoParent=0; // grandpar. is not phi
    nPiTrg++;
    //if( (PiPT>XPtPion) && (fabs(PiCosth)<XCosPion)  ) nPiTrg++;
    if( m_Event[parent-1].fFlafor == 113 ) NumbPiRho++;
  }
  // ******* incresingly restrictive definition of pi+pi- pair *******
  int TrigPiTag = 0;
  int TrigPiTwo = 0;
  int TrigPiRho = 0;
  if( jPlus*jMinu )            TrigPiTag = 1;  // Any 2 pions
  if( TrigPiTag && nPiTrg==2 ) TrigPiTwo = 1;  // Exactly 2 pions, rho + background(pi+pi-pi0)
  if( TrigPiTwo && RhoParent ) TrigPiRho = 1;  // Exactly 2 pions, only with rho parent
  int TrigPions = TrigPiRho;
  // if( (PiPT1 < XPtPion)           || (PiPT2 < XPtPion)  )           TrigPions =0;
  if( (fabs(PiCosth1) > XCosPion) || (fabs(PiCosth2) > XCosPion)  ) TrigPions =0;
  // ************ Photon part of trigger ****************************
  TLorentzVector VMissing; //
  if( TrigPiTag ) VMissing = m_pbea1 +m_pbea2 -m_PiPl -m_PiMn;
  int TrigPho  = 1;  // photon trigger
  if( VMissing.Energy()         < XEnePho  ) TrigPho =0;  // kill soft
  if( fabs(VMissing.CosTheta()) < XCosPho2 ) TrigPho =0;  // kill large angles

  // ****************************************************************
  double Q2  = -1e-5;
  if( TrigPiTag )  Q2= (m_PiPl+m_PiMn)*(m_PiPl+m_PiMn);
  if( TrigPiRho && TrigPho && (m_count1<7) && (Q2<0.3) ){
    m_count1++;
    cout<<"**********************************>>> two pions <<<************************************"<<endl;
    cout<<"nPiTrg= "<< nPiTrg<<"  "<<" PiPT=  "<<" Q2=  "<< Q2<<" s1=  "<< s1 <<endl;
    cout<<"  m_PiPl=  "; MomPrint( m_PiPl );
    cout<<"  m_PiMn=  "; MomPrint( m_PiMn );
    KKMC_generator->PyList(2);
    //PyPrint(1);
  }
  // ****************************************************************
  // ****  Two Muons trigger, background excluded, no photon tagging
  int TrigMu  = 0;
  // find muons, excluding muons from phi decays!!!
  long jMu1 =PartFindStable( 13);    // fortran numbering!!!
  long jMu2 =PartFindStable(-13);    // fortran numbering!!!
  long par1=m_Event[jMu1-1].fParent;
  long par2=m_Event[jMu2-1].fParent;
  if( (jMu1*jMu1)  && (par1 == par2) && (par1 == 3) ) TrigMu  = 1; // excluding backgr. of phi
  //**************************************************************
  if( TrigMu && (m_count1<17) ){
    m_count1++;
    cout<<"**********************************>>> two muons <<<************************************"<<endl;
    KKMC_generator->PyList(2);      
  }
  // *********************************************************************
  //          Most of histogramming starts here
  // *********************************************************************
  int SeleA = TrigPiRho;          // inclusive 2pi
  int SeleB = TrigPions*TrigPho;  // 2pi in the detector + cut on photon

  if(TrigMu==0){
    Hst_Q2hadA->Fill(s1, WtMain); // Any n pions
   }
 
  if(SeleA){
    Hst_Q2piA->Fill(Q2, WtMain);
    Hst_piCosA->Fill(PiCosth1, 0.5*WtMain);
    Hst_piCosA->Fill(PiCosth2, 0.5*WtMain);
    Hst_phCosA->Fill(VMissing.CosTheta(), WtMain);
  }
  if(SeleB){
    Hst_Q2piB->Fill(Q2, WtMain);
    Hst_piCosB->Fill(PiCosth1, 0.5*WtMain);
    Hst_piCosB->Fill(PiCosth2, 0.5*WtMain);
    Hst_phCosB->Fill(VMissing.CosTheta(), WtMain);
  }
  //double WtA; 
  //long   idw;
  if(TrigMu){
    Hst_Q2muA->Fill(s1, WtMain);
    Hst_Vmu01->Fill(vv, KKMC_generator->GetWtAlter( 1));
    Hst_Vmu71->Fill(vv, KKMC_generator->GetWtAlter(71));
    Hst_Vmu72->Fill(vv, KKMC_generator->GetWtAlter(72));
    Hst_Vmu73->Fill(vv, KKMC_generator->GetWtAlter(73));
    Hst_Vmu74->Fill(vv, KKMC_generator->GetWtAlter(74));
    Hst_Vmu10->Fill(vv, KKMC_generator->GetWtAlter(10));
    Hst_Vmu11->Fill(vv, KKMC_generator->GetWtAlter(11));
    Hst_Vmu20->Fill(vv, KKMC_generator->GetWtAlter(20));
    Hst_Vmu21->Fill(vv, KKMC_generator->GetWtAlter(21));
    Hst_Vmu22->Fill(vv, KKMC_generator->GetWtAlter(22));
    Hst_Vmu202->Fill(vv, KKMC_generator->GetWtAlter(202)); // CEEX O1, with IFI
    Hst_Vmu203->Fill(vv, KKMC_generator->GetWtAlter(203)); // CEEX O2, with IFI
    Hst_Vmu252->Fill(vv, KKMC_generator->GetWtAlter(252)); // CEEX O1, with IFI
    Hst_Vmu253->Fill(vv, KKMC_generator->GetWtAlter(253)); // CEEX O2, with IFI
  }
}//ROBOL2::Production


///////////////////////////////////////////////////////////////////////////////
void ROBOL2::Finalize()
{
//   Finalize MC  run, final printouts, cleaning etc.,
//   Plotting histograms is done independently using root file
  cout << "ROBOL2::Finalize: m_NevGen = "<< m_NevGen<<endl;
  double XsNormPb, XsErroPb;
  KKMC_generator->Finalize(XsNormPb, XsErroPb);
  cout << "ROBOL2::Finalize: XsNormPb [pb] = "<<  XsNormPb << "  +-  "<< XsErroPb <<endl;
}//ROBOL2::Finalize


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
