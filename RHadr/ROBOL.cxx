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

TH1D *ROBOL::HistoUP(const char* name, const char* title, int nbins, 
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
void ROBOL::Initialize(KKMC *KKgen, TFile *DiskFile, long &NevTot, int iBook)
{
  //////////////////////////////////////////////////////////////
  //   Initialize MC generator and analysis programs          //
  //////////////////////////////////////////////////////////////
  KKMC_generator = KKgen;
  m_DiskFile     = DiskFile;
  m_iBook        = iBook;
  m_NevGen=0;
  m_count1=0;
  //
  const int jmax =10000;
  ReaData("../../.KK2f_defaults", jmax, m_xpar);  // numbering as in input!!!
  ReaData("./pro.input",         -jmax, m_xpar);  // jmax<0 means no-zeroing
  double ypar[jmax];
  for(int j=0;j<jmax;j++) ypar[j]=m_xpar[j+1];    // ypar has c++ numbering
  //
  NevTot = (long)m_xpar[0];                       // NevTot hidden in xpar[0] !!!
  KKMC_generator->Initialize(ypar);
  cout<<"ROBOL::Initialize:  NevTot = "<<NevTot<<endl;
  //  ************* user histograms  *************
  double CMSene = m_xpar[1];
  int nbin =NevTot/100;
  if(nbin>1000) nbin=1000;

  cout<<"ROBOL iBook= "<<iBook<<endl;

  m_YSum= 0.0;
  m_YSum2=0.0;

  if(iBook){
  //  ************* special histo  *************
  HST_KKMC_NORMA = HistoUP("HST_KKMC_NORMA","KKMC normalization &xpar",jmax,0.0,10000.0);
  for(int j=1; j<=jmax; j++)
    HST_KKMC_NORMA->SetBinContent(j,m_xpar[j]);    // xpar encoded
  //  ******************************************
  hst_weight  = HistoUP("hst_weight" ,  "MC weight",      100, 0.000 , 5.0);
  hst_Mff     = HistoUP("hst_Mff"    ,  "Mass(f-fbar)",  nbin, 0.300 ,CMSene);  
  hst_Q2kloe  = HistoUP("hst_Q2kloe" ,  "Q2 as in KLOE",   43, 0.100 ,0.960); 
  // 
  hst_nPhAll  = HistoUP("hst_nPhAll" , "no. of photons, E>10MeV",  4, -0.5 ,3.5);
  hst_nPhTrg1 = HistoUP("hst_nPhTrg1" ,"no. of photons, trigger1", 4, -0.5 ,3.5);
  hst_nPhTrg2 = HistoUP("hst_nPhTrg2" ,"no. of photons, trigger2", 4, -0.5 ,3.5);
  //
  hst_thPhAll = HistoUP("hst_thPhAll" , "#theta photon, E>10MeV", 180, 0.0 ,180.0);
  hst_thPhTrg1= HistoUP("hst_thPhTrg1", "#theta photon, trigger1", 180, 0.0 ,180.0);
  //
  hst_Q2All   = HistoUP("hst_Q2All" ,   "Q2 phot. unrestricted "     , 43, 0.100 ,0.960);  
  hst_Q2Trig0 = HistoUP("hst_Q2Trig0" , "Q2 2pi, photon.unrestr"     , 43, 0.100 ,0.960);
  hst_Q2Trig1 = HistoUP("hst_Q2Trig1" , "Q2  Egamma>10MeV      "     , 43, 0.100 ,0.960);
  hst_Q2Trig2 = HistoUP("hst_Q2Trig2" , "Q2, Egam>10MeV, PTpi>200MeV", 43, 0.100 ,0.960);
  hst_Q2Trig2b= HistoUP("hst_Q2Trig2b", "Q2, trig2 including backgr.", 43, 0.100 ,0.960);
  //
  hst_vvMuTrg0= HistoUP("hst_vvMuTrg0", "dSig/dv,  muons.         ",100, 0.000 ,1.000);
  hst_vvMuTrg1= HistoUP("hst_vvMuTrg1", "dSig/dv,  muons.  trig1  ",100, 0.000 ,1.000);
  hst_Q2MuTrg0= HistoUP("hst_Q2MuTrg0", "Q2dSig/dQ2, mumu, no cut ", 43, 0.100 ,0.960);
  hst_Q2MuTrg1= HistoUP("hst_Q2MuTrg1", "Q2dSig/dQ2, muons. trig1 ", 43, 0.100 ,0.960);
  //
  hst_Q2dif1   = HistoUP("hst_Q2dif1" ,  "Q2, EEX3-CEEX2  ", 43, 0.100 ,0.960);  
  hst_Q2dif2   = HistoUP("hst_Q2dif2" ,  "Q2, EEX2-EEX3   ", 43, 0.100 ,0.960);  
  hst_Q2dif3   = HistoUP("hst_Q2dif3" ,  "Q2, CEEX1-CEEX2 ", 43, 0.100 ,0.960);
  //
  hst_PiCosthTrg1 = HistoUP("hst_PiCosthTrg1" ,"Cos Theta of pi, with photon cuts ", 40, -1.00,1.00);
  hst_PiCosthTrg2 = HistoUP("hst_PiCosthTrg2" ,"Cos Theta of pi, PTpi>200MeV trig2", 40, -1.00,1.00);
  hst_PiPtTrg1 = HistoUP("hst_PiPtTrg1" ,"PT of pi, with photon cuts ", 50, 0.00,0.50);
  hst_PiPtTrg2 = HistoUP("hst_PiPtTrg2" ,"PT of pi, PTpi>200MeV trig2", 50, 0.00,0.50);
  hst_PivvY= HistoUP("hst_PivvY", "dSig/dv, rho->2pi, trigY ", 50, 0.000 ,1.000);
  hst_Qaux0   = HistoUP("hst_Qaux0" ,   "Q auxiliary "     , 200, 3.095 ,3.098);  
  hst_Qaux1   = HistoUP("hst_Qaux1" ,   "Q auxiliary "     , 200, 3.684 ,3.687);  
  }
  HST_KKMC_NORMA->SetEntries(m_NevPrimOld);
}
///////////////////////////////////////////////////////////////////////////////
void ROBOL::KKMC_NORMA()
{
  // Transfer normalization Record of KKMC to local histogram.
  // For later use in re-normalizing histostograms
  // Code below loook a litle bit artificial but I cannot put
  // HST_KKMC_NORMA->Fill inside rejection lop inside fortran code
  // in order to fill it event per event. Hence I do it for group of events.
  //
  double XsPrim; long NevPrim, NevPrimTot;
  KKMC_generator->GetPrimaNorma(XsPrim, NevPrim);
  NevPrimTot= m_NevPrimOld+NevPrim;
  //
  HST_KKMC_NORMA->SetBinContent(0,XsPrim*NevPrimTot);
  HST_KKMC_NORMA->SetEntries(NevPrimTot);
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
  // VSumPhot  = Total photon 4-vector, from KKMC event, no FSR from Photos!
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
  double phEne,phTheta,phCosth;
  double XCosPho1 = cos(Pi*21.0/180.0); // 2001 KLOE Paper
  double XCosPho2 = cos(Pi* 5.0/180.0); // 2001 KLOE Paper
  double XCosPion = cos(Pi*55.0/180.0); // 2001 KLOE Paper
  double XPtPion  = 0.200;              // 2001 KLOE Paper
  double XEnePho  = 0.010;              // 2001 KLOE Paper
  //
  double YCosPhot = cos(Pi*15.0/180.0); // 2002 New
  double YCosPion = cos(Pi*40.0/180.0); // 2002 New
  double YEnePho  = 0.010;              // 2002 New
  // ********************************************************************
  // ***                  Trigger logics                              ***
  // ********************************************************************
  int TrigPho, nph_ene, nph_trig;
  TLorentzVector *ph_trig1[100];
  TrigPho  =0;
  nph_ene  =0;
  nph_trig =0;
  long iphot;
  for(iphot=0;iphot<m_Nphot;iphot++){
    phEne   = m_phot[iphot].Energy();
    phCosth = m_phot[iphot].CosTheta();
    phTheta = m_phot[iphot].Theta()*180/Pi;
    if(phEne>XEnePho){
      nph_ene++;
      hst_thPhAll->Fill(phTheta,WtMain);   // photon angle before cut
    }
    if( (phEne>XEnePho) && (fabs(phCosth)>XCosPho1) && (fabs(phCosth)<XCosPho2) ){
      nph_trig++;
      ph_trig1[iphot]=&(m_phot[iphot]);
      hst_thPhTrg1->Fill(phTheta,WtMain);  // photon angle after cut
    }
  }
  if(nph_trig) TrigPho = 1;                // one tagged photon is enough
  //
  int TrigNew =1;  // New trigger of Stefano di Falco, Oct 02
  // **********************************************************************
  // **** Pion trigger TrigPiRho, E>200MeV, angle 55degr, no photon tagging
  int NumbPiRho = 0; // No of pions from rho
  double PiCosth =-2;
  double Q2  = 0.0;
  double PT  = -1;
  long   parent, grandparent;
  long   RhoParent=1;
  TLorentzVector VMissing= m_pbea1+m_pbea2; //
  int nPiTrg=0;
  int nPiTot=0;
  // find pi0, works if pi0 is not decayed
  long jPize = PartFindStable(111);    // NOT USED
  if(jPize) nPiTot++;
  // find pi+
  long jPlus = PartFindStable(211);    // fortran numbering!!!
  if(jPlus){
    nPiTot++;
    parent  = m_Event[jPlus-1].fParent;
    grandparent = m_Event[parent-1].fParent;
    m_PiPl  = m_Event[jPlus-1].fMom;
    PT      = m_PiPl.Pt();
    PiCosth = m_PiPl.CosTheta();
    if( m_Event[parent-1].fFlafor      != 113) RhoParent=0; // parent has to be rho
    if( m_Event[grandparent-1].fFlafor == 333) RhoParent=0; // grandpar. is not phi
    if( (PT>XPtPion) && (fabs(PiCosth)<XCosPion) ) nPiTrg++;
    if( m_Event[parent-1].fFlafor == 113 ) NumbPiRho++;
    if( fabs(PiCosth) > YCosPion ) TrigNew =0;  // kill if pi angle outside
    VMissing -= m_PiPl;
  }else{ TrigNew =0; }
  // find pi-
  long jMinu =PartFindStable(-211);   // fortran numbering!!!
  if(jMinu){
    nPiTot++;
    parent  = m_Event[jMinu-1].fParent;
    grandparent = m_Event[parent-1].fParent;
    m_PiMn  = m_Event[jMinu-1].fMom;
    PT      = m_PiMn.Pt();
    PiCosth = m_PiMn.CosTheta();
    if( m_Event[parent-1].fFlafor      != 113) RhoParent=0; // parent has to be rho
    if( m_Event[grandparent-1].fFlafor == 333) RhoParent=0; // grandpar. is not phi
    if( (PT>XPtPion) && (fabs(PiCosth)<XCosPion)  ) nPiTrg++;
    if( m_Event[parent-1].fFlafor == 113 ) NumbPiRho++;
    if( fabs(PiCosth) > YCosPion ) TrigNew =0;  // kill if pi angle outside
    VMissing -= m_PiMn;
  }else{ TrigNew =0; }
  //
  int TrigPiTag = 0; // Any 2 pions tagged
  int TrigPiRho = 0; // Tagged 2 pions from rho
  int TrigPiTwo = 0; // Tagged 2 pions + background alowed
  if( jPlus*jMinu )            TrigPiTag = 1;  // Any 2 pions
  if( nPiTrg==2 )              TrigPiTwo = 1;  // Exactly 2 pions, rho + background(pi+pi-pi0)
  if( nPiTrg==2 && RhoParent ) TrigPiRho = 1;  // Exactly 2 pions with rho parent
  //
  if( RhoParent==0 )        TrigNew  = 0;  // eliminate omega->pi+pi-pi0
  //if( VSumPhot.Energy()         < YEnePho  ) TrigNew =0;  // kill soft
  //if( fabs(VSumPhot.CosTheta()) < YCosPhot ) TrigNew =0;  // kill large angles
  if( VMissing.Energy()         < YEnePho  ) TrigNew =0;  // kill soft
  if( fabs(VMissing.CosTheta()) < YCosPhot ) TrigNew =0;  // kill large angles

  //***********************************
  if(nPiTrg==2) Q2  = (m_PiPl+m_PiMn)*(m_PiPl+m_PiMn);
  if( TrigPiRho && TrigPho && (m_count1<7) && (Q2<0.3) ){
    m_count1++;
    cout<<"**********************************>>> two pions <<<************************************"<<endl;
    cout<<"nPiTrg= "<< nPiTrg<<"  "<<" PT=  "<<" Q2=  "<< Q2<<" s1=  "<< s1 <<endl;
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
  if( (jMu1*jMu1)  && (par1 == par2) && (par1 == 3) ) TrigMu  = 1; // excludeding backgr. of phi
  //**************************************************************
  if( TrigMu && (m_count1<17) ){
    m_count1++;
    cout<<"**********************************>>> two muons <<<************************************"<<endl;
    KKMC_generator->PyList(2);      
  }
  // *********************************************************************
  //          Most of histogramming starts here
  // *********************************************************************
  int Trg0 = TrigPiTag;
  int Trg1 =0;    if( TrigPho && TrigPiTag) Trg1 =1;
  int Trg2 =0;    if( TrigPho && TrigPiRho) Trg2 =1;
  int Trg2b=0;    if( TrigPho && TrigPiTwo) Trg2b=1;
  // Photon multiplicity (integrated x-sections [nb])
  if(Trg0)  hst_nPhAll->Fill( nph_ene,  WtMain);   // histogramming
  if(Trg1)  hst_nPhTrg1->Fill(nph_trig, WtMain);   // histogramming
  if(Trg2)  hst_nPhTrg2->Fill(nph_trig, WtMain);   // histogramming
  // distr. of pion angle and PT
  if(Trg1)  hst_PiCosthTrg1->Fill(PiCosth,WtMain); // histogramming
  if(Trg1)  hst_PiPtTrg1->Fill(   PT,     WtMain); // histogramming
  if(Trg2)  hst_PiCosthTrg2->Fill(PiCosth,WtMain); // histogramming
  if(Trg2)  hst_PiPtTrg2->Fill(   PT,     WtMain); // histogramming
  // Q^2 distrib
  if(NumbPiRho==2) 
            hst_Q2Trig0->Fill( s1, WtMain);        // 2 pi from rho
  if(Trg0)  hst_Q2All->Fill(   s1, WtMain);        // Histogramming
  if(Trg1)  hst_Q2Trig1->Fill( s1, WtMain);        // Histogramming
  if(Trg2)  hst_Q2Trig2->Fill( Q2, WtMain);        // Histogramming
  if(Trg2b) hst_Q2Trig2b->Fill(Q2, WtMain);        // Histogramming
  // muons,  vv, Q^2 distrib
  if(TrigMu           ) hst_vvMuTrg0->Fill(vv, WtMain);    // Muons, no photon trig.
  if(TrigMu && TrigPho) hst_vvMuTrg1->Fill(vv, WtMain);    // Muons, no photon trig.
  if(TrigMu           ) hst_Q2MuTrg0->Fill(s1, WtMain); // Trg0 for Muons
  if(TrigMu && TrigPho) hst_Q2MuTrg1->Fill(s1, WtMain); // Trg1 for Muons
  //if(TrigPho && TrigMu) hst_Q2TrigMu->Fill(s1, s1*WtMain);   // Trg1 for Muons
  //
  double WtEEX2  = KKMC_generator->GetWtAlter( 73);    //  Second ord. EEX2 O(alf2)
  double WtEEX3  = KKMC_generator->GetWtAlter( 74);    //  Third order EEX3 O(alf3)
  //cout<< "&&&&&& WtEEX2,3= "<<WtEEX2<<"  "<<WtEEX3<<endl;
  double WtCEEX1 = KKMC_generator->GetWtAlter(202);    //  CEEX Weight O(alf1)
  double WtCEEX2 = KKMC_generator->GetWtAlter(203);    //  CEEX Weight O(alf2)
  //cout<< "&&&&&& WtCEEX1,2= "<<WtCEEX1<<"  "<<WtCEEX2<<"  WtMain= "<< WtMain <<endl;
  if(Trg2) hst_Q2dif1->Fill(Q2, WtEEX3  -WtMain );  // EEX3  -CEEX2, O(alf^2 L), O(alf^3 L^3)
  if(Trg2) hst_Q2dif2->Fill(Q2, WtEEX2  -WtEEX3 );  // EEX2  -CEEX3, O(alf^3 L^3)
  if(Trg2) hst_Q2dif3->Fill(Q2, WtCEEX1 -WtCEEX2);  // CEEX1 -CEEX2, O(alf^2 L^2), O(alf^3 L^3)
  // New selection
  if(TrigNew){
    hst_PivvY->Fill(vv, WtMain);
    m_YSum  += WtMain;
    m_YSum2 += WtMain*WtMain;
  }
  // Miscelaneous
  hst_weight->Fill(WtMain);              // histogramming
  hst_Mff->Fill(Mff,WtMain);             // histogramming
  hst_Qaux0->Fill(Mff,WtMain);
  hst_Qaux1->Fill(Mff,WtMain);
  //hst_Qaux0->Fill(Mff,1.0);
  //hst_Qaux1->Fill(Mff,1.0);
  //
  double vvk,x1,x2;
  karlud_getvvxx_(vvk,x1,x2);
  hst_Q2kloe->Fill(s*(1-vvk),WtMain);     // histogramming
}
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
  int      nbt  = hst_Q2kloe->GetNbinsX();
  Double_t tmax = hst_Q2kloe->GetXaxis()->GetXmax();
  Double_t tmin = hst_Q2kloe->GetXaxis()->GetXmin();
  Double_t Fact = nbt*XsNormPb/1000/(tmax-tmin)/m_NevGen; // now [nb]
  // **** re-normalized histogram is a clone
  TH1D *HST_Q2kloe =(TH1D*)hst_Q2kloe->Clone();
  HST_Q2kloe->Sumw2();               // is it necessary???
  HST_Q2kloe->SetName("HST_Q2kloe"); // otherwise you get 2 histograms with same name
  HST_Q2kloe->SetTitle("dsigma/dQ2 [nb/GeV^2] ");
  HST_Q2kloe->Scale(Fact);
  // **** alternatively, re-normalized histo can be defined as a new one ****
  //TH1D *HST_Q2kloe =HistoUP("HST_Q2kloe","dsigma/dQ2 [nb/GeV^2]",nbt,tmin,tmax);
  //HST_Q2kloe->Sumw2();
  //HST_Q2kloe->Add(hst_Q2kloe, Fact);
  // *********************************************************************
  // **** alternatively HST_KKMC_NORMA is used at later stage (plotting)
  long   NevPrim = HST_KKMC_NORMA->GetEntries();
  double XsPrima = HST_KKMC_NORMA->GetBinContent(0)/NevPrim;
  cout << "HST_KKMC_NORMA: XsPrima [nb] = "<< XsPrima << " NevPrim= "<< NevPrim <<endl;
  // *********************************************************************
  // New trigger of Stefano di Falco
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
  for(long j=0; j<m_Npart;j++){
    KKMC_generator->GetPyParticle( j, m_Event[j]);  // import one particle
    //m_Event[j].Print(1);
  }
}
///////////////////////////////////////////////////////////////////////////////
long ROBOL::PartCount(const long flav){
  long jCount=0;
  for(long j=0; j<m_Npart;j++)
    if(m_Event[j].fFlafor == flav){
      jCount++;
    }
  return jCount;
}
///////////////////////////////////////////////////////////////////////////////
long ROBOL::PartFindAny(const long flav){
// fortran numbering!!!
  long jPosition=0;
  for(long j=0; j<m_Npart;j++)
    if(m_Event[j].fFlafor == flav){
      jPosition=j+1; break;
    }
  return jPosition;
}
///////////////////////////////////////////////////////////////////////////////
long ROBOL::PartFindStable(const long flav){
// fortran numbering!!!
  long jPosition=0;
  for(long j=0; j<m_Npart;j++)
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
  for(long j=0; j<m_Npart;j++){
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
void ROBOL::ReaData(char DiskFile[], int imax, double xpar[])
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
