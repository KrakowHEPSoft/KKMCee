///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                Class ROBOL                                                //
//                                                                           //
//    It contains Makers for MC generation and Analysis event per event      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
#include "RoboFoam.h"

# define sw2 setprecision(10) << setw(18)

///////////////////////////////////////////////////////////////////////////////
//      *************** temporary entries from KKMC ****************
//      SUBROUTINE KarLud_GetVVxx(vv,x1,x2)
extern "C" void  karlud_getvvxx_(double&, double&, double&);
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
void RoboFoam::Initialize(long &NevTot)
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
  cout<<" RoboFoam::Initialize:  NevTot = "<<NevTot<<endl;
  //  ************* user histograms  *************

  hst_weight3 = new TH1D("hst_weight3" ,  "MC weight",      100, -1.0, 2.0);
  hst_weight3->Sumw2();
  hst_weight5 = new TH1D("hst_weight5" ,  "MC weight",      100, -1.0, 2.0);
  hst_weight5->Sumw2();

  // scatergrams
  int nbv = 50;
  HST_xx_Ceex2  = new TH1D("HST_xx_Ceex2" ,   "dSig/dv",   nbv, 0.0, 1.0);
  HST_xx_Ceex2->Sumw2();
  HST_xx_Ceex2n = new TH1D("HST_xx_Ceex2n" ,  "dSig/dv",   nbv, 0.0, 1.0);
  HST_xx_Ceex2n->Sumw2();
  nbv = 50;
  int nbc = 50;
  SCA_xc_Ceex2 = new TH2D("SCA_xc_Ceex2",   "dSig/dc/dv ", nbv, 0.0 ,1.0, nbc, -1.0 ,1.0);
  SCA_xc_Ceex2->Sumw2();
  SCA_xc_Ceex2n= new TH2D("SCA_xc_Ceex2n",  "dSig/dc/dv ", nbv, 0.0 ,1.0, nbc, -1.0 ,1.0);
  SCA_xc_Ceex2n->Sumw2();

  //  New bigger scatergrams, restricted vmax
  int NBv =100; int NBc = 100;
  double vmx2= 0.20;
  //sct_vAcPR_Ceex2= new TH2D("sct_vAcPR_Ceex2",  "dSig/dc/dv ", NBv, 0.0 ,vmx2, NBc, -1.0 ,1.0);
  //sct_vAcPR_Ceex2->Sumw2();
  //sct_vAcPR_Ceex2n= new TH2D("sct_vAcPR_Ceex2n","dSig/dc/dv ", NBv, 0.0 ,vmx2, NBc, -1.0 ,1.0);
  //sct_vAcPL_Ceex2n->Sumw2();

  //************* special normalization histos  *************
  HST_FOAM_NORMA3 = new TH1D("HST_KKMC_NORMA3","KKMC normalization &xpar",jmax,0.0,10000.0);
  HST_FOAM_NORMA5 = new TH1D("HST_KKMC_NORMA5","KKMC normalization &xpar",jmax,0.0,10000.0);
  for(int j=1; j<=jmax; j++){
    HST_FOAM_NORMA3->SetBinContent(j,m_xpar[j]);    // xpar encoded
    HST_FOAM_NORMA5->SetBinContent(j,m_xpar[j]);    // xpar encoded
  }
  //****************** Setting up Foam generators **************************
  // Interface to KKfoam providing Foam integrand
  LibSem = new KKfoam();
  LibSem->Initialize(ypar);
  //------------------------------------------------------------------
  PseRan   = new TRandom3();  // Create random number generator
  PseRan->SetSeed(4357);
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%% FOAM simulators/integrators %%%%%%%%%%%%%%%
  // ISR+FSR+IFI
  MC_Gen5 = new TFoam("MC_Gen5");   // Create Simulator
  MC_Gen5->SetkDim(5);         // No. of dimensions, obligatory!
  MC_Gen5->SetnCells( 10000);  // No. of cells, can be omitted, default=2000
  MC_Gen5->SetnSampl(100000);  // No. of MC evts/cell in exploration, default=200
  MC_Gen5->SetRho(LibSem);
  MC_Gen5->SetPseRan(PseRan);  // Set random number generator, mandatory!
  MC_Gen5->SetOptRej(0);       // wted events (=0), default wt=1 events (=1)
  LibSem->m_Mode = 5;           // Choose Density5
  MC_Gen5->Initialize();       // Initialize simulator, may take time...
  double Xsav5, dXsav5;          //  For renormalizing histograms
  MC_Gen5->GetIntNorm(m_Xsav5,dXsav5);
  HST_FOAM_NORMA5 = new TH1D("HST_FOAM_NORMA5","MC normalization",10,0.0,10000.0);
  // ISR+FSR only
  MC_Gen3 = new TFoam("MC_Gen3");   // Create Simulator
  MC_Gen3->SetkDim(3);         // No. of dimensions, obligatory!
  MC_Gen3->SetnCells( 10000);  // No. of cells, can be omitted, default=2000
  MC_Gen3->SetnSampl(100000);  // No. of MC evts/cell in exploration, default=200
  MC_Gen3->SetRho(LibSem);
  MC_Gen3->SetPseRan(PseRan);  // Set random number generator, mandatory!
  MC_Gen3->SetOptRej(0);       // wted events (=0), default wt=1 events (=1)
  LibSem->m_Mode = 3;           // Choose Density3
  LibSem->m_count =0;           // resetting debug counter
  MC_Gen3->Initialize();       // Initialize simulator, may take time...
  double Xsav3, dXsav3;          //  For renormalizing histograms
  MC_Gen3->GetIntNorm(m_Xsav3,dXsav3);
  HST_FOAM_NORMA3 = new TH1D("HST_FOAM_NORMA3","MC normalization",10,0.0,10000.0);

}//RoboFoam::Initialize

///////////////////////////////////////////////////////////////////////////////
void RoboFoam::Production(long &iEvent)
{
/////////////////////////////////////////////////////////////////////////
//
//   GENERATE AND ANALYZE SINGLE EVENT
//
/////////////////////////////////////////////////////////////////////////
// ****************************************************************
// ************ Generate event and import it here  ****************
  m_NevGen++;
  double wt3,wt5,xx,CosTheta;
/// Generate ISR+FSR event ISR+FSR
  LibSem->m_Mode = 3;
  MC_Gen3->MakeEvent();            // generate MC event
  MC_Gen3->GetMCwt(wt3);
  xx  = LibSem->m_xx;
  CosTheta = LibSem->m_CosTheta;
  hst_weight3->Fill(wt3,1.0);
  HST_xx_Ceex2n->Fill(xx,wt3);
  SCA_xc_Ceex2n->Fill(xx,CosTheta,wt3);
  HST_FOAM_NORMA3->Fill(-1.0,m_Xsav3);  // fill normalization into underflow

/// Generate ISR+FSR+IFI event
  LibSem->m_Mode = -5;
  MC_Gen5->MakeEvent();            // generate MC event
  MC_Gen5->GetMCwt(wt5);
  xx  = LibSem->m_xx;
  CosTheta = LibSem->m_CosTheta;
  hst_weight5->Fill(wt5,1.0);
  SCA_xc_Ceex2->Fill(xx,CosTheta,wt5);
  HST_FOAM_NORMA5->Fill(-1.0,m_Xsav5);  // fill normalization into underflow

  // debug debug debug debug debug debug debug
  if(iEvent<15){
    cout<<"-------------------------------  "<<iEvent;
    cout<<"  -------------------------------"<<endl;
    cout<< " xx = "<<xx<< "      CosTheta = "<<CosTheta<<endl;
  }
} //


///////////////////////////////////////////////////////////////////////////////
void RoboFoam::KKMC_NORMA()
{
  // Transfer normalization Record of KKMC to local histogram.
  // For later use in re-normalizing histostograms
  HST_FOAM_NORMA3->SetBinContent(0,m_Xsav3*m_NevGen);
  HST_FOAM_NORMA3->SetEntries(m_NevGen);
  HST_FOAM_NORMA5->SetBinContent(0,m_Xsav5*m_NevGen);
  HST_FOAM_NORMA5->SetEntries(m_NevGen);
  cout<<" RoboFoam::KKMC_NORMA: m_Xsav5, m_NevGen ="<< m_Xsav5 <<"  "<< m_NevGen << endl;
}// RoboFoam::KKMC_NORMA



///////////////////////////////////////////////////////////////////////////////
void RoboFoam::Finalize()
{
//   Finalize MC  run, final printouts, cleaning etc., xcheck of normalization
  double XsNormPb, XsErroPb;
  cout << "///////////////////////////Finalize////////////////////////////////////"<<endl;
}



///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                             UTILITIES                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


void RoboFoam::ReaData(char DiskFile[], int imax, double xpar[])
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
