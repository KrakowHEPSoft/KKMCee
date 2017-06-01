//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//               Class   TMCgenKKMC                                               //
//                                                                          //
//              Interface (wrapper)  to MC event generator KKMC             //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
#include "TMCgenKKMC.h"


ClassImp(TMCgenKKMC);

TMCgenKKMC::TMCgenKKMC():
  TMCgen()
{
  /// This constructor is for ROOT streamers ONLY
  cout<< "----> TMCgenKKMC Default Constructor (for ROOT only) "<<endl;
  //m_Foam3 = NULL;
}


///_____________________________________________________________
TMCgenKKMC::TMCgenKKMC(const char* Name):
  TMCgen(Name)
{
//! all defaults defined here can be changed by the user
//! before calling TMCgen::Initialize
  //m_Foam3 = NULL;
///////////////////////////////////////////////////
// Physics
  m_NevTot = 0;
  m_EvenCounter = 0;

  cout<< "----> TMCgenKKMC::TMCgenFOAM USER Constructor "<<endl;
}///TMCgenKKMC

///______________________________________________________________________________________
TMCgenKKMC::~TMCgenKKMC()
{
  //!Explicit destructor
  cout<< "----> TMCgenKKMC::TMCgenKKMC !!!! DESTRUCTOR !!!! "<<endl;
}///destructor

///______________________________________________________________________________________
void TMCgenKKMC::Initialize(TRandom *RNgen, ofstream *OutFile, TH1D* h_NORMA)
{
  cout<< "----> TMCgenKKMC::Initialize, Entering "<<endl;
  ///	      SETTING UP RANDOM NUMBER GENERATOR
  TMCgen::Initialize(  RNgen, OutFile, h_NORMA);

//////////////////////////////////////////////////////////////////////////////
//void TMCgenKKMC::Initialize(double ypar[])
//{

  m_NevTot = 0;
  m_EvenCounter = 0;
  m_jmax  = 10000;

// Initialisation input data before generation
  cout << "*******************************" << endl;
  cout << "**   TMCgenKKMC   Initialize **" << endl;
  cout << "*******************************" << endl;

  /////////////////////////////////////////////////////////
  const int jmax =m_jmax;
  ReaData("../../.KK2f_defaults", jmax, m_xpar);  // numbering as in input!!!
  ReaData("./pro.input",         -jmax, m_xpar);  // jmax<0 means no-zeroing
  //double ypar[jmax];
  for(int j=0;j<jmax;j++) m_ypar[j]=m_xpar[j+1];    // ypar has c++ numbering
  //
  double CMSene  = m_xpar[ 1];
  double vvmax   = m_xpar[17];

  cout<<" TMCgen::Initialize: m_CMSene="<<CMSene<<endl;
  cout<<" TMCgen::Initialize: m_vvmax="<<vvmax<<endl;

  //=============================================================
  //   opening disk fime for fortran part of code
  m_out = m_ypar[3];
  const char *output_file = "./pro.output";
  int sl2 = strlen(output_file);
  kk2f_fort_open_(m_out,output_file,sl2);

  //*******************//
  kk2f_initialize_(m_ypar);
  //*******************//

  int NevPrim;
  //kk2f_getxsnormpb_( m_XsNormPb, m_XsErroPb);  //   To be called AFTER CALL KK2f_Finalize !!!!!!
  kk2f_getprimanorma_( m_XsNormPb, NevPrim);     //   Primary Xsection for normalization NANOBARNS ????
  cout<<"/////// TMCgen::Initialize: m_XsNormPb="<<m_XsNormPb<<endl;
  cout<<"/////// TMCgen::Initialize: m_XsErroPb="<<m_XsErroPb<<endl;


  //PyGive("MDCY(113,1)=1;");            // allow rho0 decay
  //PyGive("MDCY(333,1)=0;");            // inhibit phi decay
  //PyGive("MDCY(111,1)=0;");            // inhibit pi0 decay
  //PyGive("MSTP(41)=1;");               // allow all resoanace decay (2->1)
  //PyGive("MSTP(41)=0;");               // inhibit all resoanace decay (2->1)

  PyGive("MSTU(21)=1;");               // no stop due to errors !!!!

  for(int j=1; j<=jmax; j++)  h_NORMA->SetBinContent(j,  m_xpar[j]  );    // xpar encoded
  cout<<" TMCgen::Initialize:  xpar filled into h_NORMA  "<<endl;


  /////////////////////////////////////////////////////////
  if(f_IsInitialized != 0)
	  cout<< "----> TMCgenFOAM::Initialize, already initialized "<<endl;

}// Initialize


///////////////////////////////////////////////////////////////////////////////
void TMCgenKKMC::Generate()
{
  f_NevGen++;
  kk2f_make_();
//  f_TMCgen_NORMA->Fill(-1, m_XsNormPb/1000);    // does not work

  double XsPrimPb; int NevPrim;
  //KKMC_generator->GetPrimaNorma(XsPrim, NevPrim);
  kk2f_getprimanorma_( XsPrimPb, NevPrim);     //   Primary Xsection for normalization NANOBARNS
  f_TMCgen_NORMA->SetBinContent(0,XsPrimPb*NevPrim);  // Pb
  f_TMCgen_NORMA->SetEntries(NevPrim);

}
///////////////////////////////////////////////////////////////////////////////
void TMCgenKKMC::Finalize()
{
// XsNormPb +- XsErroPb is xsection to be readily used to normalize
// histograms both for wt=1 and wt-ed events
  TMCgen::Finalize();

  cout << "*****************************" << endl;
  cout << "**   TMCgenKKMC   Finalize **" << endl;
  cout << "*****************************" << endl;
  //
  kk2f_finalize_();
  //
//  kk2f_getxsnormpb_( m_XsNormPb, m_XsErroPb); // moved to Initialize
  //
  kk2f_fort_close_(m_out);
}//Finalize

/*[[[[
///////////////////////////////////////////////////////////////////////////////
int TMCgenKKMC::GetPyNpart()
{
// provides no. of entries in Lund/Pythia common block
  int npart = cb_PYjets.n;
  //cout<<"KKMC::Make: npart ="<<npart<<endl;
  return npart;
}

///////////////////////////////////////////////////////////////////////////////
void TMCgenKKMC::GetPyParticle( const int j, TPartLund &Event)
{
// Export one particle from /PYJETS/
  Event= TPartLund(j+1,
         cb_PYjets.k[0][j], cb_PYjets.k[1][j],
	 cb_PYjets.k[2][j], cb_PYjets.k[3][j],  cb_PYjets.k[4][j],
	 cb_PYjets.p[0][j], cb_PYjets.p[1][j],
	 cb_PYjets.p[2][j], cb_PYjets.p[3][j],  cb_PYjets.p[4][j],
	 cb_PYjets.v[0][j], cb_PYjets.v[1][j],
	 cb_PYjets.v[2][j], cb_PYjets.v[3][j],  cb_PYjets.v[4][j]);
  //cout<<" here we are!!! j= "<<j<<endl;
  //Event.Print(1);
}
*/

///////////////////////////////////////////////////////////////////////////////
void TMCgenKKMC::Print1()
{
// print event using KKMC
  kk2f_print1_(m_out);
}

/*
///////////////////////////////////////////////////////////////////////////////
void TMCgenKKMC::PyList(int lev)
{
// print event using pythia
  PyGive("MSTU(11)=16");
  pylist_(lev);
  PyGive("MSTU(11)=6");
  pylist_(lev);
}
*/
///////////////////////////////////////////////////////////////////////////////
void TMCgenKKMC::PyGive(const char *directive)
{
// set pythia directive
  int s1;
  s1 = strlen(directive);
  pygive_(directive, s1);
}
///////////////////////////////////////////////////////////////////////////////
void TMCgenKKMC::GetWt( double &WtMain,  double &WtCrude)
{
// get MC weight of event
  kk2f_getwt_(WtMain, WtCrude);
}
///////////////////////////////////////////////////////////////////////////////
void TMCgenKKMC::GetBeams( TLorentzVector &B1,  TLorentzVector &B2)
{
// get 4-momenta of beams
  double p1[4];
  double p2[4];
  kk2f_getbeams_(p1,p2);
  B1.SetPxPyPzE(p1[0],p1[1],p1[2],p1[3]);
  B2.SetPxPyPzE(p2[0],p2[1],p2[2],p2[3]);
}
///////////////////////////////////////////////////////////////////////////////
void TMCgenKKMC::GetFermions( TLorentzVector &F1,  TLorentzVector &F2)
{
// get 4-momenta of final fermions
  double p1[4];
  double p2[4];
  kk2f_getfermions_(p1,p2);
  F1.SetPxPyPzE(p1[0],p1[1],p1[2],p1[3]);
  F2.SetPxPyPzE(p2[0],p2[1],p2[2],p2[3]);
  //for(int i=0;i<4;i++) cout<<B1[i]<<"  "; cout<<endl;
  //for(int i=0;i<4;i++) cout<<B2[i]<<"  "; cout<<endl;
}
///////////////////////////////////////////////////////////////////////////////
void TMCgenKKMC::GetFermKarlud( TLorentzVector &F1,  TLorentzVector &F2)
{
// get 4-momenta of final fermions
  double p1[4];
  double p2[4];
  karlud_getfermions_(p1,p2);
  F1.SetPxPyPzE(p1[0],p1[1],p1[2],p1[3]);
  F2.SetPxPyPzE(p2[0],p2[1],p2[2],p2[3]);
  //for(int i=0;i<4;i++) cout<<B1[i]<<"  "; cout<<endl;
  //for(int i=0;i<4;i++) cout<<B2[i]<<"  "; cout<<endl;
}
///////////////////////////////////////////////////////////////////////////////
void TMCgenKKMC::GetXsecMC( double &xSecPb,  double &xErrPb)
{
// get MC final xsection PICOBARNS
  kk2f_getxsecmc_(xSecPb, xErrPb);
}
///////////////////////////////////////////////////////////////////////////////
void TMCgenKKMC::GetPrimaNorma(double &XsPrim, int &NevPrim)
{
// get normalization elements NANOBARNS
  int NevPrim1;
  kk2f_getprimanorma_( XsPrim, NevPrim1);
  NevPrim=NevPrim1;
}
///////////////////////////////////////////////////////////////////////////////
void TMCgenKKMC::GetPhoton1(const int iphot, TLorentzVector &phot)
{
// get one photon 4-vector from KKMC
  double p1[4];
  int iphot1=iphot;
  kk2f_getphoton1_(iphot1, p1);
  phot.SetPxPyPzE(p1[0],p1[1],p1[2],p1[3]);
}
///////////////////////////////////////////////////////////////////////////////
//void TMCgenKKMC::GetNphot(int &Nphot)
void TMCgenKKMC::GetNphot( int &Nphot)
{
// get photon multiplicity from KKMC
  kk2f_getnphot_( Nphot);
}
///////////////////////////////////////////////////////////////////////////////
double TMCgenKKMC::GetWtAlter(const int id)
{
  double WtAlter;
  int id1 =id;
  kk2f_getwtalter_( id1, WtAlter);
  return WtAlter;
}
///////////////////////////////////////////////////////////////////////////////
void TMCgenKKMC::GetKFfin(int &KF)
{
// get KF code of final fermion
  hepevt_getkffin_( KF);
}

void TMCgenKKMC::ReaData(const char *DiskFile, int imax, double xpar[])
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
}// ReaData


//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//              End of  Class   TMCgenKKMC                                  //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
