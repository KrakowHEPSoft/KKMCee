//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//               Class   KKplot                                             //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
// KKplot is multipurpose toolbox for KKMC testing.
//  1. Interfaces (erappers) to KKMC and KKsem F77 subrograms
//  2. Integrand for Foam
//  3. A few routines for producing tatex table out of histograms
//////////////////////////////////////////////////////////////////////////////

#include "KKplot.h"

KKplot::KKplot(const char* Name)
{
		cout<< "----> KKplot USER Constructor "<<endl;
  m_jmax =10000;
}

///////////////////////////////////////////////////////////////////////////////////
void KKplot::Initialize(TFile &DiskFileA){
//------------------------------------------------------------------------
//  For single proceser MC runs and xpar stored in normalization histogram
//------------------------------------------------------------------------
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  cout<<"==================================================================="<<endl;
  cout<<"================ KKplot::initialization1 ==========================="<<endl;
  int nbin_NORM = HST_KKMC_NORMA->GetNbinsX();
  cout<< "No fo bins in HST_KKMC_NORMA = "<< nbin_NORM <<endl;
  if( m_jmax != nbin_NORM){
     cout<<"KKplot::Initialize: wrong no. of bins in HST_KKMC_NORMA"<<endl;
     exit(9);
  }//if
  double xpar[10001];
  for(int j=1; j<=m_jmax; j++) xpar[j]=HST_KKMC_NORMA->GetBinContent(j);
  // The trick to find out no. of processors in farm
  int Nfarm = xpar[511];  // normally xpar[511]=1
  cout<< "//////////  Nfarm = "<< Nfarm <<endl;
  for(int j=1; j<=m_jmax; j++) xpar[j] = xpar[j]/Nfarm;
  //-------------------
  Initialize(xpar);
}//Initialize

///////////////////////////////////////////////////////////////////////////////////
void KKplot::Initialize( double xpar[]){
//------------------------------------------------------------------------
//  xpar importent from outside, for multi-procesor MC runs (farming)
//------------------------------------------------------------------------
  cout<<"================ KKplot::initialization2 ==========================="<<endl;
  if( m_jmax != 10000) cout<<"++++KKplot::Initialize jmax ="<<m_jmax<<endl;
  if( m_jmax != 10000) exit(9);
  for(int j=1; j<=m_jmax; j++) m_ypar[j-1]= xpar[j];
  //
  int iofset =1;   // ???
  //iofset = 0;      // new
  m_CMSene  = m_ypar[ 1 -iofset];
  m_vvmax   = m_ypar[17 -iofset];
  m_alfinv  = m_ypar[30 -iofset];     // 1/alphaQED at Q^2=0
  m_gnanob  = m_ypar[31 -iofset];     // GeV^2 -> nanobarns
  cout<< " m_CMSene = "<<m_CMSene<<endl;
  cout<< " m_alfinv = ypar[ 30 -iofset]= " <<  m_alfinv<< endl;
  cout<< " m_gnanob = ypar[ 31 -iofset]= " <<  m_gnanob<< endl;
  int     KFini = m_ypar[400 -1];
  double     MZ = m_ypar[502 -1];
  cout<< "    KFini = ypar[400 -iofset]= " <<  KFini<< endl;
  cout<< "       MZ = ypar[502 -iofset]= " <<  MZ << endl;
//[[[
  cout<< "  ircdat= ypar[ 74 -iofset]= " <<  m_ypar[ 74 -1] << endl;

  m_pi      = 3.1415926535897932;
  m_ceuler  = 0.57721566;
//
  m_alfpi   = 1/m_alfinv/m_pi;
//
  m_beam    = 0.510999e-3;  // electron
  m_chini   = 1.0;          // electron

  m_fin     = 0.105;        // final ferm. muon mass
  m_chfin   = 1.0;          // final ferm. muon charge

  m_KFini   = 11; // electron
  m_KFf     = 13; // muon

  cout<<"***********************KKsem::initialize****************************"<<endl;
  for(int j=0;j<30;j++)
    cout<<j+1<<"   "<<m_ypar[j]<<endl;
  //
  const char *output_file = "./kksem.output";
  int stl2 = strlen(output_file);
  int mout =16;
  kk2f_fort_open_(mout,output_file,stl2);
  kk2f_initialize_(m_ypar);
  kksem_initialize_(m_ypar);
  cout<<"================ KKsem initialization END   ==========================="<<endl;
  cout<<"==================================================================="<<endl;
  //int kdum; kk2f_getkeyfsr_(kdum); //<-- to avoid linker bug! (if no kk2f_initialize)
}

///////////////////////////////////////////////////////////////////////////////////
void KKplot::VVplot( TH1 *hstNew, int KF, char chak[5], int KeyDis, int KeyFob)
{
  int   nbin = hstNew->GetNbinsX();
  double xmin = hstNew->GetXaxis()->GetXmin();
  double xmax = hstNew->GetXaxis()->GetXmax();
  double Bin[nbin];
  //
  kksem_setkffin_(KF); // set m_KFfin in KKsem
  bornv_setkf_( KF );  // set SINGLE Final State
  //
  kksem_setkeyfob_(KeyFob);
  kksem_vvplot_vec_(KeyDis,chak,nbin,xmin,xmax,Bin);
  //
  hstNew->Reset();
  for(int ib=0;ib<nbin;ib++){
    //cout<<"KKplot::VVplot: ib= "<<ib<<"  Bin(ib) =  "<<Bin[ib]<<endl;
    hstNew->SetBinContent(ib+1, Bin[ib]);
    hstNew->SetBinError(  ib+1, 0.0);
  }
}// KKplot::VVplot


///////////////////////////////////////////////////////////////////////////////////
void KKplot::Ord1fill( TH1 *hstNew, int KeyDis)
// AFB(vmax)
{
  int    nbin = hstNew->GetNbinsX();
  double vmin = hstNew->GetXaxis()->GetXmin();
  double vmax = hstNew->GetXaxis()->GetXmax();
  double vv, Result;
  //
  hstNew->Reset();
 // for(int ib=1; ib <= nbin ; ib++) {
  for(int ib=1; ib <= nbin ; ib++) {
    vv = vmin + (ib*(vmax-vmin))/nbin;  // RHS of the bin
	kksem_afb_calc_( KeyDis, m_KFini, m_KFf, m_CMSene, vv, Result);
    if(ib <4) cout<<"KKplot::Ord1fill: KeyDis,ib= "<<KeyDis<<", "<<ib<<"   vv="<<vv<<"  Result =  "<<Result<<endl;
    hstNew->SetBinContent(ib, Result);
    hstNew->SetBinError(  ib, 0.0);
  }
}// KKplot::Ord1fill



///////////////////////////////////////////////////////////////////////////////////
void KKplot::Cplot( TH1 *hstNew,
		   int KF, char chak[5], int KeyDis, int KeyFob, double vmin, double vmax)
{
  int   nbc = hstNew->GetNbinsX();
  double cmin = hstNew->GetXaxis()->GetXmin();
  double cmax = hstNew->GetXaxis()->GetXmax();
  double cBin[nbc];
  int   nbv = 1;   // only one bin in v needed
  double vBin[100]; // only one bin in v needed
  //
  kksem_setkffin_(KF); // set m_KFfin in KKsem
  bornv_setkf_( KF );  // set SINGLE Final State
  //
  if( (KeyFob == 0) ||(KeyFob == -100)  ){
    kksem_setkeyfob_(KeyFob);
  }else{
    cout<<"++++ KKsem::Cplot wrong KeyFob= "<< KeyFob<<endl;
    exit(7);
  }
  //
  double c1,c2,dc,fnor,Dsig;
  dc = (cmax-cmin)/nbc;
  hstNew->Reset();
  for(int ib=0;ib<nbc;ib++){
    c1= cmin +ib*dc;
    c2= c1+dc;
    kksem_setcrange_(c1,c2);
    kksem_vvplot_vec_(KeyDis,chak,nbv,vmin,vmax,vBin); // vBin has only one bin!
    Dsig=vBin[0]/dc*(vmax-vmin); // integration over v!
    cout<<"KKplot::Cplot: ib= "<<ib<<" cBin(ib) =  "<<Dsig<<endl;
    hstNew->SetBinContent(ib+1, Dsig);
    hstNew->SetBinError(  ib+1, 0.0);
  }
  kksem_setcrange_(-1.0,1.0); // back to normal
}// KKplot::VVplot


void KKplot::ReaData(const char DiskFile[], int imax, double xpar[])
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




///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                End of Class KKplot                                        //
///////////////////////////////////////////////////////////////////////////////

