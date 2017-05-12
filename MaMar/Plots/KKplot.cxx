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
}

///////////////////////////////////////////////////////////////////////////////////
void KKplot::Initialize(TFile &DiskFileA){
  //------------------------------------------------------------------------
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  cout<<"==================================================================="<<endl;
  cout<<"================ KKsem initialization begin ==========================="<<endl;
  m_jmax =10000;
  for(int j=1; j<=m_jmax; j++)
    m_ypar[j-1]=HST_KKMC_NORMA->GetBinContent(j);    // xpar encoded
  //[[[[[

  m_CMSene  = m_ypar[ 1 -1];
  m_vvmax   = m_ypar[17 -1];

  m_alfinv   = m_ypar[30 -1];     // 1/alphaQED at Q^2=0
  m_gnanob   = m_ypar[31 -1];     // GeV^2 -> nanobarns

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
  //]]]]]
  char *output_file = "./kksem.output";
  long stl2 = strlen(output_file);
  long mout =16;
  kk2f_fort_open_(mout,output_file,stl2);
  kk2f_initialize_(m_ypar);
  kksem_initialize_(m_ypar);
  cout<<"================ KKsem initialization END   ==========================="<<endl;
  cout<<"==================================================================="<<endl;
  //long kdum; kk2f_getkeyfsr_(kdum); //<-- to avoid linker bug! (if no kk2f_initialize)
}

///////////////////////////////////////////////////////////////////////////////////
void KKplot::VVplot( TH1 *hstNew, long KF, char chak[5], long KeyDis, long KeyFob)
{
  long   nbin = hstNew->GetNbinsX();
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
    cout<<"KKplot::VVplot: ib= "<<ib<<"  Bin(ib) =  "<<Bin[ib]<<endl;
    hstNew->SetBinContent(ib+1, Bin[ib]);
    hstNew->SetBinError(  ib+1, 0.0);
  }
}// KKplot::VVplot



///////////////////////////////////////////////////////////////////////////////////
void KKplot::Cplot( TH1 *hstNew,
		   long KF, char chak[5], long KeyDis, long KeyFob, double vmin, double vmax)
{
  long   nbc = hstNew->GetNbinsX();
  double cmin = hstNew->GetXaxis()->GetXmin();
  double cmax = hstNew->GetXaxis()->GetXmax();
  double cBin[nbc];
  long   nbv = 1;   // only one bin in v needed
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




///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                End of Class KKplot                                        //
///////////////////////////////////////////////////////////////////////////////

