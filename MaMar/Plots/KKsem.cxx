///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                Class KKsem                                                //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
#include "KKsem.h"

///////////////////////////////////////////////////////////////////////////////////
void KKsem::Initialize(TFile &DiskFileA){
  //------------------------------------------------------------------------
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  cout<<"==================================================================="<<endl;
  cout<<"================ KKsem::initialize begin ==========================="<<endl;
  const int jmax =10000;
  double ypar[jmax];
  for(int j=1; j<=jmax; j++)
    ypar[j-1]=HST_KKMC_NORMA->GetBinContent(j);    // xpar encoded
  //[[[[[
  //for(int j=0;j<30;j++)
  //  cout<<j<<"   "<<ypar[j]<<endl;
  //]]]]]
  char *output_file = "./kksem.output";
  long stl2 = strlen(output_file);
  long mout =16;
  kk2f_fort_open_(mout,output_file,stl2);
  kk2f_initialize_(ypar);
  kksem_initialize_(ypar);
  cout<<"================ KKsem_initialize END   ==========================="<<endl;
  cout<<"==================================================================="<<endl;
  //long kdum; kk2f_getkeyfsr_(kdum); //<-- to avoid linker bug! (if no kk2f_initialize)
}

///////////////////////////////////////////////////////////////////////////////////
void KKsem::VVplot( TH1 *hstNew, long KF, char chak[5], long KeyDis, long KeyFob)
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
    cout<<"KKsem::VVplot: ib= "<<ib<<"  Bin(ib) =  "<<Bin[ib]<<endl;
    hstNew->SetBinContent(ib+1, Bin[ib]);
    hstNew->SetBinError(  ib+1, 0.0);
  }
}// KKsem::VVplot


///////////////////////////////////////////////////////////////////////////////////
void KKsem::Cplot( TH1 *hstNew, 
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
    cout<<"KKsem::Cplot: ib= "<<ib<<" cBin(ib) =  "<<Dsig<<endl;
    hstNew->SetBinContent(ib+1, Dsig);
    hstNew->SetBinError(  ib+1, 0.0);
  }
  kksem_setcrange_(-1.0,1.0); // back to normal
}// KKsem::VVplot

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                End of Class KKsem                                         //
//////////////////////////////////////////////////////////////////////////////

