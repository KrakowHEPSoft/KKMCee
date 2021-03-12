///----------------------------------------------------------------
/// Sept 2008, replacement for farmqsubmit c-shell
/// Arguments: qsubCMD = name of the submit command
///            nfarm   = number of jobs
/// Note: nfarm may be bigger than no. of subdirs ./farm/*/
///----------------------------------------------------------------
#include<iostream>
int SubmFarmS(TString qsubCMD = "qsub", int nfarm = 6)
{
gROOT->Reset();
gSystem->Load("../../MCdev/.libs/libMCdev.so");
///_________________________________________________________________
TString Dir1=gSystem->WorkingDirectory();
cout<<" working directiry: " << Dir1<<endl;
int isDir;
isDir = gSystem->ChangeDirectory(Dir1);
cout<<isDir<<endl;
for(int i=1; i<nfarm+1; i++){
  TString DIRi = Dir1;
  DIRi        +="/farm/";
  DIRi        += i;
  DIRi        += "/";
  isDir = gSystem->ChangeDirectory(DIRi);
  //cout<<"%%% isDir = "<<isDir<<endl;
  if(isDir){
    cout<<"%%% DIRi= "<< DIRi<<endl;
    gSystem->Exec("ls ./*.cmd.*");
    gSystem->Exec("chmod 755 ./*.cmd.*");
    gSystem->Exec(qsubCMD+" ./*.cmd.* ");
  }
}
///*****************************************************************
///TBrowser oglad("oglad");
}
