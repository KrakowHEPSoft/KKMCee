#include<iostream>
{
//_________________________________________________________________
cout<<"%%%%%%%%%%%%%%%%%%%%%%%% farm-querry %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;

gROOT->Reset();
gSystem->Load("../../MCdev/.libs/libMCdev.so");

TString Dir1=gSystem->WorkingDirectory();
cout<<Dir1<<endl;

int isDir;
isDir = gSystem->ChangeDirectory(Dir1);
double TotEvent=0;

int Nfarm =250;
for(int i=1; i<Nfarm+1; i++){
  TString DIRi = Dir1;
  DIRi        +="/farm/";
  DIRi        += i;
  DIRi        += "/";
  isDir = gSystem->ChangeDirectory(DIRi);
  //cout<<"%%% isDir = "<<isDir<<endl;
  if(isDir){
    cout<<"%%% DIRi= "<< DIRi<<endl;
    TFile SemFile("semaf.root","READ","Semaphore");
    TSemaf *CheckSem = (TSemaf*)SemFile.Get("Semafor");// read semaphore from the disk
    cout<<"%%% status  = "<<CheckSem->m_flag<<"     ";
    cout<<"%%% nevtot  = "<<CheckSem->m_nevtot<<"     ";
    cout<<"%%% nevent  = "<<CheckSem->m_nevent<<endl;
    TotEvent += (double)CheckSem->m_nevent;
    //SemFile.ls();
    //SemFile.Close();
    //TFile HstFile("histo.root","READ","r.n.generator, MCgens");
    //HstFile.ls();
    //HstFile.Close();
  }//if
}//for
//*****************************************************************
//TBrowser oglad("oglad");
cout<<"%%% Total No. of Events  = " << TotEvent <<endl;
cout<<"%%%%%%%%%%%%%%%%%%%% end farm-querry %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
return 0;
}
