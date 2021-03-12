{
gROOT->Reset();
gSystem->Load("../../MCdev/.libs/libMCdev.so");
//_________________________________________________________________

TString Dir1=gSystem->WorkingDirectory();
cout<<Dir1<<endl;

int isDir;
isDir = gSystem->ChangeDirectory(Dir1);
cout<<isDir<<endl;

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
    TFile SemFile("semaf.root","UPDATE","Semaphore");
    SemFile.cd();
    TSemaf *CheckSem = (TSemaf*)SemFile.Get("Semafor");    // read semaphore from the disk
    cout<<"%%% Stop.C>>  Old flag= "<<CheckSem->m_flag<<endl;
    if( CheckSem->m_flag != "START") CheckSem->m_flag = "STOP";                     // modify flag
    CheckSem->Write("Semafor",TObject::kOverwrite);// write modified semaphore
    SemFile.Write("",TObject::kOverwrite);
    //SemFile.ls();
    //SemFile.ShowStreamerInfo();
    SemFile.Close();
  }
}

//*****************************************************************
//TBrowser oglad("oglad");

}
