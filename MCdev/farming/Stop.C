{
gROOT->Reset();
gSystem->Load("../../lib/libMCdev.so");
//_________________________________________________________________
TFile SemFile("./semaf.root","UPDATE","Semafor");
SemFile.cd();
CheckSem = (TSemaf*)SemFile.Get("Semafor");    // read semaphore from the disk
cout<<"%%% Stop.C>>  Old flag= "<<CheckSem->m_flag<<endl;
CheckSem->m_flag = "STOP";                     // modify flag
CheckSem->Write("Semafor",TObject::kOverwrite);// write modified semaphore
SemFile.Write("",TObject::kOverwrite);
//SemFile.ls();
//SemFile.ShowStreamerInfo();
SemFile.Close();

//*****************************************************************
TBrowser oglad("oglad");

}
