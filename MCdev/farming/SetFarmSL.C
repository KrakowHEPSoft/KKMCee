///----------------------------------------------------------------
/// This is upgraded version of sept 2008
///----------------------------------------------------------------
#include<iostream>
int SetFarmSL(TString dname = "WORK", int nfarm = 6)
{
///_________________________________________________________________
  cout<<"%%% ================== SetFarm ================== %%%%"<<endl;
  gROOT->Reset();
  gSystem->Load("../../MCdev/.libs/libMCdev.so");
  TString Dir1=gSystem->WorkingDirectory();
  cout<<"%%% working directory= " <<Dir1<<endl;
//   cout<<"%%% working directory= " <<Dir1<<endl;
//   TString Dir1=gSystem->WorkingDirectory();
//   cout<<"%%% working directory= " <<Dir1<<endl;
//   TString path_to_lib = Dirl+"../../lib/libProj.so";
//   cout << "path to libProj.so: " << path_to_lib << endl;
//   gSystem->Load(Form(path_to_lib));
///
  int Nfarm     = nfarm; // no. of farm nodes
  TString DNAME = dname; // test name = work dir. name
///
  cout<<"%%% Nfarm= "<< Nfarm <<" DNAME= "<< DNAME <<endl;
///
  
///
  TRandom Rngen;
///
  int isDir = gSystem->ChangeDirectory("farm");
  if(isDir)
    cout<<"%%% SetFarm STOPPED !!! Remoove ./farm"<<endl;
  else{
    gSystem->MakeDirectory("farm");
    gSystem->ChangeDirectory("farm");
    // loop over farm nodes
    for(int i=1; i<Nfarm+1; i++){
      TString DIRnew ="./";
      DIRnew      += i;
      DIRnew      += "/";
      //cout<<"%%% DIRnew= "<< DIRnew<<endl;
      gSystem->MakeDirectory(DIRnew);   // make ./farm/i/
      int isDIRnew = gSystem->ChangeDirectory(DIRnew);
      if(isDIRnew ){
         TString DIRi=gSystem->WorkingDirectory();// local dir. for batch script
         cout<<"%%%%%%%%% New DIRi= "<< DIRi <<endl;
///--------------WRITE BATCH SCRIPT ------------------
         TString ScrName=DNAME;
         ScrName += ".cmd."; ScrName += i;
         ofstream   OutFile(ScrName,ios::out);
///----------------------------
         OutFile<< "#!/bin/bash -l   " <<endl;
         OutFile<< "#SBATCH -J "<<DNAME<<"-"<<i <<endl;
         OutFile<< "#SBATCH -o slurm.%N.%j.out # STDOUT " <<endl;
         OutFile<< "#SBATCH -e slurm.%N.%j.err # STDERR" <<endl;
///----------------------------
         OutFile<< "#######################################" <<endl;
         OutFile<< "cd "<<DIRi<<endl;
         OutFile<< "# ========= EXECUTION ========= " <<endl;
         OutFile<< "time ../"<<DNAME<<".exe" <<endl;
         OutFile<< "exit 0          " <<endl;
         OutFile.close();
         gSystem->Exec("chmod 755 ./*.cmd.*");
///------- COPY MCGEN.ROOT AND RESET RN. SEED ---------------
         gSystem->Exec("cp ../../pro.input  ./pro.input");
         gSystem->Exec("cp ../../histo.root ./histo.root");
         gSystem->Exec("cp ../../semaf.root ./semaf.root");
         gSystem->Exec("cp ../../mcgen.root ./mcgen.root");
         gSystem->Exec("ls -altr");
         TFile GenFile("./mcgen.root","UPDATE","r.n.generator, MCgens");
        //GenFile.ls();
         long iniseed = Rngen.Integer(1e8)+1;
         //iniseed=54217137;   // <--this is for tests only
         cout<<"%%% New RN gen. seed = "<< iniseed <<endl;
         TRandom  *RN_gen   = (TRandom*)GenFile.Get("RN_gen");  // read r.n. generator
         RN_gen->SetSeed(iniseed);
         RN_gen->Write("RN_gen",TObject::kOverwrite);   // RN generator status
         GenFile.Write();
         //GenFile2.ls();
         GenFile.Close(); //!!!!
///---------Redefine/increase group statistics nevgrp for farming
         TFile SemFile("./semaf.root","UPDATE","Semafor");
         SemFile.cd();
         TSemaf *CheckSem = (TSemaf*)SemFile.Get("Semafor");       //! read semaphore from the disk
         CheckSem->m_nevgrp = (Nfarm/8.0) *CheckSem->m_nevgrp; //! increase group no. of MC evts.
         cout<<"%%% New m_nevgrp    =  "<< CheckSem->m_nevgrp <<endl;
         CheckSem->Write("Semafor",TObject::kOverwrite);   //! write modified semaphore
         SemFile.Write("",TObject::kOverwrite);
         SemFile.Close();
     }else{
         cout<<"+++++ something wrong with creating DIRnew= "<<DIRnew<<endl;
      }/// if isDir
      gSystem->ChangeDirectory("../"); // IMPORTANT!!!!
    }/// for i, loop over farm nodes
  }/// if isDir
  cout<<"%%% ============== End SetFarm ================== %%%%"<<endl;
  //*****************************************************************
  //TBrowser oglad("oglad");
  return 0;
}///-------------------------------------------------------------------------

