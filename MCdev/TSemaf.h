#ifndef Semaph_H
#define Semaph_H
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                CLASS TSemaf                                               //
//                                                                           //
//    Utility for managing main MC production loop.                          //
//    It reads from semaphore file flag "semaphor" which tells running       //
//    program to continue or stop.                                           //
//    It alse initializes external random number generator from the seed     //
//    read from disk file and dumps position of the random number generator  //
//    at the end of the run.                                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
#include <iostream>
using namespace std;

#include "TROOT.h"


class TSemaf : public TObject {
public:
  TString    m_flag;        // flag = "CONTINUE","START","STOP"
  double     m_nevtot;      // no. of MC events to be generated
  double     m_nevent;      // no. of MC events to be generated
  long       m_nevgrp;      // no. of groups of MC events
public:
  TSemaf();                 // only for streamer
  TSemaf(TString flag);
  TSemaf(TString flag, double nevtot);
  TSemaf(TString flag, double nevtot, long nevgrp);
  ~TSemaf() {;}
//private:
//  TSemaf(const TSemaf &org) { } // dummy copy constructor
////////////////////////////////////////////////////////////////////////////
                      ClassDef(TSemaf,1);
};
////////////////////////////////////////////////////////////////////////////
#endif
