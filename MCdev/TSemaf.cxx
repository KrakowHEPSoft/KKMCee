///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                CLASS TSemaf                                               //
//                                                                           //
//    Utility for managing main MC production loop.                          //
//    It records semaphore member "m_flag" which tells running               //
//    program to continue or stop.                                           //
//    It alse encapsulates original random number seed and records           //
//    actual total number of generated events.                               //
//    It is a translation of similar tool from fortran (bhlumi)              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
#include "TSemaf.h"

ClassImp(TSemaf);

///////////////////////////////////////////////////////////////////////////////
TSemaf::TSemaf()
{
  // default Constructor of TSemaf, for streamer
  
  //cout<<"====> TSemaf::TSemaf DEFAULT Constructor "<<endl;
} // DEFAULT Constructor


///////////////////////////////////////////////////////////////////////////////
TSemaf::TSemaf(TString flag)
{
  // user Constructor of TSemaf, for streamer

  //cout<<"====> TSemaf::TSemaf USER Constructor "<<endl;
  m_flag     = flag;
  m_nevtot   = 2e12;
  m_nevgrp   = 1e5;
  m_nevent   = 0;
}// USER Constructor

///////////////////////////////////////////////////////////////////////////////
TSemaf::TSemaf(TString flag, double nevtot)
{
  // user Constructor of TSemaf, for streamer

  cout<<"====> TSemaf::TSemaf USER Constructor "<<endl;
  m_flag     = flag;
  m_nevtot   = nevtot;
  m_nevgrp   = 1e5;
  m_nevent   = 0;
}// USER Constructor

///////////////////////////////////////////////////////////////////////////////
TSemaf::TSemaf(TString flag, double nevtot, long nevgrp)
{
  // user Constructor of TSemaf, for streamer

  cout<<"====> TSemaf::TSemaf USER Constructor "<<endl;
  m_flag     = flag;
  m_nevtot   = nevtot;
  m_nevgrp   = nevgrp;
  m_nevent   = 0;
}// USER Constructor

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//           End of CLASS TSemaf                                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
