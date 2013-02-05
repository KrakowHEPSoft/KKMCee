///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                CLASS Semaph                                               //
//                                                                           //
//    Utility for managing main MC production loop.                          //
//    It reads from semaphore file flag "semaphor" which tells running       //
//    program to continue or stop.                                           //
//    It alse initializes random number generator RanMar from the seed       //
//    read from disk file and dumps position of the random number generator  //
//    at the end of the run.                                                 //
//    It is a direct translation of similar tool from fortran (bhlumi)       //
//    and will be hopefully replaced with something mor modern.              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <iostream.h>
#include <fstream.h>
#include <iomanip.h>

#include "TROOT.h"

#include "Semaph.h"


///////////////////////////////////////////////////////////////////////////////
//                           Fortran connections                             //
//       Contact with random number generator from Koralw generator          //
///////////////////////////////////////////////////////////////////////////////
//extern "C" void marini_(long &ijklin, long &ntotin, long &ntot2n);
//extern "C" void marout_(long &ijklin, long &ntotin, long &ntot2n);
extern "C" void pseumar_initialize_(long &ijklin, long &ntotin, long &ntot2n);
extern "C" void pseumar_out_(long &ijklin, long &ntotin, long &ntot2n);
///////////////////////////////////////////////////////////////////////////////



//ClassImp(Semaph)

Semaph::Semaph()
{
  //
  // Constructor of Semaph
  //
  ijklin = 46785717; // test (default is 54217137)
  ntotin =        0; // default
  ntot2n =        0; // default
}
///////////////////////////////////////////////////////////////////////////////
void Semaph::Initialize(TString &Semaphore)
{
  //////////////////////////////////////////////////////////////
  //   Initialization before generation of the first event    //
  //////////////////////////////////////////////////////////////
  char txtdum [200];
#define SKIPLINE(F) F.getline(txtdum,200,'\n')

  InfileSema.open("./semaphore");  //read directive from semaphore file
  InfileSema>>setw(4)>>Semaphore;

  if( Semaphore == "STAR" )
    {
      ////////////////////////////////////////////////
      //   start from scratch, read seed from disk  //
      ////////////////////////////////////////////////
      InfileSema.close();
      //
      InfileSeed.open("./iniseed");
      InfileSeed >> ijklin; SKIPLINE(InfileSeed);
      InfileSeed >> ntotin; SKIPLINE(InfileSeed);
      InfileSeed >> ntot2n; SKIPLINE(InfileSeed);
      InfileSeed.close();
      //initialization of the random number generator
      //marini_(ijklin,ntotin,ntot2n);
      pseumar_initialize_(ijklin,ntotin,ntot2n);
    }
  if( Semaphore == "CONT" )
    {
      ////////////////////////////////////////////////
      // read seed from semaphore file and continue //
      ////////////////////////////////////////////////
      InfileSema >> ijklin; SKIPLINE(InfileSema);
      InfileSema >> ntotin; SKIPLINE(InfileSema);
      InfileSema >> ntot2n; SKIPLINE(InfileSema);
      InfileSema.close();
      //initialization of the random number generator
      //marini_(ijklin,ntotin,ntot2n);
      pseumar_initialize_(ijklin,ntotin,ntot2n);
    }
  if( Semaphore == "STOP" )
    {
      InfileSema.close();
      ////////////////////////////////////////////////
      // Stop flag encountered for initialization   //
      ////////////////////////////////////////////////
      cout<< " ++++++++++++++++++++++++++++++++++++++ \nl";
      cout<< "        WRONG SEMAPHORE                 \nl";
      cout<< " ++++++++++++++++++++++++++++++++++++++ \nl";
    }
}
//////////////////////////////////////////////////////////////////////////////
void Semaph::ReadStatus(TString &Semaphore)
{
  ////////////////////////////////////////////////////////////////////////
  // Reads semaphore file during the event generation in the main loop  //
  // This routine reads closes and overwrites semaphore file            //
  // Decision whether to stop or continue done in main prog.            //
  ////////////////////////////////////////////////////////////////////////
  InfileSema.open("./semaphore");  //read directive from semaphore file
  InfileSema>>setw(4)>>Semaphore;
  InfileSema.close();

  //marout_(ijklin,ntotin,ntot2n);
  pseumar_out_(ijklin,ntotin,ntot2n);

  OufileSema.open("./semaphore");
  OufileSema<<"CONT"<<"\n";
  OufileSema<<ijklin<<"\n";
  OufileSema<<ntotin<<"\n";
  OufileSema<<ntot2n<<"\n";
  OufileSema.close();
}
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//           End of CLASS Semaph                                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
