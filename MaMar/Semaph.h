#ifndef Semaph_H
#define Semaph_H
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
#include "TROOT.h"


//class Semaph : public TNamed {
class Semaph{
public:
//
  ifstream   InfileSema;  //! Input semaphore file
  ifstream   InfileSeed;  //! Input seed file
  ofstream   OufileSema;  //! Output semaphore file
  long       ijklin;      // for ranmar
  long       ntotin;      // for ranmar
  long       ntot2n;      // for ranmar
public:
  Semaph();
  ~Semaph() {;}
// dummy copy constructor
private:
  Semaph(const Semaph &org) { }
public:
void Semaph::Initialize(TString &Semaphore);
void Semaph::ReadStatus(TString &Semaphore);
////////////////////////////////////////////////////////////////////////////
//                      ClassDef(Semaph,1)
};
////////////////////////////////////////////////////////////////////////////
#endif
