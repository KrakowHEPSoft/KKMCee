*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//                     Pseudo-CLASS  HepEvt                                 //
*//                                                                          //
*//  Purpose:  keep and serve event in HEPEVT format                         //
*//                                                                          //
*//  Output of KK2f   is encoded in double precission /d_hepevt/             //
*//  which is double precision version of /hepevt/                           //
*//  It was necessary to rename /hepevt/                                     //
*//  because older Jetset uses REAL              version of /hepevt/                    //
*//                                                                          //
*//  We introduce luhepc->luhepcd which is                                   //
*//  Double Precision version of luhepc. It translates                       //
*//  from /hepevtd/ to (from) Jetset-common-blocks                           //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////

* ----------------------------------------------------------------------
      INTEGER nmxhep         ! maximum number of particles
      PARAMETER (nmxhep=2000)
      DOUBLE PRECISION   phep, vhep
      INTEGER nevhep, nhep, isthep, idhep, jmohep, jdahep
      COMMON /d_hepevt/
     $     nevhep,           ! serial number
     $     nhep,             ! number of particles
     $     isthep(nmxhep),   ! status code
     $     idhep(nmxhep),    ! particle ident KF
     $     jmohep(2,nmxhep), ! parent particles
     $     jdahep(2,nmxhep), ! childreen particles
     $     phep(5,nmxhep),   ! four-momentum, mass [GeV]
     $     vhep(4,nmxhep)    ! vertex [mm]
      SAVE  /d_hepevt/
* ----------------------------------------------------------------------
      LOGICAL qedrad
      COMMON /phoqed/ 
     $     qedrad(nmxhep)    ! Photos flag
      SAVE   /phoqed/
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//                      End of CLASS  HepEvt                                //
*//////////////////////////////////////////////////////////////////////////////
