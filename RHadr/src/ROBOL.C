///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                Class ROBOL                                                //
//                                                                           //
//    Factory of factories                                                   //
//    It contains Makers for MC generation and Analysis event per event      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "ROBOL.h"


ClassImp(ROBOL)

///////////////////////////////////////////////////////////////////////////////
void ROBOL::Initialize(long &NevTot)
{
  //////////////////////////////////////////////////////////////
  //   Initialize MC generator and analysis programs          //
  //////////////////////////////////////////////////////////////
  KoralW.ReadData( NevTot );         // Read data from disk
  KoralW.Initialize();               // Initialize generator

  //  BEfactory: Setting Initial parameters for BE model
  /////////////
  //double range   = 0.2;  // Q^2 ragne  , Epifany97!!!
  //long   ifun    = 1;    // Gausian UA1, Epifany97!!!
  //double pp      = 0.20; // Gausian UA1, Epifany97!!!
  //double radius  = 1.0;  // Gausian UA1, Epifany97!!!
  /////////////
  //double range   = 0.2;  // Q^2 range
  //long   ifun    = 2;    // Exp UA1
  //double pp      = 0.40; // Exp UA1
  //double radius  = 1.3;  // Exp UA1
  /////////////
  //double range   = 0.2;  // Q^2 ragne
  //long   ifun    = 1;    // Gausian Aleph
  //double pp      = 0.50; // Gausian Aleph
  //double radius  = 1.0;  // Gausian Aleph
  /////////////

  //  BEfactory: Setting Initial parameters for BE model
  double range   = KoralW.xpar[4061-1];
  long   ifun    = KoralW.xpar[4062-1];
  double pp      = KoralW.xpar[4063-1];
  double radius  = KoralW.xpar[4064-1];
  BEfactory.SetModel(range,ifun,pp,radius);

  //  BEfactory: Setting renormalization constants
  double lambda  = KoralW.xpar[4065-1];
  double avewt   = KoralW.xpar[4066-1];
  double lambda2 = KoralW.xpar[4067-1];
  double avewt2  = KoralW.xpar[4068-1];
  BEfactory.SetRenorm(lambda,avewt,lambda2,avewt2);

  // ROBOL seting rejection parameters
  KeyRej =  KoralW.xpar[4069-1];
  WtMax  =  KoralW.xpar[4070-1];
}
///////////////////////////////////////////////////////////////////////////////
void ROBOL::Production(long &iEvent)
{
  //////////////////////////////////////////////////////////////
  //   Generate and analyze single event                      //
  //////////////////////////////////////////////////////////////

BackToWork:
  ////////////////////////////////////////////////
  //     Generate  Raw Koralw event             //
  ////////////////////////////////////////////////
  KoralW.Generate(current_event);

  ////////////////////////////////////////////////
  //     Set up jets in final state             //
  ////////////////////////////////////////////////
  KoralW.JetDefine( current_event );

  ////////////////////////////////////////////////
  //     Test Printouts                         //
  if(iEvent < 5 ) 
    {
      current_event.Print( 1);           // Print from class KorEvent
      KoralW.LuList( current_event, 1);  // Print from Jetset
    }
  //                                            //
  ////////////////////////////////////////////////

  ////////////////////////////////////////////////
  //     calculate BoseEinstein weights etc.    //
  ////////////////////////////////////////////////
  long ntot;  // total pion multiplicity
  double wtBE, wtBE2;
  BEfactory.MakeWeight(current_event, ntot, wtBE, wtBE2);

  hst_BEwt->Fill(wtBE2);

  ////////////////////////////////////////////////
  //     REJECTION according to weight wtBE2    //
  ////////////////////////////////////////////////
  if( KeyRej == 1 ) {
    double Wt=wtBE2/WtMax;
    double rnd;
    long   len = 1;
    KoralW.VarRan(&rnd,len);
    if ( Wt <= rnd ) goto BackToWork;
    if ( Wt <= 1.0 ) Wt = 1.0;
    wtBE2 = Wt;
  };

  hst_BEwtAR->Fill(wtBE);

  ////////////////////////////////////////////////
  //     BOOK BoseEinstein weights, pair-masses //
  ////////////////////////////////////////////////
  int kf;     // KF code of pion
  kf= 211; 
  BEfactory.BookLSP(current_event,kf,ntot,wtBE,wtBE2);
  kf=-211; 
  BEfactory.BookLSP(current_event,kf,ntot,wtBE,wtBE2);
  ////kf= 111;
  ////BEfactory.BookLSP(current_event,kf,ntot,wtBE,wtBE2);
  
  ////////////////////////////////////////////////
  //          W mass reconstruction             //
  ////////////////////////////////////////////////
  J4factory.Book(current_event,ntot,wtBE,wtBE2);
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void ROBOL::Finalize()
{
  //////////////////////////////////////////////////////////////
  //   Finalize MC  run, final printouts, cleaning etc.       //
  //////////////////////////////////////////////////////////////
  KoralW.Finalize();      // Finalize generation
  double xsect = KoralW.xpar[20];
  double ersect= KoralW.xpar[21];
  cout << " KORALW: xsect  ="<<  xsect << "    ersect ="<< ersect << "\n" ;
}
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//           End of Class ROBOL                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
