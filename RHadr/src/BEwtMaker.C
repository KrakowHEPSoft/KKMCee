///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                CLASS BEwtMaker                                            //
//                                                                           //
//    calculation of Bose-Einstein weight for a given event from Koralw      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "BEwtMaker.h"

#include <stdlib.h>
#include <math.h>


ClassImp(BEwtMaker)
/////////////////////////////////////////////////
BEwtMaker::BEwtMaker()
{
  //
  // Constructor of BEwtMaker
  //
  BEoutfile.open("BEwtMaker.output");
  //
  BEievent     = 0;
  BEprintlevel = 1;
  BEprintlast  = 20;
  //
  hst_nclu  = new TH1F("hst_nclu" ,"cluster dist",20,0,20);
  hst_unQ2  = new TH1F("hst_unQ2" ,"Q2 2-pi"     ,100,0,1);
  hst_wtQ2  = new TH1F("hst_wtQ2" ,"Q2 2-pi wted",100,0,1);
  hst_wt2Q2 = new TH1F("hst_wt2Q2","Q2 2-pi wted",100,0,1);
  hst_unQ3  = new TH1F("hst_unQ3" ,"Q2 3-pi"     ,100,0,1);
  hst_wtQ3  = new TH1F("hst_wtQ3" ,"Q2 3-pi wted",100,0,1);
  hst_wt2Q3 = new TH1F("hst_wt2Q3","Q2 3-pi wted",100,0,1);
  //
  ctuple2_counter  =0;
  ctuple2_max      =10000;
  ctuple2 = new TNtuple("ctuple2","Q2 ","Q:np:wt");
  ctuple3_counter  =0;
  ctuple3_max      =10000;
  ctuple3 = new TNtuple("ctuple3","Q3 ","Q:np:wt");
  //
}
//////////////////////////////////////////////////////////////////////////////
BEwtMaker::SetModel(double range, long FuncType, double pp, double radius)
{
// Initialize weight parameters
  BErange   = range;
  BEFuncType= FuncType;
  BEpp      = pp;
  BEradius  = radius;
}
//////////////////////////////////////////////////////////////////////////////
BEwtMaker::SetRenorm(double lambda,  double avewt, 
		     double lambda2, double avewt2)
{
// Initialize weight parameters
  BElambda  = lambda;
  BEavewt   = avewt;
  BElambda2 = lambda2;
  BEavewt2  = avewt2;
}
//////////////////////////////////////////////////////////////////////////////
double BEwtMaker::CorFun1(double Q2, int mult, double Rf, double pp)
{
  // Correlation function EXPONENTIAL UA1 fit
  double R = Rf/0.2;                   // conversion to GeV^-1
  R = R*sqrt((mult*(mult-1))/2.0);     // UA1 fit, R increasing with n
  double Q2p = Q2*2/(mult*(mult-1));
  double Qp  = sqrt(Q2p);
  double x = Qp*R;
  double zexp= 0; 
  if( x < 100.0 ) zexp = exp(-x);
  return zexp;
}
//////////////////////////////////////////////////////////////////////////////
double BEwtMaker::CorFun2(double Q2, int mult, double Rf, double pp)
{
  // Correlation function GAUSSIAN UA1 fit
  double R = Rf/0.2;                   // conversion to GeV^-1
  double Q2p = Q2*2/(mult*(mult-1));
  double x = Q2p*(R*R);
  double zexp= 0; 
  if( x < 100.0 ) zexp = exp(-x);
  return zexp;
}
//////////////////////////////////////////////////////////////////////////////
void BEwtMaker::WTcluster(int mult, PartLund *first, double &wt)
{ 
  //*************************************************************************//
  // Operates on a SINGLE cluster *first which is still in form of a 
  // LIST starting at the element *first,
  // Calculates weight wt.
  //*************************************************************************//
  int i,j,k,l,m,n;
  int ntot;
  float x;
  PartLund *clmat   = new PartLund[mult];
  float A[100][100];
  int link[100][100];
  int noli[100];
  float suli[100];
  
  // Copy etire cluster into matrix clmat
  PartLund *actual; j=0;
  for (actual = first; actual != NULL; actual = actual->previous)
    {
      clmat[j] = *actual;
      j++;
    };

  // real    matrix A containing the value of relative Q2
  for(j=0; j<mult;j++) A[j][j]=0.0;  // diagonal elements
  double Q2jk;
  double Q2 = 0;
  for(j=0; j<mult;j++)               // off-diagonal elements
      for(k=0; k<j;k++)
	{
	  double Q2jk= Q2pair( &(clmat[j]),&(clmat[k]));
	  A[j][k]=Q2jk; A[k][j]=Q2jk;
	  Q2=Q2 + Q2jk;
	}

  wt=1;
  if(mult>1)
    {
      ////double pp;   // coherence parameter pp=p=lambda
      ////double R ;   // radius in fermi <===== !!!!!
      ////int ifun=1;  // Correlation function type
      double zexp;
      int     ifun = BEFuncType;   // Correlation function type
      double  pp   = BEpp;         // coherence parameter pp=p=lambda
      double  R    = BEradius;     // radius in fermi <===== !!!!!
      if( ifun == 2 )
	{
	  // exponential
	  ////pp=0.40; //  UA1 fit
	  ////R=1.3;   //  UA1 fit
	  zexp =CorFun1(Q2, mult, R, pp);
	}
      else
	{
	  // gaussian
	  ////pp=0.20;// UA1 fit
	  ////R=1.0;  // UA1 fit
	  ////pp=0.40;  // Aleph fit
	  ////R=1.0;    // Aleph fit
	  zexp =CorFun2(Q2, mult, R, pp);
	}
      // parameters of Kacper
      double pz=pp*zexp;
      double uz=1-pp;
      double wwjap;
    //**********************
      japww(mult,pz,uz,wwjap);
    //**********************
      wt=wwjap;
      if(wt>100.0)
	{
	  double Qp = sqrt(Q2*2/(mult*(mult-1)));
	  cout<<"mult= "<<mult<<" wt= "<<wt;
	  cout<<" pz="<<pz<<" pp="<<pp<<" Qp= "<<Qp ;
	  cout<<" zexp= "<<zexp<<endl;
	}
    }

  // IMPORTANT!!! detach dynamical storage!!!
  // IMPORTANT!!! detach dynamical storage!!!
  delete [] clmat;
}
//////////////////////////////////////////////////////////////////////////////
void BEwtMaker::BEchain2(KorEvent &event, int kfP, int &nP, double &wt)
{
//==========================================================//
// Forms chain of MANY CLUSTERS of particles type kfP
//==========================================================//
  int iPrint =0;  // Printing flag for debugs
  int j,k;
  double Qp,QpMin;
  PartLund *head, *actual;
  PartLund *first, *second;
  PartLund *clhead, *cltail, *next;
  PartLund *clust[500];
  int      nclust[500];
  int mulclus,np;
  int nptot=0;
  double Range=0.20; // maximum sqrt(Q2) range for a single pair to qualify

  nP=0; // default value 
  wt=1; // default value
  
  if(iPrint) cout << "||||||||||||||| BEchain2 Begin ||||||||||||||||" << endl;
  head=NULL;  
  for (j=0; j< event.m_npart ;j++)
    if( event.m_part[j].m_kflavor ==  kfP ) 
      {
	event.m_part[j].previous = head;
	head   = &event.m_part[j];
	if(iPrint) {(*head).m_pmom.print(); cout<<endl;}
	nP++;
      };

  if(iPrint) event.Print( 1);
  if(iPrint) {cout << "ALL: "; (*head).ListPrint(); cout<<endl;}
  if(iPrint) cout << "==============================================" << endl;

  if(head == NULL) goto End;

  // I this part we form clusters of particles as LISTS which
  // begin at pointer *clust[i], i< mulclus
  // multiplicity of i-th cluster is stored in nclust[i]

  for(mulclus = 0; mulclus< 500; mulclus++)
    {
      np=1;
      cltail = head;          // last element in cluster segment
      for (actual = head; actual->previous != NULL; )
	{
	  if(iPrint) {(*actual).print(0); cout << " <-actual " << endl;}
	  next = actual->previous;
	  QpMin = 1e20;
	  for(first =head; ; first = first->previous)
	    {
	      Qp= sqrt(Q2pair(first, next));
	      if(Qp < QpMin) QpMin=Qp;
    if(iPrint){(*first).print(0);cout<<" <-first,QpMin="<<QpMin<<endl;}
	      if( first == cltail ) goto EE;
	    };
	EE:
	  if(QpMin < Range)
	    {
	      if(iPrint) {(*next).print(0);cout<<" <-detach "<< endl;}
	      actual->previous = actual->previous->previous; // detach
	      next->previous   = head;                       // attach
	      head= next;                                    // attach
	      np++;
	    }
	  else
	    {
	      actual= actual->previous;  // move further on
	    };
	};
      clhead = head;              // head of cluster
      head   = cltail->previous;  // new head of the remainder
      cltail->previous = NULL;    // grounding end of cluster
      clust[mulclus] = clhead;
      nclust[mulclus] = np;
      //
      int mlt2=np*(np-1)/2;
      hst_nclu->Fill(np+0.001,mlt2);     // cluster multiplicity
      //
      nptot = nptot+np;
     
      if(iPrint)
	{
	  cout<<"-----------> " << mulclus <<endl;
	  cout<<"np= "<<np << ": "; (*clhead).ListPrint(); cout<<endl;
	  (*head).ListPrint(); cout << endl;
	}

      if(head == NULL ) goto EEX;
    };
EEX:

  if(iPrint)
    {
      if( kfP ==  211) cout << "+++ "; 
      if( kfP == -211) cout << "--- "; 
      if( kfP ==  111) cout << "000 "; 
      cout << "|| nptot= "<<nptot << ":  ";
      cout << "mulclus= "<< mulclus << ":  ";
      for (j=0; j<mulclus; j++) cout << nclust[j] << " ";
      cout << endl;
    }

  double wtclu;
  for (j=0; j<mulclus; j++)
    if( nclust[j] > 1)
      {
	WTcluster(nclust[j], clust[j], wtclu);
	wt=wt*wtclu;
      }

End:
  if(iPrint) cout << "||||||||||||||| BEchain2 End ||||||||||||||||" << endl;
};  // END of BEchain2
//
//////////////////////////////////////////////////////////////////////////////
void BEwtMaker::MakeWeight(KorEvent &event, 
			   long &ntot, double &wtot, double &wtot2)
{
  BEievent++;
  //
  // Makes Bose-Einstein weight for a given event
  //
  int kf_PiPlus  =  211;
  int kf_PiMinus = -211;
  int kf_PiZero  =  111;

  int   nPiPlus, nPiMinus, nPiZero;
  double wtPlus,  wtMinus,  wtZero;
  BEchain2(event, kf_PiPlus,  nPiPlus,  wtPlus);
  BEchain2(event, kf_PiMinus, nPiMinus, wtMinus);
  BEchain2(event, kf_PiZero,  nPiZero,  wtZero);


  ntot= nPiPlus+nPiMinus+nPiZero;
  wtot2  = wtPlus*wtMinus*wtZero*exp(-ntot *log(BElambda2))/BEavewt2;

  double Avewt  = exp(log(BEavewt)/3); // qubic root
  wtPlus  = wtPlus *exp(-nPiPlus *log(BElambda))/Avewt;
  wtMinus = wtMinus*exp(-nPiMinus*log(BElambda))/Avewt;
  wtZero  = wtZero *exp(-nPiZero *log(BElambda))/Avewt;
  
  wtot = wtPlus*wtMinus*wtZero;

  if( (wtot>1000.0) || (BEievent <= BEprintlast ) )
    {
      cout << "(" << BEievent << "),  ";
      cout << " n = " <<" "<< event.m_npart  << ",  ";
      cout << nPiPlus <<" "<< nPiMinus <<" "<< nPiZero  << "  ";
      cout << " wtPlus  = " << wtPlus  << "   ";
      cout << " wtMinus = " << wtMinus << "   ";
      cout << " wtZero  = " << wtZero  << "   ";
      cout << " wtot = "    << wtot << endl;
      //
      BEoutfile << "(" << BEievent << "),  ";
      BEoutfile << " n = " <<" "<< event.m_npart  << ",  ";
      BEoutfile << nPiPlus <<" "<< nPiMinus <<" "<< nPiZero  << "  ";
      BEoutfile << " wtPlus  = " << wtPlus  << "   ";
      BEoutfile << " wtMinus = " << wtMinus << "   ";
      BEoutfile << " wtZero  = " << wtZero  << "   ";
      BEoutfile << " wtot = "    << wtot << endl;
    }
}
//////////////////////////////////////////////////////////////////////////////
void BEwtMaker::BookLSP(KorEvent &event, 
			int kfP, int ntot, double wt, double wt2)
{
  //
  //   booking histos/ntuples on Like-Sign-Pions and triplets
  //
  int j,k,l;
  PartLund *pion[500];

  //cout<<"***********************************************"<<endl; 
  int nP=0;
  for (j=0; j< event.m_npart ;j++)
    if( event.m_part[j].m_kflavor ==  kfP ) 
      {
	pion[nP] = &event.m_part[j];
	nP++;
      };

  //cout<<"----------------------------------------------"<<endl; 
  //cout<<" nP= "<<nP<<" wt= "<<wt<<endl;

  ///////////////////////////////
  //    like-sign-pion PAIRS   //
  ///////////////////////////////
  for(j=0; j<nP;j++)
      for(k=0; k<j;k++)
	{
	  double Q2jk= Q2pair(pion[j],pion[k]);
	  hst_unQ2->Fill(Q2jk);
	  hst_wtQ2->Fill(Q2jk,wt);
	  hst_wt2Q2->Fill(Q2jk,wt2);
	  if(ctuple2_counter < ctuple2_max)
	    if(Q2jk<1.0) {ctuple2->Fill(Q2jk,ntot,wt); ctuple2_counter++;}
	  //////////////////////////////////
	  //    like-sign-pion TRIPLETS   //
	  //////////////////////////////////
	  for(l=0; l<k;l++)
	    {
	      double Q2jl= Q2pair(pion[j],pion[l]);
	      double Q2kl= Q2pair(pion[k],pion[l]);
	      double Q3jkl = Q2jk+Q2jl+Q2kl;
	      hst_unQ3->Fill(Q3jkl);
	      hst_wtQ3->Fill(Q3jkl,wt);
	      hst_wt2Q3->Fill(Q3jkl,wt2);
	      if(ctuple3_counter < ctuple3_max)
		if(Q3jkl<1.0){ctuple3->Fill(Q3jkl,ntot,wt); ctuple3_counter++;}
	    }
	}

}
///////////////////////////////////////////////////////////////////////////////
double BEwtMaker::Q2pair(PartLund *first, PartLund *second)
{
//==========================================================================//
//                   Miscelaneous function for BE effect                    //
// calculates positive Q2 = -(q1-q2)^2 for a pair of particles              //
//==========================================================================//
  VLorenz  Q;
  double   Q2;

  Q = (*first).m_pmom - (*second).m_pmom;
  Q2 = -(Q*Q);
  return Q2;
}
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//           End of CLASS BEwtMaker                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
