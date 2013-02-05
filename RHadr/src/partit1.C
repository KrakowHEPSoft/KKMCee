//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//       Combinatorics and polynomial for single cluster weight             //
//                                                                          //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <iostream.h>
#include <fstream.h>
#include <math.h>
#include <iomanip.h>


//////////////////////////////////////////////////////////////////////////////
double dfactorial(int k)
{
  double x=1;
  for(int n = k; n>1; n--) x=x*n;
  return x;
}
//////////////////////////////////////////////////////////////////////////////
inline double dpower(double x, int k)
{
  return( exp( k*log(x)) );
}
//////////////////////////////////////////////////////////////////////////////
void printkm1( int d, int k[], int m[])
{
  int j,mm,kk;
  for(j=1; j<d+1; j++)
    {
      for(mm=1; mm<m[j]+1; mm++)
	{
	  cout<<k[j];
	}
    }
}
//////////////////////////////////////////////////////////////////////////////
void printkm2( int d, int k[], int m[])
{
  int j,mm,kk;
  for(j=1; j<d+1; j++)
    {
      for(mm=1; mm<m[j]+1; mm++)
	{
	  for(kk=1; kk<k[j]+1; kk++ ) cout<<"O";
	  cout<<endl;
	}
    }
}
//////////////////////////////////////////////////////////////////////////////
void printkm3( int d, int k[], int m[])
{
  int j;
  cout<<"multip m[j]= "; for(j=1;j<d+1;j++) cout<<"  "<<m[j]; cout<<endl;
  cout<<"       k[j]= "; for(j=1;j<d+1;j++) cout<<"  "<<k[j]; cout<<endl;
}
//////////////////////////////////////////////////////////////////////////////
double  hfun( int k, double p, double u)
{
  if(k==1)
    return( 1 );
  else
    return( exp((k-1)*log(p)) *(p/k +u));
}
//////////////////////////////////////////////////////////////////////////////
void japww(int n, double p, double u, double &ww)
{
  int k[100];
  int m[100];
  int d,i,j;
  int s;
  double x;
  
  double silnia = dfactorial(n);
  // initializations
  for(j=0; j<n+1;j++)  k[j]=0;
  for(j=0; j<n+1;j++)  m[j]=0;
  d=1; k[1]=n; m[1]=1;
  double wwsum = hfun(n,p,u);
//  cout<< "     wwsum*n! = "<<wwsum*silnia <<"\n";
//  printkm1(d,k,m); cout <<endl;
//  printkm3(d,k,m); cout <<endl;
  int ntot=1;
  while( k[1]>1 )
    {
      ntot++;
      int suma=0;
      if(k[d]==1) {suma = m[d]; d--; }
      suma = suma + k[d]; m[d]--; s=k[d]-1;
      if(m[d]>0) {d++; }
      k[d] = s; m[d]= suma/s; s = suma % s;
      if( s!=0) { d++; k[d] = s; m[d] = 1;}
      double prod=1;
      for(i=1;i<d+1;i++) 
	{
	  x= hfun(k[i],p,u);
	  prod=prod*dpower(x,m[i])/dfactorial(m[i]);
	}
//      cout <<" prod = "<<prod;
//      cout <<"    prod*n!= "<<silnia*prod<<"  "<<"\n";
//      printkm1(d,k,m); cout <<endl;
//      printkm3(d,k,m); cout <<endl;
      wwsum = wwsum + prod;
    };
//  cout << "ntot = "<<ntot<<"\n";
  ww=wwsum*silnia;
}
//////////////////////////////////////////////////////////////////////////////
