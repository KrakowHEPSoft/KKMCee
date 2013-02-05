
*=======================================================================
*=======================================================================
* attics  attics attics attics attics attics attics attics attics attics
* attics  attics attics attics attics attics attics attics attics attics
*=======================================================================
*=======================================================================

      FUNCTION Bornv_old(svari,costhe)
*     ***********************************
* this routine provides Born differential cross section
* a version without complex*16
*     ***********************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      COMMON / inout  / ninp,nout
      COMMON / weking / ene,amaz,gammz,amel,amfin,xk0,sinw2,ide,idf
      COMMON  /xpar/  xpar(1000)
      SAVE
*
      
      KeyZet = xpar(501)
      IF(abs(costhe) .GT. 1d0) WRITE(nout,*) ' Bornv: costhe=',costhe
      qe= -1d0
      qf= -1d0
      aa= 4d0*sqrt(sinw2*(1d0-sinw2))
      ve= (-1d0+4*sinw2)/aa
      ae= 1d0/aa
      IF(keyzet .LE. 0) THEN
        ve=0d0
        ae=0d0
      ENDIF
      vf= (-1d0+4*sinw2)/aa
      af= 1d0/aa
      s = svari
      chi2 =         s**2/((s-amaz**2)**2+(gammz*amaz)**2)
      rechi=(s-amaz**2)*s/((s-amaz**2)**2+(gammz*amaz)**2)
      xe= ve**2 +ae**2
      xf= vf**2 +af**2
      ye= 2*ve*ae
      yf= 2*vf*af
      ff0= qe**2*qf**2 +2*rechi*qe*qf*ve*vf +chi2*xe*xf
      ff1=             +2*rechi*qe*qf*ae*af +chi2*ye*yf
      Born    = (1d0+ costhe**2)*ff0 +2d0*costhe*ff1
*************
* this is a bit crude method of introducing threshold behaviour
      IF(    svari .LE.  4d0*amfin**2) THEN
        thresh=0d0
      ELSEIF(svari .LE. 16d0*amfin**2) THEN
        amx2=4d0*amfin**2/svari
        thresh=sqrt(1d0-amx2)*(1d0+amx2/2d0)
      ELSE
        thresh=1d0
      ENDIF
      Bornv_old= Born*thresh
      END


      SUBROUTINE xxrho3(vv,rho)
*     ************************
* !!!!!! obsolete !!!!!!!!
* muon mass distr. initial+final st.bremss.
* gauss method, change of variables with help of chbin3
*     ************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      EXTERNAL xdis3
      COMMON / cxdis3 / vvv
      SAVE   / cxdis3 /
*
      vvv =  vv
      prec= -0.5d-3
*c      prec= -1d-3
      xa= 0d0
      xb= 1d0
      CALL gausjd(xdis3,xa,xb,prec,result)
      rho=result
      END
      FUNCTION xdis3(r)
*     ******************
* !!!!!! obsolete !!!!!!!!
* integrand for yrho initial + final state bremss. convolution
*     ************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER( fleps =1d-30)
      PARAMETER(pi=3.1415926535897932d0)
      PARAMETER(gnanob = 389.385d-30*1.d33  )
      COMMON / inout  / ninp,nout
      COMMON / weking / ene,amaz,gammz,amel,amfin,xk0,sinw2,ide,idf
      COMMON  /xpar/  xpar(1000)
* indirect input from CALLing PROGRAMs
* keydis=keydis_vvrho + 1000*keydis_uurho
      COMMON / keydis / keydis
      COMMON / cxdis3 / vvv
      DIMENSION softi(5)
      SAVE
*
      KeyZet = xpar(501)
      alfinv = xpar(30)
*
      alf1 = 1d0/alfinv/pi
      keydi1= MOD(keydis,1000000)/1000
      keydi2= MOD(keydis,1000)
      svar   =  4d0*ene**2
      sig0nb =  4d0*pi/(alfinv**2*3d0*svar)  *gnanob
      beti   =  2*alf1*(dlog(svar/amel**2) -1)
      betf2  =  2*alf1*(dlog(svar/amfin**2)-1)
* importance sampling PARAMETERs for chbin3
      alf =  beti
      bet =  betf2
      rmin = fleps**alf
      rmax = 1d0 - fleps**bet
      IF(rmax .LE. rmin) GOTO 900
* translate r into t
      x = max(r,rmin)
      x = min(r,rmax)
      CALL chbin3(x,alf,bet,t,t1,rjac)
* kinematics
      v = vvv*t
      u = vvv*t1/(1d0-vvv*t)
      svar1  =  svar * (1d0  -v)
* components
      Bornb  =  Borny(svar1)/(1d0-v) *sig0nb
* constant  x-section for tests
      IF(keyzet .EQ. -2) Bornb=Bornb*(1d0-v)
      CALL  vvrho(keydi1, svar,  amel,  v, rhoini,beti,softi,hardi,ite)
      CALL  uurho(keydi2, svar1, amfin, u ,rhofin,betf,softf,hardf)
      fact1  =  vvv/(1d0-vvv*t)
* total integrand
      xdis3  =  fact1 *rhoini *rhofin *Bornb *rjac
      RETURN
 900  WRITE(nout,*) ' +++++ xdis3: something wrong!'
      WRITE(   6,*) ' +++++ xdis3: something wrong!'
      STOP
      END

*======================================================================
*======================================================================
*======================================================================
      SUBROUTINE figmis
*     *****************
*----------------------------------------------------------------------
* misterious exercise (probably from skrzypek) od legendre FUNCTION
* and zero-order final state exponentiation
*---------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      COMMON / keydis / keydis
      SAVE   / keydis /
      EXTERNAL vrh2

* books histo and fills it with a FUNCTION
      tmin =  -4d00
      tmax =  -0.001d0
      nbin = 200
*--- o(alf0)
      keydis=   10
      CALL GLK_BookFun1(2010,' o(alf0)  all  $',nbin,tmin,tmax,vrh2)
*---lege
      keydis=  510
      CALL GLK_BookFun1(2510,' o(alf0)  lege $',nbin,tmin,tmax,vrh2)

*--- v*rho(v)
*--- o(alf0)
      CALL GLK_Yminim( 0, 0.0000d0)
      CALL GLK_Ymaxim( 0, 0.1200d0)
      CALL GLK_Print(2010)
      CALL GLK_Print(2510)
      CALL GLK_Print(300)
* nonte carlo
*     CALL GLK_PlSet('dmod',1d0)
      CALL GLK_Plot(300,' ',' ',0)
* analytical
*     CALL GLK_PlSet('dmod',1d0)
      CALL GLK_Plot(2010,'s',' ',0)
*     CALL GLK_PlSet('dmod',1d0)
      CALL GLK_Plot(2510,'s',' ',0)
      END
      FUNCTION vrh2(t)
*     ****************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      vv= 10.0**t
      vrh2 = vvrh2(vv)*vv
      END
      FUNCTION vvrh2(uu)
*---------------------------------------------------------------------
* this is  a version which fits final state bramsstrahlung ll monte
* with betas in terms of sudakov variables calculated using initial
* virtual beams (not like in yfs3 where final state momenta are used)
* various types of the rho(u) distribution
*--------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER(pi= 3.1415926535897932d0)
      PARAMETER(ceuler =0.57721566d0)
      COMMON / weking / ene,amaz,gammz,amel,amfin,xk0,sinw2,ide,idf
      COMMON / inout  / ninp,nout
      COMMON / keydis / keydis
      COMMON  /xpar/  xpar(1000)
      SAVE
*
      alfinv = xpar(30)
      keyd = keydis
      alf1   = 1d0/pi/alfinv
      svar   = 4d0*ene**2
      sprim  = svar*(1-uu)
      bilg   = dlog(sprim/amfin**2)
      betf   = 2d0*alf1*(bilg-1d0)
* -------------yfs formula ------------------------------------------
      IF(keyd .GE. 10 .AND. keyd .LE. 12) THEN
*$$$     delb   = alf1*(0.5d0*bilg -1d0  +pi**2/3d0)
         delb   = alf1*(0.5d0*bilg )
         delb   = delb -betf/2 *dlog(1-uu)
         gamfac =exp(-ceuler*betf)/dpgamm(1d0+betf)
* ....zero   order
         IF(keyd   .EQ. 10)  THEN
            dels = 0d0
            delh = -betf/4 *log(1-uu)
* ....first  order
         ELSEIF(keyd   .EQ. 11)  THEN
            dels = betf/2
            delh = uu*(-1 +uu/2)
     $       +betf*(-uu**2/2 - 0.25d0*(1-uu)**2*log(1-uu) )
* ....second order
         ELSEIF(keyd   .EQ. 12)  THEN
            dels  = betf/2d0 +betf**2/8d0
            delh  = uu*(-1d0+uu/2d0)
     $        +betf*( -uu/2  +(2*uu-uu**2)/8*log(1-uu) )
         ENDIF
         distr= gamfac*exp(delb)*betf*uu**(betf-1d0)*(1 +dels +delh )
* -------------------------------------------------------------------
* -------------contributions  from various beta's -------------------
      ELSEIF(keyd .GE. 310 .AND. keyd .LE. 322)  THEN
*$$$$    delb   = alf1*(0.5d0*bilg -1d0  +pi**2/3d0)
         delb   = alf1*(0.5d0*bilg -0.5            )
         delb   = delb -betf/2 *dlog(1-uu)
         gamfac = exp(-ceuler*betf)/dpgamm(1d0+betf)
         soft   = 0d0
         delh   = 0d0
* ...beta0 first  order
         IF(keyd .EQ. 310) THEN
            soft = 1 + betf/2
            delh = -betf/4 *log(1-uu)
* ...beta1 first  order
         ELSEIF(keyd .EQ.  311)  THEN
            delh =
     $      uu*(-1d0+uu/2/(1+betf))*(1-0.5*betf*log(1-uu))
* ...beta0 second order
         ELSEIF(keyd .EQ. 320) THEN
            soft = 1 + betf/2  +betf**2/8
            delh = -betf/4 *log(1-uu)
* ...beta1 second order
         ELSEIF(keyd .EQ.  321)  THEN
            delh =
     $      uu*(1+0.5*betf)*(-1d0+uu/2/(1+betf))*(1-0.5*betf*log(1-uu))
     $       +0.25*betf*log(1-uu)*0.5*(1+(1-uu)**2)
*           delh = uu*(-1d0+uu/2d0)
*    $       +betf*(-uu/2 -uu**2/4 +(2+2*uu-uu**2)/8d0*log(1-uu))
* ...beta2 second order
         ELSEIF(keyd .EQ.  322)  THEN
            delh =    betf*uu**2/4d0
         ENDIF
         distr=  gamfac*exp(delb)*betf*uu**(betf-1d0)*(soft+delh)
* ...lege
         ELSEIF(keyd .EQ.  510)  THEN
* this not necessarily makes sense:
           gamfac = exp(-ceuler*betf)/dpgamm(1d0+betf)
           delb   = alf1*(0.5d0*bilg -0.5            )
           distr=  gamfac*exp(delb)*betf*alege(betf,uu)
      ELSE
         WRITE(nout,*) ' ===--->  wrong keydis in vvrh2'
         STOP
      ENDIF
      vvrh2 = distr
      END

      FUNCTION alege(beta,uu)
*     ***********************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      nu = beta/2
      x  = uu
      alege= dlege(nu-1d0,2d0/x-1d0)*x**(nu-1d0)*2d0**nu
      END

      FUNCTION dlege(nu,x)
*     ********************
* not normalized legendre FUNCTION d /gothic/
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER(amel=.511d-3,alfa=1d0/137d0,g=.57721566490153d0,
     1   pi=3.14159265358979d0,q=92d0)
      DOUBLE PRECISION nu

      an=1d0
      dlege=0d0
      DO 1 n=0,1000
      dlege=dlege +an
      an=an*(1d0+.5d0*nu+n)*.5d0*(1d0+nu+2d0*n)/(1.5d0+nu+n)/x**2/(n+1)
      IF(abs(an) .LT. 1d-50) go to 2
  1   CONTINUE
  2   dlege=dlege/(x)**(nu+1d0)
* normalization /not included/
*    @      sqrt(pi)/2d0**(nu+1d0)*dpgamm(nu+1d0)/dpgamm(nu+1.5d0)
* normalization /not included/
      END

*=======================================================================
*=======================================================================
*     END OF attics
*=======================================================================
*=======================================================================
