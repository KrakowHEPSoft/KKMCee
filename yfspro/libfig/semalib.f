      SUBROUTINE Semalib_MZplot(key,chak,id,title,idb)
*////////////////////////////////////////////////////////////////////////
*//                                                                    //
*//                                                                    //
*////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
* parameters
      SAVE
      CHARACTER*5 chak
      CHARACTER*80 title
      INTEGER key,id,idb
*---
      INCLUDE "semalib.h"
*
      DOUBLE PRECISION yy(400),er(400)
      DOUBLE PRECISION soft(5)
      LOGICAL GLK_Exist
*
      DOUBLE PRECISION  Mmin,Mmax,MMmin,MMmax
      DOUBLE PRECISION  xxmin,xxmax
      INTEGER nbin
      DOUBLE PRECISION  vvmin,vvmax,vvmax0,vv,amel,amfin
      DOUBLE PRECISION  CMSene,svar,sum,dist,beti,betf,hard
      INTEGER k,id1,ite
      CHARACTER*80 title2
*-----------------------------------------------------------------------
      WRITE(6,"(a,a,a,i10)") "Semalib_MZplot: ",chak," KeyDis=",Key
      CALL Semalib_SetDis(Key)

      IF(GLK_Exist(idb)) THEN
         CALL GLK_hinbo1(idb,title2,nbin,Mmin,Mmax)
      ELSE
         WRITE(*,*) "Semalib_MZplot: wrong idb= ",idb
      ENDIF
      CMSene  = m_CMSene
      vvmin   = m_vvmin
      vvmax   = m_vvmax
      amel    = m_amel
      amfin   = m_amfin
      svar   = CMSene**2
      DO k=1,nbin
         yy(k)= 0d0
         er(k)= 0d0
      ENDDO

      id1 = abs(id)
      IF(GLK_Exist(id1)) THEN
        WRITE(6,*) " ++++ Semalib_VVplot: warning deleted id1= ",id1
        CALL GLK_Delet(id1)
      ENDIF
      CALL GLK_Book1(id1,title,nbin,Mmin,Mmax)

*///////////////////////////////////////////////////////
*//       Integration over individual bins.           //
*///////////////////////////////////////////////////////
      vvmax0 = 1d0 -4*amfin**2/CMSene**2
      DO k=1,nbin
         MMmin = Mmin + (Mmax-Mmin)/nbin * (k - 1)
         MMmax = Mmin + (Mmax-Mmin)/nbin * (k - 0)
         xxmin = 1d0 -MMmax**2/svar
         xxmax = 1d0 -MMmin**2/svar
         xxmax = MIN(xxmax,vvmax)
         xxmax = MIN(xxmax,vvmax0)
         dist  = 0d0
         IF( xxmax  .GT.  xxmin ) THEN
            IF(    chak  .EQ.  "VRHO2" ) THEN
               CALL Semalib_vvchi2(xxmin,xxmax,dist)
            ELSEIF(chak  .EQ.  "URHO2") THEN
               CALL Semalib_uuchi2(xxmin,xxmax,dist)
            ELSEIF(chak  .EQ.  "XRHO2") THEN
               CALL Semalib_xxchi2(xxmin,xxmax,dist)
            ENDIF
         ENDIF
         dist = dist * nbin/(Mmax-Mmin)
         yy(k)= dist
      ENDDO
      CALL GLK_Pak (id1,yy)
      CALL GLK_Pake(id1,er)
      END


      SUBROUTINE Semalib_VVplot(key,chak,id,title,idb)
*////////////////////////////////////////////////////////////////////////
*//                                                                    //
*// Creates and fills histogram id1 with the plot                      //
*// chak= "VRHOS" rho(v)        initial miscelaneous and oldies        //
*//       "VRHO " rho(v)        initial                                //
*//       "URHO " rho(u)        final                                  //
*//       "VRHO2" rho(v)*Born   initial                                //
*//       "URHO2" rho(v)*Born   initial                                //
*//       "XRHO2" rho(u)*Born   initial+final gauss                    //
*// cumulatives                                                        //
*//       "VCHI2" sigma(xmax)   initial       gauss                    //
*//       "XCHI2" sigma(xmax)   initial+final gauss                    //
*// log10 x-scale for negative id                                      //
*//                                                                    //
*////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
* parameters
      SAVE
      CHARACTER*5 chak
      CHARACTER*80 title
      INTEGER key,id,idb
*---
      INCLUDE "semalib.h"
*
      DOUBLE PRECISION yy(400),er(400)
      DOUBLE PRECISION soft(5)
      LOGICAL GLK_Exist
*
      DOUBLE PRECISION  xmin,xmax
      INTEGER nbin
      DOUBLE PRECISION  vvmin,vvmax,vv,amel,amfin
      DOUBLE PRECISION  xxmin,xxmax
      DOUBLE PRECISION  CMSene,svar,sum,dist,beti,betf,hard
      INTEGER k,id1,ite
      CHARACTER*80 title2
*-----------------------------------------------------------------------
*
      WRITE(6,"(a,a,a,i10)") "Semalib_VVplot: ",chak," KeyDis=",Key
*
      CALL Semalib_SetDis(Key)
*
      CALL GLK_hinbo1(idb,title2,nbin,xmin,xmax)
*
      CMSene  = m_CMSene
      vvmin   = m_vvmin
      vvmax   = m_vvmax
      amel    = m_amel
      amfin   = m_amfin
*
      svar   = CMSene**2
*
      DO k=1,nbin
         yy(k)= 0d0
      ENDDO
*
      id1 = abs(id)
      IF(GLK_Exist(id1)) THEN
        WRITE(6,*) " ++++ Semalib_VVplot: warning deleted id1= ",id1
        CALL GLK_Delet(id1)
      ENDIF
      CALL GLK_Book1(id1,title,nbin,xmin,xmax)

      IF(chak .EQ. "VRHOS") THEN  !!!! historical and miscelaneous
         DO k=1,nbin
*     vv in middle of the bin !!!
            vv = xmin + (xmax-xmin)/nbin * (k-0.5d0)
            IF(id .LT. 0) vv = 10.**vv
            CALL Semalib_vvrhoS(m_KeyDis,svar,amel,vv,vvmin,dist)
*     ... log-scale induces vv-factor
            IF(id .LT. 0) dist = dist*vv
            yy(k)= dist
         ENDDO
      ELSEIF(chak .EQ. "VRHO ") THEN  !!!! historical, no Born(s(1-v))
         DO k=1,nbin
*     vv in middle of the bin !!!
            vv = xmin + (xmax-xmin)/nbin * (k-0.5d0)
            IF( id  .LT.  0) vv = 10.**vv
            CALL Semalib_vvrho(m_KeyDis,svar,amel,vv,dist,beti,soft,hard,ite)
*     ... log-scale induces vv-factor
            IF(id .LT. 0) dist = dist*vv
            yy(k)= dist
         ENDDO
      ELSEIF(chak .EQ. "URHO ") THEN  !!!!  historical
         DO k=1,nbin
*     vv in middle of the bin !!!
            vv = xmin + (xmax-xmin)/nbin * (k-0.5d0)
            IF( id  .LT.  0) vv = 10.**vv
            CALL Semalib_uurho(m_KeyDis,svar,amfin,vv,dist,betf,soft,hard)
*     ... log-scale induces vv-factor
            IF( id  .LT.  0) dist = dist*vv
            yy(k)= dist
         ENDDO
*///////////////////////////////////////////////////////
*//       NEW, integration over individual bins.      //
*///////////////////////////////////////////////////////
      ELSEIF(      (chak  .EQ.  "VRHO2")
     $       .OR.  (chak  .EQ.  "URHO2")
     $       .OR.  (chak  .EQ.  "XRHO2") ) THEN
         DO k=1,nbin
            xxmin = xmin + (xmax-xmin)/nbin * (k - 1)
            xxmax = xmin + (xmax-xmin)/nbin * (k - 0)
            IF(id .LT. 0) xxmin = 10.**xxmin
            IF(id .LT. 0) xxmax = 10.**xxmax
            IF(    chak  .EQ.  "VRHO2" ) THEN
               xxmax = min(xxmax,vvmax)
               xxmax = min(xxmax, (1-4*amfin**2/CMSene**2))
               IF( xxmax  .LE.  xxmin ) GOTO 300
               CALL Semalib_vvchi2(xxmin,xxmax,dist)
            ELSEIF(chak  .EQ.  "URHO2") THEN
               CALL Semalib_uuchi2(xxmin,xxmax,dist)
            ELSEIF(chak  .EQ.  "XRHO2") THEN
               xxmax = min(xxmax,vvmax)
               xxmax = min(xxmax, (1-4*amfin**2/CMSene**2))
               IF( xxmax  .LE.  xxmin ) GOTO 300
               CALL Semalib_xxchi2(xxmin,xxmax,dist)
            ENDIF
            dist = dist * nbin/(xmax-xmin)
            IF(id  .LT.  0) dist = dist/log(10.) !! bad convention
            yy(k)= dist
         ENDDO
 300     CONTINUE
*///////////////////////////////////////////////////////
*//                   CUMULATIVE                      //
*///////////////////////////////////////////////////////
      ELSEIF((chak .EQ. "VCHI2")
     $  .OR. (chak .EQ. "UCHI2")
     $  .OR. (chak .EQ. "XCHI2")
     $  .OR. (chak .EQ. "ZCHI1")
     $  .OR. (chak .EQ. "ZCHI2")) THEN
*     first bin separately
         xxmin = 0d0
         xxmax = xmin + (xmax-xmin)/nbin
         IF(id  .LT.  0) xxmax = 10.**xxmax
         IF(    chak  .EQ.  "VCHI2") THEN
            CALL Semalib_vvchi2(xxmin,xxmax,sum)
         ELSEIF(chak  .EQ.  "UCHI2") THEN
            CALL Semalib_uuchi2(xxmin,xxmax,sum)
         ELSEIF(chak  .EQ.  "XCHI2") THEN
            CALL Semalib_xxchi2(xxmin,xxmax,sum)
         ELSEIF(chak  .EQ.  "ZCHI1") THEN
            xxmax = xmin + (xmax-xmin)/nbin
            xxmax = 10.**xxmax/(1+10.**xxmax)
            CALL Semalib_vvchi2(xxmin,xxmax,sum)
         ELSEIF(chak  .EQ.  "ZCHI2") THEN
            xxmax = xmin + (xmax-xmin)/nbin
            xxmax = 10.**xxmax/(1+10.**xxmax)
            CALL Semalib_xxchi2(xxmin,xxmax,sum)
         ENDIF
         yy(1)= sum
         DO k=2,nbin
            xxmin = xmin + (xmax-xmin)/nbin *(k - 1)
            xxmax = xmin + (xmax-xmin)/nbin *(k - 0)
            IF(id .LT. 0) xxmin = 10d0**xxmin
            IF(id .LT. 0) xxmax = 10d0**xxmax
            dist = 0d0
            IF(    chak .EQ. "VCHI2") THEN
               xxmax = MIN(xxmax,vvmax)
               xxmax = MIN(xxmax, (1-4*amfin**2/CMSene**2))
               IF( xxmax .GT. xxmin ) CALL Semalib_vvchi2(xxmin,xxmax,dist)
            ELSEIF(chak .EQ. "UCHI2") THEN
               CALL Semalib_uuchi2(xxmin,xxmax,dist)
            ELSEIF(chak .EQ. "XCHI2") THEN
               xxmax = MIN(xxmax,vvmax)
               xxmax = MIN(xxmax, 1-4*amfin**2/CMSene**2)
               IF( xxmax .GT. xxmin ) CALL Semalib_xxchi2(xxmin,xxmax,dist)
            ELSEIF(chak .EQ. "ZCHI1") THEN
               xxmin = 10d0**xxmin/(1+10d0**xxmin)
               xxmax = 10d0**xxmax/(1+10d0**xxmax)
               xxmax = MIN(xxmax,vvmax)
               xxmax = MIN(xxmax, (1-4*amfin**2/CMSene**2))
               IF( xxmax .GT. xxmin ) CALL Semalib_vvchi2(xxmin,xxmax,dist)
            ELSEIF(chak .EQ. "ZCHI2") THEN
               xxmin = 10d0**xxmin/(1+10d0**xxmin)
               xxmax = 10d0**xxmax/(1+10d0**xxmax)
               xxmax = MIN(xxmax,vvmax)
               xxmax = MIN(xxmax, (1-4*amfin**2/CMSene**2))
               IF( xxmax .GT. xxmin ) CALL Semalib_xxchi2(xxmin,xxmax,dist)
            ENDIF
            sum = sum+dist
            yy(k)= sum
         ENDDO
      ELSE
         WRITE(6,*) "Semalib_VVplot: wrong chak=",chak
      ENDIF

      DO k=1,nbin
         er(k)= 0.d0
      ENDDO
      CALL GLK_Pak (id1,yy)
      CALL GLK_Pake(id1,er)
      END



      SUBROUTINE Semalib_Initialize(xpar_input)
*////////////////////////////////////////////////////////////////////////
*//                                                                    //
*//   Initialize xpar and other data members                           //
*//                                                                    //
*////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      DOUBLE PRECISION xpar_input(*)
*---
      INCLUDE "semalib.h"
      INTEGER i
      INTEGER KFbeam,KFfin
*-----------------------------------------------------------
      m_out = 16
      DO i=1,imax
         m_xpar(i)=xpar_input(i)
      ENDDO
      m_CMSene = m_xpar( 1)
      m_vvmin  = m_xpar(16)
      m_vvmax  = m_xpar(17)
      m_Xenph  = m_xpar(40)

      WRITE(*,*) ' ==================================================='
      WRITE(*,*) ' =================SEMALIB==========================='
      WRITE(*,*) ' ###### m_Xenph =',m_Xenph
      WRITE(*,*) ' ==================================================='

*
      KFbeam = 11           ! KF=11 is electron
      m_amel   = m_xpar(500+10*KFbeam+6)
      KFfin = 13            ! KF=13 is muon
      m_amfin  = m_xpar(500+10*KFfin+6)
*
      m_KeyZet = m_xpar(501)
      m_Zmass  = m_xpar(502)
      END





      SUBROUTINE Semalib_vvchi2(vvmin,vvmax,chi2)
*///////////////////////////////////////////////////////////////////////////
*//                                                                       //
*//                   ISR                                                 //
*//                                                                       //
*//   sigma(vmin,vmax)                                                    //
*//   Integral over vvrho(v)*Born(s*(1-v))                                //
*//                                                                       //
*///////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      SAVE
      DOUBLE PRECISION  vvmin,vvmax,chi2
*---
      INCLUDE "semalib.h"
*---
      DOUBLE PRECISION   Semalib_vvchi2_fun
      EXTERNAL Semalib_vvchi2_fun
*---
      DOUBLE PRECISION   svar,amel,v1,v2
      DOUBLE PRECISION   sf(5),z1,z2,prec,eps,chi
      DOUBLE PRECISION   Born,beti,result,dum1,dum2
      DOUBLE PRECISION   Semalib_Born
      INTEGER  KeyZet, ite
*---------------------------------------------------------------
      KeyZet = m_KeyZet
      svar   = m_CMSene**2
      amel   = m_amel
      Born   = Semalib_Born(svar)
      v1=vvmin
      v2=vvmax
*//////////////////////////////////////////////////
*//  Integration over hard part.                 //
*//  v=0 soft singular part is subtracted        //      
*//////////////////////////////////////////////////
      IF(  KeyZet.GT.0 ) THEN
*        Change of variables to get rid of 1/(1-v) from photon propagator
         z1= -log(1-v1)
         z2= -log(1-v2)
         prec  = 1d-3 *Born     ! good for vmax=1, fast!
         CALL Semalib_GausJad(Semalib_vvchi2_fun,z1,z2,prec,result)
      ELSE
         eps = 0.0001           ! good for vmax=0.9999
         eps = 0.000001         ! good for vmax=1 (slow)
         CALL Semalib_Gauss(Semalib_vvchi2_fun,v1,v2,eps,result)
         prec  = 1d-3 *Born     ! good for vmax=0.9999
         prec  = 1d-4 *Born     ! good for vmax=1 (slow)
******   CALL Semalib_GausJad(Semalib_vvchi2_fun,v1,v2,prec,result)
      ENDIF
******WRITE(6,"(a,4f20.12)") "x,error=   ",x,(resul3-result)/result
*//////////////////////////////////////////////////
*//     Combine soft and hard integrals          //
*//////////////////////////////////////////////////
      CALL Semalib_vvrho(m_KeyDis,svar,amel,1d-10,dum1,beti,sf,dum2,ite)
      IF( ite .EQ. 0) THEN
* EXPonentiated
         chi =  Born*sf(1)*v2**beti +result
         IF(v1 .NE. 0d0 ) chi = chi
     $         -Born*sf(1)*v1**beti
      ELSEIF( ite .EQ. 100) THEN
* UNEXP zero photons
         chi= Born*sf(1)
         IF(v1 .NE. 0d0 ) chi = 0d0
      ELSEIF( ite .EQ. 101) THEN
* UNEXP one photon
         chi= Born*(sf(1)*(1 +beti*log(v2)) +sf(2)) +result
         IF(v1 .NE. 0d0 ) chi = chi
     $       -Born*(sf(1)*(1 +beti*log(v1)))
      ELSEIF( ite .EQ. 102) THEN
* UNEXP Two photons
         chi= Born*(sf(1)*(1 +beti*log(v2) +1/2d0*beti**2*log(v2)**2) 
     $             +sf(2)*(1 +beti*log(v2)) +sf(3)) +result
         IF(v1 .NE. 0d0 ) chi = chi
     $       -Born*(sf(1)*(1 +beti*log(v1) +1/2d0*beti**2*log(v1)**2)
     $             +sf(2)*(1 +beti*log(v1)) +sf(3))
      ELSE
         WRITE(*,*) 'Semalib_vvchi2: STOP'
         STOP
      ENDIF
      chi2 = chi
      END

      DOUBLE PRECISION FUNCTION Semalib_vvchi2_fun(rr)
*///////////////////////////////////////////////////////////////////////////
*//   Integrand, hard part of  vvrho(v)*Born(s*(1-v))                     //
*//   Soft v=0 part subtracted according to type of bremsstrahlung        //
*///////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      DOUBLE PRECISION rr
      SAVE
*---
      INCLUDE "semalib.h"
*---
      DOUBLE PRECISION   svar,amel
      DOUBLE PRECISION   softi(5),Born,Born2,dJacob,vv,dist,hard,beti,rho
      DOUBLE PRECISION   Semalib_Born
      INTEGER  KeyZet,ite
*----------------------------------------------------------------------
      svar   = m_CMSene**2
      amel   = m_amel
      KeyZet = m_KeyZet
      IF(  KeyZet.GT.0 ) THEN
*        Change of variables to eliminate 1/(1-v)
         vv= 1d0 -exp(-rr)
         dJacob = 1-vv
      ELSE
         vv=rr
         dJacob = 1d0
      ENDIF
* Basic elements of the ISR function
      CALL  Semalib_vvrho(m_KeyDis,svar,amel,vv,rho,beti,softi,hard,ite)
*
      Born  = Semalib_Born(svar)
      Born2 = Semalib_Born(svar*(1-vv))
*
      IF( KeyZet .LT. 0) Born2 = Born
*
      IF( ite .EQ. 0) THEN
         dist = rho*Born2 -softi(1)*beti*vv**(beti-1)*Born
      ELSEIF( ite .EQ. 100) THEN
         dist= 0
      ELSEIF( ite .EQ. 101) THEN
         dist= Born2*rho  -Born*( softi(1)*beti/vv )
      ELSEIF( ite .EQ. 102) THEN
         dist= Born2*rho -Born*( softi(1)*beti/vv*(1 +beti*log(vv)) 
     $                          +softi(2)*beti/vv )
      ELSE
         WRITE(*,*) 'Semalib_vvchi2_fun: STOP'
         STOP
      ENDIF
      Semalib_vvchi2_fun = dist*dJacob
      END

      SUBROUTINE Semalib_uuchi2(uumin,uumax,chi)
*///////////////////////////////////////////////////////////////////////////
*//                                                                       //
*//                   FSR                                                 //
*//                                                                       //
*//   sigma(umin,umax)                                                    //
*//                                                                       //
*//   Integral over uurho*Born(s)                                         //
*//                                                                       //
*///////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      DOUBLE PRECISION    uumin,uumax,chi
      SAVE
*---
      DOUBLE PRECISION             Born,betf,soft
      COMMON /lc_uuchi2/ Born,betf,soft
*---
      INCLUDE "semalib.h"
*---
      DOUBLE PRECISION   Semalib_uuchi2_fun
      EXTERNAL Semalib_uuchi2_fun
*---
      DOUBLE PRECISION   svar,amfin,prec,eps,result,dum1,dum2
      DOUBLE PRECISION   Semalib_Born
*---------------------------------------------------------------------------
      svar   = m_CMSene**2
      amfin  = m_amfin
      Born   = Semalib_Born(svar)
      CALL Semalib_uurho(m_KeyDis,svar,amfin,1d-10,dum1,betf,soft,dum2)

      eps = 0.00001
      CALL Semalib_Gauss(Semalib_uuchi2_fun,uumin,uumax,eps,result)
      prec  = -1d-3
      prec  =  1d-3 *Born
******CALL Semalib_GausJad(Semalib_uuchi2_fun,uumin,uumax,prec,result)

      IF( uumin  .EQ.  0d0 ) THEN
         chi = soft*uumax**betf*Born + result
      ELSE
         chi = soft*uumax**betf*Born -soft*uumin**betf*Born + result
      ENDIF
*c      WRITE(6,"(a,4f20.12)") "x,chi=",uumax,chi,chi/Born
      END

      DOUBLE PRECISION FUNCTION Semalib_uuchi2_fun(u)
*///////////////////////////////////////////////////////////////////////////
*//                                                                       //
*///////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      DOUBLE PRECISION   u
*---
      DOUBLE PRECISION             Born,betf,soft
      COMMON /lc_uuchi2/ Born,betf,soft
*---
      INCLUDE "semalib.h"
*---
      DOUBLE PRECISION   svar,amfin
      DOUBLE PRECISION   rho,dum1,dum2,dum3
*---------------------------------------------------------------------------
      svar   = m_CMSene**2
      amfin  = m_amfin
      CALL  Semalib_uurho(m_KeyDis,svar,amfin,u,rho,dum1,dum2,dum3)
      Semalib_uuchi2_fun = rho*Born -soft*betf*u**(betf-1)*Born
      END


      SUBROUTINE Semalib_xxchi2(xxmin,xxmax,chi)
*///////////////////////////////////////////////////////////////////////////
*//                                                                       //
*//                ISR * FSR                                              //
*//                                                                       //
*//             sigma(Xmin,Xmax)                                          //
*//                                                                       //
*///////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      SAVE
      DOUBLE PRECISION xxmin,xxmax,chi
*---
      INCLUDE "semalib.h"
*---
      EXTERNAL Semalib_xxchi2_fun
      DOUBLE PRECISION   Semalib_xxchi2_fun
*---
      DOUBLE PRECISION               svar,amel,amfin
      DOUBLE PRECISION   Born,x1,x2,z1,z2
      DOUBLE PRECISION   eps,prec,beti,soft,dum1,dum2,result,fleps
      DOUBLE PRECISION   Semalib_Born
      INTEGER  KeyZet
*-------------------------------------------------------------------
      KeyZet= m_KeyZet
      svar  = m_CMSene**2
      amel  = m_amel
      amfin = m_amfin
      Born  = Semalib_Born(svar)
      x1=xxmin
      x2=xxmax
*////////////////////////////////////////////
*//        External hard integral          //
*////////////////////////////////////////////
      IF(  KeyZet.GT.0 ) THEN
*     Change of variables to get rid of 1/(1-v) from photon propagator
         z1= -log(1-x1)
         z2= -log(1-x2)
         eps = 0.0001           ! good for vmax=0.9999
         eps = 0.000001         ! good for vmax=1 
         CALL Semalib_Gauss2(Semalib_xxchi2_fun,z1,z2,eps,result)  ! Semalib_Gauss2 is clone of Semalib_Gauss
         prec  = 1d-3 *Born     ! good for vmax=0.9999
******   CALL Semalib_GausJad2(Semalib_xxchi2_fun,z1,z2,prec,result) ! Semalib_GausJad2 is clone of Semalib_GausJad
      ELSE
         eps = 0.0001           ! good for vmax=0.9999
         eps = 0.000001         ! good for vmax=1
         CALL Semalib_Gauss2(Semalib_xxchi2_fun,x1,x2,eps,result) ! Semalib_Gauss2 is clone of Semalib_Gauss
         prec  = 1d-3 *Born     ! good for vmax=0.9999
         prec  = 1d-4 *Born     ! good for vmax=1 (slow)
******   CALL Semalib_GausJad2(Semalib_xxchi2_fun,x1,x2,prec,result) ! Semalib_GausJad2 is clone of Semalib_GausJad
      ENDIF
*////////////////////////////////////////////
*//    Add soft part integrated by hand    //
*////////////////////////////////////////////
      fleps = 1d-7
      CALL Semalib_xxrho2(fleps,dum1,beti,soft,dum2)
      IF( xxmin .EQ. 0d0 ) THEN
         chi   = soft*xxmax**beti + result
      ELSE
         chi   = soft*xxmax**beti -soft*xxmin**beti + result
      ENDIF
*
******WRITE(6,"(a,4f20.12)") " Semalib_xxchi2: xxmax,chi= ",xxmax,chi,chi/Born
      END

      DOUBLE PRECISION FUNCTION Semalib_xxchi2_fun(rr)
*///////////////////////////////////////////////////////////////////////////
*//                                                                       //
*///////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      SAVE
      DOUBLE PRECISION rr
*---
      INCLUDE "semalib.h"
*---
      DOUBLE PRECISION   svar,amel,amfin
      DOUBLE PRECISION   vv,rho,beti,soft,hard,dJacob
      INTEGER  KeyZet
*------------------------------------------------------------------
      svar  = m_CMSene**2
      amel  = m_amel
      amfin = m_amfin
      KeyZet = m_KeyZet
      IF(  KeyZet.GT.0 ) THEN
*        Change of variables to eliminate 1/(1-v)
         vv= 1d0 -exp(-rr)
         dJacob = 1-vv
      ELSE
         vv=rr
         dJacob = 1d0
      ENDIF
*---
      CALL  Semalib_xxrho2(vv,rho,beti,soft,hard)
      Semalib_xxchi2_fun = Hard*dJacob
      END

      SUBROUTINE Semalib_xxrho2(vv9,rho,beti,soft,hard)
*///////////////////////////////////////////////////////////////////////////
*//                                                                       //
*//           ISR * FSR Convolution                                       //
*//                                                                       //
*//               d_sigma/d_X                                             //
*//                                                                       //
*//   Integral over  vvrho(v)*uurho(u)*Born(s*(1-v)*delta(X-u+v-u*v)      //
*//   where    v = X*t                                                    //
*//            u = X*(1-t)/(1-vv*t)                                       //
*//   t is in (0,1) range, X= u+v-u*v                                     //
*//                                                                       //
*///////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      DOUBLE PRECISION  vv9,rho,beti,soft,hard
*---
      SAVE
      INCLUDE "semalib.h"
*---
      DOUBLE PRECISION               vv,bti,btf,c1,c2
      COMMON / lc_xxrho2 / vv,bti,btf,c1,c2
*---
      DOUBLE PRECISION  softi(5)
*---
      DOUBLE PRECISION   Semalib_xxrho2_fun
      EXTERNAL Semalib_xxrho2_fun
*---
      DOUBLE PRECISION   svar,amel,amfin
      INTEGER  KeyZet,keydi1,keydi2,ite
      DOUBLE PRECISION   Born1,Born2
      DOUBLE PRECISION   dum3,eps,prec,softf0,btf0,result,rhof,rhoi
      DOUBLE PRECISION   svar1,cregul,dum1,softf
      DOUBLE PRECISION   Semalib_Born,Semalib_Gamma
*-----------------------------------------------------------------------------
      vv  = vv9
      KeyZet = m_KeyZet

      svar  = m_CMSene**2
      amel  = m_amel
      amfin = m_amfin


      keydi1= MOD(m_KeyDis,1000000)/1000
      keydi2= MOD(m_KeyDis,1000)

      svar1 = svar*(1-vv)
      Born1 = Semalib_Born(svar )
      Born2 = Semalib_Born(svar1)
*
      IF( keyzet  .LT.  0) Born2 = Born1
*
      CALL Semalib_vvrho(keydi1,svar, amel,vv,rhoi,bti,softi,dum3,ite)
      CALL Semalib_uurho(keydi2,svar,amfin,vv,rhof,btf,softf,dum3)
*///////////////////////////////////////////////////////
*//     Define Soft part and integrate analyticaly    //
*///////////////////////////////////////////////////////
* Residue at t=0 and  t=1
* Note btf, softf has to be calculated at svar1, u=0 !!!
      CALL Semalib_uurho(keydi2,svar1,amfin,1d-10,dum1,btf,softf,dum3)
      IF(  KeyZet.GT.0 ) THEN
         c1 = (vv/(1-vv))**bti    *bti*softi(1)    *rhof*Born1
         c2 =          vv**btf    *btf*softf       *rhoi*Born2
      ELSE
         c1 =          vv**bti    *bti*softi(1)    *rhof*Born1
         c2 = (vv/(1-vv))**btf    *btf*softf       *rhoi*Born2
      ENDIF
      cregul= 1/(bti*btf)*Semalib_Gamma(1+bti)*Semalib_Gamma(1+btf)
     $                   /Semalib_Gamma(1+bti+btf)
     $                  *(c1*btf+c2*bti)
*///////////////////////////////////////////////////////
*//        Integrate Hard part numericaly             //
*///////////////////////////////////////////////////////
      IF( abs(cregul) .GT. 0.0001) prec=  1d-4 * abs(cregul)
      prec=  1d-3*Born1         ! good for vmax=0.9999
******CALL Semalib_GausJad(Semalib_xxrho2_fun,0d0,1d0,prec,result) ! Infinitely slow
      eps = 0.0001              ! good for vmax=0.9999
      eps = 0.000001            ! good for vmax=1
      eps = 0.00001             ! sufficient for vmax=1, tested! twice faster
      CALL Semalib_Gauss(Semalib_xxrho2_fun,0d0,1d0,eps ,result)
*
*  Total result 
      rho = cregul + result
*///////////////////////////////////////////////////////
*//   Define Soft and Hard parts                      //
*///////////////////////////////////////////////////////
*     power and residue at vv-->0 limit
      CALL Semalib_uurho(keydi2,svar,amfin,1d-10,dum1,btf0,softf0,dum3)
      beti = bti + btf0
      soft = softi(1)*softf0*Born1
     $       *Semalib_Gamma(1+bti)*Semalib_Gamma(1+btf0)/Semalib_Gamma(1+bti+btf0)
      hard = rho - soft*beti*vv**(beti-1)
      END

      DOUBLE PRECISION  FUNCTION Semalib_xxrho2_fun(r)
*///////////////////////////////////////////////////////////////////////////
*//                                                                       //
*//  Integrand for xxrho2, ISR * FSR Convolution                          //
*//                                                                       //
*///////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      SAVE
      DOUBLE PRECISION r
*---
      INCLUDE "semalib.h"
*---
      DOUBLE PRECISION               vv,bti,btf,c1,c2
      COMMON / lc_xxrho2 / vv,bti,btf,c1,c2
*---
      DOUBLE PRECISION   svar,amel,amfin
      DOUBLE PRECISION   tt
      DOUBLE PRECISION   v,rhoi,beti,softi(5),hardi
      DOUBLE PRECISION   u,rhof,betf,softf,hardf,distr,svar1
      DOUBLE PRECISION   dJacob,dJacob1,Bornb
      DOUBLE PRECISION   Semalib_Born
      INTEGER  ite,KeyZet,KeyDi1,KeyDi2
*-----------------------------------------------------------------
      KeyZet = m_KeyZet
      svar  = m_CMSene**2
      amel  = m_amel
      amfin = m_amfin

      KeyDi1= MOD(m_KeyDis,1000000)/1000
      KeyDi2= MOD(m_KeyDis,1000)

      tt = MAX(r, 1d-9)
      tt = MIN(tt, 1d0 -1d-9)
      IF(  KeyZet.GT.0 ) THEN
*        In principle this should better eliminate 1/(1-v) in Born
         v = vv*tt/(1-vv*(1-tt))
         u = vv*(1-tt)
         dJacob = (1-v)  *vv/(1-vv)
      ELSE
         v = vv*tt
         u = vv*(1-tt)/(1-vv*tt)
         dJacob =   vv/(1-v)
      ENDIF
*
* Born
      svar1 = svar*(1-v)
      Bornb = Semalib_Born(svar1)
      IF( KeyZet .LT. 0) Bornb = Semalib_Born(svar)
*
* Bremsstrahlung functions
      CALL Semalib_vvrho(keydi1,svar, amel, v,rhoi,beti,softi,hardi,ite)
      CALL Semalib_uurho(keydi2,svar1,amfin,u,rhof,betf,softf,hardf)
*
* Basic integrand
      distr   = dJacob *rhoi*rhof *Bornb
*
* Minus subtractions. Residua c1,c2 defined in Semalib_xxrho2
      Semalib_xxrho2_fun   = distr   -c1*tt**(bti-1) *(1-tt)**btf
     $                               -c2*tt**bti     *(1-tt)**(btf-1)
      END




      SUBROUTINE Semalib_vvrho(KeyIni,svar,amel,vv,distr,beti,softi,hard,ite)
*///////////////////////////////////////////////////////////////////////////
*//                                                                       //
*//                ISR  Bremsstrahlung functions                          //
*//                                                                       //
*///////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      DOUBLE PRECISION  svar,amel,vv,distr,beti,softi(5),hard
      INTEGER KeyIni,ite
*---
      DOUBLE PRECISION  pi,zeta2,zeta3,ceuler
      PARAMETER( pi= 3.1415926535897932d0 )
      PARAMETER(zeta2  =1.6449340668482264D0)
      PARAMETER(zeta3  =1.2020569031595943D0)
*****>PARAMETER(zeta2 = 1d0/6*pi**2)
      PARAMETER(ceuler =0.57721566d0)
      SAVE
*---
      INCLUDE "semalib.h"
*---
      DOUBLE PRECISION   alfinv,vvmin,alf1,bilg,beti2,ffact,damel
      DOUBLE PRECISION   delh,delb,gamfac,delgri,etah,facgri,dels
      DOUBLE PRECISION   Semalib_Gamma,Semalib_DiLog
      DOUBLE PRECISION   DilJac0
      INTEGER  keywtm,keyd 
*------------------------------------------------------------------
      alfinv = m_xpar(30)
      vvmin  = m_xpar(16)
      keywtm = m_xpar(26)
*
      keyd = KeyIni
      ite  = 0
      alf1   = 1d0/pi/alfinv
      bilg   = dlog(svar/amel**2)
      beti   = 2d0*alf1*(bilg-1d0)
      beti2  = 2d0*alf1*bilg
      beti2  = beti2 * m_Xenph  !!! photon enhancement factor
      IF(keywtm .EQ. 1) beti2=beti
* ---------------------------------------------------------------
*     Technical tests
* ---------------------------------------------------------------
      IF(keyd .GE. 0 .AND. keyd .LE. 100) THEN
         gamfac = exp(-ceuler*beti)/Semalib_Gamma(1d0+beti)
         DilJac0 = (1d0+1d0/SQRT(1d0-vv))/2d0
         IF(keyd   .EQ. 51)  THEN
            delh = -beti/4 *log(1-vv)
     $           -1/2d0 *alf1 *log(1-vv)**2 ! NLL from phase space
            softi(1)  = gamfac
            hard  = gamfac*beti*vv**(beti-1d0) *delh
            distr = softi(1) *beti*vv**(beti-1d0) + hard
         ELSEIF(keyd   .EQ. 52)  THEN
            IF(vv .GT. vvmin) THEN
               damel=beti2/beti*(vv/vvmin)**(beti2-beti)
            ELSE
               damel=1d0
            ENDIF
            softi(1)  = 1d0
c[[[            distr = beti*vv**(beti-1d0)*damel
            distr = beti*vv**(beti-1d0)*damel*DilJac0
            hard  = distr - softi(1) *beti*vv**(beti-1d0)
         ELSEIF(keyd   .EQ. 53)  THEN
            softi(1)  = 1d0
            distr = beti*vv**(beti-1d0)
            hard  = 0d0
         ELSEIF(keyd   .EQ. 54)  THEN
            softi(1)  = gamfac
            distr = gamfac*beti*vv**(beti-1d0)
            hard  = 0d0
         ELSEIF(keyd   .EQ. 55)  THEN
            softi(1)  = 0d0
            delh = -beti/4 *log(1-vv)
     $           -1/2d0 *alf1 *log(1-vv)**2 ! NLL from phspace
            hard  = gamfac*beti*vv**(beti-1d0) *delh
            distr = hard
         ENDIF
*============================================================================
*            YFS exp. up to 3-rd order
*============================================================================
      ELSEIF(keyd .GE. 300 .AND. keyd .LE. 305) THEN
         delb   = beti/4d0 +alf1*(-0.5d0  +pi**2/3d0)
         gamfac = exp(-ceuler*beti)/Semalib_Gamma(1d0+beti)
         ffact  = gamfac*exp(delb)
*-----------------------------------------------------------------
         IF(keyd   .EQ. 300)  THEN
* ....zero   order
            dels = 0d0
            delh = -beti/4 *log(1-vv)
            softi(1) = ffact*(1 +dels)
            distr = ffact*beti*vv**(beti-1d0)*(1 +dels +delh)
*-----------------------------------------------------------------
         ELSEIF(keyd   .EQ. 301)  THEN
* ....first  order
            dels = beti/2
            delh = vv*(-1 +vv/2)
     $           -beti/2*vv**2 - beti/4*(1-vv)**2*log(1-vv)
            softi(1) = ffact*(1 +dels)
            distr = ffact*beti*vv**(beti-1d0)*(1 +dels +delh)
*-----------------------------------------------------------------
* ....second order
         ELSEIF(keyd   .EQ. 302)  THEN
            dels = beti/2d0 +beti**2/8d0
            delh = vv*(-1d0+vv/2d0)
     $           +beti*.5d0*(-0.25d0*(4d0-6d0*vv+3d0*vv**2)*dlog(1d0-vv)-vv)
            softi(1) = ffact*(1 +dels)
            distr = ffact*beti*vv**(beti-1d0)*(1 +dels +delh)
*-----------------------------------------------------------------
* O(alf3)LL added
         ELSEIF(keyd   .EQ. 303)  THEN
            dels = beti/2d0 +beti**2/8d0
            delh = vv*(-1d0+vv/2d0)
     $           +beti*.5d0*(
     $               -0.25d0*(4d0-6d0*vv+3d0*vv**2)*dlog(1d0-vv) -vv )
*!!!!!!!!! O(alf3)LL part, worth 0.5% for 200GeV!!!
     @           +beti**2*(
     @              +(3d0*vv-2d0)*vv/16*dlog(1d0-vv)
     @              +(8d0-14d0*vv+7d0*vv**2)/96*dlog(1d0-vv)**2
     @              +vv**2/8d0
     @              +(2d0-vv)*vv/8*Semalib_DiLog(vv)  
     @           )
            softi(1) = ffact*(1 +dels)
            distr = ffact*beti*vv**(beti-1d0)*(1 +dels +delh)
*-----------------------------------------------------------------
         ELSEIF(keyd   .EQ. 304)  THEN
* O(alf3)LL, exp(3/4*gam) puled out,  No O(alf2)NLL yet
            delgri = beti*3/4d0 +alf1*(-0.5d0  +pi**2/3d0)
            facgri = gamfac*exp(delgri)
            dels = 0d0
            etah = (1+(1-vv)**2)/2d0
*     O(alf2)LL part
     @           +beti*( 
     @              -(1d0+3d0*(1d0-vv)**2)/8d0*log(1d0-vv) 
     @              -.25d0*vv**2)
*     O(alf3)LL part
     @           +beti**2*(
     @              +(3d0*vv-2d0)*vv/16*dlog(1d0-vv)
     @              +(8d0-14d0*vv+7d0*vv**2)/96*dlog(1d0-vv)**2
     @              +vv**2/8d0
     @              +(2d0-vv)*vv/8*Semalib_DiLog(vv)  )
            softi(1) = facgri*(1 +dels)
            distr    = facgri*beti*vv**(beti-1d0)*(etah +dels)
*-----------------------------------------------------------------
         ELSEIF(keyd   .EQ. 305)  THEN
* The best pure photonic formula for ISR !!!!!
* O(alf3)LL, exp(3/4*gam) puled out,  O(alf2)NLL added!!!
            delgri = beti*3/4d0 +alf1*(-0.5d0  +pi**2/3d0)
            facgri = gamfac*exp(delgri)
            dels =
     @           beti*alf1*(3d0/32d0 -pi**2/8d0 +3d0/2d0*zeta3)
*
            etah = (1+(1-vv)**2)/2d0
*     O(alf2)LL part
     @           +beti*( 
     @              -(1d0+3d0*(1d0-vv)**2)/8d0*log(1d0-vv) 
     @              -.25d0*vv**2)
*     O(alf3)LL part
     @           +beti**2*(
     @              +(3d0*vv-2d0)*vv/16*dlog(1d0-vv)
     @              +(8d0-14d0*vv+7d0*vv**2)/96*dlog(1d0-vv)**2
     @              +vv**2/8d0
     @              +(2d0-vv)*vv/8*Semalib_DiLog(vv)  )
*     O(alf2)NLL part
     @           +alf1*1/8d0*(
     @               4*(1+(1-vv)**2)*(Semalib_DiLog(vv)+log(vv)*log(1-vv))
     @              -(1+3*(1-vv)**2)*(log(1-vv))**2
     @              +2*(3+2*(1-vv)+(1-vv)**2)*log(1-vv)
     @              +2*vv*(3-2*(1-vv)) )
            softi(1) = facgri*(1 +dels)
            distr    = facgri*beti*vv**(beti-1d0)*(etah +dels)
         ENDIF
*----------------------------
*     modernized hard
         hard     = distr -softi(1)*beti*vv**(beti-1d0)
*============================================================================
*  -------------- YFS exclusive exponentiation --------------------
*  -------------contributions  from various beta's ----------------
*============================================================================
      ELSEIF(keyd .GE. 310 .AND. keyd .LE. 322)  THEN
         delb   = alf1*(0.5d0*bilg -1d0  +pi**2/3d0)
         gamfac = exp(-ceuler*beti)/Semalib_Gamma(1d0+beti)
         ffact  = gamfac*exp(delb)
         softi(1)   = 0d0
         delh   = 0d0
* ...beta0 first  order
         IF(keyd .EQ. 310) THEN
            softi(1) = 1 + beti/2
            delh = -beti/4d0 *log(1-vv)
* ...beta1 first  order
         ELSEIF(keyd .EQ.  311)  THEN
            delh =
     $      vv*(-1d0+vv/2/(1+beti))*(1-0.5*beti*log(1-vv))
* ...beta0 second order
         ELSEIF(keyd .EQ. 320) THEN
            softi(1) = 1 + beti/2  +beti**2/8
            delh = -beti/4 *log(1-vv)
* ...beta1 second order
         ELSEIF(keyd .EQ.  321)  THEN
            delh = vv*(-1+vv/2)
     $      -beti*vv/2 -beti*vv**2/4+beti/8*(-2+6*vv-3*vv**2)*log(1-vv)
* ...beta2 second order
         ELSEIF(keyd .EQ.  322)  THEN
            delh =    beti*  vv**2/4d0
         ENDIF
         softi(1)  = softi(1)*ffact
         hard  = ffact*beti*vv**(beti-1d0) *delh
         distr = softi(1) *beti*vv**(beti-1d0) + hard
*----------------------------------------------------------------
* Special NLL case for beta0
      ELSEIF(keyd .GE. 400 .AND. keyd .LE. 422)  THEN
         delb   = alf1*(0.5d0*bilg -1d0  +pi**2/3d0)
         gamfac = exp(-ceuler*beti)/Semalib_Gamma(1d0+beti)
         ffact  = gamfac*exp(delb)
         softi(1)   = 0d0
         delh   = 0d0
         IF(keyd   .EQ. 400)  THEN
* ....zero   order
            dels = 0d0
            delh = -beti/4 *log(1-vv)
     $           -1/2d0 *alf1 *log(1-vv)**2 !!!! NLL from phsp
            softi(1) = ffact*(1 +dels)
            distr = ffact*beti*vv**(beti-1d0)*(1 +dels +delh)
         ENDIF
*     modernized hard
         hard     = distr -softi(1)*beti*vv**(beti-1d0)
*----------------------------------------------------------------
*     UNEXPonentiated
*----------------------------------------------------------------
      ELSEIF(keyd .GE. 600 .AND. keyd .LE. 700)  THEN
         IF(keyd .EQ.  660)  THEN
            ite = 100
            softi(1) =  1
            hard    = 0
            distr   = softi(1)
         ELSEIF(keyd .EQ.  661)  THEN
            ite = 101
            softi(1) =  1
            softi(2) =  0
            hard     = 0
            distr    = hard + softi(1)*beti/vv
         ELSEIF(keyd .EQ.  662)  THEN
            ite = 102
            softi(1) =  1
            softi(2) =  0
            softi(3) =  -1/12d0*pi**2 *beti**2
            hard     = -beti**2/4d0 *log(1-vv)/vv
     $           -1/2d0 *alf1*beti *log(1-vv)**2/vv ! NLL from phspace
            distr    = Hard
     $           +softi(1)*beti/vv*(1 +beti*log(vv)) 
     $           +softi(2)*beti/vv
         ELSEIF(keyd .EQ.  663)  THEN
            ite = 102
            softi(1) =  0
            softi(2) =  0
            softi(3) =  0
            hard     = -beti**2/4d0 *log(1-vv)/vv
     $           -1/2d0 *alf1*beti *log(1-vv)**2/vv ! nll from phspace
            distr    = Hard
     $           +softi(1)*beti/vv*(1 +beti*log(vv)) 
     $           +softi(2)*beti/vv
         ENDIF
      ELSE
         WRITE(*,*) 'Semalib_vvrho: ===--->  wrong KeyIni in vvrho',keyd
         STOP
      ENDIF
      END


      SUBROUTINE Semalib_uurho(KeyFin,svar,amfin,uu,distr,betf,soft,hard)
*///////////////////////////////////////////////////////////////////////////
*//                                                                       //
*//                FSR  Bremsstrahlung functions                          //
*//                                                                       //
*///////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INTEGER  KeyFin
      DOUBLE PRECISION   svar,amfin,uu,distr,betf,soft,hard
      SAVE
*---
      DOUBLE PRECISION  pi,ceuler
      PARAMETER(pi= 3.1415926535897932d0)
      PARAMETER(ceuler =0.57721566d0)
*---
      INCLUDE "semalib.h"
*---
      DOUBLE PRECISION  alfinv,sprim,bilg,delb,gamfac,ffact,alf1
      DOUBLE PRECISION  dels,delh
      DOUBLE PRECISION  Semalib_Gamma
      INTEGER keyd
*----------------------------------------------------------------------
      alfinv = m_xpar(30)
      keyd   = KeyFin
      alf1   = 1d0/pi/alfinv
      sprim  = svar*(1-uu)
      bilg   = dlog(sprim/amfin**2)
      betf   = 2d0*alf1*(bilg-1d0)
* -------------yfs formula ------------------------------------------
      IF(keyd .GE. 300 .AND. keyd .LE. 302) THEN
         delb   = alf1*(0.5d0*bilg -1d0  +pi**2/3d0)
         delb   = delb -betf/2 *dlog(1-uu)
         gamfac = exp(-ceuler*betf)/Semalib_Gamma(1d0+betf)
         ffact  = gamfac*exp(delb)
* ....zero   order
         IF(keyd   .EQ. 300)  THEN
            dels = 0d0
            delh = -betf/4d0 *log(1d0-uu)
* ....first  order
         ELSEIF(keyd   .EQ. 301)  THEN
            dels = betf/2d0
***            delh = uu*(-1d0 +uu/2d0)
***     $       +betf*(-uu**2/2d0 - 0.25d0*(1-uu)**2*log(1-uu) )
**** final state specific contr.
***     $       +betf*uu/2d0*(1-uu/2)*log(1-uu)
            delh = uu*(-1d0 +uu/2d0)
     $          +betf*( -0.5d0*uu**2
     $               +0.25*(-1d0+4d0*uu-2d0*uu**2)*log(1d0-uu) )
* ....second order
         ELSEIF(keyd   .EQ. 302)  THEN
            dels  = betf/2d0 +betf**2/8d0
            delh  = uu*(-1d0+uu/2d0)
     $        +betf*(-0.5d0*uu-0.25d0*uu*(-1d0+0.5d0*uu)*log(1d0-uu))
         ENDIF
         soft = ffact*(1d0+dels)
         hard = ffact*betf*uu**(betf-1d0)*delh
         distr= betf*uu**(betf-1d0)* soft +hard
* -------------------------------------------------------------------
* -------------contributions  from various beta's -------------------
      ELSEIF(keyd .GE. 310 .AND. keyd .LE. 322)  THEN
         delb   = alf1*(0.5d0*bilg -1d0  +pi**2/3d0)
         delb   = delb -betf/2 *dlog(1-uu)
         gamfac = exp(-ceuler*betf)/Semalib_Gamma(1d0+betf)
         ffact  = gamfac*exp(delb)
         soft   = 0d0
         delh   = 0d0
* ...beta0 first  order
         IF(keyd .EQ. 310) THEN
            soft = 1d0 + betf/2d0
            delh = -betf/4d0 *log(1d0-uu)
* ...beta1 first  order
         ELSEIF(keyd .EQ.  311)  THEN
            delh = uu*(-1d0 +uu/2d0)
     $       +betf*(-uu**2/2d0 -uu*(-1d0+0.5d0*uu)*log(1d0-uu))
*cccc            delh =
*cccc     $      uu*(-1d0+0.5d0*uu/(1d0+betf))*(1d0-betf*log(1d0-uu))
*cccc     $       -betf*uu/2d0*(-1d0+uu/2)*log(1d0-uu)
* ...beta0 second order
         ELSEIF(keyd .EQ. 320) THEN
            soft = 1d0 + betf/2d0  +betf**2/8d0
            delh = -betf/4d0 *log(1d0-uu)
* ...beta1 second order
         ELSEIF(keyd .EQ.  321)  THEN
*cc            delh =
*cc     $      uu*(1d0+0.5*betf)*(-1d0+0.5d0*uu/(1d0+betf))
*cc     $                     *(1d0-0.5*betf*log(1d0-uu))
*cc     $       +0.25*betf*log(1-uu)*0.5*(1+(1-uu)**2)
*ccc final state specific contrib.
*cc     $       +betf*uu/2*(1-uu/2)*log(1-uu)
            delh = uu*(-1d0+uu/2d0)
     $          +betf*(-uu/2d0 -uu**2/4d0
     $                 +(2d0+6d0*uu-3d0*uu**2)/8d0*log(1d0-uu))

* ...beta2 second order
         ELSEIF(keyd .EQ.  322)  THEN
            delh = betf*(uu**2/4d0 +uu*(-1d0+uu/2d0)/2d0*log(1d0-uu))
         ENDIF
         soft = ffact*soft
         hard = ffact*betf*uu**(betf-1d0)*delh
         distr= betf*uu**(betf-1d0)* soft +hard
      ELSE
         WRITE(m_out,*) ' ===--->  wrong KeyFin in Semalib_uurho'
         STOP
      ENDIF
      END


      DOUBLE PRECISION FUNCTION Semalib_Born(svar)
*///////////////////////////////////////////////////////////////////////////
*//                                                                       //
*///////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      DOUBLE PRECISION    svar
      DOUBLE PRECISION    pi
      PARAMETER(pi=3.1415926535897932d0)
*---
      INCLUDE "semalib.h"
*---
      DOUBLE PRECISION    gnanob,alfa,alfinv,Semalib_Borny
*---------------------------------------------------------------
      alfinv = m_xpar(30)
      gnanob = m_xpar(31)
      alfa=1/alfinv
      Semalib_Born = 4*pi*alfa**2/(3d0*svar )*Semalib_Borny(svar ) *gnanob
      END

      DOUBLE PRECISION FUNCTION Semalib_Borny(svar1)
*///////////////////////////////////////////////////////////////////////////
*//                                                                       //
*// this routine calculates total Born cross section.                     //
*// it exploits the fact that Born x. section = a + b*c + d*c**2          //
*//                                                                       //
*///////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      DOUBLE PRECISION    svar1
      DOUBLE PRECISION    svar,Semalib_Bornv
*----------------------------------------------------------------
      svar=svar1
      Semalib_Borny = Semalib_Bornv(svar, 0.d0)
      END

      DOUBLE PRECISION FUNCTION Semalib_Bornv(svar,costhe)
*///////////////////////////////////////////////////////////////////////////
*//                                                                       //
*// Similar to BornV_Simple from yfs3.f                                   //
*//                                                                       //
*// This routine provides unsophisticated Born differential cross section //
*// at the crude x-section level, with Z and gamma s-chanel exchange.     //
*// It is ising heavily xpar input inherited from the main program        //
*///////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      DOUBLE PRECISION   svar,costhe
      SAVE
*---
      INCLUDE "semalib.h"
*---
      DOUBLE PRECISION   s,amaz,sinw2,gammz
      DOUBLE PRECISION   sum,T3e,t3f,qf,Qe,deno,Ve,Ae,thresh,born
      DOUBLE PRECISION   xe,yf,xf,ye,rechi,ff0,ff1,amx2,chi2,amfin,vf,af
      INTEGER  KeyZet,HadMin,KFbeam
      INTEGER  i,ke,KFfin,ncf,kf
*--------------------------------------------------------------------
*
      s = svar

      KeyZet = m_KeyZet         ! defined in initialization
      amaz   = m_Zmass          ! defined in initialization

      sinw2  = m_xpar(503)
      gammz  = m_xpar(504)
      HadMin = m_xpar(51)

* Z and gamma couplings to beams (electrons)

      KFbeam = 11           ! KF=11 is electron
      ke = 500+10*KFbeam
      T3e = m_xpar(ke+4)/2d0  ! isospin, L-hand component
      Qe  = m_xpar(ke+3)/3d0  ! electric charge
      deno= 4d0*sqrt(sinw2*(1d0-sinw2))
      Ve= (2*T3e -4*Qe*sinw2)/deno
      Ae=  2*T3e             /deno

* Z and gamma couplings to final fermions
* Loop over all flavours defined in m_xpar(400+i)
      sum = 0d0
      DO i=1,20
         KFfin = m_xpar(400+i)
         Born =0d0
         IF(  ((KFfin .GE.  1) .AND. (KFfin .LE.  6)) .OR.
     $        ((KFfin .GE. 11) .AND. (KFfin .LE. 16))   ) THEN
*
            kf = 500+10*KFfin
            amfin = m_xpar(kf+6)
            NCf   = m_xpar(kf+2)
            T3f   = m_xpar(kf+4)/2d0 ! isospin, L-hand component
            Qf    = m_xpar(kf+3)/3d0 ! electric charge
            deno  = 4d0*sqrt(sinw2*(1d0-sinw2))
            Vf = (2*T3f -4*Qf*sinw2)/deno
            Af =  2*T3f             /deno
            
            IF(KeyZet .LE. 0) THEN
               Ve=0d0
               Ae=0d0
            ENDIF
            IF(KeyZet .EQ. 9) THEN
               Qe=0d0
               Qf=0d0
            ENDIF
c[[[[[[[!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c         Ve     = 2*T3e             /Deno
c         Vf     = 2*T3e             /Deno
c         Ae     =  0d0
c         Af     =  0d0
c]]]]]]]!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

*************
            IF(abs(costhe) .GT. 1d0) 
     $           WRITE(*,*) ' Semalib_Bornv: costhe=',costhe
*************

            chi2 =         s**2/((s-amaz**2)**2+(gammz*s/amaz)**2)
            rechi=(s-amaz**2)*s/((s-amaz**2)**2+(gammz*s/amaz)**2)
**            chi2 =         s**2/((s-amaz**2)**2+(gammz*amaz)**2)
**            rechi=(s-amaz**2)*s/((s-amaz**2)**2+(gammz*amaz)**2)


            xe= Ve**2 +Ae**2
            xf= Vf**2 +Af**2
            ye= 2*Ve*Ae
            yf= 2*Vf*Af
            ff0= qe**2*qf**2 +2*rechi*qe*qf*Ve*Vf +chi2*xe*xf
            ff1=             +2*rechi*qe*qf*Ae*Af +chi2*ye*yf
            Born    = (1d0+ costhe**2)*ff0 +2d0*costhe*ff1
* Colour factor
            Born = NCf*Born
* This is a bit crude method of introducing threshold behaviour
* cos(theta) depencence incorrect!!!
            IF(    svar .LE.  4d0*amfin**2) THEN
               thresh=0d0
            ELSEIF(svar .LE. 16d0*amfin**2) THEN
               amx2=4d0*amfin**2/svar
               thresh=sqrt(1d0-amx2)*(1d0+amx2/2d0)
            ELSE
               thresh=1d0
            ENDIF
            Born= Born*thresh
         ENDIF
* For light quarks u,d,s, special cut on mass (just in case)
         IF( (KFfin .GE. 1) .AND. (KFfin .LE. 3)) THEN
            IF( svar .LE. HadMin**2) Born=0d0
         ENDIF
         sum = sum +Born
      ENDDO
      Semalib_Bornv = sum
      END

      SUBROUTINE Semalib_vvrhoS(KeyIni,svar,amel,vv,vvmin,distr)
*///////////////////////////////////////////////////////////////////
*//                                                               //
*//  THIS IS OLDER COLLECTION OF DISTRIBUTIONS MAINLY FOR SUSSEX  //
*//                                                               //
*///////////////////////////////////////////////////////////////////
* convention for KeyIni
* pedagogical exercises
*     KeyIni   =  1      crude distribution for initial state mc
*     KeyIni   =  9      reference distr.  of yfs2 cpc paper
*     KeyIni   =  50-52  obsolete test distr. for yfs2 cpc paper
*     KeyIni   =  101    soft part yfs       first  order
*     KeyIni   =  102    soft part yfs       second order
*     KeyIni   =  105    hard non-exp.       first  order
*     KeyIni   =  106    hard non-exp.       second order
* total results
*     KeyIni   =  0 + r*100                  zero   order
*     KeyIni   =  1 + r*100                  first  order
*     KeyIni   =  2 + r*100                  second order
*     r = 200 kuraev-fadin
*     r = 400 yfs single electron ll str. funct.
*-------------------------------------------------------------
      IMPLICIT NONE
      INTEGER KeyIni
      DOUBLE PRECISION  svar,amel,vv,vvmin,distr
      SAVE
*---
      DOUBLE PRECISION  pi,ceuler,dz2
      PARAMETER(pi= 3.1415926535897932d0 )
      PARAMETER(ceuler =0.57721566d0)
      PARAMETER(dz2 = pi**2/6d0)
*---
      INCLUDE "semalib.h"
*---
      DOUBLE PRECISION  alfinv,alf1,bilg,beti,dilat,beti2
      DOUBLE PRECISION  dist2,delh2,dist1,soft,delh,delvs 
      DOUBLE PRECISION  gamfac,delb,damel,z,dels
      DOUBLE PRECISION  Semalib_Gamma
      INTEGER keyd
*-------------------------------------------------------------------
      alfinv = m_xpar(30)
      keyd = KeyIni
      alf1   = 1d0/pi/alfinv
      bilg   = dlog(svar/amel**2)
      beti   = 2d0*alf1*(bilg-1d0)
*===================================================================
* ---------------------- keyd = 1 ----------------------------------
* ---- crude distribution in yfs2 initial state monte carlo --------
* ------------------------------------------------------------------
* dilat is related to dilatation jacobian in yfsgen
* damel is responsible for modification of photon ang. distribution
* see also weight wt=wt1 in   angbre
      IF(keyd .GE. 1 .AND. keyd .LT. 100) THEN
         dilat=1d0
         IF(vv .GT. vvmin) dilat=(1d0+1d0/sqrt(1d0-vv))/2d0
         beti2  = 2d0*alf1*bilg
         beti2  = beti2 * m_Xenph !!! photon enhancement factor
         damel=1d0
         IF(vv .GT. vvmin) damel=beti2/beti*(vv/vvmin)**(beti2-beti)
*---------
         IF    (keyd .EQ. 1)  THEN
            distr= beti*vv**(beti-1d0)*dilat*damel
* ...reference distribution used in yfs2 paper --------------------
         ELSEIF(keyd .EQ.  9)  THEN
            distr= beti*vv**(beti-1d0)*(1+(1-vv)**2)/2
* basic reference distribution  xrefer=sigma-ref
         ELSEIF(keyd .EQ. 50) THEN
            distr= beti*vv**(beti-1d0)
* xrefer times damel
         ELSEIF(keyd .EQ. 51) THEN
            distr= beti*vv**(beti-1d0)*damel
* xrefer times dilatation factor dilat
         ELSEIF(keyd .EQ. 52) THEN
            distr= beti*vv**(beti-1d0)*dilat
         ENDIF
* ------------------------------------------------------------
* ------------- soft part only  -- yfs exponentiation
* ------------------------------------------------------------
      ELSEIF(keyd .GE.  101 .AND. keyd .LE. 102) THEN
         delb   = alf1*(0.5d0*bilg -1d0  +pi**2/3d0)
         gamfac =exp(-ceuler*beti)/Semalib_Gamma(1d0+beti)
* ....first order
         IF(keyd   .EQ.  101) THEN
            dels = alf1*(bilg-1d0)
* ....second order
         ELSEIF(keyd   .EQ.  102) THEN
            dels = alf1*(bilg-1d0) +0.5d0*(alf1*bilg)**2
         ENDIF
         distr=  gamfac*exp(delb)*beti*vv**(beti-1d0)*( 1d0+dels)
* ---------------------------------------------------------------
* ------------- hard part only - no exponentiation
* ---------------------------------------------------------------
      ELSEIF(keyd .GE.  105 .AND. keyd .LE. 106 ) THEN
         z=1d0-vv
* ....first  order
         IF(keyd   .EQ.  105) THEN
            distr =  beti*(1d0 +z**2)/2d0/vv
* ....second order
         ELSEIF(keyd   .EQ.  106) THEN
            dist1 =  beti*(1d0 +z**2)/2d0/vv
            delh2 =  (alf1*bilg)**2*(
     $             -(1d0+z*z)/vv*dlog(z)
     $              +(1d0+z)*(0.5d0*dlog(z)-2d0*dlog(vv))
     $              -2.5d0-0.5d0*z  )
            dist2 =  (alf1*bilg)**2*( 3 +4*dlog(vv))/vv +delh2
            distr =  dist1+dist2
         ENDIF
* ---------------------------------------------------------------
* ---------------------------------------------------------------
* ----- exponentiation following kuraiev fadin prescription
* ---------------------------------------------------------------
      ELSEIF(keyd .GE. 201 .AND. keyd .LE. 202) THEN
* first order  from gerrit burgers (polarization book)
         z=1d0-vv
         IF(keyd   .EQ. 201) THEN
            delvs= alf1*(1.5d0*bilg +2d0*dz2-2d0)
            delh=  -alf1*(1d0+z)*(bilg-1d0)
         ELSEIF(keyd   .EQ. 202) THEN
            delvs= alf1*(1.5d0*bilg +2d0*dz2-2d0)
     $            +alf1**2*(9d0/8d0-2d0*dz2)*bilg**2
            delh=  -alf1*(1d0+z)*(bilg-1d0)
     $      +alf1**2*( -(1d0+z*z)/vv     *dlog(z)
     $              +(1d0+z)*(0.5d0*dlog(z)-2d0*dlog(vv))
     $              -2.5d0-0.5d0*z)*bilg**2
         ENDIF
         distr=  beti*vv**(beti-1d0)*( 1d0+delvs) +delh
* -------------------------------------------------------------------
* -------------single fermion ll fragmentation ----------------------
* -------------------yfs formula-------------------------------------
* -------------------------------------------------------------------
      ELSEIF(keyd .GE. 400 .AND. keyd .LE. 402)  THEN
*&&&&&&  delb   = beti/4
         delb   = beti/4  +alf1*(-0.5d0  +pi**2/3d0)
         gamfac = exp(-ceuler*beti)/Semalib_Gamma(1d0+beti)
         soft   = 0d0
         delh   = 0d0
* ...zero   order
         IF(keyd .EQ. 400) THEN
            soft = 1
* ...first  order
         ELSEIF(keyd .EQ. 401) THEN
            soft = 1 + beti/2
            delh =
     $      vv*(-1d0+vv/2/(1+beti))
* ...second order
         ELSEIF(keyd .EQ. 402) THEN
            soft = 1 + beti/2  +beti**2/8
*           delh =
*    $      vv*(1+0.5*beti)*(-1d0+vv/2/(1+beti))
*    $      -0.50*beti*log(1-vv)*0.5*(1+(1-vv)**2)
*    $      + beti/4d0*vv*(vv +(1-vv/2)*dlog(1d0-vv))
            delh =
     $      vv*(-1d0+vv/2)
     $      + beti*(-vv/2 -(1+3*(1-vv)**2)/8*dlog(1d0-vv))
         ENDIF
         distr=  gamfac*exp(delb)*beti*vv**(beti-1d0)*(soft+delh)
* -------------------------------------------------------------------
* -------------single fermion ll fragmentation ----------------------
* -------------contributions  from various beta's -------------------
* -------------------------------------------------------------------
      ELSEIF(keyd .GE. 400 .AND. keyd .LE. 422)  THEN
*&&&&&   delb   = beti/4
         delb   = beti/4  +alf1*(-0.5d0  +pi**2/3d0)
         gamfac = exp(-ceuler*beti)/Semalib_Gamma(1d0+beti)
         soft   = 0d0
         delh   = 0d0
* ...beta0 zero   order
         IF(keyd .EQ. 400) THEN
            soft = 1
* ...beta0 first  order
         ELSEIF(keyd .EQ. 410) THEN
            soft = 1 + beti/2
* ...beta1 first  order
         ELSEIF(keyd .EQ.  411)  THEN
            delh =
     $      vv*(-1d0+vv/2/(1+beti))
* ...beta0 second order
         ELSEIF(keyd .EQ. 420) THEN
            soft = 1 + beti/2  +beti**2/8
* ...beta1 second order
         ELSEIF(keyd .EQ.  421)  THEN
            delh =
     $      vv*(1+0.5*beti)*(-1d0+vv/2/(1+beti))
     $      -0.50*beti*log(1-vv)*0.5*(1+(1-vv)**2)
* ...beta2 second order
         ELSEIF(keyd .EQ.  422)  THEN
            delh = beti/4d0*vv*(vv +(1-vv/2)*dlog(1d0-vv))
         ELSE
            GOTO 900
         ENDIF
         distr=  gamfac*exp(delb)*beti*vv**(beti-1d0)*(soft+delh)
* ----------------------------------------------------------------
* -------------single fermion ll fragmentation -------------------
* -------------i n f i n i t e   o r d e r -----------------------
* ----------------------------------------------------------------
      ELSEIF(keyd .EQ. 502)  THEN
         delb   = beti*0.75d0
         gamfac = exp(-ceuler*beti)/Semalib_Gamma(1d0+beti)
         soft   = 1
         delh =
     $   vv*(-1d0+vv/2)
     $   + beti*(-vv**2/4 -(1+3*(1-vv)**2)/8*dlog(1d0-vv))
         distr=  gamfac*exp(delb)*beti*vv**(beti-1d0)*(soft+delh)
      ELSE
         GOTO 900
      ENDIF
      RETURN
 900  WRITE(6,*) ' ===--->  wrong KeyIni in Semalib_vvrhoS',keyd
      STOP
      END


      DOUBLE PRECISION FUNCTION Semalib_Gamma(z)
*///////////////////////////////////////////////////////////////////////////
*//                                                                       //
*//   Douple Precission Gamma function                                    //
*//                                                                       //
*///////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
* DOUBLE PRECISION gamma FUNCTION
      DOUBLE PRECISION z,z1,x,x1,x2,d1,d2,s1,s2,s3,pi,c(20),const
      INTEGER  k,mm,nn,kk
      SAVE c,pi,const
      DATA c( 1) / 8.3333333333333333333333333332d-02/
      DATA c( 2) /-2.7777777777777777777777777777d-03/
      DATA c( 3) / 7.9365079365079365079365079364d-04/
      DATA c( 4) /-5.9523809523809523809523809523d-04/
      DATA c( 5) / 8.4175084175084175084175084175d-04/
      DATA c( 6) /-1.9175269175269175269175269175d-03/
      DATA c( 7) / 6.4102564102564102564102564102d-03/
      DATA c( 8) /-2.9550653594771241830065359477d-02/
      DATA c( 9) / 1.7964437236883057316493849001d-01/
      DATA c(10) /-1.3924322169059011164274322169d+00/
      DATA c(11) / 1.3402864044168391994478951001d+01/
      DATA c(12) /-1.5684828462600201730636513245d+02/
      DATA c(13) / 2.1931033333333333333333333333d+03/
      DATA c(14) /-3.6108771253724989357173265219d+04/
      DATA c(15) / 6.9147226885131306710839525077d+05/
      DATA c(16) /-1.5238221539407416192283364959d+07/
      DATA c(17) / 3.8290075139141414141414141414d+08/
      DATA c(18) /-1.0882266035784391089015149165d+10/
      DATA c(19) / 3.4732028376500225225225225224d+11/
      DATA c(20) /-1.2369602142269274454251710349d+13/
      DATA pi    / 3.1415926535897932384626433832d+00/
      DATA const / 9.1893853320467274178032973641d-01/
      IF(z .GT. 5.75d 1)          GOTO  6666
      nn = z
      IF (z  -  dble(float(nn)))  3,1,3
    1 IF (z     .LE.     0.d0)   GOTO 6667
      Semalib_Gamma = 1.d0
      IF (z     .LE.     2.d0)   RETURN
      z1 = z
    2 z1 = z1  -  1.d0
      Semalib_Gamma = Semalib_Gamma * z1
      IF (z1  -  2.d0)           61,61,2
    3 IF (dabs(z)  .LT.  1.d-29)  GOTO 60
      IF (z .LT. 0.d0)           GOTO 4
      x  = z
      kk = 1
      GOTO 10
    4 x  = 1.d0  -  z
      kk = 2
   10 x1 = x
      IF (x .GT. 19.d0)         GOTO 13
      d1 = x
   11 x1 = x1  +  1.d0
      IF (x1 .GE. 19.d0)        GOTO 12
      d1 = d1 * x1
      GOTO 11
   12 s3 = -dlog(d1)
      GOTO 14
   13 s3 = 0.d0
   14 d1 = x1 * x1
      s1 = (x1  -  5.d-1) * dlog(x1)  -  x1  +  const
      DO 20                  k=1,20
      s2 = s1  +  c(k)/x1
      IF (dabs(s2  -  s1) .LT. 1.d-28) GOTO 21
      x1 = x1 * d1
   20 s1 = s2
   21 s3 = s3  +  s2
      GOTO (50,22),    kk
   22 d2 = dabs(z  -  nn)
      d1 = d2 * pi
      IF (d1 .LT. 1.d-15)      GOTO 31
   30 x2 =  dlog(pi/dsin(d1))  -  s3
      GOTO 40
   31 x2 = -dlog(d2)
   40 mm = dabs(z)
      IF(x2 .GT. 1.74d2)       GOTO 6666
      Semalib_Gamma = dexp(x2)
      IF (mm  .NE.  (mm/2) * 2)   RETURN
      Semalib_Gamma = -Semalib_Gamma
      RETURN
   50 IF(s3  .GT. 1.74d2)  GOTO 6666
      Semalib_Gamma = dexp(s3)
      RETURN
 6666 WRITE(*,*) "///// Semalib_Gamma ..... argument too large.  /////"
      RETURN
 6667 WRITE(*,*) "///// pgamm ..... argument is at pole.  /////"
      RETURN
   60 Semalib_Gamma = 0.d0
      IF(dabs(z) .LT. 1.d-77)   RETURN
      Semalib_Gamma = 1.d0/z
   61 RETURN
      END


      DOUBLE PRECISION FUNCTION Semalib_DiLog(x)
*///////////////////////////////////////////////////////////////////////////
*//                                                                       //
*// Dilogarithm FUNCTION: dilog(x)=int( -ln(1-z)/z ) , 0 < z < x .        //
*// This is the CERNLIB version.                                          //
*//                                                                       //
*///////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      DOUBLE PRECISION x
      DOUBLE PRECISION z,b,a,y,t,s
*--------------------------------------------------
      z=-1.644934066848226d0
      IF(x  .LT. -1.d0)  GOTO 1
      IF(x  .LE.  0.5d0) GOTO 2
      IF(x  .EQ.  1.d0)  GOTO 3
      IF(x  .LE.  2.d0)  GOTO 4
      z=3.289868133696453d0
    1 t=1.d0/x
      s=-0.5d0
      z=z-0.5d0*dlog(dabs(x))**2
      GOTO 5
    2 t=x
      s=0.5d0
      z=0.d0
      GOTO 5
    3 Semalib_DiLog=1.644934066848226d0
      RETURN
    4 t=1.d0-x
      s=-0.5d0
      z=1.644934066848226d0-dlog(x)*dlog(dabs(t))
    5 y=2.666666666666667d0*t+0.666666666666667d0
      b=      0.000000000000001d0
      a=y*b  +0.000000000000004d0
      b=y*a-b+0.000000000000011d0
      a=y*b-a+0.000000000000037d0
      b=y*a-b+0.000000000000121d0
      a=y*b-a+0.000000000000398d0
      b=y*a-b+0.000000000001312d0
      a=y*b-a+0.000000000004342d0
      b=y*a-b+0.000000000014437d0
      a=y*b-a+0.000000000048274d0
      b=y*a-b+0.000000000162421d0
      a=y*b-a+0.000000000550291d0
      b=y*a-b+0.000000001879117d0
      a=y*b-a+0.000000006474338d0
      b=y*a-b+0.000000022536705d0
      a=y*b-a+0.000000079387055d0
      b=y*a-b+0.000000283575385d0
      a=y*b-a+0.000001029904264d0
      b=y*a-b+0.000003816329463d0
      a=y*b-a+0.000014496300557d0
      b=y*a-b+0.000056817822718d0
      a=y*b-a+0.000232002196094d0
      b=y*a-b+0.001001627496164d0
      a=y*b-a+0.004686361959447d0
      b=y*a-b+0.024879322924228d0
      a=y*b-a+0.166073032927855d0
      a=y*a-b+1.935064300869969d0
      Semalib_DiLog = s*t*(a-b)+z
      END



      SUBROUTINE Semalib_GausJad(fun,aa,bb,eeps,result)
*     ****************************************
*////////////////////////////////////////////////////////////////////////////
*//  Gauss integration by s. jadach, oct. 90.                              //      
*//  This is non-adaptive (!!!!) unoptimized (!!!) integration subprogram  //
*////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      DOUBLE PRECISION aa,bb,eeps,result
      SAVE
*---
      EXTERNAL fun
      DOUBLE PRECISION   fun
      INTEGER  i,k,ndivi,iter,itermx
      DOUBLE PRECISION   xplus,xminu,range,sum16,sum8,erabs,erela,fplus,fminu
      DOUBLE PRECISION   a,b,eps,calk8,calk16,x1,x2,xmidle,delta
*---
      DOUBLE PRECISION wg(12),xx(12)
      DATA wg
     $/0.101228536290376d0, 0.222381034453374d0, 0.313706645877887d0,
     $ 0.362683783378362d0, 0.027152459411754d0, 0.062253523938648d0,
     $ 0.095158511682493d0, 0.124628971255534d0, 0.149595988816577d0,
     $ 0.169156519395003d0, 0.182603415044924d0, 0.189450610455069d0/
      DATA xx
     $/0.960289856497536d0, 0.796666477413627d0, 0.525532409916329d0,
     $ 0.183434642495650d0, 0.989400934991650d0, 0.944575023073233d0,
     $ 0.865631202387832d0, 0.755404408355003d0, 0.617876244402644d0,
     $ 0.458016777657227d0, 0.281603550779259d0, 0.095012509837637d0/
      DATA itermx / 15/
*-----------------------------------------------------------------------------
      eps=abs(eeps)
      a=aa
      b=bb
      ndivi=1
* iteration over subdivisions terminated by precision requirement
      DO 400 iter=1,itermx
      calk8  =0d0
      calk16 =0d0
* sum over delta subintegrals
      DO 200 k = 1,ndivi
      delta = (b-a)/ndivi
      x1    =  a + (k-1)*delta
      x2    =  x1+ delta
      xmidle= 0.5d0*(x2+x1)
      range = 0.5d0*(x2-x1)
      sum8 =0d0
      sum16=0d0
* 8- and 12-point   gauss integration over single delta subinterval
      DO 100 i=1,12
      xplus= xmidle+range*xx(i)
      xminu= xmidle-range*xx(i)
      fplus=fun(xplus)
      fminu=fun(xminu)
      IF(i .LE. 4) THEN
          sum8 =sum8  +(fplus+fminu)*wg(i)/2d0
      ELSE
          sum16=sum16 +(fplus+fminu)*wg(i)/2d0
      ENDIF
  100 CONTINUE
      calk8 = calk8 + sum8 *(x2-x1)
      calk16= calk16+ sum16*(x2-x1)
  200 CONTINUE
      erabs = abs(calk16-calk8)
      erela = 0d0
      IF(calk16.ne.0d0) erela= erabs/abs(calk16)
******WRITE(6,*) 'Semalib_GausJad: calk8,calk16=',iter,calk8,calk16,erela
* PRECISION check to terminate integration
      IF(eeps .GT. 0d0) THEN
        IF(erabs .LT.  eps) GOTO 800
      ELSE
        IF(erela .LT.  eps) GOTO 800
      ENDIF
  400 ndivi=ndivi*2
      WRITE(*,*) ' +++++ Semalib_GausJad:  required precision to high!'
      WRITE(*,*) ' +++++ Semalib_GausJad:  iter,erela=',iter,erela
  800 result= calk16
      END



      SUBROUTINE Semalib_GausJad2(fun,aa,bb,eeps,result)
*     ****************************************
*////////////////////////////////////////////////////////////////////////////
*//  Gauss integration by s. jadach, oct. 90.                              //      
*//  This is non-adaptive (!!!!) unoptimized (!!!) integration subprogram  //
*////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      DOUBLE PRECISION aa,bb,eeps,result
      SAVE
*---
      EXTERNAL fun
      DOUBLE PRECISION   fun
      INTEGER  i,k,ndivi,iter,itermx
      DOUBLE PRECISION   xplus,xminu,range,sum16,sum8,erabs,erela,fplus,fminu
      DOUBLE PRECISION   a,b,eps,calk8,calk16,x1,x2,xmidle,delta
*---
      DOUBLE PRECISION wg(12),xx(12)
      DATA wg
     $/0.101228536290376d0, 0.222381034453374d0, 0.313706645877887d0,
     $ 0.362683783378362d0, 0.027152459411754d0, 0.062253523938648d0,
     $ 0.095158511682493d0, 0.124628971255534d0, 0.149595988816577d0,
     $ 0.169156519395003d0, 0.182603415044924d0, 0.189450610455069d0/
      DATA xx
     $/0.960289856497536d0, 0.796666477413627d0, 0.525532409916329d0,
     $ 0.183434642495650d0, 0.989400934991650d0, 0.944575023073233d0,
     $ 0.865631202387832d0, 0.755404408355003d0, 0.617876244402644d0,
     $ 0.458016777657227d0, 0.281603550779259d0, 0.095012509837637d0/
      DATA itermx / 15/
*-----------------------------------------------------------------------------
      eps=abs(eeps)
      a=aa
      b=bb
      ndivi=1
* iteration over subdivisions terminated by precision requirement
      DO 400 iter=1,itermx
      calk8  =0d0
      calk16 =0d0
* sum over delta subintegrals
      DO 200 k = 1,ndivi
      delta = (b-a)/ndivi
      x1    =  a + (k-1)*delta
      x2    =  x1+ delta
      xmidle= 0.5d0*(x2+x1)
      range = 0.5d0*(x2-x1)
      sum8 =0d0
      sum16=0d0
* 8- and 12-point   gauss integration over single delta subinterval
      DO 100 i=1,12
      xplus= xmidle+range*xx(i)
      xminu= xmidle-range*xx(i)
      fplus=fun(xplus)
      fminu=fun(xminu)
      IF(i .LE. 4) THEN
          sum8 =sum8  +(fplus+fminu)*wg(i)/2d0
      ELSE
          sum16=sum16 +(fplus+fminu)*wg(i)/2d0
      ENDIF
  100 CONTINUE
      calk8 = calk8 + sum8 *(x2-x1)
      calk16= calk16+ sum16*(x2-x1)
  200 CONTINUE
      erabs = abs(calk16-calk8)
      erela = 0d0
      IF(calk16.ne.0d0) erela= erabs/abs(calk16)
******WRITE(6,*) 'Semalib_GausJad2: calk8,calk16=',iter,calk8,calk16,erela
* PRECISION check to terminate integration
      IF(eeps .GT. 0d0) THEN
        IF(erabs .LT.  eps) GOTO 800
      ELSE
        IF(erela .LT.  eps) GOTO 800
      ENDIF
  400 ndivi=ndivi*2
      WRITE(*,*) ' +++++ Semalib_GausJad2:  required precision to high!'
      WRITE(*,*) ' +++++ Semalib_GausJad2:  iter,erela=',iter,erela
  800 result= calk16
      END



      SUBROUTINE Semalib_Gauss(f,a,b,eps,dgauss)
*/////////////////////////////////////////////////////////////////////////////
*//                                                                         //
*//     Adaptive DOUBLE PRECISION Gaussian quadrature from CERNLIB          //
*//                                                                         //
*//     Gauss is set equal to the approximate value of the integral of      //
*//     the function f over the interval (a,b), with accuracy parameter eps //
*//                                                                         //
*/////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      DOUBLE PRECISION  f,a,b,eps,dgauss
      DOUBLE PRECISION  w(12),x(12),aa,bb,c1,c2,u,s8,s16,const
      INTEGER i
      EXTERNAL f
***   LOGICAL mflag,rflag
*
*     REAL*16 version for future use, some compilers dont like it
*
*      DATA w / 0.10122 85362 90376 25915 25313 543d0,
*     $         0.22238 10344 53374 47054 43559 944d0,
*     $         0.31370 66458 77887 28733 79622 020d0,
*     $         0.36268 37833 78361 98296 51504 493d0,
*     $         0.27152 45941 17540 94851 78057 246d-1,
*     $         0.62253 52393 86478 92862 84383 699d-1,
*     $         0.95158 51168 24927 84809 92510 760d-1,
*     $         0.12462 89712 55533 87205 24762 822d0,
*     $         0.14959 59888 16576 73208 15017 305d0,
*     $         0.16915 65193 95002 53818 93120 790d0,
*     $         0.18260 34150 44923 58886 67636 680d0,
*     $         0.18945 06104 55068 49628 53967 232d0/
*
*      DATA x / 0.96028 98564 97536 23168 35608 686d0,
*     $         0.79666 64774 13626 73959 15539 365d0,
*     $         0.52553 24099 16328 98581 77390 492d0,
*     $         0.18343 46424 95649 80493 94761 424d0,
*     $         0.98940 09349 91649 93259 61541 735d0,
*     $         0.94457 50230 73232 57607 79884 155d0,
*     $         0.86563 12023 87831 74388 04678 977d0,
*     $         0.75540 44083 55003 03389 51011 948d0,
*     $         0.61787 62444 02643 74844 66717 640d0,
*     $         0.45801 67776 57227 38634 24194 430d0,
*     $         0.28160 35507 79258 91323 04605 015d0,
*     $         0.95012 50983 76374 40185 31933 543d-1/

      DATA w / 
     $     0.10122 85362 90376 259d0,
     $     0.22238 10344 53374 471d0,
     $     0.31370 66458 77887 287d0,
     $     0.36268 37833 78361 983d0,
     $     0.27152 45941 17540 949d-1,
     $     0.62253 52393 86478 929d-1,
     $     0.95158 51168 24927 848d-1,
     $     0.12462 89712 55533 872d0,
     $     0.14959 59888 16576 732d0,
     $     0.16915 65193 95002 538d0,
     $     0.18260 34150 44923 589d0,
     $     0.18945 06104 55068 496d0/

      DATA x / 
     $     0.96028 98564 97536 232d0,
     $     0.79666 64774 13626 740d0,
     $     0.52553 24099 16328 986d0,
     $     0.18343 46424 95649 805d0,
     $     0.98940 09349 91649 933d0,
     $     0.94457 50230 73232 576d0,
     $     0.86563 12023 87831 744d0,
     $     0.75540 44083 55003 034d0,
     $     0.61787 62444 02643 748d0,
     $     0.45801 67776 57227 386d0,
     $     0.28160 35507 79258 913d0,
     $     0.95012 50983 76374 402d-1/
*
*--------------------------------------------------------------
*  START
*
      dgauss=0.0d0
      IF(b .EQ. a) RETURN
      const=0.005d0/(b-a)
      bb=a
*
*  COMPUTATIONAL LOOP
*
 1    aa=bb
      bb=b
 2    c1=0.5d0*(bb+aa)
      c2=0.5d0*(bb-aa)
      s8=0.0d0
      DO 3 i=1,4
         u=c2*x(i)
         s8=s8+w(i)*(f(c1+u)+f(c1-u))
 3    CONTINUE
      s8=c2*s8
      s16=0.0d0
      DO 4 i=5,12
         u=c2*x(i)
         s16=s16+w(i)*(f(c1+u)+f(c1-u))
 4    CONTINUE
      s16=c2*s16
      IF( abs(s16-s8)  .LE.  eps*(1.d0+abs(s16)) ) GOTO 5
      bb=c1
      IF( 1.d0+ABS(const*c2) .NE. 1.d0) GOTO 2
      dgauss=0.0d0
***      CALL kermtr('d103.1',lgFILE,mflag,rflag)
***      IF(mflag) THEN
***         IF(lgFILE .EQ. 0) THEN
***            WRITE(*,6)
***         ELSE
            WRITE(6,6)
***         ENDIF
***      ENDIF
***      IF(.NOT. rflag) CALL abend
      STOP
***      RETURN
 5    dgauss=dgauss+s16
      IF(bb .NE. b) GOTO 1
      RETURN
*
 6    FORMAT( 4x, 'Function Semalib_Gauss ... too high accuracy required')
      END



      SUBROUTINE Semalib_Gauss2(f,a,b,eps,dgauss)
*/////////////////////////////////////////////////////////////////////////////
*//                                                                         //
*//     Adaptive DOUBLE PRECISION Gaussian quadrature from CERNLIB          //
*//                                                                         //
*//     Gauss2 is set equal to the approximate value of the integral of     //
*//     the function f over the interval (a,b), with accuracy parameter eps //
*//                                                                         //
*/////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      DOUBLE PRECISION  f,a,b,eps,dgauss
      DOUBLE PRECISION  w(12),x(12),aa,bb,c1,c2,u,s8,s16,const
      INTEGER i
      EXTERNAL f
***   LOGICAL mflag,rflag
*
*     REAL*16 version for future use, some compilers dont like it
*
*      DATA w / 0.10122 85362 90376 25915 25313 543d0,
*     $         0.22238 10344 53374 47054 43559 944d0,
*     $         0.31370 66458 77887 28733 79622 020d0,
*     $         0.36268 37833 78361 98296 51504 493d0,
*     $         0.27152 45941 17540 94851 78057 246d-1,
*     $         0.62253 52393 86478 92862 84383 699d-1,
*     $         0.95158 51168 24927 84809 92510 760d-1,
*     $         0.12462 89712 55533 87205 24762 822d0,
*     $         0.14959 59888 16576 73208 15017 305d0,
*     $         0.16915 65193 95002 53818 93120 790d0,
*     $         0.18260 34150 44923 58886 67636 680d0,
*     $         0.18945 06104 55068 49628 53967 232d0/
*
*      DATA x / 0.96028 98564 97536 23168 35608 686d0,
*     $         0.79666 64774 13626 73959 15539 365d0,
*     $         0.52553 24099 16328 98581 77390 492d0,
*     $         0.18343 46424 95649 80493 94761 424d0,
*     $         0.98940 09349 91649 93259 61541 735d0,
*     $         0.94457 50230 73232 57607 79884 155d0,
*     $         0.86563 12023 87831 74388 04678 977d0,
*     $         0.75540 44083 55003 03389 51011 948d0,
*     $         0.61787 62444 02643 74844 66717 640d0,
*     $         0.45801 67776 57227 38634 24194 430d0,
*     $         0.28160 35507 79258 91323 04605 015d0,
*     $         0.95012 50983 76374 40185 31933 543d-1/

      DATA w / 
     $     0.10122 85362 90376 259d0,
     $     0.22238 10344 53374 471d0,
     $     0.31370 66458 77887 287d0,
     $     0.36268 37833 78361 983d0,
     $     0.27152 45941 17540 949d-1,
     $     0.62253 52393 86478 929d-1,
     $     0.95158 51168 24927 848d-1,
     $     0.12462 89712 55533 872d0,
     $     0.14959 59888 16576 732d0,
     $     0.16915 65193 95002 538d0,
     $     0.18260 34150 44923 589d0,
     $     0.18945 06104 55068 496d0/

      DATA x / 
     $     0.96028 98564 97536 232d0,
     $     0.79666 64774 13626 740d0,
     $     0.52553 24099 16328 986d0,
     $     0.18343 46424 95649 805d0,
     $     0.98940 09349 91649 933d0,
     $     0.94457 50230 73232 576d0,
     $     0.86563 12023 87831 744d0,
     $     0.75540 44083 55003 034d0,
     $     0.61787 62444 02643 748d0,
     $     0.45801 67776 57227 386d0,
     $     0.28160 35507 79258 913d0,
     $     0.95012 50983 76374 402d-1/
*
*--------------------------------------------------------------
*  START
*
      dgauss=0.0d0
      IF(b .EQ. a) RETURN
      const=0.005d0/(b-a)
      bb=a
*
*  COMPUTATIONAL LOOP
*
 1    aa=bb
      bb=b
 2    c1=0.5d0*(bb+aa)
      c2=0.5d0*(bb-aa)
      s8=0.0d0
      DO 3 i=1,4
         u=c2*x(i)
         s8=s8+w(i)*(f(c1+u)+f(c1-u))
 3    CONTINUE
      s8=c2*s8
      s16=0.0d0
      DO 4 i=5,12
         u=c2*x(i)
         s16=s16+w(i)*(f(c1+u)+f(c1-u))
 4    CONTINUE
      s16=c2*s16
      IF( abs(s16-s8)  .LE.  eps*(1.d0+abs(s16)) ) GOTO 5
      bb=c1
      IF( 1.d0+ABS(const*c2) .NE. 1.d0) GOTO 2
      dgauss=0.0d0
***      CALL kermtr('d103.1',lgFILE,mflag,rflag)
***      IF(mflag) THEN
***         IF(lgFILE .EQ. 0) THEN
***            WRITE(*,6)
***         ELSE
            WRITE(6,6)
***         ENDIF
***      ENDIF
***      IF(.NOT. rflag) CALL abend
      STOP
***      RETURN
 5    dgauss=dgauss+s16
      IF(bb .NE. b) GOTO 1
      RETURN
*
 6    FORMAT( 4x, 'Function Semalib_Gauss2 ... too high accuracy required')
      END

      SUBROUTINE Semalib_GetBorn(Born)
*////////////////////////////////////////////////////////////////////////
*//                                                                    //
*//   Provides value of Born xsect in nanogarns                        //
*//                                                                    //
*////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      DOUBLE PRECISION Born
*---
      INCLUDE "semalib.h"
      DOUBLE PRECISION cmsene,Semalib_Born
*-----------------------------------------------------------
      cmsene=m_xpar(1)
      Born = Semalib_Born(cmsene**2)
      END

      SUBROUTINE Semalib_GetCMSene(CMSene)
*////////////////////////////////////////////////////////////////////////
*//                                                                    //
*//   Provides CMS energy as recorded in semalib                       //
*//                                                                    //
*////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      DOUBLE PRECISION CMSene
*---
      INCLUDE "semalib.h"
*-----------------------------------------------------------
      CMSene=m_xpar(1)
      END

      SUBROUTINE Semalib_GetZmass(Zmass)
*////////////////////////////////////////////////////////////////////////
*//                                                                    //
*//   Provides Z mass                                                  //
*//                                                                    //
*////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      DOUBLE PRECISION Zmass
*---
      INCLUDE "semalib.h"
*-----------------------------------------------------------
      Zmass = m_Zmass
      END

      SUBROUTINE Semalib_SetZmass(Zmass)
*////////////////////////////////////////////////////////////////////////
*//                                                                    //
*//   (re-)defines Z-mass                                              //
*//                                                                    //
*////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      DOUBLE PRECISION Zmass
*---
      INCLUDE "semalib.h"
*-----------------------------------------------------------
      m_Zmass = Zmass
      END


      SUBROUTINE Semalib_SetDis(Key)
*////////////////////////////////////////////////////////////////////////
*//                                                                    //
*//   Sets type of bremsstralung                                       //
*//                                                                    //
*////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INTEGER Key
*---
      INCLUDE "semalib.h"
*-----------------------------------------------------------
      m_KeyDis = Key
      END


*////////////////////////////////////////////////////////////////////////
*//                                                                    //
*//   End of Class SemaLIB                                             //
*//                                                                    //
*////////////////////////////////////////////////////////////////////////
