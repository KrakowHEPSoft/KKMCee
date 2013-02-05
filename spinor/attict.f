
      SUBROUTINE GPS_BornOld(KFi,KFf,Qs,p1,m1,p2,m2,p3,m3,p4,m4,AmpBorn)
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//   Older version, codded differently, output the same                            //
*//                                                                                 //
*//   Born spin amplitudes calculated with spinor methods.                          //
*//   Mass of the final fermion kept exactly.                                       //
*//                                                                                 //
*//   Input:                                                                        //
*//   KFi, Kff = beam and final fermion flavour codes (to define charges)           //
*//   Qs       = s-chanel momentum for gamma and Z propagators (not for spinors)    //
*//   pi,mi    are for spinors, not for gamma and Z propagators                     //
*//   p1,m1    =fermion momentum and mass (beam)                                    //
*//   p2,m2    =fermion momentum and mass (beam)                                    //
*//   p3,m3    =fermion momentum and mass final state                               //
*//   p4,m4    =fermion momentum and mass final state                               //
*//                                                                                 //
*//   Output:                                                                       //
*//   AmpBorn   = spin amplitudes                                                   //
*//                                                                                 //
*//   Notes:                                                                        //
*//   Electron mass neglected in spinors, this is why we may use Chisholm!          //
*//   Final fermion mass kept exactly.                                              //
*//   Gamma and Z in s-chanel.                                                      //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE "GPS.h"
*
      INTEGER    KFi,KFf
      REAL*8     Qs(4),p1(4),p2(4),p3(4),p4(4)
      REAL*8     m1,m2,m3,m4
      COMPLEX*16 AmpBorn(2,2,2,2)
*
      REAL*8   T3e,Qe,Ve,Ae
      REAL*8   T3f,Qf,Af,Vf
      INTEGER  KeyZet
      INTEGER  NCf,NCe
      REAL*8   Sw2,MZ,GammZ
      REAL*8   svar
*-----------------------------------------------------------------------------
      COMPLEX*16 PropGam,PropZet
      COMPLEX*16 AmpGam(2,2,2,2),AmpZet(2,2,2,2)
      INTEGER    i,j,k,l
      INTEGER    j1,j2,j3,j4
      INTEGER    Hel1,Hel2,Hel3,Hel4
      COMPLEX*16 s31,s24,s14,s32
      REAL*8     CupGam,CupZet,CupZet1,CupZet2
      COMPLEX*16 HeliFactor
      COMPLEX*16 GPS_iProd1
      COMPLEX*16 GPS_iProd2
      REAL*8     dummy
*-----------------------------------------------------------------------------
      CALL GPS_Initialize
* get charges, izospin, color
      CALL BornV_GetParticle(KFi, dummy, Qe,T3e,NCe)
      CALL BornV_GetParticle(KFf, dummy, Qf,T3f,NCf)

****{{{ set electron momenta exactly massles 
**** DANGEROUS because affects calling program!!!!!
**      p1(3) =  p1(4)
**      p2(3) = -p2(4)
****}}}
* EW parameters
      CALL BornV_GetSwsq(  Sw2 )
      CALL BornV_GetMZ(    MZ   )
      CALL BornV_GetGammZ( GammZ)
* Couplings
      Ve    = (2*T3e -4*Qe*Sw2)/DSQRT(Sw2*(1d0-Sw2))/4d0
      Ae    =  2*T3e           /DSQRT(Sw2*(1d0-Sw2))/4d0
      Vf    = (2*T3f -4*Qf*Sw2)/DSQRT(Sw2*(1d0-Sw2))/4d0
      Af    =  2*T3f           /DSQRT(Sw2*(1d0-Sw2))/4d0
* Possibility to switch off Z
      CALL BornV_GetKeyZet(KeyZet)
      IF(KeyZet .LE. 0) THEN
         Ve=0d0
         Ae=0d0
      ENDIF
*=============================================================
* Propagators
      svar=Qs(4)**2-Qs(3)**2-Qs(2)**2-Qs(1)**2
      IF(svar .LE. (ABS(m3)+ABS(m4))**2 ) RETURN
      PropGam =    DCMPLX(  1d0/svar,  0d0)
      PropZet =    1d0/DCMPLX(svar-MZ**2, GammZ*svar/MZ)
*=============================================================
* Clean
      DO j1 = 1,2
         DO j2 = 1,2
            DO j3 = 1,2
               DO j4 = 1,2
                  AmpGam( j1,j2,j3,j4) = DCMPLX(0d0,0d0)
                  AmpZet( j1,j2,j3,j4) = DCMPLX(0d0,0d0)
                  AmpBorn(j1,j2,j3,j4) = DCMPLX(0d0,0d0)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
*************************************************
*   Spinor, massive, based on Chisholm
      DO j1 = 1,2
         DO j3 = 1,2
            DO j4 = 1,2
               Hel1 = 3-2*j1
               Hel3 = 3-2*j3
               Hel4 = 3-2*j4
               Hel2 = -Hel1     ! exact helicity conservation for beams
               j2  = (3-Hel2)/2
               s31 = GPS_iProd2(  Hel3, p3, m3,   Hel1, p1, m1) ! t
               s24 = GPS_iProd2(  Hel2, p2, m2,   Hel4, p4, m4) ! t
               s14 = GPS_iProd2(  Hel1, p1,-m1,   Hel4, p4, m4) ! u
               s32 = GPS_iProd2(  Hel3, p3, m3,   Hel2, p2,-m2) ! u
               CupGam  = Qe*Qf
               CupZet1 = (Ve +Hel2*Ae)*(Vf +Hel1*Af) ! t
               CupZet2 = (Ve +Hel2*Ae)*(Vf +Hel2*Af) ! u
               AmpGam( j1,j2,j3,j4) = PropGam*(CupGam *s31*s24 +CupGam *s14*s32)
               AmpZet( j1,j2,j3,j4) = PropZet*(CupZet1*s31*s24 +CupZet2*s14*s32)
               AmpBorn(j1,j2,j3,j4) = AmpGam(j1,j2,j3,j4)+AmpZet(j1,j2,j3,j4)
            ENDDO
         ENDDO
      ENDDO
      END                       !!!GPS_BornOld!!!


      SUBROUTINE GPS_1PhotSimple(KFi,KFf,p1,p2,q1,q2,ph,y,z,Dist)
*///////////////////////////////////////////////////////////////////////////////
*//                                                                           //
*//   Simple 1-photon first order                                             //
*//   For tests only                                                          //
*//   obsoleted by direct access to QED3                                      //
*//                                                                           //
*///////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INTEGER   KFi,KFf
      REAL*8    p1(4),p2(4),q1(4),q2(4),ph(4),y,z,Dist
*
      REAL*8    PP(4),QQ(4)
      INTEGER   i,j,k
      REAL*8    cosq1p1,cosq1p2,cosq2p1,cosq2p2
      REAL*8    svar,svar1
      REAL*8    dist10,amel,amfin,delp,delq,sfacj,hfac,gf1,gf2,gi1,gi2
      REAL*8    BornV_GetMass
      REAL*8    andi11,andi12,andi21,andi22
*----------------
      REAL*8    wm0,wm1,wmd,a,b,del
      wm0(del,a,b)= 1d0 -2d0*del -del*(a/b+b/a)
      wmd(del,a,b)= 1d0 + del*(a/b+b/a)*(a**2+b**2)/((1-a)**2+(1-b)**2)
      wm1(del,a,b)= 1d0 - del*(a/b+b/a)*(1d0-a)*(1d0-b)*2d0/((1d0-a)**2+(1d0-b)**2)
*----------------
      DO k=1,4
         PP(k)= p1(k)+p2(k)
         QQ(k)= q1(k)+q2(k)
      ENDDO
      svar  = PP(4)**2-PP(3)**2-PP(2)**2-PP(1)**2
      svar1 = QQ(4)**2-QQ(3)**2-QQ(2)**2-QQ(1)**2
*
      amel   = BornV_GetMass(KFi)
      amfin  = BornV_GetMass(KFf)
      delp  =  amel**2/svar
      delq  = amfin**2/svar1
* case of ISR
      CALL KinLib_ThetaR(QQ,p1,p2,q1,q2,cosq1p1,cosq1p2,cosq2p1,cosq2p2)
      CALL GPS_BornSimple(KFi,KFf, svar1, cosq1p1, andi11)
      CALL GPS_BornSimple(KFi,KFf, svar1, cosq1p2, andi12)
      CALL GPS_BornSimple(KFi,KFf, svar1, cosq2p1, andi21)
      CALL GPS_BornSimple(KFi,KFf, svar1, cosq2p2, andi22)
* BornSimple provides distributions in R-units, i.e. without 1/svar1 factor
      andi11 = andi11*svar/svar1
      andi12 = andi12*svar/svar1
      andi21 = andi21*svar/svar1
      andi22 = andi22*svar/svar1
*
      sfacj  =  2d0/(y*z) *wm0(delp,y,z)
c[[   hfac   =  sfacj     *wmd(delp,y,z) ! as in BHLUMI
      hfac   =  2d0/(y*z) *wm1(delp,y,z) ! Original O(alf1)
      gf1 = 0.5d0
      gf2 = 0.5d0
      gi1 = ((1-y)**2)/2d0
      gi2 = ((1-z)**2)/2d0
      dist10= (  gi1*gf1*andi11   +gi1*gf2*andi12
     $          +gi2*gf1*andi21   +gi2*gf2*andi22)*hfac

      Dist = dist10

c[[[[[[
c      Write(*,*) '!!!!!!!!!!!!!!!!!!  GPS_1PhotSimple !!!!!!!!!!!!!!!!!!'
c      Write(*,'(a,6f15.7)') 'svar,svar1             = ',svar,svar1,1-svar1/svar
c      Write(*,'(a,6g20.14)') 'y,z                    = ',y,z,y*z
c      Write(*,'(a,6f15.7)') 'cq1p1,cq1p2,cq2p1,cq2p2= ',cosq1p1,cosq1p2,cosq2p1,cosq2p2
c      Write(*,'(a,6g20.12)') 'sfacj, dist10 =  ', sfacj, dist10/sfacj
c      Write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
c]]]]]]
      END                       !!!!!!  GPS_1PhotSimple



      SUBROUTINE BornTest
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//   This is just junk, essential parts are mooves somewhere else                  //
*//   nevertheless it could be kept for a couple of months, in in case              //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE

      REAL*8 pi
      PARAMETER( pi =3.1415926535897932d0)
      REAL*8  qf1(4),qf2(4)
      REAL*8  pf1(4),pf2(4)
      REAL*8  rf1(4),rf2(4)
      REAL*8  tf1(4),tf2(4)
      REAL*8  xi(4),xir(4)
      REAL*8  eta(4),etar(4)
      INTEGER nphot,KFfin,KFini
      REAL*8  svar,costhe
      REAL*8  BornV_Simple
      REAL*8  BornV_Differential
      REAL*8  Born1,BornX
      REAL*8  WtMain,WtCrud
      REAL*8  cmsene,beta,eta1,eta2,phi,theta
      REAL*8  BornV_GetMass,amfin
      REAL*8  Rot1(4,4),Rot2(4,4)

      INTEGER icont
      DATA icont /0/

*----------------------------------------------

c      KFini = 11                ! KF=11 is electron
c      CALL MBR_GetKF(KFfin)     ! Actual KFcode of final fermion
c      amfin= BornV_GetMass(KFfin)

c      CALL KarFin_GetFermions( qf1,qf2)
c      CALL KarLud_GetBeams(    pf1,pf2)
c      CALL HepEvt_GetNPhot(    nphot)
c      CALL YFS3ff_GetWt(WtMain,WtCrud)

      IF(icont.LE.10 .AND.  nphot .EQ. 0 .AND. WtCrud.NE.0d0) THEN
c         icont = icont+1
c         svar = (pf1(4) +pf2(4))**2
c         costhe=     (pf1(3)*qf1(3)+pf1(2)*qf1(2)+pf1(1)*qf1(1))
c     $          /SQRT(pf1(3)*pf1(3)+pf1(2)*pf1(2)+pf1(1)*pf1(1))
c     $          /SQRT(qf1(3)*qf1(3)+qf1(2)*qf1(2)+qf1(1)*qf1(1))
c[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
* substitute qf1,qf2
c         theta= 1d-9           ! forward
c         theta= pi/6            ! 30degr !!!!
c         theta= pi/4            ! 45degr
c         theta= pi/2           ! 90degr
c         phi= 0d0           ! x-z plane
c         phi= pi/2d0        ! y-z plane
c         phi= 1/4d0*pi      ! in between !!!!
c         cmsene=sqrt(svar)
c         costhe= COS(theta)
c         CALL KinLib_givpair(cmsene,amfin,amfin,qf1,qf2,beta,eta1,eta2)
c         CALL KinLib_RotEul(theta,phi,qf1,qf1)
c         CALL KinLib_RotEul(theta,phi,qf2,qf2)
c]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
c         WRITE(*,*) 
c     $   '//////////////////////////////////////////////////////////////////////////'
c         WRITE(*,*) 
c     $   '//////////////////////////////////////////////////////////////////////////'
c         WRITE(*,*) 'svar,costhe=',svar,costhe, (1-costhe)/2

cc         Born1 = BornV_Simple(         KFini,KFfin,svar,costhe)
cc         Born1 = BornV_Differential( 1,Kffin, svar,costhe, 0.d0,0.d0, 0.d0,0.d0 )
c         CALL GPS_BornSimple(KFini,KFfin,pf1,pf2,qf1,qf2,Born1)

c(((((
c         CALL GPS_MakeBorn(KFini,KFfin,pf1,pf2,qf1,qf2,BornX)
c         CALL GPS_BornPrint(0)
c         WRITE(*,'(a,3f20.14)') ' GPS_Born,BornSimple= ', BornX,Born1,BornX/Born1
c)))))

c(((((
*//  Born spin amplitudes, spinor method, massles, based on Chisholm          //
c         CALL GPS_MakeBorn1(KFini,KFfin,pf1,pf2,qf1,qf2,BornX)
c         CALL GPS_BornPrint(1)
c         WRITE(*,'(a,3f20.14)') ' GPS_Born1,BornSimple= ', BornX,Born1,BornX/Born1
c)))))
c         CALL KinLib_VecPrint(6,'Qf1=    ',qf1)
c         CALL KinLib_VecPrint(6,'Qf2=    ',qf2)

************************************************************************************
************************************************************************************
************************************************************************************


* Get first rotation matrix from GPS to JacobWick
c         CALL GPS_GetXi(xi,eta)
cc         CALL KinLib_VecPrint(6,'xi=       ',xi)
cc         CALL KinLib_VecPrint(6,'eta=      ',eta)
c         CALL GPS_TraJacobWick(1,qf1,xi, xir)  ! to tau rest frame
c         CALL GPS_TraJacobWick(1,qf1,eta,etar) ! to tau rest frame
c         CALL KinLib_VecPrint(6,'xir1=       ',xir)
c         CALL KinLib_VecPrint(6,'etar1=      ',etar)
c         CALL GPS_GPS(xir,etar,Rot1)
* Get second rotation matrix from GPS to JacobWick
c         CALL GPS_GetXi(xi,eta)
c         CALL GPS_TraJacobWick(1,qf2,xi, xir)  ! to tau rest frame
c         CALL GPS_TraJacobWick(1,qf2,eta,etar) ! to tau rest frame
c         CALL KinLib_VecPrint(6,'xir2=       ',xir)
c         CALL KinLib_VecPrint(6,'etar2=      ',etar)
c         CALL GPS_GPS(xir,etar,Rot2)
c         CALL KinLib_LorPrint(6,'Rot1    ',Rot1)
c         CALL KinLib_LorPrint(6,'Rot2    ',Rot2)
c         CALL GPS_RmatMult(Rot1,Rot2)


*         CALL GPS_RmatPrint(6,'JacobW frame')

c         tf1(1) = 0d0
c         tf1(2) = 0d0
c         tf1(3) = 0d0
c         tf1(4) = 1.777d0
c         CALL GPS_TraJacobWick(1,qf1,qf1, rf1)  ! from LAB to tau rest frame
c         CALL KinLib_VecPrint(6,'rf1=    ',rf1)

c         CALL GPS_TralorDoIt(1,tf1,rf1)
c         CALL GPS_TralorDoIt(2,tf1,rf2)
c         CALL KinLib_VecPrint(6,'rf1=    ',rf1)
c         CALL KinLib_VecPrint(6,'rf2=    ',rf2)

**************************************************
* Principal bootstrap test
c         CALL GPS_GetXi(xi,eta)
c         CALL GPS_TralorUnDo(2,xi,xir)
c         CALL GPS_TralorUnDo(2,eta,etar)
c         CALL KinLib_VecPrint(6,'xir=    ',xir)
c         CALL KinLib_VecPrint(6,'etar=   ',etar)
* Additional reverse xcheck
c         xir(1)= 0d0
c         xir(2)= 0d0
c         xir(3)=-1d0
c         xir(4)= 1d0
c         etar(1) = -1d0
c         etar(2) =  0d0
c         etar(3) =  0d0
c         etar(4) =  0d0
c         CALL KinLib_VecPrint(6,'xir=       ',xir)
c         CALL KinLib_VecPrint(6,'etar=      ',etar)
c         CALL GPS_TralorDoIt(2,xir,xi)
c         CALL GPS_TralorDoIt(2,etar,eta)
c         CALL KinLib_VecPrint(6,'xi=     ',xi)
c         CALL KinLib_VecPrint(6,'eta=    ',eta)
**************************************************

c         rf1(4)= 1.777d0
c         rf1(3)= 0d0
c         rf1(2)= 0d0
c         rf1(1)= 0d0
c         CALL GPS_TralorDoIt(2,rf1,rf1)
c         CALL GPS_TraJacobWick(-1,qf1,rf1,rf1)
c         CALL KinLib_VecPrint(6,'rf1=    ',rf1)

c         CALL GPS_TraJacobWick( 1,qf1,qf1,rf1)
c         CALL KinLib_VecPrint(6,'rf1=    ',rf1)

*************************************************************************************
*//   Clasical, massles, in terms of cos(theta)                               //
c         CALL GPS_MakeBorn2(KFini,KFfin,pf1,pf2,qf1,qf2,BornX)
c         CALL GPS_BornPrint(2)
c         CALL GPS_RmatMake2
c         CALL GPS_RmatPrint(6,'bricolage')
c         WRITE(*,'(a,3f20.14)') 'GPS_Born2,BornSimple= ', BornX,Born1,BornX/Born1

      ENDIF

      END
