C This FILE contains main interface for the test of DOUBLE brem.
C spin amplitudes. It consist of sections (see the name for 
C the search of beginning):
C  (1)  SECTION-INITIALIZATION
C  (2)  SECTION-CONSTRUCTION OF SPIN AMPLITUDE WITH TWO PHOTONS
C  (3)  SECTION-ELA'S-TEST-WITH-ITS-INITIALIZATION
C  (4)  SECTION-PHASE_SPACE, note it is often restricted. Check.
C  (5)  SECTION-RANDOM-NUMBERS
C##########################################################
C##########################################################
C##########################################################
C##########################################################



C  (1)  SECTION-INITIALIZATION
      SUBROUTINE GPS_brem_klu(mode,xkeyii,xkeyif,xkeyfi,xkeyff)
C in this routine keys for ini-ini,ini-fin,fin-ini,fin-fin
C are defined. Routine is CALL in many places of this FILE
C and interbrem.f to transfer these switches.
      INTEGER mode,xkeyii,xkeyif,xkeyfi,xkeyff
      INTEGER      ykeyii,ykeyif,ykeyfi,ykeyff
      SAVE         ykeyii,ykeyif,ykeyfi,ykeyff
*** Only one 1 and three 0's make sense for comparison with Ela
***      DATA   ykeyii,ykeyif,ykeyfi,ykeyff /0,0,0,1/   ! Final 0001
***      DATA   ykeyii,ykeyif,ykeyfi,ykeyff /1,0,0,0/   ! Initial 1000
***      DATA   ykeyii,ykeyif,ykeyfi,ykeyff /0,0,1,0/   ! Interf1 0010
      DATA   ykeyii,ykeyif,ykeyfi,ykeyff /0,1,0,0/   ! Interf2 0100

      IF (mode .EQ. 0) THEN
         ykeyii=xkeyii
         ykeyif=xkeyif
         ykeyfi=xkeyfi
         ykeyff=xkeyff
      ELSE
         xkeyii=ykeyii
         xkeyif=ykeyif
         xkeyfi=ykeyfi
         xkeyff=ykeyff
      ENDIF
      END

      SUBROUTINE GPS_MakeTwoPh
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//   Spin amplitudes O(alf2) with hard photons is calculated                       //
*//   optionally sample is generated and tests are CALLed                           //
*//                                                                                 //
*//   INPUT: exactly two hard photon kinematics                                     //
*//                                                                                 //
*//   OUTPUT:                                                                       //
*//   printouts, COMMON working space is avoided                                    //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE "GPS.h"
*
      INTEGER    KFi,KFf,Nphot
      DOUBLE PRECISION     p1(4),p2(4),p3(4),p4(4),Phot(100,4),PH1(4),PH2(4)
*
      DOUBLE PRECISION     ph(4),PX(4),PP(4),QQ(4)
      INTEGER              i,j,k,l,n,last,loop,loop2,Hel,j1,j2,j3,j4
      DOUBLE PRECISION     ChaIni,ChaFin
      DOUBLE PRECISION     svar,svarX,svarX1,svarQ,svar1,svar2,Ene,betaf
      DOUBLE COMPLEX       GPS_soft,GPS_softb
      DOUBLE COMPLEX       Sini(2,100),Sfin(2,100),sProd,Sactu,Cfact0
      DOUBLE COMPLEX       Amp2Phot(2,2,2,2,2,2),AmpBornII(2,2,2,2),AmpBornIF(2,2,2,2)
      DOUBLE COMPLEX       AmpBornFI(2,2,2,2),AmpBornFF(2,2,2,2)
      DOUBLE PRECISION     Fleps,Massf,Mbeam,m1,m2,m3,m4,mph,Vcut(3)
      DOUBLE PRECISION     XborSum,XboxSum,CrudSum,Exp0Sum,Exp1Sum
      DOUBLE PRECISION     Xborn,Xboxy
      DOUBLE PRECISION     DistCru,BornCru,fLLux
      DOUBLE PRECISION     CrudNorm, ExpoNorm
      DOUBLE PRECISION     RhoCru3
      DOUBLE PRECISION     alfQED,alfpini,alfpfin, alfpmix, Emin, MasPhot
      DOUBLE PRECISION     BVR_SForFac,BVR_TForFac
      DOUBLE PRECISION     Yisr,Yfsr,Yint
      DOUBLE PRECISION     YFSkonIni, YFSkonFin, YFS_IRini, YFS_IRfin, YFS_isr, YFS_fsr
      DOUBLE PRECISION     BornV_GetMass, BornV_GetCharge, BornV_Simple
      DOUBLE PRECISION     Wt0,Wt1
      DOUBLE PRECISION     dummy
      DOUBLE PRECISION     WT,eps,amtau
      INTEGER              iSAVECPU,lll
*-----------------------------------------------
      INTEGER   Icont
      SAVE      Icont
      DATA      Icont /0/
*-----------------------------------------------
      INTEGER xkeyii,xkeyif,xkeyfi,xkeyff
      INTEGER nout
*------------------------------------------------------------------------------------------
      DOUBLE PRECISION    Y_IR, N_IR
      Y_IR=1D0                  ! YES, IR included in     GPS_**Plus
      Y_IR=0D0                  ! No,  IR not included in GPS_**Plus
      N_IR=1D0-Y_IR
*--------------------------------------------------------------------------------------
      nout = 16
      CALL GPS_brem_klu(1,xkeyii,xkeyif,xkeyfi,xkeyff)
      CALL GPS_Initialize
*
      CALL KK2f_GetKFini( KFi)     ! Normaly for beam KFi=11 is electron
      CALL MBrA_GetKF(    KFf)     ! Actual KFcode of the final fermion
      Mbeam  =  BornV_GetMass(   KFi)
      Massf  =  BornV_GetMass(   KFf)
      ChaIni =  BornV_GetCharge( KFi)
      ChaFin =  BornV_GetCharge( KFf)

      Fleps =  1d-100
      m1  = Mbeam
      m2  = Mbeam
      m3  = Massf
      m4  = Massf
      mph = Fleps
      CALL KK2f_GetBeams(    p1,p2)
      CALL KK2f_GetFermions( p3,p4)
      CALL KK2f_GetPhotAll(Nphot,Phot) ! ordered in energy
      CALL KK2f_GetVcut( Vcut)
      CALL KK2f_GetEmin( Emin)
      CALL KK2f_GetMasPhot(MasPhot)
C  ###############################################
C  ###############################################
***   IF(nphot .NE. 2) RETURN
C  ###############################################
C  ###############################################
      amtau=p1(4)+p2(4)
      eps=1d-7
      DO l=1,4
         ph1(l)=PHOT(1,l)
         ph2(l)=PHOT(2,l)
         PP(l)=p3(l)+p4(l)+ph1(l)+ph2(l)
      ENDDO
*///////////////////////////////////////////////////////////
      DO lll=1,6 !< main cancerogenic loop
*///////////////////////////////////////////////////////////
         CALL KINE4C(AMTAU,m3,m4,eps,eps,P3,P4,Ph1,ph2,WT)
         DO l=1,4
            phot(1,l)=ph1(l)
            phot(2,l)=ph2(l)
         ENDDO
*///////////////////////////////////////////////////////////
         CALL GPS_Amp2Print(6,0,'test   $',Amp2Phot,p1,p2,p3,p4,ph1,ph2)
         CALL GPS_SoftPlus2(KFi,KFf,p1,m1,p2,m2,p3,m3,p4,m4,ph1,ph2,mph,Amp2Phot) !
         CALL GPS_Amp2Print(nout,1,'test2   $',Amp2Phot,p1,p2,p3,p4,ph1,ph2) !
*///////////////////////////////////////////////////////////
      ENDDO                     !< main cancerogenic loop
*///////////////////////////////////////////////////////////
      END                       !!!end of GPS_Make!!!



C  (2)  SECTION-CONSTRUCTION OF SPIN AMPLITUDE WITH TWO PHOTONS
      SUBROUTINE GPS_SoftPlus2(KFi,KFf,p1,m1,p2,m2,p3,m3,p4,m4,ph1,ph2,mph,Amp2Phot) !
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//   IR-part of 2-photon amplitudes is enriched here by Orer(alpha)                //
*//   GPS_Hini, GPS_Hfin and genuine order(alpha**2) terms                          //
*//   All photon helicity configurations are calculated                             //
*//                                                                                 //
*//   This routine is used in test but can serve as an example to start             //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE "GPS.h"
*
      INTEGER             KFi,KFf
      DOUBLE PRECISION    PX(4),p1(4),p2(4),p3(4),p4(4),ph1(4),ph2(4),ph(4)
      DOUBLE PRECISION    m1,m2,m3,m4,mph
      DOUBLE COMPLEX      Sactu,Sprod,CNorm,CFact0
      INTEGER             Hel1,Hel2
      DOUBLE COMPLEX      Sini(2,100),Sfin(2,100)
      DOUBLE COMPLEX      Amp2Phot(2,2,2,2,2,2),AmpWork(2,2,2,2)
      DOUBLE COMPLEX      AmpBorn(2,2,2,2)
      DOUBLE PRECISION    BornV_GetMass, BornV_GetCharge,ChaIni,ChaFin
      INTEGER             j,j1,j2,j3,j4,k,Sig,l
      INTEGER             xkeyii,xkeyif,xkeyfi,xkeyff
      DOUBLE COMPLEX      GPS_soft,GPS_softb
      DOUBLE COMPLEX      gI,gF
      DOUBLE PRECISION    svarX,svarQ,PP(4),QQ(4)
      DOUBLE PRECISION    Y_IR, N_IR
*------------------------------------------------------------------------------------------
      Y_IR=1D0                  ! YES, IR included in     GPS_**Plus
      Y_IR=0D0                  ! No,  IR not included in GPS_**Plus
      N_IR=1D0-Y_IR

      DO k=1,4
         PP(k) = p1(k)+p2(k)
         QQ(k) = p3(k)+p4(k)
      ENDDO
      svarX = PP(4)**2 - PP(3)**2 - PP(2)**2 - PP(1)**2
      svarQ = QQ(4)**2 - QQ(3)**2 - QQ(2)**2 - QQ(1)**2
*//////////////////////////////////////////////////////////////////
*//                       S-factors                              //
*//////////////////////////////////////////////////////////////////
      ChaIni =  BornV_GetCharge( KFi)
      ChaFin =  BornV_GetCharge( KFf)
      gI = DCMPLX(ChaIni*m_e_QED)
      gF = DCMPLX(ChaFin*m_e_QED)
      IF( m_KeyArb  .EQ.  0 ) THEN
         Sini(1,1)  =  gI *GPS_soft(  1,ph1,p1,p2)
         Sfin(1,1)  = -gF *GPS_soft(  1,ph1,p3,p4)
         Sini(1,2)  =  gI *GPS_soft(  1,ph2,p1,p2)
         Sfin(1,2)  = -gF *GPS_soft(  1,ph2,p3,p4)
      ELSE
         Sini(1,1)  =  gI *GPS_softb( 1,ph1,p1,m1,p2,m2)
         Sfin(1,1)  = -gF *GPS_softb( 1,ph1,p3,m3,p4,m4)
         Sini(1,2)  =  gI *GPS_softb( 1,ph2,p1,m1,p2,m2)
         Sfin(1,2)  = -gF *GPS_softb( 1,ph2,p3,m3,p4,m4)
      ENDIF
      DO j=1,2
         Sini(2,j) = -DCONJG(Sini(1,j))
         Sfin(2,j) = -DCONJG(Sfin(1,j))
      ENDDO

      CALL GPS_brem_klu(1,xkeyii,xkeyif,xkeyfi,xkeyff)

* remember only one of xkey is non-zero!
      DO k=1,4
         PX(k) = p1(k)+p2(k)-ph2(k)*xkeyii-ph1(k)*xkeyii-ph2(k)*xkeyfi-ph1(k)*xkeyif !
      ENDDO
      CALL GPS_Born(KFi,KFf,PX,p1,m1,p2,-m2,p3,m3,p4,-m4,AmpBorn)

      DO Hel1=1,2
      DO Hel2=1,2
      CALL GPS_BornZero(AmpWork)
      CNorm=DCMPLX(1d0,0d0)
C--->  ph1-ini ph2-ini:--------------------------------------------------------------------------------
      IF (xkeyii .EQ. 1) THEN
         Sprod =  Sini(Hel1,1)*Sini(Hel2,2) *N_IR
         CALL GPS_BornAdd(Sprod, AmpBorn, AmpWork) !
         Sactu= Sini(Hel1,1)
         CALL GPS_HiniAdd(KFi,KFf,PX,p1,m1,p2,m2,p3,m3,p4,m4, Hel1,ph1,mph, Sactu,SProd, AmpBorn,AmpWork) !
         Sactu= Sini(Hel2,2)
         CALL GPS_HiniAdd(KFi,KFf,PX,p1,m1,p2,m2,p3,m3,p4,m4, Hel2,ph2,mph, Sactu,SProd, AmpBorn,AmpWork) !
         CALL GPS_HiiAdd(CNorm, KFi,KFf,PX,p1,m1,p2,m2,p3,m3,p4,m4, Hel1,ph1, Hel2,ph2,mph,      AmpWork) !
      ENDIF
C--->  ph1-ini ph2-fin:--------------------------------------------------------------------------------
      IF (xkeyif .EQ. 1) THEN
         Sprod =  Sini(Hel1,1)*Sfin(Hel2,2) *N_IR
         CALL GPS_BornAdd(Sprod, AmpBorn, AmpWork) !
         Sactu= Sini(Hel1,1)
         CALL GPS_HiniAdd(KFi,KFf,PX,p1,m1,p2,m2,p3,m3,p4,m4, Hel1,ph1,mph, Sactu,SProd, AmpBorn,AmpWork) !
         Sactu= Sfin(Hel2,2)
         CALL GPS_HfinAdd(KFi,KFf,PX,p1,m1,p2,m2,p3,m3,p4,m4, Hel2,ph2,mph, Sactu,SProd, AmpBorn,AmpWork) !
         CALL GPS_HifAdd(CNorm,  KFi,KFf,PX,p1,m1,p2,m2,p3,m3,p4,m4, Hel1,ph1, Hel2,ph2,mph,     AmpWork) !
      ENDIF
C---> ph1-fin  ph2-ini:--------------------------------------------------------------------------------
      IF (xkeyfi .EQ. 1) THEN
         Sprod =  Sfin(Hel1,1)*Sini(Hel2,2) *N_IR
         CALL GPS_BornAdd(Sprod, AmpBorn, AmpWork) !
         Sactu= Sfin(Hel1,1)
         CALL GPS_HfinAdd(KFi,KFf,PX,p1,m1,p2,m2,p3,m3,p4,m4, Hel1,ph1,mph, Sactu,SProd, AmpBorn,AmpWork) !
         Sactu= Sini(Hel2,2)
         CALL GPS_HiniAdd(KFi,KFf,PX,p1,m1,p2,m2,p3,m3,p4,m4, Hel2,ph2,mph, Sactu,SProd, AmpBorn,AmpWork) !
         CALL GPS_HifAdd(CNorm,  KFi,KFf,PX,p1,m1,p2,m2,p3,m3,p4,m4, Hel2,ph2, Hel1,ph1,mph,     AmpWork) !
      ENDIF
C--->  ph1-fin ph2-fin:--------------------------------------------------------------------------------
      IF (xkeyff .EQ. 1) THEN
         Sprod =  Sfin(Hel1,1)*Sfin(Hel2,2) *N_IR
         Cfact0 = sProd
         Cfact0 = sProd  *(svarX/svarQ)
         CALL GPS_BornAdd(Cfact0, AmpBorn, AmpWork) !
         Sactu= Sfin(Hel1,1)
         CALL GPS_HfinAdd(KFi,KFf,PX,p1,m1,p2,m2,p3,m3,p4,m4, Hel1,ph1,mph, Sactu,SProd, AmpBorn,AmpWork) !
         Sactu= Sfin(Hel2,2)
         CALL GPS_HfinAdd(KFi,KFf,PX,p1,m1,p2,m2,p3,m3,p4,m4, Hel2,ph2,mph, Sactu,SProd, AmpBorn,AmpWork) !
         CALL GPS_HffAdd(CNorm,   KFi,KFf,PX,p1,m1,p2,m2,p3,m3,p4,m4, Hel1,ph1, Hel2,ph2,mph,    AmpWork) !
      ENDIF
      DO j1=1,2
         DO j2=1,2
            DO j3=1,2
               DO j4=1,2
                  Amp2Phot(j1,j2,j3,j4,Hel1,Hel2)=AmpWork(j1,j2,j3,j4)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      ENDDO
      ENDDO
      END



      SUBROUTINE GPS_Amp2Print(nout,iprint,word,Amp2Phot,p1,p2,p3,p4,ph1,ph2)
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//   Print 2-photon 64 spin amplitudes in a READible FORMAT on unit nout           //
*//   The above is commented out, but tests are invoked and results printed         //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      DOUBLE PRECISION     p1(4),p2(4),p3(4),p4(4),PH1(4),PH2(4),PP(4),xsect(4)
      DOUBLE COMPLEX       Amp2Phot(2,2,2,2,2,2)
      CHARACTER*8          word
      INTEGER              j1,j2,j3,j4,k,l
      DOUBLE PRECISION     Sum,STOPien,discrep
      DATA init/0/
      DATA STOPien /0d0/
      SAVE init
*
      DO l=1,4
        PP(l)=p3(l)+p4(l)+ph1(l)+ph2(l)
      ENDDO
      DO j1=1,2
         DO j2=1,2
***            WRITE(*,'(a,4(a,4i2,a))')  '     ',
***     $           (('{', 3-2*j1, 3-2*j2, 3-2*j3 , 3-2*j4 ,'}  ', j3=1,2),j4=1,2)
         ENDDO
      ENDDO
      DO l=1,2
         DO k=1,2
            DO j1=1,2
               DO j2=1,2
****               WRITE(nout,'(4(a,2f14.8,a))') 
****               WRITE(nout,'(4(a,2g14.6,a))') 
****     $              (('[',Amp2Phot(j1,j2,j3,j4,k,l),'] ', j3=1,2),j4=1,2)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      Sum=0d0
      DO l=1,2
         DO k=1,2
            DO j1=1,2
               DO j2=1,2
                  DO j3=1,2
                     DO j4=1,2
                        Sum=Sum+ CDABS(Amp2Phot(j1,j2,j3,j4,k,l))**2
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
*     **************************************
      IF (init .EQ. 0) CALL Z0_test          !!!<<<------
*     **************************************
      CALL twoela(p1,p2,p3,p4,ph2,ph1,xsect) !!!<<<------
*     **************************************
      discrep=(sum/xsect(3)-.840913E-02)**2
      ered=4*(p3(4)**2-p3(3)**2-p3(2)**2-p3(1)**2)/
     $     ( (p3(4)+p4(4))**2-(p3(3)+p4(3))**2-(p3(2)+p4(2))**2-(p3(1)+p4(1))**2) !
      ered=SQRT(ered)
      init=init+1
      IF (iprint .EQ. 1) THEN
!!!!      IF (iprint .EQ. 1 .AND. (discrep .GT. STOPien)) THEN
         STOPien=discrep
      IF (IPRINT .EQ. 1) THEN
         WRITE(nout,'(a,i7)') 'evt no= ',init
         WRITE(nout,'(4a)') '+++++++++++++++++++',
     $                     ' 2-phot amplits: ','cofiguration:',
     $                     '+++++++++++++++++++'
         CALL KinLib_VecPrint(nout,'beam1   ',p1)
         CALL KinLib_VecPrint(nout,'beam2   ',p2)
         CALL KinLib_VecPrint(nout,'f1      ',p3)
         CALL KinLib_VecPrint(nout,'f2      ',p4)
         CALL KinLib_VecPrint(nout,'ph1     ',ph1)
         CALL KinLib_VecPrint(nout,'ph2     ',ph2)
         CALL KinLib_VecPrint(nout,'PP      ',PP)
         sgn=-1
         beam1ph1=0d0
         beam1ph2=0d0
         beam2ph1=0d0
         beam2ph2=0d0
         f1ph1=0d0
         f1ph2=0d0
         f2ph1=0d0
         f2ph2=0d0
         DO k=1,4
            IF(k .EQ. 4) sgn=1d0
            beam1ph1=beam1ph1+sgn*(p1(k)-ph1(k))**2
            beam1ph2=beam1ph2+sgn*(p1(k)-ph2(k))**2
            beam2ph1=beam2ph1+sgn*(p2(k)-ph1(k))**2
            beam2ph2=beam2ph2+sgn*(p2(k)-ph2(k))**2
            f1ph1=f1ph1+sgn*(p3(k)+ph1(k))**2
            f1ph2=f1ph2+sgn*(p3(k)+ph2(k))**2
            f2ph1=f2ph1+sgn*(p4(k)+ph1(k))**2
            f2ph2=f2ph2+sgn*(p4(k)+ph2(k))**2
         ENDDO
         WRITE(nout,'(a,g20.12)') 'beam1ph1=',beam1ph1
         WRITE(nout,'(a,g20.12)') 'beam1ph2=',beam1ph2
         WRITE(nout,'(a,g20.12)') 'beam2ph1=',beam2ph1
         WRITE(nout,'(a,g20.12)') 'beam2ph2=',beam2ph2
         WRITE(nout,'(a,g20.12)') 'f1ph1=   ',f1ph1
         WRITE(nout,'(a,g20.12)') 'f1ph2=   ',f1ph2
         WRITE(nout,'(a,g20.12)') 'f2ph1=   ',f2ph1
         WRITE(nout,'(a,g20.12)') 'f2ph2=   ',f2ph2
      ENDIF
      WRITE(nout,'(4a)') '+++++++++++++++++++',
     $                    ' 2-phot amplits: ',word,
     $                   '+++++++++++++++++++'
***      WRITE(nout,'(a,f20.12)') '++++++++++ Sum= ',Sum
      WRITE(nout,'(a,g20.12)') '++++++++++      Sum= ',Sum
      IF (IPRINT .EQ. 555)      
     $WRITE(nout,'(2(a,g20.12,2x))') '++++++++++ xsect(1)= ',xsect(1),'sum/xsect(1)=',sum/xsect(1)
!!!!!      WRITE(nout,'(2(a,g20.12,2x))') '++++++++++ xsect(2)= ',xsect(2),'sum/xsect(2)=',sum/xsect(2)
      IF (IPRINT .EQ. 1)      
     $WRITE(nout,'(4(a,g17.10,2x))') '++++++++++ xsect(3)= ',xsect(3),'sum/xsect(3)=',
     $      sum/xsect(3),'normaliz=',(1-sum/xsect(3)/.840913E-02)/ered,'ered=',ered
!      WRITE(nout,'(3(a,g20.12,2x))') '++++++++++ xsect(4)= ',xsect(4),'sum/xsect(4)=',sum/xsect(4)
      ENDIF
      END
            
*/////////////////////////////////////////////////////////////////////////////////////
C  (3)  SECTION-ELA'S-TEST-WITH-ITS-INITIALIZATION
*/////////////////////////////////////////////////////////////////////////////////////
      SUBROUTINE twoela(p1x,p2x,p3,p4,PH1x,PH2x,xsect)
********************************* SUBROUTINE Z0DEC2(MODE,XPAR,NPAR)
C ***********************************
C main generating routine of TOPIK
C process   q qbar ->l lbar photon photon
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (PI=3.1415926535897932D0, ALFINV=137.03604D0) 
      PARAMETER(GNANOB = 389385D0)
      PARAMETER( ALFPI=  1D0/PI/ALFINV ,ALFA=1D0/ALFINV)
      COMMON / INOUT  / NINP,NOUT,NOUT2
      COMMON /PARFUN/ XMIN,Z1,Z2
      COMMON /BREMSTR/ KEYBRE,KEYPRO
      COMMON /GENER2/  X1,X2,WTMOD
      COMMON /WEIGHTS/ WTINF, WTYFS, WTEXA
      COMMON /MOMLAB2/ P1L(4),Q1L(4),P2L(4),Q2L(4),PHOT1L(4),PHOT2L(4)
      COMMON /MOMCMS2/ P1(4) ,Q1(4) ,P2(4) ,Q2(4) ,PHOT1(4) ,PHOT2(4)
      DIMENSION phot3(4),test(4)
      COMMON /FRAME/ XMSENE,YMSENE,AMINI,AMFIN
      COMMON / WEKINP / AMAZ,GAMMZ,SINW2
      COMMON /FLAVOUR/ IFLEV
      COMMON / RANPAR / KEYRND
      DIMENSION NPAR(99),XPAR(99)
      REAL *4 XPPR(-6:6),ULALPS 
      DIMENSION RN(1),APHOT(4)
      REAL*8     p1x(4),p2x(4),p3(4),p4(4),PH1(4),PH2(4),PP(4),xsect(4),ph1x(4),ph2x(4)
      INTEGER xkeyii,xkeyif,xkeyfi,xkeyff

      KEYRND = 1

      CALL GPS_brem_klu(1,xkeyii,xkeyif,xkeyfi,xkeyff)
      IF (xkeyfi .EQ. 1) THEN
       DO k=1,4
         ph1(k)=ph2x(k)
         ph2(k)=ph1x(k)
       ENDDO
      ELSE
       DO k=1,4
         ph2(k)=ph2x(k)
         ph1(k)=ph1x(k)
       ENDDO
      ENDIF

        svar=(p1x(4)+p2x(4))**2
        XMSENE = DSQRT(SVAR)
C generating tree-body (top,top,photon) phase space and vectors in reduced frame
C this was overritten, we get kinematics from outside
        AMPHOT=0D0
      DO k=1,4
       xsect(k)=1
       p1(k)=p1x(k)
       q1(k)= p2x(k)
       p2(k)=p3(k)
       q2(k)= p4(k)
       phot1(k)=ph1(k)
       phot2(k)=ph2(k)
      ENDDO
      wtkin=1
c.....control integral on the phase space generator
        CALL WMONIT(0,50,WTKIN,1D0,0D0)

        YMSENE=DSQRT( (P2(4)+Q2(4))**2-(P2(3)+Q2(3))**2
     #               -(P2(2)+Q2(2))**2-(P2(1)+Q2(1))**2   )
        zMSENE=DSQRT( (P2(4)+Q2(4)+ph1(4))**2-(P2(3)+Q2(3)+ph1(3))**2
     #               -(P2(2)+Q2(2)+ph1(2))**2-(P2(1)+Q2(1)+ph1(1))**2   )
         DISTRY = 1D0
C   cross section  for q qbar --> Z0 --> l lbar photon
C   according to exact matrix element calculated numeriCALLy with help of spin
C   amplitude technique
C and compact formulas from YFS3
C quark up in initial state
       CALL INTEFB

       XNORM= 1d0 !! (3d0/2d0)**4    !!!!   >>>>>>= ALFA**2*(4D0*PI*ALFA)**2 
       IFLEV=0
       xnorms= 1d0
       IF(KEYBRE .EQ. 1) xnorms=xnorms*(svar/ymsene**2)**2
       IF(KEYBRE .EQ. 2) xnorms=xnorms*(svar/zmsene**2)**2
C tutaj amplitudy w trzech wersjach 
C  DINIINF --podczerwony i bezmasowy born
C  DINIAPR --YFS3
C  DUBLINI --spinowe
       IF(KEYBRE .EQ. 1) THEN
         CALL DINIINF(SECT1)
         XCROS1 =SECT1 *XNORMs !*(2D0/3D0)**4
         CALL DINIAPR(SECT2)
         XCROS2 =SECT2 *XNORM !*(2D0/3D0)**4
         CALL DUBLINI(SECT3)
         XCROS3 =SECT3 *XNORM !*(2D0/3D0)**4
         xsect(1)=xcros1
         xsect(2)=xcros2
         xsect(3)=xcros3
       ELSEIF(KEYBRE .EQ. 3) THEN
         CALL DFININF(SECT11)
         XCROS11=SECT11*XNORM
         CALL DFINAPR(SECT12)
         XCROS12=SECT12*XNORM
         CALL DUBLFIN(SECT13)
         XCROS13=SECT13*XNORM
         xsect(1)=xcros11
         xsect(2)=xcros12
         xsect(3)=xcros13
      ELSEIF(KEYBRE .EQ. 2) THEN
c extra factor *2D0 because the conbinatorial factor 1/2 from
c the phase space shoul be cancel, the are two photons
c but matrix element in principle shoul be averaged over them
c it is not easy because of the generation presampling
c is it O.K.   !!!!!!?????????!!!!!!!!!!!
         CALL DMIXINF(SECT21)
         XCROS21=SECT21*XNORMs !*(2D0/3D0)**2 *2D0
         CALL DMIXAPR(SECT22)
         XCROS22=SECT22*XNORM !*(2D0/3D0)**2 *2D0
         CALL DUBLMIX(SECT23)
         XCROS23=SECT23*XNORM !*(2D0/3D0)**2 *2D0
         xsect(1)=xcros21
         xsect(2)=xcros22
         xsect(3)=xcros23
        
       ENDIF




C   normalization factor (comming from the Born cross section CALLibration
C   with PYTHIA)
         COMFAC=PI/SDOT * (8D0*PI) 
C...matrix element not calculated for WTKIN=0D0
 

 100     CONTINUE


      END




      SUBROUTINE Z0_test
C     *****************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON  /    / BLAN(100000) 
      COMMON / INOUT  / NINP,NOUT,NOUT2
      COMMON / RANPAR / KEYRND
      DIMENSION XPAR(99),NPAR(99)
      INTEGER xkeyii,xkeyif,xkeyfi,xkeyff
      CALL GPS_brem_klu(1,xkeyii,xkeyif,xkeyfi,xkeyff)
      CALL BornV_GetKeyZet(m_KeyZet)
      CALL BornV_GetMZ(    AMZ)
      CALL BornV_GetGammZ( GammZ)
      CALL BornV_GetSwsq(  Sw2 )
      CALL MBrA_GetKF(KFfin)
      amferm  = BornV_GetMass(KFfin)
      WRITE(nout,*) m_Keyzet,amz,gammz,sw2,amferm
      CALL GLIMIT(100000)

      NOUT =46
      NOUT2=45

      NOUT3=99

c....input PARAMETERs>>>
C------option on running PROGRAM :
C               KEYOPT=1  running PYTHIA56
C               KEYOPT=2  running simple MC for t-tbar 
C                       (to test compatibility with PYTHIA56)
C               KEYOPT=3  running option for Z-->l l +  1 photon 
C               KEYOPT=4  running option for Z-->l l +  2 photons 
C               KEYOPT=5  running simple MC for Z--> l l  + PHOTOS_20 final
      KEYOPT = 4
C final state multiplicity.
      KEYDIM  = 2+2
C----- random generator
      KEYRND = 1
C----- bremsstrahlung option
C          NPAR(1)=1 initial state
C          NPAR(1)=2 initial/final state
C          NPAR(1)=3 final state
      NPAR(1)=2
      IF (xkeyii .EQ. 1) NPAR(1)=1
      IF (xkeyff .EQ. 1) NPAR(1)=3
C----- Z0/gamma mixing switched on/off 
C          NPAR(2)=1 Z0/gamma
C          NPAR(2)=2 gamma only
C          NPAR(2)=3 Z0 only
C          NPAR(2)=4 no axial coupling
C          NPAR(2)=5 no vector coupling
C          NPAR(2)=6 only Z0 axial coupling
C          NPAR(2)=7 only Z0 vector coupling
      IF(m_KeyZet .EQ. 1) NPAR(2)=1
      IF(m_KeyZet .EQ. 0) NPAR(2)=2
C------proton-proton CMS energy in GeV
      XPAR(1)= 200D0
C------initial quark  mass in GeV
      XPAR(2)= 0.511D-3  
C------final lepton mass in GeV
      XPAR(3)= amferm  !!.10565830!  1.777 !!    0.511D-3  
c.....trigerr PARAMETERs>>>low level trigger
C------minim photon energy in laboratory  frame
      XPAR(5)= 0.01D0*XPAR(1)/2D0
C------maximum photon energy in laboratory  frame
      XPAR(6)= 0.990D0*XPAR(1)/2D0 !0.99D0*XPAR(1)/2D0
c.....trigerr PARAMETERs>>>physical trigger
c------minimu photon transverse momenta in laboratory frame
      XPAR(7)=   10D0
c------minimu lepton transverse momenta in laboratory frame
      XPAR(8)=   2D0
c------minimum lepton-lepton mass
      XPAR(9)=   2D0
 
c------GSW PARAMETERs
       XPAR(50)  =  amz     !! 91.187D0
       XPAR(51)  =  gammz   !! 2.50072032D0
       XPAR(52)  =  sw2     !! 0.222767773D0


C------number of requested events
      IF(KEYOPT .EQ. 1) NPAR(10)=       100
      IF(KEYOPT .EQ. 2) NPAR(10)=     10 
      IF(KEYOPT .EQ. 3) NPAR(10)=   1 000
      IF(KEYOPT .EQ. 4) NPAR(10)=   40 !      1 000
      IF(KEYOPT .EQ. 5) NPAR(10)=   400 000
      IF(KEYOPT .EQ. 6) NPAR(10)=      1000
 
c-------number of repetation of CALLing PHOTOS from one gevent
      NPAR(11) = 100
C-----output ident.
      IF(KEYOPT .EQ. 1) THEN
         NOUT =36
         NOUT2=35
         NOUTH=30
      ELSEIF(KEYOPT .EQ. 2) THEN
         NOUT =56
         NOUT2=55
         NOUTH=50
      ELSEIF(KEYOPT .EQ. 3) THEN
         NOUT =66
         NOUT2=65
         NOUTH=60
      ELSEIF(KEYOPT .EQ. 4) THEN
         NOUT =76
         NOUT2=75
         NOUTH=70
      ELSEIF(KEYOPT .EQ. 5) THEN
         NOUT =86
         NOUT2=85
         NOUTH=80
      ELSEIF(KEYOPT .EQ. 6) THEN
         NOUT =96
         NOUT2=95
         NOUTH=90
       ENDIF
C------initialize histo output
      CALL GOUTPU(NOUT)
C....initialization of main routines
      IF(KEYOPT .EQ. 1) THEN
c        CALL ROBOL0(-1)   
      ELSEIF(KEYOPT .EQ. 2) THEN
C**        CALL Z0DEC0(-1,XPAR,NPAR)
C**        CALL BOKER1(-1,XPAR,NPAR)
c        CALL BOKER2(-1,XPAR,NPAR)
C**        CALL BOKER3(-1,XPAR,NPAR)
C**      ELSEIF(KEYOPT .EQ. 3) THEN
C**        CALL Z0DEC1(-1,XPAR,NPAR)
C**        CALL BOK1PH(-1,XPAR,NPAR)
C**        CALL BOKER6(-1,XPAR,NPAR)
      ELSEIF(KEYOPT .EQ. 4) THEN
        IF (KEYDIM .EQ. 5) THEN
         CALL Z0DEC3(-1,XPAR,NPAR)
         CALL BOK3PH(-1,XPAR,NPAR)
        ELSE
         CALL Z0DEC2(-1,XPAR,NPAR)
         CALL BOK2PH(-1,XPAR,NPAR)
        ENDIF

c        CALL BOKER7(-1,XPAR,NPAR)
ccc        CALL TRZYNUE(-1,XPAR,NPAR)
C**      ELSEIF(KEYOPT .EQ. 5) THEN
C**        CALL Z0DEC0(-1,XPAR,NPAR)
ccc        CALL BOKER2(-1,XPAR,NPAR)
C**        CALL Z0RAD0(-1,XPAR,NPAR)
C**      ELSEIF(KEYOPT .EQ. 6) THEN
C**        CALL Z0DEC1(-1,XPAR,NPAR)
ccc        CALL BOKER2(-1,XPAR,NPAR)
C**        CALL Z0RAD1(-1,XPAR,NPAR)
      ENDIF
c.... initialization of PYTHIA56
c      CALL PREPYT(XPAR,NPAR)
C
C....generating mode
      DO 10 IEV=1,NPAR(10)
      IF(KEYOPT .EQ. 1) THEN
      IF(MOD(iev,500) .EQ. 1) WRITE(6,*)    'event no=',iev
c         CALL PYEVNT 
c         CALL LUHEPC(1)
c         CALL HEPLUJ
c         CALL ROBOL0(0)
      ELSEIF(KEYOPT .EQ. 2) THEN
      IF(MOD(iev,50 000) .EQ. 1) WRITE(6,*)    'event no=',iev
C**         CALL Z0DEC0(0,XPAR,NPAR)      
C**         CALL BOKER1(0,XPAR,NPAR)
c         CALL BOKER2(0,XPAR,NPAR)
C**         CALL BOKER3(0,XPAR,NPAR)
      ELSEIF(KEYOPT .EQ. 3) THEN
C**      IF(MOD(iev,500) .EQ. 1) WRITE(6,*)    'event no=',iev
C**         CALL Z0DEC1(0,XPAR,NPAR)
C**         CALL BOK1PH(0,XPAR,NPAR)      
C**         CALL BOKER6(0,XPAR,NPAR)
      ELSEIF(KEYOPT .EQ. 4) THEN
      IF(MOD(iev, 100 ) .EQ. 1) WRITE(6,*)    'event no=',iev
        IF (KEYDIM .EQ. 5) THEN
         CALL Z0DEC3(0,XPAR,NPAR)
         CALL BOK3PH(0,XPAR,NPAR)
        ELSE
         CALL Z0DEC2(0,XPAR,NPAR)
         CALL BOK2PH(0,XPAR,NPAR)
        ENDIF

!        CALL BOKER7(0,XPAR,NPAR)
ccc        CALL TRZYNUE(0,XPAR,NPAR)
      ELSEIF(KEYOPT .EQ. 5) THEN
C**      IF(MOD(iev, 100 ) .EQ. 1) WRITE(6,*)    'event no=',iev
C**        CALL Z0DEC0( 0,XPAR,NPAR)
ccc        CALL BOKER2( 0,XPAR,NPAR)
C**        CALL Z0RAD0( 0,XPAR,NPAR)
      ELSEIF(KEYOPT .EQ. 6) THEN
C**      IF(MOD(iev, 100 ) .EQ. 1) WRITE(6,*)    'event no=',iev
C**        CALL Z0DEC1( 0,XPAR,NPAR)
ccc        CALL BOKER2( 0,XPAR,NPAR)
C**        CALL Z0RAD1( 0,XPAR,NPAR)
      ENDIF
   10 CONTINUE
C
C....post generation mode
      IF(KEYOPT .EQ. 1) THEN
c         CALL PYSTAT(1)
c         CALL ROBOL0(1)
      ELSEIF(KEYOPT .EQ. 2) THEN  
C**         CALL Z0DEC0( 1,XPAR,NPAR)
C**         CALL BOKER1( 1,XPAR,NPAR)
c         CALL BOKER2( 1,XPAR,NPAR)
C**         CALL BOKER3( 1,XPAR,NPAR)
      ELSEIF(KEYOPT .EQ. 3) THEN  
C**         CALL Z0DEC1( 1,XPAR,NPAR)
C**         CALL BOK1PH( 1,XPAR,NPAR)
C**         CALL BOKER6( 1,XPAR,NPAR)
      ELSEIF(KEYOPT .EQ. 4) THEN  
        IF (KEYDIM .EQ. 5) THEN
         CALL Z0DEC3(1,XPAR,NPAR)
         CALL BOK3PH(1,XPAR,NPAR)
        ELSE
         CALL Z0DEC2(1,XPAR,NPAR)
         CALL BOK2PH(1,XPAR,NPAR)
        ENDIF
!         CALL BOKER7( 1,XPAR,NPAR)
ccc         CALL TRZYNUE( 1,XPAR,NPAR)
      ELSEIF(KEYOPT .EQ. 5) THEN
C**        CALL Z0DEC0( 1,XPAR,NPAR)
ccc        CALL BOKER2( 1,XPAR,NPAR)
C**        CALL Z0RAD0( 1,XPAR,NPAR)
      ELSEIF(KEYOPT .EQ. 6) THEN
C**        CALL Z0DEC1( 1,XPAR,NPAR)
ccc        CALL BOKER2( 1,XPAR,NPAR)
C**       CALL Z0RAD1( 1,XPAR,NPAR)
      ENDIF
C ------------WRITING HISTOS ON THE DISK ------------------------
      CALL GRFILE(NOUTH,' ','N')
      CALL GROUT( 0,ICY,' ')
      CALL GREND(DNAME)
C ------------THE END OF HISTO WRITING -------------------------
 
CC>>>>>>>>>>>>>>
      END

C  (4)  SECTION-PHASE_SPACE
       SUBROUTINE 
     #     KINE4C(AMTAU,AMP3,AMP2,AMP1,AMNUTA,PIM3,PIM2,PIM1,PN,WT)
C     *******************************************************************
C generator of 4 final state momenta in CMS system
C         AMTAU  -  energy of CMS system
C         AMP1,AMP2,AMP3,AMNUTA  - masses of particles
C         PIM1,PIM2,PIM3,PN  - generated momenta 
C         presampling on  infrared singularity
C         for PN, PIM1 momenta
C    factor 1/2 for two identical particles included
C         WT  - weight
C SUBROUTINE is based on SUBROUTINE DPHTRE from TAUOLA
C     *******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (PI=3.1415926535897932D0, ALFINV=137.03604D0) 
      DIMENSION PIM1(4),PIM2(4),PIM3(4),PN(4),PR(4),PAA(4),PT(4)
      DIMENSION RRR(8)
      COMMON /POMOC2/ Y,Z,DELZ,DEL1,DEL2,BETA,C,SPHI,C1,S1,CJAC

C
C FOUR BODY PHASE SPACE NORMALISED AS IN BJORKEN-DRELL
C D**3 P /2E/(2PI)**3 (2PI)**4 DELTA4(SUM P)
      PHSPAC=1.D0/2**17/PI**8
C TAU MOMENTUM
      PT(1)=0.D0
      PT(2)=0.D0
      PT(3)=0.D0
      PT(4)=AMTAU
C
      CALL VARRAN(RRR,8)
C>>> DO testow na konfiguracje podczerwone
C>>> 2 soft photons
!      RRR(1)=0.99875D0  !ph2=~(1-rrr(1)
!      RRR(2)=0.9975D0  ! final i.e. drugi foton
       rrr(1)=0.1+.9*rrr(1)
       rrr(2)=0.1+.9*rrr(2)
c>>> fix taus

!      rrr(3)=0.99999
!      rrr(4)=0.5
!       rrr(5)=1-0.99999
!      rrr(6)=0.5
c SUBROUTINE spherd
!      rrr(7)=0.005
!      rrr(8)=0.7

C>>> 1 soft photon
C>>>          RRR(1)=0.99999D0
C MASS OF (REAL/VIRTUAL) A1 -> ( RHO + PIM1) flat phase space
           AMS1=(AMP1+AMP2+AMP3)**2
           AMS2=(AMTAU-AMNUTA)**2
           AM3SQ=AMS1+   RRR(1)*(AMS2-AMS1)
           AM3 =SQRT(AM3SQ)
           PHSPAC = PHSPAC*(AMS2-AMS1)
C PHASE SPACE WITH INFRARED SINGULARITY FOR PN (y generated with density 1/y)
C        AMS1=(AMP1+AMP2+AMP3)**2       
C        Y1  = (AMTAU**2-AMS1+AMNUTA**2)/AMTAU**2
C        Y0  = 0.001D0
C        IF(Y1 .LT. Y0) Y1=Y0
C        YL1 = DLOG(Y1/Y0)
C        YL0 = DLOG(Y0)
C        Y   = DEXP(YL1*RRR(1)+YL0)
C        AM3SQ=(1D0-Y)*AMTAU**2+AMNUTA**2
C        AM3=SQRT(AM3SQ)
C        PHSPAC=PHSPAC*AMTAU**2*YL1*Y
        IF(PHSPAC .EQ. 0D0) GOTO 900
C MASS OF (REAL/VIRTUAL) RHO -> (PIM2+PIM3) flat phase space
        AMS1=(AMP2+AMP3)**2
        AMS2=(AM3-AMP1)**2
        AM2SQ=AMS1+   RRR(2)*(AMS2-AMS1)
        AM2 =SQRT(AM2SQ)
        PHSPAC=PHSPAC*(AMS2-AMS1)
C PHASE SPACE WITH INFRARED SINGULARITY FOR PIM1 (y generated with density 1/y)
!        AMS1=(AMP2+AMP3)**2
!        Z1  = (AM3**2-AMS1+AMP1**2)/AM3**2
!        Z0  = 0.001D0
!        IF(Z1 .LT. Z0) Z1=Z0
!        ZL1 = DLOG(Z1/Z0)
!        ZL0 = DLOG(Z0)
!        Z   = DEXP(ZL1*RRR(2)+ZL0)
!        AM2SQ=(1D0-Z)*AM3**2+AMP1**2
!        AM2=SQRT(AM2SQ)
!        PHSPAC=PHSPAC*AM3**2*ZL1*Z
        IF(PHSPAC .EQ. 0D0) GOTO 900
* RHO RESTFRAME, DEFINE PIPL AND PIM1
        ENQ1=(AM2SQ-AMP2**2+AMP3**2)/(2*AM2)
        ENQ2=(AM2SQ+AMP2**2-AMP3**2)/(2*AM2)
        PPI=         ENQ1**2-AMP3**2
        PPPI=SQRT(ABS(ENQ1**2-AMP3**2))
        PHSPAC=PHSPAC*(4*PI)*(2*PPPI/AM2)
* PI MINUS MOMENTUM IN RHO REST FRAME
        CALL SPHERDx(PPPI,PIM3,rrr(7),rrr(8))
        PIM3(4)=ENQ1
* PI0 1 MOMENTUM IN RHO REST FRAME
        DO 30 I=1,3
 30     PIM2(I)=-PIM3(I)
        PIM2(4)= ENQ2
* A1 REST FRAME, DEFINE PIM1
*       RHO  MOMENTUM
        PR(1)=0.D0
        PR(2)=0.D0
        PR(4)=1.D0/(2*AM3)*(AM3**2+AM2**2-AMP1**2)
        PR(3)= SQRT(ABS(PR(4)**2-AM2**2))
        PPI  =          PR(4)**2-AM2**2
*       PI0 2 MOMENTUM
        PIM1(1)=0.D0
        PIM1(2)=0.D0
        PIM1(4)=1.D0/(2*AM3)*(AM3**2-AM2**2+AMP1**2)
        PIM1(3)=-PR(3)
      PHSPAC=PHSPAC*(4*PI)*(2*PR(3)/AM3)
* OLD PIONS BOOSTED FROM RHO REST FRAME TO A1 REST FRAME
      EXE=(PR(4)+PR(3))/AM2
      CALL BOSTD3(EXE,PIM2,PIM2)
      CALL BOSTD3(EXE,PIM3,PIM3)
* ALL PIONS AND RHO ROTATED IN THE A1 REST FRAME
      THET =ACOS(-1.D0+2*RRR(3))
      PHI = 2*PI*RRR(4)
      CALL ROTPOD(THET,PHI,PIM1)
      CALL ROTPOD(THET,PHI,PIM2)
      CALL ROTPOD(THET,PHI,PIM3)
      CALL ROTPOD(THET,PHI,PR)
C
* NOW TO THE TAU REST FRAME, DEFINE A1 AND NEUTRINO MOMENTA
* A1  MOMENTUM
      PAA(1)=0.D0
      PAA(2)=0.D0
      PAA(4)=1.D0/(2*AMTAU)*(AMTAU**2-AMNUTA**2+AM3**2)
      PAA(3)= SQRT(ABS(PAA(4)**2-AM3**2))
      PPI   =          PAA(4)**2-AM3**2
      PHSPAC=PHSPAC*(4*PI)*(2*PAA(3)/AMTAU)
* TAU-NEUTRINO MOMENTUM
      PN(1)=0.D0
      PN(2)=0.D0
      PN(4)=1.D0/(2*AMTAU)*(AMTAU**2+AMNUTA**2-AM3**2)
      PN(3)=-PAA(3)
* ALL PIONS BOOSTED FROM A1  REST FRAME TO TAU REST FRAME
* Z-AXIS ANTIPARALLEL TO NEUTRINO MOMENTUM
      EXE=(PAA(4)+PAA(3))/AM3
      CALL BOSTD3(EXE,PIM1,PIM1)
      CALL BOSTD3(EXE,PIM2,PIM2)
      CALL BOSTD3(EXE,PIM3,PIM3)
      CALL BOSTD3(EXE,PR,PR)
* ALL PIONS AND NEUTRINO ROTATED IN THE TAU REST FRAME
      THET =ACOS(-1.D0+2*RRR(5))
      PHI = 2*PI*RRR(6)
      CALL ROTPOD(THET,PHI,PIM1)
      CALL ROTPOD(THET,PHI,PIM2)
      CALL ROTPOD(THET,PHI,PIM3)
      CALL ROTPOD(THET,PHI,PN)
C THE STATISTICAL FACTOR FOR IDENTICAL PI'S 
        PHSPAC=PHSPAC/2.D0
C FINAL WEIGHT
      WT = PHSPAC
      RETURN
 900  WT=0D0

      END

      SUBROUTINE SPHERDx(R,X,r1,r2)
C ----------------------------------------------------------------------
C GENERATES UNIFORMLY THREE-VECTOR X ON SPHERE  OF RADIUS R
C DOUBLE PRECISON VERSION OF SPHERA
C ----------------------------------------------------------------------
      REAL*8  R,X(4),PI,COSTH,SINTH
      REAL*8 RRR(2),r1,r2
      DATA PI /3.141592653589793238462643D0/
C
C>>>>>>>>>>>>>      CALL VARRAN(RRR,2)
      rrr(1)=r1  !1-0.005
      rrr(2)=r2  !0.7
      COSTH=-1+2*RRR(1)
      SINTH=SQRT(1 -COSTH**2)
      X(1)=R*SINTH*COS(2*PI*RRR(2))
      X(2)=R*SINTH*SIN(2*PI*RRR(2))
      X(3)=R*COSTH
      RETURN
      END


      SUBROUTINE ROTPOD(THET,PHI,PP)
C ----------------------------------------------------------------------
C
C     rotation on the sphere
C ----------------------------------------------------------------------
      REAL *8  PP(4),thet,phi
C
      CALL ROTOD2(THET,PP,PP)
      CALL ROTOD3( PHI,PP,PP)
      RETURN
      END

C  (5)  SECTION-RANDOM-NUMBERS

      SUBROUTINE VARRAN(DRVEC,LEN)
C     ***************************
C Switchable random number generator
C Translation to DOUBLE PRECISION
C     ***************************
      COMMON / RANPAR / KEYRND
      DOUBLE PRECISION DRVEC(*)
      DIMENSION RVEC(1000)

      KEYRND = 1

      IF(LEN .LT. 1 .OR. LEN .GT. 1000) GOTO 901
      KEYRND=1
   10 CONTINUE
      IF(KEYRND .EQ. 1) THEN
         CALL MARRAN(RVEC,LEN)
      ELSEIF(KEYRND .EQ. 2) THEN
         CALL RANECU(RVEC,LEN)
      ELSE
         GOTO 902
      ENDIF
C random numbers 0 and 1 not accepted
      DO 30 I=1,LEN
      IF(RVEC(I) .LE. 0E0 .OR. RVEC(I) .GE. 1E0) THEN
        WRITE(6,*) ' +++++ VARRAN: RVEC=',RVEC(I)
        GOTO 10
      ENDIF
      DRVEC(I)=RVEC(I)
   30 CONTINUE
      RETURN
  901 WRITE(6,*) ' +++++ STOP IN VARRAN: LEN=',LEN
      STOP
  902 WRITE(6,*) ' +++++ STOP IN VARRAN: WRONG KEYRND',KEYRND
      STOP
      END
      SUBROUTINE MARRAN(RVEC,LENV)
C =======================S. JADACH===================================
C == This commes from F. James, The name of RANMAR is changed to   ==
C == MARRAN in order to avoid interference with the version        ==
C == alREADy in use and the public library version (if present).   ==
C ==      THIS IS THE ONLY MODIFICATION !!!!                       ==
C ========================S. JADACH==================================
C Universal random number generator proposed by Marsaglia and Zaman
C in report FSU-SCRI-87-50
C        modified by F. James, 1988 and 1989, to generate a vector
C        of pseudorandom numbers RVEC of length LENV, and to put in
C        the COMMON block everything needed to specify currrent state,
C        and to add input and output entry points RMARIN, RMARUT.
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C!!!  CALLing sequences for RANMAR:                                  ++
C!!!      CALL RANMAR (RVEC, LEN)   RETURNs a vector RVEC of LEN     ++
C!!!                   32-bit random floating point numbers between  ++
C!!!                   zero and one.                                 ++
C!!!      CALL RMARIN(I1,N1,N2)   initializes the generator from one ++
C!!!                   32-bit INTEGER I1, and number counts N1,N2    ++
C!!!                  (for initializing, set N1=N2=0, but to restart ++
C!!!                    a previously generated sequence, use values  ++
C!!!                    output by RMARUT)                            ++
C!!!      CALL RMARUT(I1,N1,N2)   outputs the value of the original  ++
C!!!                  seed and the two number counts, to be used     ++
C!!!                  for restarting by initializing to I1 and       ++
C!!!                  skipping N2*100000000+N1 numbers.              ++
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DIMENSION RVEC(*)
      COMMON/RASET1/U(97),C,I97,J97
      PARAMETER (MODCNS=1000000000)
      SAVE CD, CM, TWOM24, NTOT, NTOT2, IJKL
      DATA NTOT,NTOT2,IJKL/-1,0,0/
C
      IF (NTOT  .GE.  0)  GO TO 50
C
C        Default initialization. User has CALLed RANMAR without RMARIN.
      IJKL = 54217137
      NTOT = 0
      NTOT2 = 0
      KALLED = 0
      GO TO 1
C
      ENTRY      RMARIN(IJKLIN, NTOTIN,NTOT2N)
C         Initializing routine for RANMAR, may be CALLed before
C         generating pseudorandom numbers with RANMAR. The input
C         values should be in the ranges:  0<=IJKLIN<=900 OOO OOO
C                                          0<=NTOTIN<=999 999 999
C                                          0<=NTOT2N<<999 999 999!
C To get the standard values in Marsaglia's paper, IJKLIN=54217137
C                                            NTOTIN,NTOT2N=0
      IJKL = IJKLIN
      NTOT = MAX(NTOTIN,0)
      NTOT2= MAX(NTOT2N,0)
      KALLED = 1
C          always come here to initialize
    1 CONTINUE
      IJ = IJKL/30082
      KL = IJKL - 30082*IJ
      I = MOD(IJ/177, 177) + 2
      J = MOD(IJ, 177)     + 2
      K = MOD(KL/169, 178) + 1
      L = MOD(KL, 169)
      WRITE(6,'(A,5I10)')
     $ ' MARran INITIALIZED: IJ,KL,IJKL,NTOT,NTOT2=',IJ,KL,IJKL,NTOT,NTO
CCC      PRINT '(A,4I10)', '   I,J,K,L= ',I,J,K,L
      DO 2 II= 1, 97
      S = 0.
      T = .5
      DO 3 JJ= 1, 24
         M = MOD(MOD(I*J,179)*K, 179)
         I = J
         J = K
         K = M
         L = MOD(53*L+1, 169)
         IF (MOD(L*M,64)  .GE.  32)  S = S+T
    3    T = 0.5*T
    2 U(II) = S
      TWOM24 = 1.0
      DO 4 I24= 1, 24
    4 TWOM24 = 0.5*TWOM24
      C  =   362436.*TWOM24
      CD =  7654321.*TWOM24
      CM = 16777213.*TWOM24
      I97 = 97
      J97 = 33
C       Complete initialization by skipping
C            (NTOT2*MODCNS + NTOT) random numbers
      DO 45 LOOP2= 1, NTOT2+1
      NOW = MODCNS
      IF (LOOP2  .EQ.  NTOT2+1)  NOW=NTOT
      IF (NOW  .GT.  0)  THEN
        WRITE(6,'(A,I15)') ' RMARIN SKIPPING OVER ',NOW
       DO 40 IDUM = 1, NTOT
       UNI = U(I97)-U(J97)
       IF (UNI  .LT.  0.)  UNI=UNI+1.
       U(I97) = UNI
       I97 = I97-1
       IF (I97  .EQ.  0)  I97=97
       J97 = J97-1
       IF (J97  .EQ.  0)  J97=97
       C = C - CD
       IF (C  .LT.  0.)  C=C+CM
   40  CONTINUE
      ENDIF
   45 CONTINUE
      IF (KALLED  .EQ.  1)  RETURN
C
C          Normal entry to generate LENV random numbers
   50 CONTINUE
      DO 100 IVEC= 1, LENV
      UNI = U(I97)-U(J97)
      IF (UNI  .LT.  0.)  UNI=UNI+1.
      U(I97) = UNI
      I97 = I97-1
      IF (I97  .EQ.  0)  I97=97
      J97 = J97-1
      IF (J97  .EQ.  0)  J97=97
      C = C - CD
      IF (C  .LT.  0.)  C=C+CM
      UNI = UNI-C
      IF (UNI  .LT.  0.) UNI=UNI+1.
      RVEC(IVEC) = UNI
C             Replace exact zeros by uniform distr. *2**-24
         IF (UNI  .EQ.  0.)  THEN
         ZUNI = TWOM24*U(2)
C             An exact zero here is very unlikely, but let's be safe.
         IF (ZUNI  .EQ.  0.) ZUNI= TWOM24*TWOM24
         RVEC(IVEC) = ZUNI
         ENDIF
  100 CONTINUE
      NTOT = NTOT + LENV
         IF (NTOT  .GE.  MODCNS)  THEN
         NTOT2 = NTOT2 + 1
         NTOT = NTOT - MODCNS
         ENDIF
      RETURN
C           Entry to output current status
      ENTRY RMARUT(IJKLUT,NTOTUT,NTOT2T)
      IJKLUT = IJKL
      NTOTUT = NTOT
      NTOT2T = NTOT2
      RETURN
      END

      SUBROUTINE RCARRY(RVEC,LENV)
C         Add-and-carry random number generator proposed by
C         Marsaglia and Zaman in SIAM J. Scientific and Statistical
C             Computing, to appear probably 1990.
C         modified with enhanced initialization by F. James, 1990
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C!!!  CALLing sequences for RCARRY:                                  ++
C!!!      CALL RCARRY (RVEC, LEN)   RETURNs a vector RVEC of LEN     ++
C!!!                   32-bit random floating point numbers between  ++
C!!!                   zero and one.                                 ++
C!!!      CALL RCARGO(INT)     initializes the generator from one    ++
C!!!                   32-bit INTEGER INT                            ++
C!!!      CALL RCARIN(IVEC)    restarts the generator from vector    ++
C!!!                   IVEC of 25 32-bit INTEGERs (see RCARUT)       ++
C!!!      CALL RCARUT(IVEC)    outputs the current values of the 25  ++
C!!!                 32-bit INTEGER seeds, to be used for restarting ++
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DIMENSION RVEC(LENV)
      DIMENSION SEEDS(24), ISEEDS(24), ISDEXT(25)
      PARAMETER (TWOP12=4096.)
      PARAMETER (ITWO24=2**24, ICONS=2147483563)
      SAVE NOTYET, I24, J24, CARRY, SEEDS, TWOM24
      LOGICAL NOTYET
      DATA NOTYET/.TRUE./
      DATA I24,J24,CARRY/24,10,0./
C
C              Default Initialization by Multiplicative Congruential
      IF (NOTYET) THEN
         NOTYET = .FALSE.
         JSEED = 314159265
         WRITE(6,'(A,I12)') ' RCARRY DEFAULT INITIALIZATION: ',JSEED
            TWOM24 = 1.
         DO 25 I= 1, 24
            TWOM24 = TWOM24 * 0.5
         K = JSEED/53668
         JSEED = 40014*(JSEED-K*53668) -K*12211
         IF (JSEED  .LT.  0)  JSEED = JSEED+ICONS
         ISEEDS(I) = MOD(JSEED,ITWO24)
   25    CONTINUE
         DO 50 I= 1,24
         SEEDS(I) = REAL(ISEEDS(I))*TWOM24
   50    CONTINUE
         I24 = 24
         J24 = 10
         CARRY = 0.
         IF (SEEDS(24)  .LT.  SEEDS(14)) CARRY = TWOM24
      ENDIF
C
C          The Generator proper: "Subtract-with-borrow",
C          as proposed by Marsaglia and Zaman,
C          Florida State University, March, 1989
C
      DO 100 IVEC= 1, LENV
      UNI = SEEDS(I24) - SEEDS(J24) - CARRY
      IF (UNI  .LT.  0.)  THEN
         UNI = UNI + 1.0
         CARRY = TWOM24
      ELSE
         CARRY = 0.
      ENDIF
      SEEDS(I24) = UNI
      I24 = I24 - 1
      IF (I24  .EQ.  0)  I24 = 24
      J24 = J24 - 1
      IF (J24  .EQ.  0)  J24 = 24
      RVEC(IVEC) = UNI
  100 CONTINUE
      RETURN
C           Entry to input and float INTEGER seeds from previous run
      ENTRY RCARIN(ISDEXT)
         TWOM24 = 1.
         DO 195 I= 1, 24
  195    TWOM24 = TWOM24 * 0.5
      WRITE(6,'(A)') ' FULL INITIALIZATION OF RCARRY WITH 25 INTEGERS:'
      WRITE(6,'(5X,5I12)') ISDEXT
      DO 200 I= 1, 24
      SEEDS(I) = REAL(ISDEXT(I))*TWOM24
  200 CONTINUE
      CARRY = REAL(MOD(ISDEXT(25),10))*TWOM24
      ISD = ISDEXT(25)/10
      I24 = MOD(ISD,100)
      ISD = ISD/100
      J24 = ISD
      RETURN
C                    Entry to ouput seeds as INTEGERs
      ENTRY RCARUT(ISDEXT)
      DO 300 I= 1, 24
         ISDEXT(I) = INT(SEEDS(I)*TWOP12*TWOP12)
  300 CONTINUE
      ICARRY = 0
      IF (CARRY  .GT.  0.)  ICARRY = 1
      ISDEXT(25) = 1000*J24 + 10*I24 + ICARRY
      RETURN
C                    Entry to initialize from one INTEGER
      ENTRY RCARGO(INSEED)
      JSEED = INSEED
      WRITE(6,'(A,I12)') ' RCARRY INITIALIZED FROM SEED ',INSEED
C      TWOM24 = 1.
         DO 325 I= 1, 24
           TWOM24 = TWOM24 * 0.5
         K = JSEED/53668
         JSEED = 40014*(JSEED-K*53668) -K*12211
         IF (JSEED  .LT.  0)  JSEED = JSEED+ICONS
         ISEEDS(I) = MOD(JSEED,ITWO24)
  325    CONTINUE
         DO 350 I= 1,24
         SEEDS(I) = REAL(ISEEDS(I))*TWOM24
  350    CONTINUE
         I24 = 24
         J24 = 10
         CARRY = 0.
         IF (SEEDS(24)  .LT.  SEEDS(14)) CARRY = TWOM24
      RETURN
      END

      SUBROUTINE RANECU(RVEC,LEN)
C         Random number generator given by L'Ecuyer in
C            Comm. ACM Vol 31, p.742, 1988
C            modified by F. James to RETURN a vector of numbers
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C!!!  CALLing sequences for RANECU:                                  ++
C!!!      CALL RANECU (RVEC, LEN)   RETURNs a vector RVEC of LEN     ++
C!!!                   32-bit random floating point numbers between  ++
C!!!                   zero and one.                                 ++
C!!!      CALL RECUIN(I1,I2)    initializes the generator from two   ++
C!!!                   32-bit INTEGERs I1 and I2                     ++
C!!!      CALL RECUUT(I1,I2)    outputs the current values of the    ++
C!!!                   two INTEGER seeds, to be used for restarting  ++
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DIMENSION RVEC(*)
      SAVE ISEED1,ISEED2
      DATA ISEED1,ISEED2 /12345,67890/
C
      DO 100 I= 1, LEN
      K = ISEED1/53668
      ISEED1 = 40014*(ISEED1 - K*53668) - K*12211
      IF (ISEED1  .LT.  0) ISEED1=ISEED1+2147483563
C
      K = ISEED2/52774
      ISEED2 = 40692*(ISEED2 - K*52774) - K* 3791
      IF (ISEED2  .LT.  0) ISEED2=ISEED2+2147483399
C
      IZ = ISEED1 - ISEED2
      IF (IZ  .LT.  1)  IZ = IZ + 2147483562
C
      RVEC(I) = REAL(IZ) * 4.656613E-10
  100 CONTINUE
      RETURN
C
      ENTRY RECUIN(IS1,IS2)
      ISEED1 = IS1
      ISEED2 = IS2
      RETURN
C
      ENTRY RECUUT(IS1,IS2)
      IS1 = ISEED1
      IS2 = ISEED2
      RETURN
      END






