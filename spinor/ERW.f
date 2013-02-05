*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//                                                                                 //
*//                         Pseudo-CLASS  ERW                                       //
*//                                                                                 //
*//       Collections of subroutines from E. Richter-Was                            //
*//       Used for testing 1 and 2 photon amplitudes                                //
*//       Only cosmetic corrections are introduced                                  //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////

      SUBROUTINE ERW_initialize
*///////////////////////////////////////////////////////////////////////////////////
*//                                                                               //
*//    Initialization at the very beginning                                       //
*//                                                                               //
*///////////////////////////////////////////////////////////////////////////////////
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / WEKINP / AMAZ,GAMMZ,SINW2,IDE    
      COMMON / COEFF  / V,A    
      COMMON / FOTON  / ARBIT(4)
      REAL*8 PX(4)
*------------------------------------
      INTEGER init
      SAVE    init
      DATA init/0/
*------------------------------------
      IF(init .EQ. 1) RETURN
      init = 1
*------------------------------------

* EW parameters
      CALL BornV_GetSwsq(  SINW2 )
      CALL BornV_GetMZ(    AMAZ   )
      CALL BornV_GetGammZ( GAMMZ)

      QE=-1D0
      AA= 4D0*DSQRT(SINW2*(1D0-SINW2))
      V = (-1D0+4*SINW2)/AA
      A =-1D0/AA

* Possibility to switch off Z
      CALL BornV_GetKeyZet(KeyZet)
      IF(KeyZet .LE. 0) THEN
         V=0d0
         A=0d0
      ENDIF

ccc      WRITE(*,*) 'AMAZ,GAMMZ,SINW2= ',AMAZ,GAMMZ,SINW2
ccc      WRITE(*,*) 'V,A             = ',V,A

*ARBITRARY VECTOR FOR PHOTON POLARIZATION IS FIXED
      PX(1)=0d0
      PX(2)=0d0
      PX(3)=1d0
      PX(4)=1d0
      ARBIT(4)=PX(4)
      ARBIT(1)=PX(1)
      ARBIT(3)=PX(3)*COS(25.)-PX(2)*SIN(25.)
      ARBIT(2)=PX(2)*COS(25.)+PX(3)*SIN(25.)
      CALL KinLib_VecPrint(6,'ARBIT=      ',ARBIT)
      END



      SUBROUTINE ERW_INTEFA
C     *********************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / MOMSET / P1(4),Q1(4),P2(4),Q2(4),PH(4)
      COMMON / BHPAR2 / CMSENE,AMEL
      COMMON / BHPAR1 / CMS,AMMI

cc[[      COMMON /MOMCMS1/ P1(4),Q1(4),P2(4),Q2(4),PHOT1(4)
cc[[      COMMON /FRAME/ XMSENE,YMSENE,AMINI,AMFIN

      SAVE
      INTEGER  KFi,KFf,nphot
*--------------------------------------------
      CALL ERW_initialize

      CALL KarLud_GetBeams(    p1,q1)
      CALL KarFin_GetFermions( p2,q2)

* one ISR photon or FSR photon
      CALL KarLud_GetNphot(nphox)     
      CALL KarFin_GetNphot(nphoy)   
      IF( nphox .EQ. 1 ) THEN
         CALL KarLud_GetPhoton1(   1,ph)
      ELSEIF( nphoy .EQ. 1 ) THEN
         CALL KarFin_GetPhoton1(   1,ph)
      ELSE
         WRITE(*,*) ' +++++ ERW_INTEFA: not a single photon !!!!'
         STOP
      ENDIF
c[[      CALL KinLib_VecPrint(6,'ph=      ',ph)
c[[      CALL KinLib_VecPrint(6,'p1=      ',p1)
c[[      CALL KinLib_VecPrint(6,'q1=      ',q1)
c[[      CALL KinLib_VecPrint(6,'p2=      ',p2)
c[[      CALL KinLib_VecPrint(6,'q2=      ',q2)
*-------------
      KFi = 11                      ! KF=11 is electron
      CALL KarLud_GetKFfin(KFf)     ! Actual KFcode of final fermion
      AMEL =  BornV_GetMass(KFi)
      AMMI =  BornV_GetMass(KFf)
      CMSENE = p1(4)+q1(4)
      CMS = DSQRT((p2(4)+q2(4))**2-(p2(3)+q2(3))**2-(p2(2)+q2(2))**2-(p2(1)+q2(1))**2)
*-------------
c[[
c      WRITE(nout,'(a,5g20.14)') ' CMSENE,CMS= ', CMSENE,CMS,(CMSENE/CMS)**2
c      WRITE(nout,'(a,5g20.14)') ' AMEL,AMMI = ', AMEL,AMMI
c]]
*---
c      AMEL  =AMINI
c      CMSENE=XMSENE
c      AMMI  =AMFIN
c      CMS   =YMSENE
c      DO K=1,4
c        PP1(K)=P1(K)
c        PP2(K)=P2(K)
c        QQ1(K)=Q1(K)
c        QQ2(K)=Q2(K)
c        PH(K) =PHOT1(K)
c      ENDDO

      END



*CROSS SECTION FOR S CHANEL INITIAL SINGL BREMSTR
*FROM KLEISS, Z.PHYS.C33.433-443 (1987)
      SUBROUTINE ERW_SKLINI(SECTION)
*     ***********************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 M1,M2,M3,M4,M9,M10,M11,M12
      COMPLEX*16 VINI,CONS1
      COMPLEX*16 XSECTION,SPLUS,SMINS,HIS

      COMMON / MOMSET / P1(4),Q1(4),P2(4),Q2(4),PH(4)
      COMMON / BHPAR2 / CMSENE,AMEL
      COMMON / BHPAR1 / CMS,AMMI        ! Not Used!!!!
      SAVE
      DATA PI,ALFINV /3.1415926535897932D0, 137.03604D0/
*---------------------------------------------------------------------
      CALL ERW_INTEFA
*     ***************
      ALFA=1D0/ALFINV
      CONS1=DCMPLX(0D0,1D0)
 
      S0=CMSENE**2
      S1= 2D0*(P2(4)*Q2(4)-P2(3)*Q2(3)-P2(2)*Q2(2)-P2(1)*Q2(1))
 
      VINI=SMINS(Q2,P2)/(SPLUS(PH,P1)*SPLUS(Q1,PH))
 
*HELICITY AMPLITUDE FOR CONFIGURATION (P1,Q1,P2,Q2,PH)=(+ + + + +)
      M1=CONS1*SPLUS(P1,Q2)**2    *(HIS(1,1,S1)*VINI)
 
*HELICITY AMPLITUDE FOR CONFIGURATION (P1,Q1,P2,Q2,PH)=(+ + + + -)
      M2=CONS1*SMINS(Q1,P2)**2    *(HIS(1,1,S1)*DCONJG(VINI))
 
*HELICITY AMPLITUDE FOR CONFIGURATION (P1,Q1,P2,Q2,PH)=(+ + - - +)
      M3=-CONS1*SPLUS(P1,P2)**2   *(HIS(1,-1,S1)*VINI)
 
*HELICITY AMPLITUDE FOR CONFIGURATION (P1,Q1,P2,Q2,PH)=(+ + - - -)
      M4=-CONS1*SMINS(Q1,Q2)**2   *(HIS(1,-1,S1)*DCONJG(VINI))
 
*HELICITY AMPLITUDE FOR CONFIGURATION (P1,Q1,P2,Q2,PH)=(- - + + +)
      M9=-CONS1*SPLUS(Q1,Q2)**2   *(HIS(-1,1,S1)*VINI)
 
*HELICITY AMPLITUDE FOR CONFIGURATION (P1,Q1,P2,Q2,PH)=(- - + + -)
      M10=-CONS1*SMINS(P1,P2)**2  *(HIS(-1,1,S1)*DCONJG(VINI))
 
*HELICITY AMPLITUDE FOR CONFIGURATION (P1,Q1,P2,Q2,PH)=(- - - - +)
      M11=CONS1*SPLUS(Q1,P2)**2  *(HIS(-1,-1,S1)*VINI)
 
*HELICITY AMPLITUDE FOR CONFIGURATION (P1,Q1,P2,Q2,PH)=(- - - - -)
      M12=CONS1*SMINS(P1,Q2)**2  *(HIS(-1,-1,S1)*DCONJG(VINI))
 
*CROSS SECTION
      XSECTION=          (
     $        +M1*DCONJG(M1)+M2*DCONJG(M2)+M3*DCONJG(M3)
     $        +M4*DCONJG(M4)
     $        +M9*DCONJG(M9)
     $        +M10*DCONJG(M10)+M11*DCONJG(M11)+M12*DCONJG(M12)
     $                    )

*CORRECTION FOR THE FINITE MASS EFFECT AND COLLINEAR EFFECT
      X1=      P1(4)*PH(4)-P1(3)*PH(3)-P1(2)*PH(2)-P1(1)*PH(1)
      X2=      Q1(4)*PH(4)-Q1(3)*PH(3)-Q1(2)*PH(2)-Q1(1)*PH(1)
      XP=PH(4)/P1(4)
      XQ=PH(4)/Q1(4)
 
*FROM KLEISS, Z.PHYS.C33.433-443 (1987)
*CORRECTION TO THE FINITE MASS AND BREMSTR FROM INITIAL STATE
*SWITCH KEY NN=1
 
      WM=SPLUS(P1,PH)*DCONJG(SPLUS(P1,PH))/2D0/X1
     #  *SPLUS(Q1,PH)*DCONJG(SPLUS(Q1,PH))/2D0/X2
     #  *(1D0-AMEL**2/X1*XP*(1D0-XP)/(2D0-2D0*XP+XP**2))
     #  *(1D0-AMEL**2/X2*XQ*(1D0-XQ)/(2D0-2D0*XQ+XQ**2))
 
c[[[  SECTION=REAL(XSECTION)
      SECTION=REAL(XSECTION*WM)
      END

c[[      FUNCTION VINI(NUMBER)
c[[*     ***************************
c[[      IMPLICIT REAL*8(A-H,O-Z)
c[[      COMPLEX*16 SPLUS,SMINS,VINI
c[[      COMMON / MOMSET / P1(4),Q1(4),P2(4),Q2(4),PH(4)
c[[ 
c[[      VINI=SMINS(Q2,P2)/(SPLUS(PH,P1)*SPLUS(Q1,PH))
c[[      END



*CROSS SECTION FOR S CHANEL FINAL SINGL BREMSTR
*FROM KLEISS, Z.PHYS.C33.433-443 (1987)
      SUBROUTINE ERW_SklFin(SECTION)
*     ***********************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 M1,M2,M3,M4,M9,M10,M11,M12
      COMPLEX*16 VFIN,CONS1
      COMPLEX*16 XSECTION,SPLUS,SMINS,HIS
      COMMON / MOMSET / P1(4),Q1(4),P2(4),Q2(4),PH(4)
      COMMON / BHPAR2 / CMSENE,AMEL
      COMMON / BHPAR1 / CMS,AMMI
      DATA PI,ALFINV /3.1415926535897932D0, 137.03604D0/
*---------------------------------------------------------------------
      CALL ERW_INTEFA
*     ***************
 
      ALFA=1D0/ALFINV
      CONS1=DCMPLX(0D0,1D0)
 
      S0=CMSENE**2
      S1= 2D0*(P2(4)*Q2(4)-P2(3)*Q2(3)-P2(2)*Q2(2)-P2(1)*Q2(1))

*HELICITY AMPLITUDE FOR CONFIGURATION (P1,Q1,P2,Q2,PH)=(+ + + + +)
      M1=CONS1*SPLUS(P1,Q2)**2   *(HIS(1,1,S0)*VFIN(NN))
 
*HELICITY AMPLITUDE FOR CONFIGURATION (P1,Q1,P2,Q2,PH)=(+ + + + -)
      M2=CONS1*SMINS(Q1,P2)**2   *(HIS(1,1,S0)*DCONJG(VFIN(NN)))
 
*HELICITY AMPLITUDE FOR CONFIGURATION (P1,Q1,P2,Q2,PH)=(+ + - - +)
      M3=-CONS1*SPLUS(P1,P2)**2  *(HIS(1,-1,S0)*VFIN(NN))
 
*HELICITY AMPLITUDE FOR CONFIGURATION (P1,Q1,P2,Q2,PH)=(+ + - - -)
      M4=-CONS1*SMINS(Q1,Q2)**2  *(HIS(1,-1,S0)*DCONJG(VFIN(NN)))
 
*HELICITY AMPLITUDE FOR CONFIGURATION (P1,Q1,P2,Q2,PH)=(- - + + +)
      M9=-CONS1*SPLUS(Q1,Q2)**2  *(HIS(-1,1,S0)*VFIN(NN))
 
*HELICITY AMPLITUDE FOR CONFIGURATION (P1,Q1,P2,Q2,PH)=(- - + + -)
      M10=-CONS1*SMINS(P1,P2)**2  *(HIS(-1,1,S0)*DCONJG(VFIN(NN)))
 
*HELICITY AMPLITUDE FOR CONFIGURATION (P1,Q1,P2,Q2,PH)=(- - - - +)
      M11=CONS1*SPLUS(Q1,P2)**2  *(HIS(-1,-1,S0)*VFIN(NN))
 
*HELICITY AMPLITUDE FOR CONFIGURATION (P1,Q1,P2,Q2,PH)=(- - - - -)
      M12=CONS1*SMINS(P1,Q2)**2  *(HIS(-1,-1,S0)*DCONJG(VFIN(NN)))
 
*CROSS SECTION
      XSECTION=          (
     $        +M1*DCONJG(M1)+M2*DCONJG(M2)+M3*DCONJG(M3)
     $        +M4*DCONJG(M4)
     $        +M9*DCONJG(M9)
     $        +M10*DCONJG(M10)+M11*DCONJG(M11)+M12*DCONJG(M12)
     $                    )
*CORRECTION FOR THE FINITE MASS EFFECT AND COLLINEAR EFFECT
      X1=      P2(4)*PH(4)-P2(3)*PH(3)-P2(2)*PH(2)-P2(1)*PH(1)
      X2=      Q2(4)*PH(4)-Q2(3)*PH(3)-Q2(2)*PH(2)-Q2(1)*PH(1)
      XP=PH(4)/P2(4)
      XQ=PH(4)/Q2(4)
 
*FROM KLEISS, Z.PHYS.C33.433-443 (1987)
*CORRECTION TO THE FINITE MASS AND BREMSTR FROM INITIAL STATE
*SWITCH KEY NN=1
      WM=SPLUS(P2,PH)*DCONJG(SPLUS(P2,PH))/2D0/X1
     #  *SPLUS(Q2,PH)*DCONJG(SPLUS(Q2,PH))/2D0/X2
     #  *(1D0-AMMI**2/X1*XP*(1D0-XP)/(2D0-2D0*XP+XP**2))
     #  *(1D0-AMMI**2/X2*XQ*(1D0-XQ)/(2D0-2D0*XQ+XQ**2))
 
*[[   SECTION=REAL(XSECTION)
      SECTION=REAL(XSECTION*WM)
 
      END

*RADIATION FACTOR TO FINAL STATE
*FROM KLEISS, Z.PHYS.C33.433-443 (1987)
      FUNCTION VFIN(NUMBER)
*     ***************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 SPLUS,SMINS,VFIN
      COMMON / MOMSET / P1(4),Q1(4),P2(4),Q2(4),PH(4)
 
      VFIN=SMINS(Q1,P1)/(SPLUS(PH,P2)*SPLUS(Q2,PH))
 
      END


C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C>>>>
C>>>>  S P I N     A M P L I T U D E S
C>>>>
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C.. INITIAL STATE BREMSTRAHLUNG
      SUBROUTINE ERW_SingIni(SECTION)
*     ***********************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 XMSINI,XSECTION,S,SUM
      COMMON / MOMSET / P1(4),Q1(4),P2(4),Q2(4),PH(4)
      COMMON / BHPAR2 / CMSENE,AMEL
      DATA PI,ALFINV /3.1415926535897932D0, 137.03604D0/
*---------------------------------------------------------------------
      CALL ERW_INTEFA
*     *************** 
      ALFA=1D0/ALFINV
      S0=CMSENE**2
      SUM=(0D0,0D0)
*SUM OVER POSSIBLE HELICITY CONFIGURATIONS
      DO 10 I=1,3,2
      LAM1=2-I
      DO 10 K=1,3,2
      LAM2=2-K
      DO 10 N=1,3,2
      LAM3=2-N
      LAM4=2-N
      DO 10 M=1,3,2
      LEPS=2-M
      S= XMSINI( LAM1,LAM2,LAM3,LAM4,LEPS)
      SUM=SUM+S*DCONJG(S)
  10  CONTINUE
 
*CROSS SECTION
      XSECTION=SUM
      SECTION=DBLE(XSECTION)
      END
 
 
*SPIN AMPLITUDE
      FUNCTION XMSINI(L1,L2,L3,L4,LE)
*     *******************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16  HIS ,XMSINI,X1,X2
      COMPLEX *16 CR1,CR2,C,C1,CPH
      COMMON / MOM / P11(4),P12(4),P21(4),P22(4),
     $               Q11(4),Q12(4),Q21(4),Q22(4)
      COMMON /FOTON/ ARBIT(4)
      COMMON / MOMSET / P1(4),Q1(4),P2(4),Q2(4),PH(4)
      COMMON /HELP1/ C,C1,CPH,ZER(4)
      COMMON / BHPAR1 / CMS,AMMI 

      DO 5 I=1,4
  5   ZER(I)=0D0
      C=DCMPLX(1D0,0D0)
      CR1=PROPIN(P1,PH)*(1D0,0D0)
      CR2=PROPIN(Q1,PH)*(1D0,0D0)
      CPH=DSQRT(-PROPIN(PH,ARBIT)*2D0)*C
      C1=C*DELTA(L3,L4)
      S1= CMS**2 
*--------------------------------------------------------
*SPIN AMPLITUDE -BREMSTRAHLUNG FROM LINE P1
      CALL AMPLIA(L1,L2,L3,L4,LE,P1 ,PH, CR1,X1)
*--------------------------------------------------------
*SPIN AMPLITUDE -BREMSTRAHLUNG FROM LINE Q1
      CALL AMPLIB(L1,L2,L3,L4,LE,PH,Q1 ,CR2,X2)
*--------------------------------------------------------
      XMSINI=(X1+X2)*HIS(L1,L3,S1)     
      END
 
 
      SUBROUTINE AMPLIA(L1,L2,L3,L4,LE,R1A,R1B,CR,X )
*     **************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX *16   CR,X
      COMPLEX *16 C1,C,CPH
      COMPLEX *16   SEL ,SPH,SFI,SPO,S,SA,SB
      COMMON / MOMSET / P1(4),Q1(4),P2(4),Q2(4),PH(4)
      COMMON /FOTON/ ARBIT(4)
      COMMON /HELP1/ C,C1,CPH,ZER(4)
      COMMON / BHPAR2 / CMSENE,AMEL
      COMMON / BHPAR1 / CMS,AMMI
      DIMENSION SEL (2,2),SPO(2,2),SFI(2,2),SPH(2,2),S(2,2),
     $           SA (2,2),SB (2,2)
      DIMENSION R1A(4),R1B(4)

      CALL MULTI( L1,P1 ,ZER,C,-C,    LE,PH,ARBIT,CPH,CPH,   SEL, C)

      CALL MULTI( LE,ARBIT,PH ,C,C,    1,R1A,R1A,C,C,   SPH, CR)
      CALL MULTI(  1,R1A,  R1A,C,C,   L3,P2 ,Q2 ,C,C,   SFI, C1)
      CALL ILOCZ(SPH ,SFI,SA)

      CALL MULTI( LE,ARBIT,PH,C,C,     1,R1B,R1B,C,C,   SPH, -CR)
      CALL MULTI(  1,R1B, R1B,C,C,    L3,P2 ,Q2 ,C,C,   SFI, C1)
      CALL ILOCZ(SPH ,SFI,SB )

      CALL DODAJ(SA,SB,S )
 
      CALL MULTI( L4,Q2,P2,C,C,        L2,Q1 ,ZER,C,C,  SPO,  C)
 
      CALL ADD1(SEL ,S,SPO,X)
 
      END
 
      SUBROUTINE AMPLIB(L1,L2,L3,L4,LE,R1A,R1B,CR,X )
*     **************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX *16   CR,X
      COMPLEX *16 C1,C,CPH
      COMPLEX *16   SEL ,SPH,SFI,SPO,SA,SB,S
      COMMON / MOMSET / P1(4),Q1(4),P2(4),Q2(4),PH(4)
      COMMON /FOTON/ ARBIT(4)
      COMMON /HELP1/ C,C1,CPH,ZER(4)
      COMMON / BHPAR2 / CMSENE,AMEL
      COMMON / BHPAR1 / CMS,AMMI
      DIMENSION SEL (2,2),SPO(2,2),SFI(2,2),SPH(2,2),S(2,2),
     $           SA (2,2),SB (2,2)
      DIMENSION R1A(4),R1B(4)
*----------------------------------------
      S1= CMS**2

      CALL MULTI( L1,P1 ,ZER, C,C,   L3, P2 , Q2 ,   C,C,      SEL , C1)

      CALL MULTI( L4,Q2 ,P2 , C,C,   1 , R1A, R1A,   C,C,      SFI , CR)
      CALL MULTI(  1,R1A,R1A, C,C,   LE, PH , ARBIT, CPH,CPH,  SPH , C)
      CALL ILOCZ( SFI, SPH, SA )

      CALL MULTI( L4,Q2 ,P2 , C,C,   1 , R1B, R1B,   C,C,      SFI ,-CR)
      CALL MULTI(  1,R1B,R1B, C,C,   LE, PH , ARBIT, CPH,CPH,  SPH , C)
      CALL ILOCZ( SFI ,SPH, SB )

      CALL DODAJ( SA, SB, S )
 
      CALL MULTI( LE,ARBIT,PH ,C,C,  L2, Q1 , ZER,   C,C,      SPO,  C)
 
      CALL ADD1( SEL, S, SPO, X)
      END


C.. FINAL STATE BREMSTRAHLUNG
      SUBROUTINE ERW_SingFin(SECTION)
*     ***********************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 XMSFIN,XSECTION,S,SUM
      COMMON / MOMSET / P1(4),Q1(4),P2(4),Q2(4),PH(4)
      COMMON / BHPAR2 / CMSENE,AMEL
      DATA PI,ALFINV /3.1415926535897932D0, 137.03604D0/
*---------------------------------------------------------------------
      CALL ERW_INTEFA
*     ***************
      ALFA=1D0/ALFINV
      S0=CMSENE**2
      SUM=(0D0,0D0)
*SUM OVER POSSIBLE HELICITY CONFIGURATIONS
      DO 10 I=1,3,2
      LAM1=2-I
      LAM2=2-I
      DO 10 K=1,3,2
      LAM3=2-K
      DO 10 N=1,3,2
      LAM4=2-N
      DO 10 M=1,3,2
      LEPS=2-M
      S= XMSFIN( LAM1,LAM2,LAM3,LAM4,LEPS)
      SUM=SUM+S*DCONJG(S)
  10  CONTINUE
 
*CROSS SECTION
      XSECTION=SUM
      SECTION=DBLE(XSECTION)
      END
 
 
*SPIN AMPLITUDE
      FUNCTION XMSFIN(L1,L2,L3,L4,LE)
*     *******************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16  HIS ,XMSFIN,X1,X2
      COMPLEX *16 CR1,CR2,C,C1,CPH
      COMMON / BHPAR2 / CMSENE,AMEL
      COMMON / MOM / P11(4),P12(4),P21(4),P22(4),
     $               Q11(4),Q12(4),Q21(4),Q22(4)
      COMMON /FOTON/ ARBIT(4)
      COMMON / MOMSET / P1(4),Q1(4),P2(4),Q2(4),PH(4)
      COMMON /HELP1/ C,C1,CPH,ZER(4)
 
      S0=CMSENE**2
      DO 5 I=1,4
  5   ZER(I)=0D0
      C=DCMPLX(1D0,0D0)
      CR1=PROPFIN(P2,PH)*(1D0,0D0)
      CR2=PROPFIN(Q2,PH)*(1D0,0D0)
      CPH=DSQRT(-PROPIN(PH,ARBIT)*2D0)*C
      C1=C*DELTA(L1,L2)
      S0=CMSENE**2
*--------------------------------------------------------
*SPIN AMPLITUDE -BREMSTRAHLUNG FROM LINE P2
      CALL FMPLIA(L1,L2,L3,L4,LE,P2 ,PH, CR1,X1)
*--------------------------------------------------------
*SPIN AMPLITUDE -BREMSTRAHLUNG FROM LINE Q2
      CALL FMPLIB(L1,L2,L3,L4,LE,PH,Q2 ,CR2,X2)
*--------------------------------------------------------
      XMSFIN=(X1+X2)*HIS(L1,L3,S0)
      END
 
 
      SUBROUTINE FMPLIA(L1,L2,L3,L4,LE,R1A,R1B,CR,X )
*     **************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX *16   CR,X
      COMPLEX *16 C1,C,CPH
      COMPLEX *16   SEL ,SPH,SFI,SPO,S,SA,SB
      COMMON / MOMSET / P1(4),Q1(4),P2(4),Q2(4),PH(4)
      COMMON /FOTON/ ARBIT(4)
      COMMON /HELP1/ C,C1,CPH,ZER(4)
      COMMON / BHPAR2 / CMSENE,AMEL
      DIMENSION SEL (2,2),SPO(2,2),SFI(2,2),SPH(2,2),S(2,2),
     $           SA (2,2),SB (2,2)
      DIMENSION R1A(4),R1B(4)

 
      CALL MULTI( L3,P2 ,ZER,C,-C, LE,PH ,ARBIT,CPH , CPH ,SEL ,C  )

      CALL MULTI( LE,ARBIT,PH ,C   , C   ,1 ,R1A,R1A,C  , C,SPH, CR)
      CALL MULTI(  1,R1A,R1A,C, C, L1,P1 ,Q1 ,C  , C,SFI,C1)
      CALL ILOCZ(SPH ,SFI,SA)
      CALL MULTI( LE,ARBIT,PH ,C   , C   ,1 ,R1B,R1B,C  , C,SPH, CR)
      CALL MULTI(  1,R1B,R1B,C, C, L1,P1 ,Q1 ,C  , C,SFI,C1)
      CALL ILOCZ(SPH ,SFI,SB )
      CALL DODAJ(SA,SB,S )
 
      CALL MULTI( L2,Q1 ,P1 ,C, C, L4,Q2 ,ZER,C  , C,SPO ,C  )
 
      CALL ADD1(SEL ,S,SPO,X)
 
      END
 
      SUBROUTINE FMPLIB(L1,L2,L3,L4,LE,R1A,R1B,CR,X )
*     **************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX *16   CR,X
      COMPLEX *16 C1,C,CPH
      COMPLEX *16   SEL ,SPH,SFI,SPO,SA,SB,S
      COMMON / MOMSET / P1(4),Q1(4),P2(4),Q2(4),PH(4)
      COMMON /FOTON/ ARBIT(4)
      COMMON /HELP1/ C,C1,CPH,ZER(4)
      DIMENSION SEL (2,2),SPO(2,2),SFI(2,2),SPH(2,2),S(2,2),
     $           SA (2,2),SB (2,2)
      DIMENSION R1A(4),R1B(4)
 

      CALL MULTI( L3,P2 ,ZER,C,C, L1,P1 ,Q1 ,C  , C,SEL ,C1)

      CALL MULTI( L2,Q1 ,P1 ,C, C, 1 ,R1A,R1A,C  , C,SFI ,-CR )
      CALL MULTI(  1,R1A,R1A,C,C,LE,PH,ARBIT,CPH,CPH ,SPH ,C)
      CALL ILOCZ(SFI ,SPH,SA )
      CALL MULTI( L2,Q1 ,P1 ,C, C, 1 ,R1B,R1B,C  , C,SFI ,-CR )
      CALL MULTI(  1,R1B,R1B,C,C,LE,PH,ARBIT,CPH,CPH ,SPH ,C)
      CALL ILOCZ(SFI ,SPH,SB )
      CALL DODAJ(SA,SB,S )
 
      CALL MULTI( LE,ARBIT,PH ,C   , C   ,L4,Q2 ,ZER,C,C,SPO,C  )
 
      CALL ADD1(SEL ,S,SPO,X)
 
      END


*////////////////////////////////////////////////////////////////////////////////
*//                                                                            //
*//      ERW  Utilities                                                        //
*//                                                                            //
*////////////////////////////////////////////////////////////////////////////////

*SPINOR PRODUCT
*FROM KLEISS, Z.PHYS.C33.433-443 (1987)
      FUNCTION SPLUS(P1,P2)
*     ********************
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 SPLUS,X,Y
      DIMENSION P1(4),P2(4)
 
      X=DCMPLX(P1(2),P1(3))
      Y=DCMPLX(P2(2),P2(3))
      SPLUS=X*DSQRT(P2(4)-P2(1))/DSQRT(P1(4)-P1(1))
     $     -Y*DSQRT(P1(4)-P1(1))/DSQRT(P2(4)-P2(1))
      END

*SPINOR PRODUCT
*FROM KLEISS, Z.PHYS.C33.433-443 (1987)
      FUNCTION SMINS(P1,P2)
*     ********************
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 SPLUS,SMINS
      DIMENSION P1(4),P2(4)
 
      SMINS=-DCONJG(SPLUS(P1,P2))
      END

      FUNCTION HIS(NI1,NI2,X)       
*     ***************************   
*FROM BOHM,DENNIER,HOLLIK,NUCL.PHYS.B304(1988),687
      IMPLICIT REAL*8(A-H,O-Z)      
      COMPLEX *16 Z,HIS     
      COMMON / WEKINP / AMAZ,GAMMZ,SINW2,IDE    
      COMMON / COEFF  / V,A
      SAVE

      BZ0RE=  (X-AMAZ**2)   /((X-AMAZ**2)**2+GAMMZ**2*AMAZ**2)  
      BZ0IM= -GAMMZ*AMAZ    /((X-AMAZ**2)**2+GAMMZ**2*AMAZ**2)  
      Z = DCMPLX(BZ0RE,BZ0IM) 
      HIS = (1D0/X + (V+NI1*A)*(V+NI2*A)*Z )     
      END       

*PROPAGATOR-INITIAL STATE
      FUNCTION PROPIN(P1,PH)
*     **************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION P1(4),PH(4)
 
      PHPH=PH(4)**2-PH(1)**2-PH(2)**2-PH(3)**2
      P1PH=P1(4)*PH(4)-P1(1)*PH(1)-P1(2)*PH(2)-P1(3)*PH(3)
      PROPIN=1D0/(-2D0*P1PH+PHPH)
 
      END

*PROPAGATOR-FINAL STATE
      FUNCTION PROPFIN(P1,PH)
*     **************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION P1(4),PH(4)
 
      PHPH=PH(4)**2-PH(1)**2-PH(2)**2-PH(3)**2
      P1PH=P1(4)*PH(4)-P1(1)*PH(1)-P1(2)*PH(2)-P1(3)*PH(3)
      PROPFIN=1D0/( 2D0*P1PH+PHPH)
 
      END


*FUNCTION DIRAC DELTA FOR INTEGER ARGUMENTS
      FUNCTION DELTA(L1,L2)
*    ****************************************
      IMPLICIT REAL*8(A-H,O-Z)
 
      N=(1+L1*L2)
      IF (N.EQ.2) DELTA=1D0
        IF (N.EQ.0) DELTA=0D0
      END

*PRODUCT OF VB(L1,P1,Q1,A1,B1)*V(L2,P2,Q2,A2,B2)
      SUBROUTINE MULTI(L1,P1,Q1,A1,B1,L2,P2,Q2,A2,B2,SS,C)
*     *******************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 A1,A2,SS,C,B1,B2,SSPM
      DIMENSION P1(4),P2(4),Q1(4),Q2(4),SS(2,2)

      A1=DCONJG(A1)
      SS(1,1)=A1*A2*SSPM( L1, L2,P1,P2)*C
      SS(1,2)=A1*B2*SSPM( L1,-L2,P1,Q2)*C
      SS(2,1)=B1*A2*SSPM(-L1, L2,Q1,P2)*C
      SS(2,2)=B1*B2*SSPM(-L1,-L2,Q1,Q2)*C
      END

      SUBROUTINE ILOCZ(SEL,SPO,S)
*     **************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX *16   SEL,SPO,S
      DIMENSION SEL(2,2),SPO(2,2),S(2,2)
 
* MULTIPLE MATRIX S =SEL*SPO
      S(1,1)=SEL(1,1)*SPO(1,1)+SEL(1,2)*SPO(2,1)
      S(1,2)=SEL(1,1)*SPO(1,2)+SEL(1,2)*SPO(2,2)
      S(2,1)=SEL(2,1)*SPO(1,1)+SEL(2,2)*SPO(2,1)
      S(2,2)=SEL(2,1)*SPO(1,2)+SEL(2,2)*SPO(2,2)
      END

      SUBROUTINE  DODAJ(SEL,SPO,S)
*     **************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX *16   SEL,SPO,S
      DIMENSION SEL(2,2),SPO(2,2),S(2,2)
 
* ADDING   MATRIX S =SEL+SPO
      S(1,1)=SEL(1,1)+SPO(1,1)
      S(1,2)=SEL(1,2)+SPO(1,2)
      S(2,1)=SEL(2,1)+SPO(2,1)
      S(2,2)=SEL(2,2)+SPO(2,2)
      END

      SUBROUTINE ADD1(SEL,SF,SPO,X)
*     **************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX *16   SEL,SPO,S,X,SF,SS
      DIMENSION SEL(2,2),SPO(2,2),S(2,2),SF(2,2),SS(2,2)
 
* MULTIPLE MATRIX S =SEL*SF
      CALL ILOCZ(SEL,SF ,S)
* MULTIPLE MATRIX SS=S  *SPO
      CALL ILOCZ(S  ,SPO,SS)
*CONSTRACTION TO THE SCALAR OBJECT
      X=SS(1,1)+SS(1,2)+SS(2,1)+SS(2,2)
      END


*FUNCTION OF SPINOR PRODUCT OF MASSIVE PARTICLES
      FUNCTION SSPM(L1,L2,P1,P2)
*     ********************
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 SPLUS,SMINS,SSPM
      DIMENSION P1(4),P2(4)
 
      IF (P1(4).LT.1D-15) THEN
      SSPM=(0D0,0D0)
      ELSE IF (P2(4).LT.1D-15) THEN
      SSPM=(0D0,0D0)
      ELSE
      R1=P1(4)-P1(3)
      R2=P2(4)-P2(3)
      AM1=(P1(4)**2-P1(1)**2-P1(2)**2-P1(3)**2)
      AM2=(P2(4)**2-P2(1)**2-P2(2)**2-P2(3)**2)
      IF (AM1.GT.0D0) THEN
      AM1=DSQRT(AM1)
      IF (ABS(R1).LT.1D-3) AM1=-AM1
      ELSE
      AM1=0D0
      ENDIF
      IF (AM2.GT.0D0) THEN
      AM2=DSQRT(AM2)
      IF (ABS(R2).LT.1D-3) AM2=-AM2
      ELSE
      AM2=0D0
      ENDIF
      ETA2=DSQRT(P2(4)-P2(1))
      ETA1=DSQRT(P1(4)-P1(1))
      SSPM=DELTA(L1, 1)*DELTA(L2,-1)*SPLUS(P1,P2)
     $   +DELTA(L1,-1)*DELTA(L2, 1)*SMINS(P1,P2)
     $   +(DELTA(L1, 1)*DELTA(L2, 1)+DELTA(L1,-1)*DELTA(L2,-1))
     $   *(AM1*ETA2/ETA1+AM2*ETA1/ETA2)
      ENDIF
      END

*////////////////////////////////////////////////////////////////////////////////
*//                                                                            //
*//                      End of CLASS  ERW                                     //
*//                                                                            //
*////////////////////////////////////////////////////////////////////////////////
