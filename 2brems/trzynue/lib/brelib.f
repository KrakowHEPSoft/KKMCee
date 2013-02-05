C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C>>>>
C>>>>  I N T E R F A C E
C>>>>
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      SUBROUTINE INTEFA
C     *******************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / MOMSET / PP1(4),QQ1(4),PP2(4),QQ2(4),PH(4)
      COMMON / BHPAR2 / CMSENE,AMEL
      COMMON / BHPAR1 / CMS,AMMI
      DIMENSION AA(4),X(4)
      COMMON /FOTON/ ARBIT(4)


      COMMON /MOMCMS1/ P1(4),Q1(4),P2(4),Q2(4),PHOT1(4)
      COMMON /FRAME/ XMSENE,YMSENE,AMINI,AMFIN


      AMEL  =AMINI
      CMSENE=XMSENE
      AMMI  =AMFIN
      CMS   =YMSENE
      DO K=1,4
        PP1(K)=P1(K)
        PP2(K)=P2(K)
        QQ1(K)=Q1(K)
        QQ2(K)=Q2(K)
        PH(K) =PHOT1(K)
      ENDDO
*===========================================INITIALIZATION FOTON
*ARBITRARY VECTOR FOR PHOTON POLARIZATION IS FIXED
*ARBITRARY VECTOR FOR PHOTON POLARIZATION IS FIXED
      CALL MOMENTA(P1,AA,X)
      ARBIT(4)=AA(4)
      ARBIT(1)=AA(1)
      ARBIT(3)=AA(3)*COS(25.)-AA(2)*SIN(25.)
      ARBIT(2)=AA(2)*COS(25.)+AA(3)*SIN(25.)

      END
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C>>>>
C>>>>  S P I N     A M P L I T U D E S
C>>>>
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 
C.. INITIAL STATE BREMSTRAHLUNG
      SUBROUTINE SINGINI(SECTION)
*     ***********************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 XMSINI,XSECTION,S,SUM
      COMMON / MOMSET / P1(4),Q1(4),P2(4),Q2(4),PH(4)
      COMMON / BHPAR2 / CMSENE,AMEL
      DATA PI,ALFINV /3.1415926535897932D0, 137.03604D0/
 
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



      CALL MULTI( L1,P1 ,ZER,C,-C, LE,PH ,ARBIT,CPH , CPH ,SEL ,C  )
*Z0 coefficients    
      CALL MULTI( LE,ARBIT,PH ,C   , C   ,1 ,R1A,R1A,C  , C,SPH, CR)
      CALL MULTI(  1,R1A,R1A,C, C, L3,P2 ,Q2 ,C  , C,SFI,C1)
      CALL ILOCZ(SPH ,SFI,SA)
      CALL MULTI( LE,ARBIT,PH ,C   , C   ,1 ,R1B,R1B,C  , C,SPH,-CR)
      CALL MULTI(  1,R1B,R1B,C, C, L3,P2 ,Q2 ,C  , C,SFI,C1)
      CALL ILOCZ(SPH ,SFI,SB )
      CALL DODAJ(SA,SB,S )
 
      CALL MULTI( L4,Q2 ,P2 ,C, C, L2,Q1 ,ZER,C  , C,SPO ,C  )
 
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

      S1= CMS**2
*Z0 coefficients    


      CALL MULTI( L1,P1 ,ZER, C ,C, L3,P2 ,Q2 ,C  , C,SEL ,C1)

      CALL MULTI( L4,Q2 ,P2 ,C, C, 1 ,R1A,R1A,C  , C,SFI , CR )
      CALL MULTI(  1,R1A,R1A,C, C, LE,PH ,ARBIT,CPH , CPH ,SPH ,C)
      CALL ILOCZ(SFI ,SPH,SA )
      CALL MULTI( L4,Q2 ,P2 ,C, C, 1 ,R1B,R1B,C  , C,SFI ,-CR )
      CALL MULTI(  1,R1B,R1B,C, C, LE,PH ,ARBIT,CPH , CPH ,SPH ,C)
      CALL ILOCZ(SFI ,SPH,SB )
      CALL DODAJ(SA,SB,S )
 
      CALL MULTI( LE,ARBIT,PH ,C   , C   ,L2,Q1 ,ZER,C,C,SPO,C  )
 
      CALL ADD1(SEL ,S,SPO,X)
 
      END
 
C.. FINAL STATE BREMSTRAHLUNG
      SUBROUTINE SINGFIN(SECTION)
*     ***********************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 XMSFIN,XSECTION,S,SUM
      COMMON / MOMSET / P1(4),Q1(4),P2(4),Q2(4),PH(4)
      COMMON / BHPAR2 / CMSENE,AMEL
      DATA PI,ALFINV /3.1415926535897932D0, 137.03604D0/
 
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
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C>>>>
C>>>>  K L E I S S   F O R M U L A S
C>>>>
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 
*CROSS SECTION FOR S CHANEL INITIAL SINGL BREMSTR
*FROM KLEISS, Z.PHYS.C33.433-443 (1987)
      SUBROUTINE SKLINI(SECTION)
*     ***********************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 M1,M2,M3,M4,M9,M10,M11,M12
      COMPLEX*16 VINI,CONS1
      COMPLEX*16 XSECTION,SPLUS,SMINS,HIS
      COMMON / MOMSET / P1(4),Q1(4),P2(4),Q2(4),PH(4)
      COMMON / BHPAR2 / CMSENE,AMEL
      COMMON / BHPAR1 / CMS,AMMI
      DATA PI,ALFINV /3.1415926535897932D0, 137.03604D0/
 
      ALFA=1D0/ALFINV
      CONS1=DCMPLX(0D0,1D0)
 
      S0=CMSENE**2
      S1= 2D0*(P2(4)*Q2(4)-P2(3)*Q2(3)-P2(2)*Q2(2)-P2(1)*Q2(1))
 
 
*HELICITY AMPLITUDE FOR CONFIGURATION (P1,Q1,P2,Q2,PH)=(+ + + + +)
      M1=CONS1*SPLUS(P1,Q2)**2
     $  *(HIS(1,1,S1)*VINI(NN))
 
*HELICITY AMPLITUDE FOR CONFIGURATION (P1,Q1,P2,Q2,PH)=(+ + + + -)
      M2=CONS1*SMINS(Q1,P2)**2
     $  *(HIS(1,1,S1)*DCONJG(VINI(NN)))
 
*HELICITY AMPLITUDE FOR CONFIGURATION (P1,Q1,P2,Q2,PH)=(+ + - - +)
      M3=-CONS1*SPLUS(P1,P2)**2
     $   *(HIS(1,-1,S1)*VINI(NN))
 
*HELICITY AMPLITUDE FOR CONFIGURATION (P1,Q1,P2,Q2,PH)=(+ + - - -)
      M4=-CONS1*SMINS(Q1,Q2)**2
     $   *(HIS(1,-1,S1)*DCONJG(VINI(NN)))
 
 
*HELICITY AMPLITUDE FOR CONFIGURATION (P1,Q1,P2,Q2,PH)=(- - + + +)
      M9=-CONS1*SPLUS(Q1,Q2)**2
     $   *(HIS(-1,1,S1)*VINI(NN))
 
*HELICITY AMPLITUDE FOR CONFIGURATION (P1,Q1,P2,Q2,PH)=(- - + + -)
      M10=-CONS1*SMINS(P1,P2)**2
     $   *(HIS(-1,1,S1)*DCONJG(VINI(NN)))
 
*HELICITY AMPLITUDE FOR CONFIGURATION (P1,Q1,P2,Q2,PH)=(- - - - +)
      M11=CONS1*SPLUS(Q1,P2)**2
     $  *(HIS(-1,-1,S1)*VINI(NN))
 
*HELICITY AMPLITUDE FOR CONFIGURATION (P1,Q1,P2,Q2,PH)=(- - - - -)
      M12=CONS1*SMINS(P1,Q2)**2
     $  *(HIS(-1,-1,S1)*DCONJG(VINI(NN)))
 
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
 
      SECTION=REAL(XSECTION)
 
      END
*CROSS SECTION FOR S CHANEL FINAL SINGL BREMSTR
*FROM KLEISS, Z.PHYS.C33.433-443 (1987)
      SUBROUTINE SKLFIN(SECTION)
*     ***********************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 M1,M2,M3,M4,M9,M10,M11,M12
      COMPLEX*16 VFIN,CONS1
      COMPLEX*16 XSECTION,SPLUS,SMINS,HIS
      COMMON / MOMSET / P1(4),Q1(4),P2(4),Q2(4),PH(4)
      COMMON / BHPAR2 / CMSENE,AMEL
      COMMON / BHPAR1 / CMS,AMMI
      DATA PI,ALFINV /3.1415926535897932D0, 137.03604D0/
 
      ALFA=1D0/ALFINV
      CONS1=DCMPLX(0D0,1D0)
 
      S0=CMSENE**2
      S1= 2D0*(P2(4)*Q2(4)-P2(3)*Q2(3)-P2(2)*Q2(2)-P2(1)*Q2(1))
 
 
*HELICITY AMPLITUDE FOR CONFIGURATION (P1,Q1,P2,Q2,PH)=(+ + + + +)
      M1=CONS1*SPLUS(P1,Q2)**2
     $  *(HIS(1,1,S0)*VFIN(NN))
 
*HELICITY AMPLITUDE FOR CONFIGURATION (P1,Q1,P2,Q2,PH)=(+ + + + -)
      M2=CONS1*SMINS(Q1,P2)**2
     $  *(HIS(1,1,S0)*DCONJG(VFIN(NN)))
 
*HELICITY AMPLITUDE FOR CONFIGURATION (P1,Q1,P2,Q2,PH)=(+ + - - +)
      M3=-CONS1*SPLUS(P1,P2)**2
     $   *(HIS(1,-1,S0)*VFIN(NN))
 
*HELICITY AMPLITUDE FOR CONFIGURATION (P1,Q1,P2,Q2,PH)=(+ + - - -)
      M4=-CONS1*SMINS(Q1,Q2)**2
     $   *(HIS(1,-1,S0)*DCONJG(VFIN(NN)))
 
 
*HELICITY AMPLITUDE FOR CONFIGURATION (P1,Q1,P2,Q2,PH)=(- - + + +)
      M9=-CONS1*SPLUS(Q1,Q2)**2
     $   *(HIS(-1,1,S0)*VFIN(NN))
 
*HELICITY AMPLITUDE FOR CONFIGURATION (P1,Q1,P2,Q2,PH)=(- - + + -)
      M10=-CONS1*SMINS(P1,P2)**2
     $   *(HIS(-1,1,S0)*DCONJG(VFIN(NN)))
 
*HELICITY AMPLITUDE FOR CONFIGURATION (P1,Q1,P2,Q2,PH)=(- - - - +)
      M11=CONS1*SPLUS(Q1,P2)**2
     $  *(HIS(-1,-1,S0)*VFIN(NN))
 
*HELICITY AMPLITUDE FOR CONFIGURATION (P1,Q1,P2,Q2,PH)=(- - - - -)
      M12=CONS1*SMINS(P1,Q2)**2
     $  *(HIS(-1,-1,S0)*DCONJG(VFIN(NN)))
 
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
 
      SECTION=REAL(XSECTION)
 
      END

C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C>>>>
C>>>>          Y F S  3
C>>>>
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 
*THIS PROGRAM CALCULATED CROSS SECTION FOR INI/FINAL SING. BREMSTR.
*USING COMPACT FORM AS IN YFS24M VERSION.MOST OF THE ROUTINES
*ARE TAKEN FROM THIS PROGRAM
 
*INFRARED LIMIT FOR SINGLE INITIAL BREMSTR CROSS SECTION
      SUBROUTINE SINIINF(SECTION)
C     ******************************************
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / MOMSET / P1(4),Q1(4),P2(4),Q2(4),PH(4)
      COMMON / BHPAR2 / CMSENE,AMEL
      COMMON / BHPAR1 / CMS,AMMI
      DIMENSION XX(4)
      DATA PI,ALFINV /3.1415926535897932D0, 137.03604D0/

 
      ALFA=1D0/ALFINV
      S0=CMSENE**2
      S1= 2D0*(P2(4)*Q2(4)-P2(3)*Q2(3)-P2(2)*Q2(2)-P2(1)*Q2(1))
     %    +AMMI**2
      DO 10 I=1,4
      XX(I)=P2(I)+Q2(I)
  10  CONTINUE
 
      CALL SFACH0(P1,Q1,PH,SFACT)
      CALL NDIST0(XX,P1,Q1,P2,Q2,BETA00)
      SUM= SFACT*BETA00
      SECTION=SUM*(S0/S1)
 
      END
      SUBROUTINE NDIST0(XX,P1,P2,Q1,Q2,ANDIS)
C     ***************************************************************
C Provides elements of beta0,             
C for transparency reasons the full reduction of momenta is done.
C     *********************************** 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      COMMON / KEYYFS / KEYGSW,KEYRAD,KEYFIX,KEYRED,KEYWGT   
      SAVE   / KEYYFS /
      DIMENSION XX(*),P1(*),P2(*),Q1(*),Q2(*)                
      DIMENSION PR1(4),PR2(4),QR1(4),QR2(4)                  
C   
      CALL REDUZ0(XX,P1,P2,PR1,PR2)       
      CALL REDUZ0(XX,Q1,Q2,QR1,QR2)       
      CALL GTHET0(PR1,QR1,COSTH1)         
      COSTH  = (PR1(1)*QR1(1) +PR1(2)*QR1(2) +PR1(3)*QR1(3)) 
     $            /SQRT((QR1(1)**2 +QR1(2)**2 +QR1(3)**2)    
     $                 *(PR1(1)**2 +PR1(2)**2 +PR1(3)**2))   
      SVAR1  = XX(4)**2-XX(3)**2-XX(2)-XX(1)**2              
      ANDIS  = BORNV (SVAR1,COSTH )       
cc      CALL BVIRT0(P1,P2,DELI1,DELI2)      
cc      CALL BVIRT0(Q1,Q2,DELF1,DELF2)      
C ...Initial/final state bremsstrahlung switches             
cc      KEYBIN  = MOD(KEYRAD,10)            
cc      KEYBFI  = MOD(KEYRAD,100)/10        
cc      DELI1   = DELI1*KEYBIN              
cc      DELI2   = DELI2*KEYBIN              
cc      DELF1   = DELF1*KEYBFI              
cc      DELF2   = DELF2*KEYBFI              
      END              
 
*COMPACT FORM FOR SINGLE INI BREMSTR CROSS SECTION
      SUBROUTINE SINIAPR(SECTION)
C     ******************************************
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / MOMSET / P1(4),Q1(4),P2(4),Q2(4),PH(4)
      COMMON / BHPAR2 / CMSENE,AMEL
      COMMON / BHPAR1 / CMS,AMMI
      DIMENSION XX(4)
      DATA PI,ALFINV /3.1415926535897932D0, 137.03604D0/
 
      ALFA=1D0/ALFINV
      S0=CMSENE**2
      S1= 2D0*(P2(4)*Q2(4)-P2(3)*Q2(3)-P2(2)*Q2(2)-P2(1)*Q2(1))
     %    +AMMI**2

      DO 10 I=1,4
      XX(I)=P2(I)+Q2(I)
  10  CONTINUE
      CALL NDIST1(XX,P1,Q1,P2,Q2,PH,DIST10,DIST11)
      SECTION=DIST10*(S0/S1)
 
      END
 
      SUBROUTINE NDIST1(QQ,P1,P2,Q1,Q2,PH,DIST,DELI1)
C     *********************************** ***********
C Provides single bremsstrahlung distribution, INITIAL STATE 
C INPUT:  P1,P2,Q1,Q2,PH, four momenta    
C OUTPUT: 
C         DIST           is first order result, exact.               
C         DISI*(1+DELI1) is second erder LL+NLL result      
C     *********************************** 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      COMMON / KEYYFS / KEYGSW,KEYRAD,KEYFIX,KEYRED,KEYWGT   
      SAVE   / KEYYFS /
      DIMENSION QQ(*),P1(*),P2(*),Q1(*),Q2(*),PH(*)          
      DIMENSION PR1(4),PR2(4),QR1(4),QR2(4),PHR(4)           
C   
      CALL REDUZ1(QQ,P1,P2,PH,PR1,PR2,PHR)                   
      CALL REDUZ0(QQ,Q1,Q2,QR1,QR2)       
C Single bremsstrahlung Xsection          
      CALL GSOFA1(P1,P2,PH,GF1,GF2)       
      CALL GTHET1(PR1,PR2,QR1,COSTH1,COSTH2)                 
      SVAR1 = QQ(4)**2-QQ(3)**2-QQ(2)**2-QQ(1)**2            
      ANDI11= BORNV(SVAR1,COSTH1)         
      ANDI12= BORNV(SVAR1,COSTH2)         
      DIST  =  GF1*ANDI11+ GF2*ANDI12     
C Virtual correction in the second order case (K-factor style)
cc      CALL BVIRT1(P1,P2,PH,DELI1)         
C ...Initial/final state bremsstrahlung switches             
cc      KEYBIN  = MOD(KEYRAD,10)            
cc      DELI1   = DELI1*KEYBIN 
        DELI1=0D0             
      END              
 
*INFRARED LIMIT FOR SINGLE FINAL BREMSTR CROSS SECTION
      SUBROUTINE SFININF(SECTION)
C     ******************************************
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / MOMSET / P1(4),Q1(4),P2(4),Q2(4),PH(4)
      COMMON / BHPAR2 / CMSENE,AMEL
      COMMON / BHPAR1 / CMS,AMMI
      DIMENSION XX(4),PHK(4)
      DATA PI,ALFINV /3.1415926535897932D0, 137.03604D0/
 
      ALFA=1D0/ALFINV
      S0=CMSENE**2
      S1= 2D0*(P2(4)*Q2(4)-P2(3)*Q2(3)-P2(2)*Q2(2)-P2(1)*Q2(1))
     %    +AMMI**2
      DO K=1,4
       PHK(K)=-PH(K)
      ENDDO
      DO 10 I=1,4
      XX(I)=P1(I)+Q1(I)
  10  CONTINUE
 
      CALL SFACH0(P2,Q2,PHK,SFACT)
      CALL NDIST0(XX,P1,Q1,P2,Q2,BETA00)
      SUM= SFACT*BETA00
      SECTION=SUM*(S0/S1)
 
      END
 
*COMPACT FORM FOR SINGLE FINAL BREMSTR CROSS SECTION
      SUBROUTINE SFINAPR(SECTION)
C     ******************************************
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / MOMSET / P1(4),Q1(4),P2(4),Q2(4),PH(4)
      COMMON / BHPAR2 / CMSENE,AMEL
      COMMON / BHPAR1 / CMS,AMMI
      DIMENSION XX(4),PPK(4)
      DATA PI,ALFINV /3.1415926535897932D0, 137.03604D0/
 
      ALFA=1D0/ALFINV
      S0=CMSENE**2
      S1= 2D0*(P2(4)*Q2(4)-P2(3)*Q2(3)-P2(2)*Q2(2)-P2(1)*Q2(1))
     %    +AMMI**2
      DO 10 I=1,4
  10  XX(I)=P1(I)+Q1(I)
C note minus sign here!!!!
      DO K=1,4
       PPK(K)=-PH(K)
      ENDDO
 
      CALL FDIST1(XX,P1,Q1,P2,Q2,PPK,DIST10,DIST11)
      SECTION=DIST10*(S0/S1)

 
      END
 
      SUBROUTINE FDIST1(QQ,P1,P2,Q1,Q2,PH,DIST10,DIST11)
C     **************************************************
C Provides FIRST ORDER single FINAL state bremsstrahlung distribution
C INPUT:  QQ,P1,P2,Q1,Q2,PH, four momenta
C OUTPUT:
C         DIST10 is first order result, exact.
C     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION QQ(*),P1(*),P2(*),Q1(*),Q2(*),PH(*)
      COMMON / KEYYFS / KEYGSW,KEYRAD,KEYFIX,KEYRED,KEYWGT
      DIMENSION PR1(4),PR2(4),QR1(4),QR2(4),PHR(4)
C
      CALL REDUZ1(QQ,Q1,Q2,PH,QR1,QR2,PHR)
      CALL REDUZ0(QQ,P1,P2,PR1,PR2)
      SVAR1 = QQ(4)**2-QQ(3)**2-QQ(2)**2-QQ(1)**2
C Infrared factor from reduced momenta
C Single bremsstrahlung Xsection
      CALL GSFIN1(Q1,Q2,PH,GF1,GF2)
C????
      CALL GTHET1(QR1,QR2,PR1,COSTH1,COSTH2)
      ANDI11= BORNV(SVAR1,COSTH1)
      ANDI12= BORNV(SVAR1,COSTH2)
      DIST  = GF1*ANDI11+ GF2*ANDI12
C Virtual correction in the second order case (K-factor style)
C     CALL BVIRT1(Q1,Q2,PH,DELF1)
C     CALL BVIRT0(P1,P2,DELI1,DUMMY)
C ...Initial/final state bremsstrahlung switches
C     KEYBIN  = MOD(KEYRAD,10)
C     KEYBFI  = MOD(KEYRAD,100)/10
C     DELI1   = DELI1*KEYBIN
C     DELF1   = DELF1*KEYBFI
C     DIST11=  DIST*(1D0+DELF1+DELI1)
      DIST11=  DIST
      DIST10=  DIST
      END
 
