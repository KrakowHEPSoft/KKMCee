
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C>>>>
C>>>>  I N T E R F A C E
C>>>>
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      SUBROUTINE INTEFB
C     *******************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / MOMS / PP1(4),QQ1(4),PP2(4),QQ2(4),PK1(4),PK2(4)
      COMMON / BHPAR2 / CMSENE,AMEL
      COMMON / BHPAR1 / CMS,AMMI
      DIMENSION AA(4),X(4)
      COMMON /FOTON/ ARBIT(4)


      COMMON /MOMCMS2/ P1(4),Q1(4),P2(4),Q2(4),PHOT1(4),PHOT2(4)
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
c...here reverse because PHOT2 from phase space is presampled for initial state
        PK1(K) =PHOT2(K)
        PK2(K) =PHOT1(K)
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


C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C>>>>>>
C>>>>>   S P I N   A M P L I T U D E S 
C>>>>>>
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      SUBROUTINE DUBLINI(SECTION)
*     ***********************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 XMDINI,S,SUM,C
      COMMON / MOMS / P1(4),Q1(4),P2(4),Q2(4),PK1(4),PK2(4)
      COMMON / BHPAR2 / CMSENE,AMEL
      COMMON /NIC/ XNORM
      DATA PI,ALFINV /3.1415926535897932D0, 137.03604D0/
 
      ALFA=1D0/ALFINV
      S0=CMSENE**2
      C=DCMPLX(1D0,0D0)
      SUM=0D0*C
      XNORM=ALFA**4/PI**4/S0/16D0
*SUM OVER POSSIBLE HELICITY CONFIGURATIONS
      DO 10 L=1,3,2
      LAM3=2-L
      LAM4=2-L
      DO 10 J=1,3,2
      LAM1=2-J
      DO 10 I=1,3,2
      LAM2=2-I
      DO 10 M=1,3,2
      LEPS1=2-M
      DO 10 N=1,3,2
      LEPS2=2-N
      S= XMDINI( LAM1,LAM2,LAM3,LAM4,LEPS1,LEPS2)
      S=S*DCONJG(S)
      SUM=SUM+S
  10  CONTINUE
*CROSS SECTION
      SECTION=DBLE(SUM)
      END
 
      SUBROUTINE DUBLFIN(SECTION)
*     ***********************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 XMDFIN,S,SUM,C
      COMMON / MOMS / P1(4),Q1(4),P2(4),Q2(4),PK1(4),PK2(4)
      COMMON / BHPAR2 / CMSENE,AMEL
      COMMON /NIC/ XNORM
      DATA PI,ALFINV /3.1415926535897932D0, 137.03604D0/
 
      ALFA=1D0/ALFINV
      S0=CMSENE**2
      C=DCMPLX(1D0,0D0)
      SUM=0D0*C
      XNORM=ALFA**4/PI**4/S0/16D0
 
*SUM OVER POSSIBLE HELICITY CONFIGURATIONS
      DO 10 L=1,3,2
      LAM1=2-L
      LAM2=2-L
      DO 10 J=1,3,2
      LAM3=2-J
      DO 10 I=1,3,2
      LAM4=2-I
      DO 10 M=1,3,2
      LEPS1=2-M
      DO 10 N=1,3,2
      LEPS2=2-N
      S= XMDFIN( LAM1,LAM2,LAM3,LAM4,LEPS1,LEPS2)
      S=S*DCONJG(S)
      SUM=SUM+S
 
  10  CONTINUE
*CROSS SECTION
      SECTION=DBLE(SUM)
      END
 
 
 
*SPIN AMPLITUDE  INITIAL DOUBLE
      FUNCTION XMDINI(L1,L2,L3,L4,LE1,LE2)
*     *******************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16  HIS ,XMDINI,XP1P1,XP1Q1,XQ1Q1
      COMPLEX *16 CPH1,CPH2,C,C1,CR1,CR2,S
      COMMON /HELP2/ C,C1,CPH1,CPH2,ZER(4)
      COMMON / MOMS / P1(4),Q1(4),P2(4),Q2(4),PK1(4),PK2(4)
      COMMON / CURR2 / PH1(4),PH2(4),PH(4)
      COMMON /FOTON/ ARBIT(4)
 
      S1= 2D0*(P2(4)*Q2(4)-P2(3)*Q2(3)-P2(2)*Q2(2)-P2(1)*Q2(1))
 
      DO 15 I=1,4
      PH(I)=PK1(I)+PK2(I)
 15   ZER(I)=0D0
      PHPH=PH(4)**2-PH(1)**2-PH(2)**2-PH(3)**2
      C=DCMPLX(1D0,0D0)
      C1=C
      S=(0D0,0D0)
 
      DO 10 KK=1,2
      IF (KK.EQ.1) THEN
      DO 5 II=1,4
      PH1(II)=PK1(II)
  5   PH2(II)=PK2(II)
      LL1=-LE1
      LL2=-LE2
      ELSE
      DO 6 II=1,4
      PH1(II)=PK2(II)
  6   PH2(II)=PK1(II)
      LL1=-LE2
      LL2=-LE1
      ENDIF
 
      CPH1=DSQRT(-PROPIN(PH1,ARBIT)*2D0)*C
      CPH2=DSQRT(-PROPIN(PH2,ARBIT)*2D0)*C
*--------------------------------------------------------
*SPIN AMPLITUDE - TWO BREMSTRAHLUNGS FROM LINE P1
      CR1=PROPIN(P1,PH1)*C
      CR2=PROPIN(P1,PH)*C
 
      CALL DUBLEA(L1,L2,L3,L4,LL1,LL2,P1,PH1,CR1,P1,PH1,PH2,
     $           CR2,XP1P1)
 
*--------------------------------------------------------
*SPIN AMPLITUDE - ONE BREMS FROM LINE P1,ONE FROM LINE Q1
      CR1=PROPIN(P1,PH1)*C
      CR2=PROPIN(Q1,PH2)*C
 
      CALL DUBLEC(L1,L2,L3,L4,LL1,LL2,P1,PH1,CR1,Q1,PH2,CR2,XP1Q1)
*--------------------------------------------------------
*SPIN AMPLITUDE - TWO BREMSTRAHLUNGS FROM LINE Q1
      CR1=PROPIN(Q1,PH)*C
      CR2=PROPIN(Q1,PH2)*C
 
      CALL DUBLEB(L1,L2,L3,L4,LL1,LL2,PH1,PH2,Q1,CR1,PH2,Q1,
     %            CR2,XQ1Q1)
 
*----------------------==--------------------------------
*TOTAL SPIN AMPLITUDE
      S=S+(XP1P1-XP1Q1+XQ1Q1)
 
 10   CONTINUE
      XMDINI=S*HIS(L1,L3,S1)
      END
 
*SPIN AMPLITUDE  FINAL  DOUBLE
      FUNCTION XMDFIN(L1,L2,L3,L4,LE1,LE2)
*     *******************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16  HIS ,XMDFIN,XP2P2,XP2Q2,XQ2Q2
      COMPLEX *16 CPH1,CPH2,C,C1,CR1,CR2,S
      COMMON /HELP2/ C,C1,CPH1,CPH2,ZER(4)
      COMMON / MOMS / P1(4),Q1(4),P2(4),Q2(4),PK1(4),PK2(4)
      COMMON / CURR2 / PH1(4),PH2(4),PH(4)
      COMMON /BHPAR2/ CMSENE, AMEL
      COMMON /FOTON/ ARBIT(4)
 
      S0= CMSENE**2
 
      DO 15 I=1,4
      PH(I)=PK1(I)+PK2(I)
 15   ZER(I)=0D0
      PHPH=PH(4)**2-PH(1)**2-PH(2)**2-PH(3)**2
      C=DCMPLX(1D0,0D0)
c     C1=C*DELTA(L3,L4)
      C1=C
      S=(0D0,0D0)
 
      DO 10 KK=1,2
      IF (KK.EQ.1) THEN
      DO 5 II=1,4
      PH1(II)=PK1(II)
  5   PH2(II)=PK2(II)
      LL1=-LE1
      LL2=-LE2
      ELSE
      DO 6 II=1,4
      PH1(II)=PK2(II)
  6   PH2(II)=PK1(II)
      LL1=-LE2
      LL2=-LE1
      ENDIF
 
      CPH1=DSQRT( PROPFIN(PH1,ARBIT)*2D0)*C
      CPH2=DSQRT( PROPFIN(PH2,ARBIT)*2D0)*C
*--------------------------------------------------------
*SPIN AMPLITUDE - TWO BREMSTRAHLUNGS FROM LINE P1
      CR1=PROPFIN(P2,PH1)*C
      CR2=PROPFIN(P2,PH)*C
 
      CALL FUBLEA(L1,L2,L3,L4,LL1,LL2,P2,PH1,CR1,P2,PH1,PH2,
     $           CR2,XP2P2)
 
*--------------------------------------------------------
*SPIN AMPLITUDE - ONE BREMS FROM LINE P1,ONE FROM LINE Q1
      CR1=PROPFIN(P2,PH1)*C
      CR2=PROPFIN(Q2,PH2)*C
 
      CALL FUBLEC(L1,L2,L3,L4,LL1,LL2,P2,PH1,CR1,Q2,PH2,CR2,XP2Q2)
*--------------------------------------------------------
*SPIN AMPLITUDE - TWO BREMSTRAHLUNGS FROM LINE Q1
      CR1=PROPFIN(Q2,PH)*C
      CR2=PROPFIN(Q2,PH2)*C
 
      CALL FUBLEB(L1,L2,L3,L4,LL1,LL2,PH1,PH2,Q2,CR1,PH2,Q2,
     %            CR2,XQ2Q2)
 
*----------------------==--------------------------------
*TOTAL SPIN AMPLITUDE
      S=S+(XP2P2-XP2Q2+XQ2Q2)
 
 10   CONTINUE
      XMDFIN=S*HIS(L1,L3,S0)
      END
 
      SUBROUTINE DUBLEA(L1,L2,L3,L4,LE1,LE2,R1A,R1B,CR1,
     %                  R2A,R2B,R2C,CR2,X)
*     **************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX *16   CR1,CR2,X
CZBW
      COMPLEX *16   Z
      COMPLEX *16 C1,C,CPH1,CPH2
      COMPLEX *16   SEL ,SPH1,SP,SPH2,SFI,SPO,SA,SB,S,SS,SSA,SSB,SC
      COMMON / MOMS / P1(4),Q1(4),P2(4),Q2(4),PK1(4),PK2(4)
      COMMON / CURR2 / PH1(4),PH2(4),PH(4)
      COMMON /FOTON/ ARBIT(4)
      COMMON / BHPAR2 / CMSENE,AMEL
      COMMON /HELP2/ C,C1,CPH1,CPH2,ZER(4)
      DIMENSION SEL (2,2),SPO(2,2),SFI(2,2),SPH1(2,2),SPH2(2,2),
     $          SP(2,2),SA(2,2),SB(2,2),S(2,2),SC(2,2),
     $          SS(2,2),SSA(2,2),SSB(2,2)
      DIMENSION R1A(4),R1B(4),R2A(4),R2B(4),R2C(4)
*P1 LINE
      CALL MULTI( L1,P1 ,ZER,C,-C,LE1,ARBIT,PH1,CPH1 , CPH1 ,SEL,C)
*FIRST PROPAGATOR
      CALL MULTI( LE1,PH1,ARBIT ,C, C,1  , R1A,R1A,C , C,SPH1,CR1)
      CALL MULTI(  1, R1A, R1A,C, C,LE2,ARBIT,PH2,CPH2 , CPH2 ,SP,C)
      CALL ILOCZ(SPH1,SP ,SSA)
      CALL MULTI( LE1,PH1,ARBIT ,C, C,1  , R1B,R1B,C , C,SPH1,-CR1)
      CALL MULTI(  1, R1B, R1B,C, C,LE2,ARBIT,PH2,CPH2 , CPH2 ,SP,C)
      CALL ILOCZ(SPH1,SP ,SSB)
      CALL DODAJ(SSA,SSB,SS)
*SECOND PROPAGATOR
      CALL MULTI( LE2,PH2,ARBIT,C, C,1  , R2A,R2A,C , C,SPH2,CR2)
      CALL MULTI(  1, R2A, R2A,C, C, L3,P2 ,Q2 ,C  , C,SFI,C1)
      CALL ILOCZ(SPH2,SFI,SA)
      CALL MULTI( LE2,PH2,ARBIT,C, C,1  , R2B,R2B,C , C,SPH2,-CR2)
      CALL MULTI(  1, R2B, R2B,C, C, L3,P2 ,Q2 ,C  , C,SFI,C1)
      CALL ILOCZ(SPH2,SFI,SB)
      CALL MULTI( LE2,PH2,ARBIT,C, C,1  , R2C,R2C,C , C,SPH2,-CR2)
      CALL MULTI(  1, R2C, R2C,C, C, L3,P2 ,Q2 ,C  , C,SFI,C1)
      CALL ILOCZ(SPH2,SFI,SC)
      CALL DODAJ(SA,SB,S)
      CALL DODAJ(S,SC,S)
*Q1 LINE
CZBW
      Z=0d0
      CALL MULTI( L4,Q2 ,P2 ,C, C, L2,Q1 ,ZER,C  , Z,SPO ,C  )
      CALL ILOCZ(S,SPO,S)
      CALL ADD3(SEL ,SS,S,X)
 
      END
 
      SUBROUTINE DUBLEC(L1,L2,L3,L4,LE1,LE2,R1A,R1B,CR1,R2A,R2B,CR2,X)
*     **************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX *16   CR1,CR2,X
      COMPLEX *16 C1,C,CPH1,CPH2
      COMPLEX *16   SEL ,SPH1,SP,SPH2,SFI,SPO,SA,SB,S,SS,SSA,SSB
      COMMON / MOMS / P1(4),Q1(4),P2(4),Q2(4),PK1(4),PK2(4)
      COMMON / CURR2 / PH1(4),PH2(4),PH(4)
      COMMON / BHPAR2 / CMSENE,AMEL
      COMMON /HELP2/ C,C1,CPH1,CPH2,ZER(4)
      COMMON /FOTON/ ARBIT(4)
      DIMENSION SEL (2,2),SPO(2,2),SFI(2,2),SPH1(2,2),SPH2(2,2),
     $          SP(2,2),SA(2,2),SB(2,2),S(2,2),
     $          SS(2,2),SSA(2,2),SSB(2,2)
      DIMENSION R1A(4),R1B(4),R2A(4),R2B(4)
*P1 LINE
      CALL MULTI( L1,P1 ,ZER,C,-C,LE1,ARBIT,PH1,CPH1 , CPH1 ,SEL,C)
*FIRST PROPAGATOR
      CALL MULTI( LE1,PH1,ARBIT ,C, C,1 , R1A,R1A,C , C,SPH1,CR1)
      CALL MULTI(  1, R1A, R1A,C, C,L3 ,P2 ,Q2,C , C,SFI,C1)
      CALL ILOCZ(SPH1,SFI,SSA)
      CALL MULTI( LE1,PH1,ARBIT ,C, C,1  , R1B,R1B,C , C,SPH1,-CR1)
      CALL MULTI(  1, R1B,R1B,C, C,L3 ,P2 ,Q2,C , C,SFI,C1)
      CALL ILOCZ(SPH1,SFI,SSB)
      CALL DODAJ(SSA,SSB,SS)
*SECOND PROPAGATOR
      CALL MULTI( L4 ,Q2 ,P2,C ,C ,1  , R2A,R2A,C , C,SP  , CR2)
      CALL MULTI(  1, R2A, R2A,C, C,LE2,ARBIT,PH2,CPH2,CPH2,SPH2,C)
      CALL ILOCZ(SP  ,SPH2,SA)
      CALL MULTI( L4 ,Q2 ,P2,C ,C, 1  , R2B,R2B,C , C,SP  ,-CR2)
      CALL MULTI(  1, R2B, R2B,C, C,LE2,ARBIT,PH2,CPH2,CPH2,SPH2,C)
      CALL ILOCZ(SP  ,SPH2,SB)
      CALL DODAJ(SA,SB,S)
*Q1 LINE
      CALL MULTI(LE2,PH2,ARBIT,C ,C,L2,Q1 ,ZER,C  , C,SPO,C )
 
      CALL ADD4(SEL ,SS,S,SPO,X)
 
      END
 
      SUBROUTINE DUBLEB(L1,L2,L3,L4,LE1,LE2,R1A,R1B,R1C,CR1,
     %                 R2A,R2B,CR2,X)
*     **************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX *16   CR1,CR2,X
      COMPLEX *16 C1,C,CPH1,CPH2
      COMPLEX *16   SEL ,SPH1,SP,SPH2,SFI,SPO,SA,SB,S,SS,SSA,SSB,SSC
      COMMON / MOMS / P1(4),Q1(4),P2(4),Q2(4),PK1(4),PK2(4)
      COMMON / CURR2 / PH1(4),PH2(4),PH(4)
      COMMON / BHPAR2 / CMSENE,AMEL
      COMMON /FOTON/ ARBIT(4)
      COMMON /HELP2/ C,C1,CPH1,CPH2,ZER(4)
      DIMENSION SEL (2,2),SPO(2,2),SFI(2,2),SPH1(2,2),SPH2(2,2),
     $          SP(2,2),SA(2,2),SB(2,2),S(2,2),
     $          SS(2,2),SSA(2,2),SSB(2,2),SSC(2,2)
      DIMENSION R1A(4),R1B(4),R2A(4),R2B(4),R1C(4)
*P1 LIN2
      CALL MULTI( L1,P1 ,ZER,C,-C,L3 ,P2 ,Q2,C , C ,SEL,C1)
*FIRST PROPAGATOR
      CALL MULTI( L4,Q2 ,P2 ,C, C,1  ,R1A,R1A ,C , C ,SFI ,CR1)
      CALL MULTI(  1, R1A, R1A,C, C,LE1,ARBIT,PH1,CPH1 , CPH1 ,SPH1,C )
      CALL ILOCZ(SFI ,SPH1,SSA)
      CALL MULTI( L4,Q2 ,P2 ,C, C,1  ,R1B ,R1B ,C , C ,SFI , CR1)
      CALL MULTI(  1, R1B, R1B,C, C,LE1,ARBIT,PH1,CPH1 , CPH1 ,SPH1,C )
      CALL ILOCZ(SFI ,SPH1,SSB)
      CALL DODAJ(SSA,SSB,SS)
      CALL MULTI( L4,Q2 ,P2 ,C, C,1  ,R1C ,R1C ,C , C ,SFI ,-CR1)
      CALL MULTI(  1, R1C, R1C,C, C,LE1,ARBIT,PH1,CPH1 , CPH1 ,SPH1,C )
      CALL ILOCZ(SFI ,SPH1,SSC)
      CALL DODAJ(SS,SSC,SS)
*SECOND PROPAGATOR
      CALL MULTI( LE1,PH1,ARBIT ,C, C,1  , R2A,R2A,C , C,SP  ,CR2)
      CALL MULTI(  1, R2A, R2A,C, C,LE2,ARBIT,PH2,CPH2 , CPH2 ,SPH2,C)
      CALL ILOCZ(SP  ,SPH2,SA)
      CALL MULTI( LE1,PH1,ARBIT ,C, C,1  , R2B,R2B,C , C,SP  ,-CR2)
      CALL MULTI(  1, R2B, R2B,C, C,LE2,ARBIT,PH2,CPH2 , CPH2 ,SPH2,C)
      CALL ILOCZ(SP  ,SPH2,SB)
      CALL DODAJ(SA,SB,S)
*Q1 LINE
      CALL MULTI(LE2, PH2,ARBIT,C, C, L2,Q1 ,ZER,C  , C,SPO,C )
 
      CALL ADD4(SEL ,SS,S,SPO,X)
!       write(*,*) 'dubleb a',R2B,LE2
!       write(*,*) 'dubleb b',C,C1,CPH1,CPH2,ZER
!       write(*,*) 'dubleb c',CPH2 
!      write(*,*) 'dubleb',sp,sph2
 !     write(*,*) 'a=',P1(4),Q1(4),P2(4),Q2(4),PK1(4),PK2(4)
 !      write(*,*) 'b=',PH1(4),PH2(4),PH(4)

      END
 
      SUBROUTINE FUBLEA(L1,L2,L3,L4,LE1,LE2,R1A,R1B,CR1,
     %                  R2A,R2B,R2C,CR2,X)
*     **************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX *16   CR1,CR2,X
      COMPLEX *16 C1,C,CPH1,CPH2
CZBW
      COMPLEX*16 Z
      COMPLEX *16   SEL ,SPH1,SP,SPH2,SFI,SPO,SA,SB,S,SS,SSA,SSB,SC
      COMMON / MOMS / P1(4),Q1(4),P2(4),Q2(4),PK1(4),PK2(4)
      COMMON / CURR2 / PH1(4),PH2(4),PH(4)
      COMMON /FOTON/ ARBIT(4)
      COMMON / BHPAR2 / CMSENE,AMEL
      COMMON /HELP2/ C,C1,CPH1,CPH2,ZER(4)
      DIMENSION SEL (2,2),SPO(2,2),SFI(2,2),SPH1(2,2),SPH2(2,2),
     $          SP(2,2),SA(2,2),SB(2,2),S(2,2),SC(2,2),
     $          SS(2,2),SSA(2,2),SSB(2,2)
      DIMENSION R1A(4),R1B(4),R2A(4),R2B(4),R2C(4)
*P1 LINE
      CALL MULTI( L3,P2 ,ZER,C,-C,LE1,ARBIT,PH1,CPH1 , CPH1 ,SEL,C)
*FIRST PROPAGATOR
      CALL MULTI( LE1,PH1,ARBIT ,C, C,1  , R1A,R1A,C , C,SPH1,CR1)
      CALL MULTI(  1, R1A, R1A,C, C,LE2,ARBIT,PH2,CPH2 , CPH2 ,SP,C)
      CALL ILOCZ(SPH1,SP ,SSA)
      CALL MULTI( LE1,PH1,ARBIT ,C, C,1  , R1B,R1B,C , C,SPH1, CR1)
      CALL MULTI(  1, R1B, R1B,C, C,LE2,ARBIT,PH2,CPH2 , CPH2 ,SP,C)
      CALL ILOCZ(SPH1,SP ,SSB)
      CALL DODAJ(SSA,SSB,SS)
*SECOND PROPAGATOR
      CALL MULTI( LE2,PH2,ARBIT,C, C,1  , R2A,R2A,C , C,SPH2, CR2)
      CALL MULTI(  1, R2A, R2A,C, C, L1,P1 ,Q1 ,C  , C,SFI,C1)
      CALL ILOCZ(SPH2,SFI,SA)
      CALL MULTI( LE2,PH2,ARBIT,C, C,1  , R2B,R2B,C , C,SPH2, CR2)
      CALL MULTI(  1, R2B, R2B,C, C, L1,P1 ,Q1 ,C  , C,SFI,C1)
      CALL ILOCZ(SPH2,SFI,SB)
      CALL MULTI( LE2,PH2,ARBIT,C, C,1  , R2C,R2C,C , C,SPH2, CR2)
      CALL MULTI(  1, R2C, R2C,C, C, L1,P1 ,Q1 ,C  , C,SFI,C1)
      CALL ILOCZ(SPH2,SFI,SC)
      CALL DODAJ(SA,SB,S)
      CALL DODAJ(S,SC,S)
*Q1 LINE
CZBW
      z=0d0
      CALL MULTI( L2,Q1 ,P1 ,C, C, L4,Q2 ,ZER,C  , Z,SPO ,C  )
      CALL ILOCZ(S,SPO,S)
      CALL ADD3(SEL ,SS,S,X)
 
      END
 
      SUBROUTINE FUBLEC(L1,L2,L3,L4,LE1,LE2,R1A,R1B,CR1,R2A,R2B,CR2,X)
*     **************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX *16   CR1,CR2,X
      COMPLEX *16 C1,C,CPH1,CPH2
      COMPLEX *16   SEL ,SPH1,SP,SPH2,SFI,SPO,SA,SB,S,SS,SSA,SSB
      COMMON / MOMS / P1(4),Q1(4),P2(4),Q2(4),PK1(4),PK2(4)
      COMMON / CURR2 / PH1(4),PH2(4),PH(4)
      COMMON / BHPAR2 / CMSENE,AMEL
      COMMON /HELP2/ C,C1,CPH1,CPH2,ZER(4)
      COMMON /FOTON/ ARBIT(4)
      DIMENSION SEL (2,2),SPO(2,2),SFI(2,2),SPH1(2,2),SPH2(2,2),
     $          SP(2,2),SA(2,2),SB(2,2),S(2,2),
     $          SS(2,2),SSA(2,2),SSB(2,2)
      DIMENSION R1A(4),R1B(4),R2A(4),R2B(4)
*P1 LINE
      CALL MULTI( L3,P2 ,ZER,C,-C,LE1,ARBIT,PH1,CPH1 , CPH1 ,SEL,C)
*FIRST PROPAGATOR
      CALL MULTI( LE1,PH1,ARBIT ,C, C,1 , R1A,R1A,C , C,SPH1,CR1)
      CALL MULTI(  1, R1A, R1A,C, C,L1 ,P1 ,Q1,C , C,SFI,C1)
      CALL ILOCZ(SPH1,SFI,SSA)
      CALL MULTI( LE1,PH1,ARBIT ,C, C,1  , R1B,R1B,C , C,SPH1, CR1)
      CALL MULTI(  1, R1B,R1B,C, C,L1 ,P1 ,Q1,C , C,SFI,C1)
      CALL ILOCZ(SPH1,SFI,SSB)
      CALL DODAJ(SSA,SSB,SS)
*SECOND PROPAGATOR
      CALL MULTI( L2 ,Q1 ,P1,C ,C ,1  , R2A,R2A,C , C,SP  , CR2)
      CALL MULTI(  1, R2A, R2A,C, C,LE2,ARBIT,PH2,CPH2,CPH2,SPH2,C)
      CALL ILOCZ(SP  ,SPH2,SA)
      CALL MULTI( L2 ,Q1 ,P1,C ,C, 1  , R2B,R2B,C , C,SP  , CR2)
      CALL MULTI(  1, R2B, R2B,C, C,LE2,ARBIT,PH2,CPH2,CPH2,SPH2,C)
      CALL ILOCZ(SP  ,SPH2,SB)
      CALL DODAJ(SA,SB,S)
*Q1 LINE
      CALL MULTI(LE2,PH2,ARBIT,C ,C,L4,Q2 ,ZER,C  , C,SPO,C )
 
      CALL ADD4(SEL ,SS,S,SPO,X)
 
      END
 
      SUBROUTINE FUBLEB(L1,L2,L3,L4,LE1,LE2,R1A,R1B,R1C,CR1,
     %                 R2A,R2B,CR2,X)
*     **************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX *16   CR1,CR2,X
      COMPLEX *16 C1,C,CPH1,CPH2
      COMPLEX *16   SEL ,SPH1,SP,SPH2,SFI,SPO,SA,SB,S,SS,SSA,SSB,SSC
      COMMON / MOMS / P1(4),Q1(4),P2(4),Q2(4),PK1(4),PK2(4)
      COMMON / CURR2 / PH1(4),PH2(4),PH(4)
      COMMON / BHPAR2 / CMSENE,AMEL
      COMMON /FOTON/ ARBIT(4)
      COMMON /HELP2/ C,C1,CPH1,CPH2,ZER(4)
      DIMENSION SEL (2,2),SPO(2,2),SFI(2,2),SPH1(2,2),SPH2(2,2),
     $          SP(2,2),SA(2,2),SB(2,2),S(2,2),
     $          SS(2,2),SSA(2,2),SSB(2,2),SSC(2,2)
      DIMENSION R1A(4),R1B(4),R2A(4),R2B(4),R1C(4)
*P1 LIN2
      CALL MULTI( L3,P2 ,ZER,C,-C,L1 ,P1 ,Q1,C , C ,SEL,C1)
*FIRST PROPAGATOR
      CALL MULTI( L2,Q1 ,P1 ,C, C,1  ,R1A,R1A ,C , C ,SFI ,-CR1)
      CALL MULTI(  1, R1A, R1A,C, C,LE1,ARBIT,PH1,CPH1 , CPH1 ,SPH1,C )
      CALL ILOCZ(SFI ,SPH1,SSA)
      CALL MULTI( L2,Q1 ,P1 ,C, C,1  ,R1B ,R1B ,C , C ,SFI ,-CR1)
      CALL MULTI(  1, R1B, R1B,C, C,LE1,ARBIT,PH1,CPH1 , CPH1 ,SPH1,C )
      CALL ILOCZ(SFI ,SPH1,SSB)
      CALL DODAJ(SSA,SSB,SS)
      CALL MULTI( L2,Q1 ,P1 ,C, C,1  ,R1C ,R1C ,C , C ,SFI ,-CR1)
      CALL MULTI(  1, R1C, R1C,C, C,LE1,ARBIT,PH1,CPH1 , CPH1 ,SPH1,C )
      CALL ILOCZ(SFI ,SPH1,SSC)
      CALL DODAJ(SS,SSC,SS)
*SECOND PROPAGATOR
      CALL MULTI( LE1,PH1,ARBIT ,C, C,1  , R2A,R2A,C , C,SP  ,-CR2)
      CALL MULTI(  1, R2A, R2A,C, C,LE2,ARBIT,PH2,CPH2 , CPH2 ,SPH2,C)
      CALL ILOCZ(SP  ,SPH2,SA)
      CALL MULTI( LE1,PH1,ARBIT ,C, C,1  , R2B,R2B,C , C,SP  ,-CR2)
      CALL MULTI(  1, R2B, R2B,C, C,LE2,ARBIT,PH2,CPH2 , CPH2 ,SPH2,C)
      CALL ILOCZ(SP  ,SPH2,SB)
      CALL DODAJ(SA,SB,S)
*Q1 LINE
      CALL MULTI(LE2, PH2,ARBIT,C, C, L4,Q2 ,ZER,C  , C,SPO,C )
 
      CALL ADD4(SEL ,SS,S,SPO,X)
 
      END
 
C.. SIMULTAMEUS INITIAL/FINAL BREMSTRAHLUNG
      SUBROUTINE DUBLMIX(SECTION)
*     ***********************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 XMSMIX,XSECTION,S,SUM
      COMMON / MOMS / P1(4),Q1(4),P2(4),Q2(4),PH1(4),PH2(4)
      COMMON / BHPAR2 / CMSENE,AMEL
      DATA PI,ALFINV /3.1415926535897932D0, 137.03604D0/
 
      ALFA=1D0/ALFINV
      S0=CMSENE**2
      SUM=(0D0,0D0)
*SUM OVER POSSIBLE HELICITY CONFIGURATIONS
      DO 10 KL=1,3,2
      LAM1=2-KL
      DO 10 L=1,3,2
      LAM2=2-L
      DO 10 K=1,3,2
      LAM3=2-K
      DO 10 N=1,3,2
      LAM4=2-N
      DO 10 M=1,3,2
      LEPS1=2-M
      DO 10 J=1,3,2
      LEPS2=2-J
      S= XMSMIX( LAM1,LAM2,LAM3,LAM4,LEPS1,LEPS2)
C     U=2D0*S*DCONJG(S)
C     PRINT *,U
      SUM=SUM+S*DCONJG(S)
  10  CONTINUE
 
      XSECTION=SUM
      SECTION=DBLE(XSECTION)
      END
 
 
*SPIN AMPLITUDE
      FUNCTION XMSMIX(L1,L2,L3,L4,LE1,LE2)
*     *******************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16  HIS ,XMSMIX,X1,X2,X3,X4
      COMPLEX *16 CR1,CR2,CR3,CR4,C,C1,CPH1,CPH2
      COMMON / MOMS / P1(4),Q1(4),P2(4),Q2(4),PH1(4),PH2(4)
      COMMON /VIRTUAL/ S00
      COMMON /FOTON/ ARBIT(4)
      COMMON /HELP/ C,C1,CPH1,CPH2,ZER(4)
      DIMENSION XX(4)
 
      DO 5 I=1,4
      XX(I)=P1(I)+Q1(I)-PH1(I)
  5   ZER(I)=0D0
      S00=XX(4)**2-XX(1)**2-XX(2)**2-XX(3)**2
      C=DCMPLX(1D0,0D0)
      CR1=PROPIN(P1,PH1)*(1D0,0D0)
      CR2=PROPIN(Q1,PH1)*(1D0,0D0)
      CR3=-PROPFIN(P2,PH2)*(1D0,0D0)
      CR4=-PROPFIN(Q2,PH2)*(1D0,0D0)
      CPH1=DSQRT(-PROPIN(PH1,ARBIT)*2D0)*C
      CPH2=DSQRT(-PROPIN(PH2,ARBIT)*2D0)*C
      C1=C
*--------------------------------------------------------
*SPIN AMPLITUDE -PH1 FORM P1 LINE  PH2 FROM P2 LINE
      CALL AMIXA(L1,L2,L3,L4,LE1,P1 ,PH1,CR1,LE2,P2,PH2,CR3,X1)
*--------------------------------------------------------
*SPIN AMPLITUDE -PH1 FORM P1 LINE  PH2 FROM Q2 LINE
      CALL AMIXB(L1,L2,L3,L4,LE1,P1 ,PH1,CR1,LE2,Q2,PH2,CR4,X2)
*--------------------------------------------------------
*SPIN AMPLITUDE -PH1 FORM Q1 LINE  PH2 FROM P2 LINE
      CALL AMIXC(L1,L2,L3,L4,LE1,Q1 ,PH1,CR2,LE2,P2,PH2,CR3,X3)
*--------------------------------------------------------
*SPIN AMPLITUDE -PH1 FORM Q1 LINE  PH2 FROM Q2 LINE
      CALL AMIXD(L1,L2,L3,L4,LE1,Q1 ,PH1,CR2,LE2,Q2,PH2,CR4,X4)
*--------------------------------------------------------
c....extra minus sign for final helicity state this is because in this
c....amplitudes is the reverse convention for the initial-final vertex
c....this you can check with AMPLI0 code for Born cross section

      XMSMIX=X1*HIS(L2,-L4,S00)-X2*HIS(L2,-L3,S00)
     #      -X3*HIS(L1,-L4,S00)+X4*HIS(L1,-L3,S00)


      END
*SPIN AMPLITUDE -PH1 FORM P1 LINE  PH2 FROM P2 LINE
      SUBROUTINE AMIXA(L1,L2,L3,L4,LE1,R1A,R1B,CR1,LE2,R2A,R2B,CR2,X )
*     **************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX *16   CR1,CR2,X1,X2,X3,X4,X,SSFT,C2
      COMPLEX *16 C1,C,CPH1,CPH2,CX
      COMPLEX *16   SEL ,SPH,SFI,SPO,S,SA,SB,SC
      COMMON / MOMS / P1(4),Q1(4),P2(4),Q2(4),PH1(4),PH2(4)
      COMMON /FOTON/ ARBIT(4)
      COMMON /HELP/ C,CX,CPH1,CPH2,ZER(4)
      COMMON / BHPAR2 / CMSENE,AMEL
      DIMENSION SEL (2,2),SPO(2,2),SFI(2,2),SPH(2,2),S(2,2),
     $           SA (2,2),SB (2,2),SC(2,2)
      DIMENSION R1A(4),R1B(4),R2A(4),R2B(4)
C.....P1-P2
      CALL MULTI( L3,P2 ,ZER,C,-C, LE2,ARBIT,PH2,CPH2,CPH2,SFI ,C  )
      CALL MULTI( LE2,PH2,ARBIT,C   , C   , L4,R2A,ZER,C  , C,SPH,CR2)
      CALL ILOCZ(SFI ,SPH,SB)
C..............
      CALL MULTI( L1,P1 ,ZER,C,-C, LE1,ARBIT,PH1,CPH1,CPH1,SEL ,C  )
      CALL MULTI( LE1,PH1,ARBIT,C   , C   ,1,R1A,R1A,C  , C,SPH,CR1)
      CALL ILOCZ(SEL ,SPH,SA)
      CALL MULTI( 1,R1A,R1A,C, C, L4,Q2 ,R2A,C  , C,SFI,C )
      CALL MULTI( L4,R2A,Q2 ,C, C, L2,Q1 ,ZER,C  , C,SPO ,C  )
      CALL ILOCZ(SFI,SPO,SC)
      CALL ADD1(SA,SC,SB,X1)

C.....PH1-P2
      CALL MULTI( L1,P1 ,ZER,C,-C, LE1,ARBIT,PH1,CPH1,CPH1,SEL ,C  )
      CALL MULTI( LE1,PH1,ARBIT,C   , C   , 1,R1B,R1B,C  , C,SPH,-CR1)
      CALL ILOCZ(SEL ,SPH,SA)
      CALL MULTI(  1,R1B,R1B,C, C, L4,Q2 ,R2A,C  , C,SFI,C )
      CALL MULTI( L4,R2A,Q2 ,C, C, L2,Q1 ,ZER,C  , C,SPO ,C  )
      CALL ILOCZ(SFI,SPO,SC)
      CALL ADD1(SA,SC,SB,X2)

C.....P1-PH2
      CALL MULTI( L3,P2 ,ZER,C,-C, LE2,ARBIT,PH2,CPH2,CPH2,SFI ,C  )
      CALL MULTI( LE2,PH2,ARBIT,C   , C   ,L4,R2B,ZER,C  , C,SPH,CR2)
      CALL ILOCZ(SFI ,SPH,SB)

C............
      CALL MULTI( L1,P1 ,ZER,C,-C, LE1,ARBIT,PH1,CPH1,CPH1,SEL ,C  )
      CALL MULTI( LE1,PH1,ARBIT,C   , C   , 1,R1A,R1A,C  , C,SPH, CR1)
      CALL ILOCZ(SEL ,SPH,SA)
      CALL MULTI(  1,R1A,R1A,C, C, L4,Q2 ,R2B,C  , C,SFI,C )
      CALL MULTI( L4,R2B,Q2 ,C, C, L2,Q1 ,ZER,C  , C,SPO ,C  )
      CALL ILOCZ(SFI,SPO,SC)
      CALL ADD1(SA,SC,SB,X3)
C.....PH1-PH2
      CALL MULTI( L1,P1 ,ZER,C,-C, LE1,ARBIT,PH1,CPH1,CPH1,SEL ,C  )
      CALL MULTI( LE1,PH1,ARBIT,C   , C   ,1,R1B,R1B,C  , C,SPH,-CR1)
      CALL ILOCZ(SEL ,SPH,SA)
CCC ZW: here was an error 8.03.1999
C      CALL MULTI( 1,R1B,R1B,C, C, L4,Q2 ,R2B,C  , C,SFI,C1)
      CALL MULTI( 1,R1B,R1B,C, C, L4,Q2 ,R2B,C  , C,SFI,C)
      CALL MULTI( L4,R2B,Q2 ,C, C, L2,Q1 ,ZER,C  , C,SPO ,C  )
      CALL ILOCZ(SFI,SPO,SC)
      CALL ADD1(SA,SC,SB,X4)
C
      X=X1+X2+X3+X4

 
      END
 
*SPIN AMPLITUDE -PH1 FORM Q1 LINE  PH2 FROM Q2 LINE
      SUBROUTINE AMIXC(L1,L2,L3,L4,LE1,R1A,R1B,CR1,LE2,R2A,R2B,CR2,X )
*     **************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX *16   CR1,CR2,X1,X2,X3,X4,X
      COMPLEX *16 C1,C,CPH1,CPH2,CX
      COMPLEX *16   SEL ,SPH,SFI,SPO,S,SA,SB,SC
      COMMON / MOMS / P1(4),Q1(4),P2(4),Q2(4),PH1(4),PH2(4)
      COMMON /FOTON/ ARBIT(4)
      COMMON /HELP/ C,CX,CPH1,CPH2,ZER(4)
      COMMON / BHPAR2 / CMSENE,AMEL
      DIMENSION SEL (2,2),SPO(2,2),SFI(2,2),SPH(2,2),S(2,2),
     $           SA (2,2),SB (2,2),SC(2,2)
      DIMENSION R1A(4),R1B(4),R2A(4),R2B(4)
C.....P1-P2
      CALL MULTI( L3,P2 ,ZER,C,-C, LE2,ARBIT,PH2,CPH2,CPH2,SFI ,C  )
      CALL MULTI( LE2,PH2,ARBIT,C   , C   ,L4 ,R2A,ZER,C  , C,SPH,CR2)
      CALL ILOCZ(SFI ,SPH,SB)

C......................
      CALL MULTI( L1,P1 ,ZER,C, C, L4,Q2 ,R2A,C  , C,SFI,C)
      CALL MULTI( L4,R2A,Q2 ,C, C, 1,R1A,R1A,C  , C,SPO ,C  )
      CALL ILOCZ(SFI,SPO,SA)
      CALL MULTI(1,R1A,R1A, C, C, LE1,ARBIT,PH1,CPH1,CPH1,SEL ,C  )
      CALL MULTI( LE1,PH1,ARBIT,C   , C   ,L2,Q1,ZER ,C  , C,SPH,CR1)
      CALL ILOCZ(SEL ,SPH,SC)
      CALL ADD1(SA,SC,SB,X1)
C.....PH1-P2
      CALL MULTI( L1,P1 ,ZER,C, C, L4,Q2 ,R2A,C  , C,SFI,C)
      CALL MULTI( L4,R2A,Q2 ,C, C, 1,R1B,R1B,C  , C,SPO ,C  )
      CALL ILOCZ(SFI,SPO,SA)
      CALL MULTI(1,R1B,R1B, C, C, LE1,ARBIT,PH1,CPH1,CPH1,SEL ,C  )
      CALL MULTI( LE1,PH1,ARBIT,C   , C   ,L2,Q1,ZER ,C  , C,SPH,-CR1)
      CALL ILOCZ(SEL ,SPH,SC)
      CALL ADD1(SA,SC,SB,X2)
C.....P1-PH2
      CALL MULTI( L3,P2 ,ZER,C,-C, LE2,ARBIT,PH2,CPH2,CPH2,SFI ,C  )
      CALL MULTI( LE2,PH2,ARBIT,C   , C   ,L4 ,R2B,ZER,C  , C,SPH,CR2)
      CALL ILOCZ(SFI ,SPH,SB)

C.....................
      CALL MULTI( L1,P1 ,ZER,C, C, L4,Q2 ,R2B,C  , C,SFI,C)
      CALL MULTI( L4,R2B,Q2 ,C, C, 1,R1A,R1A,C  , C,SPO ,C  )
      CALL ILOCZ(SFI,SPO,SA)
      CALL MULTI(1,R1A,R1A, C, C, LE1,ARBIT,PH1,CPH1,CPH1,SEL ,C  )
      CALL MULTI( LE1,PH1,ARBIT,C   , C   ,L2,Q1,ZER ,C  , C,SPH,CR1)
      CALL ILOCZ(SEL ,SPH,SC)
      CALL ADD1(SA,SC,SB,X3)
C.....PH1-PH2
      CALL MULTI( L1,P1 ,ZER,C, C, L4,Q2 ,R2B,C  , C,SFI,C)
      CALL MULTI( L4,R2B,Q2 ,C, C, 1,R1B,R1B,C  , C,SPO ,C  )
      CALL ILOCZ(SFI,SPO,SA)
      CALL MULTI(1,R1B,R1B, C, C, LE1,ARBIT,PH1,CPH1,CPH1,SEL ,C  )
      CALL MULTI( LE1,PH1,ARBIT,C   , C   ,L2,Q1,ZER ,C  , C,SPH,-CR1)
      CALL ILOCZ(SEL ,SPH,SC)
      CALL ADD1(SA,SC,SB,X4)
 
      X=X1+X2+X3+X4
 
      END
 
*SPIN AMPLITUDE -PH1 FORM Q1 LINE  PH2 FROM P2 LINE
      SUBROUTINE AMIXB(L1,L2,L3,L4,LE1,R1A,R1B,CR1,LE2,R2A,R2B,CR2,X )
*     **************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX *16   CR1,CR2,X1,X2,X3,X4,X
      COMPLEX *16 C1,C,CPH1,CPH2,CX
      COMPLEX *16   SEL ,SPH,SFI,SPO,S,SA,SB,SC
      COMMON / MOMS / P1(4),Q1(4),P2(4),Q2(4),PH1(4),PH2(4)
      COMMON /FOTON/ ARBIT(4)
      COMMON /HELP/ C,CX,CPH1,CPH2,ZER(4)
      COMMON / BHPAR2 / CMSENE,AMEL
      DIMENSION SEL (2,2),SPO(2,2),SFI(2,2),SPH(2,2),S(2,2),
     $           SA (2,2),SB (2,2),SC(2,2)
      DIMENSION R1A(4),R1B(4),R2A(4),R2B(4)
C.....P1-P2
      CALL MULTI( L3,R2A,ZER, C, C, LE2,ARBIT,PH2,CPH2,CPH2,SFI ,C  )
      CALL MULTI( LE2,PH2,ARBIT,C , C ,L4,Q2,ZER,C,-C,SPH,CR2)
      CALL ILOCZ(SFI ,SPH,SB)

C...........
      CALL MULTI( L1,P1 ,ZER,C,-C, LE1,ARBIT,PH1,CPH1,CPH1,SEL ,C  )
      CALL MULTI( LE1,PH1,ARBIT,C   , C   , 1,R1A,R1A,C  , C,SPH,CR1)
      CALL ILOCZ(SEL ,SPH,SA)
      CALL MULTI(  1,R1A,R1A,C, C, L3,R2A,P2 ,C  , C,SFI,C )
      CALL MULTI( L3,P2 ,R2A,C, C, L2,Q1 ,ZER,C  , C,SPO ,C  )
      CALL ILOCZ(SFI ,SPO,SC)
      CALL ADD1(SA,SC,SB,X1)


C.....PH1-P2
      CALL MULTI( L1,P1 ,ZER,C,-C, LE1,ARBIT,PH1,CPH1,CPH1,SEL ,C  )
      CALL MULTI( LE1,PH1,ARBIT,C   , C   ,1,R1B,R1B,C  , C,SPH,-CR1)
      CALL ILOCZ(SEL ,SPH,SA)
      CALL MULTI( 1,R1B,R1B,C, C, L3,R2A,P2 ,C  , C,SFI,C )
      CALL MULTI( L3,P2,R2A,C, C, L2,Q1 ,ZER,C  , C,SPO ,C  )
      CALL ILOCZ(SFI ,SPO,SC)
      CALL ADD1(SA,SC,SB,X2)
C.....P1-PH2
      CALL MULTI( L3,R2B,ZER, C, C, LE2,ARBIT,PH2,CPH2,CPH2,SFI ,C  )
      CALL MULTI( LE2,PH2,ARBIT,C   , C   ,L4,Q2,ZER,C,-C,SPH,CR2)
      CALL ILOCZ(SFI ,SPH,SB)

C................
      CALL MULTI( L1,P1 ,ZER,C,-C, LE1,ARBIT,PH1,CPH1,CPH1,SEL ,C  )
      CALL MULTI( LE1,PH1,ARBIT,C   , C   ,1,R1A,R1A,C  , C,SPH, CR1)
      CALL ILOCZ(SEL ,SPH,SA)
      CALL MULTI( 1,R1A,R1A,C, C, L3,R2B,P2 ,C  , C,SFI,C)
      CALL MULTI( L3,P2,R2B ,C, C, L2,Q1 ,ZER,C  , C,SPO ,C  )
      CALL ILOCZ(SFI ,SPO,SC)
      CALL ADD1(SA,SC,SB,X3)
C.....PH1-PH2
      CALL MULTI( L1,P1 ,ZER,C,-C, LE1,ARBIT,PH1,CPH1,CPH1,SEL ,C  )
      CALL MULTI( LE1,PH1,ARBIT,C   , C   ,1,R1B,R1B,C  , C,SPH,-CR1)
      CALL ILOCZ(SEL ,SPH,SA)
      CALL MULTI( 1,R1B,R1B,C, C, L3,R2B,P2 ,C  , C,SFI,C)
      CALL MULTI( L4,P2 ,R2B,C, C, L2,Q1 ,ZER,C  , C,SPO ,C  )
      CALL ILOCZ(SFI ,SPO,SC)
      CALL ADD1(SA,SC,SB,X4)
 
      X=X1+X2+X3+X4


      END
 
*SPIN AMPLITUDE -PH1 FORM Q1 LINE  PH2 FROM Q2 LINE
      SUBROUTINE AMIXD(L1,L2,L3,L4,LE1,R1A,R1B,CR1,LE2,R2A,R2B,CR2,X )
*     **************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX *16   CR1,CR2,X1,X2,X3,X4,X
      COMPLEX *16 C1,C,CPH1,CPH2,CX
      COMPLEX *16   SEL ,SPH,SFI,SPO,S,SA,SB,SC
      COMMON / MOMS / P1(4),Q1(4),P2(4),Q2(4),PH1(4),PH2(4)
      COMMON /FOTON/ ARBIT(4)
      COMMON /HELP/ C,CX,CPH1,CPH2,ZER(4)
      COMMON / BHPAR2 / CMSENE,AMEL
      DIMENSION SEL (2,2),SPO(2,2),SFI(2,2),SPH(2,2),S(2,2),
     $           SA (2,2),SB (2,2),SC(2,2)
      DIMENSION R1A(4),R1B(4),R2A(4),R2B(4)
C.....P1-P2
      CALL MULTI( L3,R2A,ZER, C, C, LE2,ARBIT,PH2,CPH2,CPH2,SFI ,C  )
      CALL MULTI( LE2,PH2,ARBIT,C   , C   ,L4,Q2,ZER,C,-C,SPH,CR2)
      CALL ILOCZ(SFI ,SPH,SB)


C..................
      CALL MULTI( L1,P1 ,ZER,C, C, L3,R2A,P2 ,C  , C,SFI,C)
      CALL MULTI( L3,P2 ,R2A,C, C, 1,R1A,R1A,C  , C,SPO ,C  )
      CALL ILOCZ(SFI,SPO,SA)
      CALL MULTI(1,R1A,R1A, C, C, LE1,ARBIT,PH1,CPH1,CPH1,SEL ,C  )
      CALL MULTI( LE1,PH1,ARBIT,C   , C   ,L2,Q1,ZER ,C  , C,SPH,CR1)
      CALL ILOCZ(SEL ,SPH,SC)
      CALL ADD1(SA,SC,SB,X1)

C.....PH1-P2
      CALL MULTI( L1,P1 ,ZER,C, C, L3,R2A,P2 ,C  , C,SFI,C)
      CALL MULTI( L3,P2 ,R2A,C, C, 1,R1B,R1B,C  , C,SPO ,C  )
      CALL ILOCZ(SFI,SPO,SA)
      CALL MULTI(1,R1B,R1B, C, C, LE1,ARBIT,PH1,CPH1,CPH1,SEL ,C  )
      CALL MULTI( LE1,PH1,ARBIT,C   , C   ,L2,Q1,ZER ,C  , C,SPH,-CR1)
      CALL ILOCZ(SEL ,SPH,SC)
      CALL ADD1(SA,SC,SB,X2)
C.....P1-PH2
      CALL MULTI( L3,R2B,ZER, C, C, LE2,ARBIT,PH2,CPH2,CPH2,SFI ,C  )
      CALL MULTI( LE2,PH2,ARBIT,C   , C   ,L4,Q2,ZER,C,-C,SPH,CR2)
      CALL ILOCZ(SFI ,SPH,SB)

C..............
      CALL MULTI( L1,P1 ,ZER,C, C, L3,R2B,P2 ,C  , C,SFI,C)
      CALL MULTI( L3,P2 ,R2B,C, C, 1,R1A,R1A,C  , C,SPO ,C  )
      CALL ILOCZ(SFI,SPO,SA)
      CALL MULTI(1,R1A,R1A, C, C, LE1,ARBIT,PH1,CPH1,CPH1,SEL ,C  )
      CALL MULTI( LE1,PH1,ARBIT,C   , C   ,L2,Q1,ZER ,C  , C,SPH,CR1)
      CALL ILOCZ(SEL ,SPH,SC)
      CALL ADD1(SA,SC,SB,X3)

C.....PH1-PH2
      CALL MULTI( L1,P1 ,ZER,C, C, L3,R2B,P2 ,C  , C,SFI,C)
      CALL MULTI( L3,P2 ,R2B,C, C, 1,R1B,R1B,C  , C,SPO ,C  )
      CALL ILOCZ(SFI,SPO,SA)
      CALL MULTI(1,R1B,R1B, C, C, LE1,ARBIT,PH1,CPH1,CPH1,SEL ,C  )
      CALL MULTI( LE1,PH1,ARBIT,C   , C   ,L2,Q1,ZER ,C  , C,SPH,-CR1)
      CALL ILOCZ(SEL ,SPH,SC)
      CALL ADD1(SA,SC,SB,X4) 

      X=X1+X2+X3+X4


      END
 
 
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C>>>>>>
C>>>>>    Y F S  3
C>>>>>>
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 
 
*THIS PROGRAM CALCULATED CROSS SECTION FOR INI DOUBLE BREMSTRAHLUNG
*USING COMPACT FORM AS IN YFS24M VERSION.MOST OF THE ROUTINES
*ARE TAKEN FROM THIS PROGRAM
 
*INFRARED LIMIT FOR CROSS SECTION- INITIAL STATE
      SUBROUTINE DINIINF(SECTION)
C     ******************************************
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / MOMS / P1(4),Q1(4),P2(4),Q2(4),PK1(4),PK2(4)
      COMMON / BHPAR2 / CMSENE,AMEL
      DIMENSION XX(4)
      DATA PI,ALFINV /3.1415926535897932D0, 137.03604D0/
 
      ALFA=1D0/ALFINV
      S0=CMSENE**2
      S1=(P2(4)+Q2(4))**2-(P2(1)+Q2(1))**2-(P2(2)+Q2(2))**2
     $   -(P2(3)+Q2(3))**2
      DO 10 I=1,4
      XX(I)=P2(I)+Q2(I)
  10  CONTINUE
      CALL NDIST0(XX,P1,Q1,P2,Q2,BETA00)
      CALL SFACH0(P1,Q1,PK1,SFACT1)
      CALL SFACH0(P1,Q1,PK2,SFACT2)
      SUM= SFACT1*SFACT2*BETA00
      SECTION=SUM*S0/S1
 
      END
 
*INFRARED LIMIT FOR CROSS SECTION  -FINAL STATE
      SUBROUTINE DFININF(SECTION)
C     ******************************************
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / MOMS / P1(4),Q1(4),P2(4),Q2(4),PK1(4),PK2(4)
      COMMON / BHPAR2 / CMSENE,AMEL
      DIMENSION XX(4),PPK1(4),PPK2(4)
      DATA PI,ALFINV /3.1415926535897932D0, 137.03604D0/
 
      ALFA=1D0/ALFINV
      S0=CMSENE**2
      S1=(P2(4)+Q2(4))**2-(P2(1)+Q2(1))**2-(P2(2)+Q2(2))**2
     $   -(P2(3)+Q2(3))**2
      DO 10 I=1,4
      XX(I)=P1(I)+Q1(I)
  10  CONTINUE
C note minus sign here!!!!
      DO K=1,4
       PPK1(K)=-PK1(K)
       PPK2(K)=-PK2(K)
      ENDDO
 
      CALL NDIST0(XX,P1,Q1,P2,Q2,BETA00)
      CALL SFACH0(P2,Q2,PPK1,SFACT1)
      CALL SFACH0(P2,Q2,PPK2,SFACT2)
      SUM= SFACT1*SFACT2*BETA00
      SECTION=SUM*S0/S1
 
      END
 
*COMPACT FORM FOR DOUBLE INI BREMSTR CROSS SECTION
      SUBROUTINE DINIAPR(SECTION)
C     ******************************************
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / MOMS / P1(4),Q1(4),P2(4),Q2(4),PK1(4),PK2(4)
      COMMON / BHPAR2 / CMSENE,AMEL
      DIMENSION XX(4)
      DATA PI,ALFINV /3.1415926535897932D0, 137.03604D0/
 
      ALFA=1D0/ALFINV
      S0=CMSENE**2
      S1=(P2(4)+Q2(4))**2-(P2(1)+Q2(1))**2-(P2(2)+Q2(2))**2
     $   -(P2(3)+Q2(3))**2
      DO 10 I=1,4
      XX(I)=P2(I)+Q2(I)
  10  CONTINUE

      CALL NDIST2(XX,P1,Q1,P2,Q2,PK1,PK2,DIST2)
      SECTION=DIST2*S0/S1
 
      END
      SUBROUTINE NDIST2(QQ,P1,P2,Q1,Q2,PH1,PH2,DIST2)        
C     ***********************************************        
C Provides double bremsstrahlung distribution - INITIAL state brem.
C INPUT:  P1,P2,Q1,Q2,PH1,PH2, four momenta                  
C OUTPUT: DIST2     double bremsstrahlung distribution       
C     *********************************** 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      DIMENSION QQ(*),P1(*),P2(*),Q1(*),Q2(*),PH1(*),PH2(*)  
      DIMENSION PR1(4),PR2(4),PH1R(4),PH2R(4),QR1(4),QR2(4)  
    
      CALL REDUZ2(QQ,P1,P2,PH1,PH2,PR1,PR2,PH1R,PH2R)        
      CALL REDUZ0(QQ,Q1,Q2,QR1,QR2)       
      SVAR1 = QQ(4)**2-QQ(3)**2-QQ(2)**2-QQ(1)**2            
C infrared factors from reduced momenta   
C double bremsstrahlung Xsect in next-to-leading log approx. 
      CALL GSOFA2(P1,P2,PH1,PH2,GF1,GF2)  
      CALL GTHET1(PR1,PR2,QR1,COSTH1,COSTH2)                 
      ANDI11= BORNV(SVAR1,COSTH1)         
      ANDI12= BORNV(SVAR1,COSTH2)         
      DIST2 =   GF1*ANDI11+   GF2*ANDI12  
      END              
 
      SUBROUTINE GSOFA2(P1,P2,PH1,PH2,F1,F2)
C     **************************************
C CALCULATES INGREDIENTS FOR REAL DOUBLE PHOTON DIFF. XSECTION
C     *****************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PH1(*),PH2(*),P1(*),P2(*)
C
      WM1(A,B)=     (1D0-A)**2
      WM2(A,B)=     (1D0-B)**2
      WMS(A,B)=     ((1D0-A)**2+(1D0-B)**2)
      WWM(A,B)=
     $   1D0-AM*2D0*(1D0-A)*(1D0-B)/((1D0-A)**2+(1D0-B)**2)*(A/B+B/A)
C
      PP = P1(4)*P2(4)-P1(1)*P2(1)-P1(2)*P2(2)-P1(3)*P2(3)
      AM2= P1(4)**2-P1(1)**2-P1(2)**2-P1(3)**2
      AM = AM2/(2D0*PP)
C     A1 = (P2(4)*PH1(4)-P2(1)*PH1(1)-P2(2)*PH1(2)-P2(3)*PH1(3))/PP
C     B1 = (P1(4)*PH1(4)-P1(1)*PH1(1)-P1(2)*PH1(2)-P1(3)*PH1(3))/PP
      B1 = (P2(4)*PH1(4)-P2(1)*PH1(1)-P2(2)*PH1(2)-P2(3)*PH1(3))/PP
      A1 = (P1(4)*PH1(4)-P1(1)*PH1(1)-P1(2)*PH1(2)-P1(3)*PH1(3))/PP
      SFAC1  =  2D0/(PP*A1*B1)*WWM(A1,B1)
C     A2 = (P2(4)*PH2(4)-P2(1)*PH2(1)-P2(2)*PH2(2)-P2(3)*PH2(3))/PP
C     B2 = (P1(4)*PH2(4)-P1(1)*PH2(1)-P1(2)*PH2(2)-P1(3)*PH2(3))/PP
      B2 = (P2(4)*PH2(4)-P2(1)*PH2(1)-P2(2)*PH2(2)-P2(3)*PH2(3))/PP
      A2 = (P1(4)*PH2(4)-P1(1)*PH2(1)-P1(2)*PH2(2)-P1(3)*PH2(3))/PP
      SFAC2  =  2D0/(PP*A2*B2)*WWM(A2,B2)
      A1P= A1/(1D0-A2)
      B1P= B1/(1D0-B2)
      A2P= A2/(1D0-A1)
      B2P= B2/(1D0-B1)
      IF((A1+B1).GT.(A2+B2)) THEN
        X1=WM1(A1,B1)*WMS(A2P,B2P) +WM1(A1P,B1P)*WMS(A2,B2)
        X2=WM2(A1,B1)*WMS(A2P,B2P) +WM2(A1P,B1P)*WMS(A2,B2)
      ELSE
        X1=WM1(A2,B2)*WMS(A1P,B1P) +WM1(A2P,B2P)*WMS(A1,B1)
        X2=WM2(A2,B2)*WMS(A1P,B1P) +WM2(A2P,B2P)*WMS(A1,B1)
      ENDIF
      F1 = X1*SFAC1*SFAC2/8D0
      F2 = X2*SFAC1*SFAC2/8D0
      END
 
 
*COMPACT FORM FOR DOUBLE FINAL BREMSTR CROSS SECTION
      SUBROUTINE DFINAPR(SECTION)
C     ******************************************
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / MOMS / P1(4),Q1(4),P2(4),Q2(4),PK1(4),PK2(4)
      COMMON / BHPAR2 / CMSENE,AMEL
      DIMENSION XX(4),PPK1(4),PPK2(4)
      DATA PI,ALFINV /3.1415926535897932D0, 137.03604D0/
 
      ALFA=1D0/ALFINV
      S0=CMSENE**2
      S1=(P2(4)+Q2(4))**2-(P2(1)+Q2(1))**2-(P2(2)+Q2(2))**2
     $   -(P2(3)+Q2(3))**2
      DO 10 I=1,4
      XX(I)=P1(I)+Q1(I)
  10  CONTINUE
C note minus sign here!!!!
      DO K=1,4
       PPK1(K)=-PK1(K)
       PPK2(K)=-PK2(K)
      ENDDO


      CALL FDIST2(XX,P1,Q1,P2,Q2,PPK1,PPK2,DIST2)
      
      SECTION=DIST2*S0/S1

 
      END
 
      SUBROUTINE FDIST2(QQ,P1,P2,Q1,Q2,PH1,PH2,DIST2)
C     ***********************************************
C Provides double bremsstrahlung distribution - FINAL state brem.
C INPUT:  P1,P2,Q1,Q2,PH1,PH2, four momenta
C OUTPUT: DIST2     double bremsstrahlung distribution
C     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION QQ(*),P1(*),P2(*),Q1(*),Q2(*),PH1(*),PH2(*)
      DIMENSION PR1(4),PR2(4),PH1R(4),PH2R(4),QR1(4),QR2(4)
 
      CALL REDUZ2(QQ,Q1,Q2,PH1,PH2,QR1,QR2,PH1R,PH2R)
      CALL REDUZ0(QQ,P1,P2,PR1,PR2)
      SVAR1 = QQ(4)**2-QQ(3)**2-QQ(2)**2-QQ(1)**2
C infrared factors from reduced momenta
C double bremsstrahlung Xsect in next-to-leading log approx.
      CALL GSFIN2(Q1,Q2,PH1,PH2,GF1,GF2)
C??????
      CALL GTHET1(QR1,QR2,PR1,COSTH1,COSTH2)
      ANDI11= BORNV(SVAR1,COSTH1)
      ANDI12= BORNV(SVAR1,COSTH2)
      DIST2 =   GF1*ANDI11+   GF2*ANDI12
      END
 
      SUBROUTINE GSFIN2(P1,P2,PH1,PH2,F1,F2)
C     **************************************
C CALCULATES INGREDIENTS FOR REAL DOUBLE PHOTON DIFF. XSECTION
C     *****************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PH1(*),PH2(*),P1(*),P2(*)
C
      WM (A  )=     (1D0-A)**2
      WMS(A,B)=     ((1D0-A)**2+(1D0-B)**2)
      WWM(A,B)=
     $   1D0-AM*2D0*(1D0-A)*(1D0-B)/((1D0-A)**2+(1D0-B)**2)*(A/B+B/A)
C
      PP = P1(4)*P2(4)-P1(1)*P2(1)-P1(2)*P2(2)-P1(3)*P2(3)
      AM2= P1(4)**2-P1(1)**2-P1(2)**2-P1(3)**2
      AM = AM2/(2D0*PP)
      BB1=ABS(P2(4)*PH1(4)-P2(1)*PH1(1)-P2(2)*PH1(2)-P2(3)*PH1(3))/PP
      AA1=ABS(P1(4)*PH1(4)-P1(1)*PH1(1)-P1(2)*PH1(2)-P1(3)*PH1(3))/PP
      BB2=ABS(P2(4)*PH2(4)-P2(1)*PH2(1)-P2(2)*PH2(2)-P2(3)*PH2(3))/PP
      AA2=ABS(P1(4)*PH2(4)-P1(1)*PH2(1)-P1(2)*PH2(2)-P1(3)*PH2(3))/PP
      AA1P= AA1/(1D0+AA2)
      BB1P= BB1/(1D0+BB2)
      AA2P= AA2/(1D0+AA1)
      BB2P= BB2/(1D0+BB1)
      A1  = AA1/(1+AA1+BB1)
      A2  = AA2/(1+AA2+BB2)
      B1  = BB1/(1+AA1+BB1)
      B2  = BB2/(1+AA2+BB2)
      A1P = AA1P/(1+AA1P+BB1P)
      A2P = AA2P/(1+AA2P+BB2P)
      B1P = BB1P/(1+AA1P+BB1P)
      B2P = BB2P/(1+AA2P+BB2P)
C[[[[[[[[
C     WRITE(2,'(A,4F20.10)') '++++',A1,A2,A1P,B1
      SFAC1  =  2D0/(PP*AA1*BB1)*WWM(A1,B1)
      SFAC2  =  2D0/(PP*AA2*BB2)*WWM(A2,B2)
      IF((A1+B1).GT.(A2+B2)) THEN
        X1=WM (A1   )*WMS(A2P,B2P) +WM (A1P    )*WMS(A2,B2)
        X2=WM (   B1)*WMS(A2P,B2P) +WM (    B1P)*WMS(A2,B2)
      ELSE
        X1=WM (A2   )*WMS(A1P,B1P) +WM (A2P    )*WMS(A1,B1)
        X2=WM (   B2)*WMS(A1P,B1P) +WM (    B2P)*WMS(A1,B1)
      ENDIF
      F1 = X1*SFAC1*SFAC2/8D0
      F2 = X2*SFAC1*SFAC2/8D0
C.. correction E.W. november 1989................................
C.. temporary solution
C.. this correction reconstructs double collinear limit up 50% error
C.. and affects below photon-fermion angle  <0.1 ammi/ene
C.. without correction error in this limit 1000%
      SFAC1  =  2D0/(PP*AA1*BB1)
      SFAC2  =  2D0/(PP*AA2*BB2)
      DELT=(AM2/(2D0*PP))**2*(B2**2*A1**2+A2**2*B1**2)*
     #  ( B1*B2/(A1*A2)/(A1+A2)**2
     #   +A1*A2/(B1*B2)/(B1+B2)**2  )
      WMINF=2D0*DELT/(X1+X2)
      WMM=WWM(A1,B1)*WWM(A2,B2)+WMINF
      F1 = X1*SFAC1*SFAC2/8D0*WMM
      F2 = X2*SFAC1*SFAC2/8D0*WMM
C...end of correction............................................
      END
 
      SUBROUTINE REDUZ2(QQ,P1,P2,PH1,PH2,PR1,PR2,PH1R,PH2R)  
C     *****************************************************  
C Reduction for beta2  
C           P1,P2,PH1,PH2 ==--> PR1,PR2,PH1R,PH2R            
C such that  PR1+PR2 = PH1R+PH2R+QQ       
C Input:  QQ,P1,P2,PH1,PH2                
C Output: PR1,PR2,PH1R,PH2R               
C     *********************************** 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      PARAMETER( EPS1 =1D-15)                   
      COMMON / INOUT  / NINP,NOUT   
      SAVE   / INOUT  /
      DIMENSION QQ(*), P1(*),  P2(*),  PH1(*),  PH2(*)       
      DIMENSION        PR1(*), PR2(*), PH1R(*), PH2R(*)      
      DIMENSION PP(4),QQK(4),PPK(4),PPX(4)       
      DIMENSION PX1(4),PX2(4),PH1X(4),PH2X(4),SPH(4)         
C   
      DO 20 K=1,4      
      PP(K)   = P1(K)+P2(K)               
      PPK(K)  = P1(K)+P2(K)-PH1(K)-PH2(K) 
 20   QQK(K)  = QQ(K)+PH1(K)+PH2(K)       
      SVAR  =  PP(4)**2 -PP(3)**2 -PP(2)**2 -PP(1)**2        
      SVAR1 =  QQ(4)**2 -QQ(3)**2 -QQ(2)**2 -QQ(1)**2        
      SS1   = PPK(4)**2-PPK(3)**2-PPK(2)**2-PPK(1)**2        
      SS2   = QQK(4)**2-QQK(3)**2-QQK(2)**2-QQK(1)**2  
      IF((PP(1)**2+PP(2)**2+PP(3)**2)/PP(4)**2 .GT. EPS1) THEN 
C transform all momenta to QQ rest-frame  
         CALL BOSTDQ( 1,QQ,P1 ,PX1)       
         CALL BOSTDQ( 1,QQ,P2 ,PX2)       
         CALL BOSTDQ( 1,QQ,PH1,PH1X)      
         CALL BOSTDQ( 1,QQ,PH2,PH2X)      
         CALL BOSTDQ( 1,QQ,PP ,PPX)       
C transform all momenta to PP rest-frame  
         CALL BOSTDQ( 1,PPX,PX1,PX1)      
         CALL BOSTDQ( 1,PPX,PX2,PX2)      
         CALL BOSTDQ( 1,PPX,PH1X,PH1X)    
         CALL BOSTDQ( 1,PPX,PH2X,PH2X)    
      ELSE             
C do nothing if we are already in PP rest-frame              
         DO 23 K=1,4   
            PH1X(K)=PH1(K)                
            PH2X(K)=PH2(K)                
            PX1(K)=P1(K)                  
   23       PX2(K)=P2(K)                  
      ENDIF 
C construct reduced beam momenta PR1,PR2  
C note: they are understood to be in QQ rest-frame           
      VV2   = 1D0 - SS2/SVAR              
      IF(ABS(VV2).GT.1D-6) THEN           
C construct reduced beam momenta PR1,PR2  
C start with dilatation of beams          
         DO 24 K=1,4   
         PP(K)  =  PX1(K)+PX2(K)          
  24     SPH(K) =  PH1(K)+PH2(K)          
         PK     =  PP(4)*SPH(4)           
         SK2    =  SPH(4)**2 -SPH(3)**2 -SPH(2)**2 -SPH(1)**2
CCCC     XLAM   =  SQRT((SVAR1-SK2)/SVAR+(PK/SVAR)**2)+PK/SVAR
         XLAM   =  SQRT(SVAR1/SS1) 
         AMEL2  =  P1(4)**2-P1(3)**2-P1(2)**2-P1(1)**2       
         PXMOD  =  SQRT(PX1(1)**2+PX1(2)**2+PX1(3)**2)       
         PX1(4) =  PX1(4)*XLAM            
         PX2(4) =  PX2(4)*XLAM            
CCCC     PRMOD  =  SQRT(PX1(4)**2-AMEL2)  
         PRMOD  =      PX1(4)**2-AMEL2 
         IF(PRMOD.LE.0D0) WRITE(NOUT,*) ' REDUZ2: PRMOD=', PRMOD
         IF(PRMOD.LE.0D0) WRITE(   6,*) ' REDUZ2: PRMOD=', PRMOD 
         PRMOD  = SQRT(ABS(PRMOD))        
         DO 30 K=1,3   
         PX1(K) = PX1(K)/PXMOD*PRMOD      
 30      PX2(K) = PX2(K)/PXMOD*PRMOD      
         DO 31 K=1,4   
         PH1X(K)= PH1X(K)*XLAM            
 31      PH2X(K)= PH2X(K)*XLAM            
      ENDIF            
C then, boost away the three-vector part of P1+P2-PH1-PH2    
C that is transform to QQ rest frame      
      DO 35 K=1,4      
 35   PP(K)= PX1(K)+PX2(K)-PH1X(K)-PH2X(K)                   
      CALL BOSTDQ( 1,PP,PX1,PR1)          
      CALL BOSTDQ( 1,PP,PX2,PR2)          
      CALL BOSTDQ( 1,PP,PH1X,PH1R)        
      CALL BOSTDQ( 1,PP,PH2X,PH2R)        
      END                        
 

 
*THIS PROGRAM CALCULATED CROSS SECTION FOR SIMULTANEUS INITIAL
*AND FINAL STATE SINGLE BREMSTRHALUNG
*USING COMPACT FORM AS IN YFS251 VERSION.MOST OF THE ROUTINES
*ARE TAKEN FROM THIS PROGRAM
*IS EXPLICITLY ASSUMED FOR THE POUURPOSE OF TESTING THIS VERSION OF
*YFS THAT PHOTON PK1 IS FROM INITIAL STATE, PHOTON PK2 IS FROM FINAL
*STATE
 
*INFRARED LIMIT FOR CROSS SECTION  INITIAL-FINAL STATE
      SUBROUTINE DMIXINF(SECTION)
C     ******************************************
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / MOMS / P1(4),Q1(4),P2(4),Q2(4),PK1(4),PK2(4)
      COMMON / BHPAR2 / CMSENE,AMEL
      DIMENSION XX(4),PPK2(4)
      DATA PI,ALFINV /3.1415926535897932D0, 137.03604D0/
 
      ALFA=1D0/ALFINV
      S0=CMSENE**2
      S1=(P2(4)+Q2(4))**2-(P2(1)+Q2(1))**2-(P2(2)+Q2(2))**2
     $   -(P2(3)+Q2(3))**2
      DO 10 I=1,4
      XX(I)=P1(I)+Q1(I)-PK1(I)
  10  CONTINUE
C note minus sign here!!!!
      DO K=1,4
       PPK2(K)=-PK2(K)
      ENDDO
 
      CALL NDIST0(XX,P1,Q1,P2,Q2,BETA00)
      CALL SFACH0(P1,Q1,PK1,SFACT1)
      CALL SFACH0(P2,Q2,PPK2,SFACT2)
      SUM= SFACT1*SFACT2*BETA00
      SECTION=SUM*S0/S1
 
      END
 
*COMPACT FORM FOR INITIAL/FINAL BREMSTR CROSS SECTION
      SUBROUTINE DMIXAPR(SECTION)
C     ******************************************
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / MOMS / P1(4),Q1(4),P2(4),Q2(4),PK1(4),PK2(4)
      COMMON / BHPAR2 / CMSENE,AMEL
      DIMENSION XX(4),PPK2(4)
      DATA PI,ALFINV /3.1415926535897932D0, 137.03604D0/
 
      ALFA=1D0/ALFINV
      S0=CMSENE**2
      S1=(P2(4)+Q2(4))**2-(P2(1)+Q2(1))**2-(P2(2)+Q2(2))**2
     $   -(P2(3)+Q2(3))**2
      DO 10 I=1,4
      XX(I)=P1(I)+Q1(I)-PK1(I)
C note minus sign here!!!!
      PPK2(I)=-PK2(I)
  10  CONTINUE
 
      CALL NFDIST(XX,P1,Q1,P2,Q2,PK1,PPK2,DIST2)
      SECTION=DIST2*S0/S1
 
      END
 
      SUBROUTINE NFDIST(QQ,P1,P2,Q1,Q2,PH1,PH2,DIST2)
C     ***********************************************
C Provides distribution for simultaneous initial and final state
C single bremsstrahlung.
C INPUT:  P1,P2,Q1,Q2,PH1,PH2 four momenta
C OUTPUT: DIST2 is second order result, leading+subleading log. appr.
C     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION QQ(*),P1(*),P2(*),Q1(*),Q2(*),PH1(*),PH2(*)
      DIMENSION PR1(4),PR2(4),QR1(4),QR2(4),PHR1(4),PHR2(4)
C
      CALL REDUZ1(QQ,P1,P2,PH1,PR1,PR2,PHR1)
      CALL REDUZ1(QQ,Q1,Q2,PH2,QR1,QR2,PHR2)
C Single bremsstrahlung Xsection
      CALL GSOFA1(P1,P2,PH1,GI1,GI2)
      CALL GSFIN1(Q1,Q2,PH2,GF1,GF2)
      CALL GTHET3(PR1,PR2,QR1,QR2,CTH11,CTH12,CTH21,CTH22)
      SVAR1 = QQ(4)**2-QQ(3)**2-QQ(2)**2-QQ(1)**2
      ANDI11= BORNV(SVAR1,CTH11)
      ANDI12= BORNV(SVAR1,CTH12)
      ANDI21= BORNV(SVAR1,CTH21)
      ANDI22= BORNV(SVAR1,CTH22)
      DIST2 =  GI1*GF1*ANDI11+ GI1*GF2*ANDI12
     &        +GI2*GF1*ANDI21+ GI2*GF2*ANDI22
      END
      SUBROUTINE GTHET3(P1,P2,Q1,Q2,CTH11,CTH12,CTH21,CTH22)
C     ***************************************************
C Calculates CosTh1 and CosTh2 between BEAM amd FINAL
C fermion momenta in Z RESONANCE rest frame Q1(4)+Q2(4)=0
C     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION P1(*),P2(*),Q1(*),Q2(*)
      Q1D=        SQRT(Q1(1)**2 +Q1(2)**2 +Q1(3)**2)
      Q2D=        SQRT(Q2(1)**2 +Q2(2)**2 +Q2(3)**2)
      P1D=        SQRT(P1(1)**2 +P1(2)**2 +P1(3)**2)
      P2D=        SQRT(P2(1)**2 +P2(2)**2 +P2(3)**2)
      CTH11 = (Q1(1)*P1(1) +Q1(2)*P1(2) +Q1(3)*P1(3))/Q1D/P1D
      CTH12 =-(Q1(1)*P2(1) +Q1(2)*P2(2) +Q1(3)*P2(3))/Q1D/P2D
      CTH21 =-(Q2(1)*P1(1) +Q2(2)*P1(2) +Q2(3)*P1(3))/Q2D/P1D
      CTH22 = (Q2(1)*P2(1) +Q2(2)*P2(2) +Q2(3)*P2(3))/Q2D/P2D
      END
 
