
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C>>>>
C>>>>  I N T E R F A C E
C>>>>
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      SUBROUTINE INTEFC
C     *******************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / MOMS3 / 
     #       PP1(4),QQ1(4),PP2(4),QQ2(4),PK1(4),PK2(4),PK3(4)
      COMMON / BHPAR2 / CMSENE,AMEL
      COMMON / BHPAR1 / CMS,AMMI
      DIMENSION AA(4),X(4)
      COMMON /FOTON/ ARBIT(4)


      COMMON /MOMCMS3/ 
     #       P1(4),Q1(4),P2(4),Q2(4),PHOT1(4),PHOT2(4),PHOT3(4)
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
        PK3(K) =PHOT3(K)
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
      SUBROUTINE TRZYINI(SECTION)
*     ***********************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 XMTINI,S,SUM,C
      COMMON / MOMS3 / P1(4),Q1(4),P2(4),Q2(4),PK1(4),PK2(4),PK3(4)
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
      DO 10 NN=1,3,2
      LEPS3=2-NN
      S= XMTINI( LAM1,LAM2,LAM3,LAM4,LEPS1,LEPS2,LEPS3)
      S=S*DCONJG(S)
      SUM=SUM+S 
 
  10  CONTINUE
*CROSS SECTION  

      SECTION=DBLE(SUM)
      END
 
 
 
*SPIN AMPLITUDE  INITIAL TRIPLE
      FUNCTION XMTINI(L1,L2,L3,L4,LE1,LE2,LE3)
*     ******************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16  HIS ,XMTINI,XP1P1P1,XP1P1Q1,XP1Q1Q1,XQ1Q1Q1
      COMPLEX *16 CPH1,CPH2,CPH3,C,C1,CR1,CR2,CR3,S
      COMMON /HELP3/ C,C1,CPH1,CPH2,CPH3,ZER(4)
      COMMON / MOMS3 / P1(4),Q1(4),P2(4),Q2(4),PK1(4),PK2(4),PK3(4)
      COMMON / CURR3 / PH1(4),PH2(4),PH(4),PH3(4),PH123(4),PH23(4)
      COMMON /FOTON/ ARBIT(4)
      COMPLEX X1PK1,X2PK1,X1PK2,X2PK2,X1PK3,X2PK3
      COMPLEX X1, X2
 
      S1= 2D0*(P2(4)*Q2(4)-P2(3)*Q2(3)-P2(2)*Q2(2)-P2(1)*Q2(1))
 
      C=DCMPLX(1D0,0D0)
      C1=C
      S=(0D0,0D0)
      X1=(0D0,0D0)
      X2=(0D0,0D0)
      X1PK1=(0D0,0D0)
      X2PK1=(0D0,0D0)
      X1PK2=(0D0,0D0)
      X2PK2=(0D0,0D0)
      X1PK3=(0D0,0D0)
      X2PK3=(0D0,0D0)

      DO 10  KK=1,6
      IF (KK.EQ.1) THEN
      npk1=1.
      npk2=2.
      npk3=3.
      DO II=1,4
       PH1(II)=PK1(II)
       PH2(II)=PK2(II)
       PH3(II)=PK3(II)
      ENDDO
      LL1=-LE1
      LL2=-LE2     
      LL3=-LE3
      ELSEIF(KK.EQ.2) THEN
      npk1=1.
      npk2=3.
      npk3=2.
      DO II=1,4
       PH1(II)=PK1(II)
       PH2(II)=PK3(II)
       PH3(II)=PK2(II)
      ENDDO
      LL1=-LE1
      LL2=-LE3
      LL3=-LE2
      ELSEIF(KK.EQ.3) THEN
      npk1=3.
      npk2=2.
      npk3=1.
      DO II=1,4
       PH1(II)=PK3(II)
       PH2(II)=PK2(II)
       PH3(II)=PK1(II)
      ENDDO
      LL1=-LE3
      LL2=-LE2
      LL3=-LE1
      ELSEIF(KK.EQ.4) THEN
      npk1=2.
      npk2=3.
      npk3=1.
      DO II=1,4
       PH1(II)=PK3(II)
       PH2(II)=PK1(II)
       PH3(II)=PK2(II)
      ENDDO
      LL1=-LE3
      LL2=-LE1
      LL3=-LE2
      ELSEIF(KK.EQ.5) THEN
      npk1=3.
      npk2=1.
      npk3=2.
      DO II=1,4
       PH1(II)=PK2(II)
       PH2(II)=PK3(II)
       PH3(II)=PK1(II)
      ENDDO
      LL1=-LE2
      LL2=-LE3
      LL3=-LE1
      ELSEIF(KK.EQ.6) THEN
      npk1=2.
      npk2=1.
      npk3=3.
      DO II=1,4
       PH1(II)=PK2(II)
       PH2(II)=PK1(II)
       PH3(II)=PK3(II)
      ENDDO
      LL1=-LE2
      LL2=-LE1
      LL3=-LE3
      ENDIF



      DO 15 I=1,4
      PH(I)   =PH1(I)+PH2(I)
      PH123(I)=PH1(I)+PH2(I) +PH3(I)
      PH23(I) =PH2(I)+PH3(I)
 15   ZER(I)  =0D0
      PHPH=PH(4)**2-PH(1)**2-PH(2)**2-PH(3)**2
      PHPH123=PH123(4)**2-PH123(1)**2-PH123(2)**2-PH123(3)**2
 
      CPH1=DSQRT(-PROPIN(PH1,ARBIT)*2D0)*C
      CPH2=DSQRT(-PROPIN(PH2,ARBIT)*2D0)*C
      CPH3=DSQRT(-PROPIN(PH3,ARBIT)*2D0)*C


*--------------------------------------------------------
*SPIN AMPLITUDE - THREE  BREMSTRAHLUNGS FROM LINE P1
      CR1=PROPIN(P1,PH1)*C
      CR2=PROPIN(P1,PH)*C
      CR3=PROPIN(P1,PH123)*C
 
      CALL TRZYLEA(L1,L2,L3,L4,LL1,LL2,LL3,P1,PH1,CR1,P1,PH1,PH2,
     $           CR2,P1,PH1,PH2,PH3,CR3,XP1P1P1)

      if(npk1.eq.1)X1PK1=X1PK1+XP1P1P1
      if(npk1.eq.2)X1PK1=X1PK1+XP1P1P1
      if(npk1.eq.3)X1PK1=X1PK1+XP1P1P1
      if(npk2.eq.1)X1PK2=X1PK2+XP1P1P1
      if(npk2.eq.2)X1PK2=X1PK2+XP1P1P1
      if(npk2.eq.3)X1PK2=X1PK2+XP1P1P1
      if(npk3.eq.1)X1PK3=X1PK3+XP1P1P1
      if(npk3.eq.2)X1PK3=X1PK3+XP1P1P1
      if(npk3.eq.3)X1PK3=X1PK3+XP1P1P1

*--------------------------------------------------------
*SPIN AMPLITUDE - TWO BREMSTRAHLUNGS FROM LINE P1, ONE FROM LINE Q1
      CR1=PROPIN(P1,PH1)*C
      CR2=PROPIN(P1,PH) *C
      CR3=PROPIN(Q1,PH3)*C
 
      CALL TRZYLEC(L1,L2,L3,L4,LL1,LL2,LL3,P1,PH1,CR1,P1,PH1,PH2,
     $           CR2,PH3,Q1,CR3,XP1P1Q1)

      if(npk1.eq.1)X1PK1=X1PK1-XP1P1Q1
      if(npk1.eq.2)X1PK1=X1PK1-XP1P1Q1
      if(npk1.eq.3)X2PK1=X2PK1-XP1P1Q1
      if(npk2.eq.1)X1PK2=X1PK2-XP1P1Q1
      if(npk2.eq.2)X1PK2=X1PK2-XP1P1Q1
      if(npk2.eq.3)X2PK2=X2PK2-XP1P1Q1
      if(npk3.eq.1)X1PK3=X1PK3-XP1P1Q1
      if(npk3.eq.2)X1PK3=X1PK3-XP1P1Q1
      if(npk3.eq.3)X2PK3=X2PK3-XP1P1Q1

*--------------------------------------------------------
*SPIN AMPLITUDE - ONE BREMS FROM LINE P1,TWO FROM LINE Q1
      CR1=PROPIN(P1,PH1)*C
      CR2=PROPIN(Q1,PH23)*C
      CR3=PROPIN(Q1,PH3)*C
 
      CALL TRZYLED(L1,L2,L3,L4,LL1,LL2,LL3,P1,PH1,CR1,PH2,PH3,Q1,
     $             CR2,PH3,Q1,CR3,XP1Q1Q1)

      if(npk1.eq.1)X1PK1=X1PK1+XP1Q1Q1
      if(npk1.eq.2)X2PK1=X2PK1+XP1Q1Q1
      if(npk1.eq.3)X2PK1=X2PK1+XP1Q1Q1
      if(npk2.eq.1)X1PK2=X1PK2+XP1Q1Q1
      if(npk2.eq.2)X2PK2=X2PK2+XP1Q1Q1
      if(npk2.eq.3)X2PK2=X2PK2+XP1Q1Q1
      if(npk3.eq.1)X1PK3=X1PK3+XP1Q1Q1
      if(npk3.eq.2)X2PK3=X2PK3+XP1Q1Q1
      if(npk3.eq.3)X2PK3=X2PK3+XP1Q1Q1

*------------------------------------------------------    
*SPIN AMPLITUDE - THREE BREMSTRAHLUNGS FROM LINE Q1
      CR1=PROPIN(Q1,PH123)*C
      CR2=PROPIN(Q1,PH23)*C
      CR3=PROPIN(Q1,PH3)*C
 
      CALL TRZYLEB(L1,L2,L3,L4,LL1,LL2,LL3,PH1,PH2,PH3,Q1,CR1,PH2,
     $            PH3,Q1,CR2,PH3,Q1,CR3,XQ1Q1Q1)

      if(npk1.eq.1)X2PK1=X2PK1-XQ1Q1Q1
      if(npk1.eq.2)X2PK1=X2PK1-XQ1Q1Q1
      if(npk1.eq.3)X2PK1=X2PK1-XQ1Q1Q1
      if(npk2.eq.1)X2PK2=X2PK2-XQ1Q1Q1
      if(npk2.eq.2)X2PK2=X2PK2-XQ1Q1Q1
      if(npk2.eq.3)X2PK2=X2PK2-XQ1Q1Q1
      if(npk3.eq.1)X2PK3=X2PK3-XQ1Q1Q1
      if(npk3.eq.2)X2PK3=X2PK3-XQ1Q1Q1
      if(npk3.eq.3)X2PK3=X2PK3-XQ1Q1Q1

*----------------------==--------------------------------
*TOTAL SPIN AMPLITUDE
      S=S+ (XP1P1P1-XP1P1Q1+XP1Q1Q1-XQ1Q1Q1)

 
 10   CONTINUE

      XMTINI=S*HIS(L1,L3,S1)
cc        write(6,*)'>>>>>>>>>'
cc        write(6,*)'section',l1,l2,l3,l4,ll1,ll2,ll3,s*dconjg(s)
cc        write(6,*)'x1-pk1',x1pk1
cc        write(6,*)'x2-pk1',x2pk1
cc        write(6,*)'x1-pk2',x1pk2
cc        write(6,*)'x2-pk2',x2pk2
cc        write(6,*)'x1-pk3',x1pk3
cc        write(6,*)'x2-pk3',x2pk3
cc        write(6,*)'>>>>>>>>>'

 
      END
 
      SUBROUTINE TRZYLEA(L1,L2,L3,L4,LE1,LE2,LE3,R1A,R1B,CR1,
     %             R2A,R2B,R2C,CR2,R3A,R3B,R3C,R3D,CR3,X)
*     **************************************************
      IMPLICIT NONE
      INTEGER L1, L2, L3, L4, LE1, LE2, LE3
      COMPLEX *16   CR1,CR2,CR3,X,X1,X2
      COMPLEX *16 C1,C,CPH1,CPH2,CPH3
      COMPLEX *16   SEL ,SPH1,SP,SPH2,SFI,SPO,SA,SB,S,SS,
     $              SSA,SSB,SC,SSS,SPH3,SSSA,SSSB,SSSC,SSSD
      REAL *8 P1,Q1,P2,Q2,PK1,PK2,PK3
      COMMON / MOMS3 / P1(4),Q1(4),P2(4),Q2(4),PK1(4),PK2(4),PK3(4)
      REAL *8 PH1,PH2,PH,PH3,PH123,PH23
      COMMON / CURR3 / PH1(4),PH2(4),PH(4),PH3(4),PH123(4),PH23(4)
      REAl *8 ARBIT
      COMMON /FOTON/ ARBIT(4)
      REAL *8  CMSENE,AMEL
      COMMON / BHPAR2 / CMSENE,AMEL
      REAL *8 ZER
      COMMON /HELP3/ C,C1,CPH1,CPH2,CPH3,ZER(4)
      DIMENSION SEL (2,2),SPO(2,2),SFI(2,2),SPH1(2,2),SPH2(2,2),
     $          SP(2,2),SA(2,2),SB(2,2),S(2,2),SC(2,2),SPH3(2,2),
     $          SS(2,2),SSA(2,2),SSB(2,2),SSS(2,2),
     $          SSSA(2,2),SSSB(2,2),SSSC(2,2),SSSD(2,2)
      REAL *8 R1A,R1B,R2A,R2B,R2C,R3A,R3B,R3C,R3D
      DIMENSION R1A(4),R1B(4),R2A(4),R2B(4),R2C(4),R3A(4),R3B(4),
     $          R3C(4),R3D(4)
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
      CALL MULTI(  1, R2A, R2A,C, C, LE3,ARBIT ,PH3 ,CPH3 , CPH3,SFI,C)
      CALL ILOCZ(SPH2,SFI,SA)
      CALL MULTI( LE2,PH2,ARBIT,C, C,1  , R2B,R2B,C , C,SPH2,-CR2)
      CALL MULTI(  1, R2B, R2B,C, C, LE3,ARBIT ,PH3 ,CPH3 , CPH3,SFI,C)
      CALL ILOCZ(SPH2,SFI,SB)
      CALL MULTI( LE2,PH2,ARBIT,C, C,1  , R2C,R2C,C , C,SPH2,-CR2)
      CALL MULTI(  1, R2C, R2C,C, C, LE3,ARBIT ,PH3 ,CPH3 , CPH3,SFI,C)
      CALL ILOCZ(SPH2,SFI,SC)
      CALL DODAJ(SA,SB,S)
      CALL DODAJ(S,SC,S)
*THIRD  PROPAGATOR
      CALL MULTI( LE3,PH3,ARBIT,C, C,1  , R3A,R3A,C , C,SPH3,CR3)
      CALL MULTI(  1, R3A, R3A,C, C, L3,P2 ,Q2 ,C  , C,SFI,C)
      CALL ILOCZ(SPH3,SFI,SSSA)
      CALL MULTI( LE3,PH3,ARBIT,C, C,1  , R3B,R3B,C , C,SPH3,-CR3)     
      CALL MULTI(  1, R3B, R3B,C, C, L3,P2 ,Q2 ,C  , C,SFI,C)
      CALL ILOCZ(SPH3,SFI,SSSB)
      CALL MULTI( LE3,PH3,ARBIT,C, C,1  , R3C,R3C,C , C,SPH3,-CR3)
      CALL MULTI(  1, R3C, R3C,C, C, L3,P2 ,Q2 ,C  , C,SFI,C)
      CALL ILOCZ(SPH3,SFI,SSSC)
      CALL MULTI( LE3,PH3,ARBIT,C, C,1  , R3D,R3D,C , C,SPH3,-CR3)
      CALL MULTI(  1, R3D, R3D,C, C, L3,P2 ,Q2 ,C  , C,SFI,C)
      CALL ILOCZ(SPH3,SFI,SSSD)
      CALL DODAJ(SSSA,SSSB,SSS)
      CALL DODAJ(SSS,SSSC,SSS)
      CALL DODAJ(SSS,SSSD,SSS)
*Q1 LINE
      CALL MULTI( L4,Q2 ,P2 ,C, C, L2,Q1 ,ZER,C ,C , SPO ,C  )
      CALL ADD5(SEL,SS,S,SSS,SPO,X)

ccc   only infrared part
c      CALL ADD5(SEL,SSA,SA,SSSA,SPO,X)

      END
      SUBROUTINE TRZYLEC(L1,L2,L3,L4,LE1,LE2,LE3,R1A,R1B,CR1,
     %                  R2A,R2B,R2C,CR2,R3A,R3B,CR3,X)
*     **************************************************
      IMPLICIT NONE
      INTEGER L1, L2, L3, L4, LE1, LE2, LE3
      COMPLEX *16   CR1,CR2,CR3,X,X1,X2
      COMPLEX *16 C1,C,CPH1,CPH2,CPH3
      COMPLEX *16   SEL ,SPH1,SP,SPH2,SFI,SPO,SA,SB,S,SS,SSA,SSB,
     $              SPH3, SSS,SC,SSSA,SSSB,SSSC
      REAL *8 P1,Q1,P2,Q2,PK1,PK2,PK3
      COMMON / MOMS3 / P1(4),Q1(4),P2(4),Q2(4),PK1(4),PK2(4),PK3(4)
      REAL *8 PH1,PH2,PH,PH3,PH123,PH23
      COMMON / CURR3 / PH1(4),PH2(4),PH(4),PH3(4),PH123(4),PH23(4)
      REAl *8 ARBIT
      COMMON /FOTON/ ARBIT(4)
      REAL *8  CMSENE,AMEL
      COMMON / BHPAR2 / CMSENE,AMEL
      REAL *8 ZER
      COMMON /HELP3/ C,C1,CPH1,CPH2,CPH3,ZER(4)
      DIMENSION SEL (2,2),SPO(2,2),SFI(2,2),SPH1(2,2),SPH2(2,2),
     $          SP(2,2),SA(2,2),SB(2,2),S(2,2),SC(2,2),
     $          SS(2,2),SSA(2,2),SSB(2,2),SPH3(2,2),SSS(2,2),
     $          SSSA(2,2),SSSB(2,2),SSSC(2,2)
      REAL *8 R1A,R1B,R2A,R2B,R2C,R3A,R3B
      DIMENSION R1A(4),R1B(4),R2A(4),R2B(4),R2C(4),R3A(4),R3B(4)

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
      CALL ILOCZ(SPH2,SFI,SSSA)
      CALL MULTI( LE2,PH2,ARBIT,C, C,1  , R2B,R2B,C , C,SPH2,-CR2)
      CALL MULTI(  1, R2B, R2B,C, C, L3,P2 ,Q2 ,C  , C,SFI,C1)
      CALL ILOCZ(SPH2,SFI,SSSB)
      CALL MULTI( LE2,PH2,ARBIT,C, C,1  , R2C,R2C,C , C,SPH2,-CR2)
      CALL MULTI(  1, R2C, R2C,C, C, L3,P2 ,Q2 ,C  , C,SFI,C1)
      CALL ILOCZ(SPH2,SFI,SSSC)
      CALL DODAJ(SSSA,SSSB,SSS)
      CALL DODAJ(SSS,SSSC,SSS)
*THIRD PROPAGATOR
      CALL MULTI( L4 ,Q2 ,P2,C ,C ,1  , R3A,R3A,C , C,SP  ,-CR3)
      CALL MULTI(  1, R3A, R3A,C, C,LE3,ARBIT,PH3,CPH3,CPH3,SPH3,C)
      CALL ILOCZ(SP  ,SPH3,SA)
      CALL MULTI( L4 ,Q2 ,P2,C ,C, 1  , R3B,R3B,C , C,SP  , CR3)
      CALL MULTI(  1, R3B, R3B,C, C,LE3,ARBIT,PH3,CPH3,CPH3,SPH3,C)
      CALL ILOCZ(SP  ,SPH3,SB)
      CALL DODAJ(SA,SB,S)
*Q1 LINE
      CALL MULTI(LE3,PH3,ARBIT,C ,C,L2,Q1 ,ZER,C  , C,SPO,C )

      CALL ADD5(SEL ,SS,SSS,S,SPO,X)

ccc only infrared part
c     CALL ADD5(SEL ,SSA,SSSA,SB,SPO,X)
 
      END
      SUBROUTINE TRZYLED(L1,L2,L3,L4,LE1,LE2,LE3,R1A,R1B,CR1,
     #                    R2A,R2B,R2C,CR2,R3A,R3B,CR3,X)
*     *******************************************************
      IMPLICIT NONE
      INTEGER L1, L2, L3, L4, LE1, LE2, LE3
      COMPLEX *16   CR1,CR2,CR3,X,X1,X2
      COMPLEX *16   C1,C,CPH1,CPH2,CPH3
      COMPLEX *16   SEL ,SPH1,SP,SPH2,SFI,SPO,SA,SB,S,SS,
     $              SSA,SSB,SSS,SC,SPH3,SSSA,SSSB
      REAL *8 P1,Q1,P2,Q2,PK1,PK2,PK3
      COMMON / MOMS3 / P1(4),Q1(4),P2(4),Q2(4),PK1(4),PK2(4),PK3(4)
      REAL *8 PH1,PH2,PH,PH3,PH123,PH23
      COMMON / CURR3 / PH1(4),PH2(4),PH(4),PH3(4),PH123(4),PH23(4)
      REAL *8  CMSENE,AMEL
      COMMON / BHPAR2 / CMSENE,AMEL
      REAL *8 ZER
      COMMON /HELP3/ C,C1,CPH1,CPH2,CPH3,ZER(4)
      REAl *8 ARBIT
      COMMON /FOTON/ ARBIT(4)
      DIMENSION SEL (2,2),SPO(2,2),SFI(2,2),SPH1(2,2),SPH2(2,2),
     $          SP(2,2),SA(2,2),SB(2,2),S(2,2),SC(2,2),SPH3(2,2),
     $          SS(2,2),SSA(2,2),SSB(2,2),SSS(2,2),SSSA(2,2),
     $          SSSB(2,2)
      REAL *8 R1A,R1B,R2A,R2B,R2C,R3A,R3B
      DIMENSION R1A(4),R1B(4),R2A(4),R2B(4),R2C(4),R3A(4),R3B(4)
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
      CALL MULTI( L4 ,Q2 ,P2,C ,C ,1  , R2A,R2A,C , C,SP  ,-CR2)
      CALL MULTI(  1, R2A, R2A,C, C,LE2,ARBIT,PH2,CPH2,CPH2,SPH2,C)
      CALL ILOCZ(SP  ,SPH2,SA)
      CALL MULTI( L4 ,Q2 ,P2,C ,C, 1  , R2B,R2B,C , C,SP  ,-CR2)
      CALL MULTI(  1, R2B, R2B,C, C,LE2,ARBIT,PH2,CPH2,CPH2,SPH2,C)
      CALL ILOCZ(SP  ,SPH2,SB)
      CALL MULTI( L4 ,Q2 ,P2,C ,C, 1  , R2C,R2C,C , C,SP  , CR2)
      CALL MULTI(  1, R2C, R2C,C, C,LE2,ARBIT,PH2,CPH2,CPH2,SPH2,C)
      CALL ILOCZ(SP  ,SPH2,SC)
      CALL DODAJ(SA,SB,S)
      CALL DODAJ(S,SC,S)
*THIRD PROPAGATOR
      CALL MULTI( LE2,PH2,ARBIT,C ,C ,1  , R3A,R3A,C , C,SP  ,-CR3)
      CALL MULTI(  1, R3A, R3A,C, C,LE3,ARBIT,PH3,CPH3,CPH3,SPH3,C)
      CALL ILOCZ(SP  ,SPH3,SSSA)
      CALL MULTI( LE2,PH2,ARBIT,C ,C, 1  , R3B,R3B,C , C,SP  , CR3)
      CALL MULTI(  1, R3B, R3B,C, C,LE3,ARBIT,PH3,CPH3,CPH3,SPH3,C)
      CALL ILOCZ(SP  ,SPH3,SSSB)
      CALL DODAJ(SSSA,SSSB,SSS)
*Q1 LINE
      CALL MULTI(LE3,PH3,ARBIT,C ,C,L2,Q1 ,ZER,C  , C,SPO,C )
 
      CALL ADD5(SEL ,SS,S,SSS,SPO,X)

ccc only infrared part
c      CALL ADD5(SEL ,SSA,SC,SSSB,SPO,X)
 
      END
      SUBROUTINE TRZYLEB(L1,L2,L3,L4,LE1,LE2,LE3,R1A,R1B,R1C,R1D,CR1,
     %                 R2A,R2B,R2C,CR2,R3A,R3B,CR3,X)
*     **************************************************
      IMPLICIT NONE
      INTEGER L1, L2, L3, L4, LE1, LE2, LE3
      COMPLEX *16   CR1,CR2,CR3,X,X1,X2
      COMPLEX *16 C1,C,CPH1,CPH2,CPH3
      COMPLEX *16   SEL ,SPH1,SP,SPH2,SFI,SPO,SA,SB,S,SS,SSA,
     $              SPH3,SSS,SSB,SSC,SSSA,SSSB,SSSC,SSSD
      REAL *8 P1,Q1,P2,Q2,PK1,PK2,PK3
      COMMON / MOMS3 / P1(4),Q1(4),P2(4),Q2(4),PK1(4),PK2(4),PK3(4)
      REAL *8 PH1,PH2,PH,PH3,PH123,PH23
      COMMON / CURR3 / PH1(4),PH2(4),PH(4),PH3(4),PH123(4),PH23(4)
      REAL *8  CMSENE,AMEL
      COMMON / BHPAR2 / CMSENE,AMEL
      REAl *8 ARBIT
      COMMON /FOTON/ ARBIT(4)
      REAL *8 ZER
      COMMON /HELP3/ C,C1,CPH1,CPH2,CPH3,ZER(4)
      DIMENSION SEL (2,2),SPO(2,2),SFI(2,2),SPH1(2,2),SPH2(2,2),
     $          SPH3(2,2),SP(2,2),SA(2,2),SB(2,2),S(2,2),
     $          SS(2,2),SSA(2,2),SSB(2,2),SSC(2,2),SSS(2,2),
     $          SSSA(2,2),SSSB(2,2),SSSC(2,2),SSSD(2,2)
      REAL *8 R1A,R1B,R1C,R1D,R2A,R2B,R2C,R3A,R3B
      DIMENSION R1A(4),R1B(4),R1C(4),R1D(4),R2A(4),R2B(4),R2C(4),
     $          R3A(4),R3B(4)
*P1 LIN2
      CALL MULTI( L1,P1 ,ZER,C,-C,L3 ,P2 ,Q2,C , C ,SEL,C1)
*FIRST PROPAGATOR
      CALL MULTI( L4,Q2 ,P2 ,C, C,1  ,R1A,R1A ,C , C ,SFI ,-CR1)
      CALL MULTI(  1, R1A, R1A,C, C,LE1,ARBIT,PH1,CPH1 , CPH1 ,SPH1,C )
      CALL ILOCZ(SFI ,SPH1,SSSA)
      CALL MULTI( L4,Q2 ,P2 ,C, C,1  ,R1B ,R1B ,C , C ,SFI ,-CR1)
      CALL MULTI(  1, R1B, R1B,C, C,LE1,ARBIT,PH1,CPH1 , CPH1 ,SPH1,C )
      CALL ILOCZ(SFI ,SPH1,SSSB)
      CALL DODAJ(SSSA,SSSB,SSS)
      CALL MULTI( L4,Q2 ,P2 ,C, C,1  ,R1C ,R1C ,C , C ,SFI ,-CR1)
      CALL MULTI(  1, R1C, R1C,C, C,LE1,ARBIT,PH1,CPH1 , CPH1 ,SPH1,C )
      CALL ILOCZ(SFI ,SPH1,SSSC)
      CALL DODAJ(SSS,SSSC,SSS)
      CALL MULTI( L4,Q2 ,P2 ,C, C,1  ,R1D ,R1D ,C , C ,SFI , CR1)
      CALL MULTI(  1, R1D, R1D,C, C,LE1,ARBIT,PH1,CPH1 , CPH1 ,SPH1,C )
      CALL ILOCZ(SFI ,SPH1,SSSD)
      CALL DODAJ(SSS,SSSD,SSS)
*SECOND PROPAGATOR
      CALL MULTI( LE1,PH1,ARBIT ,C, C,1  ,R2A,R2A ,C , C ,SFI ,-CR2)
      CALL MULTI(  1, R2A, R2A,C, C,LE2,ARBIT,PH2,CPH2 , CPH2 ,SPH2,C )
      CALL ILOCZ(SFI ,SPH2,SSA)
      CALL MULTI( LE1,PH1,ARBIT ,C, C,1  ,R2B ,R2B ,C , C ,SFI , -CR2)
      CALL MULTI(  1, R2B, R2B,C, C,LE2,ARBIT,PH2,CPH2 , CPH2 ,SPH2,C )
      CALL ILOCZ(SFI ,SPH2,SSB)
      CALL DODAJ(SSA,SSB,SS)
      CALL MULTI( LE1,PH1,ARBIT ,C, C,1  ,R2C ,R2C ,C , C ,SFI , CR2)
      CALL MULTI(  1, R2C, R2C,C, C,LE2,ARBIT,PH2,CPH2 , CPH2 ,SPH2,C )
      CALL ILOCZ(SFI ,SPH2,SSC)
      CALL DODAJ(SS,SSC,SS)
*THIRD PROPAGATOR
      CALL MULTI( LE2,PH2,ARBIT ,C, C,1  , R3A,R3A,C , C,SP  ,-CR3)
      CALL MULTI(  1, R3A, R3A,C, C,LE3,ARBIT,PH3,CPH3 , CPH3 ,SPH3,C) 
      CALL ILOCZ(SP  ,SPH3,SA)
      CALL MULTI( LE2,PH2,ARBIT ,C, C,1  , R3B,R3B,C , C,SP  , CR3)
      CALL MULTI(  1, R3B, R3B,C, C,LE3,ARBIT,PH3,CPH3 , CPH3 ,SPH3,C)
      CALL ILOCZ(SP  ,SPH3,SB)
      CALL DODAJ(SA,SB,S)
*Q1 LINE
      CALL MULTI(LE3, PH3,ARBIT,C, C, L2,Q1 ,ZER,C  , C,SPO,C )
 
      CALL ADD5(SEL,SSS,SS,S,SPO,X)

ccc  infrared part only
c      CALL ADD5(SEL,SSSD,SSC,SB,SPO,X)
 
      END
      SUBROUTINE TRZIINIINF(SECTION)
*     ***********************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 XMTINIINF,S,SUM,C
      COMMON / MOMS3 / P1(4),Q1(4),P2(4),Q2(4),PK1(4),PK2(4),PK3(4)
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
      DO 10 NN=1,3,2
      LEPS3=2-NN
      S= XMTINIINF( LAM1,LAM2,LAM3,LAM4,LEPS1,LEPS2,LEPS3)
      S=S*DCONJG(S)
      SUM=SUM+S 

  10  CONTINUE
*CROSS SECTION  
      SECTION=DBLE(SUM)
      END
 
 

*SPIN AMPLITUDE  INITIAL TRIPLE
      FUNCTION XMTINIINF(L1,L2,L3,L4,LE1,LE2,LE3)
*     ******************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16  HIS ,XMTINIINF,XP1P1P1,XP1P1Q1,XP1Q1Q1,XQ1Q1Q1
      COMPLEX *16 CPH1,CPH2,CPH3,C,C1,CR1,CR2,CR3,S
      COMMON /HELP3/ C,C1,CPH1,CPH2,CPH3,ZER(4)
      COMMON / MOMS3 / P1(4),Q1(4),P2(4),Q2(4),PK1(4),PK2(4),PK3(4)
      COMMON / CURR3 / PH1(4),PH2(4),PH(4),PH3(4),PH123(4),PH23(4)
      COMMON /FOTON/ ARBIT(4)
      COMPLEX X1PK1,X2PK1,X1PK2,X2PK2,X1PK3,X2PK3
      COMPLEX X1, X2
 
      S1= 2D0*(P2(4)*Q2(4)-P2(3)*Q2(3)-P2(2)*Q2(2)-P2(1)*Q2(1))
 
      C=DCMPLX(1D0,0D0)
      C1=C
      S=(0D0,0D0)
      X1=(0D0,0D0)
      X2=(0D0,0D0)
      X1PK1=(0D0,0D0)
      X2PK1=(0D0,0D0)
      X1PK2=(0D0,0D0)
      X2PK2=(0D0,0D0)
      X1PK3=(0D0,0D0)
      X2PK3=(0D0,0D0)

      DO 10  KK=1,6
      IF (KK.EQ.1) THEN
      npk1=1.
      npk2=2.
      npk3=3.
      DO II=1,4
       PH1(II)=PK1(II)
       PH2(II)=PK2(II)
       PH3(II)=PK3(II)
      ENDDO
      LL1=-LE1
      LL2=-LE2     
      LL3=-LE3
      ELSEIF(KK.EQ.2) THEN
      npk1=1.
      npk2=3.
      npk3=2.
      DO II=1,4
       PH1(II)=PK1(II)
       PH2(II)=PK3(II)
       PH3(II)=PK2(II)
      ENDDO
      LL1=-LE1
      LL2=-LE3
      LL3=-LE2
      ELSEIF(KK.EQ.3) THEN
      npk1=3.
      npk2=2.
      npk3=1.
      DO II=1,4
       PH1(II)=PK3(II)
       PH2(II)=PK2(II)
       PH3(II)=PK1(II)
      ENDDO
      LL1=-LE3
      LL2=-LE2
      LL3=-LE1
      ELSEIF(KK.EQ.4) THEN
      npk1=2.
      npk2=3.
      npk3=1.
      DO II=1,4
       PH1(II)=PK3(II)
       PH2(II)=PK1(II)
       PH3(II)=PK2(II)
      ENDDO
      LL1=-LE3
      LL2=-LE1
      LL3=-LE2
      ELSEIF(KK.EQ.5) THEN
      npk1=3.
      npk2=1.
      npk3=2.
      DO II=1,4
       PH1(II)=PK2(II)
       PH2(II)=PK3(II)
       PH3(II)=PK1(II)
      ENDDO
      LL1=-LE2
      LL2=-LE3
      LL3=-LE1
      ELSEIF(KK.EQ.6) THEN
      npk1=2.
      npk2=1.
      npk3=3.
      DO II=1,4
       PH1(II)=PK2(II)
       PH2(II)=PK1(II)
       PH3(II)=PK3(II)
      ENDDO
      LL1=-LE2
      LL2=-LE1
      LL3=-LE3
      ENDIF



      DO 15 I=1,4
      PH(I)   =PH1(I)+PH2(I)
      PH123(I)=PH1(I)+PH2(I) +PH3(I)
      PH23(I) =PH2(I)+PH3(I)
 15   ZER(I)  =0D0
      PHPH=PH(4)**2-PH(1)**2-PH(2)**2-PH(3)**2
      PHPH123=PH123(4)**2-PH123(1)**2-PH123(2)**2-PH123(3)**2
 
      CPH1=DSQRT(-PROPIN(PH1,ARBIT)*2D0)*C
      CPH2=DSQRT(-PROPIN(PH2,ARBIT)*2D0)*C
      CPH3=DSQRT(-PROPIN(PH3,ARBIT)*2D0)*C
*--------------------------------------------------------
*SPIN AMPLITUDE - THREE  BREMSTRAHLUNGS FROM LINE P1
      CR1=PROPIN(P1,PH1)*C
      CR2=PROPININFRA(P1,PH)*C
      CR3=PROPININFRA(P1,PH123)*C
 
      CALL TRZYLEA(L1,L2,L3,L4,LL1,LL2,LL3,P1,zer,CR1,P1,zer,zer,
     $           CR2,P1,zer,zer,zer,CR3,XP1P1P1)

      if(npk1.eq.1)X1PK1=X1PK1+XP1P1P1
      if(npk1.eq.2)X1PK1=X1PK1+XP1P1P1
      if(npk1.eq.3)X1PK1=X1PK1+XP1P1P1
      if(npk2.eq.1)X1PK2=X1PK2+XP1P1P1
      if(npk2.eq.2)X1PK2=X1PK2+XP1P1P1
      if(npk2.eq.3)X1PK2=X1PK2+XP1P1P1
      if(npk3.eq.1)X1PK3=X1PK3+XP1P1P1
      if(npk3.eq.2)X1PK3=X1PK3+XP1P1P1
      if(npk3.eq.3)X1PK3=X1PK3+XP1P1P1

*--------------------------------------------------------
*SPIN AMPLITUDE - TWO BREMSTRAHLUNGS FROM LINE P1, ONE FROM LINE Q1
      CR1=PROPIN(P1,PH1)*C
      CR2=PROPININFRA(P1,PH) *C
      CR3=PROPIN(Q1,PH3)*C
 
      CALL TRZYLEC(L1,L2,L3,L4,LL1,LL2,LL3,P1,zer,CR1,P1,zer,zer,
     $           CR2,zer,Q1,CR3,XP1P1Q1)

      if(npk1.eq.1)X1PK1=X1PK1-XP1P1Q1
      if(npk1.eq.2)X1PK1=X1PK1-XP1P1Q1
      if(npk1.eq.3)X2PK1=X2PK1-XP1P1Q1
      if(npk2.eq.1)X1PK2=X1PK2-XP1P1Q1
      if(npk2.eq.2)X1PK2=X1PK2-XP1P1Q1
      if(npk2.eq.3)X2PK2=X2PK2-XP1P1Q1
      if(npk3.eq.1)X1PK3=X1PK3-XP1P1Q1
      if(npk3.eq.2)X1PK3=X1PK3-XP1P1Q1
      if(npk3.eq.3)X2PK3=X2PK3-XP1P1Q1

*--------------------------------------------------------
*SPIN AMPLITUDE - ONE BREMS FROM LINE P1,TWO FROM LINE Q1
      CR1=PROPIN(P1,PH1)*C
      CR2=PROPININFRA(Q1,PH23)*C
      CR3=PROPIN(Q1,PH3)*C
 
      CALL TRZYLED(L1,L2,L3,L4,LL1,LL2,LL3,P1,zer,CR1,zer,zer,Q1,
     $             CR2,zer,Q1,CR3,XP1Q1Q1)

      if(npk1.eq.1)X1PK1=X1PK1+XP1Q1Q1
      if(npk1.eq.2)X2PK1=X2PK1+XP1Q1Q1
      if(npk1.eq.3)X2PK1=X2PK1+XP1Q1Q1
      if(npk2.eq.1)X1PK2=X1PK2+XP1Q1Q1
      if(npk2.eq.2)X2PK2=X2PK2+XP1Q1Q1
      if(npk2.eq.3)X2PK2=X2PK2+XP1Q1Q1
      if(npk3.eq.1)X1PK3=X1PK3+XP1Q1Q1
      if(npk3.eq.2)X2PK3=X2PK3+XP1Q1Q1
      if(npk3.eq.3)X2PK3=X2PK3+XP1Q1Q1

*------------------------------------------------------    
*SPIN AMPLITUDE - THREE BREMSTRAHLUNGS FROM LINE Q1
      CR1=PROPININFRA(Q1,PH123)*C
      CR2=PROPININFRA(Q1,PH23)*C
      CR3=PROPIN(Q1,PH3)*C
 
      CALL TRZYLEB(L1,L2,L3,L4,LL1,LL2,LL3,zer,zer,zer,Q1,CR1,zer,
     $            zer,Q1,CR2,zer,Q1,CR3,XQ1Q1Q1)

      if(npk1.eq.1)X2PK1=X2PK1-XQ1Q1Q1
      if(npk1.eq.2)X2PK1=X2PK1-XQ1Q1Q1
      if(npk1.eq.3)X2PK1=X2PK1-XQ1Q1Q1
      if(npk2.eq.1)X2PK2=X2PK2-XQ1Q1Q1
      if(npk2.eq.2)X2PK2=X2PK2-XQ1Q1Q1
      if(npk2.eq.3)X2PK2=X2PK2-XQ1Q1Q1
      if(npk3.eq.1)X2PK3=X2PK3-XQ1Q1Q1
      if(npk3.eq.2)X2PK3=X2PK3-XQ1Q1Q1
      if(npk3.eq.3)X2PK3=X2PK3-XQ1Q1Q1

*----------------------==--------------------------------
*TOTAL SPIN AMPLITUDE
      S=S+ (XP1P1P1-XP1P1Q1+XP1Q1Q1-XQ1Q1Q1)
 
 10   CONTINUE

      XMTINIINF=S*HIS(L1,L3,S1)
cc        write(6,*)'>>>>>>>>>'
cc        write(6,*)'section',l1,l2,l3,l4,ll1,ll2,ll3,s*dconjg(s)
cc        write(6,*)'x1-pk1',x1pk1
cc        write(6,*)'x2-pk1',x2pk1
cc        write(6,*)'x1-pk2',x1pk2
cc        write(6,*)'x2-pk2',x2pk2
cc        write(6,*)'x1-pk3',x1pk3
cc        write(6,*)'x2-pk3',x2pk3
cc        write(6,*)'>>>>>>>>>'

 
      END

*INFRARED LIMIT FOR CROSS SECTION- INITIAL STATE 3 PHOTONS YFS FORMULA
      SUBROUTINE TRZIYFSINF(SECTION)
C     ******************************************
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / MOMS3 / P1(4),Q1(4),P2(4),Q2(4),PK1(4),PK2(4),PK3(4)
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
      CALL SFACH0(P1,Q1,PK3,SFACT3)
      SUM= SFACT1*SFACT2*SFACT3*BETA00
      SECTION=SUM*S0/S1
 
      END

*INFRARED LIMIT FOR CROSS SECTION- INITIAL STATE 2 PHOTNS SPIN AMPL.
*TIMES THIRD PHOTON INFRARED
      SUBROUTINE TRZIDBL(SECTION)
C     ******************************************
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / MOMS3 / P1(4),Q1(4),P2(4),Q2(4),PK1(4),PK2(4),PK3(4)
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
      CALL SFACH0(P1,Q1,PK3,SFACT3)
      CALL DUBLINI(SECT3)
      SUM= SFACT3*SECT3
      SECTION=SUM
 
      END

*INFRARED LIMIT FOR CROSS SECTION- INITIAL STATE 2 PHOTONS YFS TIMES
*THIRD PHOTON INFRARED
      SUBROUTINE TRZIYFS(SECTION)
C     ******************************************
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / MOMS3 / P1(4),Q1(4),P2(4),Q2(4),PK1(4),PK2(4),PK3(4)
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
      CALL SFACH0(P1,Q1,PK3,SFACT3)
      CALL DINIAPR(SECT3)
      SUM= SFACT3*SECT3
      SECTION=SUM
 
      END


