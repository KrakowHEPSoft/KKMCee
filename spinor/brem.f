*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//  Collection of old subprograms from KORALB and KORALZ                           //
*//                                                                                 //
*//  WARNING:                                                                       //
*//  Due to short common block names this code overwrites variables in other parts  //
*//  of the code.                                                                   //
*//  One should add suffix c_Brem_ to names of all commons as soon as possible!     //
*//  I have done it for /c_Brem_INOUT/  SJ.                                         //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////


      SUBROUTINE Brem_Porownanie(KFi,KFf,p1,p2,p3,p4,ph,oldist)
**    TEN PROGRAM WYLICZA ROZNICKOWY PRZEKROJ CZYNNY NA KILKA
**    SPOSOBOW. -4- TAK JAK KORALB KAWALEK NIESPOLARYZOWANY.
**              -2- Z AMPLITUD SPINOWYCH KORALB.
**              -1- Z MUSTRAALA FUNTIZ
**              -3- Z AMPLITUD SPINOWYCH Z DOKTORKI
!      PRINT *,'mustraal                              FUNTIZ= ', X  1
!      PRINT *,'spin ampl                             KORALB= ', Y  2
!      PRINT *,'doktorat  spin ampl - density m koralb ZAMPL= ', U  3
!      PRINT *,'KORALB                                 WINTH= ', T  4
**    JAKO PARAMETR JEST TEZ POLARYZACJA BEAMU ORAZ Z0 ON OFF
**    WYLICZAJA SIE ZAWSZE WSZYSKIE PRZEKROJE NA TYLE NA ILE MOGA
**    NP 2,3 UWZGLEDNIA CALE WEKTORY POLARYZACJI
**    1 TYLKO POLARYZACJE PDLOZNA
**    TEN PROGRAM REKONSTRUUJE KONTY Z KORALA-B Z CZTEROPEDOW
**    (W PROCEDURZE PEDYPR)
**    UWAGA ! W AMPLITUDACH Z KORALA B NANIESIONA POPRAWKA !!!!!!!!!!!!
**    DOTYCZY ONA WSPOLCZYNNIKOW W CZESCI SPIN FLIP INITIAL BREMSSTRAHLUNG
**    NIE MAJA ONE WPLYWU DOPUKI NIE UWZGLEDNIA SIE POLARYZACJI POPRZECZNEJ
**    BEAMU. POPRAWKI OZNACZONE:
**
*$ POPRAWKA
**
***************************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
CC    ************************
c[[[      DIMENSION E1(3),E2(3),h1(3),h2(3)
      DIMENSION E1(4),E2(4),h1(4),h2(4)
      dimension p1(4),p2(4),p3(4),p4(4),ph(4),oldist(4)
      real*4 rvec(10)
      common /breminprm/  xENE,xAMZ,xGAMM,xAMFIN,xsinw2,ixzonon
      common /ifinifin/  ifini,iffin
      common /ifzrun/ ifrunz
      DATA IFINIT /0/
**    ************************
      DO 10 K=1,3
      h1(k)=0
      h2(k)=0
      E1(K)=0.D0
 10   E2(K)=0.D0
      e1(1)=.2
       e1(2)=.3
       e1(3)=.5
      e2(1)=.5
      e2(2)=.3
      e2(3)=-.55
         CALL PseuMar_MakeVec(rvec,6)
       h1(1)=rvec(1)-.5
       h2(1)=rvec(2)-.5
       h1(2)=rvec(3)-.5
       h2(2)=rvec(4)-.5
       h1(3)=(rvec(5)-.5)                    !     *  .5
       h2(3)=(rvec(6)-.5)                    !     * 2
       hh1=sqrt(h1(1)**2+h1(2)**2+h1(3)**2)  !     * 1.1
       hh2=sqrt(h2(1)**2+h2(2)**2+h2(3)**2)  !     * 1.1
       do k=1,3
        h1(k)=h1(k)/hh1                      !      * 0
        h2(k)=h2(k)/hh2                      !      * 0
       enddo
      IF (IFINIT.EQ.0) THEN
       WRITE(16,*) 'POLARIZATION         E1               E2'
       WRITE(16,*) '1             ',E1(1),'   ',E2(1)
       WRITE(16,*) '2             ',E1(2),'   ',E2(2)
       WRITE(16,*) '3             ',E1(3),'   ',E2(3)
!       STOP
      ENDIF
      IFINIT=1
* Possibility to switch off Z
      CALL BornV_GetKeyZet(KeyZet)
      IxZONON=KeyZet
      ifini=1
      iffin=1
      ifrunz=1
c[[[[[[[[
c      CALL Brem_SpinStore(0,h1,h2)
c      CALL Brem_SpinBeam( 0,e1,e2)
c]]]]]]]]
      CALL STest_SetHvectors(h1,h2)
      CALL STest_SetPolBeams(e1,e2)

      CALL BornV_GetSwsq(  xSINW2 )
      CALL BornV_GetMZ(    xAMZ   )
      CALL BornV_GetGammZ( xGAMM)
      xamfin =  BornV_GetMass(KFf)
      xENE  =(p1(4)+p2(4))/2

      CALL BEGIN(E1,E2)
      CALL PEDYPR(p1,p2,p3,p4,ph)
      CALL PEDY
      CALL TEST1(E1,E2,oldist)
      do k=1,4
       oldist(k)=2*oldist(k)/xene**2
      enddo
      END


      SUBROUTINE BEGIN(E1,E2)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON / CONST / PI,ALFA,ALF1,QE2,QF2,QEF
      COMMON/c_Brem_INOUT/ NINP,NOUT
      COMMON / KEYGSW / KEYGSW
      COMMON / INIFIN/ ENE,SVAR,AEL2,AMU2,ALGEL,ALGMU,BETI,BETF
      COMMON / WEAK1 / ZZ,AMZ,GAMM,GM,GG,BW,ABW,ABW1,FI,DLR
      COMMON / WEAK  / QCE,QCF,CVE,CVF,CAE,CAF,ZPROP
      COMMON /INSPIN / SEPS1,SEPS2
      common /breminprm/  xENE,xAMZ,xGAMM,xAMFIN,xsinw2,ixzonon
      DIMENSION E1(3),E2(3)
**    ************************
**    SET BEAM POLARIZATIONS A. S. O.
!      PRINT *, 'E1?'
!      READ *, E1(1),E1(2),E1(3)
!      PRINT *, 'E2?'
!      READ *, E2(1),E2(2),E2(3)
      ENE=xENE
      AMZ=xAMZ
      GAMM=xGAMM
      AMFIN=xAMFIN
      sinw2=xsinw2
      izonon=ixzonon
c[[[      SEPS1=-E1(3)
c[[[      SEPS2=-E2(3)
      SEPS1= E1(3)
      SEPS2= E2(3)
**    ************************
      KEYGSW=0
 !     PRINT *, 'ENE AMZ GAMM AMFIN ?'
  !    READ  *,  ENE,AMZ,GAMM,AMFIN
 !     PRINT *, 'IZONON = ?        '
 !     READ  *,  IZONON
 !     PRINT *, 'SINW2  = ?        '
!      READ  *,  SINW2
**    SET COUPLINGS............
      IDE=-2
      IDF=-2
      CALL SETCUP(IDE,SINW2,QCE,CVE,CAE)
      CALL SETCUP(IDF,SINW2,QCF,CVF,CAF)
      NINP=5
      NOUT=16
      PI=4.D0*ATAN(1.D0)
      ALFA=1.D0/137.03604D0
      ALF1=ALFA/PI
      IF(IZONON.EQ.0) THEN
        CVE=0
        CAE=0
        CVF=0
        CAF=0
      ENDIF
C;;
      QE2=QCE**2
      QF2=QCF**2
***   QCF=0.D0
      QEF=QCE*QCF
      CALL COUPLER(IDF,CVE,CAE,CVF,CAF)
***   SOME PARAMETERS  CALCULATION
      SVAR=4.*ENE*ENE
      GM=GAMM*AMZ/SVAR
      GG=GM*GM
      ZZ=1.-AMZ**2/SVAR
      BW=ZZ*ZZ+GG
      AEL=.5111E-3/ENE
***   AMFIN=0.10566
***   AMFIN=AMFI1
      AMU=AMFIN/ENE
      AEL2=AEL**2
      AMU2=AMU**2
  !    PRINT *,'************************************************* '
  !    PRINT *, 'POLARIZATION VECTORS OF E+ AND E-'
  !    PRINT *, 'E1=', E1(1),' ',E1(2),' ',E1(3)
  !    PRINT *, 'E2=', E2(1),' ',E2(2),' ',E2(3)
  !    PRINT *,'************************************************* '
  !    PRINT *,'ENE=   ', ENE
  !    PRINT *,'AMFIN= ', AMFIN
  !    PRINT *,'AMZ=   ', AMZ
  !    PRINT *,'GAMM=  ', GAMM
  !    PRINT *,'SINW2= ', SINW2
  !    PRINT *,'IZONON=', IZONON
  !    PRINT *,'************************************************* '
      RETURN
      END

      SUBROUTINE TEST1(E1,E2,oldist)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /INSPIN / SEPS1,SEPS2
      DIMENSION E1(3),E2(3),oldist(4)
**    ************************
      X=FUNTIZ(2)/8.D0
      DUMM=0.D0
      Y=XRAL(E1,E2)/2.D0
      Z=Y-X

      T= WINTH(DUM)
      U= ZAMPL(E1,E2)/4.D0
!      PRINT *,'mustraal                              FUNTIZ= ', X
!      PRINT *,'spin ampl                             KORALB= ', Y
!      PRINT *,'doktorat  spin ampl - density m koralb ZAMPL= ', U
!      PRINT *,'KORALB                                 WINTH= ', T
      oldist(1)=x
      oldist(2)=y
      oldist(3)=u
      oldist(4)=t
      RETURN
      END



      SUBROUTINE PEDY
      IMPLICIT REAL*8(A-H,O-Z)
**    ************************
      COMMON / UTIL1 / XK,C1,S1,C2,S2,CF,SF,CG,SG,V
      COMMON / UTIL  / QP(4),QM(4),PH(4)
      COMMON / INIFIN/ ENE,SVAR,AEL2,AMU2,ALGEL,ALGMU,BETI,BETF
      COMMON / ENERG / EN1,AEL1,AMF1,AM1,ALGE1,ALGM1,BET1,BT1,ATH2
      REAL*8 PI1(4),PI2(4),PIX(4),PIY(4)
      EN1=ENE
      AEL1=AEL2
      AMF1=AMU2
      AM1=SQRT(AMU2)


!      TET1=0.3 D0/2.D0
!      TET2=0.09D0
!      FIJE=0.700    D0
!      FIDW=0.600    D0
!      xk=0.4
!      PRINT *, ' TET1,TET2 '
    !  READ *, TET1,TET2
 !     PRINT *, ' FIJE FIDW '
   !   READ *, FIJE,FIDW
 !     PRINT *, ' XK '
   !   READ *, XK



 !     C1=COS(TET1)
 !     S1=SIN(TET1)

  !    C2=COS(TET2)
  !    S2=SIN(TET2)

  !    CF=COS(FIJE)
  !    SF=SIN(FIJE)

!      CG=COS(FIDW)
!      SG=SIN(FIDW)


!      PRINT *, '************************************************'
!      PRINT *, 'XK= ',XK
!      PRINT *, 'C1= ',C1,' C2= ',C2
!      PRINT *, 'CF= ',CF,' SF= ',SF
!      PRINT *, 'CG= ',CG,' SG= ',SG
!      PRINT *, '************************************************'


      DO 420 I=1,3
      QP(I)=0.D0
  420 QM(I)=0.D0
      QP(4)=AM1
      CALL Brem_Tralor(1,QP)
      QM(4)=AM1
      CALL Brem_Tralor(2,QM)
      DO 430 I=1,4
  430 PH(I)=-QP(I)-QM(I)
      PH(4)=2.D0+PH(4)
      IF(XK.EQ.0.) PH(4)=0.D0
***********************************************************************

      RETURN
      END
      SUBROUTINE PEDYPR(p1,p2,p3,p4,phx)
      IMPLICIT REAL*8(A-H,O-Z)
**    ************************
      COMMON / UTIL1 / XK,C1,S1,C2,S2,CF,SF,CG,SG,V
      COMMON / UTIL  / QP(4),QM(4),PH(4)
      COMMON / INIFIN/ ENE,SVAR,AEL2,AMU2,ALGEL,ALGMU,BETI,BETF
      COMMON / ENERG / EN1,AEL1,AMF1,AM1,ALGE1,ALGM1,BET1,BT1,ATH2
      dimension p1(4),p2(4),p3(4),p4(4),phx(4)
      REAL*8 PI1(4),PI2(4),PIX(4),PIY(4)
!      PRINT *, '************************************************'
!      PRINT *, 'XK= ',XK
!      PRINT *, 'C1= ',C1,' C2= ',C2
!      PRINT *, 'CF= ',CF,' SF= ',SF
!      PRINT *, 'CG= ',CG,' SG= ',SG
!      PRINT *, '************************************************'
      do k=1,4
       qp(k)=p3(k)/ene
       qm(k)=p4(k)/ene
       ph(k)=phx(k)/ene
      enddo
!      write(*,*) 'pedypr'
!      write(*,*) qp
!      write(*,*) qm
!      write(*,*) ph
!      write(*,*) qp
***********************************************************************
!      PRINT *, 'KONTY Z PEDOW '
      CC1=PH(3)/PH(4)
***   CALCULATION OF C2 ---------------------------
      CALL MULSK(PH,QM,R)
      CALL MULSK(PH,QP,R1)
      am1=BornV_GetMass(13)/ene
      CC2=(R-R1)/SQRT(1.D0-PH(4)-AM1*AM1)/2.D0/PH(4)*SQRT(1.D0-PH(4))
***   CALCULATION OF CG , SG -----------------------
      CCG= PH(2)/SQRT(PH(1)**2+PH(2)**2)
      CSG=-PH(1)/SQRT(PH(1)**2+PH(2)**2)
***   CALCULATION OF CF _---------------------------
      PIX(1)=0.D0
      PIX(2)=0.D0
      PIX(3)=1.D0
      PIX(4)=0.D0

      CALL MULVX(PH,QP,PI1)
      CALL MULVX(PH,PIX,PI2)
      CALL MULSK(PI1,PI2,R)
      CALL MULSK(PI1,PI1,R1)
      CALL MULSK(PI2,PI2,R2)
      CCF=R/SQRT(R1*R2)
****  CALCULATION OF SF ---------------------------
      CALL MULSK3(PH,QP,R)
      CALL MULSK3(PH,PIX,R1)

      CALL MULSK3(PH,PH,R3)
      DO 9 KKK=1,4
      PI1(KKK)=QP(KKK)-R/R3*PH(KKK)
   9  PI2(KKK)=PIX(KKK)-R1/R3*PH(KKK)

      CALL MULVX(PI1,PI2,PIY)
      CALL MULSK3(PH,PIY,R)
      CALL MULSK3(PI1,PI1,R1)
      CALL MULSK3(PI2,PI2,R2)
      CSF=-R/SQRT(-R1*R2*R3)
 !     PRINT *, 'KONTY Z PEDOW '
 !     PRINT *, '************************************************'
 !     PRINT *, 'XK= ',PH(4)
 !     PRINT *, 'C1= ',CC1,' C2= ',CC2
 !     PRINT *, 'CF= ',CCF,' SF= ',CSF
 !     PRINT *, 'CG= ',CCG,' SG= ',CSG
 !     PRINT *, '************************************************'
C
      c1=cc1
      s1=sqrt(1d0-c1**2)
      c2=cc2
      s2=sqrt(1d0-c2**2)
      cf=ccf
      sf=csf
      cg=ccg
      sg=csg
      xk=ph(4)
      RETURN
      END

      SUBROUTINE MULVX(X,Y,R)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(4),Y(4),R(4)
      R(4)=0.D0
      R(1)=X(2)*Y(3)-X(3)*Y(2)
      R(2)=X(3)*Y(1)-X(1)*Y(3)
      R(3)=X(1)*Y(2)-X(2)*Y(1)
      RETURN
      END
      SUBROUTINE MULSK(X,Y,R)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(4),Y(4)
      R=0.D0
      R=R+X(4)*Y(4)-X(3)*Y(3)
      R=R-X(2)*Y(2)-X(1)*Y(1)
      RETURN
      END
      SUBROUTINE MULSK3(X,Y,R)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(4),Y(4)
      R=0.D0
      R=R-X(3)*Y(3)
      R=R-X(2)*Y(2)-X(1)*Y(1)
      RETURN
      END

      FUNCTION FUNTIZ(IT)
C     *******************
      IMPLICIT REAL*8(A-H,O-Z)
      common /ifinifin/  ifini,iffin
      COMMON / INIFIN/ ENE,SVAR,AEL2,AMU2,ALGEL,ALGMU,BETI,BETF
      COMMON / CONST / PI,ALFA,ALF1,QE2,QF2,QEF
      COMMON / WEAK1 / ZZ,AMZ,GAMM,GM,GG,BW,ABW,ABW1,FI,DLR
      COMMON / WEAK2 / C0,C1,C2,D0,D1,D2,C3,C4,D3,D4
      COMMON / UTIL  / QP(4),QM(4),PH(4)
      COMMON / UTIL2 / XK0,XKMIN,XKMAX
      COMMON / UTIL3 / XK,C,S,CG,SG
      COMPLEX*16 ZAT1,ZAT0,PAT0,PAT1,CPRZ0,PIAA,CVPAA
      DATA ICONT/0/
      ICONT=ICONT+1
C

C

      T=2.*(QP(4)-QP(3))
      U=2.*(QM(4)-QM(3))
      T1=2.*(QM(4)+QM(3))
      U1=2.*(QP(4)+QP(3))
      Y1=PH(4)-PH(3)*(1.D0-0.5D0*AEL2)
      Y2=PH(4)+PH(3)*(1.D0-0.5D0*AEL2)
      Y3=2.*(1.-QM(4))
      Y4=2.*(1.-QP(4))
CC    INTERFERENCE CONTRIBUTION
      W=T*Y2*Y4+T1*Y1*Y3-U*Y2*Y3-U1*Y1*Y4
      W=W/(Y1*Y2*Y3*Y4)
      XK=PH(4)
      XKM=1.-XK
      TT=T*T+T1*T1+U*U+U1*U1
      VV=.5*(U*U+U1*U1-T*T-T1*T1)
      TT0=C0*TT+D0*VV
      TT1=C1*TT+D1*VV
      TT2=C2*TT+D2*VV
      D12=1.-.5*AEL2/(XKM+1./XKM)*(Y1/Y2+Y2/Y1)
      D34=1.-.5*AMU2/XKM/(XKM+1./XKM)*(Y3/Y4+Y4/Y3)

      BWK=(ZZ-XK)*(ZZ-XK)+GG
      BW0=BW
      SVA1=SVAR*XKM
      PAT0=SVAR/DCMPLX(SVAR,0D0)
      PAT1=SVAR/DCMPLX(SVA1,0D0)
      ZAT0=SVAR*CPRZ0(SVAR)
      ZAT1=SVAR*CPRZ0(SVA1)
      A1A1  = CDABS(PAT1)**2
      A1Z1RE=  REAL(PAT1*DCONJG(ZAT1))
      Z1Z1  = CDABS(ZAT1)**2
      A0A0  = CDABS(PAT0)**2
      A0Z0RE=  REAL(PAT0*DCONJG(ZAT0))
      Z0Z0  = CDABS(ZAT0)**2
      A0A1RE=  REAL(PAT0*DCONJG(PAT1))
      Z0A1RE=  REAL(ZAT0*DCONJG(PAT1))
      A0Z1RE=  REAL(PAT0*DCONJG(ZAT1))
      Z0Z1RE=  REAL(ZAT0*DCONJG(ZAT1))
      A0Z1IM= DIMAG(PAT0*DCONJG(ZAT1))
      Z0A1IM= DIMAG(ZAT0*DCONJG(PAT1))
      Z0Z1IM= DIMAG(ZAT0*DCONJG(ZAT1))
      ASQI=D12/Y1/Y2*XKM*(TT0*A1A1  +TT1*A1Z1RE         +TT2*Z1Z1)
      ASQF=D34/Y3/Y4*    (TT0*A0A0  +TT1*A0Z0RE         +TT2*Z0Z0)
      AIN=.25*W*(TT0*A0A1RE+.5*TT1*(Z0A1RE+A0Z1RE)+TT2*Z0Z1RE)
      XIN=(T*T-T1*T1)*(C3*(A0Z1IM+Z0A1IM)+2.*C4*Z0Z1IM)
     $   -(U*U-U1*U1)*(D3*(A0Z1IM+Z0A1IM)+2.*D4*Z0Z1IM)
      XIN=XIN*2.*XK/(Y1*Y2*Y3*Y4)*(QP(1)*QM(2)-QP(2)*QM(1))
      FUNTI= QE2*ifini*ASQI+QF2*iffin*ASQF +QEF*(AIN-XIN)*ifini*iffin
      FUNTIZ=FUNTI
      RETURN
      END



      FUNCTION CPRZ0(SS)
C     ******************
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 CPRZ0,SIGMA,CINTZZ
      COMMON / WEAK1 / ZZ,AMZ,GAMM,GM,GG,BW,ABW,ABW1,FI,DLR
      COMMON / WEAK2 / C0,C1,C2,D0,D1,D2,C3,C4,D3,D4
      COMMON / KEYGSW / KEYGSW
      common /ifzrun/ ifrunz
      S=SS
      IF(KEYGSW.EQ.0) THEN
        SIGMA=DCMPLX(0D0,AMZ*GAMM)
      ELSE
***     SIGMA=CINTZZ(0,S)
        SIGMA=DCMPLX(0D0,AMZ*GAMM)
      ENDIF
      if (ifrunz.eq.1) SIGMA=DCMPLX(0D0,s*GAMM/amz)
      CPRZ0=DCMPLX((S-AMZ**2),0D0)
      CPRZ0=1/(CPRZ0+SIGMA)
      RETURN
      END


      SUBROUTINE COUPLER(IDF,CVE,CAE,CVF,CAF)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON / WEAK2 / C0,C1,C2,D0,D1,D2,C3,C4,D3,D4
      COMMON / CONST / PI,ALFA,ALF1,QE2,QF2,QEF
      COMMON / WEAK3 /    XC0(4),XC1(4),XC2(4),XD1(4),XD2(4),
     $ XD0(4),XC3(4),XC4(4),XD3(4),XD4(4)
      COMMON /INSPIN / SEPS1,SEPS2
C     COLOUR FACTOR
      COLR=1.
      IF(IABS(IDF).GT.2) COLR=3.

      C0=0.
      C1=0.
      C2=0.
      C3=0.
      C4=0.

      D0=0.
      D1=0.
      D2=0.
      D3=0.
      D4=0.
C$ZWAS
C----- PART FOR SPINER
      DO 777 K=1,4
C     AT FIRST C0....D4 ARE DEFINED FOR DIFFERENT TAU HELICITIES
C     LATER APPOPRIATE QUANTITIES ARE CALCULATED AND STORED IN WEAK3
      IF(K.EQ.1) T1= 1.
      IF(K.EQ.1) T2= 1.
      IF(K.EQ.2) T1=-1.
      IF(K.EQ.3) T1= 1.
      IF(K.EQ.3) T2=-1.
      IF(K.EQ.4) T1=-1.
      PIN=1.+SEPS1*SEPS2
      POUT=1.+T1*T2
      SIN=SEPS1+SEPS2
      SOUT=T1+T2
      VIN =CVE*PIN +CAE*SIN
      VOUT=CVF*POUT+CAF*SOUT
      AIN =CAE*PIN +CVE*SIN
      AOUT=CAF*POUT+CVF*SOUT
      V2IN =(CVE**2+CAE**2)*PIN +2.*CVE*CAE*SIN
      V2OUT=(CVF**2+CAF**2)*POUT+2.*CVF*CAF*SOUT
      A2IN =(CVE**2+CAE**2)*SIN +2.*CVE*CAE*PIN
      A2OUT=(CVF**2+CAF**2)*SOUT+2.*CVF*CAF*POUT
      XC0(K)=COLR*QE2*QF2*PIN*POUT
      XD0(K)=2.*COLR*QE2*QF2*SIN*SOUT
      XC1(K)=2.*COLR*QEF*VIN*VOUT
      XD1(K)=4.*COLR*QEF*AIN*AOUT
      XC2(K)=COLR*V2IN*V2OUT
      XD2(K)=2.*COLR*A2IN*A2OUT
      XC3(K)=QEF*COLR*(VIN*AOUT-AIN*VOUT)
      XD3(K)=QEF*COLR*(VIN*AOUT+AIN*VOUT)
      XC4(K)=0.5*COLR*(V2IN*A2OUT-A2IN*V2OUT)
      XD4(K)=0.5*COLR*(V2IN*A2OUT+A2IN*V2OUT)
C     COUPLING CONSTANS FOR UNPOLARIZED CROSS SECTION (FINAL STATE)
C     SUMMED FROM CONTRIBUTIONS OF DIFFERENT SPIN CONFIGURATIONS
C     OF COURSE C0... CONTAIN EFFECTS OF ELECTRON BEAM POLARIZATION
      C0=C0+XC0(K)/4.
      C1=C1+XC1(K)/4.
      C2=C2+XC2(K)/4.
      C3=C3+XC3(K)/4.
      C4=C4+XC4(K)/4.
      D0=D0+XD0(K)/4.
      D1=D1+XD1(K)/4.
      D2=D2+XD2(K)/4.
      D3=D3+XD3(K)/4.
      D4=D4+XD4(K)/4.
  777 CONTINUE
      RETURN
      END
      SUBROUTINE SETCUP(IDF,SIN2,QC,CV,CA)
C     ************************************
      IMPLICIT REAL*8(A-H,O-Z)

      B=1./SQRT(16.*SIN2*(1.-SIN2))
C;;   AMZ=2.*B*74.584
C;;   GAMM=2.5

      IDFERM=IABS(IDF)
      GO TO (1,2,3,4), IDFERM
C     NOT USED  ( NEUTRINO )
    1 QC=0.
      CV=B
      CA=B
      GO TO 500
C     ELECTRON, MUON ...
    2 QC=-1.
      CV=(-1.+4.*SIN2)*B
      CA=-B
      GO TO 500
C     UP QUARKS
    3 QC=2./3.
      CV=(1.-8./3.*SIN2)*B
      CA=B
      GO TO 500
C     DOWN QUARKS
    4 QC=-1./3.
      CV=(-1.+4./3.*SIN2)*B
      CA=-B
  500 IF(IDF.GT.0) RETURN
C     ANTIPARTICLES
      QC=-QC
      CV=-CV
      CA=+CA
      RETURN
      END


CC    ##############################################################







      FUNCTION ZAMPL(E1,E2)
C
C     ****************************************************************
C     *          *************************************               *
C     *          ***           KORAL               ***               *
C     *          *************************************               *
C     *                                                              *
C     *     S. JADACH,  Z.WAS,  JAGELLONIAN UNIVERSITY, KRAKOW       *
C     *     THE MONTE CARLO PROGRAM SIMULATING   THE   PROCESS       *
C     *         E+  E-   INTO   TAU+  TAU-  ( PHOTON )               *
C     *     IN QED TO ORDER  ALPHA**3  INCLUDING SPIN EFFECTS        *
C     *     AND MASS IN THE FINAL STATE.      Z0   CONTRIBUTION      *
C     *     INCLUDED   IN   THE  HIGH   ENERGY  APPROXIMATION.       *
C     *                  VERSION  ????????  ??                       *
C     ****************************************************************
C     *  THIS IS CENTRAL MENAGEMENT ROUTINE FOR TAUPAIR PRODUCTION   *
C     *  PROCESS IT CALLS OTHER ROUTINES CALCULATES THE SPIN WEIGHT  *
C     *  WHICH IS NEXT USED TO DECIDE WHETHER EVENT IS ACCEPTED OR   *
C     *  REJECTED. THIS WEIGHING AND REJECTING PROCEDURE INTRODUCES  *
C     *  CORRELATION IN THE DECAYS OF TWO TAUS AND OTHER SPIN        *
C     *  EFFECTS. E1 AND E2 ARE POLARISATION VECTORS OF E+ AND E-    *
C     ****************************************************************
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON / UTIL1 / XK,C1,S1,C2,S2,CF,SF,CG,SG,V
      COMMON / CONTRL/ SWT(6),KSPIN
c[[[[      DIMENSION E1(3),E2(3),SA(4),SB(4),H1(4),H2(4),HI1(3),HI2(3)
      DIMENSION E1(3),E2(3),SA(4),SB(4),H1(4),H2(4),HI1(4),HI2(4)
      COMPLEX*16 T1(4,4),T2(4,4),S(4,4),E(4,4)
      COMPLEX*16 ZERO,ROT,TOT
      LOGICAL LSOF,LHAR
      DATA ICONT/0/
      ZERO=DCMPLX(0.D0,0.D0)
C-----IF ELECTRONS ARE UNPOLARIZED WE TAKE EASIER PATH
      SUM=0.
      DO 30 I=1,3
   30 SUM=SUM+E1(I)**2+E2(I)**2
      IF(SUM.LT.0.0001) GO TO 200
C-----
!      PRINT *,'Z POLARYZACJA BEAMU'
C-----POLARIZED ELECTRONS

      CALL AMPLIZ(T1,T2)

      KTO=1

      KTO=2

c[[[[[[[
cc      CALL Brem_SpinStore(1,hI1,hI2)
c]]]]]]]
      CALL STest_GetHvectors(hI1,hI2)

      DO 33 KK=1,3
      H1(KK)=HI1(KK)
33    H2(KK)=HI2(KK)

      H1(4)=1.D0
      H2(4)=1.D0
      CALL SPIN(E,H1,H2)
C-----ROTATION OF THE INITIAL POLARIZATION TO THE DYNAMIC COORDINATES
      SA(1)= CG*E1(1)+SG*E1(2)
      SA(2)=-SG*E1(1)+CG*E1(2)
      SA(3)=E1(3)
      SA(4)=1.
      SB(1)= CG*E2(1)+SG*E2(2)
      SB(2)=-SG*E2(1)+CG*E2(2)
      SB(3)=E2(3)
      SB(4)=1.
      CALL SPIN(S,SA,SB)
C---- CALCULATION OF THE SPIN WEIGHT FOR POLARIZED ELECTRONS
      LSOF= XK.EQ.0.
      LHAR= .NOT.LSOF

      ROT=ZERO
      TOT=ZERO
      DO 100 I=1,4
      DO 100 J=1,4
      IF(LSOF) TOT=TOT+DCONJG(T1(I,J))*T2(I,J)+DCONJG(T2(I,J))*T1(I,J)
      IF(LHAR) TOT=TOT+DCONJG(T1(I,J))*T1(I,J)+DCONJG(T2(I,J))*T2(I,J)
      DO 100 K=1,4
      DO 100 L=1,4
      IF(LSOF) ROT=ROT+
     +S(K,L)*(DCONJG(T1(K,I))*T2(L,J)+DCONJG(T2(K,I))*T1(L,J))*E(J,I)
      IF(LHAR) ROT=ROT+
     +S(K,L)*(DCONJG(T1(K,I))*T1(L,J)+DCONJG(T2(K,I))*T2(L,J))*E(J,I)
  100 CONTINUE


      XLSP=8.
      GOTO 400
C-----
C-----UNPOLARIZED ELECTRONS
  200 CONTINUE
!      PRINT *,'BEZ POLARYZACJI BEAMU'
      CALL AMPLIZ(T1,T2)

      KTO=1

      KTO=2

c[[[[
cc      CALL Brem_SpinStore(1,hI1,hI2)
c]]]]
      CALL STest_GetHvectors(hI1,hI2)


      DO 34 KK=1,3
      H1(KK)=HI1(KK)
 34   H2(KK)=HI2(KK)

      H1(4)=1.
      H2(4)=1.
      CALL SPIN(E,H1,H2)
C-----CALCULATE SPIN WEIGHT FOR UNPOLARIZED ELECTRONS
      LSOF= XK.EQ.0.
      LHAR= .NOT.LSOF
      ROT=ZERO
      TOT=ZERO
      DO 300 I=1,4
      DO 300 J=1,4
      IF(LSOF) TOT=TOT+DCONJG(T1(I,J))*T2(I,J)+DCONJG(T2(I,J))*T1(I,J)
      IF(LHAR) TOT=TOT+DCONJG(T1(I,J))*T1(I,J)+DCONJG(T2(I,J))*T2(I,J)
      DO 300 K=1,4
      IF(LSOF) ROT=ROT+
     +(DCONJG(T1(K,I))*T2(K,J)+DCONJG(T2(K,I))*T1(K,J))*E(J,I)
      IF(LHAR) ROT=ROT+
     +(DCONJG(T1(K,I))*T1(K,J)+DCONJG(T2(K,I))*T2(K,J))*E(J,I)
  300 CONTINUE
      XLSP=2.
  400 CONTINUE
C-----
      ICONT=ICONT+1
      XAX=REAL(ROT)

ccc      IF(ICONT.LE.2) PRINT 3300,XAX
      IF(ICONT.LE.2) WRITE(16, 3300) XAX
 3300 FORMAT(25H ########## SPIN WEIGHT2 ,E20.9)
      ZAMPL=REAL(ROT)
!      DO 160 K=1,4
!      DO 160 I=1,4
!      PRINT *, T1(I,K),',T1(I,K),',I,K
! 160  PRINT *, T2(I,K),',T2(I,K),',I,K
      RETURN
      END



      SUBROUTINE AMPLIZ(T1,T2)
C********************************************************************
C  IN THIS ROUTINE THE COMPLEX SPIN AMPLITUDES FOR THE TAU          *
C  PRODUCTION PROCESS ARE CALCULATED                                *
C********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      common /ifinifin/  ifini,iffin
      COMMON / CONST / PI,ALFA,ALF1,QE2,QF2,QEF
      COMMON / WEAK  / QCE,QCF,CVE,CVF,CAE,CAF,ZPROP
      COMMON / UTIL1 / XK,C1,S1,C2,S2,CF,SF,CG,SG,V
      COMMON / ENERG / ENE,AEL2,AMF2,AMF,ALGEL,ALGMF,BETI,BT1,ATH2
      common /cyk/ cykus
      COMPLEX*16 T1(4,4),T2(4,4),ZERO,FAC,TCI(4,4),TCF(4,4)
      COMPLEX*16 ONEC,IMAG,F2,ZH,DE,CK,CL,CN,CX,SX,dex
      ZERO=DCMPLX(0.D0,0.D0)
      IMAG=DCMPLX(0.D0,1.D0)
      ONEC=DCMPLX(1.D0,0.D0)
C
C  *******************************************
C  *                                         *
C  *       HARD PHOTON CASE                  *
C  *                                         *
C  *******************************************
      Y=SQRT(1.-XK)
      BB=.5*XK/Y
      V=SQRT(ABS(1.-AMF2/Y/Y))
      GB=(1.-.5*XK)/Y
      AM=AMF/Y
!      PRINT *, 'AMPLITUDY AMPLIZ Z RACHUNKU'

      HINI= QCE/(AEL2+S1**2)*ifini
      HFIN= QCF/(AM**2+V**2*S2**2)*iffin
C-----

      DO 140 K=1,4
      DO 140 I=1,4

      T1 (I  ,K)= ZERO
  140 T2 (I  ,K)= ZERO
C-----
      DO 150 K=1,2
      DO 150 I=1,2
      EPS=-(3.D0-2.D0*I)
      TAU=-(3.D0-2.D0*K)
      GE=1.D0

      CX=GB*CF*ONEC+BB*SF*GE*IMAG
      SX=BB*CF*ONEC+GB*SF*GE*IMAG
      CK=-(EPS*TAU+C1*C2)*CX+S1*S2-GE*(TAU*C1+EPS*C2)*SX
      CN=(C1+GE*EPS)*(C2-C2*GE*TAU)*ONEC
      CL=(C1-GE*EPS*C1)*(C2+GE*TAU)*ONEC
      T1(I,K)= IMAG/2.D0/BB*((GB-BB)*S1*HINI*DE(XK,EPS,TAU)
     &+S2*HFIN*(CF*ONEC-SF*GE*IMAG)*DE(0.D0,EPS,TAU))*CK

      T1(I+2,K)=0.5D0*Y*SQRT(AEL2)*C1*HINI*(CF*ONEC+SF*GE*IMAG)
     &         *DE(XK,-C1*EPS,TAU)*CL

       T1(I,K+2)=0.5D0*AMF/Y*C2*HFIN*DE(0.D0,EPS,-C2*TAU)*CN  
!       T1(I,K+2)=0.5D0*AMF/Y*C2*HFIN*DE(0.D0,EPS,    TAU)*CN  !  *0 


C-----
      GE=-1.D0

      CX=GB*CF*ONEC+BB*SF*GE*IMAG
      SX=BB*CF*ONEC+GB*SF*GE*IMAG
      CK=-(EPS*TAU+C1*C2)*CX+S1*S2-GE*(TAU*C1+EPS*C2)*SX
      CN=(C1+GE*EPS)*(C2-C2*GE*TAU)*ONEC
      CL=(C1-GE*EPS*C1)*(C2+GE*TAU)*ONEC
**    PRINT *, 'CL=(C1-GE*EPS*C1)*(C2+GE*TAU)'
**    PRINT *, 'TYLE CL POWINNO BYC W AMPLIZ SMIEC DOPASOW. DO KORAL-B'
**    CL=(C1+GE*EPS)*(-C2+GE*TAU)*ONEC
      T2(I,K)= IMAG/2.D0/BB*((GB-BB)*S1*HINI*DE(XK,EPS,TAU)
     &+S2*HFIN*(CF*ONEC-SF*GE*IMAG)*DE(0.D0,EPS,TAU))*CK

      T2(I+2,K)=0.5D0*Y*SQRT(AEL2)*C1*HINI*(CF*ONEC+SF*GE*IMAG)
     &         *DE(XK,-C1*EPS,TAU)*CL

      T2(I,K+2)=0.5D0*AMF/Y*C2*HFIN*DE(0.D0,EPS,-C2*TAU)*CN   
!      T2(I,K+2)=0.5D0*AMF/Y*C2*HFIN*DE(0.D0,EPS,    TAU)*CN  ! *0 



 150  CONTINUE
      cykus=c2
      CALL CHANGE(T1)
      CALL CHANGE(T2)

      RETURN
      END

      SUBROUTINE CHANGE(T)
      COMPLEX*16 T(4,4),TB(4,4),ONEC,IMAG
      IMAG=DCMPLX(0.D0,1.D0)
      ONEC=DCMPLX(1.D0,0.D0)
      DO 7 K=1,4
      TB(1,K)= (T(1,K)-T(2,K))*(-IMAG)
      TB(2,K)=-(T(1,K)+T(2,K))*(-IMAG)
      TB(3,K)= (T(3,K)-T(4,K))*(-IMAG)
   7  TB(4,K)=-(T(3,K)+T(4,K))*(-IMAG)

      DO 8 K=1,4
      T(K,1)= (TB(K,1)-TB(K,2))/2.D0
      T(K,2)=-(TB(K,1)+TB(K,2))/2.D0
      T(K,3)= (TB(K,3)-TB(K,4))/2.D0
   8  T(K,4)=-(TB(K,3)+TB(K,4))/2.D0

      RETURN
      END

      FUNCTION DE(XK,EPS,TAU)
C********************************************************************
C  IN THIS ROUTINE THE COMPLEX DOUBLE PROPAGATOR                    *
C                     IS  CALCULATED                                *
C********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON / WEAK  / QCE,QCF,CVE,CVF,CAE,CAF,ZPROP
      COMMON / ENERG / ENE,AEL2,AMF2,AMF,ALGEL,ALGMF,BETI,BT1,ATH2
      common /breminprm/  xENE,xAMZ,xGAMM,xAMFIN,xsinw2,ixzonon

      COMPLEX*16 ZERO
      COMPLEX*16 ONEC,IMAG,F2,ZH,DE,CPRZ0
      ZERO=DCMPLX(0.D0,0.D0)
      IMAG=DCMPLX(0.D0,1.D0)
      ONEC=DCMPLX(1.D0,0.D0)
C
      if (ixzonon.eq.3) onec=zero
      Y=SQRT(1.-XK)
      BB=.5*XK/Y
      GB=(1.-.5*XK)/Y
      SVAR=ENE*ENE*4.D0
      SVAR1=SVAR*(1.D0-XK)
      FAC=(CVE+EPS*CAE)*(CVF+TAU*CAF)
      DE=QCE*QCF/(1.D0-XK)*ONEC+FAC*SVAR*CPRZ0(SVAR1)

      RETURN
      END

      FUNCTION DEx(XK,EPS,TAU)
C********************************************************************
C  IN THIS ROUTINE THE COMPLEX DOUBLE PROPAGATOR                    *
C                     IS  CALCULATED                                *
C********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON / WEAK  / QCE,QCF,CVE,CVF,CAE,CAF,ZPROP
      COMMON / ENERG / ENE,AEL2,AMF2,AMF,ALGEL,ALGMF,BETI,BT1,ATH2
      common /breminprm/  xENE,xAMZ,xGAMM,xAMFIN,xsinw2,ixzonon

      COMPLEX*16 ZERO
      COMPLEX*16 ONEC,IMAG,F2,ZH,DE,CPRZ0,dex
      ZERO=DCMPLX(0.D0,0.D0)
      IMAG=DCMPLX(0.D0,1.D0)
      ONEC=DCMPLX(1.D0,0.D0)
C
      if (ixzonon.eq.3) onec=zero
      Y=SQRT(1.-XK)
      BB=.5*XK/Y
      GB=(1.-.5*XK)/Y
      SVAR=ENE*ENE*4.D0
      SVAR1=SVAR*(1.D0-XK)
      FAC=(CVE+EPS*CAE)*(CVF-TAU*CAF)
      DEx=QCE*QCF/(1.D0-XK)*ONEC+FAC*SVAR*CPRZ0(SVAR1)

      RETURN
      END
 
      FUNCTION XRAL(E1,E2)
C
C     ****************************************************************
C     *          *************************************               *
C     *          ***           KORAL               ***               *
C     *          *************************************               *
C     *                                                              *
C     *     S. JADACH,  Z.WAS,  JAGELLONIAN UNIVERSITY, KRAKOW       *
C     *     THE MONTE CARLO PROGRAM SIMULATING   THE   PROCESS       *
C     *         E+  E-   INTO   TAU+  TAU-  ( PHOTON )               *
C     *     IN QED TO ORDER  ALPHA**3  INCLUDING SPIN EFFECTS        *
C     *     AND MASS IN THE FINAL STATE.      Z0   CONTRIBUTION      *
C     *     INCLUDED   IN   THE   LOW   ENERGY  APPROXIMATION.       *
C     *                  VERSION  DECEMBER  83                       *
C     ****************************************************************
C     *  THIS IS CENTRAL MENAGEMENT ROUTINE FOR TAUPAIR PRODUCTION   *
C     *  PROCESS IT CALLS OTHER ROUTINES CALCULATES THE SPIN WEIGHT  *
C     *  WHICH IS NEXT USED TO DECIDE WHETHER EVENT IS ACCEPTED OR   *
C     *  REJECTED. THIS WEIGHING AND REJECTING PROCEDURE INTRODUCES  *
C     *  CORRELATION IN THE DECAYS OF TWO TAUS AND OTHER SPIN        *
C     *  EFFECTS. E1 AND E2 ARE POLARISATION VECTORS OF E+ AND E-    *
C     ****************************************************************
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON / UTIL1 / XK,C1,S1,C2,S2,CF,SF,CG,SG,V
      COMMON / CONTRL/ SWT(6),KSPIN
c[[[[      DIMENSION E1(3),E2(3),SA(4),SB(4),H1(4),H2(4),HI1(3),HI2(3)
      DIMENSION E1(3),E2(3),SA(4),SB(4),H1(4),H2(4),HI1(4),HI2(4)
      COMPLEX*16 T1(4,4),T2(4,4),S(4,4),E(4,4),TA(4,4),TB(4,4)
      COMPLEX*16 ZERO,ROT,TOT,IMAG
      LOGICAL LSOF,LHAR
      DATA ICONT/0/
      ZERO=DCMPLX(0.D0,0.D0)
      IMAG=DCMPLX(0.D0,1.D0)
C-----IF ELECTRONS ARE UNPOLARIZED WE TAKE EASIER PATH
      SUM=0.
      DO 30 I=1,3
   30 SUM=SUM+E1(I)**2+E2(I)**2
      IF(SUM.LT.0.0001) GO TO 200
C-----
!      PRINT *,'Z POLARYZACJA BEAMU'
C-----POLARIZED ELECTRONS

      CALL AMPLIT(T1,T2)

      KTO=1

      KTO=2

c[[[[
cc      CALL Brem_SpinStore(1,hI1,hI2)
c]]]]
      CALL STest_GetHvectors(hI1,hI2)

      DO 33 KK=1,3
      H1(KK)=HI1(KK)
 33   H2(KK)=HI2(KK)

      H1(4)=1.D0
      H2(4)=1.D0
      CALL SPIN(E,H1,H2)
C-----ROTATION OF THE INITIAL POLARIZATION TO THE DYNAMIC COORDINATES
      SA(1)= CG*E1(1)+SG*E1(2)
      SA(2)=-SG*E1(1)+CG*E1(2)
      SA(3)=E1(3)
      SA(4)=1.
      SB(1)= CG*E2(1)+SG*E2(2)
      SB(2)=-SG*E2(1)+CG*E2(2)
      SB(3)=E2(3)
      SB(4)=1.
      CALL SPIN(S,SA,SB)
C---- CALCULATION OF THE SPIN WEIGHT FOR POLARIZED ELECTRONS
      LSOF= XK.EQ.0.
      LHAR= .NOT.LSOF

      ROT=ZERO
      TOT=ZERO
      DO 100 I=1,4
      DO 100 J=1,4
      IF(LSOF) TOT=TOT+DCONJG(T1(I,J))*T2(I,J)+DCONJG(T2(I,J))*T1(I,J)
      IF(LHAR) TOT=TOT+DCONJG(T1(I,J))*T1(I,J)+DCONJG(T2(I,J))*T2(I,J)
      DO 100 K=1,4
      DO 100 L=1,4
      IF(LSOF) ROT=ROT+
     +S(K,L)*(DCONJG(T1(K,I))*T2(L,J)+DCONJG(T2(K,I))*T1(L,J))*E(J,I)
      IF(LHAR) ROT=ROT+
     +S(K,L)*(DCONJG(T1(K,I))*T1(L,J)+DCONJG(T2(K,I))*T2(L,J))*E(J,I)
  100 CONTINUE


      XLSP=8.
      GOTO 400
C-----
C-----UNPOLARIZED ELECTRONS
  200 CONTINUE
!      PRINT *,'BEZ POLARYZACJI BEAMU'
      CALL AMPLIT(T1,T2)

      KTO=1

      KTO=2
c[[[[[[
cc      CALL Brem_SpinStore(1,hI1,hI2)
c]]]]]]
      CALL STest_GetHvectors(hI1,hI2)

      DO 34 KK=1,3
      H1(KK)=HI1(KK)
 34   H2(KK)=HI2(KK)

      H1(4)=1.
      H2(4)=1.
      CALL SPIN(E,H1,H2)
C-----CALCULATE SPIN WEIGHT FOR UNPOLARIZED ELECTRONS
      LSOF= XK.EQ.0.
      LHAR= .NOT.LSOF
!      IF(LHAR) PRINT *, 'HARD WAY'
      ROT=ZERO
      TOT=ZERO
      DO 300 I=1,4
      DO 300 J=1,4
      IF(LSOF) TOT=TOT+DCONJG(T1(I,J))*T2(I,J)+DCONJG(T2(I,J))*T1(I,J)
      IF(LHAR) TOT=TOT+DCONJG(T1(I,J))*T1(I,J)+DCONJG(T2(I,J))*T2(I,J)
      DO 300 K=1,4
      IF(LSOF) ROT=ROT+
     +(DCONJG(T1(K,I))*T2(K,J)+DCONJG(T2(K,I))*T1(K,J))*E(J,I)
      IF(LHAR) ROT=ROT+
     +(DCONJG(T1(K,I))*T1(K,J)+DCONJG(T2(K,I))*T2(K,J))*E(J,I)
  300 CONTINUE
      XLSP=2.
  400 CONTINUE
C-----
      ICONT=ICONT+1
      XAX=REAL(ROT)
ccc      IF(ICONT.LE.2) PRINT 3300,XAX
      IF(ICONT.LE.2) WRITE(16, 3300) XAX
 3300 FORMAT(25H $$$$$$$$$$ SPIN WEIGHT  ,F20.9)
      XRAL=REAL(ROT)
      DO 160 K=1,4
      DO 160 I=1,4
      TA(I,K)=T2(I,K)+IMAG*T1(I,K)
      TB(I,K)=T2(I,K)-IMAG*T1(I,K)

!      PRINT *, TA(I,K),',T1(I,K),',I,K
!      PRINT *, TB(I,K),',T2(I,K),',I,K
  160  continue

      RETURN
      END


      SUBROUTINE AMPLIT(T1,T2)
C********************************************************************
C  IN THIS ROUTINE THE COMPLEX SPIN AMPLITUDES FOR THE TAU          *
C  PRODUCTION PROCESS ARE CALCULATED                                *
C********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      common /ifinifin/  ifini,iffin
      COMMON / CONST / PI,ALFA,ALF1,QE2,QF2,QEF
      COMMON / WEAK  / QCE,QCF,CVE,CVF,CAE,CAF,ZPROP
      COMMON / UTIL1 / XK,C1,S1,C2,S2,CF,SF,CG,SG,V
      COMMON / BOXY  / Z1,Z2,Z3
      COMMON / ENERG / ENE,AEL2,AMF2,AMF,ALGEL,ALGMF,BETI,BT1,ATH2
      COMPLEX*16 T1(4,4),T2(4,4),ZERO,FAC,TCI(4,4),TCF(4,4)
      COMPLEX*16 ONEC,IMAG,F2,ZH,Z1,Z2,Z3
      ZERO=DCMPLX(0.D0,0.D0)
      IMAG=DCMPLX(0.D0,1.D0)
      ONEC=DCMPLX(1.D0,0.D0)
C
C  *******************************************
C  *                                         *
C  *       HARD PHOTON CASE                  *
C  *                                         *
C  *******************************************
!      PRINT *, 'AMPLITUDY AMPLIT POPRAWIONE WZGLEDEM KORAL-B '
      Y=SQRT(1.-XK)
      BB=.5*XK/Y
      V=SQRT(ABS(1.-AMF2/Y/Y))
      GB=(1.-.5*XK)/Y
      AM=AMF/Y
      HINI= QCE/Y/BB/(AEL2+S1**2)*ifini
      HFIN= QCF/BB/(AM**2+V**2*S2**2)*iffin
C-----

      TCI(1,1)= (-BB*SF                              )*HINI*ONEC
      TCI(1,2)= ( BB*C2*CF                           )*HINI*IMAG
      TCI(1,3)=                                             ZERO
      TCI(1,4)= ( AM*BB*S2*CF                        )*HINI*ONEC
      TCI(2,1)= ( BB*C1*CF                           )*HINI*ONEC
      TCI(2,2)= ( BB*C1*C2*SF                        )*HINI*IMAG
      TCI(2,3)=                                             ZERO
      TCI(2,4)= ( AM*BB*C1*S2*SF                     )*HINI*ONEC
      TCI(3,1)= (-GB*CF                              )*HINI*ONEC
      TCI(3,2)= (-GB*C2*SF                           )*HINI*IMAG
      TCI(3,3)=                                             ZERO
      TCI(3,4)= (-AM*GB*S2*SF                        )*HINI*ONEC
      TCI(4,1)= (-GB*C1*SF                           )*HINI*ONEC
      TCI(4,2)= ( GB*C1*C2*CF-S1*S2                  )*HINI*IMAG
      TCI(4,3)=                                             ZERO
      TCI(4,4)= ( AM*(GB*C1*S2*CF+S1*C2)             )*HINI*ONEC
C-----
      TCF(1,1)= (-BB*V*S2*SF                         )*HFIN*ONEC
      TCF(1,2)= ( BB*V*S2*C2*CF                      )*HFIN*IMAG
      TCF(1,3)= ( AM*BB*SF                           )*HFIN*IMAG
      TCF(1,4)= (-AM*BB*V*C2*C2*CF                   )*HFIN*ONEC
      TCF(2,1)= ( BB*V*S2*C1*CF                      )*HFIN*ONEC
      TCF(2,2)= ( BB*V*S2*C2*C1*SF                   )*HFIN*IMAG
      TCF(2,3)= (-AM*BB*C1*CF                        )*HFIN*IMAG
      TCF(2,4)= (-AM*BB*V*C2*C2*C1*SF                )*HFIN*ONEC
      TCF(3,1)= (-GB*V*S2*CF                         )*HFIN*ONEC
      TCF(3,2)= (-GB*V*S2*C2*SF                      )*HFIN*IMAG
      TCF(3,3)= ( AM*BB*CF                           )*HFIN*IMAG
      TCF(3,4)= (-AM*V*(GB*S2*S2-BB)*SF              )*HFIN*ONEC
      TCF(4,1)= (-GB*V*S2*C1*SF                      )*HFIN*ONEC
      TCF(4,2)= ( V*S2*(GB*C1*C2*CF-S1*S2)           )*HFIN*IMAG
      TCF(4,3)= ( AM*BB*C1*SF                        )*HFIN*IMAG
      TCF(4,4)= ( AM*V*((GB*S2*S2-BB)*C1*CF+S2*C2*S1))*HFIN*ONEC
      DO 140 K=1,4
      DO 140 I=1,2
      ZH=TCF(I,K)
      TCF(I  ,K)=    ZH*CF-TCF(I+2,K)*SF
  140 TCF(I+2,K)=    ZH*SF+TCF(I+2,K)*CF
C-----
*$ POPRAWKA
      SGX=C1*SQRT(AEL2)
*$ POPRAWKA
      SGX1= SGX*BB/GB
      DO 150 K=1,4
      T1(1,K)= S1*TCI(1,K)+TCF(1,K)
      T1(2,K)=(S1*TCI(2,K)+TCF(2,K))*IMAG
*$ POPRAWKA
      T1(3,K)=-SGX*TCI(1,K)*IMAG*C1
      T1(4,K)=-SGX*TCI(2,K)
C-----
      T2(1,K)= S1*TCI(3,K)+TCF(3,K)
      T2(2,K)=(S1*TCI(4,K)+TCF(4,K))*IMAG
*$ POPRAWKA
      T2(3,K)=-SGX1*TCI(3,K)*IMAG*C1
 150  T2(4,K)=-SGX1*TCI(4,K)
      RETURN
      END
      SUBROUTINE SPIN(S,E,F)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 S(4,4)
      DIMENSION E(4),F(4)
      a=1d0
      S(1,1)=DCMPLX(1.D0-E(1)*F(1)+E(2)*F(2)+E(3)*F(3),  0.D0)
      S(2,2)=DCMPLX(1.D0+E(1)*F(1)-E(2)*F(2)+E(3)*F(3),  0.D0)
      S(3,3)=DCMPLX(1.D0-E(1)*F(1)-E(2)*F(2)-E(3)*F(3),  0.D0) *a
      S(4,4)=DCMPLX(1.D0+E(1)*F(1)+E(2)*F(2)-E(3)*F(3),  0.D0) *a
      S(1,2)=DCMPLX(   F(3)+E(3)        ,+E(1)*F(2)+E(2)*F(1))
      S(1,3)=DCMPLX(   F(1)-E(1)        ,-E(2)*F(3)+E(3)*F(2)) *a
      S(1,4)=DCMPLX( E(1)*F(3)+E(3)*F(1),    E(2)  +  F(2)   ) *a
      S(2,3)=DCMPLX(-E(1)*F(3)+E(3)*F(1),   -E(2)  +  F(2)   ) *a
      S(2,4)=DCMPLX(   F(1)+E(1)        , E(2)*F(3)+E(3)*F(2)) *a
      S(3,4)=DCMPLX(  -F(3)+E(3)        ,-E(1)*F(2)+E(2)*F(1)) *a
      S(4,1)=DCONJG(S(1,4))
      S(4,2)=DCONJG(S(2,4))
      S(4,3)=DCONJG(S(3,4))
      S(3,1)=DCONJG(S(1,3))
      S(3,2)=DCONJG(S(2,3))
      S(2,1)=DCONJG(S(1,2))
      
      RETURN
      END
      SUBROUTINE Brem_Tralor(KTO,VEC)
      IMPLICIT REAL*8(A-H,O-Z)
C
C     ******************************************************************
C     *   TRALOR TRANSFORMES FOUR-VECTOR VEC FROM TAU REST SYSTEM      *
C     *   TO LAB SYSTEM. RECOMMENDED TO USE  FOR DECAY PRODUCTS.       *
C     *   KTO=1,2 DENOTES TAU+ AND TAU- CORRESPONDINGLY                *
C     ******************************************************************
C
      COMMON / UTIL1 / XK,C1,S1,C2,S2,CF,SF,CG,SG,V

      COMMON / ENERG / ENE,AEL2,AMF2,AMF,ALGEL,ALGMF,BETI,BT1,ATH2
      DIMENSION VEC(4),VEC1(4),TL(4,4)
C
      Y=SQRT(1.-XK)
      A=Y/AMF
      B=SQRT((A-1.)*(A+1.))
      IF(KTO.EQ.2) B=-B
      BB=.5*XK/Y
      GB=(1.-.5*XK)/Y
C
      H1= C1*S2*CF+GB*S1*C2
      H2=-S1*S2*CF+GB*C1*C2
      TL(1,1)= CF
      TL(1,2)=-C2*SF
      TL(1,3)=-A*S2*SF
      TL(1,4)=-B*S2*SF
      TL(2,1)= C1*SF
      TL(2,2)= C1*C2*CF-GB*S1*S2
      TL(2,3)=-B*BB*S1+A*H1
      TL(2,4)=-A*BB*S1+B*H1
      TL(3,1)=-S1*SF
      TL(3,2)=-S1*C2*CF-GB*C1*S2
      TL(3,3)=-B*BB*C1+A*H2
      TL(3,4)=-A*BB*C1+B*H2
      TL(4,1)= 0.
      TL(4,2)= BB*S2
      TL(4,3)= B*GB-A*BB*C2
      TL(4,4)= A*GB-B*BB*C2
C
      DO 100 I=1,4
  100 VEC1(I)=VEC(I)
      DO 110 I=1,4
      SUM=0.
      DO 105 J=1,4
  105 SUM=SUM+TL(I,J)*VEC1(J)
  110 VEC(I)=SUM
C     ROTATION ARROUND THE BEAM
      A=VEC(1)
      VEC(1)=CG*A-SG*VEC(2)
      VEC(2)=SG*A+CG*VEC(2)
C
      RETURN
      END
      FUNCTION WINTH(DUM)
      IMPLICIT REAL*8(A-H,O-Z)
C     ****************************************************************
C     * HERE THE CONTRIBUTION FROM THE INTERFERENCE OF THE INITIAL   *
C     * AND FINAL STATE BREMSSTRAHLUNG IS CALCULATED FOR USE IN      *
C     * SUBROUTINE EVENT  ( HARD BREMSSTARHLUNG ONLY)                *
C     ****************************************************************
C
      common /ifinifin/  ifini,iffin
      COMMON / CONST / PI,ALFA,ALF1,QE2,QF2,QEF
      COMMON / UTIL1 / XK,CC1,SS1,CC2,SS2,CF,SF,CG,SG,V
      COMMON / UTIL2 / XK0,XKMIN,XKMAX
      COMMON / ENERG / ENE,AEL2,AMF2,AMF,ALGEL,ALGMF,BETI,BT1,ATH2
      EQUIVALENCE (AMF2,AM),(AEL2,AL)
C
      XINI=SQRT(1.D0-AEL2)
      V=SQRT(ABS(1.-AMF2/(1.D0-XK)))
      C1=CC1*XINI
      S1=SS1*XINI
      C2=CC2*V
      S2=SS2*V
      Z1=SQRT(1.-XK)*S1*S2*CF-(1.-.5*XK)*C1*C2
      Z2=1.-.5*XK
      T = Z1+Z2+.5*XK*( C1-C2)
      U1=-Z1+Z2+.5*XK*(-C1-C2)
      U =-Z1+Z2+.5*XK*( C1+C2)
      T1= Z1+Z2+.5*XK*(-C1+C2)
      TT=T*T
      UU=U*U
      TT1=T1*T1
      UU1=U1*U1
      SS=2.
      SP=2.*(1.-XK)
      X1=XK*(1.-C1)
      X2=XK*(1.+C1)
      IF(C1.GT.0.9) X1=XK*(AEL2+S1**2)/(1.+C1)
      IF(C1.LT.-.9) X2=XK*(AEL2+S1**2)/(1.-C1)
      Y1=XK*(1.-C2)
      Y2=XK*(1.+C2)
      AINI=(TT+UU+AM/SP*(T+U)**2)*(1.-AL/SP*X1/X2)
     $   +(TT1+UU1+AM/SP*(T1+U1)**2)*(1.-AL/SP*X2/X1)
      AINI=AINI*QE2/(SP*X1*X2)* ifini
      AFIN=(TT+UU1+AM*SS)*(1.-AM*(Y1+Y2)/Y2/SS)
     $    +(UU+TT1+AM*SS)*(1.-AM*(Y1+Y2)/Y1/SS)
     $    +AM*(X1*X1+X2*X2)/SS-4.*AM*(SS-SP)
      AFIN=AFIN*QF2/(SS*Y1*Y2)*iffin
      ANTR=(TT+TT1+UU+UU1+AM*SS+AM*SP)
     $     *(T*X2*Y2+T1*X1*Y1-U*X2*Y1-U1*X1*Y2)
     $      +AM*X1*X2*((SS-SP)*(T+T1-U-U1)-(X1-X2)*(Y1-Y2))
      ANTR=ANTR*QEF/(SS*SP*X1*X2*Y1*Y2)* ifini*iffin
      WINTH=ANTR+(AINI+AFIN)
      RETURN
      END

      FUNCTION RNDM(DUMM)
      IMPLICIT REAL*8 (A-H,O-Z)
      RNDM=0.5D0
      RETURN
      END



      SUBROUTINE Brem_SpinStore(mode,h1,h2)
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//   stores, mode=0 or gives out mode=1 polarimetric 3-vectors of tau decay        //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INTEGER mode,k
      REAL*8 h1(3),h2(3),Sh1(3),Sh2(3)
      SAVE Sh1,Sh2
      Do k=1,3
        If (mode.eq.0) then
          sh1(k)=h1(k)
          sh2(k)=h2(k)
        else
          h1(k)=sh1(k)
          h2(k)=sh2(k)
        endif
      enddo
      end

      SUBROUTINE Brem_SpinBeam(mode,e1,e2)
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//   stores, mode=0 or gives out mode=1 polarzation 3-vectors of the beams         //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INTEGER mode,k
      REAL*8 e1(3),e2(3),Se1(3),Se2(3)
      SAVE Se1,Se2
      Do k=1,3
        If (mode.eq.0) then
          se1(k)=e1(k)
          se2(k)=e2(k)
        else
          e1(k)=se1(k)
          e2(k)=se2(k)
        endif
      enddo
      end
