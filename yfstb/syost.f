* % This file can be printed with TeX. The following line sets the format.
* \def\alphat{\tilde\alpha} \def\betat{\tilde\beta} \input fortex


***********************************************************************
***********************************************************************
*************        YFS3 approximate distr    ************************
***********************************************************************
***********************************************************************

 
      SUBROUTINE sfdist2(QQ,P1,P2,Q1,Q2,PH1,PH2,DIST2)
C     ***********************************************
C Provides double bremsstrahlung distribution - FINAL state brem.
C INPUT:  P1,P2,Q1,Q2,PH1,PH2, four momenta
C OUTPUT: DIST2     double bremsstrahlung distribution
C     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DCSkey, LLkey, Zkey
      DIMENSION QQ(*),P1(*),P2(*),Q1(*),Q2(*),PH1(*),PH2(*)
      DIMENSION PR1(4),PR2(4),PH1R(4),PH2R(4),QR1(4),QR2(4)
      COMMON /trmkey/ DCSkey, LLkey, Zkey
      SAVE /trmkey/
 
      CALL sreduz2(QQ,Q1,Q2,PH1,PH2,QR1,QR2,PH1R,PH2R)
      CALL sreduz0(QQ,P1,P2,PR1,PR2)
      SVAR1 = QQ(4)**2-QQ(3)**2-QQ(2)**2-QQ(1)**2
C infrared factors from reduced momenta
C double bremsstrahlung Xsect in next-to-leading log approx.
      CALL sgsfin2(Q1,Q2,PH1,PH2,GF1,GF2)
      CALL sgthet1(QR1,QR2,PR1,COSTH1,COSTH2)
      ANDI11= bornsv(SVAR1,COSTH1)
      ANDI12= bornsv(SVAR1,COSTH2)
      DIST2 =   GF1*ANDI11+   GF2*ANDI12
      END
 
 
      SUBROUTINE sgthet1(P1,P2,Q1,COSTH1,COSTH2)
C     *****************************************
C Calculates CosTh1 and CosTh2 between BEAM amd FINAL
C fermion momenta in final fermion rest frame Q1(4)+Q2(4)=0
C     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION P1(*),P2(*),Q1(*)
      COSTH1 = (P1(1)*Q1(1) +P1(2)*Q1(2) +P1(3)*Q1(3))
     $    /SQRT((Q1(1)**2 +Q1(2)**2 +Q1(3)**2)
     $ *(P1(1)**2 +P1(2)**2 +P1(3)**2))
      COSTH2 =-(P2(1)*Q1(1) +P2(2)*Q1(2) +P2(3)*Q1(3))
     $    /SQRT((Q1(1)**2 +Q1(2)**2 +Q1(3)**2)
     $ *(P2(1)**2 +P2(2)**2 +P2(3)**2))
      END
      SUBROUTINE sreduz0(QQ,P1,P2,PR1,PR2)
C     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION QQ(4),P1(4),P2(4),PR1(4),PR2(4)
      DO 20 K=1,4
      PR1(K)=P1(K)
 20   PR2(K)=P2(K)
      END
      SUBROUTINE sreduz2(QQ,P1,P2,PH1,PH2,PR1,PR2,PH1R,PH2R)
C     *****************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION QQ(*), P1(*),  P2(*),  PH1(*),  PH2(*)
      DIMENSION        PR1(*), PR2(*), PH1R(*), PH2R(*)
      DO 20 K=1,4
      PH1R(K)=PH1(K)
      PH2R(K)=PH2(K)
      PR1(K)=P1(K)
 20   PR2(K)=P2(K)
      END
 
      SUBROUTINE sgsfin2(P1,P2,PH1,PH2,F1,F2)
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
c       call dumpt(6,"===p1===",p1)
c       call dumpt(6,"   p2   ",p2)
c       call dumpt(6,"--ph1---",ph1)
c       call dumpt(6,"  ph2   ",ph2)
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
C.. this correction reconstructs double collinear limit with 50%
C.. and affects below photon-fermion angle  <0.1 ammi/ene
C.. without correction error in this limit 1000%
c      SFAC1  =  2D0/(PP*AA1*BB1)
c      SFAC2  =  2D0/(PP*AA2*BB2)
c      DELT=(AM2/(2D0*PP))**2*(B2**2*A1**2+A2**2*B1**2)*
c     #  ( B1*B2/(A1*A2)/(A1+A2)**2
c     #   +A1*A2/(B1*B2)/(B1+B2)**2  )
c      WMINF=2D0*DELT/(X1+X2)
c      WMM=WWM(A1,B1)*WWM(A2,B2)+WMINF
c      F1 = X1*SFAC1*SFAC2/8D0*WMM
c      F2 = X2*SFAC1*SFAC2/8D0*WMM
C...end of correction............................................
      END
 
      FUNCTION bornsv(SVARI,COSTHE)
C     ***********************************
C THIS ROUTINE PROVIDES BORN DIFFERENTIAL CROSS SECTION
C a version without COMPLEX*16
C     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DCSkey, LLkey, Zkey
C     COMMON / WEKING / ENE,AMAZ,GAMMZ,AMEL,AMFIN,XK0,SINW2,IDE,IDF
      COMMON / WEKINP / AMAZ,GAMMZ,SINW2,IDE
      COMMON / BHPAR1 / CMS,AMFIN
      COMMON / BHPAR2 / CMSENE,AMEL
      COMMON /COEFF/ VE,AE
      COMMON /trmkey/ DCSkey, LLkey, Zkey
      SAVE / WEKING /,/ BHPAR1 /,/ BHPAR2 /,/ WEKINP /,/COEFF/,
     $     /trmkey/
 
C@@   Z switch added by S.Y.
![[[[[ gamma only
      IF (Zkey .EQ. 0) THEN
          BORN    = 1D0+ COSTHE**2
      ELSE
![[[[[ include Z exchange
          SINW2=0.2259D0
          AMAZ=91.161D0
          GAMMZ=2.534D0
          AMEL = 0.0D0
          AMFIN = 0.0D0
C@@
C         VF=VE
C         AF=AE
          QE= -1D0
          QF= -1D0
          AA= 4D0*SQRT(SINW2*(1D0-SINW2))
          VE= (-1D0+4*SINW2)/AA
          AE= 1D0/AA
          VF= (-1D0+4*SINW2)/AA
          AF= 1D0/AA
          S = SVARI
          CHI2 = S**2/((S-AMAZ**2)**2+(GAMMZ*AMAZ)**2)
          RECHI=(S-AMAZ**2)*S/((S-AMAZ**2)**2+(GAMMZ*AMAZ)**2)
          XE= VE**2 +AE**2
          XF= VF**2 +AF**2
          YE= 2*VE*AE
          YF= 2*VF*AF
          FF0= QE**2*QF**2 +2*RECHI*QE*QF*VE*VF +CHI2*XE*XF
          FF1=     +2*RECHI*QE*QF*AE*AF +CHI2*YE*YF
          BORN    = (1D0+ COSTHE**2)*FF0 +2D0*COSTHE*FF1
c         BORN    = (1D0+ COSTHE**2
c    1       +4D0*(AMEL**2+AMFIN**2)/S*(1D0-COSTHE**2)
c    1       +16*AMEL**2*AMFIN**2/S**2*COSTHE**2)*FF0
c    1       +2D0*COSTHE*FF1
      ENDIF
C************
C THIS IS A BIT CRUDE METHOD OF INTRODUCING THRESHOLD BEHAVIOUR
      IF(    SVARI.LE. 4D0*AMFIN**2) THEN
        THRESH=0D0
      ELSEIF(SVARI.LE.16D0*AMFIN**2) THEN
        AMX2=4D0*AMFIN**2/SVARI
        THRESH=SQRT(1D0-AMX2)*(1D0+AMX2/2D0)
      ELSE
        THRESH=1D0
      ENDIF
      bornsv= BORN*THRESH
      END !\end




 
      !--------------------------------------------------------------
      !     BHABHA SCATTERING WITH TWO PHOTON EMISSION
      !         For interface with external program
      !     Scott Yost                              January 26, 1993
      !
      ! This program calculates the cross-section for Bhabha scattering
      ! with two photon bremsstrahlung, including Z boson exchange,
      !           $ e^{+} e^{-} \rightarrow e^{+} e^{-} + 2\gamma$,
      ! both exactly and in the leading log approximation. All calculations
      ! assume the massless limit for the electrons. However, electron
      ! masses are incorporated in the spinor and vector products, as well as
      ! in the kinematic computations. This allows a check that the mass
      ! effects are in fact not too large. The process of interest
      ! is an order $\alpha^2$ amplitude. Comments use \TeX notation.
      !-------------------------------------------------------------
 
 
******************* MASTER INTERFACE SUBROUTINE *********************
 
       SUBROUTINE TWOPHO(DCS, DCSLL, DCSS, P1, P2, Q1, Q2, K1, K2)
      !--------------------------------------------------------------
      ! Calculates the exact and leading log matrix elements
      ! for given momenta, using BHLUMI notation:
      ! P1 = incoming positron, P2 = outgoing positron
      ! Q1 = incoming electron, Q2 = outgoing electron
      ! K1, K2 = photons. (The subroutine is symmetric in these.)
      ! The components are $(x, y, z, t)$, with $z$ = beam direction.
      ! The results are DCS   = exact cross-section,
      !                 DCSLL = leading log approximation,
      !                 DCSS  = soft approximation.
      ! ``Term keys'', which determine which terms are included in the
      ! exact and leading log calculations, are passed through a
      ! common block `trmkey'.
      !    DCSkey and LLkey are 7-digit strings of 0's and 1's, which
      ! determine, in order from most-significant to least significant
      ! digit, the presence or absence of the following terms:
      ! digit     DCSkey              LLkey
      !   7   explicit mass terms  explicit mass terms
      !   6   positron line        positron line      ($t$ channel)
      !   5   electron line        electron line      ($t$ channel)
      !   4   mixed $t$ channel    mixed term $m$
      !   3   final state          final state        ($s$ channel)
      !   2   initial state        initial state      ($s$ channel)
      !   1   mixed $s$ channel    mixed term $m'$
      ! where leading log terms $m$ and $m'$ can be found in the paper.
      ! The correspondance between similarly-named terms in DCSkey and
      ! LLkey is not precise. The explicit mass terms are the $W$ factors
      ! taken from the LL calculation, in either case. Exact mass
      ! corrections are not included in the exact calculation.
      ! For example, to include all terms, but mass terms only in the
      ! LL case, use DCSkey = 111111 and LLkey = 1111111.
      !    Zkey is simply 0 or 1, and determines whether $Z^0$
      ! exchange terms are to be included.
      !
      ! A data file data.-1 for use with final.f is created.
      !-------------------------------------------------------------
       IMPLICIT NONE
       INTEGER DCSkey, LLkey, Zkey, data
       PARAMETER (data=2)
       COMMON /trmkey/ DCSkey, LLkey, Zkey
       REAL*8 ampsqr, ldlog, soft, DCS, DCSLL, DCSS, theta, phi,
     &        thetaP, thetaE, theta1, theta2, phi1, phi2, phiP,
     &        E1rat, E2rat, degree
      ! Four-momenta of particles, in BHLUMI notation:
       REAL*8 P1(4), P2(4), Q1(4), Q2(4), K1(4), K2(4)
      ! Spinor representations of the momenta, in my notation:
       COMPLEX*16     SP1(4), SP2(4), SQ1(4), SQ2(4), SK1(4), SK2(4)
       COMMON /spnrs/ SP1(4), SP2(4), SQ1(4), SQ2(4), SK1(4), SK2(4)
       EXTERNAL spnvec, ampsqr, ldlog, soft
       SAVE /trmkey/, /spinor/
 
      ! one degree, in radians:
       degree = 90.0D0/DATAN(1.0D0)
 
      ! translation of vector variables into spinors, in my notation:
       CALL spnvec('P1', P1, SP1, theta, phi)
       CALL spnvec('Q1', Q1, SQ1, theta, phi)
       CALL spnvec('P2', P2, SP2, theta, phi)
           thetaP = theta*degree
           phiP = phi*degree
       CALL spnvec('Q2', Q2, SQ2, theta, phi)
           thetaE = theta*degree
       CALL spnvec('K1', K1, SK1, theta, phi)
           E1rat = K1(4)/P1(4)
           theta1 = theta*degree
           phi1 = phi*degree
       CALL spnvec('K2', K2, SK2, theta, phi)
           E2rat = K2(4)/P1(4)
           theta2 = theta*degree
           phi2 = phi*degree
 
      ! Make data file for baba.f:
       OPEN(unit=data, file='data.-1', status='unknown')
       WRITE(data,*) '# term keys:'
       WRITE(data,*) DCSkey, LLkey, Zkey
       WRITE(data,*) '# E1rat = ', E1rat, ', E2rat value:'
       WRITE(data,*)  E2rat
       WRITE(data,*) '# theta1, theta2, thetaP, thetaE:'
       WRITE(data,100) theta1, theta2, thetaP, thetaE
       WRITE(data,*) '# phi1, phi2, phiP:'
       WRITE(data,200) phi1, phi2, phiP
       CLOSE(data)
 
       DCS   = ampsqr(SQ1, SP1, SQ2, SP2, SK1, SK2)
       DCSLL =  ldlog(SQ1, SP1, SQ2, SP2, SK1, SK2)
       DCSS  =   soft(SQ1, SP1, SQ2, SP2, SK1, SK2)
       RETURN
 
100    FORMAT(G20.13, 1X, G20.13, 1X, F17.13, 1X, F17.13)
200    FORMAT(G20.13, 1X, G20.13, 1X, G20.13)
       END ! SUBROUTINE TWOPHO
 
 
******* EVALUATION OF AVERAGED/SUMMED SQUARED MATRIX ELEMENT ********
 
 
       FUNCTION ampsqr(P, Pp, Q, Qp, K1, K2)
      !-----------------------------------------------------------------
      ! The exact calculation of the averaged/summed squared matrix element.
      ! An overall factor of $e^8$ is omitted. The helicities $\lambda$ and
      ! $\rho_i$ of the incoming electron and photons are summed in loops.
      ! The helicities $\lambda', \mu, \mu'$ of the incoming positron and
      ! outgoing electron and positron, respectively, are summed explicitly
      ! or fixed by helicity conservation. All explicit amplitudes are $t$
      ! channel expressions from my notes, and the $s$ channel is obtained
      ! by crossing $(q,\mu) \leftrightarrow\ (-p',-\lambda')$.
      !
      ! Depending on the electron and positron helicities, the total
      ! amplitude can be pure $t$ channel, pure $s$ channel, or mixed:
      !          ampT  = pure $t$ channel amplitude
      !          ampS  = pure $s$ channel amplitude
      !          ampST = mixed $s+t$ channel amplitude
      !
      ! The helicity cases are:
      !     Pure $t$ channel when $\lambda = \mu = \lambda' = \mu'$.
      !     Pure $s$ channel when $\lambda = -\mu =-\lambda' = \mu'$.
      !     Both channels appear if $\lambda = \mu = -\lambda' =-\mu'$.
      !     Both channels vanish if $\lambda = -\mu = \lambda' = -\mu'$.
      !-----------------------------------------------------------------
       IMPLICIT NONE
       INTEGER lam, rho1, rho2      ! helicities $\lambda$, $\rho_i$
       INTEGER DCSkey, LLkey, Zkey
       COMMON /trmkey/ DCSkey, LLkey, Zkey
       REAL*8 ampsqr, absqr
       COMPLEX*16 P(4), Pp(4), Q(4), Qp(4), K1(4), K2(4),
     &            PpX(4), QX(4), spr, ampS, ampT, ampST, amp
       LOGICAL getkey
       EXTERNAL spr, cross, amp, absqr, getkey
       SAVE /trmkey/
 
       ampsqr = 0.0D0
       ampT  = DCMPLX(0.0D0)
       ampS  = DCMPLX(0.0D0)
       ampST = DCMPLX(0.0D0)
       CALL cross(QX, PpX, Q, Pp)
 
       DO  lam = -1, 1, 2          ! begin loop over $\lambda$
       DO rho2 = -1, 1, 2          ! begin loop over $\rho_2$
       DO rho1 = -1, 1, 2          ! begin loop over $\rho_1$
         ! $t$ channel
          IF (getkey(DCSkey, 6) .OR. getkey(DCSkey, 5)
     &        .OR. getkey(DCSkey, 4)) THEN
           ampT  = amp(P, Pp, Q, Qp, K1, K2, lam,  lam, rho1, rho2,'t')
           ampST = amp(P, Pp, Q, Qp, K1, K2, lam, -lam, rho1, rho2,'t')
          ELSE
           ampST = DCMPLX(0.0D0)
          END IF
 
         ! $s$ channel (by crossing)
          IF (getkey(DCSkey, 3) .OR. getkey(DCSkey, 2)
     &        .OR. getkey(DCSkey, 1)) THEN
           ampS  = amp(P, PpX, QX, Qp, K1, K2, lam, lam, rho1, rho2,'s')
           ampST = ampST
     &           + amp(P, PpX, QX, Qp, K1, K2, lam,-lam, rho1, rho2,'s')
          END IF
 
          ampsqr = ampsqr + absqr(ampS) + absqr(ampT) + absqr(ampST)
       END DO  ! end loop over $\rho_1$
       END DO  ! end loop over $\rho_2$
       END DO  ! end loop over $\lambda$
 
      ! Average over four incoming helicity assignments.
       ampsqr = ampsqr/4.0D0
       RETURN
       END ! FUNCTION ampsqr
 
 
       FUNCTION amp(P, Pp, Q, Qp, K1, K2, lam, lamp, rho1, rho2, chan)
      !----------------------------------------------------------------
      ! Evaluation of the total amplitude for the $T$ channel. $S$ channel
      ! amplitudes may be obtained by crossing $p'$ and $q$ and replacing
      ! `lamp' = $\lambda'$ by $-\mu$. A `double crossing' transformation
      ! $(p, q, \lambda) \leftrightarrow (-q', -p', -\lambda')$ simplifies the
      ! calculation by relating electron line and positron line radiation.
      ! A common overall phase of $i$ is omitted from the amplitude.
      !----------------------------------------------------------------
       IMPLICIT NONE
       INTEGER lam, lamp, rho1, rho2    ! helicities
       INTEGER DCSkey, LLkey, Zkey, shift
       COMMON /trmkey/ DCSkey, LLkey, Zkey
       COMPLEX*16 amp, P(4), Pp(4), Q(4), Qp(4), K1(4),
     &            K2(4), PX(4), PpX(4), QX(4), QpX(4),
     &            poslin, elelin, diflin, same, diff
       LOGICAL getkey
       CHARACTER*1 chan                 ! channel, 's' or 't'
       EXTERNAL cross, same, diff, getkey
       SAVE /trmkey/
 
       elelin = DCMPLX(0.0D0)
       poslin = DCMPLX(0.0D0)
       diflin = DCMPLX(0.0D0)
       CALL cross(QX, PpX, Q, Pp)
       CALL cross(PX, QpX, P, Qp)
 
      ! The terms represented by positions 1 -- 3 in the $S$ channel
      ! are related by crossing to corresponding $T$ channel terms
      ! in positions 4 -- 6 of DCSkey. See subroutine TWOPHO.
       IF (chan .EQ. 's') THEN
          shift = 0
       ELSE
          shift = 3
       END IF
 
      ! electron line emission:
       IF (getkey(DCSkey, 2 + shift)) elelin =
     &    same( P,  Pp,  Q,  Qp, K1, K2,  lam,  lamp, rho1, rho2, chan)
      ! positron line emission:   (uses double crossing)
       IF (getkey(DCSkey, 3 + shift)) poslin =
     &    same(PX, PpX, QX, QpX, K1, K2, -lamp, -lam, rho1, rho2, chan)
      ! emission from different lines:  (two photon orderings)
       IF (getkey(DCSkey, 1 + shift)) diflin =
     &    diff(P, Pp, Q, Qp, K1, K2, lam, lamp, rho1, rho2, chan)
     &  + diff(P, Pp, Q, Qp, K2, K1, lam, lamp, rho2, rho1, chan)
 
       amp = elelin + poslin + diflin
       RETURN
       END ! FUNCTION amp
 
 
       FUNCTION same(P, Pp, Q, Qp, K1, K2, lam, lamp, rho1, rho2, chan)
      !----------------------------------------------------------
      ! Evaluation of the $t$ channel amplitude with both photons
      ! emitted from the same line. The electron line expressions
      ! from my notes are used here. The positron line case may be
      ! obtained from the same function using double crossing.
      ! `FUNCTION llcor' is a leading log mass correction factor.
      !----------------------------------------------------------
       IMPLICIT NONE
       INTEGER lam, lamp, rho1, rho2     ! helicities
       REAL*8 z, Del, Delhat, vpr, llcor
      ! z = argument of propagator ($t'$ in my notes). Del, Delhat are my
      !       $\Delta = (p - k_1 - k_2)^2, \hat\Delta = (q + k_1 + k_2)^2.$
       COMPLEX*16 same, P(4), Pp(4), Q(4), Qp(4), K1(4), K2(4),
     &            Ki(4), Kj(4), H(4), Hp(4), Hhat(4), Hphat(4),
     &            term1, term2, term3, denom1, denom2,
     &            prop, spr, braket
      ! H, Hp, Hhat, Hphat correspond to $h_i, h'_i, \hat{h}_i, \hat{h}'_i$
      !      in my notes, and Ki, Kj correspond to my $k_i, k_j$.
       CHARACTER*1 chan                  ! channel, 's' or 't'
       EXTERNAL prop, hspnr, spr, vpr, braket, llcor
 
       z = vpr(Pp, Pp) - vpr(Pp, Qp)
       same = 4.0D0 * prop(z, lam, -lamp, chan)
       same = same * llcor(P, Q, K1, K2)
 
       IF (rho1 .EQ. rho2) THEN
              ! equal photon helicity case:
 
               CALL hspnr( H, Hhat,  P,  Q,   lam, rho1)
               CALL hspnr(Hp, Hphat, Pp, Qp, lamp, rho2)
 
               same = same * (spr(H, Hp, rho1))**2
               same = same * spr(P, Q, rho1) * spr(Qp, Pp, -rho1)
               same = same/(spr(P, K1, rho1) * spr(Q, K1, rho1))
               same = same/(spr(P, K2, rho2) * spr(Q, K2, rho2))
           ELSE
              ! opposite photon helicity case:
 
               CALL hspnr(Ki,    Kj, K1, K2,  lam, rho1)
               CALL hspnr(Hp, Hphat, Pp, Qp, lamp,  lam)  ! ($\rho_i = \lambda$)
 
               Del    = vpr(K1, K2) - vpr(P, K1) - vpr(P, K2)
               Delhat = vpr(K1, K2) + vpr(Q, K1) + vpr(Q, K2)
               Del    = Del    + 0.5D0 * vpr(Pp, Pp)
               Delhat = Delhat + 0.5D0 * vpr(Pp, Pp)
               denom1 = spr(Kj, P, -lam) * spr(P, Ki, lam)
               denom2 = spr(Kj, Q, -lam) * spr(Q, Ki, lam)
 
              ! separate calculations of the three terms:
 
               term1 = spr(P, Kj, lam) * spr(Q, Hphat, -lam)
               term1 = term1 * braket(Ki, P, -1, Kj, Hp, lam)
               term1 = term1/(Del * denom1)
 
               term2 = spr(P, Hp, lam) * spr(Q, Ki, -lam)
               term2 = term2 * braket(Hphat, Q, 1, Ki, Kj, lam)
               term2 = term2/(Delhat * denom2)
 
               term3 = braket(Q, P, -1, Kj, Hp, lam)
               term3 = term3 * braket(Hphat, Q, 1, Ki, P, lam)
               term3 = term3/(denom1 * denom2)
 
               same = -same * (term1 + term2 + term3)
           END IF
       RETURN
       END ! FUNCTION same
 
 
       FUNCTION diff(P, Pp, Q, Qp, K1, K2, lam, lamp, rho1, rho2, chan)
      !---------------------------------------------------------
      ! Evaluation of the $t$ channel amplitude with the photons
      ! emitted from different lines. The electron line expressions
      ! from my notes are used here. The positron line case may be
      ! obtained from the same function using double crossing. Only
      ! one ordering of the two photons is included, so this function
      ! must be summed explicitly over their order.
      !---------------------------------------------------------
       IMPLICIT NONE
       INTEGER lam, lamp, rho1, rho2      ! helicities
       REAL*8 z, vpr, llcor     ! z = argument of propagator ($t_1$ in my notes)
       COMPLEX*16 diff, P(4), Pp(4), Q(4), Qp(4), K1(4),
     &            K2(4), H(4), Hp(4), Hhat(4), Hphat(4),
     &            prop, spr, braket
      !  H, Hp, Hhat, Hphat are $h_1, h'_1, \hat{h}_1, \hat{h}'_1$ in my notes.
       CHARACTER*1 chan          ! The channel, 's' or 't'
       EXTERNAL prop, hspnr, spr, vpr, braket, llcor
 
       CALL hspnr(H, Hhat, P, Q, lam, rho1)
       CALL hspnr(Hp, Hphat, Pp, Qp, lamp, rho1)
       z = vpr(Q, K1) - vpr(P, K1) - vpr(P, Q) + vpr(Pp, Pp)
       diff = 4.0D0 * prop(z, lam, -lamp, chan)
       diff = diff * llcor(P, Q, K1, K2)
 
       IF (rho1 .EQ. rho2) THEN
              ! equal photon helicity case:
               diff = diff * z * (spr(H, Hp, rho1))**2
               diff = diff/(spr( P, K1, rho1) * spr( Q, K1, rho1))
               diff = diff/(spr(Pp, K2, rho2) * spr(Qp, K2, rho2))
           ELSE
              ! opposite photon helicity case:
               diff = -diff *
     &              (braket(Hphat, K1, (lam * rho1), Hhat, H, rho1))**2
               diff = diff/(spr( P, K1, rho1) * spr( Q, K1, rho1))
               diff = diff/(spr(Pp, K2, rho2) * spr(Qp, K2, rho2))
           END IF
       RETURN
       END ! FUNCTION diff
 
 
       FUNCTION llcor(P, Q, K1, K2)
      !------------------------------------------------------------
      ! The leading log mass correction factor, translated into the
      ! electron line case by double crossing. A square root is
      ! taken since this is included in the amplitudes. Empirically,
      ! this restores the ratio of the exact and leading log results
      ! to what it would be if neither of them included explicit
      ! mass correction ($W$) factors.
      !------------------------------------------------------------
       IMPLICIT NONE
       INTEGER DCSkey, LLkey, Zkey
       COMMON /trmkey/ DCSkey, LLkey, Zkey
       REAL*8 llcor, T, alpha(2), beta(2), massqr, vpr, masfac
       COMPLEX*16 P(4), Q(4), K1(4), K2(4)
       LOGICAL getkey
       EXTERNAL vpr, masfac, getkey
       SAVE /trmkey/
 
       IF (getkey(DCSkey, 7)) THEN
          T = -vpr(P, Q)
          massqr = vpr(P, P)
          alpha(1) = vpr(P, K1)/T
          alpha(2) = vpr(P, K2)/T
          beta(1) = -vpr(Q, K1)/T
          beta(2) = -vpr(Q, K2)/T
          T = T + vpr(P, P)
          llcor = DSQRT(masfac(massqr, T, alpha, beta))
       ELSE
          llcor = 1.0D0
       END IF
       RETURN
       END ! FUNCTION llcor
 
 
***********************   SOFT PHOTON LIMIT  ****************************
 
 
       FUNCTION soft(P, Pp, Q, Qp, K1, K2)
      !----------------------------------------------------------------
      ! Soft photon approximation to the averaged/summed squared amplitude,
      ! The two terms correspond to positron line and electron line emission.
      ! The positron line case is from my notes and the electron line case
      ! is obtained by double crossing. An overall factor $e^8$ is omitted.
      !----------------------------------------------------------------
       IMPLICIT NONE
       INTEGER DCSkey, LLkey, Zkey
       COMMON /trmkey/ DCSkey, LLkey, Zkey
       REAL*8  soft, sofpos
       COMPLEX*16 P(4), Pp(4), Q(4), Qp(4), K1(4), K2(4),
     &            PX(4), PpX(4), QX(4), QpX(4)
       LOGICAL getkey
       EXTERNAL sofpos, cross, getkey
 
       CALL cross(QX, PpX, Q, Pp)
       CALL cross(PX, QpX, P, Qp)
 
       soft = 0.0D0
       IF (getkey(LLkey, 6))
          ! positron line
     &     soft = soft + sofpos( P,  Pp,  Q,  Qp, K1, K2, 't')
       IF (getkey(LLkey, 5))
          ! electron line
     &     soft = soft + sofpos(PX, PpX, QX, QpX, K1, K2, 't')
       IF (getkey(LLkey, 3))
          ! final state
     &     soft = soft + sofpos( P, PpX, QX,  Qp, K1, K2, 's')
       IF (getkey(LLkey, 2))
          ! initial state
     &     soft = soft + sofpos(PX,  Pp,  Q, QpX, K1, K2, 's')
       IF (getkey(LLkey, 4))
          ! mixed
 
     &     soft = soft - sofpos(Pp,   P,  Q,  Qp, K1, K2, 'u')
       IF (getkey(LLkey, 1))
          ! mixed
     &     soft = soft - sofpos( P,  Pp, Qp,   Q, K1, K2, 'u')
 
       soft = DABS(soft)
       RETURN
       END ! FUNCTION soft
 
 
       FUNCTION sofpos(P, Pp, Q, Qp, K1, K2, chan)
      !----------------------------------------------------------------
      ! The positron line contribution to the soft photon averaged/summed
      ! squared amplitude. The invariants are as follows:
      !    T     = $(p-q)^2 = t $
      !    Tprm  = $(p'-q')^2 = t' $
      !    U     = $(p-q')^2 = u $
      !    Uprm  = $(p'-q)^2 = u' $
      !    Tprm0 = $t' - 2m_e^2 = t'_0 $
      !    alpha(i) = $\alphat_i = -2(q'\cdot k_i)/t'_0 $
      !    beta(i)  = $\betat_i = 2(p'\cdot k_i)/t'_0 $
      !
      ! 'chan' is a channel argument in the Born amplitude
      ! evaluation: 't', 's', or 'u', describing what channel the
      ! radiation is in. See FUNCTION M2born.
      !----------------------------------------------------------------
       IMPLICIT NONE
       INTEGER DCSkey, LLkey, Zkey
       COMMON /trmkey/ DCSkey, LLkey, Zkey
       REAL*8 sofpos, S, Shat, T, Tprm, Tprm0, U, Uprm,
     &        alpha(2), beta(2), massqr,
     &        denom, M2born, vpr, func, masfac, BornS
       COMPLEX*16 P(4), Pp(4), Q(4), Qp(4), K1(4), K2(4)
       LOGICAL getkey
       CHARACTER*1 chan
       EXTERNAL M2born, vpr, func, masfac, BornS, getkey
 
       massqr = vpr(Pp, Pp)   ! $ = 2 m_e^2 $
       S     = massqr + vpr( P, Pp)
       Shat  = massqr + vpr( Q, Qp)
       T     = massqr - vpr( P,  Q)
       Tprm0 = -vpr(Pp, Qp)
       Tprm  = massqr + Tprm0
       U     = massqr - vpr( P, Qp)
       Uprm  = massqr - vpr(Pp,  Q)
 
       alpha(1) = -vpr(Qp, K1)/Tprm0
       alpha(2) = -vpr(Qp, K2)/Tprm0
        beta(1) =  vpr(Pp, K1)/Tprm0
        beta(2) =  vpr(Pp, K2)/Tprm0
 
       sofpos = 8.0D0 * M2born(S, Tprm, massqr, chan)
       denom = Tprm0 * alpha(2) * beta(1)
       sofpos = sofpos/denom
       denom = Tprm * alpha(1) * beta(2)
       sofpos = 2.0D0 * sofpos/denom
       sofpos = sofpos * (Tprm/Tprm0)
       IF (getkey(LLkey, 7))
     &    sofpos = sofpos * masfac(massqr, Tprm, alpha, beta)
       RETURN
       END ! FUNCTION sofpos
 
 
 
*********************** LEADING LOG APPROXIMATION ***********************
 
 
       FUNCTION ldlog(P, Pp, Q, Qp, K1, K2)
      !----------------------------------------------------------------
      ! Leading log approximation to the averaged/summed squared amplitude,
      ! The two terms correspond to positron line and electron line emission.
      ! The positron line case is from my notes and the electron line case
      ! is obtained by double crossing. An overall factor $e^8$ is omitted.
      !----------------------------------------------------------------
       IMPLICIT NONE
       INTEGER DCSkey, LLkey, Zkey
       COMMON /trmkey/ DCSkey, LLkey, Zkey
       REAL*8  ldlog, llpos
       COMPLEX*16 P(4), Pp(4), Q(4), Qp(4), K1(4), K2(4),
     &            PX(4), PpX(4), QX(4), QpX(4)
       LOGICAL getkey
       EXTERNAL llpos, cross, getkey
       SAVE /trmkey/
 
       CALL cross(QX, PpX, Q, Pp)
       CALL cross(PX, QpX, P, Qp)
 
       ldlog = 0.0D0
       IF (getkey(LLkey, 6))
     &    ldlog = ldlog + llpos( P,  Pp,  Q,  Qp, K1, K2, 't')
       IF (getkey(LLkey, 5))
     &    ldlog = ldlog + llpos(PX, PpX, QX, QpX, K1, K2, 't')
       IF (getkey(LLkey, 3))
     &    ldlog = ldlog + llpos( P, PpX, QX,  Qp, K1, K2, 's')
       IF (getkey(LLkey, 2))
     &    ldlog = ldlog + llpos(PX,  Pp,  Q, QpX, K1, K2, 's')
       IF (getkey(LLkey, 4))
     &    ldlog = ldlog - llpos(Pp,   P,  Q,  Qp, K1, K2, 'u')
       IF (getkey(LLkey, 1))
     &    ldlog = ldlog - llpos( P,  Pp, Qp,   Q, K1, K2, 'u')
       ldlog = DABS(ldlog)
       RETURN
       END ! FUNCTION ldlog
 
 
       FUNCTION llpos(P, Pp, Q, Qp, K1, K2, chan)
      !----------------------------------------------------------------
      ! The positron line contribution to the leading log averaged/summed
      ! squared amplitude. The invariants are as follows:
      !    T     = $(p-q)^2 = t $
      !    Tprm  = $(p'-q')^2 = t' $
      !    Tprm0 = $t' - 2m_e^2 = t'_0 $
      !    alpha(i) = $\alphat_i = -2(q'\cdot k_i)/t'_0 $
      !    beta(i)  = $\betat_i = 2(p'\cdot k_i)/t'_0 $
      !    Sborn1 = ${\tilde{s}}_1^B $
      !    Sborn2 = ${\tilde{s}}_2^B $
      !
      ! 'chan' is a channel argument in the Born amplitude
      ! evaluation: 't', 's', or 'u', describing what channel the
      ! radiation is in. See FUNCTION M2born.
      !----------------------------------------------------------------
       IMPLICIT NONE
       INTEGER DCSkey, LLkey, Zkey
       COMMON /trmkey/ DCSkey, LLkey, Zkey
       REAL*8 llpos, S, Shat, T, Tprm, Tprm0, U, Uprm,
     &        alpha(2), beta(2), massqr, Sborn1, Sborn2,
     &        denom, M2born, vpr, func, masfac, BornS
       COMPLEX*16 P(4), Pp(4), Q(4), Qp(4), K1(4), K2(4)
       CHARACTER*1 chan
       LOGICAL getkey
       EXTERNAL M2born, vpr, func, masfac, BornS, getkey
       SAVE /trmkey/
 
       massqr = vpr(Pp, Pp)   ! $ = 2 m_e^2 $
       S     = massqr + vpr( P, Pp)
       Shat  = massqr + vpr( Q, Qp)
       T     = massqr - vpr( P,  Q)
       Tprm0 = -vpr(Pp, Qp)
       Tprm  = massqr + Tprm0
       U     = massqr - vpr( P, Qp)
       Uprm  = massqr - vpr(Pp,  Q)
 
       alpha(1) = -vpr(Qp, K1)/Tprm0
       alpha(2) = -vpr(Qp, K2)/Tprm0
        beta(1) =  vpr(Pp, K1)/Tprm0
        beta(2) =  vpr(Pp, K2)/Tprm0
 
       Sborn1 = BornS(S   , T, Uprm, massqr)
       Sborn2 = BornS(Shat, T, U   , massqr)
 
       llpos = 1D0
       llpos = func(alpha,beta,alpha) * M2born(Sborn1, T, massqr, chan)
     &       + func(alpha,beta, beta) * M2born(Sborn2, T, massqr, chan)
       denom = Tprm0 * alpha(2) * beta(1)
       llpos = llpos/denom
       denom = T * alpha(1) * beta(2)
       llpos = 2.0D0 * llpos/denom
       llpos = llpos * (Tprm/Tprm0)
       IF (getkey(LLkey, 7))
     &    llpos = llpos * masfac(massqr, Tprm, alpha, beta)
       RETURN
       END ! FUNCTION llpos
 
 
       FUNCTION func(alpha, beta, gamma)
      !--------------------------------------------
      ! The function func is related to functions appearing
      ! in the leading log expressions as follows:
      ! $ f_1(\alpha, \beta) = $func$(\alpha, \beta, \alpha),$
      ! $ f_2(\alpha, \beta) = $func$(\alpha, \beta, \beta). $
      !--------------------------------------------
       IMPLICIT NONE
       INTEGER i, j
       REAL*8 func, alpha(2), beta(2), gamma(2), u, v, x, Y
      ! Define the function $Y$ in the Leading Log notes:
       Y(x, u, v) = (1-x)**2 * ((1-u)**2 + (1-v)**2)
 
       IF (alpha(1) + beta(1) .GE.
     &     alpha(2) + beta(2)) THEN
           i = 1
           j = 2
       ELSE
           i = 2
           j = 1
       END IF
 
       u = alpha(j)/(1.0D0 - alpha(i))
       v =  beta(j)/(1.0D0 -  beta(i))
       x = gamma(i)/(1.0D0 - gamma(j))
       func = Y(gamma(i), u, v) + Y(x, alpha(j), beta(j))
       RETURN
       END ! FUNCTION func
 
 
       FUNCTION masfac(massqr, Tprm, alpha, beta)
      !------------------------------------------------------
      ! The mass correction factor
      !    $ W(\alphat_1, \betat_1) W(\alphat_2, \betat_2) $
      ! in the leading log expressions of BFLW et al, crossed
      ! to the $t$ channel.
      !------------------------------------------------------
       IMPLICIT NONE
       INTEGER index
       REAL*8 masfac, massqr, Tprm, alpha(2), beta(2), ratio, Wfac
 
       masfac = 1.0D0
       DO index = 1, 2
          ratio = beta(index)/alpha(index)
          ratio = (1.0D0 - alpha(index))/(1.0D0 - beta(index))
          Wfac = ratio + 1.0D0/ratio
          ratio = beta(index)/alpha(index)
          Wfac = (ratio + 1.0D0/ratio)/Wfac
          Wfac = 1.0D0 - Wfac * massqr/Tprm
          masfac = masfac * Wfac
       END DO
       RETURN
       END ! FUNCTION masfac
 
 
******************* BASIC CROSS SECTION FUNCTIONS ***********************
 
 
       FUNCTION csnorm(P, Pp, Q, Qp, K1, K2)
      !-----------------------------------------------------------
      ! The normalization factor for the cross-section as defined by
      ! B.F.L.W. et al.  This is the coefficient of the averaged/summed
      ! squared matrix element. In my notes, $F_n$ represents the cross
      ! section for $n$ photons, up to some energy factors. A naive
      ! approximation $ F_{naive} =  \beta^2 F_0/(E_B (t')^2), $ to $F_2$
      ! is divided out, giving the normalized differential cross section
      !          DCS = $ F_2/F_{naive}. $
      ! The Born kinematic parameters for $F_0$ are the same as those in
      ! `FUNCTION llpos', and $\beta = \pi^{-2}\log[(|t'|/m_e^2) -1].$
      ! As always, factors of $e^2 = 4\pi\alpha$ are omitted. The final
      ! normalized cross-section is dimensionless and independent of $e$.
      !
      ! The parameter `chan' is 's' or 't', depending on which channel
      ! is used in the big logarithm. 's' is used only if all of the
      ! $T$ channel terms in DCSkey are zero.
      !-----------------------------------------------------------
       IMPLICIT NONE
       INTEGER DCSkey, LLkey, Zkey
       COMMON /trmkey/ DCSkey, LLkey, Zkey
       REAL*8 csnorm, EB, EP, S, Shat, Sborn, T, Tprm, Uprm,
     &        F2coef, F0, Fnaive, pi, massqr, betaf,
     &        Lfunc, M2born, energy, vpr, BornS
       COMPLEX*16 P(4), Pp(4), Q(4), Qp(4), K1(4), K2(4)
       LOGICAL getkey
       CHARACTER*1 chan
       EXTERNAL Lfunc, M2born, energy, vpr, BornS, getkey
       SAVE /trmkey/
 
       pi = 4.0D0 * DATAN(1.0D0)
       massqr = vpr(Pp, Pp)   !  $ = 2 m_e^2 $
       EP = energy(Qp)
       EB = energy(Pp)
       S    = massqr + vpr( P, Pp)
       Shat = massqr + vpr( Q, Qp)
       T    = massqr - vpr( P,  Q)
       Tprm = massqr - vpr(Pp, Qp)
       Uprm = massqr - vpr(Pp,  Q)
 
       IF (getkey(DCSkey, 6) .OR. getkey(DCSkey, 5)
     &     .OR. getkey(DCSkey, 4)) THEN
            chan = 't'
            betaf = Lfunc(Pp, Qp)/(pi**2)
            Fnaive = betaf**2/(EB * Tprm**2)
       ELSE
            chan = 's'
            betaf = Lfunc(Pp, P)/(pi**2)
            Fnaive = betaf**2/(EB * S**2)
       END IF
 
       Sborn = BornS(S, T, Uprm, massqr)
       F0 = M2born(Sborn, T, massqr, chan)
       F0 = F0/(32.0D0 * (pi**2) * Sborn)
       Fnaive = Fnaive * F0
       F2coef = EP/((2.0D0**13) * (pi**8) * S * Shat)
       csnorm = F2coef/Fnaive
       RETURN
       END ! FUNCTION csnorm
 
 
       FUNCTION Lfunc(spnr1, spnr2)
      !------------------------------------------------------
      ! The logarithm function which appears in the leading
      ! log expansion. It is evaluated using the vector product
      ! of two spinors to find the appropriate invariant, $t$ or $t'$.
      !------------------------------------------------------
       IMPLICIT NONE
       REAL*8 Lfunc, emass, vpr, invar
       COMPLEX*16 spnr1(4), spnr2(4)
       PARAMETER (emass = 0.51099906D-3)
       EXTERNAL vpr
 
       invar = vpr(spnr1, spnr1) + vpr(spnr2, spnr2)
       invar = DABS(0.5D0 * invar  - vpr(spnr1, spnr2))
       Lfunc = DLOG(invar/(emass**2)) - 1.0D0
       RETURN
       END ! FUNCTION Lfunc
 
 
       FUNCTION BornS(S, T, U, massqr)
      !---------------------------------------------------------
      ! The effective Born $s$ parameter appearing in leading log
      ! and normalization factors. It is expressed in terms of
      ! parameters S, T, U which may be specified to give $s^B_1$
      ! or $s^B_2$. See FUNCTION llpos. Here, massqr is $2m_e^2.$
      ! If the arguments S, T are interchanged, this function gives
      ! the Born $t$ parameter in the $s$ channel.
      !---------------------------------------------------------
       IMPLICIT NONE
       REAL*8 BornS, S, T, U, massqr, ratio, EPsqr, vborn2
      ! EPsqr is twice the square of the incoming positron or electron
      !    energy in the CMS frame of the outgoing fermions.
      ! The following are mass-correction velocity factors:
       EPsqr = (S + U - 2.0D0 * massqr)**2/(2.0D0 * T)
       vborn2 = 1.0D0 - 2.0D0 * massqr/T
       ratio = DSQRT(vborn2/(1.0D0 - massqr/EPsqr))
 
       BornS = (U - S)/(U + S - 2.0D0 * massqr)
       BornS = 0.5D0 * T * (BornS * ratio - vborn2)
       RETURN
       END ! FUNCTION BornS
 
 
       FUNCTION M2born(Sborn, Tborn, massqr, chan)
      !----------------------------------------------------
      ! This function returns the averaged and summed squared
      ! Born amplitude for Bhabha scattering in the massless
      ! limit, with the overall factor of $e^4$ omitted. Helicity
      ! conservation determines which channels contribute to
      ! the amplitude, as described in `FUNCTION ampsqr'.
      ! The third parameter is a character indicating which
      ! channel the second parameter is in. If 'chan' = 's', the
      ! Sborn and Tborn parameters must be interchanged. If
      ! 'chan' = 'u', then Tborn is really the 'u' parameter.
      ! `massqr' is $2m_e^2$ and is used purely for kinematics.
      !----------------------------------------------------
       IMPLICIT NONE
       INTEGER lam          ! helicity $\lambda$
       REAL*8 M2born, Sborn, Tborn, pureS, pureT, mixed, absqr,
     &        S, T, massqr
       COMPLEX*16 STamp, prop
       CHARACTER*1 chan     ! channel: 's', 't', or 'u'
       EXTERNAL prop, absqr
 
       pureS = 0.0D0
       pureT = 0.0D0
       mixed = 0.0D0
 
c[[[  Old version: Whole born amplitude calculated everywhere:
c      IF (chan .EQ. 't') THEN
c          T = Tborn
c          S = Sborn
c      ELSE IF (chan .EQ. 's') THEN
c          T = Sborn
c          S = Tborn
c      ELSE  !  IF (chan .EQ. 'u') THEN
c          T = 2.0D0 * massqr - Sborn - Tborn
c          S = Sborn
c      END IF
c          STamp = (0.0D0, 0.0D0)   ! STamp includes both channels
c
c          pureS = pureS + absqr(prop(S, lam, -lam, 's'))
c          STamp = STamp + prop(S, lam, lam, 's')
c
c          pureT = pureT + absqr(prop(T, lam, -lam, 't'))
c          STamp = STamp + prop(T, lam, lam, 't')
c
c          mixed = mixed + absqr(STamp)
c      END DO ! end loop over helicity $\lambda$
c]]]
 
c[[[  New version: only specified channel of Born amplitude included:
c          (That channel is labelled `T' , but in fact depends on
c           the arguments given. Crossing may have been used.)
       T = Tborn
       S = Sborn
       STamp = (0.0D0, 0.0D0)
 
       DO lam = -1, 1, 2  ! loop over helicity $\lambda$:
           pureT = pureT + absqr(prop(T, lam, -lam, chan))
           mixed = mixed + absqr(prop(T, lam,  lam, chan))
       END DO ! end loop over helicity $\lambda$
c]]]
 
       pureS = pureS * (T)**2
       pureT = pureT * (S)**2
       mixed = mixed * (S + T - 2.0D0 * massqr)**2
       M2born = pureS + pureT + mixed
 
       RETURN
       END ! FUNCTION M2born
 
 
       FUNCTION prop(z, hel1, hel2, chan)
      !------------------------------------------------------------
      ! The photon/Z propagator function $G_{\lambda,\mu}^{(\prime)}(z)$
      ! in my notes. The width is taken to be the $s$-dependent width
      ! $s\Gamma_Z/M_Z$ in the $s$ channel, and 0 in the $t$ channel, to
      ! give $G$ or $G'$, respectively. Here, $\lambda, \mu$ = hel1, hel2.
      ! The propagator is symmetric in the helicities.
      !
      ! fundamental parameters:
      !       Zmass  = mass of Z ($M_Z$)
      !       Zwidth = width of Z ($\Gamma_Z$)
      !       sinW2  = $ \sin^2(\theta_W) $
      !------------------------------------------------------------
       IMPLICIT NONE
       INTEGER DCSkey, LLkey, Zkey, hel1, hel2
       COMMON /trmkey/ DCSkey, LLkey, Zkey
       REAL*8 z, width, numer, Zmass, Zwidth, sinW2
       PARAMETER (Zmass=91.161D0, Zwidth=2.534D0, sinW2=0.2259D0)
       CHARACTER*1 chan         ! channel: 's' or 't'
       COMPLEX*16 prop, denom, i
       LOGICAL getkey
       PARAMETER (i = (0.0D0, 1.0D0))
       EXTERNAL getkey
       SAVE /trmkey/
 
       prop = 1.0D0/z
 
       IF (getkey(Zkey, 1)) THEN
         ! Include $Z^0$ exchange --
         ! The width is nonzero only in the $s$ channel:
          IF (chan .EQ. 's') THEN
              width = Zwidth/Zmass
          ELSE
              width = 0.0D0
          END IF
 
          numer = ((1 - hel1) - 4.0D0 * sinW2) *
     &            ((1 - hel2) - 4.0D0 * sinW2)
          numer = numer/(16.0D0 * sinW2 * (1.0D0 - sinW2))
          denom = z * (1.0D0 + i * width) - Zmass**2
          prop = prop + numer/denom
       END IF
       RETURN
       END ! FUNCTION prop
 
 
*********************** BASIC SPINOR OPERATIONS *************************
 
 
       FUNCTION angle(spnr1, spnr2)
      !----------------------------------------
      ! This function returns the angle
      ! between two spinors. Here, `sinsqr' is
      ! the squared sine of the angle. The spinors
      ! are not necessarily assumed to have $E$ > $m$.
      ! This function is for diagnostic use only.
      ! The result is single precision, in degrees.
      ! Intermediate calculations are double precision.
      !----------------------------------------
       IMPLICIT NONE
       REAL*4 angle   ! SINGLE PRECISION!
       REAL*8 vpr, pi
       COMPLEX*16 spnr1(4), spnr2(4), mcor1, mcor2, p1, p2,
     &            shift, sinsqr, momtum
       EXTERNAL vpr, momtum
 
       p1 = momtum(spnr1)
       p2 = momtum(spnr2)
       mcor1 = spnr1(4)
       mcor2 = spnr2(4)
 
       shift = p1 * mcor2 + p2 * mcor1 + mcor1 * mcor2
       sinsqr = vpr(spnr1, spnr2) - 2.0D0 * shift
       sinsqr = 0.25D0 * sinsqr/(p1 * p2)
       angle = 2.0D0 * DASIN(DSQRT(DREAL(sinsqr)))
       pi = 4.0D0 * DATAN(1.0D0)
       angle = angle * 180.0D0/pi
       RETURN
       END ! FUNCTION angle
 
 
       FUNCTION vpr(vec1, vec2)
      !--------------------------------------------------
      ! Defines the doubled vector product of vec1 and vec2,
      ! i.e. $2 v_1 \cdot v_2 = <v_1,-|v_2,+> <v_2,+|v_1,->.$
      ! for massless vectors $v_1, v_2$ with Weyl spinor
      ! representations vec1, vec2. This subroutine is valid
      ! for vectors on the forward or backward light cone.
      !--------------------------------------------------
       IMPLICIT NONE
       REAL*8 vpr
       COMPLEX*16 vec1(4), vec2(4), spr
       EXTERNAL spr
 
       vpr = DREAL(spr(vec1, vec2, 1) * spr(vec2, vec1, -1))
       RETURN
       END ! FUNCTION vpr
 
 
       FUNCTION braket(P1, P2, c, P3, P4, mu)
      !-----------------------------------------------------
      ! A function appearing in the opposite helicity amplitudes,
      ! depending on four spinors, an integer $c$, and a helicity
      ! $\mu$. The function is the `triple product' in my notes:
      ! $ <p_1, (p_2 + c p_3), p_4>_\mu. $
      !-----------------------------------------------------
       IMPLICIT NONE
       INTEGER c, mu
       COMPLEX*16 braket, P1(4), P2(4), P3(4), P4(4), spr
       EXTERNAL spr
 
       braket =      c * spr(P1, P3, -mu) * spr(P3, P4, mu)
       braket = braket + spr(P1, P2, -mu) * spr(P2, P4, mu)
       RETURN
       END ! FUNCTION braket
 
 
       FUNCTION spr(spnr1, spnr2, hel)
      !---------------------------------------------------------
      ! Defines the spinor product $<p, -\lambda | q, \lambda>,$ where
      ! hel = $\lambda$ and spnr1 = $P$, spnr2 = $Q$ are the Weyl
      ! spinors related to $p$ and $q$ as in `SUBROUTINE spinor'.
      ! Including a mass such that FUNCTION vpr reproduces a massive
      ! vector product, but the algebraic properties of spr are
      ! preserved, I define
      !   $ <p,-|q,+> = \sigma(1 + \epsilon/\sqrt{\sigma\sigma'}), $
      !   $ <p,+|q,-> = \sigma'(1 - \epsilon/\sqrt{\sigma\sigma'}) $
      ! where
      ! $ \sigma = P_2 Q_1 - P_1 Q_2,   \sigma' = Q_3 P_1 - Q_1 P_3, $
      ! $ \epsilon = \sqrt{2 E_1 E_2 (1 - v_1 v_2)}. $
      ! When $p$ and $q$ are parallel, a phase singularity requires
      ! special care. The parallel limit is, by definition,
      ! $ <p,-|q,+> = ie^{i\phi}\epsilon, <p,+|q,-> = -ie^{-i\phi}\epsilon. $
      ! This is continuous to the limiting value with equal $\phi$'s
      ! or small $\theta$'s.
      !---------------------------------------------------------
       IMPLICIT NONE
       INTEGER hel
       COMPLEX*16 spr, spnr1(4), spnr2(4), mcor1, mcor2, i,
     &            epsiln, sigma, sigmaP,   ! $ (\epsilon, \sigma, \sigma') $
     &            p1, p2, momtum
       PARAMETER (i = (0.0D0, 1.0D0))
       EXTERNAL momtum
 
       p1 = momtum(spnr1)
       p2 = momtum(spnr2)
       mcor1 = spnr1(4)
       mcor2 = spnr2(4)
 
       epsiln = p1 * mcor2 + p2 * mcor1 + mcor1 * mcor2
       epsiln = CDSQRT(-2.0D0 * i * epsiln) * CDSQRT(i)
      ! The branch cut is moved to the imaginary axis to stabilize
      ! the phase when evaluated for negative energy spinors.
       sigma  = spnr1(2) * spnr2(1) - spnr1(1) * spnr2(2)
       sigmaP = spnr1(1) * spnr2(3) - spnr1(3) * spnr2(1)
 
       IF (CDABS(sigma) .EQ. 0.0D0) THEN  ! parallel case
          IF (CDABS(spnr1(2)) .EQ. 0.0D0) THEN
             spr = i
          ELSE   ! The following is $ ie^{i\phi}$. The definition
                 ! is chosen to give the correct branch in all cases.
             spr = i * CDSQRT(spnr1(2))/CDSQRT(spnr1(3))
          END IF
          IF (hel .EQ. -1) spr = DCONJG(spr)
          spr = spr * epsiln
       ELSE           ! normal case
          epsiln = epsiln/CDSQRT(sigma * sigmaP)
          IF (hel .EQ. 1) THEN
             spr = sigma * (1.0D0 + epsiln)
          ELSE
             spr = sigmaP * (1.0D0 - epsiln)
          END IF
       END IF
       RETURN
       END ! FUNCTION spr
 
 
       SUBROUTINE spnvec(var, vec, spnr, theta, phi)
      !---------------------------------------------------------
      ! Defines a spinor representation of a given four-vector.
      ! The algorithm is chosen for stability at small angles and masses.
      ! For this reason, it is convenient to define the spinor by first
      ! deriving the polar representation of the vector, and then using
      ! SUBROUTINE spinor below. This subroutine is useful in the
      ! interface between the spinor functions above and a Monte-Carlo
      ! program such as BHLUMI, which generates four-vectors.
      ! The polar angles are also returned, for possible
      ! external use. The vector vec is to be given in components
      ! (1, 2, 3, 4 = time).
      !---------------------------------------------------------
       IMPLICIT NONE
       REAL*8 vec(4), perp, theta, phi, massqr
       COMPLEX*16 spnr(4)
       CHARACTER*2 var
       EXTERNAL spinor, compar
 
       massqr = vec(4)**2 - vec(3)**2 - vec(2)**2 - vec(1)**2
       perp = DSQRT(vec(1)**2 + vec(2)**2)
       theta = DATAN2(perp, vec(3))        ! always between 0 and $\pi$.
       IF (perp .EQ. 0.0D0) THEN
          phi = 0.0D0                      ! by convention
       ELSE
          phi = DATAN2(vec(2), vec(1))     ! between $-\pi$ and $\pi$.
       END IF
 
       CALL spinor(spnr, vec(4), theta, phi, massqr)
       CALL compar(var, vec, spnr)
       RETURN
       END ! SUBROUTINE spnvec
 
 
       SUBROUTINE compar(var, vec, spnr)
      !--------------------------------------------------------
      ! Diagnostic subroutine: Compares given vector to that
      ! derived from the associated spinor, prints difference
      ! if it exceeds a specified amount, "small".
      !--------------------------------------------------------
       IMPLICIT NONE
       INTEGER index
       REAL*8 vec(4), diff(4), small
       PARAMETER (small = 1.0D-30)
       COMPLEX*16 spnr(4), vec2(3)
       CHARACTER*2 var
       EXTERNAL vector
 
       CALL vector(vec2, spnr)
       diff(4) = 0.5D0 * DREAL(vec2(1) + vec2(2)) - vec(4)
       diff(3) = 0.5D0 * DREAL(vec2(1) - vec2(2)) - vec(3)
       diff(1) = DREAL(vec2(3)) - vec(1)
       diff(2) = DIMAG(vec2(3)) - vec(2)
cc       DO index = 1, 4
cc          IF (DABS(diff(index)) .GT. small)
cc     &      WRITE(6,*) var, ': diff(', index, ') = ', diff(index)
cc       END DO
       RETURN
       END ! SUBROUTINE compar
 
 
       SUBROUTINE spinor(spnr, E, theta, phi, massqr)
      !---------------------------------------------------------
      ! Defines a Weyl spinor representation of a
      ! momentum, following the conventions of my notes:
      ! Let $P_1$ and $P_2$ be the nonzero components of a massless
      ! Dirac spinor $|p,+>$. Explicitly, they are
      !          $ P_1 = \sqrt{2vE} \cos(\theta/2), $
      !          $ P_2 = \sqrt{2vE} \sin(\theta/2) e^{i\phi}. $
      ! These form the first two entries of a complex array $P$.
      ! Here, $v = p/E$ incorporates mass effects, in a purely
      ! kinematic sense. (Dynamically, zero mass is assumed.)
      ! The energy may be either positive or negative. The third component
      !          $ P_3 = \sqrt{2Ev} \sin(\theta/2) e^{-i\phi} $
      ! for convenience. The mass correction $E(1 - v)$ is included as $P_4$.
      !---------------------------------------------------------
       IMPLICIT NONE
       REAL*8 E, theta, phi, massqr
       COMPLEX*16 spnr(4), i, mascor, mterm, root
       PARAMETER (i = (0.0D0, 1.0D0))
       EXTERNAL mascor
 
       mterm = mascor(E, massqr)
       root = 2.0D0 * (E - mterm)
       root = CDSQRT(root)
       spnr(1) = root * DCOS(0.5D0 * theta)
       spnr(2) = root * DSIN(0.5D0 * theta) * CDEXP( i * phi)
       spnr(3) = root * DSIN(0.5D0 * theta) * CDEXP(-i * phi)
       spnr(4) = mterm
       RETURN
       END ! SUBROUTINE spinor
 
 
       SUBROUTINE vector(vctr, spnr)
      !------------------------------------------------
      ! Defines the lightcone and complex transverse
      ! vector components `vctr' of a spinor `spnr'.
      ! A spinor $P$ is related to a vector $p$ whose
      ! components satisfy
      ! $ p_+ = (P_1)^2 + P_4,     p_- = P_2 P_3 + P_4, $
      !           $  p_{\perp} = P_1 P_2.$
      ! Here, $P_4 = E(1 - v).$
      ! There is no assumption about the energy.
      ! The components are:
      ! vctr(1) = $p_+$, vctr(2) = $p_-$, vctr(3) = $p_{\perp}$.
      !------------------------------------------------
       IMPLICIT NONE
       COMPLEX*16 spnr(4), vctr(3)
 
       vctr(1) = spnr(4) + spnr(1)**2            ! $p_+$
       vctr(2) = spnr(4) + spnr(2) * spnr(3)     ! $p_-$
       vctr(3) = spnr(1) * spnr(2)               ! $p_{\perp}$
       RETURN
       END ! SUBROUTINE vector
 
 
       SUBROUTINE cross(spnr1X, spnr2X, spnr1, spnr2)
      !--------------------------------------------------
      ! Carries out a crossing transformation on two spinors
      ! spnr1, spnr2 giving crossed spinors spnr1X, spnr2X.
      ! Either positive or negative energy input spinors are
      ! allowed, and the output has the opposite sign energy.
      ! This causes the spinor components to be multiplied
      ! by a phase $i$ or $-i$.
      !--------------------------------------------------
       IMPLICIT NONE
       INTEGER index
       REAL*8 Esign1, Esign2, energy
       COMPLEX*16 spnr1(4), spnr2(4), spnr1X(4), spnr2X(4), i
       PARAMETER (i = (0.0D0, 1.0D0))
       EXTERNAL energy
 
       Esign1 = energy(spnr1)
       Esign2 = energy(spnr2)
       Esign1 = Esign1/DABS(Esign1)
       Esign2 = Esign2/DABS(Esign2)
 
       DO index = 1, 3
             spnr2X(index) = i * Esign1 * spnr1(index)
             spnr1X(index) = i * Esign2 * spnr2(index)
       END DO
       spnr2X(4) = -spnr1(4)
       spnr1X(4) = -spnr2(4)
       RETURN
       END ! SUBROUTINE cross
 
 
       SUBROUTINE hspnr(spnr1H, spnr2H, spnr1, spnr2, hel1, hel2)
      !---------------------------------------------------------------------
      ! defines helicity-dependent pair of Weyl spinors (spnr1H, spnr2H)
      ! to be (spnr1, spnr2) if the helicities hel1 and hel2 are equal, and
      ! to be (spnr2, spnr1) if they are opposite.
      !---------------------------------------------------------------------
       IMPLICIT NONE
       INTEGER hel1, hel2, index
       COMPLEX*16 spnr1H(4), spnr2H(4), spnr1(4), spnr2(4)
 
       IF (hel1 .EQ. hel2) THEN
          DO index = 1, 4
               spnr1H(index) = spnr1(index)
               spnr2H(index) = spnr2(index)
          END DO
       ELSE
          DO index = 1, 4
               spnr1H(index) = spnr2(index)
               spnr2H(index) = spnr1(index)
          END DO
       END IF
       RETURN
       END ! SUBROUTINE hspnr
 
 
       FUNCTION energy(spnr)
      !-------------------------------------
      ! Returns the energy of a spinor.
      !-------------------------------------
       IMPLICIT NONE
       REAL*8 energy
       COMPLEX*16 spnr(4), momtum
       EXTERNAL momtum
 
       energy = DREAL(momtum(spnr) + spnr(4))
       RETURN
       END ! FUNCTION energy
 
 
       FUNCTION momtum(spnr)
      !-----------------------------------------
      ! Returns the momentum $p$ of a spinor:
      ! $ p = vE = (P_1^2 + P_2 P_3). $
      !-----------------------------------------
       IMPLICIT NONE
       COMPLEX*16 momtum, spnr(4)
 
       momtum = 0.5D0 * (spnr(1)**2 + spnr(2) * spnr(3))
       RETURN
       END ! FUNCTION momtum
 
 
       FUNCTION mascor(E, massqr)
      !---------------------------------------------
      ! The difference $E - p = E(1 - v)$ for a
      ! massive particle as a function of $E$ and $m^2$.
      ! In general, this may be complex. The algorithm
      ! is chosen for numerical stability when $m$
      ! is much smaller than $E$.
      !---------------------------------------------
       IMPLICIT NONE
       REAL*8 E, massqr
       COMPLEX*16 mascor, mass, ratio, i
       PARAMETER (i = (0.0D0, 1.0D0))
 
       mass = CDSQRT(DCMPLX(massqr))
       ratio = i * (mass/E)
       mascor = CDLOG(ratio + CDSQRT(1.0D0 + ratio**2))
       mascor = CDSIN(0.5D0 * i * mascor)
       mascor = (2.0D0 * E) * (mascor**2)
       RETURN
       END ! FUNCTION mascor
 
 
       FUNCTION absqr(z)
      !-----------------------------------
      ! The absolute square function.
      !-----------------------------------
       IMPLICIT NONE
       REAL*8 absqr
       COMPLEX*16 z
 
       absqr = DREAL(z * DCONJG(z))
       RETURN
       END ! FUNCTION absqr
 
 
       FUNCTION getkey(key, place)
      !---------------------------------------
      ! Returns the logical value assigned to
      ! the specified digit place of a key:
      ! .TRUE. if nonzero, otherwise .FALSE.
      ! Places are counted from the least
      ! significant digit, which is place 1.
      !---------------------------------------
       IMPLICIT NONE
       INTEGER key, place, digit, expon
       LOGICAL getkey
 
       expon = 10**place
       digit = key - (key/expon) * expon
       expon = expon/10
       digit = digit/expon
       getkey = (digit .NE. 0)
       RETURN
       END! FUNCTION getkey            !\end
 

