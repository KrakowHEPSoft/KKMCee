*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//                     Pseudo-CLASS  BornV                                  //
*//                                                                          //
*//  Purpose:                                                                //
*//  Provide Born angular distribution and integrated x-section              //
*//  as a function of s.                                                     //
*//                                                                          //
*//  NOTES:                                                                  //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
*
      SUBROUTINE BornV_Initialize(xpar_input)
*//////////////////////////////////////////////////////////////////////////////
*//                    Class initializator                                   //
*// Notes:                                                                   //
*// This initializator should be called before any other routine of the class//
*// It defines (mostly static) class members using input from xpar matrix    //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BXformat.h'
      INCLUDE 'BornV.h'
      DOUBLE PRECISION  xpar_input(*)
*------------------------------------------------------------------------------
      DOUBLE PRECISION    amuon
      INTEGER             k,j,kxpa,KF
      DOUBLE PRECISION    vvmax
*------------------------------------------------------------------------------
      m_QCDcor = 0d0
      DO k=1,m_poinQ
         m_QCDcorR(k)=0d0
      ENDDO
      m_CMSene = xpar_input( 1)         ! Central value of CMS energy, do not change!
      m_XXXene = m_CMSene               ! Just initialization only
      m_KFini = xpar_input( 400)        ! KFcode of beam, POSITIVE!!!
*                      <<<  ff-pair spectrum >>>
      m_vvmin  = xpar_input(16)         ! minimum v, infrared cut
      vvmax    = xpar_input(17)         ! maximum v
      amuon  = 0.1056583d0
      m_vvmax  = min(vvmax, 1d0-(2*amuon/m_CMSene)**2)
      m_HadMin = xpar_input(51)         ! minimum hadronization mass
*                        <<< Basic QED >>>
      m_AlfInv = xpar_input(30)         ! Alpha_QED at Thomson limit
      m_alfpi  = 1d0/m_pi/m_AlfInv
*                  <<< Electroweak parameters >>>
      m_Gmu    = xpar_input(32)         ! Fermi constant
      m_MZ     = xpar_input(502)        ! Z mass [GeV]
      m_amh    = xpar_input(805)        ! Higgs mass, Input for Dizet
      m_amtop  = xpar_input(806)        ! Top mass,   Input for Dizet
* Note that gammz and swsq will be redefined in the case of EW corrs. are on
      m_swsq   = xpar_input(503)        ! Electroweak mixing angle
      m_gammz  = xpar_input(504)        ! Z width
      m_MW     = xpar_input(505)        ! W mass [GeV]
      m_GammW  = xpar_input(506)        ! W width[GeV]

*               <<< Static Table of ALL fermion parameters >>>
      DO j=1,20
         m_IsGenerated(j) = xpar_input(400+j)   ! Generation flag
         kxpa = 500+10*j
         m_KFferm(j)= xpar_input(kxpa+1)        ! fermion flavour code
         m_NCf(j)   = xpar_input(kxpa+2)        ! number of colours
         m_Qf(j)    = xpar_input(kxpa+3)/3d0    ! electric charge
         m_T3f(j)   = xpar_input(kxpa+4)/2d0    ! isospin, L-hand component
         m_helic(j) = xpar_input(kxpa+5)        ! helicity, polarization
         m_amferm(j)= xpar_input(kxpa+6)        ! fermion mass
         m_AuxPar(j)= xpar_input(kxpa+8)        ! auxiliary parameter
      ENDDO
*                       <<< Test switches >>>
      m_KeyElw = xpar_input(12)         ! ElectroWeak library on/off
      m_KeyZet = xpar_input(501)        ! Z-boson on/off
      m_KeyWtm = xpar_input(26)         ! Photon emission without mass terms
      m_KeyRes = xpar_input(13)         ! Exper. R for gamma*
*                       <<<  Other        >>>
      m_KeyQCD = xpar_input(53)         ! QCD FSR
      m_KeyINT = xpar_input(27)         ! This is realy copy from KK2f
      m_Xenph  = xpar_input(40)         ! This is realy copy from KK2f
      IF(m_KeyINT .EQ. 0)  m_Xenph  = 1D0
*                       <<< Miscelaneous >>>
      m_gnanob = xpar_input(31)         ! GeV^(-2) to nanobarns
*
      m_out    = xpar_input(4)
*
      WRITE(m_out,bxope)
      WRITE(m_out,bxtxt) '  BornV  Initializator                '
      WRITE(m_out,bxl1f) m_MZ    ,   'Z mass     [GeV]   ','amz   ','a1'
      WRITE(m_out,bxl1f) m_amh   ,   'Higgs mass [GeV]   ','amh   ','a2'
      WRITE(m_out,bxl1f) m_amtop ,   'Top mass   [GeV]   ','amtop ','a3'
      WRITE(m_out,bxl1f) m_gammz,    'Z width    [GeV]   ','gammz ','a4'
      WRITE(m_out,bxl1f) m_swsq,     'sin(theta_w)**2    ','sinw2 ','a5'
      WRITE(m_out,bxl1f) m_AlfInv,   '1/alfa_QED  at  Q=0','AlfInv','a6'
      WRITE(m_out,bxl1f) m_HadMin,   'MassCut light qqbar','HadMin','a6'
      WRITE(m_out,bxl1i) m_KFini ,   'KF code of beam    ','KFini ','a7'
      WRITE(m_out,bxl1g) vvmax,      'Input vvmax        ','vvmax ','a8'
      WRITE(m_out,bxl1g) m_vvmax,    'reduced vvmax in MC','vvmax ','a9'
      WRITE(m_out,bxtxt) 'Test switches:                         '
      WRITE(m_out,bxl1i) m_KeyElw,   'Electroweak lib.   ','KeyElw','10'
      WRITE(m_out,bxl1i) m_KeyZet,   'Z on/off   switch  ','KeyZet','11'
      WRITE(m_out,bxl1i) m_KeyWtm,   'mass terms on/off  ','KeyWtm','12'
      WRITE(m_out,bxl1i) m_KeyRes,   'R for gamma* on/off','KeyRes','12'
      WRITE(m_out,bxclo)
      IF( m_KeyElw .NE. 0 ) CALL BornV_StartEW(xpar_input)
      END


      SUBROUTINE BornV_ReBin1(RR,alf,bet,xmax,x,djac)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//  This mapps variable r into x.                                           //
*//  to be employed in the integration (either ordinary or Monte Carlo)      //
*//  of distributions resambling                                             //
*//  the binomial distribution x**(alf-1)*(1-x)**(bet-1)                     //
*//  with alf > 0 and  bet arbitrary.                                        //
*//  variable r is in (0,1) range and x is within (0,xmax) range.            //
*//  djac is jacobian factor d(x)/d(r).                                      //
*//  mapping is such that 1/djac is very CLOSE to                            //
*//  binomial distribution x**(alf-1)*(1-x)**(bet-1).                        //
*//  WARNING: mapping may fail very CLOSE to r=0. Practically, one is        //
*//  recommended to obey: fleps**alf < r, where fleps = 1d-100.              //
*//  problems may also arise for very small xmax ( below 1.d-12 ).           //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      DOUBLE PRECISION RR,R,alf,bet,xmax,x,djac
      DOUBLE PRECISION x0,dist,r1,p1,q1,q2
*------------------------------------------------------------------------------
      IF( alf .LE. 0d0 ) GOTO 900
      R = MAX(RR,m_fleps**alf)
      x0=(alf-1d0)/(alf+bet-2d0)
      IF(x0 .GT. xmax) x0=xmax
      x0= max(x0, 0d0)
      q1= 1d0/alf *x0**alf  *(1d0-x0)**(bet-1d0)
      q2= x0**(alf-1d0)     *1d0/bet *( (1d0-x0)**bet -(1d0-xmax)**bet )
      p1= q1/(q1+q2)
      IF( r .LE. p1 ) THEN
         x     = x0*(r/p1)**(1d0/alf)
         dist  = x**(alf-1d0)  *(1d0-x0)**(bet-1d0)
      ELSE
         r1    = (1d0-r)/(1d0-p1)
         x     = (1d0-xmax)**bet + ((1d0-x0)**bet-(1d0-xmax)**bet)*r1
         x     = 1d0 - x**(1d0/bet)
         dist  = x0**(alf-1d0) *(1d0-x)**(bet-1d0)
      ENDIF
      djac=(q1+q2)/dist
      RETURN
  900 WRITE(*,*) ' ========= STOP in BornV_ReBin1: wrong params'
      STOP
      END


      SUBROUTINE BornV_BinBack1(xmax,x,alf,bet,RR)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//  This mapps variable x back into r (inverse of BonrV_Rebin1).            //
*//  Use : if some peaks in the distribution of x are a priori known, they   //
*//  are here converted into r so that the presampler can be informed where  //
*//  they lie (presampling is done in terms of r).                           //
*//  I dnd't bother computing the jacobian (not needed)...                   //
*//                                                                          //
*//  Maarten Boonekamp sept 2001                                             //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      DOUBLE PRECISION RR,R,alf,bet,xmax,x,djac
      DOUBLE PRECISION x0,dist,r1,p1,q1,q2
*------------------------------------------------------------------------------
      IF( alf .LE. 0d0 ) GOTO 900
      x0=(alf-1d0)/(alf+bet-2d0)
      IF(x0 .GT. xmax) x0=xmax
      x0= max(x0, 0d0)
      q1= 1d0/alf *x0**alf  *(1d0-x0)**(bet-1d0)
      q2= x0**(alf-1d0)     *1d0/bet *( (1d0-x0)**bet -(1d0-xmax)**bet )
      p1= q1/(q1+q2)
      IF( x .LE. x0 ) THEN
        RR = p1*(x/x0)**alf
      ELSE
        RR =  ((1d0- x)**bet - (1d0-xmax)**bet)
     &      / ((1d0-x0)**bet - (1d0-xmax)**bet)
        RR = 1d0 - (1d0-p1)*RR
      ENDIF
      R = MAX(RR,m_fleps**alf)
      RETURN
  900 WRITE(*,*) ' ========= STOP in BornV_BinBack1: wrong params'
      STOP
      END


      SUBROUTINE BornV_BinPeaks(Npeaks,Xpeaks)
*//////////////////////////////////////////////////////////////////////////////////
*//                                                                              //
*// Converts the peaks positions according to the KK2f conventions.              //
*//                                                                              //
*// FOR USE WITH KK2F ONLY                                                       //
*//                                                                              //
*// Arguments : Npeaks    = nb of peaks                                          //
*//             Xpeaks(i) = position of peak i in units of com energy on input,  //
*//                         and converted on output.                             //
*//                                                                              //
*// CRUCIAL :                                                                    //
*// The resulting Xpeaks array should be ordered in increasing values !!!        //
*//                                                                              //
*// Maarten Boonekamp, sept. 2001                                                //
*//////////////////////////////////////////////////////////////////////////////////
C Declarations
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Xpeaks(*)
c ... KK2f common
      INCLUDE 'BornV.h'
C Conversion
      CALL BornV_MakeGami(m_CMSene,gamiCR,gami,alfi)
      DO I = 1, Npeaks
        xtemp = 1d0 - Xpeaks(I)**2
        beta = -0.5d0
        CALL BornV_BinBack1(m_vvmax,xtemp,gamiCR,beta,RR)
        Xpeaks(I) = RR
      ENDDO
C Ordering (a CERN routine)
      CALL Mathlib_FLPSOR(Xpeaks,Npeaks)
C end
      RETURN
      END



      SUBROUTINE BornV_ReBin1a(RR,alf,bet,xmax,x,djac)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//   the same as BornV_ReBin1 but pole approximation used                   //
*//                                                                          //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      DOUBLE PRECISION RR,R,alf,bet,xmax,x,djac
      DOUBLE PRECISION x0,dist,r1,p1,q1,q2
*------------------------------------------------------------------------------
      IF( alf .LE. 0d0 ) GOTO 900
      R = MAX(RR,m_fleps**alf)
      x0=(alf-1d0)/(alf+bet-2d0)
      IF(x0 .GT. xmax) x0=xmax
      x0= max(x0, 0d0)
      q1= 1d0/alf *x0**alf
      q2= 1d0/bet *( (1d0-x0)**bet -(1d0-xmax)**bet )
      p1= q1/(q1+q2)
      IF( r .LE. p1 ) THEN
         x     = x0*(r/p1)**(1d0/alf)
         dist  = x**(alf-1d0)
      ELSE
         r1    = (1d0-r)/(1d0-p1)
         x     = (1d0-xmax)**bet + ((1d0-x0)**bet-(1d0-xmax)**bet)*r1
         x     = 1d0 - x**(1d0/bet)
         dist  = (1d0-x)**(bet-1d0)
      ENDIF
      djac=(q1+q2)/dist
      RETURN
  900 WRITE(*,*) ' ========= STOP in BornV_ReBin1: wrong params'
      STOP
      END

      SUBROUTINE BornV_ReBin2(RR,alf,bet,x,xm1,djac)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//  This mapps variable r into x. xm1=1-x kept because of rounding errors   //
*//  The same as BornV_ReBin1, but xmax=1 and bet>0                          //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      DOUBLE PRECISION RR,R,alf,bet,x,xm1,djac
      DOUBLE PRECISION x0,dist,r1,p1,q1,q2
*------------------------------------------------------------------------------
      IF( alf .LE. 0d0 ) GOTO 900
      IF( bet .LE. 0d0 ) GOTO 900
      R = MAX(RR,m_fleps**alf)
      x0=(alf-1d0)/(alf+bet-2d0)
      IF( (x0 .GT. 1d0) .OR. (x0 .LT. 0d0) ) GOTO 900
      x0= max(x0, 0d0)
      q1= 1d0/alf *x0**alf  *(1d0-x0)**(bet-1d0)
      q2= x0**(alf-1d0)     *1d0/bet *(1d0-x0)**bet
      p1= q1/(q1+q2)
      IF( r .LE. p1 ) THEN
         x    = x0*(r/p1)**(1d0/alf)
         dist = x**(alf-1d0)  *(1d0-x0)**(bet-1d0)
         xm1  = 1d0-x
      ELSE
         r1   = (1d0-r)/(1d0-p1)
         r1   = MAX(r1,m_fleps**bet)
         xm1  =(1d0-x0) *r1**(1d0/bet)
         dist = x0**(alf-1d0) *xm1**(bet-1d0)
         x    = 1d0-xm1
      ENDIF
      djac=(q1+q2)/dist
      RETURN
  900 WRITE(*,*) ' ========= STOP in BornV_ReBin2: wrong params'
      STOP
      END

      SUBROUTINE BornV_ReBin2a(RR,alf,bet,x,xm1,djac)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//  This mapps variable r into x. xm1=1-x kept because of rounding errors   //
*//  The same as BornV_ReBin2, but xmax=1 and bet>0                          //
*//  and pole approximation is used for crude/simplified distribution.       //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      DOUBLE PRECISION RR,R,alf,bet,x,xm1,djac
      DOUBLE PRECISION x0,dist,r1,p1,q1,q2
*------------------------------------------------------------------------------
      IF( alf .LE. 0d0 ) GOTO 900
      IF( bet .LE. 0d0 ) GOTO 900
      R = MAX(RR,m_fleps**alf)
      x0=(alf-1d0)/(alf+bet-2d0)
      IF( (x0 .GT. 1d0) .OR. (x0 .LT. 0d0) ) GOTO 900
      x0= max(x0, 0d0)
      q1= 1d0/alf *x0**alf
      q2= 1d0/bet *(1d0-x0)**bet
      p1= q1/(q1+q2)
      IF( r .LE. p1 ) THEN
         x    = x0*(r/p1)**(1d0/alf)
         dist = x**(alf-1d0)
         xm1  = 1d0-x
      ELSE
         r1   = (1d0-r)/(1d0-p1)
         r1   = MAX(r1,m_fleps**bet)
         xm1  = (1d0-x0) *r1**(1d0/bet)
         dist = xm1**(bet-1d0)
         x    = 1d0-xm1
      ENDIF
      djac=(q1+q2)/dist
      RETURN
  900 WRITE(*,*) ' ========= STOP in BornV_ReBin2a: wrong params'
      STOP
      END

      DOUBLE PRECISION FUNCTION BornV_RhoFoamC(xarg)
*//////////////////////////////////////////////////////////////////////////////////
*//                                                                              //
*//   Integrand for FoamC in 3-dim mode for beamstrahlung                        //
*//                                                                              //
*//   Remember that BornV_Crude and BornV_MakeRho use hidden input  m_XXXene!!   //
*//   BornV_Crude is in R-units (poitnlike xsection at  sqrt(s)=m_XXXene!        //
*//                                                                              //
*//////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      DOUBLE PRECISION  xarg(10)
      DOUBLE PRECISION  R,r1,r2
      DOUBLE PRECISION  Power,Jacob,sf12
      DOUBLE PRECISION  Rho,BornV_Crude,IRC_circee
      DOUBLE PRECISION  z1, z2, XX, RhoISR, gamiCR, gami, alfi, beta, GamBig, alpha, alpha2
      DOUBLE PRECISION  Rjac0, Rjac1, Rjac2
      DOUBLE PRECISION  zbms, zisr, y1,y2, ybms,yisr, xbms,xisr
      DOUBLE PRECISION  Par(0:3)
      INTEGER Option
      INTEGER           Icont
      DATA              Icont/0/
      Icont = Icont+1
*//////////////////////////////////////////////////////////////////////////////////////
*//  gamiCR, alpha, beta are dummy parameters in variable transformations
*//  they can be varied by 25% or so and weight distribution will be exactly the same!
*//  grid will absorb their variations.
      CALL IRC_GetParamee (Par) ! dee(z) = Par(1) *z1**Par(2)  *(1-z1)**Par(3), z1=1-x1
****  IF(Icont.LE.1) WRITE(*,*) ' Par(i)= ', Par(0), Par(1),Par(2),Par(3)
      alpha =  0.40d0           ! beamsstrahl: x1**(alpha-1), alpha manualy adjusted
      alpha =  Par(3)+1d0
      beta  = -0.50d0           ! ISR crude is as (1-vv)**(-1.5)=(1-vv)**(beta-1)
*[[[
      alpha = 3.0+1.0
*]]]
*//////////////////////////////////////////////////////////////////////////////////////
      R    = xarg(1)
      r1   = xarg(2)
      r2   = xarg(3)
      Rho  = 1d0
      Option = 1                ! Option = 1, sophisticated, a litle bit better
      Option = 2                ! Option = 2, simple, good enough
      IF( Option .EQ. 1 ) THEN
*//////////////////////////////////////////////////////////////////////////////////////
*//   R --> XX,   ZZ=1-XX=(1-vv)*(1-x1)= total loss factor, ISR and beamsstrahlung
         CALL BornV_MakeGami(m_CMSene,gamiCR,gami,alfi)          ! make gamiCR at CMSene
         IF( gami .LE. 0d0 ) GOTO 800
         GamBig = gami+2d0*alpha                            ! total singularity at XX=0
         CALL BornV_ReBin1a(R,GamBig,beta,m_vvmax,XX,RJac0) ! Mapping  R => XX=1-ZZ
         Rho = Rho *RJac0
*//   r1 --> m_vv
         alpha2 = 2d0*alpha
         CALL BornV_ReBin2a(r1, gami, alpha2, yisr, ybms, RJac1) ! Mapping  r1 => m_vv
         xisr = yisr *XX
         xbms = ybms *XX/(1d0-yisr*XX)
         Rho  = Rho  *XX/(1d0-yisr*XX) *RJac1
         zisr = 1d0-xisr
         zbms = 1d0-xbms
         m_vv = xisr
*//   r2 --> m_x2
         CALL BornV_ReBin2a(r2, alpha, alpha, y1, y2, RJac2) ! Mapping  r2 => m_x2
         m_x1 =   y1*xbms
         m_x2 =   y2*xbms/(1d0-y1*xbms)
         Rho  = Rho *xbms/(1d0-y1*xbms) *RJac2 
*//////////////////////////////////////////////////////////////////////////////////////
*//   simplified analytical importance sampling transformations
      ELSEIF( Option .EQ. 2 ) THEN
         CALL BornV_MakeGami(m_CMSene,gamiCR,gami,alfi)           ! make gamiCR at CMSene
         IF( gami .LE. 0d0 ) GOTO 800
         m_vv  = R**(1d0/gami)*m_vvmax
         Rho   = Rho* m_vv/R/gami*m_vvmax
         m_x1  = r1**(1d0/alpha)                             ! Mapping  r1 => x1
         Rho   = Rho   *m_x1/r1/alpha
         m_x2  = r2**(1d0/alpha)                             ! Mapping  r2 => x2
         Rho = Rho   *m_x2/r2/alpha
***         IF( (1d0-m_vv)*(1d0-m_x1)*(1d0-m_x2) .LT. (1d0-m_vvmax) ) GOTO 800  ! WRONG
         IF( m_vv .GT. m_vvmax ) GOTO 800  ! vmax from input, the best to keep vmax=1
         IF( m_x1 .GT. 0.99999 ) GOTO 800  ! cutting off extremely perigheral
         IF( m_x2 .GT. 0.99999 ) GOTO 800  ! cutting off extremely perigheral
         IF( m_CMSene*SQRT((1d0-m_vv)*(1d0-m_x1)*(1d0-m_x2)) .LT. 1d0 ) GOTO 800  ! mass(Z)>1GeV
      ENDIF
      z1 = 1d0-m_x1
      z2 = 1d0-m_x2
*//////////////////////////////////////////////////////////////////////////////////////
*//   Calculate ISR crude structure function (the same as in Karlud)
      m_XXXene =  m_CMSene*SQRT(z1*z2)                   ! hidden input for BornV_Crude
      CALL BornV_MakeISR(RhoISR)                         !<-- uses m_XXXene and m_vv
      Rho = Rho *RhoISR
*//////////////////////////////////////////////////////////////////////////////////////
*//   Beamsstrahlung structure function, singular as m_x1**(alpha-1)
      IF( (z1.EQ.1d0) .OR. (z2.EQ.1d0) ) THEN ! rounding errors may cause problems
         SF12 = 0d0
      ELSE
***      SF12 = IRC_circee( z1, z2 )
         SF12 = Par(1) *m_x1**Par(3) *z1**Par(2)   *Par(1) *m_x2**Par(3) *z2**Par(2)
***      SF12 = Par(1) *m_x1**Par(3)               *Par(1) *m_x2**Par(3)    ! Truncated
*[[[
* Valence 2*x*u(x):   XUPV = 2.18 * X**0.5D0 * (1.D0-X)**3.D0
* Valence x*d(x):     XDNV = 1.23 * X**0.5D0 * (1.D0-X)**4.D0
* Sea     x*s(x):    XSEA = 0.6733 * X**(-0.2D0) * (1.D0-X)**7.D0 
        SF12 = 2.18 *m_x1**3.D0 *z1**0.5D0   *0.6733 *m_x2**7.D0 *z2**(-0.2D0)/z1/z2
*]]]
      ENDIF
      Rho = Rho *SF12

* Born Xsection at s' = m_XXXene**2 *(1-vv)
      BornV_RhoFoamC = Rho*BornV_Crude(m_vv)/(1d0-m_vv)
      RETURN
 800  CONTINUE
      BornV_RhoFoamC =0d0
      RETURN
 900  CONTINUE
      WRITE(*,*) ' STOP in BornV_RhoFoamC, m_x1 = ', m_x1
      WRITE(*,*) ' XX, m_vv= ', XX, m_vv
      STOP
      END

      DOUBLE PRECISION FUNCTION BornV_RhoFoamB(xarg)
*//////////////////////////////////////////////////////////////////////////////////
*//                                                                              //
*//   Integrand for FoamB in 2-dim mode for beamstrahlung                        //
*//                                                                              //
*//   Remember that BornV_Crude and BornV_MakeRho use hidden input  m_XXXene!!   //
*//   BornV_Crude is in R-units (poitnlike xsection at  sqrt(s)=m_XXXene!        //
*//                                                                              //
*//                                                                              //
*//                                                                              //
*//////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      DOUBLE PRECISION  xarg(10)
      DOUBLE PRECISION  R,r1
      DOUBLE PRECISION  Rho,BornV_Crude,IRC_circee
      DOUBLE PRECISION  gamiCR, gami, alfi, beta, RJac0, RJac1, RJacS
      DOUBLE PRECISION  alpha,  alpha2,  eps,  eps2
      DOUBLE PRECISION  GamBig, RhoISR, SF1, XX, yisr, ybms, z1, aa
      DOUBLE PRECISION  anor, xnor
      DOUBLE PRECISION  Par(0:3)
      INTEGER Option
      INTEGER           Icont
      DATA              Icont/0/
      Icont = Icont+1
*//////////////////////////////////////////////////////////////////////////////////////
*//  gamiCR, alpha, beta are dummy parameters in variable transformations
*//  they can be varied by 25% or so and weight distribution will be exactly the same!
*//  grid will absorb their variations.
      CALL IRC_GetParamee (Par) ! dee(z) = Par(1) *z1**Par(2)  *(1-z1)**Par(3), z1=1-x1
****  IF(Icont.LE.1) WRITE(*,*) ' Par(i)= ', Par(0), Par(1),Par(2),Par(3)
      alpha =  0.40d0           ! beamsstrahl: x1**(alpha-1), alpha manualy adjusted
      alpha =  Par(3)+1d0
      beta  = -0.50d0           ! ISR crude is as (1-vv)**(-1.5)=(1-vv)**(beta-1)
*//////////////////////////////////////////////////////////////////////////////////////
      R    = xarg(1)
      r1   = xarg(2)
      m_x2 = 0d0
      Rho  = 1d0
      Option = 3                ! Option = 3 for tests of normalization only
      Option = 1                ! Option = 1 sophisticated
      Option = 2                ! Option = 2 simple, good enough
      IF( Option .EQ. 1 ) THEN
*//////////////////////////////////////////////////////////////////////////////////////
*//   (over)complicated analytical importance sampling transformations
*//   R --> XX,   ZZ=1-XX=(1-vv)*(1-x1)= total loss factor, ISR and beamsstrahlung
         CALL BornV_MakeGami(m_CMSene,gamiCR,gami,alfi)          ! make gamiCR at CMSene
         IF( gami .LE. 0d0 ) GOTO 800
         GamBig = gami+alpha                                ! total singularity at XX=0
         CALL BornV_ReBin1a(R,GamBig,beta,m_vvmax,XX,RJac0) ! Mapping  R => XX
         Rho = Rho *RJac0
*//   r1 --> m_vv
         IF( gami .LE. 0d0 ) GOTO 800
         CALL BornV_ReBin2a(r1, gami, alpha, yisr, ybms, RJac1) ! Mapping  r1 => m_vv
         m_vv =  yisr* XX
         m_x1 =  ybms* XX/(1d0-yisr*XX)
         Rho  =  Rho  *XX/(1d0-yisr*XX) *RJac1
         IF( m_x1 .GE. 1d0) GOTO 900
*//////////////////////////////////////////////////////////////////////////////////////
*//   simplified analytical importance sampling transformations
      ELSEIF( Option .EQ. 2 ) THEN
         CALL BornV_MakeGami(m_CMSene,gamiCR,gami,alfi)           ! make gamiCR at CMSene
         IF( gami .LE. 0d0 ) GOTO 800
         m_vv  = R**(1d0/gami)*m_vvmax
         Rho   = Rho* m_vv/R/gami*m_vvmax
         m_x1  = r1**(1d0/alpha)                             ! Mapping  r1 => x1
         Rho   = Rho *m_x1/r1/alpha
         IF( (1d0-m_vv)*(1d0-m_x1) .LT. (1d0-m_vvmax) ) GOTO 800
      ELSEIF( Option .EQ. 3 ) THEN
*//   primitive test options, usefull for checking normalization
         m_vv = R*m_vvmax
         Rho = Rho *m_vvmax
         m_x1  = r1
         IF( (1d0-m_vv)*(1d0-m_x1) .LT. (1d0-m_vvmax) ) GOTO 800
      ENDIF
*//////////////////////////////////////////////////////////////////////////////////////
*//   Calculate ISR crude structure function (the same as in Karlud)
      m_XXXene =  m_CMSene*SQRT(1d0-m_x1)              ! hidden input for BornV_Crude
      CALL BornV_MakeISR(RhoISR)                       !<-- uses m_XXXene and m_vv
      Rho = Rho *RhoISR
*//////////////////////////////////////////////////////////////////////////////////////
*//   Beamsstrahlung structure function, singular as m_x1**(alpha-1)
      z1 = 1d0-m_x1
      IF( (z1.EQ.1d0) .OR. (m_x1.EQ.0d0) ) THEN ! rounding errors may cause problems
         SF1 = 0d0
      ELSE
*****    SF1 = 2d0 *IRC_circee( z1, 1d0 )   ! factor 2 due to implicit symmetrization x1<-->x2
*****    SF1 = 2d0 *Par(0) *Par(1) *m_x1**Par(3)               ! truncated
         SF1 = 2d0 *Par(0) *Par(1) *m_x1**Par(3) *z1**Par(2)   ! the same as circee
      ENDIF
      Rho = Rho *SF1
*[[[
      Rho = 1d-20
*]]]
*//////////////////////////////////////////////////////////////////////////////////////
*//   Born Xsection at s' =m_XXXene**2 *(1-vv) =m_CMSene**2 *(1-XX)
      BornV_RhoFoamB = Rho* BornV_Crude(m_vv)/(1d0-m_vv)
      RETURN
 800  CONTINUE
      BornV_RhoFoamB =0d0
      RETURN
 900  CONTINUE
      WRITE(*,*) ' STOP in BornV_RhoFoamB, m_x1 = ', m_x1
      WRITE(*,*) ' XX, m_vv= ', XX, m_vv
      STOP
      END                       ! BornV_RhoFoamB


      DOUBLE PRECISION FUNCTION BornV_RhoFoamA(xarg)
*//////////////////////////////////////////////////////////////////////////////////
*//                                                                              //
*//   Integrand for FoamA in 1-dim mode beamstrahlung off/on                     //
*//   !!! DEFINES m_vv !!!!                                                      //
*//                                                                              //
*//   Remember that BornV_Crude and BornV_MakeRho use hidden input  m_XXXene!!   //
*//   BornV_Crude is in R-units (poitnlike xsection at  sqrt(s)=m_XXXene!        //
*//                                                                              //
*//////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      DOUBLE PRECISION  xarg(10)
      DOUBLE PRECISION  R
      DOUBLE PRECISION  Rho,BornV_Crude
      DOUBLE PRECISION  gamiCR, gami, alfi, beta, RJac
      DOUBLE PRECISION  IRC_circee
*----------------------------------------
      R  = xarg(1)
      m_x1 = 0d0
      m_x2 = 0d0
*-----------------------------------------------
      m_XXXene =  m_CMSene        ! hidden input for BornV_Crude
      CALL BornV_MakeGami(m_XXXene,gamiCR,gami,alfi)
      IF( gami .LE. 0d0 ) GOTO 800
*     Mapping  r => vv change  to improve on efficiency
      m_vv  = R**(1d0/gami)*m_vvmax
      RJac  = m_vv/R/gami*m_vvmax
      CALL BornV_MakeISR(Rho)                  !<-- uses m_XXXene and m_vv
      Rho = Rho*RJac
*----------------------------------------
      Rho = Rho *IRC_circee(1d0,1d0)           !<-- implicit factor from circee 
*[[[
      Rho = 1d-20
*]]]
*----------------------------------------
* Born Xsection at s' = m_XXXene**2 *(1-vv)
      IF(m_KeyZet .EQ. -2) THEN   ! Artificial constant x-section for test runs
         BornV_RhoFoamA = Rho* BornV_Crude(0d0)
      ELSE                        ! 1/(1-vv) because BornV_Crude is in R-units
         BornV_RhoFoamA = Rho* BornV_Crude(m_vv)/(1d0-m_vv)
      ENDIF
      RETURN
 800  CONTINUE
      BornV_RhoFoamA =0d0
      END

      DOUBLE PRECISION  FUNCTION BornV_RhoVesko1(R)
*//////////////////////////////////////////////////////////////////////////////////
*//                                                                              //
*//   Integrand of Vesko, dSigma/dV(V)*Jacobians function of R                   //
*//   !!!! DEFINES m_vv !!!!                                                     //
*//                                                                              //
*//   Remember that BornV_Crude and BornV_MakeRho use hidden input  m_XXXene!!   //
*//   BornV_Crude is in R-units (poitnlike xsection at  sqrt(s)=m_XXXene!        //
*//                                                                              //
*//   In  the case of beamsstrahlung additional normalization                    //
*//   factor circee(1d0,1d0) is added in BStra_Initialize                        //
*//////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      DOUBLE PRECISION   R,Rho, BornV_Crude
      DOUBLE PRECISION   gamiCR, gami, alfi, beta, RJac
*-----------------------------------------------
      m_XXXene =  m_CMSene                 ! hidden input for BornV_Crude
*-----------------------------------------------
      CALL BornV_MakeGami(m_XXXene,gamiCR,gami,alfi)
      IF( gami .LE. 0d0 ) GOTO 800
*     Mapping  r => vv change  to improve on efficiency
*     IMPORTANT: alpha=gami MUST coincide with vv->0 limit in rho(vv)!!!
*     If not, then problems with m_fleps for weak bremsstrahlung! 
      beta = -0.5d0
      CALL BornV_ReBin1(R,gami,beta,m_vvmax,m_vv,RJac)
      CALL BornV_MakeISR(Rho)              ! uses m_XXXene and m_vv
      Rho = Rho*RJac
*-----------------------------------------------
* Translate R into m_vv and get QED (crude) density Rho
***** CALL BornV_MakeRho(R,Rho)
*-----------------------------------------------
* Born Xsection at s' = m_XXXene**2 *(1-vv)
      IF(m_KeyZet .EQ. -2) THEN   ! Artificial constant x-section for test runs
         BornV_RhoVesko1 = Rho* BornV_Crude(0d0)
      ELSE                        ! 1/(1-vv) because BornV_Crude is in R-units
         BornV_RhoVesko1 = Rho* BornV_Crude(m_vv)/(1d0-m_vv)
      ENDIF
      RETURN
 800  CONTINUE
      BornV_RhoVesko1 =0d0
      END

      SUBROUTINE BornV_MakeGami(CMSene,gamiCR,gami,alfi)
*//////////////////////////////////////////////////////////////////////////////
*//   Crude Gami as a function of CMSene                                     //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      DOUBLE PRECISION  CMSene, gamiCR, gami, alfi
      DOUBLE PRECISION  amel, svar, am2, beta, chini2
      INTEGER           KFbeam
*---------------------------------
      KFbeam = m_KFini      ! from the input
      amel   = m_amferm(KFbeam)
      am2  = (2d0*amel/CMSene)**2
      chini2 = m_Qf(ABS(KFbeam))**2  ! electriccharge of beam fermion
      alfi   = chini2 *m_alfpi
      IF( am2 .GT. 1d0 ) GOTO 800
      beta = SQRT(1d0-am2)
      gami    = 2d0*chini2*m_alfpi *( DLOG((1+beta)**2/am2) -1d0)
      gamiCR  = 2d0*chini2*m_alfpi *  DLOG((1+beta)**2/am2)
      gamiCR  = gamiCR *m_Xenph         !!! enhancement of crude photon multiplicity
      IF(m_KeyWtm .EQ. 1) gamiCR=gami   !!! new, for very special tests
*-------------
      RETURN
 800  CONTINUE
      gamiCR = 0d0
      gami   = 0d0
      END

      SUBROUTINE BornV_MakeISR(Rho)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//   This procedure is tightly related to ISR photon generation in Karlud   //
*//   It calculates Rho(m_vv, m_XXXene) QED crude Structure Function         //
*//                                                                          //
*//   m_AvMult is later used in KarLud_YFSini                                //
*//   m_YFSkon ,m_YFS_IR are later used in GPS_Make  and QED3_Make           //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      DOUBLE PRECISION   Rho
      DOUBLE PRECISION   gami,  gamiCR,  alfi, BornV_Crude
      DOUBLE PRECISION   xBorn, DilJac0, beta, VoluMC
*-----------------------------
      CALL BornV_MakeGami(m_XXXene,gamiCR,gami, alfi)
      IF(m_vv .GT. m_vvmin) THEN
         DilJac0   = (1d0+1d0/SQRT(1d0-m_vv))/2d0
         m_AvMult  = gamiCR*DLOG(m_vv/m_vvmin)
         VoluMC    = gamiCR/m_vv *EXP( m_AvMult )    !!! Phase space Volume CRUDE
         m_YFS_IR  = -gami*DLOG(1d0/m_vvmin)         !!! IR part of YFS formfactor
         Rho       = VoluMC *EXP(m_YFS_IR)
      ELSE
         DilJac0   = 1d0
         m_AvMult  = 0d0
         VoluMC    = 1d0
* IMPORTANT:     The integral over Rho(v<vvmin) = YFS_IR = EXP(-gami*LOG(1/vvmin))
         m_YFS_IR  = -gami*DLOG(1d0/m_vvmin)         !!! IR part of YFS formfactor
         Rho       = 1d0/m_vv *gami*m_vv**gami
      ENDIF
      Rho =  Rho * DilJac0
* YFS formfactor, finite part, YFS_form_Factor = EXP(YFS_IR + YFSkon)
* YFSkon is delegated/exported to QED3 and GPS (not used here).
      m_YFSkon =  EXP(1/4d0 *gami + alfi*( -.5d0  +m_pi**2/3d0) )
      m_YFS_IR =  EXP(m_YFS_IR)
      END

      SUBROUTINE BornV_MakeRho(R,Rho)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//   This function is tightly related to ISR photon generation in Karlud    //
*//   It calculates Rho(vv) QED crude Structure Function at XXXene           //
*//   Translates R into m_vv                                                 //
*//                                                                          //
*//   m_vv is used in Karlud and other places                                //
*//   m_AvMult is later used in KarLud_YFSini                                //
*//   m_YFSkon ,m_YFS_IR are later used in GPS_Make  and QED3_Make           //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      DOUBLE PRECISION   R,Rho
      DOUBLE PRECISION   gami,  gamiCR,  alfi, BornV_Crude
      DOUBLE PRECISION   xBorn, DilJac0, beta, RJac,  VoluMC
*-----------------------------
      CALL BornV_MakeGami(m_XXXene,gamiCR,gami,alfi)
      IF( gamiCR .LE. 0d0 ) GOTO 800
*     Mapping  r => vv change  to improve on efficiency
      beta = -0.5d0
      CALL BornV_ReBin1(R,gami,beta,m_vvmax,m_vv,RJac)
      IF(m_vv .GT. m_vvmin) THEN
         DilJac0   = (1d0+1d0/SQRT(1d0-m_vv))/2d0
         m_AvMult  = gamiCR*DLOG(m_vv/m_vvmin)
         VoluMC    = gamiCR/m_vv *EXP( m_AvMult )    !!! Phase space Volume CRUDE
         m_YFS_IR  = -gami*DLOG(1d0/m_vvmin)         !!! IR part of YFS formfactor
         Rho       = VoluMC *EXP(m_YFS_IR)
      ELSE
         DilJac0   = 1d0
         m_AvMult  = 0d0
         VoluMC    = 1d0
* IMPORTANT:     The integral over Rho(v<vvmin) = YFS_IR = EXP(-gami*LOG(1/vvmin))
         m_YFS_IR  = -gami*DLOG(1d0/m_vvmin)         !!! IR part of YFS formfactor
         Rho       = 1d0/m_vv *gami*m_vv**gami
      ENDIF
      Rho =  Rho * DilJac0*RJac
* YFS formfactor, finite part, YFS_form_Factor = EXP(YFS_IR + YFSkon)
* YFSkon is delegated/exported to QED3 and GPS (not used here).
      m_YFSkon =  EXP(1/4d0 *gami + alfi*( -.5d0  +m_pi**2/3d0) )
      m_YFS_IR =  EXP(m_YFS_IR)
      RETURN
 800  CONTINUE
      Rho  = 0d0
      m_vv = 0d0
      END


      DOUBLE PRECISION  FUNCTION BornV_Crude(vv)
*///////////////////////////////////////////////////////////////////////////////
*//                                                                           //
*// This routine calculates CRUDE total Born cross section  SUMMED OVER KF.   //
*// It exploits the fact that born x. section = a + b*c + d*c**2              //
*// Hidden input is m_XXXene                                                  //
*//                                                                           //
*///////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
*
      INTEGER KFf
      DOUBLE PRECISION   vv, svar1
      DOUBLE PRECISION   BornV_Differential
      DOUBLE PRECISION   Born, sum
* for brancher
      DOUBLE PRECISION   WMList(200),XsList(200)
      INTEGER            j, KFlist(200), Nbranch, KeyRes
      DOUBLE PRECISION   RRes_CorQQ, xmf, BornV_GetMass
*-----------------------------------------------------------------------
      svar1  = (1-vv)*m_XXXene**2
* get from brancher list of KF's and of enhancement factors
      CALL MBrA_GetKFlist(Nbranch,KFlist)
      CALL MBrA_GetWMList(Nbranch,WMList)
      sum = 0d0
      DO j=1,Nbranch
         Born =0d0
         KFf=KFlist(j)
*///////////////////////////////////////////////////////////////////////////////
         Born= BornV_Differential( 0, Kff, svar1, 0.d0, 0.d0,0.d0, 0.d0,0.d0 ) !
* Below is the modification of hadronic cross sercion according to experiemntal R
         CALL BornV_GetKeyRes(KeyRes)
         IF( KeyRes.NE.0 .AND. KFf.LE.5 )  THEN 
            xmf = BornV_GetMass(KFf)
            Born= Born * RRes_CorQQ(DSQRT(svar1),KFf,xmf)
         ENDIF
*///////////////////////////////////////////////////////////////////////////////
* For light quarks u,d,s, special cut on mass (just in case)
         IF( (ABS(KFf) .GE. 1) .AND. (ABS(KFf) .LE. 3)) THEN
            IF( svar1 .LE. m_HadMin**2) Born=0d0
         ENDIF
* The amplification factor WM goes into crude normalization
* It is countered later on by the weight from MBrA_GenKF
         sum = sum +Born*WMList(j)
         XsList(j) = Born              !<---  WtMax=WMList(j) NOT included!!!
      ENDDO
* send back to bracher xsections for generation of KF
      CALL MBrA_SetXSList(XsList)
*------
      BornV_Crude =sum                 !<---  WtMax=WMList(j) IS included!!!
      END

      DOUBLE PRECISION  FUNCTION BornV_Differential(Mode,KFf,svar,CosThe,eps1,eps2,ta,tb)
*///////////////////////////////////////////////////////////////////////////////
*// Mode=0 it is CRUDE version for pure Born, no spin, no EW corrs.           //
*//                                                                           //
*// Mode=1 full result with electroweak corrs. spin etc. added.               //
*//        used in QED3 and all kind of tests                                 //
*//                                                                           //
*// Mode=3 for tests of pretabulation, GSW(s,theta) has to be provided from   //
*//        outside with help of BornV_SetGSW                                  //
*//                                                                           //
*// Note that in the test mode KeyEwl=0 and Mode=1 we use BornV_Simple        //
*// which perhaps will have to be changed in future besause lack of spin eff. //
*// At this stage however we are bound to use it because the KeyZet etc.      //
*// are implemented only in BornV_Simple and not in BornV_Dizet.              //
*//                                                                           //
*///////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE

      INTEGER Mode
      DOUBLE PRECISION   svar,CosThe,eps1,eps2,ta,tb
      INTEGER KFf
      INCLUDE 'BornV.h'
      SAVE
      DOUBLE PRECISION   Born
      DOUBLE PRECISION   BornV_Dizet, BornV_Simple
*-----------------------------------------------------------------------------
      IF(     Mode .EQ. 1 ) THEN
         IF( m_KeyElw .EQ. 0 ) THEN
            Born= BornV_Simple( m_KFini,KFf,svar,CosThe)
         ELSE
*           Linear interpolation from tables, only for Mode=1 
            CALL BornV_InterpoGSW( ABS(KFf),  svar, CosThe)
            Born= BornV_Dizet( 1,m_KFini,KFf, svar, CosThe, eps1,eps2,ta,tb)
         ENDIF
      ELSEIF( Mode .EQ. 3 ) THEN
*           For test of pretabulation, BornV_SetGSW has to be invoked prior
            Born= BornV_Dizet( 1,m_KFini,KFf, svar, CosThe, eps1,eps2,ta,tb)
      ELSEIF( Mode .EQ. 0 ) THEN
         Born= BornV_Simple( m_KFini,KFf,svar,CosThe)
*        Another potential possibility, with a different threshold behavior is:
*        Born= BornV_Dizet( 0,m_KFini,KFf,svar,CosThe,0d0,0d0,0d0,0d0)
      ELSE
         WRITE(*,*) 'STOP in BornV_Differential: Mode =',Mode
         STOP
      ENDIF
      BornV_Differential = Born
      END


      DOUBLE PRECISION  FUNCTION BornV_Simple(KFi,KFf,svar,costhe)
*///////////////////////////////////////////////////////////////////////////////
*//                                                                           //
*// This routine provides unsophisticated Born differential cross section     //
*// at the crude x-section level, with Z and gamma s-chanel exchange.         //
*// Note that it uses m_swsq from tables (Dizet) and not from xpar(503)       //
*//                                                                           //
*///////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
*
      INTEGER KFi,KFf
      DOUBLE PRECISION   svar,costhe
*
      DOUBLE PRECISION   ss,T3e,Qe,deno,Ve,Ae
      DOUBLE PRECISION   ye,yf,xf,rechi,xe,amx2
      DOUBLE PRECISION   thresh,ff0,ff1,chi2
      DOUBLE PRECISION   t3f,amfin,af,vf,born,sum,qf
      DOUBLE PRECISION   BWD,Coef
      DOUBLE COMPLEX     PropW,WVPi
      INTEGER NCF
*----------------------
      ss = svar
* Z and gamma couplings to beams (electrons)
      T3e = m_T3f(KFi)  ! isospin, L-hand component
      Qe  = m_Qf( KFi)  ! electric charge

      deno= 4d0*sqrt(m_swsq*(1d0-m_swsq))
      Ve  = (2*T3e -4*Qe*m_swsq)/deno
      Ae  =  2*T3e              /deno

      NCf   = m_NCf(KFf)        ! number of colours
      T3f   = m_T3f(KFf)        ! isospin, L-hand component
      Qf    = m_Qf( KFf)        ! electric charge
      deno  = 4d0*sqrt(m_swsq*(1d0-m_swsq))
      Vf    = (2*T3f -4*Qf*m_swsq)/deno
      Af    =  2*T3f              /deno
* Switch off Z or gamma
      IF(m_KeyZet .EQ. 0) THEN
         Ve=0d0
         Ae=0d0
      ENDIF
      IF(m_KeyZet .EQ. 9) THEN
         Qe=0d0
         Qf=0d0
      ENDIF
c[[   BWD = (ss-m_MZ**2)**2 + (m_gammz*m_MZ)**2   !!! <--! fixed width
      BWD = (ss-m_MZ**2)**2 + (m_gammz*ss/m_MZ)**2
      chi2 = ss**2          /BWD
      rechi=(ss-m_MZ**2)*ss /BWD
      xe= Ve**2 +Ae**2
      xf= Vf**2 +Af**2
      ye= 2*Ve*Ae
      yf= 2*Vf*Af
      ff0= qe**2*qf**2 +2*rechi*qe*qf*Ve*Vf +chi2*xe*xf
      ff1=             +2*rechi*qe*qf*Ae*Af +chi2*ye*yf
      Born    = (1d0+ costhe**2)*ff0 +2d0*costhe*ff1
      IF(ABS(KFf).EQ.12) THEN
* Electron neutrino, approximate, actually matters costhe=0 for primary spectrum ....
        Coef  =1.d0/2.d0/m_swsq
        Born  = Born+0.25*(0.75D0*(1d0-costhe)**2+0.75*(1d0+costhe)**2)
     $                   *Coef**2*ss**2/(m_MZ**2*(1-m_swsq))**2
      ENDIF
*     Colour factor
      Born = NCf*Born
      IF( ABS(costhe) .GT. 1d0) WRITE(*,*) '----> BornV: costhe=',costhe
* This is a bit crude method of introducing threshold behaviour
* cos(theta) depencence incorrect!!!
      amfin = m_amferm(KFf)     ! mass
      IF(    svar .LE.  4d0*amfin**2) THEN
         thresh=0d0
      ELSEIF(svar .LE. 160d0*amfin**2) THEN
         amx2=4d0*amfin**2/svar
         thresh=sqrt(1d0-amx2)*(1d0+amx2/2d0)
      ELSE
         thresh=1d0
      ENDIF
      Born= Born*thresh
      BornV_Simple = Born
      END

      DOUBLE PRECISION  FUNCTION BornV_Integrated(KFfin,svar)
*///////////////////////////////////////////////////////////////////////////////
*//                                                                           //
*//   !!!!!!!!!!!!!! USED ONLY IN SEMIANALYTICAL programs!!!!!!!!!!           //
*//                                                                           //
*// This routine calculates total Born cross section.                         //
*// It is NOT used in MC any more                                             //
*//                                                                           //
*// It exploits the fact that born x. section = a + b*c + d*c**2              //
*//                                                                           //
*// For KFfin = 0 we sum over all alowed flavours, otherwise,                 //
*// for KFfin.NE.0 we calculate xsect for the actual value of m_KFfin         //
*//                                                                           //
*///////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
*
      INTEGER KFfin
      DOUBLE PRECISION   svar
      DOUBLE PRECISION   BornV_Differential
      DOUBLE PRECISION   Born,Sum
      INTEGER KFf
*-----------------------------------------------------------------------
* Selective Inclusive/Exclusive Loop over all final fermions
      Sum = 0d0
      DO KFf=1,20
         Born =0d0
         IF( m_IsGenerated(KFf) .NE.  0) THEN
            IF((KFfin .EQ. 0  )  .OR. ! Inclusive
     $         (KFfin .EQ. KFf)) THEN ! Exclusive
               Born= BornV_Differential( 0,Kff,svar, 0.d0, 0.d0,0.d0, 0.d0,0.d0 )
            ENDIF
         ENDIF
         Sum = Sum +Born
      ENDDO
      BornV_Integrated =Sum
      END



      DOUBLE PRECISION FUNCTION BornV_Dizet(Mode,KFi,KFf,svar,CosThe,eps1,eps2,ta,tb)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//   Calculates differential born cross section.                            //
*//   For Mode=0 pure Born and for Mode=1 electroweak corrs. are added.      //
*//   KFi,KFf can be also negative for antiparticle, in this case it is      //
*//   important to produce tables with correct input KFini, KFfin !!!        //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      INTEGER          Mode,KFi,KFf
      DOUBLE PRECISION svar,CosThe,eps1,eps2,ta,tb
*-------------------------------------------------------------------------------------
      DOUBLE PRECISION  pi
      PARAMETER( pi =3.141592653589793238462643d0 )
      INTEGER    nneut
      PARAMETER( nneut = 1)
*-------------------------------------------------------------------------------------
      DOUBLE PRECISION 
     $  xupgi(2),                   ! Left/Right coupling gamma initial
     $  xupzi(2),                   ! Left/Right coupling Z     initial
     $  xupgf(2),                   ! Left/Right coupling gamma final
     $  xupzf(2)                    ! Left/Right coupling Z     final
*-------------------------------------------------------------------------------------
      DOUBLE PRECISION 
     $  t3e,                        ! Left izospin initial
     $  qe,                         ! Charge       initial
     $  t3f,                        ! Left izospin final
     $  qf                          ! Charge       final
      INTEGER
     $  kolor                       ! Color final fermion
*-------------------------------------------------------------------------------------
      DOUBLE COMPLEX   aborn(2,2),aphot(2,2),azett(2,2)
      DOUBLE COMPLEX   xupzfp(2),xupzip(2)
      DOUBLE COMPLEX   abornm(2,2),aphotm(2,2),azettm(2,2)
      DOUBLE COMPLEX   propa,propz
      DOUBLE COMPLEX   xr,xi,propw,aw(2,2)
      DOUBLE COMPLEX   xupf,xupi,xff(4),xfem,xfota,xrho,xke,xkf,xkef
      DOUBLE COMPLEX   xthing,xve,xvf,xvef,DelW
*
      INTEGER          j,i,ivini,kdumm,kff0,mode0,ivfin
      DOUBLE PRECISION xm2,xp2,xmw,regulm,regula,del1,xef,del0,factom,factor,thresh
      DOUBLE PRECISION xm3,helit,polar2,polar1,helic,born,amin,aizor,xgw,cost0,svar0
      DOUBLE PRECISION qem,qfm,xf,xe,beta,amfin,aizol,sinthe,xcoup
      DOUBLE PRECISION RSQV,RSQA
      DOUBLE PRECISION t,s
*-------------------------------------------------------------------------------------
      INTEGER    icont
      DATA       icont /0/
*/////////////////////////////////////////////////////////////////////////////
* Translation table KF-->IV
      INTEGER IV(-16:16)
      DATA IV / -1, -2, -1, -2, -1, -2, 4*0, -3, -4, -3, -4, -3, -4,  0,  
     $           4,  3,  4,  3,  4,  3, 4*0,  2,  1,  2,  1,  2,  1    /
*/////////////////////////////////////////////////////////////////////////////
      DATA xi/(0.d0,1.d0)/,xr/(1.d0,0.d0)/
      DATA xgw/2.5d0/
* To be sure that initialization starts properly
      DATA Mode0,svar0,cost0,KFf0 /-155,-155.d0,-156.d0,-99/
*-------------------------------------------------------------------------------------

      icont=icont+1
*////////////////////////////////////////////////////////////////////////
*//    Save CPU for the same svar, CosThe  and varying spins ta, tb    //
*////////////////////////////////////////////////////////////////////////
      IF (Mode.NE.Mode0 .OR. svar.NE.svar0 .OR. CosThe.NE.cost0 
     $                  .OR. KFf.NE.KFf0  ) THEN
         Mode0  = Mode
         svar0  = svar
         cost0  = CosThe
         KFf0   = KFf
*////////////////////////////////////////////////////////////////////////
*//               Coupling constants                                   //
*////////////////////////////////////////////////////////////////////////
c[[[[[[[[[[[[[[[[[[[[[[!!!!!!!!!!!!!!!!!!!!!!
c      icont=icont+1
c      IF(icont.LE.200) THEN
c         write(*,*) '|||||||||||||||||||||BornV|||||||||||||||||||||||||||||||||||||'
c         write(*,'(a,8g22.14)') 'sqrt(svar),costhe=',sqrt(svar),costhe
c         write(*,'(a,8g22.14)') 'QCDcor=',m_QCDcor
c      ENDIF
c]]]]]]]]]]]]]]]]]]]]]]!!!!!!!!!!!!!!!!!!!!!!
         amin  = m_amferm(ABS(KFi))
         IVini = IV(KFi)
         CALL BornV_givizo( IVini, 1,aizor,qe,kdumm)
         CALL BornV_givizo( IVini,-1,aizol,qe,kdumm)
         xupgi(1)=qe
         xupgi(2)=qe
         t3e    = aizol+aizor
         xupzi(1)=(aizor-qe*m_swsq)/sqrt(m_swsq*(1-m_swsq))
         xupzi(2)=(aizol-qe*m_swsq)/sqrt(m_swsq*(1-m_swsq))
*
         amfin = m_amferm(ABS( KFf ))
         IVfin = IV( KFf )
         CALL BornV_givizo( IVfin, 1,aizor,qf,kolor)
         CALL BornV_givizo( IVfin,-1,aizol,qf,kolor)
         xupgf(1)=qf
         xupgf(2)=qf
         t3f    =  aizol+aizor
         xupzf(1)=(aizor -qf*m_swsq)/sqrt(m_swsq*(1-m_swsq))
         xupzf(2)=(aizol -qf*m_swsq)/sqrt(m_swsq*(1-m_swsq))
*
         sinthe = sqrt(1.d0-CosThe**2)
         beta   = SQRT(MAX(0d0,1d0-4d0*amfin**2/svar))
c[[[[[!!!!!!!!!!!!!!!!!!!!!!!
ccc         beta=1d0
c]]]]]!!!!!!!!!!!!!!!!!!!!!!!

* Multiply axial coupling by beta factor.
         xupzfp(1)= 0.5d0*(xupzf(1)+xupzf(2))+0.5*beta*(xupzf(1)-xupzf(2))
         xupzfp(2)= 0.5d0*(xupzf(1)+xupzf(2))-0.5*beta*(xupzf(1)-xupzf(2))
         xupzip(1)= 0.5d0*(xupzi(1)+xupzi(2))     +0.5*(xupzi(1)-xupzi(2))
         xupzip(2)= 0.5d0*(xupzi(1)+xupzi(2))     -0.5*(xupzi(1)-xupzi(2))
* Final state vector coupling
         xupf     = 0.5d0*(xupzf(1)+xupzf(2))
         xupi     = 0.5d0*(xupzi(1)+xupzi(2))
         xthing   = 0d0
*////////////////////////////////////////////////////////////////////////
*                          Propagators                                 //
*////////////////////////////////////////////////////////////////////////
         IF (Mode .EQ. 0 ) THEN
            propa =1d0/svar
            propz =1d0/dcmplx( svar -m_MZ**2, svar/m_MZ *m_gammz )
            RSQV=1d0
            RSQA=1d0
         ELSE
* Multiply axial coupling by beta factor. 
* Add formfactors initialisation of s-dependent electro-weak form factors and
* photonic vacuum polarisation 
* (electro-weak box contributions left out here, they depend on acos)
            CALL BornV_GetQCDcor2(KFf,RSQV,RSQA)
            xff(1)=m_GSW(1)
            xff(2)=m_GSW(2)
            xff(3)=m_GSW(3)
            xff(4)=m_GSW(4)
***         xffa  =UNDEFINED !!!!
            xfem  =m_GSW(6)
            xfota =m_GSW(7)
c[[[[[[[[[[[[[[[[[[[[[[!!!!!!!!!!!!!!!!!!!!!!
c            IF(icont.LE.200) THEN
c               write(*,*) '|||||||||||||||||||||BornV|||||||||||||||||||||||||||||||||||||'
c               write(*,'(a,8g22.14)') 'sqrt(svar),costhe=',sqrt(svar),costhe
c               write(*,'(a,8g22.14)') 'xfem=',xfem
c               write(*,'(a,8g22.14)') 'RSQV,RSQA=',RSQV,RSQA
c               xfem  =0d0       !!!!!!!!!!!!!!!!!!!!
c            ELSE
c               STOP
c            ENDIF
c]]]]]]]]]]]]]]]]]]]]]]!!!!!!!!!!!!!!!!!!!!!!
*-------------------------------------------------
            xrho =xff(1)
            xke  =xff(2)
            xkf  =xff(3)
            xkef =xff(4)
            qfm =dabs(qf)
            qem =dabs(qe)
            xe   =  1.d0 -4.d0*m_swsq*qem
            xf   =  1.d0 -4.d0*m_swsq*qfm
            xef  = -1.d0 +xe +xf +16.d0*qem*qfm*m_swsq*m_swsq ! xef=xe*xf !!!
            xve  =  1.d0 -4.d0*m_swsq*qem*xke
            xvf  =  1.d0 -4.d0*m_swsq*qfm*xkf
            xvef = -1.d0 +xve +xvf +16.d0*qem*qfm*m_swsq*m_swsq*xkef
* Multiply axial  coupling by beta factor.
* Multiply vector coupling by form-factor.
* Multiply final vector by RSQV and final axial by RSQA (QCD corrections)
            xupgf(1)=xupgf(1)*RSQV
            xupgf(2)=xupgf(2)*RSQV
            xupzfp(1)=0.5d0*(xupzf(1)+xupzf(2))*xvf/xf*RSQV  +0.5*(xupzf(1)-xupzf(2))*beta*RSQA !
            xupzfp(2)=0.5d0*(xupzf(1)+xupzf(2))*xvf/xf*RSQV  -0.5*(xupzf(1)-xupzf(2))*beta*RSQA !
            xupzip(1)=0.5d0*(xupzi(1)+xupzi(2))*xve/xe  +0.5*(xupzi(1)-xupzi(2)) !
            xupzip(2)=0.5d0*(xupzi(1)+xupzi(2))*xve/xe  -0.5*(xupzi(1)-xupzi(2)) !
* Final state vector coupling
            xupf     =0.5d0*(xupzf(1)+xupzf(2))*xvf/xf*RSQV
* Double vector formfactor thing
            xthing=0.25d0*(xupzf(1)+xupzf(2))*(xupzi(1)+xupzi(2))*(xvef/xef-xvf*xve/xe/xf)*RSQV !
            propa =1d0/svar/(2d0-xfem)
            propz =1d0/dcmplx(svar-m_MZ**2,svar/m_MZ*m_gammz)
* Replace Born normalization of Z propagator by the better one
            del1 =m_Gmu *m_MZ**2 *m_AlfInv/(DSQRT(2.d0)*8.d0*pi)
            del0 =1.d0/(m_swsq*(1.d0-m_swsq))/16.d0
            propz = propz*del1/del0*xrho
c[[[[[[[[[[[[[
c      IF(icont.LE.20) THEN
c         write(*,'(a,5g22.14)') '   propa= ', propa
c         write(*,'(a,5g22.14)') '   propz= ', propz
c         write(*,'(a,5g22.14)') '   xrho = ', xrho
c         write(*,'(a,5g22.14)') '   xke  = ', xke
c         write(*,'(a,5g22.14)') '   xkf  = ', xkf
c         WRITE(*,'(a,5g22.14)') '   xkef = ', xkef/(xke*xkf)
c         write(*,'(a,5g22.14)') '    swsq= ', m_swsq
c      ENDIF
c]]]]]]]]]]]]]
         ENDIF ! (Mode .EQ. 0)
*////////////////////////////////////////////////////////////////////////
*//             Additional Spin amplitudes in neutrino case            //
*////////////////////////////////////////////////////////////////////////
         DO i=1,2
            DO j=1,2
               aw(i,j)=(0.d0,0.d0)
            ENDDO
         ENDDO
         IF (iabs(IVfin) .EQ. 1.and.abs(kff).eq.12) THEN
            IF(Mode .EQ. 0) THEN
               xmw=m_MZ*dsqrt(1d0-m_swsq)
               xmw=m_MW         ! new!!!
               xcoup=1.d0/2.d0/m_swsq
               IF (IVini .LT. 0) THEN
                  aw(1,1)= -DCMPLX(xcoup*(1.d0+CosThe))/xmw/xmw
               ELSE
                  aw(2,2)= -DCMPLX(xcoup*(1.d0+CosThe))/xmw/xmw
               ENDIF
            ELSE
               t= -svar*(1.d0-CosThe)/2d0
               s= svar
cc               xmw=m_MZ*dsqrt(1d0-m_swsq)
cc               xmw=m_MW         ! new!!!
cc               xgw=0d0          ! new!!!
cc               xp2=(svar*(1.d0-CosThe)/2.+xmw*xmw)**2+(xmw*xgw)**2
cc               propw=dcmplx(-(svar*(1.d0-CosThe)/2+xmw*xmw)/xp2)
cc               propw=propw-xi*dcmplx(xmw*xgw/xp2)
               propw=DCMPLX( 1d0/(t-m_MW**2) )
               xcoup=1.d0/2.d0/m_swsq
               del1 =m_Gmu *m_MW**2 *m_AlfInv/(DSQRT(2d0)*pi) ! corrected
               del0 =1.d0/m_swsq/2.d0
               propw=propw*DCMPLX(del1/del0)
c[[[[
cccc           DelW= 1D0/m_AlfInv/m_pi/2*(-3D0/2*LOG(s/m_MW**2)+1D0/2*(LOG(-t/s))**2-m_pi**2/6+2D0)
cccc           ROW=ROW+AL1PI/2*(-3D0/2*ALSMW+1D0/2*(LOG(-TT/S))**2   ! modified
cccc        &                   -2D0/3*PI2+2D0)                      !   row           
c]]]]
               DelW=+( +3d0/2d0*LOG(m_MW**2/amin**2) +0.5d0*(LOG(-t/s))**2 )
     $              -( +3d0/2d0*LOG(s/amin**2) +4d0/6d0*m_pi**2 -2d0) ! (b) true F_1 s-chanel
cccc $              -( +3d0/2d0*LOG(s/amin**2) +1d0/6d0*m_pi**2 -2d0) ! (a) from Zbyszek&Tord
               DelW=DelW *0.5d0/(m_AlfInv*m_pi)
               propw=propw*(m_GSW(5) +DelW)
               IF (IVini .LT. 0) THEN
                  aw(1,1)= propw*dcmplx(xcoup*(1.d0+CosThe)) !orig. 1+CosThe
               ELSE
                  aw(2,2)= propw*dcmplx(xcoup*(1.d0+CosThe)) !orig. 1+CosThe
               ENDIF
            ENDIF
         ENDIF
*////////////////////////////////////////////////////////////////////////
*//             Spin amplitudes   Z+gamma case                         //
*////////////////////////////////////////////////////////////////////////
         DO i=1,2
            DO j=1,2
               regula= (3-2*i)*(3-2*j) + CosThe
               regulm=-(3-2*i)*(3-2*j) * sinthe *2.d0*amfin/sqrt(svar)
               aphot(i,j)=propa*(xupgi(i) *xupgf(j)*regula)
               azett(i,j)=propz*(xupzip(i)*xupzfp(j)+xthing)*regula
               aborn(i,j)=aphot(i,j)+azett(i,j)+aw(i,j)
               aphotm(i,j)= propa*dcmplx(0d0,1d0)  *xupgi(i)*xupgf(j)    *regulm !
               azettm(i,j)= propz*dcmplx(0d0,1d0)*(xupzip(i)*xupf+xthing)*regulm !
               abornm(i,j)=aphotm(i,j)+azettm(i,j)
c[[[[[[[[[[[[[
c               IF(icont.LE.20) THEN
c                  write(*,'(a,2i5,5g22.14)') 'amplit= ',i,j, 
c     $                 propa*xupgi(i) *xupgf(j) + propz*(xupzip(i)*xupzfp(j)+xthing)
c               ENDIF
c]]]]]]]]]]]]]
            ENDDO
         ENDDO
      ENDIF
*////////////////////////////////////////////////////////////////////////
*//    Saving CPU trick ENDs here                                      //
*////////////////////////////////////////////////////////////////////////
*////////////////////////////////////////////////////////////////////////
*//           Differential X-section out of spin amplituds             //
*//  Helicity conservation explicitly obeyed:                          //
*//  Only diagonal elements of the spin density matrices.              //
*//  (Only longitudinal polarizations)                                 //
*////////////////////////////////////////////////////////////////////////
      polar1 =  (eps1)
      polar2 = (-eps2)
      Born   =  0d0
      DO i=1,2
         helic= 3-2*i
         DO j=1,2
            helit=3-2*j
            factor=kolor*(1d0+helic*polar1)*(1d0-helic*polar2)/4d0
            factom=factor*(1+helit*ta)*(1-helit*tb)
            factor=factor*(1+helit*ta)*(1+helit*tb)
            IF(iabs(IVfin) .NE. 1) THEN
*     Normal case
*     (mass terms included in Born. is it better ??????)
               Born=Born+cdabs(aborn(i,j))**2*factor
               IF (Mode .NE. 0) THEN
                  Born=Born+CDABS(abornm(i,j))**2*factom
               ENDIF
            ELSE
*     Neutrino case
               xm2=cdabs( aborn(i,j))**2 !  +(nneut-1)*cdabs(azett(i,j))**2
               xm3=cdabs(abornm(i,j))**2 ! +(nneut-1)*cdabs(azettm(i,j))**2
               Born=Born+(xm2+xm3)*factor
            ENDIF
         ENDDO
      ENDDO
* phase space threshold factor, and multiply by svar**2 to get R-units!
      IF (svar .GT. 4d0*amfin**2) THEN
         thresh=sqrt(1-4d0*amfin**2/svar)
         Born = Born*svar**2*thresh
      ELSE
         Born=0.d0
      ENDIF
      BornV_Dizet = Born
      END


      SUBROUTINE BornV_InterpoGSW(KFf,svar,CosThe)
*//////////////////////////////////////////////////////////////////////////
*//
*//  Calculates GSW formfactors from tables using linear interpolation
*//
*//////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      DOUBLE PRECISION   svar,CosThe
      DOUBLE COMPLEX     xff(4),xfem,xfota
      INTEGER            kk,i,j, KFf
      DOUBLE PRECISION   x,y,h,hy,ww
*----------------------------------------------------------------------------
      ww = SQRT(svar)
      m_WminZ = m_MZ-m_WdelZ
      m_WmaxZ = m_MZ+m_WdelZ

      IF(  (ww.GE.m_WminZ) .AND. (ww.LE.m_WmaxZ)    ) THEN
* LEP1 near Z0 resonance 
         x= (ww-m_WminZ)/(m_WmaxZ-m_WminZ)
         i= MIN( INT(m_poin2*x)+1, m_poin2)
         h= x*m_poin2-DFLOAT(i-1)
         IF(m_poTh2.EQ.0) THEN
            DO kk=1,m_poinG
               m_GSW(kk) = m_czz( i,1,kk,KFf)*(1-h) +m_czz(i+1,1,kk,KFf)*h !
            ENDDO
         ELSE
            y= (1d0+CosThe)/2d0
            j= MIN( INT(m_poTh2*y)+1, m_poTh2)
            hy= y*m_poTh2-DFLOAT(j-1)
            DO kk=1,m_poinG
               m_GSW(kk)=m_czz( i,j,  kk,KFf)*(1-h)*(1d0-hy) +m_czz(i+1,j,  kk,KFf)*h*(1d0-hy) !
     $                  +m_czz( i,j+1,kk,KFf)*(1-h)*hy       +m_czz(i+1,j+1,kk,KFf)*h*hy !
            ENDDO
         ENDIF
         DO kk=1,m_poinQ
            m_QCDcorR(kk)=m_szz( i,kk,KFf)*(1-h) +m_szz(i+1,kk,KFf)*h
         ENDDO
      ELSEIF(  (ww.GE.m_WminLEP1) .AND. (ww.LE.m_WmaxLEP1) ) THEN
* LEP1 outside Z0 and low energies
         x= LOG( ww/m_WminLEP1) / LOG(m_WmaxLEP1/m_WminLEP1)
         i= MIN( INT(m_poin1*x)+1, m_poin1)
         h= x*m_poin1-DFLOAT(i-1)
         DO kk=1,m_poinG
            m_GSW(kk)    =m_cyy( i,kk,KFf)*(1-h) +m_cyy(i+1,kk,KFf)*h
         ENDDO
         DO kk=1,m_poinQ
            m_QCDcorR(kk)=m_syy( i,kk,KFf)*(1-h) +m_syy(i+1,kk,KFf)*h
         ENDDO
      ELSEIF(  (ww.GE.m_WmaxLEP1) .AND. (ww.LE.m_WmaxLEP2) ) THEN
* in the LEP2 region
         x= (ww-m_WmaxLEP1)/(m_WmaxLEP2-m_WmaxLEP1)
         i= MIN( INT(m_poin3*x)+1, m_poin3) 
         h= x*m_poin3-DFLOAT(i-1)
         y= (1d0+CosThe)/2d0
         j= MIN( INT(m_poTh3*y)+1, m_poTh3)
         hy= y*m_poTh3-DFLOAT(j-1)
* EW complex form-factors
         DO  kk=1,m_poinG
            m_GSW(kk)=
     $          m_ctt(i,j  , kk,KFf)*(1-h) *(1d0-hy)  +m_ctt(i+1,  j, kk,KFf) *h*(1d0-hy) !
     $         +m_ctt(i,j+1, kk,KFf)*(1-h) *hy        +m_ctt(i+1,j+1, kk,KFf) *h*hy       !
         ENDDO
* QCD correction
         DO kk=1,m_poinQ
            m_QCDcorR(kk)=m_stt( i,kk,KFf)*(1-h) +m_stt(i+1,kk,KFf)*h
         ENDDO
      ELSEIF(  (ww.GE.m_WmaxLEP2) .AND. (ww.LE.m_WmaxNLC) ) THEN
* in the NLC region
         x= (ww-m_WmaxLEP2)/(m_WmaxNLC-m_WmaxLEP2)
         i= MIN( INT(m_poin4*x)+1, m_poin4)
         h= x*m_poin4-DFLOAT(i-1)
         y= (1d0+CosThe)/2d0
         j= MIN( INT(m_poTh4*y)+1, m_poTh4)
         hy= y*m_poTh4-DFLOAT(j-1)
* EW complex form-factors
         DO  kk=1,m_poinG
            m_GSW(kk)=
     $         m_clc(i,  j, kk,KFf) *(1-h) *(1d0-hy)  +m_clc(i+1,  j, kk,KFf) *h *(1d0-hy)  !
     $        +m_clc(i,j+1, kk,KFf) *(1-h) *hy        +m_clc(i+1,j+1, kk,KFf) *h *hy        !
         ENDDO
* QCD correction
         DO kk=1,m_poinQ
            m_QCDcorR(kk)=m_slc( i,kk,KFf)*(1-h) +m_slc(i+1,kk,KFf)*h
         ENDDO
      ELSE
***         PRINT *,'STOP in BornV_InterpoGSW: s out of predefined range, ww=', ww
         PRINT *,'BornV_InterpoGSW WARNING: s out of predefined range, ww=', ww
         DO kk=1,m_poinG
            m_GSW(kk) = 0d0
         ENDDO
         DO kk=1,m_poinQ
            m_QCDcorR(kk)= 0d0
         ENDDO
***         STOP
      ENDIF
      m_QCDcor = m_QCDcorR(1)-1d0 ! <--- obsolete!!!
      END                        !BornV_GetGSW


      SUBROUTINE BornV_givizo(idferm,ihelic,sizo3,charge,kolor)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*// Provides electric charge and weak izospin of a family fermion where      //
*// idferm =           1,        2,        3,         4,                     //
*// denotes:    neutrino,   lepton,       up,      down   (quark)            //
*// negative idferm=-1,-2,-3,-4, denotes corresponding antiparticle          //
*// ihelic =     +1,  -1   denotes  right and left handednes ( chirality)    //
*// sizo3 is third projection of weak izospin (plus minus half)              //
*// and charge is electric charge in units of electron charge                //
*// kolor is a qcd colour, 1 for lepton, 3 for quarks                        //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INTEGER          idferm,ihelic,kolor
      DOUBLE PRECISION sizo3,charge
*
      INTEGER          lepqua,iupdow,ic,ih,idtype
*------------------------------------------------------------------------------
      IF(idferm  .EQ.  0  .OR.  iabs(idferm)  .GT.  4) GOTO 901
      IF(iabs(ihelic)  .NE.  1)                GOTO 901
      ih  =ihelic
      idtype =iabs(idferm)
      ic  =idferm/idtype
      lepqua=INT(idtype*0.4999999d0)
      iupdow=idtype-2*lepqua-1
      charge  =(-iupdow+2d0/3d0*lepqua)*ic
      sizo3   =0.25d0*(ic-ih)*(1-2*iupdow)
      kolor=1+2*lepqua
* note that conventionaly z0 coupling is
* xoupz=(sizo3-charge*swsq)/sqrt(swsq*(1-swsq))
      RETURN
 901  print *,' STOP in BornV_givizo: wrong params.'
      STOP
      END


      DOUBLE PRECISION  FUNCTION BornV_Sig0nb(CMSene)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//  provides pointlike muon x-section in nanobarns                          //
*//  for normalization purpose                                               //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      DOUBLE PRECISION   pi
      PARAMETER (pi =3.1415926535897932d0)
      DOUBLE PRECISION  CMSene
*---------------------------------
      BornV_Sig0nb =  4d0*pi/(m_AlfInv**2*3d0*CMSene**2)*m_gnanob
      END ! BornV_Sig0nb

*//////////////////////////////////////////////////////////////////////////////
*//                  Getters and setters                                     //
*//////////////////////////////////////////////////////////////////////////////

      SUBROUTINE BornV_GetParticle(KFferm, mass, Qf, T3f, NCf)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      INTEGER KFferm, NCf
      DOUBLE PRECISION   mass, Qf, T3f
      INTEGER KF
      KF = ABS(KFferm)
      mass  = m_amferm( KF)
      Qf    = m_Qf(     KF)
      T3f   = m_T3f(    KF)
      NCf   = m_NCf(    KF)
      END                       ! BornV_GetParticle

      SUBROUTINE BornV_GetIsGenerated(KFferm,IsGenerated)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      INTEGER KFferm, IsGenerated
      IsGenerated   = m_IsGenerated(KFferm)
      END                       ! BornV_GetIsGenerated

      SUBROUTINE BornV_SetIsGenerated(KFferm,IsGenerated)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      INTEGER KFferm, IsGenerated
      m_IsGenerated(KFferm) = IsGenerated
      END                       ! BornV_SetIsGenerated

      SUBROUTINE BornV_SetKF(KFferm)
*//////////////////////////////////////////////////////////////////////////////
*//   set just one fermion                                                   //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      INTEGER KFferm,k
      DO k=1,20
         m_IsGenerated(k)=0
      ENDDO
      m_IsGenerated(KFferm) =1d0
      END                       ! BornV_GetIsGenerated

      DOUBLE PRECISION  FUNCTION BornV_GetMass(KFferm)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      INTEGER KFferm
      BornV_GetMass = m_amferm(ABS(KFferm))
      END ! BornV_GetMass

      DOUBLE PRECISION  FUNCTION BornV_GetAuxPar(KFferm)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      INTEGER KFferm
*------
      BornV_GetAuxPar = m_AuxPar(ABS(KFferm))
      END ! BornV_GetAuxPar

      DOUBLE PRECISION  FUNCTION BornV_GetCharge(KFferm)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      INTEGER KFferm
*------
      BornV_GetCharge = m_Qf(ABS(KFferm))
      END ! BornV_GetCharge

      INTEGER FUNCTION BornV_GetColor(KFferm)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      INTEGER KFferm
*------
      BornV_GetColor = m_NCf(ABS(KFferm))
      END ! BornV_GetColor


      SUBROUTINE BornV_SetKeyQCD(KeyQCD)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      INTEGER KeyQCD
      m_KeyQCD = KeyQCD
      END

      SUBROUTINE BornV_GetKeyQCD(KeyQCD)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      INTEGER KeyQCD
      KeyQCD = m_KeyQCD
      END

      SUBROUTINE BornV_SetKeyElw(KeyElw)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      INTEGER KeyElw
      m_KeyElw = KeyElw
      END

      SUBROUTINE BornV_GetKeyElw(KeyElw)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      INTEGER KeyElw
      KeyElw = m_KeyElw
      END

      SUBROUTINE BornV_GetKeyZet(KeyZet)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      INTEGER KeyZet
*
      KeyZet = m_KeyZet
      END

      SUBROUTINE BornV_SetKeyZet(KeyZet)
*//////////////////////////////////////////////////////////////////////////////
*//   for tests only                                                         //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      INTEGER KeyZet
*
      m_KeyZet = KeyZet
      END

      SUBROUTINE BornV_SetCMSene(CMSene)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      DOUBLE PRECISION  CMSene
*
      m_CMSene = CMSene
      END

      SUBROUTINE BornV_SetMZ(MZ)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      DOUBLE PRECISION  MZ
*------------------
      m_MZ = MZ
      END

      SUBROUTINE BornV_GetMZ(MZ)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      DOUBLE PRECISION  MZ
*------------------
      MZ = m_MZ
      END

      SUBROUTINE BornV_GetGammZ(GammZ)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      DOUBLE PRECISION  GammZ
*------------------
      GammZ = m_GammZ
      END

      SUBROUTINE BornV_GetGmu(Gmu)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      DOUBLE PRECISION  Gmu
*------------------
      Gmu = m_Gmu
      END

      SUBROUTINE BornV_GetSwsq(Swsq)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      DOUBLE PRECISION  Swsq
*------------------
      Swsq = m_swsq
      END

      SUBROUTINE BornV_GetAlfInv(AlfInv)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      DOUBLE PRECISION  AlfInv
*------------------
      AlfInv = m_AlfInv
      END

      SUBROUTINE BornV_GetAvMult(AvMult)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      DOUBLE PRECISION    AvMult
*------------------
      AvMult = m_AvMult
      END

      SUBROUTINE BornV_GetYFSkon(YFSkon)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//   Get finite part of YFS form-factor                                     //
*//   Used in QED3.f                                                         //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      DOUBLE PRECISION    YFSkon
*------------------
      YFSkon = m_YFSkon
      END

      SUBROUTINE BornV_GetYFS_IR(YFS_IR)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//   Get IR (cut-off dependend) part of ISR YFS form-factor                 //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      DOUBLE PRECISION    YFS_IR
*------------------
      YFS_IR = m_YFS_IR

      END

      SUBROUTINE BornV_GetVV(vv)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      DOUBLE PRECISION  vv
*
      vv = m_vv
      END ! BornV_GetVV

      SUBROUTINE BornV_GetVXX(vv,x1,x2)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      DOUBLE PRECISION  vv,x1,x2
*
      vv = m_vv
      x1 = m_x1
      x2 = m_x2
      END ! BornV_GetVXX

      SUBROUTINE BornV_GetGSW(GSW)
*//////////////////////////////////////////////////////////////////////////
*//
*//  Exports  GSW formfactors, called in GPS/CEEX
*//  BornV_InterpoGSW has to be called before
*//
*//////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      DOUBLE COMPLEX    GSW(*)
      INTEGER      k
*     ---------------------------------------------------------------------
      DO k=1,m_poinG
         GSW(k) = m_GSW(k)
      ENDDO
      END

      SUBROUTINE BornV_SetGSW(GSW)
*//////////////////////////////////////////////////////////////////////////
*//
*//  For special tests of pretabulation GSW values can be set form outside
*//
*//////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      DOUBLE COMPLEX    GSW(*)
      INTEGER      k
*     ---------------------------------------------------------------------
      DO k=1,m_poinG
          m_GSW(k) =GSW(k)
      ENDDO
      END

      SUBROUTINE BornV_GetQCDcor2(KFf,RSQV,RSQA)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//   Get QCD correction factors for vercor and axial couplings              //
*//   QED corrections has to be removed from the Dizet QCD formfactor.       //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      INTEGER   KFf
      DOUBLE PRECISION   RSQV,RSQA, QEDcor
      INTEGER            KeyFSR

*     Standard case KeyQCD=1:
*     QCD included through K-factor of Dizet,
*     FSR O(alpha) QED corr. from Dizet is excluded (but alpha*alpha_s is included!)
      IF( ABS(KFf) .LE. 6 ) THEN
         QEDcor = m_Qf(KFf)**2 *3d0/4d0 *m_alfpi
         RSQV = SQRT(m_QCDcorR(1) -QEDcor)        ! Quarks
         RSQA = SQRT(m_QCDcorR(2) -QEDcor)
      ELSE
         RSQV = 1d0             ! Leptons
         RSQA = 1d0
      ENDIF
      IF( m_KeyQCD .EQ. 0) THEN
*     Special KeyQCD=0: for excluding QCD correction whatsoever
         RSQV = 1d0
         RSQA = 1d0
      ENDIF
*
      CALL KK2f_GetKeyFSR(KeyFSR)
      IF( (m_KeyQCD .EQ. 2) .AND. (KeyFSR .EQ .0) ) THEN
*     Special KeyQCD=2: for including QED FSR correction "by-hand"
         IF( ABS(KFf) .LE. 6 ) THEN
            RSQV = SQRT(m_QCDcorR(1)) ! Quarks, QEDcor provided by Dizet
            RSQA = SQRT(m_QCDcorR(2))
         ELSE
            QEDcor = m_Qf(KFf)**2 *3d0/4d0 *m_alfpi
            RSQV = SQRT(1d0+QEDcor) ! Leptons, explicit QEDcor
            RSQA = SQRT(1d0+QEDcor)
         ENDIF
      ENDIF
      END


      SUBROUTINE BornV_GetQCDcorR(QCDcorR)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//   Get QCD correction factor, provided by Dizet                           //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      DOUBLE PRECISION    QCDcorR(m_poinQ)
      INTEGER i
      DO i=1,m_poinQ
         QCDcorR(i) = m_QCDcorR(i)
      ENDDO
      END

      SUBROUTINE BornV_GetQCDcor(QCDcor)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//   Get QCD correction factor, provided by Dizet                           //
*//   !!!!! Obsolete !!!!!!!!!!!!!!
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      DOUBLE PRECISION    QCDcor
      QCDcor = m_QCDcor
      END

      SUBROUTINE BornV_SetQCDcor(QCDcor)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//   Set QCD correction factor, for tests of pretabulation                  //
*//   !!!!! Obsolete !!!!!!!!!!!!!!
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      DOUBLE PRECISION    QCDcor
      m_QCDcor    = QCDcor
      END

      SUBROUTINE BornV_GetMW(MW)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//   Get Mass of W
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      DOUBLE PRECISION    MW
      MW = m_MW
      END

      SUBROUTINE BornV_GetGammW(GammW)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//   Get width of W
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      DOUBLE PRECISION    GammW
      GammW = m_GammW
      END

      SUBROUTINE BornV_GetKeyRes(KeyRes)
*/////////////////////////////////////////////////////////////////////////////////////
*//   KF of beams                                                                   //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      INTEGER KeyRes
*
      KeyRes = m_KeyRes
      END

*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//                      End of CLASS  BornV                                 //
*//////////////////////////////////////////////////////////////////////////////

