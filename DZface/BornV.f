c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*//   Presently ONLY BornV_Initialize() is used in MainTabC.xx
*//   The rest is kept only for reference as a history record
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*//////////////////////////////////////////////////////////////////////////////
*//      This is truncated version of original BornV.f
*//      Only BornV_Initialize is really used in TabMain.f
*//      Other entries kept for hypothetic debug
*//////////////////////////////////////////////////////////////////////////////
**//                     Pseudo-CLASS  BornV                                  //
*//   Provide Born angular distribution and integrated x-section             //
*//   NOTES: Modified for KKMC-hh by S. Yost and S.J.                        //
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
      INCLUDE 'hhDizet.h'
      DOUBLE PRECISION  xpar_input(*)
      DOUBLE PRECISION  Npar(30), alfQCDMZ, AlStrZ, AlStrT, AlfinvMZ,
     &     AlQedZ, DAL5H, zpard(30), partz(0:11), partw(3)
*------------------------------------------------------------------------------
      INTEGER             k,j,kxpa,KF
*------------------------------------------------------------------------------
      write(*,*) '***********************************************'
      write(*,*) '************BornV_Initialize*******************'

      DO j=1,10000
         m_xpar_input(j)=xpar_input(j)
      ENDDO


      m_QCDcor = 0d0
      DO k=1,m_poinQ
         m_QCDcorR(k)=0d0
      ENDDO
      m_CMSene = xpar_input( 1)         ! Central value of CMS energy, do not change!
      m_XXXene = m_CMSene               ! Just initialization only
      m_KFini = xpar_input( 400)        ! KFcode of beam, POSITIVE!!!
      m_KeyFSR= xpar_input(  21)
*                      <<<  ff-pair spectrum >>>
      m_vvmin  = xpar_input(16)         ! minimum v, infrared cut
      m_vvmax  = xpar_input(17)         ! maximum v, infrared cut (may redefine later)
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
      m_Gammz  = xpar_input(504)        ! Z width
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
      m_out = 16
*
cc      IF (xpar_input(400).EQ.0d0) xpar_input(400) = 1d0
cc      IF (m_KeyElw .GE. 1) CALL BornV_StartEW(xpar_input)
cc      xpar_input(400) = m_KFini ! (reset to what it was)

      WRITE(m_out,bxope)
      WRITE(m_out,bxtxt) '  BornV  Initialization                '
      WRITE(m_out,bxl1f) m_MZ    ,   'Z mass     [GeV]   ','amz   ','a1'
      WRITE(m_out,bxl1f) m_amh   ,   'Higgs mass [GeV]   ','amh   ','a2'
      WRITE(m_out,bxl1f) m_amtop ,   'Top mass   [GeV]   ','amtop ','a3'
      WRITE(m_out,bxl1f) m_gammz,    'Z width    [GeV]   ','gammz ','a4'
      WRITE(m_out,bxl1f) m_swsq,     'sin(theta_w)**2    ','sinw2 ','a5'
      WRITE(m_out,bxl1f) m_AlfInv,   '1/alfa_QED  at  Q=0','AlfInv','a6'
      WRITE(m_out,bxl1f) m_HadMin,   'MassCut light qqbar','HadMin','a6'
      WRITE(m_out,bxl1i) m_KFini ,   'KF code of beam    ','KFini ','a7'
      WRITE(m_out,bxl1g) m_vvmax,    'Input vvmax        ','vvmax ','a8'
      WRITE(m_out,bxtxt) 'Test switches:                         '
      WRITE(m_out,bxl1i) m_KeyElw,   'Electroweak lib.   ','KeyElw','10'
      WRITE(m_out,bxl1i) m_KeyZet,   'Z switch           ','KeyZet','11'
      WRITE(m_out,bxl1i) m_KeyWtm,   'mass terms on/off  ','KeyWtm','12'
      WRITE(m_out,bxl1i) m_KeyRes,   'R for gamma* on/off','KeyRes','12'
      WRITE(m_out,bxclo)
      END

      DOUBLE PRECISION FUNCTION BornV_Dizet(Mode,KFi,KFf,svar,CosThe,eps1,eps2,ta,tb)
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c!!!! translated to c++ but keep it for debug for some time
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//   Calculates differential born cross section.                            //
*//   For Mode=0 pure Born and for Mode=1 electroweak corrs. are added.      //
*//   KFi,KFf can be also negative ??? for antiparticle, in this case it is  //
*//   important to produce tables with correct input KFini, KFfin !!!        //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      INCLUDE 'hhDizet.h'
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
      INTEGER          j,i,ivini,kdumm,kfi0,kff0,mode0,ivfin
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
      DATA Mode0,svar0,cost0,KFi0,KFf0 /-155,-155.d0,-156.d0,-99,-99/
*-------------------------------------------------------------------------------------

      icont=icont+1
*////////////////////////////////////////////////////////////////////////
*//    Save CPU for the same svar, CosThe  and varying spins ta, tb    //
*////////////////////////////////////////////////////////////////////////
      IF (Mode.NE.Mode0 .OR. svar.NE.svar0 .OR. CosThe.NE.cost0 
     $   .OR. KFi.NE.KFi0 .OR. KFf.NE.KFf0  ) THEN
         Mode0  = Mode
         svar0  = svar
         cost0  = CosThe
         KFi0   = KFi
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
c[[[[[!!!!!!!!!!!!!!!!!!!!!!!
c      IF(icont .EQ. 1) write(*,*) " m_swsq=",m_swsq
c      IF(icont .EQ. 1) write(*,*) " xupzi =",xupzi(1),xupzi(2)
c      IF(icont .EQ. 1) write(*,*) " xupzf =",xupzf(1),xupzf(2)
c]]]]]!!!!!!!!!!!!!!!!!!!!!!!
*
         sinthe = sqrt(1.d0-CosThe**2)
         beta   = SQRT(MAX(0d0,1d0-4d0*amfin**2/svar))

* Multiply axial coupling by beta factor.
         xupzfp(1)= 0.5d0*(xupzf(1)+xupzf(2))+0.5*beta*(xupzf(1)-xupzf(2))
         xupzfp(2)= 0.5d0*(xupzf(1)+xupzf(2))-0.5*beta*(xupzf(1)-xupzf(2))
         xupzip(1)= 0.5d0*(xupzi(1)+xupzi(2))     +0.5*(xupzi(1)-xupzi(2))
         xupzip(2)= 0.5d0*(xupzi(1)+xupzi(2))     -0.5*(xupzi(1)-xupzi(2))
* Final state vector coupling
         xupf     = 0.5d0*(xupzf(1)+xupzf(2))
         xupi     = 0.5d0*(xupzi(1)+xupzi(2))
         xthing   = 0d0
c[[[[[!!!!!!!!!!!!!!!!!!!!!!!
c      IF(icont .EQ. 1) write(*,*) " xupzfp =",xupzfp(1),xupzfp(2)
c      IF(icont .EQ. 1) write(*,*) " xupzip =",xupzip(1),xupzip(2)
c]]]]]!!!!!!!!!!!!!!!!!!!!!!!

*////////////////////////////////////////////////////////////////////////
*                          Propagators                                 //
*////////////////////////////////////////////////////////////////////////
         IF (Mode .EQ. 0 ) THEN
            propa =1d0/svar
            IF ((m_KeyZet.EQ.-1).OR.(m_KeyZet.EQ.-2)) THEN 
               ! fixed width option for KKMC-hh
               propz =1d0/dcmplx( svar -m_MZ**2, m_gammz*m_MZ)
            ELSE
               propz =1d0/dcmplx( svar -m_MZ**2, svar/m_MZ *m_gammz )
            ENDIF
            IF (m_KeyZet.EQ.-1) THEN 
               ! fixed width with redefined parameters
               propz = propz/dcmplx(1d0,m_gammz/m_MZ)
            END IF 
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
c             IF(icont.EQ.1)  write(*,'(a,8g22.14)') 'RSQV,RSQA=',RSQV,RSQA
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
            IF ((m_KeyZet.EQ.-1).OR.(m_KeyZet.EQ.-2)) THEN 
              !!! KKMC-hh: implement fixed width option...
              !!! but note: fixed width isn't really compatible with Dizet.
               propz =1d0/dcmplx(svar-m_MZ**2,m_MZ*m_gammz)
            ELSE
               propz =1d0/dcmplx(svar-m_MZ**2,svar/m_MZ*m_gammz)
            END IF
            IF (m_KeyZet.EQ.-1) THEN !!! fixed width parameter redefinition
               propz = propz/dcmplx(1d0,m_gammz/m_MZ)
            END IF
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
* KKMC-hh: electron neutrino is not special
*        IF (iabs(IVfin) .EQ. 1.and.abs(kff).eq.12) THEN
*           IF(Mode .EQ. 0) THEN
*              xmw=m_MZ*dsqrt(1d0-m_swsq)
*              xmw=m_MW         ! new!!!
*              xcoup=1.d0/2.d0/m_swsq
*              IF (IVini .LT. 0) THEN
*                 aw(1,1)= -DCMPLX(xcoup*(1.d0+CosThe))/xmw/xmw
*              ELSE
*                 aw(2,2)= -DCMPLX(xcoup*(1.d0+CosThe))/xmw/xmw
*              ENDIF
*           ELSE
*              t= -svar*(1.d0-CosThe)/2d0
*              s= svar
cc               xmw=m_MZ*dsqrt(1d0-m_swsq)
cc               xmw=m_MW         ! new!!!
cc               xgw=0d0          ! new!!!
cc               xp2=(svar*(1.d0-CosThe)/2.+xmw*xmw)**2+(xmw*xgw)**2
cc               propw=dcmplx(-(svar*(1.d0-CosThe)/2+xmw*xmw)/xp2)
cc               propw=propw-xi*dcmplx(xmw*xgw/xp2)
*              propw=DCMPLX( 1d0/(t-m_MW**2) )
*              xcoup=1.d0/2.d0/m_swsq
*              del1 =m_Gmu *m_MW**2 *m_AlfInv/(DSQRT(2d0)*pi) ! corrected
*              del0 =1.d0/m_swsq/2.d0
*              propw=propw*DCMPLX(del1/del0)
c[[[[
cccc           DelW= 1D0/m_AlfInv/m_pi/2*(-3D0/2*LOG(s/m_MW**2)+1D0/2*(LOG(-t/s))**2-m_pi**2/6+2D0)
cccc           ROW=ROW+AL1PI/2*(-3D0/2*ALSMW+1D0/2*(LOG(-TT/S))**2   ! modified
cccc        &                   -2D0/3*PI2+2D0)                      !   row           
c]]]]
*              DelW=+( +3d0/2d0*LOG(m_MW**2/amin**2) +0.5d0*(LOG(-t/s))**2 )
*    $              -( +3d0/2d0*LOG(s/amin**2) +4d0/6d0*m_pi**2 -2d0) ! (b) true F_1 s-chanel
cccc $              -( +3d0/2d0*LOG(s/amin**2) +1d0/6d0*m_pi**2 -2d0) ! (a) from Zbyszek&Tord
*              DelW=DelW *0.5d0/(m_AlfInv*m_pi)
*              propw=propw*(m_GSW(5) +DelW)
*              IF (IVini .LT. 0) THEN
*                 aw(1,1)= propw*dcmplx(xcoup*(1.d0+CosThe)) !orig. 1+CosThe
*              ELSE
*                 aw(2,2)= propw*dcmplx(xcoup*(1.d0+CosThe)) !orig. 1+CosThe
*              ENDIF
*           ENDIF
*        ENDIF
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
c[[[[[!!!!!!!!!!!!!!!!!!!!!!!
c      IF(icont .LE. 200) THEN
c        write(*,*) "####### KFi,KFf,svar,CosThe=",KFi,KFf,svar,CosThe
c        write(*,*) "####### xff(1,2,3)",xff(1),xff(2),xff(3)
c        write(*,*) "####### Born =",Born
c      ENDIF
c]]]]]!!!!!!!!!!!!!!!!!!!!!!!

      BornV_Dizet = Born
      END

      SUBROUTINE BornV_PrintGSW(nout,KFi,KFf,ww,CosThe)
*//////////////////////////////////////////////////////////////////////////
*//   This is for debug purpose, to be kept
*//////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      INCLUDE 'hhDizet.h'
      DOUBLE PRECISION   svar,CosThe,ww
      INTEGER nout,KFi,KFf

cc      CALL BornV_SetKFini(KFi)
      m_KFini =  KFi
      svar = ww**2
      CALL BornV_InterpoGSW(KFf,svar,CosThe)
      write(nout,*) '%%%%%% BornV_PrintGSW: KFi,KFf,ww,CosThe=', KFi,KFf,ww,CosThe
      write(nout,*) '%%%%%% GSW(1,2,3)',m_GSW(1),m_GSW(2),m_GSW(3)

      END ! BornV_PrintGSW

      SUBROUTINE BornV_InterpoGSW(KFf,svar,CosThe)
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c  translated to c++ but keep it for debug
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*//////////////////////////////////////////////////////////////////////////
*//  KKMC-hh Version. - SY
*//  Calculates GSW formfactors from tables using linear interpolation
*//  of the matrix in KFi and KFf constructed in hhDizet.
*//  For compatibility with KKMC, KFi is passed in the BornV common block.
*//  Modified for KKMC-hh so it behaves gracefully beyond table limits.
*//  A new calculation is done if the request is outside the table.
*//
*//////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      INCLUDE 'hhDizet.h'
      DOUBLE PRECISION   svar,CosThe
      DOUBLE COMPLEX     xff(4),xfem,xfota
      INTEGER            kk,i,j, KFi, KFf, KFl
      DOUBLE PRECISION   x,y,h,hy,ww
*----------------------------------------------------------------------------
      ww = SQRT(svar)
      m_WminZ = m_MZ-m_WdelZ
      m_WmaxZ = m_MZ+m_WdelZ
      KFi = ABS(m_KFini)
* First generation used in place of 2nd or 3rd, since only charges matter.
      IF ((KFi.EQ.3).OR.(KFi.EQ.5)) KFi = 1
      IF (KFi.EQ.4) KFi = 2
* Lepton prototypes use 2nd generation
      IF ((KFf.EQ.11).OR.(KFf.EQ.13).OR.(KFf.EQ.15)) THEN
         KFl = 13
      ELSE ! neutrino: must match initialization in hhDizet.f
         KFl = 14
      END IF

* Defaults, without FSR QCD. No need to interpolate if KeyQCD = 0.
      m_QCDcorR(1) = 1d0
      m_QCDcorR(2) = 1d0
      m_QCDcorR(3) = 0d0
      m_QCDcorR(4) = 0d0

      IF(  (ww.GE.m_WminZ) .AND. (ww.LE.m_WmaxZ)    ) THEN
* LEP1 near Z0 resonance 
         x= (ww-m_WminZ)/(m_WmaxZ-m_WminZ)
         i= MIN( INT(m_poin2*x)+1, m_poin2)
         h= x*m_poin2-DFLOAT(i-1)
         IF(m_poTh2.EQ.0) THEN
            DO kk=1,m_poinG
               m_GSW(kk) =m_czzs( i,1,kk,KFi,KFl)*(1-h)
     $                   +m_czzs(i+1,1,kk,KFi,KFl)*h !
            ENDDO
         ELSE
            y= (1d0+CosThe)/2d0
            j= MIN( INT(m_poTh2*y)+1, m_poTh2)
            hy= y*m_poTh2-DFLOAT(j-1)
            DO kk=1,m_poinG
               m_GSW(kk)=m_czzs( i,j,   kk,KFi,KFl)*(1-h)*(1d0-hy)
     $                  +m_czzs(i+1,j,  kk,KFi,KFl)*h*(1d0-hy)
     $                  +m_czzs( i,j+1, kk,KFi,KFl)*(1-h)*hy
     $                  +m_czzs(i+1,j+1,kk,KFi,KFl)*h*hy !
            ENDDO
         ENDIF
         IF (m_KeyQCD.NE.0) THEN
           DO kk=1,m_poinQ
               m_QCDcorR(kk)=m_szzs( i,kk,KFi,KFl)*(1-h)
     $                      +m_szzs(i+1,kk,KFi,KFl)*h
           END DO
        ENDIF
      ELSEIF(  (ww.GE.m_WminLEP1) .AND. (ww.LE.m_WmaxLEP1) ) THEN
* LEP1 outside Z0 and low energies
         x= LOG( ww/m_WminLEP1) / LOG(m_WmaxLEP1/m_WminLEP1)
         i= MIN( INT(m_poin1*x)+1, m_poin1)
         h= x*m_poin1-DFLOAT(i-1)
c[[[[[[[[[[[[[
c         write(*,*) '########## BornV_InterpoGSW: ww, KFi,KFl =',ww, KFi,KFl
c         write(*,*) '########## BornV_InterpoGSW: x,i,h       =',x,i,h
c]]]]]]]]]]]]]
         DO kk=1,m_poinG
            m_GSW(kk)    =m_cyys(  i,kk,KFi,KFl)*DCMPLX(1-h)
     $                   +m_cyys(i+1,kk,KFi,KFl)*DCMPLX(h)
         ENDDO
c[[[[[[[[[[[[[
c         write(*,*) '########## BornV_InterpoGSW: m_cyys(  i, 3 ,KFi,KFl) =', m_cyys(  i, 3 ,KFi,KFl)
c         write(*,*) '########## BornV_InterpoGSW: m_cyys(i+1, 3 ,KFi,KFl) =', m_cyys(i+1, 3 ,KFi,KFl)
c         write(*,*) '########## BornV_InterpoGSW: m_GSW(3)                =', m_GSW(3)
c]]]]]]]]]]]]]
         IF (m_KeyQCD.NE.0) THEN
            DO kk=1,m_poinQ
               m_QCDcorR(kk)=m_syys( i,kk,KFi,KFl)*(1-h)
     $                      +m_syys(i+1,kk,KFi,KFl)*h
            ENDDO
         END IF
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
     $          m_ctts(i,j  , kk,KFi,KFl)*(1-h) *(1d0-hy)
     $         +m_ctts(i+1,  j, kk,KFi,KFl) *h*(1d0-hy)
     $         +m_ctts(i,j+1, kk,KFi,KFl)*(1-h) *hy
     $         +m_ctts(i+1,j+1, kk,KFi,KFl) *h*hy       !
         ENDDO
* QCD correction
         IF (m_KeyQCD.NE.0) THEN
            DO kk=1,m_poinQ
               m_QCDcorR(kk)=m_stts( i,kk,KFi,KFl)*(1-h)
     $                      +m_stts(i+1,kk,KFi,KFl)*h
            ENDDO
         END If
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
     $         m_clcs(i,  j, kk,KFi,KFl) *(1-h) *(1d0-hy)
     $        +m_clcs(i+1,  j, kk,KFi,KFl) *h *(1d0-hy)  !
     $        +m_clcs(i,j+1, kk,KFi,KFl) *(1-h) *hy
     $        +m_clcs(i+1,j+1, kk,KFi,KFl) *h *hy        !
         ENDDO
* QCD correction
         IF (m_KeyQCD.NE.0) THEN
            DO kk=1,m_poinQ
               m_QCDcorR(kk)=m_slcs( i,kk,KFi,KFl)*(1-h)
     $                      +m_slcs(i+1,kk,KFi,KFl)*h
            ENDDO
         END IF
      ELSE
*!! New for KKMC-hh: Calculate values that cannot be interpolated.
         write(*,*) '@@@@@ STOP in BornV_InterpoGSW: out of interpolation range. '
         STOP
      ENDIF
c[[[[[[[[[[[[[
c         write(*,*) '########## BornV_InterpoGSW: ww, KFi,KFl =',ww, KFi,KFl
c         write(*,*) '########## BornV_InterpoGSW: x,i,h       =',x,i,h
c         write(*,*) '########## BornV_InterpoGSW: y,j,hy      =',y,j,hy
c]]]]]]]]]]]]]
      m_QCDcor = m_QCDcorR(1)-1d0 ! <--- obsolete!!!
      END                        !BornV_GetGSW



      SUBROUTINE BornV_givizo(idferm,ihelic,sizo3,charge,kolor)
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
* keep it for debug
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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


*//////////////////////////////////////////////////////////////////////////////
*//                  Getters and setters                                     //
*//////////////////////////////////////////////////////////////////////////////

      SUBROUTINE BornV_GetQCDcor2(KFf,RSQV,RSQA)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//   Get QCD correction factors for vercor and axial couplings              //
*//   QED corrections has to be removed from the Dizet QCD formfactor.       //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      INCLUDE 'hhDizet.h'
      INTEGER   KFf
      DOUBLE PRECISION   RSQV,RSQA, QEDcor
cc      INTEGER            KeyFSR

*     Standard case KeyQCD=1:
*     QCD included through K-factor of Dizet,
*     FSR O(alpha) QED corr. from Dizet is excluded (but alpha*alpha_s is included!)
*** Reorganized for efficiency with KKMC-hh, since QCD FSR is turned off...
*** (final state quarks unsupported in KKMC-hh so KFf > 6 always)
      RSQV = 1d0
      RSQA = 1d0
      IF ((m_KeyQCD.EQ.1).AND.(ABS(KFf).LE.6 )) THEN
         QEDcor = m_Qf(KFf)**2 *3d0/4d0 *m_alfpi
         RSQV = SQRT(m_QCDcorR(1) -QEDcor)        ! Quarks
         RSQA = SQRT(m_QCDcorR(2) -QEDcor)
      ELSE IF (m_KeyQCD.EQ.2) THEN
*     Special KeyQCD=2: for including QED FSR correction "by-hand"
cc         CALL KK2f_GetKeyFSR(KeyFSR)
         IF (m_KeyFSR.EQ.0) THEN
            IF( ABS(KFf) .LE. 6 ) THEN
               RSQV = SQRT(m_QCDcorR(1)) ! Quarks, QEDcor provided by Dizet
               RSQA = SQRT(m_QCDcorR(2))
            ELSE
               QEDcor = m_Qf(KFf)**2 *3d0/4d0 *m_alfpi
               RSQV = SQRT(1d0+QEDcor) ! Leptons, explicit QEDcor
               RSQA = SQRT(1d0+QEDcor)
            END IF
         ENDIF
      ENDIF
      END

*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//                      End of CLASS  BornV                                 //
*//////////////////////////////////////////////////////////////////////////////

