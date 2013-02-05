*/////////////////////////////////////////////////////////////////////////////////
*//   make tabcheck
*/////////////////////////////////////////////////////////////////////////////////
      PROGRAM tabtest
*/////////////////////////////////////////////////////////////////////////////////
*//                                                                             //
*/////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'TabCheck.h'
      INTEGER KFdown, KFup, KFstran, KFcharm, KFbotom, KFtop
      PARAMETER( 
     $     KFdown  = 1,   KFup    = 2,
     $     KFstran = 3,   KFcharm = 4,
     $     KFbotom = 5,   KFtop   = 6)
      INTEGER KFel,KFelnu,KFmu,KFmunu,KFtau,KFtaunu
      PARAMETER(
     $     KFel    = 11,  KFelnu  = 12,
     $     KFmu    = 13,  KFmunu  = 14,
     $     KFtau   = 15,  KFtaunu = 16)
      DOUBLE PRECISION    ene
      DOUBLE PRECISION    BornV_Make
      DOUBLE PRECISION    CosTheta,bor1,bor2
      INTEGER             ninp2, nout
      INTEGER             i,ibox,iloop
      CHARACTER*1         chr1
      CHARACTER*40        TableFile, BenchFile
*
      INTEGER            jmax
      PARAMETER         (jmax = 10000)
      DOUBLE PRECISION   xpar(jmax)
*----------------------------------------------------------------------------
      ninp2 =9
      nout =16
      m_Mode  = 0               ! no EWK
      m_Mode  = 1
      CALL ReaDataY('./.KK2f_defaults',1,0,jmax,xpar)
      OPEN(nout,FILE='./TabTest.output')
      m_CMSene = xpar( 1)
      m_CMSene =  15d0
      m_CMSene = 200d0
      m_MZ     = xpar(502)        ! Z mass [GeV]
      m_amh    = xpar(505)        ! Higgs mass, Input for Dizet
      m_amtop  = xpar(506)        ! Top mass,   Input for Dizet
      m_KFini  = KFel
      WRITE(*,*) 'Initialize BornV'
      CALL BornV_Initialize(xpar)
      CALL BornV_SetKeyElw(1)     ! reset KeyElw ?????
****  CALL BornV_SetKeyZet(0)     ! switch off Z
      CALL BornV_GetGammZ(m_GammZ)
      WRITE(*,*) ' -------------- test1 --------------'
      m_KFfin  = KFbotom
      m_KFfin  = KFdown
      m_KFfin  = KFup
      m_KFfin  = KFtau
      m_KFfin  = KFmu
*----------------------------------------------------------------------------
* DIZET initialization, set EW params
      WRITE(*,*) 'Enter DZface'
      ibox  = 1       !!!
      CALL DZface_Initialize( m_KFini,  m_KFfin, m_MZ, m_amh, m_amtop, ibox,  nout)
      m_seps1 = 0d0
      m_seps2 = 0d0
      m_ta    = 0d0
      m_tb    = 0d0
*----------------------------------------------------------------------------
* test of look-up tables
      CALL TabCheck_Lep1       ! SigDif=0.36E-03, AfbDif=0.11E-04 for mu
      CALL TabCheck_NearZ      ! SigDif=.120E-03, AfbDif= .964E-04 for mu
      CALL TabCheck_LEP2       ! SigDif=.187E-04, AfbDif= .144E-03 for mu
      CALL TabCheck_NLC        ! SigDif=.549E-03, AfbDif= .964E-04 for mu
      CLOSE(ninp2)
      END


      SUBROUTINE TabCheck_Lep1
*/////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'TabCheck.h'
      DOUBLE PRECISION   SigDifMax, AfbDifMax, SigDif, AfbDif
      DOUBLE PRECISION   born0,svar,BornT,BornD,AfbT,AfbD
      DOUBLE PRECISION   BornV_Sig0nb,Sig0nb,X
      DOUBLE PRECISION   BornV_Differential
      INTEGER            i
*----------------------------------------------------------------------------
      SigDifMax =-1d0
      AfbDifMax =-1d0
      DO i=m_poin1,1,-1
         IF( MOD(i,10).EQ.0) WRITE(*,*) 'Lep1 i= ',i,'/',m_poin1
         x  =DFLOAT(i)/DFLOAT(m_poin1) ! weak test
         x  =(i-0.5d0)/DFLOAT(m_poin1) ! strong test
         m_CMSene =m_WminLEP1*(m_WmaxLEP1/m_WminLEP1)**x
         Sig0nb  = BornV_Sig0nb(m_CMSene)
         svar    = m_CMSene**2
         born0   = Sig0nb*BornV_Differential( m_Mode,m_KFfin,svar,0d0,m_seps1,m_seps2,m_ta,m_tb)
         CALL TabCheck_MakeXTab(BornT,AfbT)
         CALL TabCheck_MakeXDiz(BornD,AfbD)
         SigDif =0d0
         IF(BornD.NE.0d0) SigDif = ABS(Bornt/BornD-1)
         AfbDif = ABS(AfbT-AfbD)
         IF(SigDif.GT.SigDifMax) THEN
            SigDifMax=SigDif
            WRITE( *,'(a,f15.7,9g19.10)')'CMSene, BornT, BornD, SigDif',m_CMSene,BornT, BornD, SigDif
         ENDIF
         IF(AfbDif.GT.AfbDifMax) THEN
            AfbDifMax=AfbDif
            WRITE( *,'(a,f15.7,9g17.10)')'CMSene,   AfbT, AfbD, AfbDif',m_CMSene,AfbT, AfbD, AfbDif
         ENDIF
      ENDDO
*----------------------------------------------------------------------------
      END

      SUBROUTINE TabCheck_NearZ
*/////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'TabCheck.h'
      DOUBLE PRECISION   SigDifMax, AfbDifMax, SigDif, AfbDif
      DOUBLE PRECISION   BornT,BornD,AfbT,AfbD
      DOUBLE PRECISION   BornV_Sig0nb,Sig0nb,X
      DOUBLE PRECISION   BornV_Differential
      INTEGER            i
*----------------------------------------------------------------------------
      SigDifMax =-1d0
      AfbDifMax =-1d0
      m_WminZ =m_MZ-2d0*m_GammZ
      m_WmaxZ =m_MZ+2d0*m_GammZ
      DO i=1,m_poin2
         IF( MOD(i,10).EQ.1) WRITE(*,*) 'NearZ i= ',i,'/',m_poin2
         X  =DFLOAT(i)/DFLOAT(m_poin2) !  weak test,   end of bins
         X  =(i-0.5d0)/DFLOAT(m_poin2) !  strong test, middle of
         m_CMSene =m_WminZ+(m_WmaxZ-m_WminZ)*X
         Sig0nb  = BornV_Sig0nb(m_CMSene)
         CALL TabCheck_MakeXTab(BornT,AfbT)
         CALL TabCheck_MakeXDiz(BornD,AfbD)
         SigDif =0d0
         IF(BornD.NE.0d0) SigDif = ABS(Bornt/BornD-1)
         AfbDif = ABS(AfbT-AfbD)
         IF(SigDif.GT.SigDifMax) THEN
            SigDifMax=SigDif
            WRITE( *,'(a,f15.7,9g19.10)')'CMSene, BornT, BornD, SigDif',m_CMSene,BornT, BornD, SigDif
         ENDIF
         IF(AfbDif.GT.AfbDifMax) THEN
            AfbDifMax=AfbDif
            WRITE( *,'(a,f15.7,9g17.10)')'CMSene,   AfbT, AfbD, AfbDif',m_CMSene,AfbT, AfbD, AfbDif
         ENDIF
      ENDDO
*----------------------------------------------------------------------------
      END

      SUBROUTINE TabCheck_LEP2
*/////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'TabCheck.h'
      DOUBLE PRECISION   SigDifMax, AfbDifMax, SigDif, AfbDif
      DOUBLE PRECISION   BornT,BornD,AfbT,AfbD
      DOUBLE PRECISION   BornV_Sig0nb,Sig0nb,X
      DOUBLE PRECISION   BornV_Differential
      INTEGER            i
*----------------------------------------------------------------------------
      SigDifMax =-1d0
      AfbDifMax =-1d0
      DO  i=1,m_poin3
         IF( MOD(i,10).EQ.1) WRITE(*,*) 'LEP2 i= ',i,'/',m_poin3
         X  = DFLOAT(i)/DFLOAT(m_poin3)
         X  = (i-0.5d0)/DFLOAT(m_poin3) !  strong test, middle of bin
         m_CMSene =m_WmaxLEP1+(m_WmaxLEP2-m_WmaxLEP1)*X
         Sig0nb  = BornV_Sig0nb(m_CMSene)
         CALL TabCheck_MakeXTab(BornT,AfbT)
         CALL TabCheck_MakeXDiz(BornD,AfbD)
         SigDif =0d0
         IF(BornD.NE.0d0) SigDif = ABS(Bornt/BornD-1)
         AfbDif = ABS(AfbT-AfbD)
         IF(SigDif.GT.SigDifMax) THEN
            SigDifMax=SigDif
            WRITE( *,'(a,f15.7,9g19.10)')'CMSene, BornT, BornD, SigDif',m_CMSene,BornT, BornD, SigDif
         ENDIF
         IF(AfbDif.GT.AfbDifMax) THEN
            AfbDifMax=AfbDif
            WRITE( *,'(a,f15.7,9g17.10)')'CMSene,   AfbT, AfbD, AfbDif',m_CMSene,AfbT, AfbD, AfbDif
         ENDIF
      ENDDO
*----------------------------------------------------------------------------
      END

      SUBROUTINE TabCheck_NLC
*/////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'TabCheck.h'
      DOUBLE PRECISION   SigDifMax, AfbDifMax, SigDif, AfbDif
      DOUBLE PRECISION   BornT,BornD,AfbT,AfbD
      DOUBLE PRECISION   BornV_Sig0nb,Sig0nb,X
      DOUBLE PRECISION   BornV_Differential
      INTEGER            i
*----------------------------------------------------------------------------
      SigDifMax =-1d0
      AfbDifMax =-1d0
      DO  i=1,m_poin4
         IF( MOD(i,10).EQ.1) WRITE(*,*) 'NLC i= ',i,'/',m_poin3
         X   = DFLOAT(i)/DFLOAT(m_poin4)
         X   = (i-0.5d0)/DFLOAT(m_poin4) !  strong test, middle of bin
         m_CMSene  =m_WmaxLEP2+(m_WmaxNLC-m_WmaxLEP2)*X
         Sig0nb  = BornV_Sig0nb(m_CMSene)
         CALL TabCheck_MakeXTab(BornT,AfbT)
         CALL TabCheck_MakeXDiz(BornD,AfbD)
         SigDif =0d0
         IF(BornD.NE.0d0) SigDif = ABS(Bornt/BornD-1)
         AfbDif = ABS(AfbT-AfbD)
         IF(SigDif.GT.SigDifMax) THEN
            SigDifMax=SigDif
            WRITE( *,'(a,f15.7,9g19.10)')'CMSene, BornT, BornD, SigDif',m_CMSene,BornT, BornD, SigDif
         ENDIF
         IF(AfbDif.GT.AfbDifMax) THEN
            AfbDifMax=AfbDif
            WRITE( *,'(a,f15.7,9g17.10)')'CMSene,   AfbT, AfbD, AfbDif',m_CMSene,AfbT, AfbD, AfbDif
         ENDIF
      ENDDO
*----------------------------------------------------------------------------
      END

      SUBROUTINE TabCheck_MakeXTab(Xsec,Afb)
*/////////////////////////////////////////////////////////////////////////////////
*//   xsection and Afb from pretabulation                                       //
*/////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'TabCheck.h'
      DOUBLE PRECISION   Xsec,Afb
      DOUBLE PRECISION   BornV_Sig0nb,Sig0nb,cmax,xsB,xsF
      DOUBLE PRECISION   Mathlib_Gauss
      DOUBLE PRECISION   TabCheck_DTheTab
      EXTERNAL           TabCheck_DTheTab
*     ----------------------------------------------------------
      Sig0nb   = BornV_Sig0nb(m_CMSene)
      cmax=  0.99999999999d0
      CALL TabCheck_Gaus16(TabCheck_DTheTab, -cmax,  0d0, xsB) !!! 16-point gauss
      CALL TabCheck_Gaus16(TabCheck_DTheTab,   0d0, cmax, xsF) !!! 16-point gauss
      Xsec=Sig0nb*(xsB+xsF)
      Afb =0d0
      IF(Xsec.NE.0d0) Afb = (xsF-xsB)/(xsF+xsB)
      END

      SUBROUTINE TabCheck_MakeXDiz(Xsec,Afb)
*/////////////////////////////////////////////////////////////////////////////////
*//   xsection and Afb from directly from dizet                                 //
*/////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'TabCheck.h'
      DOUBLE PRECISION   Xsec,Afb
      DOUBLE PRECISION   BornV_Sig0nb,Sig0nb,dummy,cmax,xsB,xsF
      DOUBLE PRECISION   Mathlib_Gauss
      DOUBLE PRECISION   TabCheck_DTheDiz
      EXTERNAL           TabCheck_DTheDiz
*     ----------------------------------------------------------
      CALL DZface_MakeGSW(-1, m_CMSene, 0d0, m_GSW, m_QCDcor) !initialize QCDcor calculation.
      CALL DZface_MakeGSW( 0, m_CMSene, 0d0, m_GSW, m_QCDcor) !calculate, just in case we need
      CALL BornV_SetGSW(    m_GSW   )
      CALL BornV_SetQCDcor( m_QCDcor)
      Sig0nb   = BornV_Sig0nb(m_CMSene)
      cmax=  0.99999999999d0
      CALL TabCheck_Gaus8(TabCheck_DTheDiz, -cmax,  0d0, xsB) !!! 8-point gauss
      CALL TabCheck_Gaus8(TabCheck_DTheDiz,   0d0, cmax, xsF) !!! 8-point gauss
      Xsec=Sig0nb*(xsB+xsF)
      Afb =0d0
      IF(Xsec.NE.0d0) Afb = (xsF-xsB)/(xsF+xsB)
      END

      DOUBLE PRECISION FUNCTION TabCheck_DTheTab(CosThe)
*/////////////////////////////////////////////////////////////////////////////////
*//                                                                             //
*/////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'TabCheck.h'
      DOUBLE PRECISION   CosThe,svar
      DOUBLE PRECISION   BornV_Differential
*     ---------------------------------------------------------
      svar     = m_CMSene**2
      TabCheck_DTheTab = 
     $    3d0/8d0* BornV_Differential( 1,m_KFfin,svar,CosThe,m_seps1,m_seps2,m_ta,m_tb)
      END


      DOUBLE PRECISION FUNCTION TabCheck_DTheDiz(CosThe)
*/////////////////////////////////////////////////////////////////////////////////
*//                                                                             //
*/////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'TabCheck.h'
      DOUBLE PRECISION   CosThe,svar
      DOUBLE PRECISION   BornV_Differential
*     ---------------------------------------------------------
*     we get GSW formfactors from Dizet and pump them into BornV
      IF(  m_CMSene .GT. m_WmaxLEP1) THEN     ! costheta dependence above 120GeV only
         CALL DZface_MakeGSW( 0, m_CMSene, CosThe, m_GSW, m_QCDcor)
         CALL BornV_SetGSW(    m_GSW   )
         CALL BornV_SetQCDcor( m_QCDcor)
      ENDIF
      svar     = m_CMSene**2
      TabCheck_DTheDiz = 
     $    3d0/8d0* BornV_Differential( 3,m_KFfin,svar,CosThe,m_seps1,m_seps2,m_ta,m_tb)    
      END

      SUBROUTINE TabCheck_Gaus8(fun,aa,bb,result)
*//////////////////////////////////////////////////////////////////////////////
*//   8-point Gauss                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      DOUBLE PRECISION  fun,aa,bb,result
      EXTERNAL fun
      DOUBLE PRECISION  a,b,sum8,xmidle,range,xplus,xminu
      INTEGER           k,i
      DOUBLE PRECISION  wg(4),xx(4)
      DATA wg /0.101228536290376d0, 0.222381034453374d0, 0.313706645877887d0, 0.362683783378362d0/
      DATA xx /0.960289856497536d0, 0.796666477413627d0, 0.525532409916329d0, 0.183434642495650d0/
*-------------------------------------------------------------------------------
      a  = aa
      b  = bb
      xmidle= 0.5d0*(a+b)
      range = 0.5d0*(b-a)
      sum8 =0d0
      DO i=1,4
         xplus= xmidle+range*xx(i)
         xminu= xmidle-range*xx(i)
         sum8 =sum8  +(fun(xplus)+fun(xminu))*wg(i)/2d0
      ENDDO
      result = sum8*(b-a)
      END

      SUBROUTINE TabCheck_Gaus16(fun,aa,bb,result)
*//////////////////////////////////////////////////////////////////////////////
*//   12-point Gauss                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      DOUBLE PRECISION  fun,aa,bb,result
      EXTERNAL fun
      DOUBLE PRECISION  a,b,sum16,xmidle,range,xplus,xminu
      INTEGER           k,i
      DOUBLE PRECISION  wg(8),xx(8)
      DATA wg              /0.027152459411754d0, 0.062253523938648d0,
     $ 0.095158511682493d0, 0.124628971255534d0, 0.149595988816577d0,
     $ 0.169156519395003d0, 0.182603415044924d0, 0.189450610455069d0/
      DATA xx              /0.989400934991650d0, 0.944575023073233d0,
     $ 0.865631202387832d0, 0.755404408355003d0, 0.617876244402644d0,
     $ 0.458016777657227d0, 0.281603550779259d0, 0.095012509837637d0/
*-------------------------------------------------------------------------------
      a  = aa
      b  = bb
      xmidle= 0.5d0*(a+b)
      range = 0.5d0*(b-a)
      sum16 =0d0
      DO i=1,8
         xplus= xmidle+range*xx(i)
         xminu= xmidle-range*xx(i)
         sum16 =sum16  +(fun(xplus)+fun(xminu))*wg(i)/2d0
      ENDDO
      result = sum16*(b-a)
      END
