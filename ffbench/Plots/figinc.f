*------------------------------------------------------
*     make figinc-dvi
*     make figinc-ps
*------------------------------------------------------
      PROGRAM MAIN
*     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / inout  / ninp,nout
*
      ninp=  5
      nout= 16
      CALL GLK_SetNout(nout)
*=====================================================
* Initial + Final state plots
      CALL ploter1
      END

      SUBROUTINE ploter1
*     ******************
      COMMON / INOUT  / NINP,NOUT
      SAVE
      CHARACTER*60  Tesnam, Dname
      CHARACTER*60  Hname, Hname2, DumpFile
*----------------------------------------------
      Tesnam    = 'figinc'
      OPEN( nout, file='output-'//Tesnam)

* PRINCIPAL, without hadronization
c      Dname='../Inclusive/Inclusive.input.10k'! jan 99
c      Hname='../Inclusive/pro.hst.10k'        ! jan 99
      Dname  = '../Inclusive/Inclusive.input' ! current
      Hname  = '../Inclusive/pro.hst'         ! current
* Beamsstrahlung
c      Dname  = '../Beast/Beast.input.BstON'    ! beamstrahlung
c      Hname  = '../Beast/pro.hst.BstOFF20k'    ! beamstrahlung OFF
c      Hname2 = '../Beast/pro.hst.BstON20k'     ! beamstrahlung ON

* Read data, the same as in MC run
      CALL ReaDat(Dname)
* Read histograms from MC run
      CALL GLK_ReadFile(Hname)
*      CALL GLK_ReadFile(Hname2)
      CALL GLK_ListPrint(6)     ! debug
*=========================================================================
**** This is added to fool linker!!!
      CALL KK2f_GetKeyFSR(KeyFSR)
*=========================================================================
*
      CALL figchi
c      CALL figbeast
*
*=========================================================================
* Write all histograms into dump file, for control
      DumpFile = './dump.hst'
      CALL GLK_WriteFile(DumpFile)
*=========================================================================
      CLOSE(nout)
      END


      SUBROUTINE figchi
*     ************************************************************
* MC<-->SAN, Total O(alf2),  O(alf2), linear scale, sigma(vmax)
*     ************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / inout  / ninp,nout
      PARAMETER (imax=10000)
      COMMON  /xpar/  xpar(imax)
      SAVE
*------------------------------------------------------------------
      CHARACTER*80 title
*------------------------------------------------------------------
! Parameters for tables
      DIMENSION    idl(5)
      CHARACTER*32 capt(6)
      CHARACTER*16  fmt(3), fmtx,fmty
      CHARACTER*80 mcapt
*------------------------------------------------------------------
      CHARACTER*80 CaptV(8)
      DATA CaptV /
     $'\\large{',
     $'Final fermion effective mass distribution from {\\sl KK} MC at 500GeV. ',
     $'Sharp Z ratiative return peak is seen.',
     $'Included are muons, taus and all quarks except top.',
     $'Properly normalized in nanobs.',
     $'}',
     $'\\label{fig:Figs-V}',
     $'% end-of-caption'/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
      CHARACTER*80 CaptV2(6)
      DATA CaptV2 /
     $'\\large{',
     $'Dependence of cross section on cut $v_{\\max}$, obtained from {\\sl KK} MC.',
     $'It is divided by the Born cross-section.',
     $'}',
     $'\\label{fig:Figs-V2}',
     $'% end-of-caption'/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
      CHARACTER*80 CaptKF(6)
      DATA CaptKF /
     $'\\large{',
     $'This is KF code distribution from {\\sl KK} MC at 500GeV.',
     $'Properly normalized in nanobs.',
     $'}',
     $'\\label{fig:Figs-KF}',
     $'% end-of-caption'/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
      CHARACTER*80 LabelV(6)
      DATA LabelV /
     $'\\put(300,250){\\begin{picture}( 1200,1200) ',
     $'\\put( 60, 1070){\\makebox(0,0)[l]{\\Huge ',
     $'                     $ {d\\sigma \\over d v}  $ }}',
     $'\\put( 500,  40){\\makebox(0,0)[b]{\\Large $v = 1-s^\\prime /s$ }}',
     $'\\end{picture}}',
     $'% end-of-label'/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
      CHARACTER*80 LabelV2(6)
      DATA LabelV2 /
     $'\\put(300,250){\\begin{picture}( 1200,1200) ',
     $'\\put( 60, 1070){\\makebox(0,0)[l]{\\Huge ',
     $' $ {\\sigma \\over \\sigma_{Born} }  $ }}',
     $'\\put( 500,  40){\\makebox(0,0)[b]{\\Large $v_{\\max}$ }}',
     $'\\end{picture}}',
     $'% end-of-label'/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
      CHARACTER*80 LabelKF(5)
      DATA LabelKF /
     $'\\put(300,250){\\begin{picture}( 1200,1200) ',
     $'\\put(-50, 1000){\\makebox(0,0)[r]{\\LARGE  $\\sigma_{KF}$ }}',
     $'\\put( 500,  40){\\makebox(0,0)[b]{\\Large  KF }}',
     $'\\end{picture}}',
     $'% end-of-label'/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
! Mark plots for plots
      CHARACTER*32 star,diamond,circle,times,disc,dot
      PARAMETER (diamond ='\\makebox(0,0){\\LARGE $\\diamond$}')
      PARAMETER (star    ='\\makebox(0,0){\\LARGE $\\star$}')
      PARAMETER (circle  ='\\circle{30}')
      PARAMETER (times   ='\\makebox(0,0){\\LARGE $\\times$}')
      PARAMETER (disc    ='\\circle*{20}')
      PARAMETER (dot     ='\\circle*{10}')
!-----------------------------------
      CHARACTER*80 captab(4)
      DATA captab /
     $'Total cross sections for various final fermions.',
     $'Center of mass energy $\\protect\\sqrt{s}=500 GeV$.',
     $'Not calculated cross sections are set to zero.',
     $'% end-of-caption'/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
      CHARACTER*60  TeXfile

      WRITE(NOUT,*) ' ====================================='
      WRITE(NOUT,*) ' ==========    figchi   =============='
      WRITE(NOUT,*) ' ====================================='

      CMSene=xpar(1)
      Sig0nb = BornV_Sig0nb(CMSene)
      BornR  = BornV_Integrated(0,CMSene**2)
      BornNb = BornR*Sig0nb
      BornPb = BornNb*1000
      WRITE(nout,*) 'Figchi:: Born [R]=  ',BornR
      WRITE(nout,*) 'Figchi:: Born [pb]= ',BornPb
      WRITE(*,   *) 'Figchi:: Born [R]=  ',BornR
      WRITE(*,   *) 'Figchi:: Born [pb]= ',BornPb

      idv =  50000
      kdv = 150000
      jdv = 250000

      CALL GLK_hinbo1(idv+1,title,nbiv,vmin,vmax)

      CALL GLK_Print(idv+1)
      CALL GLK_Print(idv+2)
      CALL GLK_Print(idv+3)
      CALL GLK_Print(idv+4)
      CALL GLK_Print(idv+7)
      CALL GLK_Print(idv+9)

* renormalize properly histograms
      IdGenYFS3 = 6
      CALL GLK_RenHst("NB  ",IdGenYFS3,idv+1,kdv+1)
      CALL GLK_RenHst("NB  ",IdGenYFS3,idv+9,kdv+9)
      CALL GLK_CumHis(       IdGenYFS3,idv+1,jdv+1)
      CALL GLK_CumHis(       IdGenYFS3,idv+9,jdv+9)
      CALL GLK_Operat(jdv+1,'+',jdv+1,jdv+1,0d0,1d0/BornNb)
      CALL GLK_idopt(kdv+1,'ERRO')
      CALL GLK_idopt(kdv+9,'ERRO')
      CALL GLK_idopt(jdv+1,'ERRO')
      CALL GLK_idopt(jdv+9,'ERRO')
      CALL GLK_Yminim(jdv+1, 0d0)

      CALL GLK_Print(kdv+1)
      CALL GLK_Print(kdv+9)
      CALL GLK_Print(jdv+1)
      CALL GLK_Print(jdv+9)


*-------------------------------------------
* Initialize GLK_Plot
      TeXfile   = 'figinc.tex'
      Lint =0
      CALL GLK_PlInitialize(Lint,TeXfile)
*===========================================
*          Plot of d_sigma/d_v
*===========================================
      fmtx='f10.1'
      fmty='f10.2'
* Normal
      CALL GLK_plcapt(CaptV)
      CALL GLK_plot2(kdv+1,' ',' ',dot    ,fmtx,fmty)
      CALL GLK_PlLabel(LabelV)
* Cumulative
      CALL GLK_plcapt(CaptV2)
      CALL GLK_plot2( jdv+1,' ',' ',dot    ,fmtx,fmty)
      CALL GLK_PlLabel(LabelV2)
*===========================================
*          Plot of sigma_KF
*===========================================
      fmtx='f10.0'
      fmty='f10.3'
      CALL GLK_plcapt(CaptKF)
      CALL GLK_plot2(kdv+9,' ',' ',dot    ,fmtx,fmty)
      CALL GLK_PlLabel(LabelKF)
*===================================================================
*                Tables
*===================================================================
! BARE1 WW trigger O(alf2)exp
      capt(1)='KF'
      capt(2)='$\\sigma_{KF}$'
      capt(3)='$\\sigma_{KF}^C$'
      fmt(1)='F10.0'
      fmt(2)='F10.6'
      fmt(3)='F8.6'
      idl(1)=kdv+9
      idl(2)=jdv+9
!------------------------------------------------
      CALL GLK_PlCapt(CapTab)
      Mcapt = ' {\\sl KK} at 500GeV'
      CALL GLK_PlTable2(2,idl,capt,Mcapt,fmt,' ','L',' ')
*===================================================================
Close GLK_Plot and close its file
      CALL GLK_PlEnd
      END

      SUBROUTINE figbeast
*     ************************************************************
* MC<-->SAN, Total O(alf2),  O(alf2), linear scale, sigma(vmax)
*     ************************************************************
      IMPLICIT NONE
      INTEGER           ninp,nout
      COMMON / inout  / ninp,nout
      INTEGER    imax
      PARAMETER (imax=10000)
      DOUBLE PRECISION xpar
      COMMON  /xpar/   xpar(imax)
*------------------------------------------------------------------
      CHARACTER*80  title
      CHARACTER*16  fmt(3), fmtx,fmty
      INTEGER       lint, IdGenYFS3
      INTEGER       idv,kdv,jdv, idv2
      CHARACTER*60  TeXfile
*------------------------------------------------------------------
! Mark plots for plots
      CHARACTER*32 star,diamond,circle,times,disc,dot
      PARAMETER (diamond ='\\makebox(0,0){\\LARGE $\\diamond$}')
      PARAMETER (star    ='\\makebox(0,0){\\LARGE $\\star$}')
      PARAMETER (circle  ='\\circle{30}')
      PARAMETER (times   ='\\makebox(0,0){\\LARGE $\\times$}')
      PARAMETER (disc    ='\\circle*{20}')
      PARAMETER (dot     ='\\circle*{10}')
*------------------------------------------------------------------
*    $_________|_________|_________|_________|_________|_________|_________|_________|
      CHARACTER*80 CaptE(7)
      DATA CaptE /
     $'\\large{',
     $'Distribution of energy loss due to beamstrahlung in units of beam',
     $'From KK run at 500GeV. Beamstralung spectrum according to CIRCE ',
     $'of T. Ohl, TESLA option.  Total sample 20k events. ',
     $'}',
     $'\\label{fig:Figs-E}',
     $'% end-of-caption'/
      CHARACTER*80 LabelE(5)
      DATA LabelE /
*    $_________|_________|_________|_________|_________|_________|_________|_________|
     $'\\put(300,250){\\begin{picture}( 1200,1200) ',
     $'\\put( 60, 1070){\\makebox(0,0)[l]{\\Huge $ {dN  \\over d x}  $ }}',
     $'\\put( 600,-100){\\makebox(0,0)[t]{\\Large $x = log10(E/Ebeam) $ }}',
     $'\\end{picture}}',
     $'% end-of-label'/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
      CHARACTER*80 CaptX(9)
      DATA CaptX /
     $'\\large{',
     $'Mass distribution of the final fermion pair in the Z resonance region',
     $'(radiative return) from {\\sl KK} MC.',
     $'Circles: beamstrahlung OFF.',
     $'Stars: beamstrahlung ON. Total sample 20k events.',
     $'Included are muons, taus and all quarks except top.',
     $'}',
     $'\\label{fig:Figs-X}',
     $'% end-of-caption'/
*    $_________|_________|_________|_________|_________|_________|_________|_________|
      CHARACTER*80 LabelX(5)
      DATA LabelX /
     $'\\put(300,250){\\begin{picture}( 1200,1200) ',
     $'\\put( 60, 1070){\\makebox(0,0)[l]{\\Huge $ {dN  \\over d M}  $ }}',
     $'\\put( 500,  40){\\makebox(0,0)[b]{\\Large M [GeV] }}',
     $'\\end{picture}}',
     $'% end-of-label'/
*------------------------------------------------------------------
      idv  =  50000
      kdv  = 150000
      jdv  = 250000
      idv2 =  50000 +1000000

      CALL GLK_Print(idv2 +5)
      CALL GLK_idopt(idv2 +5,'ERRO')
      CALL GLK_Print(idv+20)
      CALL GLK_idopt(idv+20,'ERRO')
      CALL GLK_Print(idv2+20)
      CALL GLK_idopt(idv2+20,'ERRO')

      IdGenYFS3 = 6
cc      CALL GLK_RenHst("NB  ",IdGenYFS3,idv+5,kdv+5)
cc      CALL GLK_idopt(kdv+5,'ERRO')
cc      CALL GLK_Print(kdv+5)
*-------------------------------------------
* Initialize GLK_Plot
      TeXfile   = 'figbeast.tex'
      Lint =0
      CALL GLK_PlInitialize(Lint,TeXfile)
*===========================================
      fmtx='f10.1'
      fmty='f10.2'
      CALL GLK_plcapt(CaptE)
      CALL GLK_plot2(idv2  +5,' ',' ',dot    ,fmtx,fmty)
      CALL GLK_PlLabel(LabelE)
*===========================================
      CALL GLK_plcapt(CaptX)
      CALL GLK_plot2(idv2 +20,' ','*',star   ,fmtx,fmty)
      CALL GLK_plot2(idv  +20,'S','*',circle ,fmtx,fmty)
      CALL GLK_PlLabel(LabelX)
*===================================================================
Close GLK_Plot and close its file
      CALL GLK_PlEnd

      END
