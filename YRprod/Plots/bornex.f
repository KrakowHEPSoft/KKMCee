*----------------------------------------------------------------
*     gmake bornex
*----------------------------------------------------------------
      PROGRAM MAIN
*     ************
      IMPLICIT NONE
      INTEGER   ninp,nout
      CHARACTER*60      Tesnam
*--------------------------------------------------------------------
      ninp=  5
      nout= 16
      Tesnam    = 'bornex'
      OPEN( nout, file='output-'//Tesnam)
      CALL GLK_SetNout(nout)

      CALL TestBorSimp          !  Calibrating KKsem wwith simplistic Born
***      CALL TestBorDres         ! Comparison for Dressed Born with EW corrections
***      CALL TestQCDCOR          ! Test of FSR QCD factor QCDCOR

      CLOSE(nout)
      END

      SUBROUTINE  TestBorSimp
*///////////////////////////////////////////////////////////////////////////////////////////////
*//   Calibrating KKsem wwith simplistic Born
*//   make bornex
*///////////////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
*------------------------------------------------------------------------------
      INTEGER    imax
      PARAMETER (imax=10000)
      DOUBLE PRECISION     xpar(imax)
      DOUBLE PRECISION  BornPb,Born,CMSene,svar
      INTEGER           ninp,nout,KeyFoB,iEne,iKF,KF
      CHARACTER*60      Dname
*------------------------------------------------------------------
      ninp=  5
      nout= 16
      Dname  = './bornex-Mu.input'
      CALL KK2f_ReaDataX('../../.KK2f_defaults', 1,imax,xpar)  ! Read data, the same as in MC run
      CALL KK2f_ReaDataX(                 Dname, 0,imax,xpar)  ! Read user input
      CALL KK2f_Initialize( xpar)                    ! Initialize generator with the production data
      CALL KKsem_Initialize(xpar)                    ! Initialize semianalytical package
*=================================================================================================
      WRITE(   *,*) '*************************************************************************'!
      WRITE(   *,*) '***************** Bornex:: TestBorDres **********************************'!
      WRITE(   *,*) '*************************************************************************'!
      CMSene=xpar(1)
      CALL KKsem_GetBorn(Born)
      BornPb = Born*1000d0
      WRITE(   *,'(a,f10.7,a)') ' Born [Pb]= ',BornPb, '  BornV_Dizet Gauss'!
      WRITE(nout,'(a,f10.7,a)') ' Born [Pb]= ',BornPb, '  BornV_Dizet Gauss'!
*------------
      KeyFoB= 10
      CALL KKsem_SetKeyFoB(KeyFoB)
      CALL KKsem_GetBorn(Born)
      BornPb = Born*1000d0
      WRITE(   *,'(a,f10.7,a,i5,a)') ' Born [Pb]= ',BornPb,' KeyFoB=',KeyFoB, '  BornV_Dizet, costhe=0'!
*------------
      KeyFoB=-11
      CALL KKsem_SetKeyFoB(KeyFoB)
      CALL KKsem_GetBorn(Born)
      BornPb = Born*1000d0
      WRITE(   *,'(a,f10.7,a,i5,a)') ' Born [Pb]= ',BornPb,' KeyFoB=',KeyFoB, '  BornV_Simple'!
*------------
      KeyFoB=-10
      CALL KKsem_SetKeyFoB(KeyFoB)
      CALL KKsem_GetBorn(Born)
      WRITE(   *,'(a,f10.7,a,i5,a)') ' Born [Pb]= ',BornPb,' KeyFoB=',KeyFoB, '  KKsem_BornV'!
*------------------------------------------------------------------------------------------------
      WRITE(*,*) '******************************************************************************'!
      WRITE(*,*) '         Calibrating KKsem wwith simplistic Born                              '!
      WRITE(*,*) '         Compare with zf_2f_ForSJ_primBORN_ZG.log.dima                        '!
      WRITE(*,*) '******************************************************************************'!
      KeyFoB=-10
      CALL KKsem_SetKeyFoB(KeyFoB)
      CALL BornV_SetKeyQCD(0)
      DO iKF=0,5
         KF=iKF
         IF(iKF.EQ.0)  KF=13
         DO iEne =1,3
            IF(iEne.EQ.1)  CMSene = 189d0
            IF(iEne.EQ.2)  CMSene = 200d0
            IF(iEne.EQ.3)  CMSene = 206d0
            svar = CMSene**2
            CALL BornV_SetKF(KF) ! set just one fermion
            CALL KKsem_MakeBorn(svar,Born)
            BornPb = Born*1000d0
            WRITE(   *,'(a,i4,a,f7.1,a,f13.6)') 'KF=',KF,' CMSene=', CMSene ,'  Born[Pb]= ',BornPb !
         ENDDO
      ENDDO
      WRITE(*,*) '******************************************************************************'!
      END


      SUBROUTINE  TestBorDres
*///////////////////////////////////////////////////////////////////////////////////////////////
*     make bornex
*     ./txp_to_eps TabBornex1
*     ./txp_to_eps TabBornex2
*     ./txp_to_eps TabBornex3
*///////////////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      DOUBLE PRECISION  BornPb,QCDcor,Born
      DOUBLE PRECISION  QCDcorR(10)
      DOUBLE PRECISION  rEne,svar,MZ
      INTEGER           ninp,nout
      CHARACTER*60      Dname
*------------------------------------------------------------------------------
      INTEGER    imax
      PARAMETER (imax=10000)
      DOUBLE PRECISION     xpar(imax)
      DOUBLE PRECISION     ymin,ymax
      INTEGER              KFferm,i,j,k,KeyFob
      INTEGER              idDown,      idUp,       idStran,      idCharm,      idBott !
      PARAMETER(           idDown=1001, idUp=1002,  idStran=1003, idCharm=1004, idBott=1005 )!
      INTEGER              jdDown,      jdUp,       jdStran,      jdCharm,      jdBott !
      PARAMETER(           jdDown=2001, jdUp=2002,  jdStran=2003, jdCharm=2004, jdBott=2005 )!
      INTEGER              kdDown,      kdUp,       kdStran,      kdCharm,      kdBott !
      PARAMETER(           kdDown=3001, kdUp=3002,  kdStran=3003, kdCharm=3004, kdBott=3005 )!
*------------------------------------------------------------------
! Parameters for tables
      CHARACTER*60  TeXfile
      INTEGER       nCols
      INTEGER       idl(6)
      CHARACTER*80  capt(7)
      CHARACTER*16  fmt(3), fmtx,fmty
      CHARACTER*80  mcapt
*------------------------------------------------------------------
      ninp=  5
      nout= 16
      Dname  = './bornex-QCD.input'
      CALL KK2f_ReaDataX('../../.KK2f_defaults', 1,imax,xpar)  ! Read data, the same as in MC run
      CALL KK2f_ReaDataX(                 Dname, 0,imax,xpar)  ! Read user input
      CALL KK2f_Initialize( xpar)                    ! Initialize generator with the production data
      CALL KKsem_Initialize(xpar)                    ! Initialize semianalytical package
*---------
      WRITE(   *,*) '*************************************************************************'!
      WRITE(nout,*) '*************************************************************************'!
*=================================================================================================
      WRITE(   *,*) '*************************************************************************'!
      WRITE(   *,*) '***************** Bornex:: TestBorDres **********************************'!
      WRITE(   *,*) '*************************************************************************'!
*--------------------------------------------------------------------------
      CALL  DefEne(1000)
      KFferm = 1
      CALL DefXsec(KFferm,idDown)
      CALL RDatZF2('zf_2f_ForSJ_primBORN.log.dima',KFferm,jdDown)
      KFferm = 2
      CALL DefXsec(KFferm,idUp)
      CALL RDatZF2('zf_2f_ForSJ_primBORN.log.dima',KFferm,jdUp)
      KFferm = 3
      CALL DefXsec(KFferm,idStran)
      CALL RDatZF2('zf_2f_ForSJ_primBORN.log.dima',KFferm,jdStran)
      KFferm = 4
      CALL DefXsec(KFferm,idCharm)
      CALL RDatZF2('zf_2f_ForSJ_primBORN.log.dima',KFferm,jdCharm)
      KFferm = 5
      CALL DefXsec(KFferm,idBott)
      CALL RDatZF2('zf_2f_ForSJ_primBORN.log.dima',KFferm,jdBott)
*
      CALL GLK_Operat(jdDown,  '-', idDown,  kdDown,  1d0, 1d0) !
      CALL GLK_Operat(kdDown , '/', idDown,  kdDown,  1d0, 1d0) !
      CALL GLK_Operat(jdUp,    '-', idUp,    kdUp,    1d0, 1d0) !
      CALL GLK_Operat(kdUp ,   '/', idUp,    kdUp,    1d0, 1d0) !
      CALL GLK_Operat(jdStran, '-', idStran, kdStran, 1d0, 1d0) !
      CALL GLK_Operat(kdStran, '/', idStran, kdStran, 1d0, 1d0) !
      CALL GLK_Operat(jdCharm, '-', idCharm, kdCharm, 1d0, 1d0) !
      CALL GLK_Operat(kdCharm, '/', idCharm, kdCharm, 1d0, 1d0) !
      CALL GLK_Operat(jdBott,  '-', idBott,  kdBott,  1d0, 1d0) !
      CALL GLK_Operat(kdBott , '/', idBott,  kdBott,  1d0, 1d0) !
*===================================================================
      fmt(1)='F10.0'
      fmt(2)='F10.4'
      fmt(3)= 'F8.4'
      capt(1)='No.'
      capt(2)='{\\color{blue}$\\sqrt{s}$}'
      capt(3)='{\\color{blue} d }'
      capt(4)='{\\color{blue} u }'!
      capt(5)='{\\color{blue} s }'!
      capt(6)='{\\color{blue} c }'!
      capt(7)='{\\color{blue} b }'!
      nCols  = 6
      idl(1) = 1000
      idl(2) = idDown
      idl(3) = idUP
      idl(4) = idStran
      idl(5) = idCharm
      idl(6) = idBott
*-----------
      TeXfile   = 'TabBornex1.txp'
      CALL GLK_PlInitialize(2,TeXfile) !Initialize GLK_Plot
      Mcapt = 'Born KKsem LEP2'
      CALL GLK_SetTabRan(2,19,1)
      CALL GLK_PlTable2(nCols,idl,capt,Mcapt,fmt,' ','R',' ')
      Mcapt = 'Born KKsem LEP1'
      CALL GLK_SetTabRan(20,40,2)
      CALL GLK_PlTable2(nCols,idl,capt,Mcapt,fmt,'S','R',' ')
      Mcapt = 'Born KKsem Low energy'
      CALL GLK_SetTabRan(41,44,1)
      CALL GLK_PlTable2(nCols,idl,capt,Mcapt,fmt,'S','R',' ')
      CALL GLK_PlEnd                      !Close GLK_Plot and close its file
*-----------
      idl(2) = jdDown
      idl(3) = jdUP
      idl(4) = jdStran
      idl(5) = jdCharm
      idl(6) = jdBott
      TeXfile   = 'TabBornex2.txp'
      CALL GLK_PlInitialize(2,TeXfile) !Initialize GLK_Plot
       Mcapt = 'Born ZFITTER LEP2'
      CALL GLK_SetTabRan(2,19,1)
      CALL GLK_PlTable2(nCols,idl,capt,Mcapt,fmt,' ','R',' ')
      Mcapt = 'Born ZFITTER LEP1'
      CALL GLK_SetTabRan(20,40,2)
      CALL GLK_PlTable2(nCols,idl,capt,Mcapt,fmt,'S','R',' ')
      Mcapt = 'Born ZFITTER Low energy'
      CALL GLK_SetTabRan(41,44,1)
      CALL GLK_PlTable2(nCols,idl,capt,Mcapt,fmt,'S','R',' ')
      CALL GLK_PlEnd                      !Close GLK_Plot and close its file
*-----------
      idl(2) = kdDown
      idl(3) = kdUP
      idl(4) = kdStran
      idl(5) = kdCharm
      idl(6) = kdBott
      TeXfile   = 'TabBornex3.txp'
      CALL GLK_PlInitialize(2,TeXfile) !Initialize GLK_Plot
       Mcapt = 'Born (ZF-KKsem)/KKsem LEP2'
      CALL GLK_SetTabRan(2,19,1)
      CALL GLK_PlTable2(nCols,idl,capt,Mcapt,fmt,' ','R',' ')
      Mcapt = 'Born (ZF-KKsem)/KKsem LEP1'
      CALL GLK_SetTabRan(20,40,2)
      CALL GLK_PlTable2(nCols,idl,capt,Mcapt,fmt,'S','R',' ')
      Mcapt = 'Born (ZF-KKsem)/KKsem Low energy'
      CALL GLK_SetTabRan(41,44,1)
      CALL GLK_PlTable2(nCols,idl,capt,Mcapt,fmt,'S','R',' ')
      CALL GLK_PlEnd                      !Close GLK_Plot and close its file
*--------------------------------------------------------------------------
      WRITE(nout,*) '*************************************************************************'!
      END

      SUBROUTINE  DefEne(idHist)
*********************************************************
      IMPLICIT NONE
      INTEGER              idHist
      DOUBLE PRECISION     rEne, MZ
      INTEGER              nEntry
      PARAMETER(           nEntry=44 )
      DOUBLE PRECISION     tabEne(nEntry)
      INTEGER              icont
*
      iCont=0
      DO rEne = 188d0,206d0,1d0
         iCont=iCont+1
         tabEne(iCont)=rEne
      ENDDO
      MZ  = 91.187D0
      DO rEne = -5d0,5d0, 0.5d0
         iCont=iCont+1
         tabEne(iCont)=MZ+rEne
      ENDDO
      DO rEne = 16d0,29d0,4d0
         iCont=iCont+1
         tabEne(iCont)=rEne
      ENDDO
      CALL GLK_Book1(idHist,'Ene  $', nEntry, 0d0,1d0*nEntry) !
      CALL GLK_Pak(idHist,tabEne)
      END


      SUBROUTINE  DefXsec(KFferm,idHist)
*********************************************************
      IMPLICIT NONE
      INTEGER              KFferm, idHist
      INTEGER              icont
      INTEGER              i,j,k,KeyFoB,nout
      DOUBLE PRECISION     rEne,svar,Born,MZ
      INTEGER              nEntry
      PARAMETER(           nEntry=44 )
      DOUBLE PRECISION     tabEne(nEntry),tabXec(nEntry)
*----------
      CALL GLK_UnPak(1000,tabEne)
      CALL GLK_Book1(idHist,'Sigma  $', nEntry, 0d0,1d0*nEntry) !
      nout =16
      CALL BornV_SetKF(KFferm)  ! set just one fermion
      CALL BornV_SetKeyQCD(0)
      KeyFoB= 0 ! Normal
      CALL KKsem_SetKeyFoB(KeyFoB)
      WRITE(   *,*) '---------------------------------------------------------------------------'!
      DO i=1,nEntry
         svar=tabEne(i)**2
         CALL KKsem_MakeBorn(svar,Born)
         tabXec(i)=Born*1d3
         WRITE(   *,'(a,i5,f7.1,a,f13.6)') 'i,CMSene=',i, tabEne(i) ,'  Born[Pb]= ',tabXec(i)!
         WRITE(nout,'(a,i5,f7.1,a,f13.6)') 'i,CMSene=',i, tabEne(i) ,'  Born[Pb]= ',tabXec(i)!
      ENDDO
      CALL GLK_Pak(idHist,tabXec)
      END

      SUBROUTINE  TestQCDCOR
*///////////////////////////////////////////////////////////////////////////////////////////////
*///////////////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      DOUBLE PRECISION  BornPb,QCDcor,Born
      DOUBLE PRECISION  QCDcorR(10)
      INTEGER   ninp,nout
      CHARACTER*60      Dname
*------------------------------------------------------------------------------
      INTEGER    imax,i,j,k
      PARAMETER (imax=10000)
      REAL*8     xpar(imax)
      REAL*8     ymin,ymax
*---------------------------------------
      ninp=  5
      nout= 16
      Dname  = './bornex-QCD.input'
      CALL KK2f_ReaDataX('../../.KK2f_defaults', 1,imax,xpar)  ! Read data, the same as in MC run
      CALL KK2f_ReaDataX(                 Dname, 0,imax,xpar)  ! Read user input
      CALL KK2f_Initialize( xpar)                    ! Initialize generator with the production data
      CALL KKsem_Initialize(xpar)                    ! Initialize semianalytical package
*=================================================================================================
      WRITE(   *,*) '*************************************************************************'!
      WRITE(   *,*) '************************** TestQCDCOR************************************'!
      WRITE(   *,*) '*************************************************************************'!
      CALL KKsem_GetBorn(Born)  ! This re-sets QCDcor
      BornPb = Born*1000d0
      WRITE(   *,'(a,f10.7)') 'Bornex:: Born [Pb], QCD on = ',BornPb
      WRITE(nout,'(a,f10.7)') 'Bornex:: Born [Pb], QCD on = ',BornPb
      CALL BornV_GetQCDcor(QCDcor)
      WRITE(   *,'(a,f10.7)') 'Bornex:: QCDcor= ',QCDcor
      WRITE(nout,'(a,f10.7)') 'Bornex:: QCDcor= ',QCDcor
      CALL BornV_SetKeyQCD(0)
      CALL KKsem_GetBorn(Born)
      BornPb = Born*1000d0
      WRITE(   *,'(a,f10.7)') 'Bornex:: Born [Pb], QCD off= ',BornPb
      WRITE(nout,'(a,f10.7)') 'Bornex:: Born [Pb], QCD off= ',BornPb
      WRITE(   *,*) '-------------------------------------------------------------' !
      CALL KKsem_GetBorn(Born)  ! This re-sets QCDcor
      CALL BornV_GetQCDcorR(QCDcorR)
      WRITE(   *,'(a,9f10.7)') 'Bornex:: QCDcorR= ',(QCDcorR(i),i=1,4)
      WRITE(nout,'(a,9f10.7)') 'Bornex:: QCDcorR= ',(QCDcorR(i),i=1,4)
      WRITE(   *,*) '-------------------------------------------------------------' !
      END

      SUBROUTINE   RDatZF2(DataFile1,iKFf,idHist)
*/////////////////////////////////////////////////////////////////////////////
*// reading sigma and  AFB with error from diskfile, KORALZ results         //
*/////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      CHARACTER*(*)     DataFile1
      INTEGER           ikFf,KF,idHist
      CHARACTER*4       cSAB
      CHARACTER*33      ch2
      DOUBLE PRECISION     ene,vmax,xs3,as3
      INTEGER              iSel,line,indf,iEne,i
      INTEGER              nEntry
      PARAMETER(           nEntry=44 )
      DOUBLE PRECISION     tabEne(nEntry),tabXec(nEntry)
* Translation matrix from INDF to KF
      INTEGER     m_KFfromINDF(0:11)
*                         nu,  e, mu, tau,  u,  d,  c,  s,  t,  b, q's,Bhabha 
      DATA m_KFfromINDF / -1, -1, 13,  15,  2,  1,  4,  3, -1,  5,   7,    -1/!
      SAVE m_KFfromINDF
*-------------
      OPEN(20,File=DataFile1)
* read header
      WRITE(*,*) '========== SAB_RDatZF2 Enter ================'
      DO line=1,10000
         READ(20,'(a,a,I2)') cSAB
         IF(cSAB .EQ. '<be>' ) GOTO 200
      ENDDO
      WRITE(*,*) ' STOP in SAB_RDatZF1, end of data file',DataFile1
      STOP
 200  CONTINUE

      CALL GLK_UnPak(1000,tabEne)
* read lines
      DO line=1,10000
*********************************
         READ(20,'(   a,   a,   I2,  f7.1, f8.3, 9f12.4)') 
     $             cSAB, ch2, indf,   ene, vmax, xs3!
c         WRITE(*,'(   a,   a,   I2,  f7.1, f8.3, 9f12.4)') 
c     $             cSAB, ch2, indf,   ene, vmax, xs3!
*********************************
         IF(cSAB .EQ. '<en>') GOTO 900
         KF= m_KFfromINDF(indf)
         IF( KF .EQ. iKFf ) THEN
            iSel= 0
            DO iEne=1,nEntry
               IF( ABS(Ene-tabEne(iEne)) .LT.0.1) iSel=iEne
            ENDDO
            tabXec(iSel)=xs3
         ENDIF
      ENDDO
 900  CONTINUE
      CLOSE(20)
      WRITE(   *,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'!
      DO i=1,nEntry
         WRITE(   *,'(a,i5,f7.1,a,f13.6)') 'i,CMSene=',i, tabEne(i) ,'  Born[Pb]= ',tabXec(i)!
      ENDDO
      WRITE(*,*) '========== SAB_RDatZF2 Exit  ================'
      CALL GLK_Book1(idHist,'Ene  $', nEntry, 0d0,1d0*nEntry) !
      CALL GLK_Pak(  idHist, tabXec )
      CALL GLK_Print(idHist)
      END
