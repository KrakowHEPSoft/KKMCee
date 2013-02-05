*///////////////////////////////////////////////////////////////////////////////
*//                   main program for MC production                          //
*///////////////////////////////////////////////////////////////////////////////
*//   make Test_SfacISR
*//   make Test_SfacFSR
*//   make Test_Rmat
*//   make Test_UV
*//   make Test_BornXsec
*//   make Test_IsrSingle
*//   make Test_FsrSingle
*//   make Test_NewSingle
*//   make Test_MultiPhot
*//   make Test_Virtual
*//   make Test_Rules
*//   make Test_DsigOverDtau
*//-----------------------------------------------------------------------------
*//   make Tau-start
*//   make Tau-start-debug
*//
*///////////////////////////////////////////////////////////////////////////////
      PROGRAM Main
*     *************
      IMPLICIT NONE
* input/output files
      INTEGER              ninp,nout
      COMMON / c_yfspro  / ninp,nout
      CHARACTER*4 semaph
      INTEGER     ntot2n, ninph, ninp2, ninp3, ijklin, ntotin
*-------------------------------------------------------------------------
* general output for everybody including glibk
      ninp =5
      nout =16  ! general output for everybody including glibk
      OPEN(nout,FILE='./pro.output')
      REWIND(nout)
      CALL GLK_SetNout(nout)
*
      WRITE(nout,*) '   '
      WRITE(nout,*) '=============================================='
      WRITE(nout,*) '==========*********************==============='
      WRITE(nout,*) '==========***    MainPro    ***==============='
      WRITE(nout,*) '==========*********************==============='
      WRITE(nout,*) '=============================================='
      WRITE(nout,*) '   '
*---------------

* READ semaphore flag
      CALL givsem(semaph)

*---------------
      IF(semaph  .EQ.  'STAR') THEN
         WRITE(6,*) ' ------- starting from the scratch ----------'
* READ initial (root) random number seed
         ninp3=3
         OPEN(ninp3,FILE='./iniseed')
         READ(ninp3,'(i10)') ijklin
         READ(ninp3,'(i10)') ntotin
         READ(ninp3,'(i10)') ntot2n
         CALL PseuMar_Initialize(ijklin,ntotin,ntot2n)
      ELSEIF(semaph  .EQ.  'CONT') THEN
         WRITE(6,*) ' ------- restoring from the disk   ----------'
* restore histograms from the disk
         ninph=10
         OPEN(ninph,FILE='./'//'pro.hst')
         CALL GLK_hrfile(ninph,' ',' ')      !transfer FILE number
         CALL GLK_hrin(   0,9999,0)          !READ from the disk
         CALL GLK_hrend(' ')                 !CLOSE FILE
* READ random number seed stored in semaphore FILE
         ninp2=2
         OPEN(ninp2,FILE='./semaphore')
         READ(ninp2,'(a4)') semaph
         READ(ninp2,'(i10)') ijklin
         READ(ninp2,'(i10)') ntotin
         READ(ninp2,'(i10)') ntot2n
         CALL PseuMar_Initialize(ijklin,ntotin,ntot2n)
         CLOSE(ninp2)
      ELSEIF(semaph  .EQ.  'STOP') THEN
         WRITE(6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++'
         WRITE(6,*) '++++ STOP: please change semaph to cont or start !'
         WRITE(6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++'
         STOP
      ELSE
         WRITE(6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++'
         WRITE(6,*) '++++ STOP: wrong key semaph = ', semaph
         WRITE(6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++'
         STOP
      ENDIF
* *******************************************
      CALL splot
* ********************************************
      END


      SUBROUTINE splot
*     ****************
      IMPLICIT NONE
*
      INTEGER              ninp,nout
      COMMON / c_yfspro  / ninp,nout
*
      INTEGER    imax
      PARAMETER (imax = 10000)               ! length ox xpar
      DOUBLE PRECISION        xpar(imax)
*
      CHARACTER*4 semaph,chdum
      INTEGER     Ie1Pri,Ie2Pri
      INTEGER     igroup, loop, nevt, iev, ngroup
      SAVE
*---------------------------------------------------------------------
      WRITE(nout,*) '   '
      WRITE(nout,*) '=============================================='
      WRITE(nout,*) '==========*********************==============='
      WRITE(nout,*) '==========***    MainPro    ***==============='
      WRITE(nout,*) '==========*********************==============='
      WRITE(nout,*) '=============================================='
      WRITE(nout,*) '   '

* Read data for main program
      OPEN( ninp,FILE='./pro.input')
      READ( ninp,'(/,i10)') nevt
      CLOSE(ninp)
      WRITE(   6,*)   nevt,' requested events '
      WRITE(nout,*)   nevt,' requested events '
*
      CALL KK2f_ReaDataX('../../.KK2f_defaults', 1,imax,xpar)  ! reading general defaults
      CALL KK2f_ReaDataX(         './pro.input', 0,imax,xpar)  ! reading user input
      CALL KK2f_Initialize(xpar)
*
      CALL robol1(-1,xpar)
*-------------------------------------------------------!
*                 main MC loop                          !
*-------------------------------------------------------!
      Ie1Pri =xpar( 6)
      Ie2Pri =xpar( 7)
      ngroup = 2000
      iev=0
      DO loop=1,10000000
        DO igroup =1,ngroup
          iev=iev+1
          IF(MOD(iev, ngroup) .EQ. 1) WRITE( 6,*)  'iev= ',iev
          CALL KK2f_Make
*         **************
          IF(iev .LE. Ie2Pri) THEN      ! Control printouts
             CALL lugive("MSTU(11)=16");
             CALL lulist(1)
          ENDIF
          IF(iev .LE. Ie2Pri) THEN
             CALL lugive("MSTU(11)=6");
             CALL lulist(1)
          ENDIF
*         ============================================
*         histograming
          CALL robol1( 0,xpar)
*         ============================================
*         check on requested no. of events
          IF(iev  .EQ.  nevt)     GOTO 300
        ENDDO
        CALL givsem(semaph)     ! check on semaphore flag
        IF(semaph  .EQ.  'STOP') GOTO 300
        CALL  dumpeh(iev)       ! dump partial results after every ngroup
      ENDDO
 300  CONTINUE
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
      WRITE(6,*) ' generation finished '
      CALL KK2f_Finalize
*------------------------------------------------
      CALL robol1(1,xpar)
*
      CALL  dumpeh(iev)
*     +++++++++++++++++
      END

      SUBROUTINE givsem(semaph)
*     ************************
      IMPLICIT NONE
      CHARACTER*4 semaph
      INTEGER     ninp2
* ------------------------------------------------------
* READ semaphore flag
* ------------------------------------------------------
      ninp2=2
      OPEN(ninp2,FILE='./semaphore')
      READ(ninp2,'(a4)') semaph
      CLOSE(ninp2)
      END

      SUBROUTINE dumpeh(nev)
*     ************************
      IMPLICIT NONE
      INTEGER ninp2, nouth, ntotin, ntot2n, ijklin, nev, icy
* ------------------------------------------------------
* WRITE histos on the disk
* ------------------------------------------------------
      nouth=11
      OPEN(nouth,FILE='./pro.hst')
      CALL GLK_hrfile(nouth,' ','n')   !transfer FILE number
      CALL GLK_hrout( 0,icy,' ')       !WRITE on the disk
      CALL GLK_hrend(' ')              !CLOSE FILE
* ------------------------------------------------------
* overWRITE semaphore FILE flag
      ninp2=2
      OPEN(ninp2,FILE='./semaphore')
      WRITE(ninp2,'(a4)') 'CONT'
* append semaphore FILE with new random number seed in
      CALL PseuMar_Out(ijklin,ntotin,ntot2n)
      WRITE(ninp2,'(i10,a)') ijklin, ' = ijklin '
      WRITE(ninp2,'(i10,a)') ntotin, ' = ntotin '
      WRITE(ninp2,'(i10,a)') ntot2n, ' = ntot2n '
      WRITE(ninp2,'(i10,a)') nev,    ' =    nev '
      CLOSE(ninp2)
* ------------------------------------------------------
      END


      SUBROUTINE Robol1(mode,xpar)
*     *********************************
*  histogramming is done here
*     *********************************
      IMPLICIT NONE
      DOUBLE PRECISION xpar(*)
* input/output files
      COMMON / c_yfspro  / ninp,nout
      INTEGER            ninp,nout
*
      INTEGER            mode, nbv, nbl, idv, ntest, itest, kffin
      INTEGER            i,j,k,l
      DOUBLE PRECISION   ss2, ss, wtcrud, vv, enetot, cmsene, wtmain, vmin, vmax
      DOUBLE PRECISION   zmin, wtz, ebeam,  xkf, eneini, enefin
      DOUBLE PRECISION   zenl5, zenl7, zenl1, zenl3, zenl9
      DOUBLE PRECISION   wtexp0, wtexp1, wtexp2
      SAVE
************************************************************************
      DOUBLE PRECISION  p1(4),p2(4),p3(4),p4(4)
      INTEGER           NphAll,NphIni,NphFin
      DOUBLE PRECISION  PhoAll(100,4)
      DOUBLE PRECISION  PhoIni(100,4),PhoFin(100,4),Ph(4)
      DOUBLE PRECISION  Rho
      DOUBLE PRECISION  WtBest, WtSet(1000),WtSet2(1000)
      SAVE
*=======================================================================
      IF(mode .EQ. -1 ) THEN
         CALL STest_Initialize
*     book histograms
         nbv = 400
         nbv = 50
         vmin= 0
         vmax= 1
         idv  = 50000
         CALL GLK_Book1(idv+1,'d_sigma/dv(v)  $',nbv,vmin,vmax)
         CALL GLK_Book1(idv+2,'Ene total ALL  $',nbv,vmin,vmax)
         CALL GLK_Book1(idv+3,'Ene total ISR  $',nbv,vmin,vmax)
         CALL GLK_Book1(idv+4,'Ene total FSR  $',nbv,vmin,vmax)
         CALL GLK_Book1(idv+7,'NphAll         $',16,0d0,16d0)
         CALL GLK_Book1(idv+9,'KFcode         $',16,1d0,17d0)
         CALL GLK_Book1(idv+10,'NphAll crude  $',16,0d0,16d0)
         CALL GLK_Book1(idv+11,'NphISR crude  $',16,0d0,16d0)
         CALL GLK_Book1(idv+12,'NphFSR crude  $',16,0d0,16d0)
         nbl = 25
         zmin = -5d0
         CALL GLK_Book1(idv+14,'log10(ene) 1-st $',nbl,zmin,0d0)
         CALL GLK_Book1(idv+15,'log10(ene) 3-st $',nbl,zmin,0d0)
         CALL GLK_Book1(idv+16,'log10(ene) 5-st $',nbl,zmin,0d0)
         CALL GLK_Book1(idv+17,'log10(ene) 7-nd $',nbl,zmin,0d0)
         CALL GLK_Book1(idv+18,'log10(ene) 9-nd $',nbl,zmin,0d0)
*
         DO k=201,203
            CALL GLK_Book1(idv+1000+k,'wt CEEX $',25,0d0, 5d0)
            CALL GLK_Book1(idv+2000+k,'wt CEEX $',25,0d0,50d0)
         ENDDO
         DO k=71,74
            CALL GLK_Book1(idv+1000+k,'wt EEX $',25,0d0, 5d0)
         ENDDO
*=======================================================================
      ELSEIF(mode .EQ. 0 ) THEN
* Histograming
         CALL KK2f_GetWtAll(WtMain,WtCrud,WtSet)

         IF( WtCrud .EQ. 0d0) RETURN

         CALL HepEvt_GetBeams(p1,p2)
         CALL HepEvt_GetFfins(p3,p4)

         CALL HepEvt_GetPhotIni(NphIni,PhoIni)  ! ISR
         CALL HepEvt_GetPhotFin(NphFin,PhoFin)  ! FSR

*[[      CALL HepEvt_GetPhotAll(NphAll,PhoAll)  ! not ordered
         CALL KK2f_GetPhotAll(NphAll,PhoAll)  ! ordered in energy

         ss  =  (p1(4)+p2(4))**2 -(p1(3)+p2(3))**2 -(p1(2)+p2(2))**2 -(p1(1)+p2(1))**2
         ss2 =  (p3(4)+p4(4))**2 -(p3(3)+p4(3))**2 -(p3(2)+p4(2))**2 -(p3(1)+p4(1))**2
         vv  = 1d0-ss2/ss
         CALL GLK_Fil1(idv+1, vv, WtMain)

         CMSene = p1(4)+p2(4)
* Sum of energies of ALL photons
         EneTot = 0d0
         DO i=1,NphAll
            EneTot=EneTot+PhoAll(i,4)
         ENDDO
         EneTot = 2d0*EneTot/CMSene
         CALL GLK_Fil1(idv+2, EneTot, WtMain)

* Sum of energies of ISR photons
         EneIni = 0d0
         DO i=1,NphIni
            EneIni=EneIni+PhoIni(i,4)
         ENDDO
         EneIni = 2d0*EneIni/CMSene
         CALL GLK_Fil1(idv+3, EneIni, WtMain)

* Sum of energies of FSR photons
         EneFin = 0d0
         DO i=1,NphFin
            EneFin=EneFin+PhoFin(i,4)
         ENDDO
         EneFin = 2d0*EneFin/CMSene
         CALL GLK_Fil1(idv+4, EneFin, WtMain)

* Photon multiplicity (total)
         CALL GLK_Fil1(idv+7, NphAll*1.00001d0, WtMain)

* Flavour distribution (KF code)
         CALL HepEvt_GetKFfin(KFfin)
         xKF = KFfin+1d-6
         CALL GLK_Fil1(idv+9, xKF, WtMain)
* Photon multiplicity (CRUDE)
         CALL GLK_Fil1(idv+10, NphAll*1.00001d0, 1d0)
         CALL GLK_Fil1(idv+11, NphIni*1.00001d0, 1d0)
         CALL GLK_Fil1(idv+12, NphFin*1.00001d0, 1d0)
* Photon energy ordered
         zenl1 = -50d0
         zenl3 = -50d0
         zenl5 = -50d0
         zenl7 = -50d0
         zenl9 = -50d0
         Ebeam = p1(4)
         IF(NphAll .GE. 1 ) zenl1= LOG10(PhoAll(1,4)/Ebeam)
         IF(NphAll .GE. 3 ) zenl3= LOG10(PhoAll(3,4)/Ebeam)
         IF(NphAll .GE. 5 ) zenl5= LOG10(PhoAll(5,4)/Ebeam)
         IF(NphAll .GE. 7 ) zenl7= LOG10(PhoAll(7,4)/Ebeam)
         IF(NphAll .GE. 9 ) zenl9= LOG10(PhoAll(9,4)/Ebeam)
         Wtz =WtMain
         Wtz =1d0
         CALL GLK_Fil1(idv+14, zenl1, Wtz)
         CALL GLK_Fil1(idv+15, zenl3, Wtz)
         CALL GLK_Fil1(idv+16, zenl5, Wtz)
         CALL GLK_Fil1(idv+17, zenl7, Wtz)
         CALL GLK_Fil1(idv+18, zenl9, Wtz)
*=======================================================================
***         IF( vv .GT. 0.999d0 ) RETURN
*=======================================================================
*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
***         CALL KK2f_Print1(6)
*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
* Exponentiation weight
         WtExp0 = WtSet(201)*WtCrud
         WtExp1 = WtSet(202)*WtCrud
         WtExp2 = WtSet(203)*WtCrud
***   WtExp0 = WtSet(251)*WtCrud  ! Interf OFF
***   WtExp1 = WtSet(252)*WtCrud  ! Interf OFF
         DO k=201,203
            CALL GLK_Fil1(idv+1000+k, WtSet(k)*WtCrud, 1d0)
            CALL GLK_Fil1(idv+2000+k, WtSet(k)*WtCrud, 1d0)
         ENDDO
         DO k=71,74
            CALL GLK_Fil1(idv+1000+k, WtSet(k)*WtCrud, 1d0)
         ENDDO
*=======================================================================
         iTest = NINT(xpar(901))
         nTest = NINT(xpar(902))
         IF(    iTest .EQ. 1) THEN
            CALL STest_SfacISR(  nout,5)
         ELSEIF(iTest .EQ. 2) THEN
            CALL STest_SfacFSR(  nout,5)
         ELSEIF(iTest .EQ. 3) THEN
            CALL STest_BornRmat( nout,5)
         ELSEIF(iTest .EQ. 4) THEN
            CALL STest_UV(       nout,5)
         ELSEIF(iTest .EQ. 5) THEN
            CALL STest_BornXsec( nout,nTest)
         ELSEIF(iTest .EQ. 6) THEN
            CALL STest_IsrSingle(nout,nTest)
         ELSEIF(iTest .EQ. 7) THEN
            CALL STest_FsrSingle(nout,nTest)
         ELSEIF(iTest .EQ. 8) THEN
            CALL STest_MultiPhot(nout,nTest)
         ELSEIF(iTest .EQ. 9) THEN
            CALL STest_Virtual(  nout,nTest)
         ELSEIF(iTest .EQ. 10) THEN
            CALL STest_Rules(    nout,1)
         ELSEIF(iTest .EQ. 11) THEN
*           X-check with DsigOverDtau makes sense ONLY for KeyINT=0
            CALL KK2f_DsigOverDtau(nout,Rho)
         ELSEIF(iTest .EQ. 12) THEN
            CALL STest_NewSingle( nout ,nTest)
         ENDIF
*=======================================================================
      ELSE
         IF( iTest .EQ. 8 ) THEN
*     finalization, printouts
            CALL GLK_Print(idv+1)
            CALL GLK_Print(idv+2)
            CALL GLK_Print(idv+3)
            CALL GLK_Print(idv+4)
            CALL GLK_Print(idv+9)
            CALL GLK_Print(idv+7)
            CALL GLK_Print(idv+10)
            CALL GLK_Print(idv+11)
            CALL GLK_Print(idv+12)
            CALL GLK_Print(idv+14)
            CALL GLK_Print(idv+15)
            CALL GLK_Print(idv+16)
            CALL GLK_Print(idv+17)
            CALL GLK_Print(idv+18)
            DO k=201,203
               CALL GLK_Print( idv+1000+k)
               CALL GLK_Print( idv+2000+k)
            ENDDO
            DO k=71,74
               CALL GLK_Print( idv+1000+k)
            ENDDO
         ENDIF
      ENDIF

      END

