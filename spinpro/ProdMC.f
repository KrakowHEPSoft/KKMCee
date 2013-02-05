** Test of double bremsstrahlung
**	make Mu-start
*///////////////////////////////////////////////////////////////////////////////
      PROGRAM Main
*///////////////////////////////////////////////////////////////////////////////
*//                   main program for MC production                          //
*///////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
*
      INTEGER            ninp, nout
      COMMON /c_MainPro/ ninp, nout
      CHARACTER*4 semaph
      INTEGER     ninp3, ntot2n, ninph, ninp2, ntotin, ijklin
      SAVE
*--------------------------------------------------------------------------------
      ninp =5   ! standard input
      nout =16  ! general output for everybody including glibk
      OPEN(nout,FILE='./pro.output')
      REWIND(nout)
      CALL GLK_SetNout(nout)
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
* *******************
      CALL ProdMC
* *******************
      END


      SUBROUTINE ProdMC
*///////////////////////////////////////////////////////////////////////////////
*//                   MC production                                           //
*///////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
*
      INTEGER            ninp, nout
      COMMON /c_MainPro/ ninp, nout
      INTEGER    imax
*
      PARAMETER (imax = 10000)               ! length ox xpar
      REAL*8        xpar(imax)
*
      CHARACTER*4   semaph,chdum
      INTEGER       kat1,kat2,kat3,kat4
      INTEGER       igroup, ngroup, nevt, loop, iev
      REAL*8        xSecPb, xErrPb
      SAVE
*-------------------------------------------------------------------------------
* Read data for main program
      OPEN( ninp,FILE='./pro.input')
      READ( ninp,'(4i2)') kat1,kat2,kat3,kat4
      WRITE(nout,'(4a6/4i6)')
     $     'kat1','kat2','kat3','kat4',
     $      kat1 , kat2 , kat3 , kat4 
      READ(ninp,'(i10)') nevt
      CLOSE(ninp)
      WRITE(   6,*)   nevt,' requested events '
      WRITE(nout,*)   nevt,' requested events '
*
      CALL KK2f_ReaDataX('../../.KK2f_defaults', 1,imax,xpar)  ! reading general defaults
      CALL KK2f_ReaDataX(         './pro.input', 0,imax,xpar)  ! reading user input
*
      CALL KK2f_Initialize(xpar)                  ! initialize generator
*
      IF(kat1 .EQ. 1) CALL Robol1(-1,xpar)
*-------------------------------------------------------!
*                 main MC loop                          !
*-------------------------------------------------------!
      ngroup = 5000
      iev=0
      DO loop=1,10000000
         DO igroup =1,ngroup
            iev=iev+1
            IF(MOD(iev, ngroup) .EQ. 1) WRITE( 6,*)  'iev= ',iev
            CALL KK2f_Make      ! make single event
*   Control printouts
*          CALL momprt(' YFSPRO ', 6,iev,1,10,pf1,pf2,qf1,qf2,nphot,sphot,KFfin)
*          CALL dumpri('*momini*', 6,iev,1,10,xf1,xf2,nphox,xphot)
            IF(iev .LE. 10) THEN
               CALL lugive('MSTU(11)=16')
               CALL lulist(1)
               CALL lugive('MSTU(11)=6')
               CALL lulist(1)
            ENDIF
            IF(kat1 .EQ. 1) CALL Robol1( 0,xpar) ! histograming
            IF(iev  .EQ.  nevt)     GOTO 300
         ENDDO
         CALL givsem(semaph)    ! check on semaphore flag
         IF(semaph  .EQ.  'STOP') GOTO 300
         CALL  dumpeh(iev)      ! dump partial results on the disk
      ENDDO
 300  CONTINUE
*-------------------------------------------------------!
*            End  main MC loop                          !
*-------------------------------------------------------!
      WRITE(6,*) ' generation finished '
      CALL KK2f_Finalize                          ! final bookkeping, printouts etc.
      CALL KK2f_GetXsecMC(xSecPb, xErrPb)         ! get MC x-section
*
      IF(kat1 .EQ. 1) CALL Robol1(1,xpar)
      CALL  dumpeh(iev)    ! dump final results on the disk
      END

      SUBROUTINE givsem(semaph)
*///////////////////////////////////////////////////////////////////////////////
*//             READ semaphore flag                                           //
*///////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      CHARACTER*4 semaph
      INTEGER     ninp2
* ------------------------------------------------------
      ninp2=2
      OPEN(ninp2,FILE='./semaphore')
      READ(ninp2,'(a4)') semaph
      CLOSE(ninp2)
      END

      SUBROUTINE dumpeh(nev)
*///////////////////////////////////////////////////////////////////////////////
*//             WRITE histos on the disk                                      //
*///////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INTEGER    ijklin, ntotin, ntot2n, ninp2, nev, nouth, icy
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
*////////////////////////////////////////////////////////////////////////////////////////
*//                                                                                    //
*//   Example of histogramming, with acces to events through getters                   //
*//                                                                                    //
*//                                                                                    //
*////////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      REAL*8    xpar(*)
      INTEGER   mode
*
      REAL*8    p1(4),p2(4),p3(4),p4(4)
      INTEGER   NphAll,NphIni,NphFin
      REAL*8    PhoAll(100,4),PhoIni(100,4),PhoFin(100,4)
      REAL*8    WtSet(1000), WtMain, WtCrud
      REAL*8    vv, vmin, vmax, ss, ss2, CMSene, EneTot, EneIni, EneFin
      INTEGER   i, nbv, idv, KFfin
      SAVE
*----------------------------------------------------------------------------------------
      IF(mode .EQ. -1 ) THEN
*     book histograms
         nbv = 400
         vmin= 0
         vmax= 1
         idv  = 50000
         CALL GLK_Book1(idv+1,'d_sigma/dv(v)  $',nbv,vmin,vmax)
         CALL GLK_Book1(idv+2,'Ene total ALL  $',nbv,vmin,vmax)
         CALL GLK_Book1(idv+3,'Ene total ISR  $',nbv,vmin,vmax)
         CALL GLK_Book1(idv+4,'Ene total FSR  $',nbv,vmin,vmax)
         CALL GLK_Book1(idv+7,'NphAll         $',16,0d0,16d0)
         CALL GLK_Book1(idv+9,'KFcode         $',16,1d0,17d0)
*=======================================================================
      ELSEIF(mode .EQ. 0 ) THEN
* Histograming
         CALL HepEvt_GetBeams(p1,p2)
         CALL HepEvt_GetFfins(p3,p4)
         CALL HepEvt_GetPhotAll(NphAll,PhoAll)
         CALL HepEvt_GetPhotIni(NphIni,PhoIni)
         CALL HepEvt_GetPhotFin(NphFin,PhoFin)
         CALL KK2f_GetWtAll(WtMain,WtCrud,WtSet)
         
         ss =  (p1(4)+p2(4))**2 -(p1(3)+p2(3))**2
     $        -(p1(2)+p2(2))**2 -(p1(1)+p2(1))**2
         ss2=  (p3(4)+p4(4))**2 -(p3(3)+p4(3))**2
     $        -(p3(2)+p4(2))**2 -(p3(1)+p4(1))**2
         vv=1d0-ss2/ss
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
         CALL GLK_Fil1(idv+9, KFfin*1.00001d0, WtMain)

c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
* may be this should go to robol2 ????
         IF(WtCrud.NE.0d0)   CALL GPS_MakeTwoPh    !<-- 2 photon bremsstrahlung 
         write(*,*) '!!!!!!!!!!!!!!!'
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
*=======================================================================
      ELSE
*     finalization, printouts
         CALL GLK_Print(idv+1)
         CALL GLK_Print(idv+2)
         CALL GLK_Print(idv+3)
         CALL GLK_Print(idv+4)
         CALL GLK_Print(idv+7)
         CALL GLK_Print(idv+9)
      ENDIF
      END

