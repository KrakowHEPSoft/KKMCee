      PROGRAM MAIN
!     ***********************************
!
! To execute:  make BeMaTabs-dvi
!              make BeMaTabs-ps
!
! BenchMark Tables for BHLUMI 4.0x, to be used in CPC article
!
!     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common / cglib / b(50000)
      COMMON / INOUT  / NINP,NOUT
!----------------------------------------------------------------------
! Communicates with READAT
      CHARACTER*60             Tesnam,TeXfile,Dname    ,Hname
      COMMON / PLflag / lendan,Tesnam,TeXfile,Dname(50),Hname(50)
!----------------------------------------------------------------------
      CALL GLIMIT(50000)
      NINP=  5
      NOUT= 16
      Tesnam    = 'BeMaTabs'
      OPEN( NOUT, file='output-'//Tesnam)
      CALL GOUTPU(NOUT)
!--------------------------------------------------------
! Stored histograms and corresponding histograms 
!--------------------------------------------------------
      lendan = 0
      lendan = lendan+1
      Dname(lendan)  = '../prod2/prod2.data.2037M'  ! April 96
      Hname(lendan)  = '../prod2/prod2.hst.2037M'   ! April 96
      Dname(lendan)  = '../prod2/prod2.data'        ! Current
      Hname(lendan)  = '../prod2/bhl.hst'           ! Current
      lendan = lendan+1
      Dname(lendan)  = '../obis2/obis2.data.2118M'   ! April 96
      Hname(lendan)  = '../obis2/obis2.hst.2118M'    ! April 96
      Dname(lendan)  = '../obis2/obis2.data'      ! Current
      Hname(lendan)  = '../obis2/bhl.hst'         ! Current
      lendan = lendan+1
      Dname(lendan)  = '../llog2/llog2.data.2114M'  ! April 96
      Hname(lendan)  = '../llog2/llog2.hst.2114M'   ! April 96
      Dname(lendan)  = '../llog2/llog2.data'        ! Current
      Hname(lendan)  = '../llog2/bhl.hst'           ! Current
!==========================================================
! Table for workshop
      CALL wshop
!==========================================================
!--------------------------------------------------------
!  dumping histogram for control
!--------------------------------------------------------
      NOUTH=20
      OPEN(NOUTH,file='dump.hst')
      CALL grfile(nouth,' ','N')
      CALL grout( 0,ICY,' ')
      CALL grend(tname)
      CLOSE(nout)
      END




      SUBROUTINE wshop
!     ************************
! tables for wshop proceedings
!     *********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(PI     =  3.1415926535897932D0)
      CHARACTER*80      BXOPE,BXCLO,BXTXT,BXL1I,BXL1F,BXL2F,BXL1G,BXL2G
      PARAMETER(
     $BXOPE =  '(//1X,15(5H=====)    )',
     $BXTXT =  '(1X,1H=,                  A48,25X,    1H=)',
     $BXL1I =  '(1X,1H=,I17,                 16X, A20,A12,A7, 1X,1H=)',
     $BXL1F =  '(1X,1H=,F17.8,               16X, A20,A12,A7, 1X,1H=)',
     $BXL2F =  '(1X,1H=,F17.8, 4H  +-, F11.8, 1X, A20,A12,A7, 1X,1H=)',
     $BXL1G =  '(1X,1H=,G17.8,               16X, A20,A12,A7, 1X,1H=)',
     $BXL2G =  '(1X,1H=,G17.8, 4H  +-, F11.8, 1X, A20,A12,A7, 1X,1H=)',
     $BXCLO =  '(1X,15(5H=====)/   )'    )
      SAVE
      COMMON / INOUT  / NINP,NOUT
!-----------------------------------------------------------------------
! Communicates with MAIN
      CHARACTER*60             Tesnam,TeXfile,Dname    ,Hname
      COMMON / PLflag / lendan,Tesnam,TeXfile,Dname(50),Hname(50)
!-----------------------------------------------------------------------
! Communicates with READAT
      COMMON / PARGEN / CMSENE,TMING,TMAXG,VMAXG,XK0,KEYOPT,KEYRAD
      COMMON / PAROBL / NPHI,NTHE,TMINW,TMAXW,TMINN,TMAXN,VMAXE,KEYTRI
!-----------------------------------------------------------------------
      DIMENSION cut(100)
      DIMENSION bin1(100),err1(100)
      CHARACTER*80 title
      DIMENSION BorW(10),BorN(10)
!-------------------------
! Parameters for tables
      DIMENSION    idl(5)
      CHARACTER*16 capt(6)
      CHARACTER*8  fmt(3), fmtx,fmty
!-----------------------------------
      CHARACTER*64 mcapt
!-----------------------------------
      CHARACTER*64 captab(50)
      DATA captab /
     $'This is from the output file 4.x-cpc/figs/BeMaTabs.tex.',
     $'Total cross sections for various symmetric Wide-Wide,', 
     $'event selections from BHLUMI, in nanobarn units.', 
     $'BHLUM4$_{pht}$ and BHLUM4$_{ZVP}$ denote standard multiphoton', 
     $'BHLUMI event generator.', 
     $'OBI+LMG denotes so-called OLDBIS+LUMLOG$_{h.o.}$ ', 
     $'additive recipe while ', 
     $'BHLPV follows SABSPV multiplicative recipe', 
     $'realized using OLDBIS and LUMLOG cross section.', 
     $'Only in BHLUM4$_{ZVP}$ Z exchange, up-down', 
     $'interference and vacuum polarization are switched on.',
     $'Center of mass energy $\\protect\\sqrt{s}=92.3 GeV$.', 
     $'Not calculated cross sections set to zero.',
     $'% end-of-caption'/
!-----------------------------------
!==================================================================

! ------------------------------------------------
! ----------------- BHLUMI -----------------------
! ------------------------------------------------
      jene=1
      WRITE(   6,BXOPE)
      WRITE(   6,BXTXT) Dname(jene)
      WRITE(   6,BXTXT) Hname(jene)
!-- Reading input-data file used for MC
      CALL ReaDat(Dname(jene))
      WRITE(   6,BXL1F) cmsene, 'total CMS energy  ','<<<---','=='
!-- Restore histograms from j-th  directory
      NINPH=10
      OPEN(NINPH,file=Hname(jene))
      CALL GRFILE(NINPH,' ',' ')
      CALL GRIN(0,0,0)
      CLOSE(NINPH)
!-- End of restoring
      JBHL= 30000 +2000
      kbhl= JBHL +10000000
!------- normalize in nanobarns
      KeyGen = 3
      DO itr=1,1
! BHLUMI O(alf2)exp
        CALL cumhis(KeyGen,JBHL+340+itr,kbhl+340+itr) ! BARE1
        CALL cumhis(KeyGen,JBHL+350+itr,kbhl+350+itr) ! CALO2
        CALL cumhis(KeyGen,JBHL+360+itr,kbhl+360+itr) ! SICAL2
        CALL cumhis(KeyGen,JBHL+300+itr,kbhl+300+itr) ! SICAL
! VP+Z included
        CALL cumhis(KeyGen,JBHL+380+itr,kbhl+380+itr) ! CALO2
        CALL cumhis(KeyGen,JBHL+390+itr,kbhl+390+itr) ! SICAL2
        CALL cumhis(KeyGen,JBHL+370+itr,kbhl+370+itr) ! SICAL
! Dummy BARE1
      CALL gopera(kbhl+371,'+',kbhl+371,kbhl+400+itr,0d0,0d0)! BARE1
      ENDDO
! ------------------------------------------------
! ----------------- OLDBIS -----------------------
! ------------------------------------------------
      jene=2
      WRITE(   6,BXOPE)
      WRITE(   6,BXTXT) Dname(jene)
      WRITE(   6,BXTXT) Hname(jene)
!-- Reading input-data file used for MC
      CALL ReaDat(Dname(jene))
      WRITE(   6,BXL1F) cmsene, 'total CMS energy  ','<<<---','=='
!-- Restore histograms from j-th  directory
      NINPH=10
      OPEN(NINPH,file=Hname(jene))
      CALL GRFILE(NINPH,' ',' ')
      CALL GRIN(0,0,0)
      CLOSE(NINPH)
!-- End of restoring
      JBIS= 10000 +2000
      kbis= JBis +10000000
!------- normalize in nanobarns
      KeyGen = 1
      DO itr=1,1
        CALL cumhis(KeyGen,JBIS+300+itr,kbis+300+itr) ! SICAL
        CALL cumhis(KeyGen,JBIS+340+itr,kbis+340+itr) ! BARE1
        CALL cumhis(KeyGen,JBIS+350+itr,kbis+350+itr) ! CALO2
        CALL cumhis(KeyGen,JBIS+360+itr,kbis+360+itr) ! SICAL2
      ENDDO
! ------------------------------------------------
! ----------------- LUMLOG -----------------------
! ------------------------------------------------
      jene=3
      WRITE(   6,BXOPE)
      WRITE(   6,BXTXT) Dname(jene)
      WRITE(   6,BXTXT) Hname(jene)
!-- Reading input-data file used for MC
      CALL ReaDat(Dname(jene))
      WRITE(   6,BXL1F) cmsene, 'total CMS energy  ','<<<---','=='
!-- Restore histograms from j-th  directory
      NINPH=10
      OPEN(NINPH,file=Hname(jene))
      CALL GRFILE(NINPH,' ',' ')
      CALL GRIN(0,0,0)
      CLOSE(NINPH)
!-- End of restoring
!------- normalize in nanobarns
      JLOG= 20000 +4000     ! from JDU  O(alf3)exp -O(alf1) 
      klog= JLOG +10000000
      KeyGen = 2
      DO itr=1,1
        CALL cumhis(KeyGen,JLOG+340+itr,klog+340+itr) ! BARE1
        CALL cumhis(KeyGen,JLOG+350+itr,klog+350+itr) ! CALO2
        CALL cumhis(KeyGen,JLOG+360+itr,klog+360+itr) ! SICAL2
        CALL cumhis(KeyGen,JLOG+300+itr,klog+300+itr) ! SICAL
      ENDDO
      LLOG= 20000 +2000     ! from JDA  O(alf1) for SABSPV
      mlog= LLOG +10000000
      DO itr=1,1
        CALL cumhis(KeyGen,LLOG+340+itr,mlog+340+itr) ! BARE1
        CALL cumhis(KeyGen,LLOG+350+itr,mlog+350+itr) ! CALO2
        CALL cumhis(KeyGen,LLOG+360+itr,mlog+360+itr) ! SICAL2
        CALL cumhis(KeyGen,LLOG+300+itr,mlog+300+itr) ! SICAL
      ENDDO
!==================================================================
!==================================================================
!==================================================================
!==================================================================
!-----------------------------------------------------------------------
!          BORN BORN BORN Born Wide and Narrow
!-----------------------------------------------------------------------
!----------------------------------
!              BARE1              !
!----------------------------------
      THbar1 = 0.024d0
      THbar2 = 0.058d0
      Nthe = 16
      PAD = (THbar2-THbar1)/Nthe
      BorW(1)  = BORNB(CMSENE,THbar1       , THbar2       )
      BorN(1)  = BORNB(CMSENE,THbar1 +1*PAD, THbar2 -1*PAD)
      write(6,*) 'BARE1  BorW(1),BorN(1) ---->' ,BorW(1),BorN(1)
!---------------------------------!
!              CALO1              !
!---------------------------------!
! The same as in BARE1
      write(6,*) 'CALO1  BorW(1),BorN(1) ---->' ,BorW(1),BorN(1)
!---------------------------------!
!              CALO2              !
!---------------------------------!
      THbar1 = 0.024d0
      THbar2 = 0.058d0
      Nthe = 16
      PAD = (THbar2-THbar1)/Nthe
      BorW(2)  = BORNB(CMSENE,THbar1 +1*PAD, THbar2 -1*PAD)
      BorN(2)  = BORNB(CMSENE,THbar1 +2*PAD, THbar2 -4*PAD)
      write(6,*) 'CALO2  BorW(1),BorN(2) ---->' ,BorW(2),BorN(2)
!---------------------------------!
!              SICAL2             !
!---------------------------------!
! The same as CALO2
      write(6,*) 'SICAL2 BorW(2),BorN(2) ---->' ,BorW(2),BorN(2)
!---------------------------------!
!              SICAL              !
!---------------------------------!
! The same as CALO2
      write(6,*) 'SICAL  BorW(2),BorN(2) ---->' ,BorW(2),BorN(2)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      khyb= 50000 +2000
      kbpv= 60000 +2000
      kdif= khyb +10000000
      nlog= llog +20000000
      kdip= llog +30000000
! ######################
      itr=1  ! Wide-Wide
! ######################
! BARE1   O+L and ((O+L)-B)/Born
      fac = 1/BorN(1)
      IF(itr .EQ. 1) fac = 1/BorW(1)
      CALL gopera(kbis+340+itr,'+',klog+340+itr,khyb+340+itr,1d0,1d0)
!>      CALL gopera(khyb+340+itr,'-',kbhl+340+itr,kdif+340+itr,fac,fac)
!     BhPV and (BhPV-B)/Born
!     (O(alf3)LL-O(alf1))LL +O(alf1)LL => O(alf3)LL
!     (Oldbis   -O(alf1)LL) /Born      => NLL/Born
!     O(alf3)LLexp * (NLL/Born)        => Corr
!     O(alf3)LLexp + Corr              => BhPV
!     (BhPV  - Bhlum4) /Born           => BhPdiff
      CALL gopera(klog+340+itr,'+',mlog+340+itr,nlog+340+itr,1d0,1d0)
      CALL gopera(kbis+340+itr,'-',mlog+340+itr,kbpv+340+itr,fac,fac)
      CALL gopera(nlog+340+itr,'*',kbpv+340+itr,kbpv+340+itr,1d0,1d0)
      CALL gopera(nlog+340+itr,'+',kbpv+340+itr,kbpv+340+itr,1d0,1d0)
!>      CALL gopera(kbpv+340+itr,'-',kbhl+340+itr,kdip+340+itr,fac,fac)
! CALO2   O+L and ((O+L)-B)/Born
      fac = 1/BorN(2)
      IF(itr .EQ. 1) fac = 1/BorW(2)
      CALL gopera(kbis+350+itr,'+',klog+350+itr,khyb+350+itr,1d0,1d0)
!>      CALL gopera(khyb+350+itr,'-',kbhl+350+itr,kdif+350+itr,fac,fac)
!    BhPV and (BhPV-B)/Born
      CALL gopera(klog+350+itr,'+',mlog+350+itr,nlog+350+itr,1d0,1d0)
      CALL gopera(kbis+350+itr,'-',mlog+350+itr,kbpv+350+itr,fac,fac)
      CALL gopera(nlog+350+itr,'*',kbpv+350+itr,kbpv+350+itr,1d0,1d0)
      CALL gopera(nlog+350+itr,'+',kbpv+350+itr,kbpv+350+itr,1d0,1d0)
!>      CALL gopera(kbpv+350+itr,'-',kbhl+350+itr,kdip+350+itr,fac,fac)
! SICAL2  O+L and ((O+L)-B)/Born
      fac = 1/BorN(2)
      IF(itr .EQ. 1) fac = 1/BorW(2)
      CALL gopera(kbis+360+itr,'+',klog+360+itr,khyb+360+itr,1d0,1d0)
!>      CALL gopera(khyb+360+itr,'-',kbhl+360+itr,kdif+360+itr,fac,fac)
!    BhPV and (BhPV-B)/Born
      CALL gopera(klog+360+itr,'+',mlog+360+itr,nlog+360+itr,1d0,1d0)
      CALL gopera(kbis+360+itr,'-',mlog+360+itr,kbpv+360+itr,fac,fac)
      CALL gopera(nlog+360+itr,'*',kbpv+360+itr,kbpv+360+itr,1d0,1d0)
      CALL gopera(nlog+360+itr,'+',kbpv+360+itr,kbpv+360+itr,1d0,1d0)
!>      CALL gopera(kbpv+360+itr,'-',kbhl+360+itr,kdip+360+itr,fac,fac)
! SICAL   O+L and ((O+L)-B)/Born
      fac = 1/BorN(2)
      IF(itr .EQ. 1) fac = 1/BorW(2)
      CALL gopera(kbis+300+itr,'+',klog+300+itr,khyb+300+itr,1d0,1d0)
!>      CALL gopera(khyb+300+itr,'-',kbhl+300+itr,kdif+300+itr,fac,fac)
!    BhPV and (BhPV-B)/Born
      CALL gopera(klog+300+itr,'+',mlog+300+itr,nlog+300+itr,1d0,1d0)
      CALL gopera(kbis+300+itr,'-',mlog+300+itr,kbpv+300+itr,fac,fac)
      CALL gopera(nlog+300+itr,'*',kbpv+300+itr,kbpv+300+itr,1d0,1d0)
      CALL gopera(nlog+300+itr,'+',kbpv+300+itr,kbpv+300+itr,1d0,1d0)
!>      CALL gopera(kbpv+300+itr,'-',kbhl+300+itr,kdip+300+itr,fac,fac)
!==================================================================

!********************************************************************
!                             TABLES                                !
!********************************************************************
!--------------------------------------------------------
!----    initialize GPLOT
!--------------------------------------------------------
      CALL GPLINT(0)
      NOUFIG=11
      TeXfile   = 'BeMaTabs.tex'
      OPEN(NOUFIG,file=TeXfile)
      CALL GPLCAP(-NOUFIG)
!========================================================
!---------------------------------!
!              BARE1              !
!---------------------------------!
      nbare1=1400241
      jbare1=2400241
      kbare1=3400241
      mbare1=4400241
      ibare1=5400241
! Repack MC results
      CALL RePack5(kbhl+340+itr,nbare1)
      CALL RePack5(khyb+340+itr,mbare1)
      CALL RePack5(kbpv+340+itr,jbare1)
      CALL RePack5(kbhl+400+itr,kbare1)
!----------------------------------
! BARE1 WW trigger O(alf2)exp
      capt(1)='$z_{min}$'
      capt(2)=' BHLUM4$_{pht}$'
      capt(3)=' OBI+LMG '
      capt(4)=' BHLPV '
      capt(5)=' BHLUM4$_{ZVP}$'
      fmt(1)='F10.3'
      fmt(2)='F10.3'
      fmt(3)='F8.3'
      idl(1)=nbare1
      idl(2)=mbare1
      idl(3)=jbare1
      idl(4)=kbare1
!------------------------------------------------
      CALL gplcapt(captab)
      mcapt = ' (a) BARE1, $Born = 175.977 nb$'
      CALL gpltab2(4,idl,capt,mcapt,fmt,' ',' ',' ')
!========================================================
!---------------------------------!
!              CALO2              !
!---------------------------------!
      ncalo2=1500241
      jcalo2=2500241
      kcalo2=3500241
      mcalo2=4500241
      icalo2=5500241
! Repack MC results
      CALL RePack5(kbhl+350+itr,ncalo2)
      CALL RePack5(khyb+350+itr,mcalo2)
      CALL RePack5(kbpv+350+itr,jcalo2)
      CALL RePack5(kbhl+380+itr,kcalo2)
!----------------------------------
! CALO2 WW trigger O(alf2)exp
      fmt(1)='F10.3'
      fmt(2)='F10.3'
      fmt(3)='F8.3'
      idl(1)=ncalo2
      idl(2)=mcalo2
      idl(3)=jcalo2
      idl(4)=kcalo2
!------------------------------------------------
      mcapt = ' (c) CALO2, $Born = 140.018 nb$'
      CALL gpltab2(4,idl,capt,mcapt,fmt,'S',' ',' ')
!========================================================
!---------------------------------!
!              SICAL2             !
!---------------------------------!
      nsical2=1600241
      jsical2=2600241
      ksical2=3600241
      msical2=4600241
      isical2=5600241
! Repack MC results
      CALL RePack5(kbhl+360+itr,nsical2)
      CALL RePack5(khyb+360+itr,msical2)
      CALL RePack5(kbpv+360+itr,jsical2)
      CALL RePack5(kbhl+390+itr,ksical2)
!----------------------------------
! SICAL2 WW trigger O(alf2)exp
      fmt(1)='F10.3'
      fmt(2)='F10.3'
      fmt(3)='F8.3'
      idl(1)=nsical2
      idl(2)=msical2
      idl(3)=jsical2
      idl(4)=ksical2
!------------------------------------------------
      mcapt = ' (d) SICAL2, $Born = 140.018 nb$'
      CALL gpltab2(4,idl,capt,mcapt,fmt,'S',' ',' ')
!========================================================
!---------------------------------!
!              SICAL              !
!---------------------------------!
      nsical=1000241
      jsical=2000241
      ksical=3000241
      msical=4000241
      isical=5000241
! Repack MC results
      CALL RePack5(kbhl+300+itr,nsical)
      CALL RePack5(khyb+300+itr,msical)
      CALL RePack5(kbpv+300+itr,jsical)
      CALL RePack5(kbhl+370+itr,ksical)
!----------------------------------
! SICAL WW trigger O(alf2)exp
      fmt(1)='F10.3'
      fmt(2)='F10.3'
      fmt(3)='F8.3'
      idl(1)=nsical
      idl(2)=msical
      idl(3)=jsical
      idl(4)=ksical
!------------------------------------------------
      mcapt = ' (e) SICAL, $Born = 140.018 nb$'
      CALL gpltab2(4,idl,capt,mcapt,fmt,'S',' ','E')
!========================================================
!   The end of plots writing
      CALL gplend
!--------------------------------------------------------
!--------------------------------------------------------
      END

      SUBROUTINE RePack5(id,id2)
!     ***********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION bin1(100),err1(100)
      id1=id
      nbin=5
      zmin= 0.0d0
      zmax= 1.0d0 
      DO i=1,nbin
         k=44-i*8
         bin1(i)= gi(id1,k)
         err1(i)=gie(id1,k)
      ENDDO
      CALL gbook1(id2,'clone    $',nbin,zmin,zmax)
      CALL gidopt(id2,'ERRO')
      CALL   gpak(id2,bin1)
      CALL  gpake(id2,err1)
!--------------------------------------------------------
      END
