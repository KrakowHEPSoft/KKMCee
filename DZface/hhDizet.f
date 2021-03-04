*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//                                                                          //
*//                     Pseudo-CLASS  hhDizet                                //
*//  Purpose:                                                                //
*//  Provides an interface to DIZET to import EW form factors.               //
*//                                                                          //
*//  hhDizet_Initialize:                                                     //
*//  A table is constructed depending on both incoming and outgoing          //
*//  flavors KFini, KFfin. This must be done once at the start of the run.   //
*//  KKMC's BornV must be initialized before using this.                     //
*//                                                                          //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
 
      SUBROUTINE hhDizet_Initialize(xpar)
*-------------------------------------------------------------------------
*     It requires BornV_Initialize(xpar) to be called before
*-------------------------------------------------------------------------
*--   Initializes the Dizet electroweak form factor tables for each
*--   required KFini and KFfin, making a set of tables indexed on both
*--   KFini and KFfin. Except for the b quark, the tables depend only on
*--   whether the quark is up-like or down-like, and whether the final
*--   lepton is charged or neutral, so only one generation of quarks and
*--   one generation of leptons is needed, apart from the b quark.
*-------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      INCLUDE 'hhDizet.h'
      DOUBLE PRECISION xpar(*)
      INTEGER KFini, KFfin, i
*------------------------------------------------------------------
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

      DO i =1,1000
         m_xpar(i) = xpar(i)
      ENDDO
*********************************************************************
      m_ndisk = 17
      write(*,*) '============ hh_InitializeDizet:', ' ============'
      write(*,*) "Saving EW table data to file: DIZET-table1"
      OPEN (m_ndisk, FILE='DIZET-table1',STATUS='REPLACE')
*------------------------------------------------------------------
* Find active chanels
      KFini = KFel
      DO i=401,416
         IF( xpar(i) .EQ. 1d0 ) THEN
            KFfin= i-400
            IF(    KFfin .EQ. KFdown  ) THEN
               WRITE(    *,*) '=========== BornV_StartDZ: DOWN quark ==========='
               WRITE(m_out,*) '=========== BornV_StartDZ: DOWN quark ==========='
               write(*,*) '--------------------------------------------'
               CALL hhDizet_InitDizet(KFini ,KFfin, xpar)
               CALL hhDizet_Tabluj
               CALL hhDizet_WriteTable(KFini, KFfin)
             ELSEIF(KFfin .EQ. KFup    )  THEN
               WRITE(    *,*) '=========== BornV_StartDZ: UP quark ==========='
               WRITE(m_out,*) '=========== BornV_StartDZ: UP quark ==========='
               CALL hhDizet_InitDizet(KFini ,KFfin, xpar)
               CALL hhDizet_Tabluj
               CALL hhDizet_WriteTable(KFini, KFfin)
            ELSE
               WRITE(*,    *) '#### STOP in dizet/BornV_StartEW, wrong KFfin=',KFfin
               WRITE(m_out,*) '#### STOP in dizet/BornV_StartEW, wrong KFfin=',KFfin
               STOP
            ENDIF
         ENDIF
      ENDDO
********************************************************************
      CLOSE(m_ndisk)
      write(*,*) 'KKMC-hh Multiflavor Dizet initialization completed.'
      END ! hhDizet_Initialize

      SUBROUTINE hhDizet_InitDizet( KFini, KFfin,  xpar)
*/////////////////////////////////////////////////////////////////////////////////////////
*//                                                                                     //
*//  Interface to DIZET of the Dubna-Zeuthen EWRC group                                 //
*//  Based on that in KORALZ 4.x                                                        //
*//  Notes:                                                                             //
*//    QED alfinv is separate from alfinv used bremsstrahlung part of KK2f.             //
*//    Note that fermion masses in Dizet are isolated from these in KK2f.               //
*//                                                                                     //
*/////////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      INCLUDE 'BXformat.h'
      DOUBLE PRECISION  xpar(*)
      INTEGER           KFfin, KFini
*
      DOUBLE PRECISION  partz(0:11),partw(3)
      DOUBLE PRECISION  zpard(30)
      INTEGER           Npar(25)
      INTEGER           ihvp,iamt4,Iqcd,imoms,iscre,ialem,imask
      INTEGER           iscal,ibarb,iftjr,ifacr,ifact,ihigs,iafmt
      INTEGER           imass,ii,i,kdumm,kolor
      INTEGER           iewlc,iczak,ihig2,iale2,igfer,iddzz
      DOUBLE PRECISION  amfin,wmass_start,qe,aizor,xolor,qf
      DOUBLE PRECISION  AlStrZ,AlQedZ,DAL5H,AlStrT
      DOUBLE PRECISION  v_tba
*/////////////////////////////////////////////////////////////////////////////
      m_KeyQCD = xpar(53)  ! QCD FSR always=0 for lepton in final state!!!
      m_KeyQCD = 0         ! QCD FSR always=0 for lepton in final state!!!
      m_out   = 16
      m_KFfin = KFfin
cc      m_KFini = xpar(400) !!!
      m_KFini = KFini
*
      m_MZ      = xpar(502) !!!
*
      m_ibox    = xpar(801) !!!
      m_amh     = xpar(805) !!!
      m_amtop   = xpar(806) !!!
*
      m_alfinvMZ  = xpar(808) !!! (128.86674175d0)
      m_alfQCDMZ  = xpar(809) !!! (0.125d0)
      DAL5H  = 0              ! just in case it is input (see Ialem param.)
      m_MW   = 80d0      ! input, rededined by Dizet
      v_tba=1    !the Vtb mixing matrix element
      WRITE(m_out,bxope)
      WRITE(m_out,bxtxt) 'hhDizet_InitDizet, Interface to Dizet 6.xx   '
      WRITE(m_out,bxl1f) m_MZ    ,   'Z mass             ','amz   ','a1'
      WRITE(m_out,bxl1f) m_amh   ,   'Higgs mass         ','amh   ','a2'
      WRITE(m_out,bxl1f) m_amtop ,   'Top mass           ','amtop ','a3'
      WRITE(m_out,bxl1i) m_KFini ,   'KF code of beam    ','KFini ','a5'
      WRITE(m_out,bxl1i) m_KFfin ,   'KF of final fermion','KFfin ','a6'
      WRITE(m_out,bxl1i) m_ibox  ,   'EW box switch      ','ibox  ','a9'
      WRITE(m_out,bxl1f) m_alfinvMZ, 'QED alfa inv. at Z ','alfinv','a1'
      WRITE(m_out,bxl1f) m_alfQCDMZ, 'QCD alfa at Z mass ','alfQCD','a2'
      WRITE(m_out,bxl1f) v_tba     , 'mixing matr. elem. ','v_tba ','a2'
      WRITE(m_out,bxl1f) DAL5H     , 'Delt_ALPH^5_had(MZ)','DAL5H ','a2'
      WRITE(m_out,bxclo)

      qe = m_Qf( m_KFini) ! from BornV.h
      qf = m_Qf( m_KFfin) ! from BornV.h

*  Default values
*      Ihvp  =  1  ! =1,2,3  (Jegerlehner/Eidelman, Jegerlehner(1988), Burkhardt etal.)
*      Iamt4 =  4  ! =0,1,2,3,4 (=4 the best, Degrassi/Gambino)
*      Iqcd  =  3  ! =1,2,3  (approx/fast/lep1, exact/Slow!/Bardin/, exact/fast/Kniehl)
*      Imoms =  1  ! =0,1    (=1 W mass recalculated)
*      Imass =  0  ! =0,1    (=1 test only, effective quark masses)
*      Iscre =  0  ! =0,1,2  ( Remainder terms,
*      Ialem =  3  ! =1,3 or 0,2, (for 1,3 DALH5 not input)
*      Imask =  0  ! =0,1 (=0: Quark masses everywhere; =1 Phys. threshold in the ph.sp.)
*      Iscal =  0  ! =0,1,2,3  ( Kniehl=1,2,3, Sirlin=4)
*      Ibarb =  2  ! =-1,0,1,2 ( Barbieri???)
*      Iftjr =  1  ! =0,1      ( FTJR corrections)
*      Ifacr =  0  ! =0,1,2,3  ( Expansion of delta_r; =0 none; =3 fully, unrecommed.)
*      Ifact =  0  ! =0,1,2,3,4,5 (Expansion of kappa; =0 none )
*      Ihigs =  0  ! =0,1      ( Leading Higgs contribution resummation)
*      Iafmt =  1  ! =0,1      (=0 for old ZF)
* new parameters of 6.x version
*      Iewlc =  1  ! =0,1   (???)
*      Iczak =  1  ! =0,1   (Czarnecki/Kuehn corrections)
*      Ihig2 =  1  ! =0,1   (Two-loop higgs  corrections off,on)
*      Iale2 =  3  ! =1,2,3 (Two-loop constant corrections in delta_alpha)
*      Igfer =  2  ! =0,1,2 (QED corrections for fermi constant)
*      Iddzz =  1  ! =0,1   (??? DD-ZZ game, internal flag)
*
* Input flags in NPAR
      DO i=1,25
         Npar(i) = xpar(900+i)
      ENDDO
      WRITE(m_out,'(a/(a8,i2,a8,i2))')
     $ ' DIZET flags, see routine Dizet for explanation:',
     $ '  Ihvp =',Npar( 1),  ' Iamt4 =',Npar( 2),
     $ '  Iqcd =',Npar( 3),  ' Imoms =',Npar( 4),
     $ ' Imass =',Npar( 5),  ' Iscre =',Npar( 6),
     $ ' Ialem =',Npar( 7),  ' Imask =',Npar( 8),
     $ ' Iscal =',Npar( 9),  ' Ibarb =',Npar(10),
     $ ' IFtjr =',Npar(11),  ' Ifacr =',Npar(12),
     $ ' IFact =',Npar(13),  ' Ihigs =',Npar(14),
     $ ' Iafmt =',Npar(15),  ' Iewlc =',Npar(16),
     $ ' Iczak =',Npar(17),  ' Ihig2 =',Npar(18),
     $ ' Iale2 =',Npar(19),  ' Igfer =',Npar(20),
     $ ' IDDZZ =',Npar(21),  ' IAMW2 =',Npar(22),
     $ ' ISFSR =',Npar(23),  ' IDMWW =',Npar(24),
     $ ' IDSWW =',Npar(25)
*     =================================================================
* Input which is not in Npar
      AlStrZ      =  m_alfQCDMZ          ! input at MZ
      AlQedZ      =  1d0/m_alfinvMZ      ! will be redefibed by dizet6
*     =================================================================
******      !!!!! Dizet 5.x, obsolete !!!!!
******      CALL dizet( Npar, m_MW, m_MZ, m_amtop, m_amh,
******     $            AlQedZ, AlStrZ,
******     $            zpard, partz, partw)
*     =================================================================
*     Dizet 6.x

      CALL DIZET(
     $  Npar,    ! Inp. Integer switches
     $  m_MW,    ! I/O. W mass (Out. for Imoms=Npar(4)=1,3)
     $  m_MZ,    ! I/O. Z mass (Out. for Imoms=Npar(4)=2,4)
     $  m_amtop, ! Inp. t-quark mass
     $  m_amh,   ! Inp. Higgs boson mass
     $  DAL5H,   ! I/O. \Delta_Alpha^5_{had}(MZ), (Inp. Ialem=0,2)(Out. Ialem=1,3)
     $  V_TBA,   !  ???
     $  ALSTRZ,  ! Inp. Alpha_strong(MZ) ???
     $  AlQedZ,  ! Out. Alpha_QED
     $  AlStrT,  ! Out. Alpha_strong(MT)
     $  zpard,   ! Out. zpar(1) = del_r, zpar(2) = del_r_rem, zpar(3) = sw2, ... etc
     $  partz,   ! Out. Z partial decay widths
     $  partw)   ! Out. W partial decay widths
*     =================================================================

      WRITE(m_out,*) '   '
      WRITE(m_out,'(a,  f12.10)') ' Alpha-QED   (MZ)    =',AlQedZ
      WRITE(m_out,'(a,  f15.10)') ' 1/Alpha-QED (MZ)    =',1./AlQedZ
      WRITE(m_out,'(a,  f8.4)')   ' Alfa strong (MZ)    =',AlStrZ
      WRITE(m_out,'(a,  f8.4)')   ' Alfa strong (Mt)    =',AlStrT
      WRITE(m_out,'(a,f15.10)')   'zpard(20): QCD corr.fact. to Z-width (no b)  =',zpard(20)
      WRITE(m_out,'(a,f15.10)')   'zpard(25): QCD corr.fact. to Z-width (into b)=',zpard(25)
      WRITE(m_out,*) '   '
      WRITE(m_out,*) 'zpar-matrix: standard output of dizet:'
      WRITE(m_out,'(a,i2,a,f12.8)') ('    zpar(',ii,')=',zpard(ii),ii=1,30)

      WRITE(m_out,bxope)
      WRITE(m_out,bxtxt) 'DZface_Initializion ended  '
      WRITE(m_out,bxclo)
      END


      SUBROUTINE hhDizet_ReaDataX(DiskFile,iReset,imax,xpar)
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//   Clone of KK2f_ReaDataX                                                        //
*//                                                                                 //
*//   DiskFile  = input file to read                                                //
*//   imax   = maximum index in xpar                                                //
*//   iReset = 1, resets xpar to 0d0                                                //
*//   iTalk=1,     prints echo into standard input                                  //
*//                                                                                 //
*//   Single data card is:    (a1,i4,d15.0,a60)                                     //
*//   First data card: BeginX                                                       //
*//   Last  data card: EndX                                                         //
*//   First character * defines comment card!                                       //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      CHARACTER*(*)     DiskFile
      DOUBLE PRECISION  xpar(*)
      CHARACTER*6       beg6
      CHARACTER*4       end4
      CHARACTER*1       mark1
      CHARACTER*60      comm60
      CHARACTER*80      comm80
      INTEGER           imax,iReset,iTalk
      INTEGER           ninp,i,line,index
      DOUBLE PRECISION  value
*////////////////////////////////////////
*//  Clear xpar and read default Umask //
*////////////////////////////////////////
      iTalk = 1
      IF(iReset .EQ. 1 ) THEN
         iTalk = 0
         DO i=1,imax
            xpar(i)=0d0
         ENDDO
      ENDIF
      ninp = 13
      OPEN(ninp,file=DiskFile)
      IF(iTalk .EQ. 1) THEN
         WRITE(  *,*) '****************************'
         WRITE(  *,*) '*  hhDizet_ReaDataX Starts  *'
         WRITE(  *,*) '****************************'
      ENDIF
* Search for 'BeginX'
      DO line =1,100000
         READ(ninp,'(a6,a)') beg6,comm60
         IF(beg6 .EQ. 'BeginX') THEN
            IF(iTalk .EQ. 1)   WRITE( *,'(a6,a)') beg6,comm60
            GOTO 200
         ENDIF
      ENDDO
 200  CONTINUE
* Read data, 'EndX' terminates data, '*' marks comment
      DO line =1,100000
         READ(ninp,'(a)') mark1
         IF(mark1 .EQ. ' ') THEN
            BACKSPACE(ninp)
            READ(ninp,'(a1,i4,d15.0,a60)') mark1,index,value,comm60
            IF(iTalk .EQ. 1)
     $           WRITE( *,'(a1,i4,g15.6,a60)') mark1,index,value,comm60
            IF( (index .LE. 0) .OR. (index .GE. imax)) GOTO 990
            xpar(index) = value
         ELSEIF(mark1 .EQ. 'E') THEN
            BACKSPACE(ninp)
            READ(  ninp,'(a4,a)') end4,comm60
            IF(iTalk .EQ. 1)   WRITE( *,'(a4,a)') end4,comm60
            IF(end4 .EQ. 'EndX') GOTO 300
            GOTO 991
         ELSEIF(mark1 .EQ. '*') THEN
            BACKSPACE(ninp)
            READ(  ninp,'(a)') comm80
            IF(iTalk .EQ. 1)    WRITE( *,'(a)') comm80
         ENDIF
      ENDDO
 300  CONTINUE
      IF(iTalk .EQ. 1)  THEN
         WRITE(  *,*) '**************************'
         WRITE(  *,*) '*   hhDizet_ReaDataX Ends   *'
         WRITE(  *,*) '**************************'
      ENDIF
      CLOSE(ninp)
      RETURN
*-----------
 990  WRITE(    *,*) '+++ hhDizet_ReaDataX: wrong index= ',index
      STOP
      RETURN
 991  WRITE(    *,*) '+++ hhDizet_ReaDataX: wrong end of data '
      STOP
      END

      SUBROUTINE hhDizet_Tabluj()
*-------------------------------------------------------------------------------
*-- Version of DZface_Tabluj with an extra dimension in the tables for KFini.
*-------------------------------------------------------------------------------
      IMPLICIT NONE 
      INCLUDE 'BornV.h'
      INCLUDE 'hhDizet.h'
      DOUBLE PRECISION ww, x, cosi, amz, amh, amtop, QCDcorR(20)
      DOUBLE PRECISION amzD, gammz1, gammw1
      DOUBLE COMPLEX GSW(100)
      INTEGER i, j, k, kk
*-------------------------------------------------------------------------------
* Initialization
      CALL hhDizet_GetPrm(amzD, m_gammz, gammz1, m_MW, m_GammW, gammw1, m_swsq)
      WRITE(6,'(a)') 'amzD, amh, amtop, swsq, gammz, amw, gammw = '
      WRITE(6,'(a,10f12.7)') '     =',amzD, m_amh, m_amtop, m_swsq, m_gammz, m_MW, m_GammW

c     CALL hhDizet_MakeGSW(-1, 0d0, 0d0, GSW, QCDcorR) ! just sets parameters
*-------------------------------------------------------------------------------
* basic s range LEP1 and below 
      write(*,*) 'hhDizet_Tabluj: basic LEP1 range pretabulation.'
      CALL hhDizet_MakeGSW(-1, ww, 0d0, GSW, QCDcorR) ! just sets parameters
      DO i=0,m_poin1
         ww = m_WminLEP1 *(m_WmaxLEP1/m_WminLEP1)**(DFLOAT(i)/DFLOAT(m_poin1)) !
         IF(MOD(i,10)  .EQ.  0) WRITE(    6,*) 'a: i,ww= ',i,ww
         CALL hhDizet_MakeGSW( 0, ww, 0d0, GSW, QCDcorR) ! at theta=pi, ibox=0
         DO kk=1,m_poinG
            m_cyys(i+1, kk, m_KFini, m_KFfin) = GSW(kk)
         ENDDO
         IF (m_KeyQCD.NE.0) THEN  ! a bit OCD... we don't have QCD FSR
            DO kk=1,m_poinQ
               m_syys(i+1, kk, m_KFini, m_KFfin) = QCDcorR(kk)
            ENDDO
         END IF 
      ENDDO
*-------------------------------------------------------------------------------
* near Z0 resonance 
      write(*,*) 'hhDizet_Tabluj: pretabulation, near Z0.'
      CALL hhDizet_MakeGSW(-1, ww, 0d0, GSW, QCDcorR) ! just sets parameters
      m_WminZ =m_MZ-m_WdelZ
      m_WmaxZ =m_MZ+m_WdelZ
      DO i=0,m_poin2
         ww = m_WminZ +(m_WmaxZ-m_WminZ)*(DFLOAT(i)/DFLOAT(m_poin2))
         IF(MOD(i,10)  .EQ.  0) WRITE(    6,*) 'b: i,ww= ',i,ww
         IF(m_poTh2.EQ.0) THEN
            cosi=0d0
            CALL hhDizet_MakeGSW( 0, ww, cosi, GSW, QCDcorR) ! at theta=pi, ibox=0
            DO kk=1,m_poinG
               m_czzs(i+1,1,kk,m_KFini,m_KFfin)=GSW(kk)
            ENDDO
         ELSE
            DO  j=0,m_poTh2
               cosi = ( -1d0+2d0*DFLOAT(j)/DFLOAT(m_poTh2) )*(1d0-1d-8)
               IF(j  .EQ.        0) cosi=cosi+.3d0/DFLOAT(m_poTh2)
               IF(j  .EQ.  m_poTh2) cosi=cosi-.3d0/DFLOAT(m_poTh2)
               CALL hhDizet_MakeGSW( m_ibox, ww, cosi, GSW, QCDcorR) ! ibox from input
               DO kk=1,m_poinG
                  m_czzs(i+1,j+1,kk,m_KFini,m_KFfin)=GSW(kk)
               ENDDO
            ENDDO
         ENDIF
         IF (m_KeyQCD.NE.0) THEN 
            DO kk=1,m_poinQ
               m_szzs(i+1,kk,m_KFini,m_KFfin) = QCDcorR(kk)
            ENDDO
         ENDIF 
      ENDDO
*-------------------------------------------------------------------------------
* The region of boxes, LEP2
      write(*,*) 'hhDizet_Tabluj: LEP2 energy zone pretabulation.'
      CALL hhDizet_MakeGSW(-1, ww, 0d0, GSW, QCDcorR) ! just sets parameters
      DO  i=0,m_poin3
         ww = m_WmaxLEP1 +(m_WmaxLEP2-m_WmaxLEP1)*(DFLOAT(i)/DFLOAT(m_poin3)) 
         IF(MOD(i,10)  .EQ.  0) WRITE(    6,*) 'c: i,ww= ',i,ww
         DO  j=0,m_poTh3
            cosi = ( -1d0+2d0*DFLOAT(j)/DFLOAT(m_poTh3) )*(1d0-1d-8)
            IF(j  .EQ.        0) cosi=cosi+.3d0/DFLOAT(m_poTh3)
            IF(j  .EQ.  m_poTh3) cosi=cosi-.3d0/DFLOAT(m_poTh3)
            CALL hhDizet_MakeGSW( m_ibox, ww, cosi, GSW, QCDcorR) ! ibox from input
            DO  kk=1,m_poinG
               m_ctts(i+1,j+1,kk,m_KFini,m_KFfin)=GSW(kk)
            ENDDO
         ENDDO
         IF (m_keyQCD.GT.0) THEN 
            DO kk=1,m_poinQ
               m_stts(i+1,kk,m_KFini,m_KFfin) = QCDcorR(kk)
            ENDDO
         END IF 
      END DO
*-------------------------------------------------------------------------------
* The region of boxes, NLC
      write(*,*) 'hhDizet_Tabluj: NLC energy range pretabulation.'
      CALL hhDizet_MakeGSW(-1, ww, 0d0, GSW, QCDcorR) ! just sets parameters
      DO  i=0,m_poin4
         ww  = m_WmaxLEP2+(m_WmaxNLC-m_WmaxLEP2)*(DFLOAT(i)/DFLOAT(m_poin4)) !
         IF(MOD(i,10)  .EQ.  0) WRITE(    6,*) 'd: i,ww= ',i,ww
         DO  j=0,m_poTh4
            cosi = ( -1d0+2d0*DFLOAT(j)/DFLOAT(m_poTh4) )*(1d0-1d-8)
            IF(j  .EQ.        0) cosi=cosi+.3d0/DFLOAT(m_poTh4)
            IF(j  .EQ.  m_poTh4) cosi=cosi-.3d0/DFLOAT(m_poTh4)
            CALL hhDizet_MakeGSW( m_ibox, ww, cosi, GSW, QCDcorR) ! ibox from input
            DO  kk=1,m_poinG
               m_clcs(i+1,j+1,kk,m_KFini,m_KFfin) = GSW(kk)
            ENDDO
         ENDDO
         IF (m_KeyQCD.NE.0) THEN 
            DO kk=1,m_poinQ
               m_slcs(i+1,kk,m_KFini,m_KFfin) = QCDcorR(kk)
            ENDDO
         ENDIF 
      ENDDO
      write(*,*) 'hhDizet_Tabluj: pretabulation finished.'
      END ! hhDizet_Tabluj 


      SUBROUTINE hhDizet_GetPrm( zmass,gamz0,gamzf,wmass,gamw0,gamwf,sin2w)
*/////////////////////////////////////////////////////////////////////////
*//     Gets out params from DIZET common blocks
*/////////////////////////////////////////////////////////////////////////
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      COMMON /cdzzwg/camz,camh,gmu,a0,gamz,gamw,calsz,calst,calxi,calqed
      COMMON /cdzwsm/amw2,amz2,r,r1,r12,r2,amh2,rw,rw1,rw12,rw2,rz,rz1,
     $               rz12,rz2,alr,alrw,alrz,sw2m,cw2m,aksx,r1w,r1w2
*     -------------------------------------------------------------
      sin2w = r1
      wmass = SQRT(amw2)
      zmass = wmass/SQRT(1d0-sin2w)
      gamw0 = gamw
      gamz0 = gamz
      END


      SUBROUTINE hhDizet_MakeGSW(iBox,ww,cosi,GSW,QCDcorR)
*/////////////////////////////////////////////////////////////////////////
*//   gets EWK FFactors and QCD corr. out of Dizet, at ww and theta     //
*//   Prior Dizet initialization is necesary!!!                         //
*//   Used in Tabluj and also in program testing tables                 //
*/////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      INCLUDE 'hhDizet.h'
      INTEGER           iBox
      DOUBLE COMPLEX    GSW(*)
      DOUBLE PRECISION  ww,cosi,QCDcorR(*),ROW,PIRE
      DOUBLE PRECISION  QCDcorD(0:14)
      DOUBLE COMPLEX    xff(4),xfem,xfota
      INTEGER           i,j,kk,k,kdumm,kolor,ibfla,KFf
      DOUBLE PRECISION  ss,uu,tt,aizor,qe,qf,dum1,dum2
c[[[[[[[[[[[
      INTEGER icont
      DATA icont /0/
c]]]]]]]]]]]
      SAVE
*     ---------------------------------------------------
      IF(    iBox.EQ.-1) THEN
         IF ( ABS(m_KFfin)  .EQ.  5) THEN
            ibfla=1
         ELSE
            ibfla=0
         ENDIF
         qe = m_Qf( m_KFini) ! from BornV.h
         qf = m_Qf( m_KFfin) ! from BornV.h

         DO k=1,m_poinG
            m_GSW(k) =0d0
         ENDDO
      ELSEIF( iBox.EQ.0 .OR. iBox.EQ.1 ) THEN
         ss  = ww*ww
         tt  = ss*(1d0-cosi)/2d0
         uu  = ss*(1d0+cosi)/2d0
         CALL rokanc( iBox,ibfla,-uu,-ss,-tt,qe,qf,xff,xfem,xfota) ! check signs of s,u,t!!!
         IF (ABS(m_KFfin).EQ.12) THEN
c[[[[[[[[[[[[
c   CALL RHOCC(-s-t,-t,S,0D0,-1D0,-1D0,0D0,ROW)
c   ROW=ROW+1D0/m_AlfInv/m_pi/2*(-3D0/2*LOG(s/m_MW**2)+1D0/2*(LOG(-t/s))**2-m_pi**2/6+2D0)
c]]]]]]]]]]]]
ccc            CALL rhocc(ss,tt,-ss+tt,-qe,qe,0D0,0D0,ROW) ! this is for nunu electron chanel
            CALL RHOCC(-ss+tt,+tt,ss,   0D0,-1D0,-1D0,0D0,ROW)
         ELSE
            row=777.5555         ! UNUSED!!!!
         ENDIF
         CALL hhDizet_MakeQCDcor(ww,QCDcorD)
         DO  kk=1,4
            GSW(kk)=xff(kk)
         ENDDO
         GSW(5 )=ROW  !!  123.456d0 ! UNUSED!!!!
         GSW(6 )=xfem
ccc         pire=-0.0351079942
ccc         GSW(6)=1 -PIRE
ccc         write(*,*) 'xfem= ',xfem, ' ss= ',ss
         GSW(7 )=xfota
         GSW(8 )=QCDcorD(1)-1d0     !  obsolete!!! test of pretabulation,may be replace with alpha_s???
* translate QCD R-factors and corrections
         QCDcorR(1) = 1d0
         QCDcorR(2) = 1d0
         QCDcorR(3) = 0d0
         QCDcorR(4) = 0d0
         KFf = ABS(m_KFfin)
         IF(     KFf .EQ. 1 ) THEN
            QCDcorR(1) = QCDcorD(3)  ! D quark, R Vector
            QCDcorR(2) = QCDcorD(4)  ! D quark, R Axial
            QCDcorR(3) = QCDcorD(13) ! singlet
            QCDcorR(4) = QCDcorD(14) ! f_1
         ELSEIF( KFf .EQ. 2 ) THEN
            QCDcorR(1) = QCDcorD(1)  ! U quark, R Vector
            QCDcorR(2) = QCDcorD(2)  ! U quark, R Axial
            QCDcorR(3) = QCDcorD(13) ! singlet
            QCDcorR(4) = QCDcorD(14) ! f_1
         ELSEIF( KFf .EQ. 3 ) THEN
            QCDcorR(1) = QCDcorD(7)  ! S quark, R Vector
            QCDcorR(2) = QCDcorD(8)  ! S quark, R Axial
            QCDcorR(3) = QCDcorD(13) ! singlet
            QCDcorR(4) = QCDcorD(14) ! f_1
         ELSEIF( KFf .EQ. 4 ) THEN
            QCDcorR(1) = QCDcorD(5)  ! C quark, R Vector
            QCDcorR(2) = QCDcorD(6)  ! C quark, R Axial
            QCDcorR(3) = QCDcorD(13) ! singlet
            QCDcorR(4) = QCDcorD(14) ! f_1
         ELSEIF( KFf .EQ. 5 ) THEN
            QCDcorR(1) = QCDcorD(11) ! B quark, R Vector
            QCDcorR(2) = QCDcorD(12) ! B quark, R Axial
            QCDcorR(3) = QCDcorD(13) ! singlet
            QCDcorR(4) = QCDcorD(14) ! f_1
         ENDIF
      ELSE
         WRITE(*,*) '++++ STOP in hhDizet_MakeGSW, wrong iBox =',iBox
         STOP
      ENDIF
      END

      SUBROUTINE hhDizet_MakeQCDcor(ww,QCDcor)
*/////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      DOUBLE PRECISION ww, ene, AlfQED             ! Input
      DOUBLE PRECISION ALPHTT,ALPHXI, QCDcor(0:14) ! Output
      INTEGER  i,j,k
*     -------------------------------------------------------------
      ene     = MAX(ww,20d0)
      AlfQED  = 1d0/m_alfinvMZ
      CALL qcdcof(ene, m_amtop, m_swsq, AlfQED, m_alfQCDMZ, ALPHTT, ALPHXI, QCDcor) !
      END
*/


      DOUBLE PRECISION FUNCTION hhDizet_AlfCor(svar)
*------------------------------------------------------------------------------
*- Calculates running alpha correction factor (relative to Thomson limit).
*------------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      DOUBLE PRECISION svar, ww
      COMPLEX*16 Xrok(4), Xfot, Xfot5
      IF (svar.EQ.0d0) THEN
          hhDizet_AlfCor = 1d0
      ELSE
        ww = SQRT(svar)
        CALL ROKANC(0,0,-svar/2d0,-svar,-svar/2d0,1d0,1d0,Xrok,Xfot,Xfot5)
        hhDizet_AlfCor = REAL(Xfot)
      END IF
      END ! hhDizet_AlfCor

      SUBROUTINE hhDizet_WriteTable(KFi, KFf)
*------------------------------------------------------------------------------
*-- Writes electroweak table for given initial and final fermions. 
*-- KFi, KFf > 0 is assumed. KFi must be a quark, KFf must be a lepton.
*------------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      INCLUDE 'hhDizet.h'
      INTEGER KFi, KFf, hhIO_OpenFile, i, j, k
      DOUBLE PRECISION ww, cosi, amzD, gammz1, gammw1
      CHARACTER(LEN=60) tableFile

      IF ((KFi.GT.5).OR.(KFi.LT.1)) THEN
         write(*,*) 'hhDizet_WriteTable: Invalid quark code ', KFi
      END IF
      IF ((KFf.GT.16).OR.(KFf.LT.11)) THEN
         write(*,*)'hhDizet_WriteTable: Invalid lepton code ', KFf
      END IF

      write(m_ndisk,*) "hhDizet_WriteTable: KFi, KFf =", KFi, KFf
*-- params out of dizet 
      CALL hhDizet_GetPrm( amzD,m_gammz,gammz1,m_MW,m_GammW,gammw1,m_swsq)
      WRITE(m_ndisk, m_fmt0) m_MZ, m_amh, m_amtop, m_swsq, m_gammz, m_MW, M_GammW
*/////////////////////////////////////////////////////////////////
*//   basic s range LEP1 and below        m_cyy(m_poin1+1, 7)   //
*/////////////////////////////////////////////////////////////////
      DO i=0, m_poin1
         ww = m_WminLEP1 *(m_WmaxLEP1/m_WminLEP1)**(DFLOAT(i)/DFLOAT(m_poin1)) !
         WRITE(m_ndisk,m_fmt1) 'a',i,ww
         WRITE(m_ndisk,     *)                    (m_cyys(i+1,    k,KFi,KFf),k=1,m_poinG) ! EW
         IF (m_KeyQCD.NE.0) WRITE(m_ndisk,m_fmt2) (m_syys(i+1,    k,KFi,KFf),k=1,m_poinQ) ! QCD
      ENDDO
*/////////////////////////////////////////////////////////////////
*/              near Z0 resonance    m_czz(m_poin2+1, 7)        //
*/////////////////////////////////////////////////////////////////
      m_WminZ =m_MZ-m_WdelZ
      m_WmaxZ =m_MZ+m_WdelZ
      DO i=0,m_poin2
         ww = m_WminZ +(m_WmaxZ-m_WminZ)*(DFLOAT(i)/DFLOAT(m_poin2))
         DO  j=0,m_poTh2
            cosi = ( -1d0+2d0*DFLOAT(j)/DFLOAT(m_poTh2) )*(1d0-1d-8)
            IF(j  .EQ.        0) cosi=cosi+.3d0/DFLOAT(m_poTh2)
            IF(j  .EQ.  m_poTh2) cosi=cosi-.3d0/DFLOAT(m_poTh2)
            WRITE(m_ndisk,m_fmt1) 'b',i,ww,j,cosi
            WRITE(m_ndisk,     *)                 (m_czzs(i+1,j+1,k,KFi,KFf),k=1,m_poinG) ! EW
         ENDDO
         IF (m_KeyQCD.NE.0) WRITE(m_ndisk,m_fmt2) (m_szzs(i+1,    k,KFi,KFf),k=1,m_poinQ) ! QCD
      ENDDO
*/////////////////////////////////////////////////////////////////////
*//   the region of boxes, LEP2,   m_ctt(m_poin3+1, m_poTh3+1, 7)   //
*/////////////////////////////////////////////////////////////////////
*     write(*,*) 'DZface_WriteFile: m_poTh3=', m_poTh3
      DO  i=0,m_poin3
         ww = m_WmaxLEP1 +(m_WmaxLEP2-m_WmaxLEP1)*(DFLOAT(i)/DFLOAT(m_poin3)) !
         DO  j=0,m_poTh3
            cosi = ( -1d0+2d0*DFLOAT(j)/DFLOAT(m_poTh3) )*(1d0-1d-8)
            IF(j  .EQ.        0) cosi=cosi+.3d0/DFLOAT(m_poTh3)
            IF(j  .EQ.  m_poTh3) cosi=cosi-.3d0/DFLOAT(m_poTh3)
            WRITE(m_ndisk,m_fmt1) 'c',i,ww,j,cosi
            WRITE(m_ndisk,    *)                  (m_ctts(i+1,j+1,k,KFi,KFf),k=1,m_poinG) ! EW
         ENDDO
         IF (m_KeyQCD.NE.0) WRITE(m_ndisk,m_fmt2) (m_stts(i+1,    k,KFi,KFf),k=1,m_poinQ) ! QCD
      ENDDO
*/////////////////////////////////////////////////////////////////////
*//   the region of boxes, NLC,    m_clc(m_poin4+1, m_poTh4+1, 7)   //
*/////////////////////////////////////////////////////////////////////
      DO  i=0,m_poin4
         ww  = m_WmaxLEP2+(m_WmaxNLC-m_WmaxLEP2)*(DFLOAT(i)/DFLOAT(m_poin4)) !
         DO  j=0,m_poTh4
            cosi = ( -1d0+2d0*DFLOAT(j)/DFLOAT(m_poTh4) )*(1d0-1d-8)
            IF(j  .EQ.        0) cosi=cosi+.3d0/DFLOAT(m_poTh4)
            IF(j  .EQ.  m_poTh4) cosi=cosi-.3d0/DFLOAT(m_poTh4)
            WRITE(m_ndisk,m_fmt1) 'd',i,ww,j,cosi
            WRITE(m_ndisk,     *)                 (m_clcs(i+1,j+1,k,KFi,KFf),k=1,m_poinG) ! EW
         ENDDO
         IF (m_KeyQCD.NE.0) WRITE(m_ndisk,m_fmt2) (m_slcs(i+1,    k,KFi,KFf),k=1,m_poinQ) ! QCD
      ENDDO
      END ! hhDizet_WriteTable


      LOGICAL FUNCTION hhDizet_ReadTable(KFi, KFf)
*------------------------------------------------------------------------------
*      OBSOLETE !!!!   OBSOLETE !!!!   OBSOLETE !!!!   OBSOLETE !!!!
*-- Reads an EW table file. Returns 1 if the parameters are consistent, else 0
*------------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      INCLUDE 'hhDizet.h'
      INTEGER KFi, KFf, hhIO_OpenRead, i, j, k, n, n1, n2
      DOUBLE PRECISION ww, cosi, amz, amh, amtop
      CHARACTER chr
      CHARACTER(LEN=60) tableFile
      LOGICAL ex
*
      WRITE(*,*) 'hhDizet_ReadTable: START!!!! KFi, KFf=,',KFi, KFf
*-- See if a table file has been saved...
      IF( KFi .EQ. 1) THEN
          tableFile ='DIZET-table1'
      ELSE
          tableFile ='DIZET-table2'
      ENDIF
      INQUIRE(FILE=tableFile, EXIST=ex)
      IF (ex) THEN 
         OPEN (27, FILE=tableFile)
         write(*,*) 'Reading saved EW table:',tableFile
      ELSE 
         hhDizet_ReadTable = .FALSE.
         write(*,*) 'No compatible saved EW tables. Generating...  ', tableFile
         RETURN
      END IF 

* Follows BornV_ReadFile in BornV_StartEW.f......
* But consistency check has already been done. 
      READ(27,m_fmt0) amz, amh, amtop, m_swsq, m_gammz, m_MW, m_GammW
      WRITE(6,'(a)') 'amz, amh, amtop, swsq, gammz, amw, gammw = '
      WRITE(6,'(a,10f12.7)') '     =',amz, amh, amtop, m_swsq, m_gammz, m_MW, m_GammW

* basic range
      DO i=0, m_poin1
         READ(27,m_fmt1) chr,n,ww
         READ(27,m_fmt2) (m_cyys(i+1,k,KFi,KFf),k=1,m_poinG) ! EW
         IF (m_KeyQCD.NE.0) READ(27,m_fmt2) (m_syys(i+1,k,KFi,KFf),k=1,m_poinQ) ! QCD
      ENDDO
* Z pole range
      DO i=0,m_poin2
         DO  j=0,m_poTh2
            READ(27,m_fmt1) chr,n,ww
            READ(27,m_fmt2) (m_czzs(i+1,j+1,k,KFi,KFf),k=1,m_poinG) ! EW
         ENDDO
         IF (m_KeyQCD.NE.0) READ(27,m_fmt2) (m_szzs(i+1,k,KFi,KFf),k=1,m_poinQ) ! QCD
      ENDDO
* LEP2 range
      DO  i=0,m_poin3
         DO  j=0,m_poTh3
            READ(27,m_fmt1) chr,n1,ww,n2,cosi
            READ(27,m_fmt2) (m_ctts(i+1,j+1,k,KFi,KFf),k=1,m_poinG) ! EW
         ENDDO
         IF (m_KeyQCD.NE.0) READ(27,m_fmt2) (m_stts(i+1,k,KFi,KFf),k=1,m_poinQ) ! QCD
      ENDDO
* NLC range
      DO  i=0,m_poin4
         DO  j=0,m_poTh4
            READ(27,m_fmt1) chr,n1,ww,n2,cosi
            READ(27,m_fmt2) (m_clcs(i+1,j+1,k,KFi,KFf),k=1,m_poinG) ! EW
         ENDDO
         IF (m_KeyQCD.NE.0) READ(27,m_fmt2) (m_slcs(i+1,k,KFi,KFf),k=1,m_poinQ) ! QCD
      ENDDO
      CLOSE(27)
      hhDizet_ReadTable = .TRUE.
      WRITE(*,*) 'hhDizet_ReadTable: ENDED!!!!'
      END ! hhDizet_ReadTable

********************************************************************************
********************************************************************************
********************************************************************************

      SUBROUTINE hhDizet_GetSMpar(swsq, gammz, MW, GammW)
*------------------------------------------------------------------------------
*     Temporary inteface for transfering EW formfactors to c++ code
*     SM parameters calculated by DIZET
*------------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      INCLUDE 'hhDizet.h'
      DOUBLE PRECISION swsq, gammz, MW, GammW
      swsq   = m_swsq
      gammz  = m_gammz
      MW     = m_MW
      GammW  = m_GammW
      END

      DOUBLE PRECISION FUNCTION hhDizet_GetELW(it,md,KFi,KFf,i,j,k)
*------------------------------------------------------------------------------
*     Temporary inteface for transfering EW formfactors to c++ code
*------------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      INCLUDE 'hhDizet.h'
      INTEGER KFi,KFf,i,j,k, it, md
      DOUBLE COMPLEX z

      IF( it .EQ. 1) THEN
         z= m_cyys(i,    k,KFi,KFf)  ! single point in theta
      ELSE IF( it .EQ. 2) THEN
         z= m_czzs(i,j,k,KFi,KFf)    ! i,j,k f77 indexing!!!
      ELSE IF( it .EQ. 3) THEN
         z= m_ctts(i,j,k,KFi,KFf)
      ELSE IF( it .EQ. 4) THEN
         z= m_clcs(i,j,k,KFi,KFf)
      ELSE
         write(*,*) '+++++ hhDizet_GetELW: STOP, wrong it =',it
         STOP
      ENDIF
      IF( md .EQ. 0) hhDizet_GetELW= DREAL(z)
      IF( md .EQ. 1) hhDizet_GetELW= DIMAG(z)
      END ! hhDizet_GetELW


      DOUBLE PRECISION FUNCTION hhDizet_GetQCD(it,KFi,KFf,i,k)
**------------------------------------------------------------------------------
*     Temporary inteface for transfering QCD formfactors to c++ code
*------------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      INCLUDE 'hhDizet.h'
      INTEGER KFi,KFf,i,k, it, md

      IF( it .EQ. 1) THEN
         hhDizet_GetQCD= m_syys(i,k,KFi,KFf)  ! single point in theta
      ELSE IF( it .EQ. 2) THEN
         hhDizet_GetQCD= m_szzs(i,k,KFi,KFf)  ! i,k, f77 indexing
      ELSE IF( it .EQ. 3) THEN
         hhDizet_GetQCD= m_stts(i,k,KFi,KFf)
      ELSE IF( it .EQ. 4) THEN
         hhDizet_GetQCD= m_slcs(i,k,KFi,KFf)
      ELSE
         write(*,*) '+++++ hhDizet_GetELW: STOP, wrong it =',it
         STOP
      ENDIF
      END ! hhDizet_GetQCD
