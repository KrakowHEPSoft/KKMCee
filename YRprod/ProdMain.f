*///////////////////////////////////////////////////////////////////////////////
*//
*//   make Quark-start
*//   make Mu-start
*//   make Tau-start
*//
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
      DOUBLE PRECISION        xpar(imax)
*
      CHARACTER*4   semaph,chdum
      INTEGER       kat1,kat2,kat3,kat4
      INTEGER       igroup, ngroup, nevt, loop, iev
      DOUBLE PRECISION        xSecPb, xErrPb
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
      CALL RobAll_Initialize(xpar)
*-------------------------------------------------------!
*                 main MC loop                          !
*-------------------------------------------------------!
      ngroup = 5000
      iev=0
      DO loop=1,10000000
        DO igroup =1,ngroup
          iev=iev+1
          IF(MOD(iev, ngroup) .EQ. 1) WRITE( 6,*)  'iev= ',iev
          CALL KK2f_Make                          ! make single event
*   Control printouts
*          CALL momprt(' YFSPRO ', 6,iev,1,10,pf1,pf2,qf1,qf2,nphot,sphot,KFfin)
*          CALL dumpri('*momini*', 6,iev,1,10,xf1,xf2,nphox,xphot)
          IF(iev .LE. 10) THEN
             CALL pygive('MSTU(11)=16')
             CALL pylist(1)
             CALL pygive('MSTU(11)=6')
             CALL pylist(1)
          ENDIF
          CALL RobAll_Make    ! histograming
          IF(iev  .EQ.  nevt)     GOTO 300
        ENDDO
        CALL givsem(semaph)  ! check on semaphore flag
        IF(semaph  .EQ.  'STOP') GOTO 300
        CALL  dumpeh(iev)  ! dump partial results on the disk
      ENDDO
 300  CONTINUE
*-------------------------------------------------------!
*            End  main MC loop                          !
*-------------------------------------------------------!
      WRITE(6,*) ' generation finished '
      CALL KK2f_Finalize                          ! final bookkeping, printouts etc.
      CALL KK2f_GetXsecMC(xSecPb, xErrPb)         ! get MC x-section
*
      CALL RobAll_Finalize
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


      SUBROUTINE RobAll_Initialize(xpar)
*////////////////////////////////////////////////////////////////////////////////////////
*//                                                                                    //
*////////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'RobAll.h'
      DOUBLE PRECISION    xpar(*)
      INTEGER i,j,k
      WRITE(*,*) 'm_nBinSAB,m_nBinANG=',m_nBinSAB,m_nBinANG
*----------------------------------------------------------------------------------------
      CALL GLK_Book1(m_kSAB,'Sigma and AFB grand storage $', m_nBinSAB, 0d0,1d0*m_nBinSAB) !
      CALL GLK_Book1(m_kANG,'dSigma/dTheta grand storage $', m_nBinANG, 0d0,1d0*m_nBinANG) !
*------------------------------------------------------
      DO i=1,100
         m_sum1(i)=0d0
         m_sum2(i)=0d0
      ENDDO
      END

      SUBROUTINE RobAll_Make
*////////////////////////////////////////////////////////////////////////////////////////
*//                                                                                    //
*//   Example of histogramming, with acces to events through getters                   //
*//                                                                                    //
*////////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'RobAll.h'
      DOUBLE PRECISION    PhoAll(100,4),PhoIni(100,4),PhoFin(100,4)
      DOUBLE PRECISION    WtList(1000), WtMain, WtCrud
      DOUBLE PRECISION    vv
      INTEGER             NphAll,NphIni,NphFin
      INTEGER             i,j,k, KFfin
      INTEGER             iVmax,iSel,iThe,iMod,jThe
      DOUBLE PRECISION    WtBest,Wt
      DOUBLE PRECISION    p1(4),p2(4),p3(4),p4(4)
      DOUBLE PRECISION    CosThe, vvPR, mProp, mPropPlus, mPropMinus, xup,xlo, dCos !
      DOUBLE PRECISION    ECM,mFin,mPropAl,CosThe1,CosThe2,CosThePL,CosPRD !
      DOUBLE PRECISION    xFlag(100),Sign, WTtest1, WTtest2
      INTEGER icont
      DATA    icont/0/
      icont=icont+1
* Get momenta and weights
      CALL HepEvt_GetBeams(p1,p2)
      CALL HepEvt_GetFfins(p3,p4)
***   CALL  KK2f_GetBeams(    p1,p2)
***   CALL  KK2f_GetFermions( p3,p4)
***   CALL KarFin_GetFermions(p3,p4)
      CALL HepEvt_GetPhotAll(NphAll,PhoAll)
      CALL HepEvt_GetPhotIni(NphIni,PhoIni)
      CALL HepEvt_GetPhotFin(NphFin,PhoFin)
*
      CALL KK2f_GetWtList(WtMain,WtList)
      CALL KK2f_GetWtCrud(WtCrud)

      CALL HepEvt_GetKFfin(KFfin) ! Flavour distribution (KF code)
*============================================================================================
      IF(WtCrud.EQ.0d0) RETURN   !!!! PROTECT against bad kinenatics
*============================================================================================
      CALL AngMake(p1,p2,p3,p4,ECM,mFin,mPropAl,CosThe1,CosThe2,CosThePL,CosPRD)    !
* definition of the angle and other cut-off params
      CosThe=CosThe1            ! for comparison with Zfitter
****  CosThe=CosPRD             ! for comparison with KKsem
      vv=1d0-(mFin/ECM)**2
* vvP of abstract propagator, from inside of MC generator, makes sense only for IFI OFF!!!
      vvPR = 0d0
      IF( (NphAll .NE. 0) .AND. (WtCrud .NE. 0d0) ) THEN
         CALL BornV_GetVV(vvPR)
      ENDIF
      mProp = ECM * SQRT(1d0-vvPR)
*============================================================================================
*   ------------------ histogramming starts here --------------------
*                         SAB = Sigma+AFB
*============================================================================================
* First primitive cuts, the same for all fermions
      DO iVmax=1,m_nVmax
         m_WtSel(iVmax)=0d0
         IF( vv.LT.m_VmaxList(iVmax) ) m_WtSel(iVmax)=1d0
      ENDDO
      m_WtThe(1)                  =  1d0          ! Total
      m_WtThe(2)                  =  m_WtThe(1)   ! Forward
      IF(CosThe.LT.0d0) m_WtThe(2)= -m_WtThe(1)   ! Backwrd
* Four model weights
      m_WtMod(1)=WtMain                  ! WtList(203) Main Weight O(alf2),  IFI ON
      m_WtMod(2)=WtList(253)             ! Aux. Weight O(alf2),  IFI OFF
      m_WtMod(3)=WtList(202)             ! Phys. Prec  O(alf1),  IFIon
      m_WtMod(4)=WtList( 73)             ! Second ord. EEX2 O(alf2)
      m_WtMod(5)=WtList( 74)             ! Third order EEX3 O(alf3)
      m_WtMod(6)=WtList(213)-WtList(203) ! Virtual PAIRS alone (IFI on)
***   m_WtMod(6)=WtList(263)-WtList(253) ! Virtual PAIRS alone (IFI off)
      m_WtMod(7)=WtList(203)-WtList(253) ! IFI alone
      m_WtMod(8)=WtList(203)-WtList(202) ! PhysPrec. O(alf2)-O(alf1)  IFIon
***   IF(icont.LE.100) WRITE(*,'(a,10F12.5)') '=====>',(m_WtMod(i), i=1,8) !
* Primitive selection, no CosTheta cut
      DO iSel=1,m_nVmax
         DO iThe=1,2
            DO iMod=1,8
               Wt= m_WtSel(iSel)*m_WtThe(iThe) *m_WtMod(iMod)
               IF( Wt.NE.0d0) CALL RobAll_FilSAB(KFfin,iSel,iThe,iMod,Wt)
            ENDDO
         ENDDO
      ENDDO
      m_sum1(1) = MAX(m_sum1(1),ABS(m_WtSel(5)*m_WtThe(1) *m_WtMod(7))) ! TEST
      m_sum1(5) = m_sum1(5) + m_WtSel(5)*m_WtThe(1) *m_WtMod(7) ! TEST
      m_sum2(5) = m_sum2(5) + m_WtSel(5)*m_WtThe(2) *m_WtMod(7) ! TEST
cc      IF(m_WtSel(5).NE.0d0)
cc     $ WRITE(*,'(a,5g20.7)') '::::: ',m_WtThe(1)*m_WtMod(1), m_WtThe(1)*m_WtMod(7) !
cc      IF(m_WtSel(5).NE.0d0)
cc     $ WRITE(*,'(a,i5,5g20.7)') '***** ', NphAll,CosThe, WtList(203), WtList(253), m_WtMod(7)/m_WtMod(1)!
*============================================================================================
* Now YR selection, wiht CosTheta cuts, flavour dependent
      DO i=1,100
         m_WtSel(i)=0d0
      ENDDO
      IF(KFfin.EQ.13) THEN      ! Muons
*     Muons IAleph5,IAleph6
         IF((mFin.GT.0.10d0*ECM)  .AND. (ABS(CosThe).LT.0.95d0)) m_WtSel(m_nYR+ 1)=1d0 !
         IF((mFin.GT.0.90d0*ECM)  .AND. (ABS(CosThe).LT.0.95d0)) m_WtSel(m_nYR+ 2)=1d0 !
*     Muons, IDelphi5, IDelphi6
         IF((mFin.GT.      75d0)  .AND. (ABS(CosThe).LT.0.95d0)) m_WtSel(m_nYR+ 3)=1d0 !
         IF((mFin.GT.0.85d0*ECM)  .AND. (ABS(CosThe).LT.0.95d0)) m_WtSel(m_nYR+ 4)=1d0 !
*     Muons, ILT9,ILT10,ILT11
         IF((mFin.GT.      75d0)  .AND. (ABS(CosThe).LT.0.90d0)) m_WtSel(m_nYR+ 5)=1d0 !
         IF((mFin.GT.0.85d0*ECM)  .AND. (ABS(CosThe).LT.0.90d0)) m_WtSel(m_nYR+ 6)=1d0 !
         IF((mFin.GT.0.85d0*ECM)  .AND. (ABS(CosThe).LT.1.00d0)) m_WtSel(m_nYR+ 7)=1d0 !
*     Muons, IOpal6,IOpal7,IOpal8,IOpal9
         IF((mProp.GT.0.10d0*ECM) .AND. (ABS(CosThe).LT.0.95d0)) m_WtSel(m_nYR+ 8)=1d0 !
         IF((mProp.GT.0.10d0*ECM) .AND. (ABS(CosThe).LT.1.00d0)) m_WtSel(m_nYR+ 9)=1d0 !
         IF((mProp.GT.0.85d0*ECM) .AND. (ABS(CosThe).LT.0.95d0)) m_WtSel(m_nYR+10)=1d0 !
         IF((mProp.GT.0.85d0*ECM) .AND. (ABS(CosThe).LT.1.00d0)) m_WtSel(m_nYR+11)=1d0 !
      ELSEIF(KFfin.EQ.15) THEN      ! Taus
*     Taus,  IAleph7,IAleph8 
         IF((mFin.GT.0.10d0*ECM)  .AND. (ABS(CosThe).LT.0.95d0)) m_WtSel(m_nYR+ 1)=1d0 !
         IF((mFin.GT.0.90d0*ECM)  .AND. (ABS(CosThe).LT.0.95d0)) m_WtSel(m_nYR+ 2)=1d0 !
*     Taus,   IDelphi7  IDelphi8
         IF((mFin.GT.      75d0)  .AND. (ABS(CosThe).LT.0.95d0)) m_WtSel(m_nYR+ 3)=1d0 !
         IF((mFin.GT.0.85d0*ECM)  .AND. (ABS(CosThe).LT.0.95d0)) m_WtSel(m_nYR+ 4)=1d0 !
*     Taus,   ILT12, ILT13,ILT14
         IF((mFin.GT.      75d0)  .AND. (ABS(CosThe).LT.0.92d0)) m_WtSel(m_nYR+ 5)=1d0 !
         IF((mFin.GT.0.85d0*ECM)  .AND. (ABS(CosThe).LT.0.92d0)) m_WtSel(m_nYR+ 6)=1d0 !
         IF((mFin.GT.0.85d0*ECM)  .AND. (ABS(CosThe).LT.1.00d0)) m_WtSel(m_nYR+ 7)=1d0 !
*     Taus,    IOpal10, IOpal11, IOpal12, IOpal13
         IF((mProp.GT.0.10d0*ECM) .AND. (ABS(CosThe).LT.0.90d0)) m_WtSel(m_nYR+ 8)=1d0 !
         IF((mProp.GT.0.10d0*ECM) .AND. (ABS(CosThe).LT.1.00d0)) m_WtSel(m_nYR+ 9)=1d0 !
         IF((mProp.GT.0.85d0*ECM) .AND. (ABS(CosThe).LT.0.90d0)) m_WtSel(m_nYR+10)=1d0 !
         IF((mProp.GT.0.85d0*ECM) .AND. (ABS(CosThe).LT.1.00d0)) m_WtSel(m_nYR+11)=1d0 !
      ELSE                      ! Quarks
         mPropPlus  = mProp     !!!!!!!!
         mPropMinus = mProp     !!!!!!!!
*     Quarks, IAleph1,  IAleph2
         IF((mPropPlus.GT.0.10d0*ECM)                             ) m_WtSel(m_nYR+ 1)=1d0 !
         IF((mPropPlus.GT.0.90d0*ECM) .AND.(ABS(CosThe).LT.0.95d0)) m_WtSel(m_nYR+ 2)=1d0 !
*     Quarks, IDelphi1, IDelphi2
         IF((mPropPlus.GT.0.10d0*ECM)                             ) m_WtSel(m_nYR+ 3)=1d0 !
         IF((mPropPlus.GT.0.85d0*ECM)                             ) m_WtSel(m_nYR+ 4)=1d0 !
*     Quarks, ILT1, ILT2, ILT3
         IF((mPropMinus.GT.0.10d0*ECM)                            ) m_WtSel(m_nYR+ 5)=1d0 !
         IF((mPropMinus.GT.      60d0)                            ) m_WtSel(m_nYR+ 6)=1d0 !
         IF((mPropMinus.GT.0.85d0*ECM)                            ) m_WtSel(m_nYR+ 7)=1d0 !
*     Quarks,  IOpal1, IOpal2
         IF((mPropMinus.GT.0.10d0*ECM)                            ) m_WtSel(m_nYR+ 8)=1d0 !
         IF((mPropMinus.GT.0.85d0*ECM)                            ) m_WtSel(m_nYR+ 9)=1d0 !
      ENDIF
      DO iSel = m_nYR+1, m_nYR+11
         DO iThe=1,2
            DO iMod=1,8
               Wt= m_WtSel(iSel)*m_WtThe(iThe) *m_WtMod(iMod)
               IF( Wt.NE.0d0) CALL RobAll_FilSAB(KFfin,iSel,iThe,iMod,Wt)
            ENDDO
         ENDDO
      ENDDO
      IF( KFfin .EQ. 13 ) THEN
*-------------------------------------------------------------------
* Realistic from Zbyszek, muons only
*-------------------------------------------------------------------
         CALL buker4(xFlag,Sign)
ccc         IF(icont.LE.100) WRITE(*,'(a,20f5.1)') '---->', (xFlag(i), i=1,20)!
* normal
         m_WtSel(m_nYR+12)= xFlag(11)! Aleph5
         m_WtSel(m_nYR+13)= xFlag(12)! Aleph6
         m_WtSel(m_nYR+14)= xFlag(13)! Delphi4
         m_WtSel(m_nYR+15)= xFlag(14)! Delphi5
         m_WtSel(m_nYR+16)= xFlag(15)! LT9
         m_WtSel(m_nYR+17)= xFlag(16)! LT10
         m_WtSel(m_nYR+18)= xFlag(17)! Opal6
         m_WtSel(m_nYR+19)= xFlag(18)! Opal7
         m_WtThe(1)                  =  1d0 ! Total
         m_WtThe(2)                  =  m_WtThe(1) ! Forward
         IF(Sign.LT.0d0) m_WtThe(2)  = -m_WtThe(1) ! Backwrd
         DO iSel = m_nYR+12, m_nYR+19
            DO iThe=1,2
               DO iMod=1,8
                  Wt= m_WtSel(iSel)*m_WtThe(iThe) *m_WtMod(iMod)
                  IF( Wt.NE.0d0) CALL RobAll_FilSAB(KFfin,iSel,iThe,iMod,Wt) !
               ENDDO
            ENDDO
         ENDDO
* photonic
         m_WtSel(m_nYR+20)= xFlag( 1) ! ALEPH-12
         m_WtSel(m_nYR+21)= xFlag( 2) ! ALEPH-15
         m_WtSel(m_nYR+22)= xFlag( 3) ! DELPHI9
         m_WtSel(m_nYR+23)= xFlag( 4) ! DELPHI12
         m_WtSel(m_nYR+24)= xFlag( 5) ! LT15
         m_WtSel(m_nYR+25)= xFlag( 6) ! LT18
         m_WtSel(m_nYR+26)= xFlag( 7) ! Opal14
         m_WtSel(m_nYR+27)= xFlag( 8) ! Opal15
         m_WtSel(m_nYR+28)= xFlag( 9) ! Opal20
         m_WtSel(m_nYR+29)= xFlag(10) ! Opal21
         m_WtThe(1)                  =  1d0 ! Total
         m_WtThe(2)                  =  0d0 ! AFB is null
         DO iSel = m_nYR+20, m_nYR+29
            DO iThe=1,2
               DO iMod=1,8
                  Wt= m_WtSel(iSel)*m_WtThe(iThe) *m_WtMod(iMod)
                  IF( Wt.NE.0d0) CALL RobAll_FilSAB(KFfin,iSel,iThe,iMod,Wt) !
               ENDDO
            ENDDO
         ENDDO
*-------------------------------------------------------------------
*              ANG = Dsigma/dCosTheta, muons only
*-------------------------------------------------------------------
         DO jThe= m_jThe1, m_jThe2
            m_WtThe(jThe)=0d0
         ENDDO
         dCos=2d0/m_nBinThe
         DO jThe= m_jThe1, m_jThe2
            xlo= -1d0+dCos*(jThe- m_jThe1)
            xup= -1d0+dCos*(jThe- m_jThe1+1)
            IF( CosThe.GE.xlo .AND. CosThe.LE.xup) m_WtThe(jThe)=1d0/dCos !
         ENDDO
         DO iSel = 1,m_nVmax
            DO jThe = m_jThe1, m_jThe2
               DO iMod = 1,8
                  Wt= m_WtSel(iSel)*m_WtThe(jThe) *m_WtMod(iMod)
                  IF( Wt.NE.0d0) CALL RobAll_FilANG(KFfin,iSel,jThe,iMod,Wt)
               ENDDO
            ENDDO
         ENDDO
         DO iSel = m_nYR+1, m_nYR+2 ! Aleph5,6 for muons
            DO jThe = m_jThe1, m_jThe2
               DO iMod = 1,8
                  Wt= m_WtSel(iSel)*m_WtThe(jThe) *m_WtMod(iMod)
                  IF( Wt.NE.0d0) CALL RobAll_FilANG(KFfin,iSel,jThe,iMod,Wt)
               ENDDO
            ENDDO
         ENDDO
      ENDIF
      END

      SUBROUTINE RobAll_Finalize
*////////////////////////////////////////////////////////////////////////////////////////
*//   printouts                                                                        //
*////////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'RobAll.h'
      INTEGER ikFf,iSel,iThe,iMod
      DOUBLE PRECISION xTot,eTot
      iKFf =13
      iThe =1
      DO iSel =  m_iSel1, m_iSel2
         iMod = 1               !Best
         CALL GetSabBin(ikFf,iSel,iThe,iMod, xTot,eTot)
         WRITE(*,'(a,i3,5g20.7)') '===iSel,xTot,eTot: ',iSel,xTot,eTot
         iMod = 7               !IFI
         CALL GetSabBin(ikFf,iSel,iThe,iMod, xTot,eTot)
         WRITE(*,'(a,i3,5g20.7)') '===iSel,xTot,eTot: ',iSel,xTot,eTot
c         iMod =2                !Pairs off
c         CALL GetSabBin(ikFf,iSel,iThe,iMod, xTot,eTot)
c         WRITE(*,*) '===iSel,xTot,eTot: ',iSel,xTot,eTot
c         iMod =6                !Pairs ON
c         CALL GetSabBin(ikFf,iSel,iThe,iMod, xTot,eTot)
c         WRITE(*,*) '---iSel,xTot,eTot: ',iSel,xTot,eTot
      ENDDO
      iSel = 5                  ! iSel = 5 is vmax=0.20
      iThe = 1
      WRITE(*,*) '------- iSel,iThe = ',iSel,iThe
      DO iMod= 1,8
         CALL GetSabBin(ikFf,iSel,iThe,iMod, xTot,eTot)
         WRITE(*,'(a,i3,5g20.7)') '////iMod,xTot,eTot: ',iMod,xTot,eTot        
      ENDDO
      iThe = 2
      WRITE(*,*) '------- iSel,iThe = ',iSel,iThe
      DO iMod= 1,8
         CALL GetSabBin(ikFf,iSel,iThe,iMod, xTot,eTot)
         WRITE(*,'(a,i3,5g20.7)') '////iMod,xTot,eTot: ',iMod,xTot,eTot        
      ENDDO
      WRITE(*,*) '------- test of histogramming --------'
      WRITE(*,'(a,5g20.7)') '*****m_sum1(5),m_sum2(5)=  ',m_sum1(5),m_sum2(5) !
      WRITE(*,'(a,5g20.7)') '*****m_sum1(1) max WT7  =  ',m_sum1(1)!
      END

      SUBROUTINE GetSabBin(ikFf,iSel,iThe,iMod, xbin,ebin)
*////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'RobAll.h'
      INTEGER           ikFf,iSel,iThe,iMod,iBin
      DOUBLE PRECISION  xbin,ebin
      iBin= (iKFf-m_iKFf1  )*(m_iSel2-m_iSel1+1)*(m_iThe2-m_iThe1+1)*(m_iMod2-m_iMod1+1)  !
     $                        +(iSel -m_iSel1  )*(m_iThe2-m_iThe1+1)*(m_iMod2-m_iMod1+1)  !
     $                                            +(iThe -m_iThe1  )*(m_iMod2-m_iMod1+1)  !
     $                                                                +(iMod -m_iMod1  )+1!
      IF( iBin .LT.1 .OR. iBin.GT.m_nBinSAB ) THEN
         WRITE(*,*) ' ----STOP in GetSABBin: iBin,m_nBinSAB=',iBin,m_nBinSAB !
         WRITE(*,*) ' ----ikFf,iSel,iThe,iMod=', ikFf,iSel,iThe,iMod !
         STOP
      ENDIF
      CALL  GLK_GetBin(m_kSAB,iBin,xbin)
      CALL  GLK_GetErr(m_kSAB,iBin,ebin)
      END

      SUBROUTINE RobAll_FilSAB(iKFf,iSel,iThe,iMod,WT)
*////////////////////////////////////////////////////////////////////////////////////////
*//                                                                                    //
*////////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'RobAll.h'
      INTEGER           iKFf,iSel,iThe,iMod
      DOUBLE PRECISION  WT
      INTEGER           iBin
      iBin= (iKFf-m_iKFf1  )*(m_iSel2-m_iSel1+1)*(m_iThe2-m_iThe1+1)*(m_iMod2-m_iMod1+1)  !
     $                        +(iSel -m_iSel1  )*(m_iThe2-m_iThe1+1)*(m_iMod2-m_iMod1+1)  !
     $                                            +(iThe -m_iThe1  )*(m_iMod2-m_iMod1+1)  !
     $                                                                +(iMod -m_iMod1  )+1!
      CALL GLK_Fil1(m_kSAB, iBin-0.5d0,  Wt)
      END

      SUBROUTINE RobAll_FilANG(jKFf,jSel,jThe,jMod,WT)
*////////////////////////////////////////////////////////////////////////////////////////
*//                                                                                    //
*////////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'RobAll.h'
      INTEGER           jKFf,jSel,jThe,jMod
      DOUBLE PRECISION  WT
      INTEGER           iBin
      iBin= (jKFf-m_jKFf1  )*(m_jSel2-m_jSel1+1)*(m_jThe2-m_jThe1+1)*(m_jMod2-m_jMod1+1)  !
     $                        +(jSel -m_jSel1  )*(m_jThe2-m_jThe1+1)*(m_jMod2-m_jMod1+1)  !
     $                                            +(jThe -m_jThe1  )*(m_jMod2-m_jMod1+1)  !
     $                                                                +(jMod -m_jMod1  )+1!
      CALL GLK_Fil1(m_kANG, iBin-0.5d0,  Wt)
      END

      SUBROUTINE  AngMake(p1,p2,p3,p4,ECM,mFin,mPropAl,CosThe1,CosThe2,CosThePL,CosPRD) !
*////////////////////////////////////////////////////////////////////////////////////////
*//                                                                                    //
*////////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      DOUBLE PRECISION    p1(4),p2(4),p3(4),p4(4)
      DOUBLE PRECISION    ECM,mFin,mPropAl,CosThe1,CosThe2,CosThePL,CosPRD !
*
      DOUBLE PRECISION    zAleph,sPropAl,sIni,sFin,CMSene
      DOUBLE PRECISION    yy1,yy2,theta1,theta2,SinThe1,SinThe2
*
      sIni =  (p1(4)+p2(4))**2 -(p1(3)+p2(3))**2 -(p1(2)+p2(2))**2 -(p1(1)+p2(1))**2 !
      sFin =  (p3(4)+p4(4))**2 -(p3(3)+p4(3))**2 -(p3(2)+p4(2))**2 -(p3(1)+p4(1))**2 !
      mFin = DSQRT(sFin)
      ECM  = DSQRT(sIni)
      CosThe1 = p3(3)/DSQRT(p3(1)**2+p3(2)**2+p3(3)**2)
      CosThe2 = p4(3)/DSQRT(p4(1)**2+p4(2)**2+p4(3)**2)
*--------------------------------------------------------------------
* Various definitions of Theta and s-propagator
*--------------------------------------------------------------------
** definition of P.L. B219, 103 (1989)
      CosThePL = ( p3(4)*CosThe1 -p4(4)*CosThe2)/(p3(4)+p4(4))
* definition of P.R. D41, 1425 (1990)
      SinThe1 = DSQRT(DABS((1d0-CosThe1)*(1d0+CosThe1)))
      SinThe2 = DSQRT(DABS((1d0-CosThe2)*(1d0+CosThe2)))
      yy1 = SinThe2/(SinThe1+SinThe2)
      yy2 = SinThe1/(SinThe1+SinThe2)
      CosPRD = yy1*CosThe1 - yy2*CosThe2
*-------------------------------
* LL formula for s'/s from angles according to ALEPH note 1996
      theta1= ACOS( CosThe1 )
      theta2= ACOS( CosThe2 )
      zAleph =  (SIN(theta1)+SIN(theta2) -ABS(SIN(theta1+theta2)))
     $         /(SIN(theta1)+SIN(theta2) +ABS(SIN(theta1+theta2)))
      sPropAl  = sIni*zAleph
      mPropAl  = DSQRT(ABS(sPropAl))
      END




      SUBROUTINE buker4(xFlag,Sign)
! observables for muons
! routine calculates whether for observable of the given index KK
! it is in or out the cuts.
! INPUT: HEPEVT COMMON block
! 
! OUTPUT: for every observable numbered KK CALL on storing routine is executed
!         (if in cuts) or is not executed ELSEwhere

      IMPLICIT REAL*8 (A-H,O-Z)
      DOUBLE PRECISION xFlag(100),Sign
      PARAMETER (NMXHEP=2000)
      COMMON/d_HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      REAL*8 phep,vhep,pol,sig,err
       SAVE  /d_HEPEVT/
      COMMON / INOUT / INUT,IOUT
      COMMON /kalinout/ wtkal(6)
      REAL*8 wtkal,wtsum(6),wtsum2(6),wto,wto2,wta,wta2,wt
      REAL*8  phsum(4),PH(4),PI,PHTR(4),qp(4),qm(4),xmp(10),xmm(10),ph1(4)
      LOGICAL IFacc,ifphot,ifpart,ifphot1
      CHARACTER*1 sorb
      CHARACTER*9 unitobsa,uniterra
      DATA pi /3.141592653589793238462643D0/
*-------------------------
      DO i=1,100
         xFlag(i)=0d0
      ENDDO

         cmsene=phep(4,1)*2
cc         wt=FHISTWT(ID)
          DO L=1,4
           QP(l)=PHEP(l,3)
           QM(l)=PHEP(l,4)
          ENDDO

         NEVTES=NEVTES+1
cc         CALL memorizer4(0,30,1D0,0D0,0D0,0D0)
C===============================================================
C     'ALEPH-12'
      kk=1
      NPH=0
      ITR=0
      IFacc=(xCOS(qp) .LT. 0.95 .AND. xCOS(qm) .LT. 0.95)
      IFacc=ifacc .AND. (enrg(qp) .GT. 5.0 .OR. enrg(qm) .GT. 5.0)
      IFacc=ifacc .AND. (xacol(qp,qm,3) .GT. 30.0)
      IFphot=.false.
      DO K=5,NHEP
         DO L=1,4
            PH(l)=PHEP(l,k)
         ENDDO
         IFpart=(XCOS(PH) .LT. 0.95)
         IFpart=ifpart .AND. (ENRG(PH) .GT. 15D0)
         IFpart=ifpart .AND. (xacol(PH,qm,3) .GT. 10.0) .AND. (xacol(PH,qp,3) .GT. 10.0)
         IFphot=ifphot .OR. ifpart
      ENDDO
      IFacc=ifacc .AND. ifphot .AND. ((qp(4)+QM(4)) .GT. 0.6*cmsene)
      IF (ifacc) THEN
cc          CALL memorizer4(0,KK,wt,0D0,0D0,0D0)
         xFlag(kk)=1d0
      ENDIF
C===============================================================
C     'ALEPH-15'
      kk=2
      NPH=0
      ITR=0
      IFacc=(xCOS(qp) .LT. 0.95 .AND. xCOS(qm) .LT. 0.95)
      IFacc=ifacc .AND. (enrg(qp) .GT. 5.0 .OR. enrg(qm) .GT. 5.0)
      IFacc=ifacc .AND. (xacol(qp,qm,3) .GT. 30.0)
      DO K=5,NHEP
         DO L=1,4
            PH(l)=PHEP(l,k)
         ENDDO
         IFpart=(XCOS(PH) .LT. 0.95)
         IFpart=ifpart .AND. (ENRG(PH) .GT. 15D0)
         IFpart=ifpart .AND. (xacol(PH,qm,3) .GT. 20.0) .AND. (xacol(PH,qp,3) .GT. 20.0)
         IF (ifpart) THEN
            nph=nph+1
            CALL pMsum(QP,PH,PHSUM)
            xmp(nph)=XMAS(phsum)
            CALL pMsum(QM,PH,PHSUM)
            xmm(nph)=XMAS(phsum)
         ENDIF
      ENDDO
      IFacc=ifacc .AND. (nph .GE. 2)
      IFpart=.false.
      DO k=1,nph
         DO l=1,nph
            IF (k .NE. l) THEN
               IFpart=ifpart .OR. (ABS(xmp(l)-xmm(k)) .LT. 5D0)
            ENDIF
         ENDDO
      ENDDO
      IF (ifacc .AND. ifpart) THEN
cc          CALL memorizer4(0,KK,wt,0D0,0D0,0D0)
         xFlag(kk)=1d0
      ENDIF
C===============================================================
C     'DELPHI9'
      kk=3
      NPH=0
      ITR=0
      esum=0
      ejetmax=0
      IFacc=(xCOS(qp) .LT. 0.9063 .AND. (enrg(qp) .GT. 5.0))
      IF (ifacc)                                   ejetmax=qp(4)
      IFacc=ifacc .OR. (xCOS(qm) .LT. 0.9063 .AND. (enrg(qm) .GT. 5.0))
      IF (xCOS(qm) .LT. 0.9063 .AND. (enrg(qm) .GT. 5.0)) ejetmax=max(ejetmax,qm(4))
      IF (xCOS(qp) .LT. 0.9397) esum=esum+qp(4)
      IF (xCOS(qm) .LT. 0.9397) esum=esum+qm(4)
      IFphot=.false.
      DO K=5,NHEP
         DO L=1,4
            PH(l)=PHEP(l,k)
         ENDDO
         IFpart=(XCOS(PH) .LT. 0.9848)
         IFpart=ifpart .AND. (ENRG(PH) .GT. 5D0)
         IF (ifpart) ejetmax=max(ejetmax,ph(4)) 
         IF (xCOS(PH) .LT. 0.9397) esum=esum+PH(4)
         IF(xacol(ph,qp,3) .LT. 15 .AND. qp(4) .GT. 1d0) IFpart=.false.
         IF(xacol(ph,qm,3) .LT. 15 .AND. qm(4) .GT. 1d0) IFpart=.false.
         erec=0
         DO i=5,NHEP
            DO j=1,4
               ph1(j)=phep(j,i)
            ENDDO
            IF (xacol(ph,ph1,3) .LT. 15 .AND. k .NE. i) erec=erec+ph1(4)
         ENDDO
         IFpart=ifpart .AND. erec .LT. 2
         IFphot=ifphot .OR. ifpart
      ENDDO
      IFacc=ifacc .AND. ifphot .AND. esum .GT. 0.2*cmsene .AND. ejetmax .GT. 10D0
      IF (ifacc) THEN
cc          CALL memorizer4(0,KK,wt,0D0,0D0,0D0)
         xFlag(kk)=1d0
      ENDIF
C===============================================================
C     'DELPHI12'
      kk=4
      NPH=0
      ITR=0
      esum=0
      IF (xCOS(qp) .LT. 0.9397) esum=esum+qp(4)
      IF (xCOS(qm) .LT. 0.9397) esum=esum+qm(4)
      ejetmax=0
      IFacc=(xCOS(qp) .LT. 0.9063 .AND. (enrg(qp) .GT. 5.0))
      IFacc=ifacc .AND. (xCOS(qm) .LT. 0.9063 .AND. (enrg(qm) .GT. 5.0))
      IFphot=.false.
      DO K=5,NHEP
         DO L=1,4
            PH(l)=PHEP(l,k)
         ENDDO
         IFpart=(XCOS(PH) .LT. 0.9848)
         IFpart=ifpart .AND. (ENRG(PH) .GT. 5D0)
         IF (ifpart) ejetmax=max(ejetmax,ph(4)) 
         IF (xCOS(PH) .LT. 0.9397) esum=esum+PH(4)
         IF(xacol(ph,qp,3) .LT. 15 .AND. qp(4) .GT. 1d0) IFpart=.false.
         IF(xacol(ph,qm,3) .LT. 15 .AND. qm(4) .GT. 1d0) IFpart=.false.
         erec=0
         DO i=5,NHEP
            DO j=1,4
               ph1(j)=phep(j,i)
            ENDDO
            IF (xacol(ph,ph1,3) .LT. 15 .AND. k .NE. i) erec=erec+ph1(4)
         ENDDO
         IFpart=ifpart .AND. erec .LT. 2
         IF (ifpart) THEN
            nph=nph+1
            CALL pMsum(QP,PH,PHSUM)
            xmp(nph)=XMAS(phsum)
            CALL pMsum(QM,PH,PHSUM)
            xmm(nph)=XMAS(phsum)
         ENDIF
      ENDDO
      IFacc=ifacc .AND. (nph .GE. 2)
      IFpart=.false.
      DO k=1,min(2,nph)
         DO l=1,min(2,nph)
            IF (k .NE. l) THEN
               IFpart=ifpart .OR. (ABS(xmp(l)-xmm(k)) .LT. 10D0)
            ENDIF
         ENDDO
      ENDDO
      IF (ifacc .AND. ifpart .AND. esum .GT. 0.2*cmsene) THEN
cc          CALL memorizer4(0,KK,wt,0D0,0D0,0D0)
         xFlag(kk)=1d0
      ENDIF
C===============================================================
C     'LT15'
      kk=5
      NPH=0
      ITR=0
      esum=0
      ejetmax=0
      IFacc=(xCOS(qp) .LT. 0.94)
      IFacc=ifacc .OR. (xCOS(qm) .LT. 0.94)
      IFphot=.false.
      DO K=5,NHEP
         DO L=1,4
            PH(l)=PHEP(l,k)
         ENDDO
         IFpart=(XCOS(PH) .LT. 0.97)
         IFpart=ifpart .AND. (ENRG(PH) .GT. 15D0)
         IF (ifpart) THEN
            nph=nph+1
            CALL pMsum(QP,PH,PHSUM)
            xmp(nph)=XMAS(phsum)
            IF(xCOS(qp) .GT. 0.94) xmp(nph)=0
            CALL pMsum(QM,PH,PHSUM)
            xmm(nph)=XMAS(phsum)
            IF(xCOS(qm) .GT. 0.94) xmm(nph)=0
         ENDIF
         IFphot=ifphot .OR. (ENRG(PH) .GT. 20D0 .AND. XCOS(PH) .LT. 0.75 .AND. max(xmp(nph),xmm(nph)) .GT. 70D0)
!          IFphot=ifphot .AND. ifpart
      ENDDO
      IFacc=ifacc .AND. ifphot .AND. nph .LE. 2
      IF (ifacc) THEN
cc          CALL memorizer4(0,KK,wt,0D0,0D0,0D0)
         xFlag(kk)=1d0
      ENDIF
c################################################################
C===============================================================
C     'LT18'
      kk=6
      NPH=0
      ITR=0
      esum=0
      ejetmax=0
      IFacc=(xCOS(qp) .LT. 0.94)
      IFacc=ifacc .AND. (xCOS(qm) .LT. 0.94)
      IFphot=.false.
      IFphot1=.true.
      DO K=5,NHEP
         DO L=1,4
            PH(l)=PHEP(l,k)
         ENDDO
         IFpart=(XCOS(PH) .LT. 0.97)
         IFpart=ifpart .AND. (ENRG(PH) .GT. 15D0)
         IF (ifpart) THEN
            nph=nph+1
            CALL pMsum(QP,PH,PHSUM)
            xmp(nph)=XMAS(phsum)
            CALL pMsum(QM,PH,PHSUM)
            xmm(nph)=XMAS(phsum)
         ENDIF
         IFphot=ifphot .OR. (ENRG(PH) .GT. 15D0 .AND. XCOS(PH) .LT. 0.75)
!          IFphot1=ifphot1 .AND. (ENRG(PH) .GT. 15D0 .AND. XCOS(PH) .LT. 0.94)
      ENDDO
      del1=ABS(xmp(1)-xmm(2))
      del2=ABS(xmm(1)-xmp(2))
      sum1=xmp(1)+xmm(2)
      sum2=xmp(2)+xmm(1)
      IFacc=ifacc .AND. ifphot .AND. ifphot1 .AND. nph .EQ. 2
      IFacc=ifacc .AND. ((del2 .LT. 10d0 .AND. sum2 .GT. 100) .OR. (del1 .LT. 10d0 .AND. sum1 .GT. 100))
      IF (ifacc) THEN
cc          CALL memorizer4(0,KK,wt,0D0,0D0,0D0)
         xFlag(kk)=1d0
      ENDIF
C===============================================================
C     'Opal14'
      kk=7
      NPH=0
      ITR=0
      evis=0
      egam=0
      IFacc=(xCOS(qp) .LT. 0.95 .AND. xCOS(qm) .LT. 0.95)
      IFacc=ifacc .AND. (pete(qp) .GT. 1.0 .AND. pete(qm) .GT. 1.0)
      IFacc=ifacc .AND. (xacol(qp,qm,3) .GT. 20.0)
      evis=qp(4)+qm(4)
      IFphot=.false.
      DO K=5,NHEP
         DO L=1,4
            PH(l)=PHEP(l,k)
         ENDDO
         IFpart=(XCOS(PH) .LT. 0.95)
         IFpart=ifpart .AND. (ENRG(PH) .GT. 0.05*cmsene/2)
         IFpart=ifpart .AND. (xacol(PH,qm,3) .GT. 20.0) .AND. (xacol(PH,qp,3) .GT. 20.0)
         IF (ifpart) egam=max(egam,ph(4))
         IFphot=ifphot .OR. ifpart
      ENDDO
      evis=evis+egam
      IFacc=ifacc .AND. ifphot .AND. (evis .GT. 1.6*cmsene/2)
      IF (ifacc) THEN
cc          CALL memorizer4(0,KK,wt,0D0,0D0,0D0)
         xFlag(kk)=1d0
      ENDIF
C===============================================================
C     'Opal15'
      kk=8
      NPH=0
      ITR=0
      evis=0
      egam=0
      IFacc=(xCOS(qp) .LT. 0.95 .AND. xCOS(qm) .LT. 0.95)
      IFacc=ifacc .AND. (pete(qp) .GT. 1.0 .AND. pete(qm) .GT. 1.0)
      IFacc=ifacc .AND. (xacol(qp,qm,3) .GT. 20.0)
      evis=qp(4)+qm(4)
      IFphot=.false.
      DO K=5,NHEP
         DO L=1,4
            PH(l)=PHEP(l,k)
         ENDDO
         IFpart=(XCOS(PH) .LT. 0.95)
         IFpart=ifpart .AND. (ENRG(PH) .GT. 0.05*cmsene/2)
         IFpart=ifpart .AND. (xacol(PH,qm,3) .GT. 20.0) .AND. (xacol(PH,qp,3) .GT. 20.0)
         IF (ifpart) egam=max(egam,ph(4))
         IFphot=ifphot .OR. ifpart
      ENDDO
      evis=evis+egam
      CALL pmsum(qp,qm,phsum)
      IFacc=ifacc .AND. (85D0 .GT. xmas(phsum) .OR. xmas(phsum) .GT. 95d0)
      IFacc=ifacc .AND. ifphot .AND. (evis .GT. 1.6*cmsene/2)
      IF (ifacc) THEN
cc          CALL memorizer4(0,KK,wt,0D0,0D0,0D0)
         xFlag(kk)=1d0
      ENDIF
C===============================================================
C     'Opal20'
      kk=9
      NPH=0
      ITR=0
      evis=0
      egam=0
      egam1=0
      IFacc=(xCOS(qp) .LT. 0.95 .AND. xCOS(qm) .LT. 0.95)
      IFacc=ifacc .AND. (pete(qp) .GT. 1.0 .AND. pete(qm) .GT. 1.0)
      IFacc=ifacc .AND. (xacol(qp,qm,3) .GT. 20.0)
      evis=qp(4)+qm(4)
      IFphot=.false.
      DO K=5,NHEP
         DO L=1,4
            PH(l)=PHEP(l,k)
         ENDDO
         IFpart=(XCOS(PH) .LT. 0.95)
         IFpart=ifpart .AND. (ENRG(PH) .GT. 0.05*cmsene/2)
         IFpart=ifpart .AND. (xacol(PH,qm,3) .GT. 20.0) .AND. (xacol(PH,qp,3) .GT. 20.0)
         IF (ifpart) egam=max(egam,ph(4))
         IF (ifpart) nph=nph+1
      ENDDO
      DO K=5,NHEP
         DO L=1,4
            PH(l)=PHEP(l,k)
         ENDDO
         IFpart=(XCOS(PH) .LT. 0.95)
         IFpart=ifpart .AND. (ENRG(PH) .GT. 0.05*cmsene/2) .AND. ph(4) .LT. egam
         IFpart=ifpart .AND. (xacol(PH,qm,3) .GT. 20.0) .AND. (xacol(PH,qp,3) .GT. 20.0)
         IF (ifpart) egam1=max(egam1,ph(4))
      ENDDO
      evis=evis+egam+egam1
      CALL pmsum(qp,qm,phsum)
      IFacc=ifacc .AND. nph .GT. 1 .AND. (evis .GT. 1.6*cmsene/2)
      IF (ifacc) THEN
cc          CALL memorizer4(0,KK,wt,0D0,0D0,0D0)
         xFlag(kk)=1d0
      ENDIF
C===============================================================
C     'Opal21'
      kk=10
      NPH=0
      ITR=0
      evis=0
      egam=0
      egam1=0
      IFacc=(xCOS(qp) .LT. 0.95 .AND. xCOS(qm) .LT. 0.95)
      IFacc=ifacc .AND. (pete(qp) .GT. 1.0 .AND. pete(qm) .GT. 1.0)
      IFacc=ifacc .AND. (xacol(qp,qm,3) .GT. 20.0)
      evis=qp(4)+qm(4)
      IFphot=.false.
      DO K=5,NHEP
         DO L=1,4
            PH(l)=PHEP(l,k)
         ENDDO
         IFpart=(XCOS(PH) .LT. 0.95)
         IFpart=ifpart .AND. (ENRG(PH) .GT. 0.05*cmsene/2)
         IFpart=ifpart .AND. (xacol(PH,qm,3) .GT. 20.0) .AND. (xacol(PH,qp,3) .GT. 20.0)
         IF (ifpart) egam=max(egam,ph(4))
         IF (ifpart) nph=nph+1
      ENDDO
      DO K=5,NHEP
         DO L=1,4
            PH(l)=PHEP(l,k)
         ENDDO
         IFpart=(XCOS(PH) .LT. 0.95)
         IFpart=ifpart .AND. (ENRG(PH) .GT. 0.05*cmsene/2) .AND. ph(4) .LT. egam
         IFpart=ifpart .AND. (xacol(PH,qm,3) .GT. 20.0) .AND. (xacol(PH,qp,3) .GT. 20.0)
         IF (ifpart) egam1=max(egam1,ph(4))
      ENDDO
      evis=evis+egam+egam1
      CALL pmsum(qp,qm,phsum)
      IFacc=ifacc .AND. (85D0 .GT. xmas(phsum) .OR. xmas(phsum) .GT. 95d0)
      IFacc=ifacc .AND. nph .GT. 1 .AND. (evis .GT. 1.6*cmsene/2)
      IF (ifacc) THEN
cc          CALL memorizer4(0,KK,wt,0D0,0D0,0D0)
         xFlag(kk)=1d0
      ENDIF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C general variables for observables of the mumu (eventually gamma)
      CALL pmsum(qp,qm,phsum)
      xminv=xmas(phsum)
      CALL xprop_view(1,phsum)
      xmprop=xmas(phsum)
      xmprop=xmas1(NHEP,PHEP)
      xmang=xmasang(cmsene,qp,qm)
C     xmang1=xmasang1(cmsene,qp,qm,PH)
      cosmu=xCOS(qm)
      sign=1
      IF(qm(3) .GT. 0d0) sign=-1
C-------------------------------------------------------------------------------------------
C-------------------------------------------------------------------------------------------
C-------------------------------------------------------------------------------------------
C-------------------------------------------------------------------------------------------
*     Realist muon sigma and AFB
C-------------------------------------------------------------------------------------------
C     Aleph5  !! need clarifications before checking what DOes it mean photon
      KK=11
      IFacc=(xCOS(qp) .LT. 0.95 .AND. xCOS(qm) .LT. 0.95)
      IFacc=ifacc .AND. (enrg(qp) .GT. 6.0 .AND. enrg(qm) .GT. 6.0)
      IFacc=ifacc .AND. (enrg(qp)+enrg(qm) .GT. 60D0)
      xmangu=xmang
      phmax=0
      DO K=5,NHEP
         DO L=1,4
            PH(l)=PHEP(l,k)
         ENDDO
         IF (xCOS(ph) .LT. 0.97 .AND. enrg(ph) .GT. 5d0 .AND. enrg(pH) .GT. phmax) THEN
            phmax=enrg(ph)
            xmangu=xmasang1(cmsene,qp,qm,PH)
         ENDIF
      ENDDO
      IF (ifacc .AND. xmangu .GT. 0.1*cmsene) THEN
cc          CALL memorizer4(0,KK,wt,0D0,wt*sign,0D0)
         xFlag(kk)=1d0
      ENDIF
C----------------------------------------------------
C Aleph6 !! need clarification
      KK=12
      IFacc=(xCOS(qp) .LT. 0.95 .AND. xCOS(qm) .LT. 0.95)
      IFacc=ifacc .AND. (enrg(qp) .GT. 6.0 .AND. enrg(qm) .GT. 6.0)
      IFacc=ifacc .AND. (enrg(qp)+enrg(qm) .GT. 60D0)
      xmangu=xmang
      phmax=0
      DO K=5,NHEP
         DO L=1,4
            PH(l)=PHEP(l,k)
         ENDDO
         IF (xCOS(ph) .LT. 0.97 .AND. enrg(ph) .GT. 5d0 .AND. enrg(pH) .GT. phmax) THEN
            phmax=enrg(ph)
            xmangu=xmasang1(cmsene,qp,qm,PH)
         ENDIF
      ENDDO
      IF (ifacc .AND. xmangu .GT. 0.9*cmsene .AND. xminv .GT. 0.74*cmsene) THEN
cc          CALL memorizer4(0,KK,wt,0D0,wt*sign,0D0)
         xFlag(kk)=1d0
      ENDIF
C----------------------------------------------------
C     Delphi4
      KK=13
      IFacc=(xCOS(qp) .LT. COS(20./180.*pi) .AND. xCOS(qm) .LT. COS(20./180.*pi))
      IFacc=ifacc .AND. (enrg(qp) .GT. 30.0 .OR. enrg(qm) .GT. 30.0)
      xmangu=xmang
      phmax=0
      IF (ifacc .AND. xminv .GT. 75) THEN
cc          CALL memorizer4(0,KK,wt,0D0,wt*sign,0D0)
         xFlag(kk)=1d0
      ENDIF
C----------------------------------------------------
C Delphi5
      KK=14
      IFacc=(xCOS(qp) .LT. COS(20./180.*pi) .AND. xCOS(qm) .LT. COS(20./180.*pi))
      IFacc=ifacc .AND. (enrg(qp) .GT. 30.0 .OR. enrg(qm) .GT. 30.0)
      xmangu=xmang
      phmax=0
      IF (ifacc .AND. xminv .GT. 0.85*cmsene) THEN
cc          CALL memorizer4(0,KK,wt,0D0,wt*sign,0D0)
         xFlag(kk)=1d0
      ENDIF
C----------------------------------------------------
C LT9
      KK=15
      IFacc=(xCOS(qp) .LT. COS(20./180.*pi) .AND. xCOS(qm) .LT. COS(20./180.*pi))
      IFacc=ifacc .AND. (enrg(qp) .GT. 35.0 .OR. enrg(qm) .GT. 35.0)
      xmangu=xmang
      phmax=0
      DO K=5,NHEP
         DO L=1,4
            PH(l)=PHEP(l,k)
         ENDDO
         IF (xCOS(ph) .LT. 0.985 .AND. enrg(ph) .GT. 15d0 .AND. enrg(pH) .GT. phmax
     $        .AND. (xacol(PH,qm,3) .GT. 10.0) .AND. (xacol(PH,qp,3) .GT. 10.)) THEN
            phmax=enrg(ph)
            xmangu=xmasang1(cmsene,qp,qm,PH)
         ENDIF
      ENDDO
      IF (ifacc .AND. xmangu .GT. 75D0) THEN
cc          CALL memorizer4(0,KK,wt,0D0,wt*sign,0D0)
         xFlag(kk)=1d0
      ENDIF
C----------------------------------------------------
C     LT10
      KK=16
      IFacc=(xCOS(qp) .LT. COS(20./180.*pi) .AND. xCOS(qm) .LT. COS(20./180.*pi))
      IFacc=ifacc .AND. (enrg(qp) .GT. 35.0 .OR. enrg(qm) .GT. 35.0)
      xmangu=xmang
      phmax=0
      DO K=5,NHEP
         DO L=1,4
            PH(l)=PHEP(l,k)
         ENDDO
         IF (xCOS(ph) .LT. 0.985 .AND. enrg(ph) .GT. 15d0 .AND. enrg(pH) .GT. phmax
     $        .AND. (xacol(PH,qm,3) .GT. 10.0) .AND. (xacol(PH,qp,3) .GT. 10.)) THEN
            phmax=enrg(ph)
            xmangu=xmasang1(cmsene,qp,qm,PH)
         ENDIF
      ENDDO
      IF (ifacc .AND. xmangu .GT. 0.85*cmsene) THEN
cc          CALL memorizer4(0,KK,wt,0D0,wt*sign,0D0)
         xFlag(kk)=1d0
      ENDIF
C----------------------------------------------------
C     Opal6
      KK=17
      IFacc=(xCOS(qp) .LT. 0.95 .AND. xCOS(qm) .LT. 0.95)
      IFacc=ifacc .AND. (enrg(qp) .GT. 6.0 .AND. enrg(qm) .GT. 6.0)
      IFacc=ifacc .AND. (xacol(qp,qm,2) .GT. 0.32/pi*180)
      xmangu=xmang
      phmax=0
      esum=qp(4)+qm(4)
      DO K=5,NHEP
         DO L=1,4
            PH(l)=PHEP(l,k)
         ENDDO
         IF (xCOS(ph) .LT. 0.985 .AND. enrg(ph) .GT. 0.8d0) THEN
!!!!!!!!!!!!!!!!             esum=esum+ph(4)
            IF      (enrg(pH) .GT. phmax
     $           .AND. (xacol(PH,qm,3) .GT. 0.2/pi*180)
     $           .AND. (xacol(PH,qp,3) .GT. 0.2/pi*180)) THEN
               xmangu=xmasang1(cmsene,qp,qm,PH)
            ENDIF
            IF (enrg(pH) .GT. phmax) THEN
               phmax=enrg(ph)
            ENDIF
            
         ENDIF
      ENDDO
      esum=esum+phmax
      IFacc=ifacc .AND. xmangu .GT. 0.1*cmsene
      IF (ifacc .AND. (
     $     (esum .GT. (0.35*cmsene+0.5*91.17**2/cmsene) .AND. esum .LT. (0.75*cmsene+0.5*91.17**2/cmsene) .AND. xminv .GT. 70D0)) 
     $     .OR. esum .LT. (0.35*cmsene+0.5*91.17**2/cmsene) .OR. esum .GT. (0.75*cmsene+0.5*91.17**2/cmsene)) THEN
cc          CALL memorizer4(0,KK,wt,0D0,wt*sign,0D0)
         xFlag(kk)=1d0
      ENDIF
C----------------------------------------------------
C     Opal7
      KK=18
      IFacc=(xCOS(qp) .LT. 0.95 .AND. xCOS(qm) .LT. 0.95)
      IFacc=ifacc .AND. (enrg(qp) .GT. 6.0 .AND. enrg(qm) .GT. 6.0)
      IFacc=ifacc .AND. (xacol(qp,qm,2) .GT. 0.32/pi*180)
      xmangu=xmang
      phmax=0
      esum=qp(4)+qm(4)
      DO K=5,NHEP
         DO L=1,4
            PH(l)=PHEP(l,k)
         ENDDO
         IF (xCOS(ph) .LT. 0.985 .AND. enrg(ph) .GT. 0.8d0) THEN
!!!!!!!!!!!!!!!             esum=esum+ph(4)
            IF      (enrg(pH) .GT. phmax
     $           .AND. (xacol(PH,qm,3) .GT. 0.2/pi*180)
     $           .AND. (xacol(PH,qp,3) .GT. 0.2/pi*180)) THEN
               xmangu=xmasang1(cmsene,qp,qm,PH)
            ENDIF
            IF (enrg(pH) .GT. phmax) THEN
               phmax=enrg(ph)
            ENDIF
         ENDIF
      ENDDO
      esum=esum+phmax
c        WRITE(*,*) '>>>>>>>>>>>>>>>>>>>>>>'
c        DO K=3,NHEP
c         WRITE(*,1000) PHEP(1,k),PHEP(2,k),PHEP(3,k),PHEP(4,k)
c        ENDDO
c 1000    FORMAT(1x,4(f10.4,2x))
C        WRITE(*,*) cmsene,'>>',xminv, xmang, xmangu
      IF (ifacc .AND. esum .GT. (0.35*cmsene+0.5*91.17**2/cmsene) 
     $     .AND. xmangu .GT. 0.85*cmsene .AND. xminv .GT. SQRT(0.1*cmsene**2+91.17**2) ) THEN
cc          CALL memorizer4(0,KK,wt,0D0,wt*sign,0D0)
         xFlag(kk)=1d0
      ENDIF
C----------------------------------------------------
C----------------------------------------------------
C----------------------------------------------------
C idealized observables !!!!!!!!!!!!!!!!!
C----------------------------------------------------
C IAleph5
        KK=19

        IF (cosmu .LT. 0.95D0 .AND. xminv .GT. 0.1*cmsene) THEN
cc          CALL memorizer4(0,KK,wt,0D0,wt*sign,0D0)
          xFlag(kk)=1d0
        ENDIF

C----------------------------------------------------
C IAleph6
        KK=20

        IF (cosmu .LT. 0.95D0 .AND. xminv .GT. 0.9*cmsene) THEN
cc          CALL memorizer4(0,KK,wt,0D0,wt*sign,0D0)
          xFlag(kk)=1d0
        ENDIF

C----------------------------------------------------
C IDelphi5
        KK=21

        IF (cosmu .LT. 0.95D0 .AND. xminv .GT. 75D0) THEN
cc          CALL memorizer4(0,KK,wt,0D0,wt*sign,0D0)
          xFlag(kk)=1d0
        ENDIF

C----------------------------------------------------
C IDelphi6
        KK=22

        IF (cosmu .LT. 0.95D0 .AND. xminv .GT. 0.85*cmsene) THEN
cc          CALL memorizer4(0,KK,wt,0D0,wt*sign,0D0)
          xFlag(kk)=1d0
        ENDIF

C----------------------------------------------------
C ILT9
        KK=23

        IF (cosmu .LT. 0.9D0 .AND. xminv .GT. 75) THEN
cc          CALL memorizer4(0,KK,wt,0D0,wt*sign,0D0)
          xFlag(kk)=1d0
        ENDIF

C----------------------------------------------------
C ILT10
        KK=24

        IF (cosmu .LT. 0.9D0 .AND. xminv .GT. 0.85*cmsene) THEN
cc          CALL memorizer4(0,KK,wt,0D0,wt*sign,0D0)
          xFlag(kk)=1d0
        ENDIF

C----------------------------------------------------
C ILT11
        KK=25

        IF (cosmu .LT. 1D0 .AND. xminv .GT. 0.85*cmsene) THEN
cc          CALL memorizer4(0,KK,wt,0D0,wt*sign,0D0)
          xFlag(kk)=1d0
        ENDIF

C----------------------------------------------------
C IOpal6
        KK=26

        IF (cosmu .LT. 0.95D0 .AND. xmprop .GT. 0.1*cmsene) THEN
cc          CALL memorizer4(0,KK,wt,0D0,wt*sign,0D0)
          xFlag(kk)=1d0
        ENDIF

C----------------------------------------------------
C IOpal7
        KK=27

        IF (cosmu .LT. 1D0 .AND. xmprop .GT. 0.1*cmsene) THEN
cc          CALL memorizer4(0,KK,wt,0D0,wt*sign,0D0)
          xFlag(kk)=1d0
        ENDIF

C----------------------------------------------------
C IOpal8
        KK=28

        IF (cosmu .LT. 0.95D0 .AND. xmprop .GT. 0.85*cmsene) THEN
cc          CALL memorizer4(0,KK,wt,0D0,wt*sign,0D0)
          xFlag(kk)=1d0
        ENDIF
C----------------------------------------------------
C IOpal9
        KK=29
        IF (cosmu .LT. 1D0 .AND. xmprop .GT. 0.85*cmsene) THEN
cc          CALL memorizer4(0,KK,wt,0D0,wt*sign,0D0)
          xFlag(kk)=1d0
        ENDIF
c################################################################
      END


C ================================================================
C ================================================================
C universal library for calculating observables for any generator
C storing final states in HEPEVT COMMON block
C in HEPEVT: FILEds 1,2 must be beams
C            fields 3,4 must be fs fermions.
C ================================================================
C ================================================================

      FUNCTION ENRG(PH)
c takes energy of the 4-vector
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PH(4)
      ENRG=PH(4)
      END

      FUNCTION XCOS(PH)
c for 4-vector calculates module of the cosine with beam 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PH(4)
      XCOS=ABS(PH(3))/SQRT(PH(1)**2+PH(2)**2+PH(3)**2)
      END

      FUNCTION PETE(PH)
c for 4-vector calculates p_t
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PH(4)
      PETE=SQRT(PH(1)**2+PH(2)**2)
      END

      FUNCTION XMAS(PH)
c for 4-vector calculates mass
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PH(4)
      XMAS=SQRT(ABS(PH(4)**2-PH(3)**2-PH(2)**2-PH(1)**2))
      XMAS1=XMAS
      END
      FUNCTION XMAS1(N,PHEP)
c calculate mass of virtual Z assuming leading log approx.
c event must have beams first fs fermions later and photons at END
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NMXHEP=2000)
      REAL*8 PHEP(5,NMXHEP)
      DIMENSION qp(4),qm(4),psum(4),pp(4),pm(4),ph(4)
      AMZ=91.17
      GAMZ=2.5
      DO K=1,4
       PP(k)=PHEP(K,1)
       PM(k)=PHEP(K,2)
       QP(k)=PHEP(K,3)
       QM(k)=PHEP(K,4)
       PSUM(K)=QP(k)+QM(K)
      ENDDO
      DO L=5,N
        DO K=1,4
          PH(k)=PHEP(K,1)
        ENDDO
        XM1=XINV(PP,PH)
        XM1=MAX(XM1,XINV(PM,PH))
        XM2=XINV(QP,PH)
        XM2=MAX(XM2,XINV(QM,PH))
        XMV1=(XMAS(PSUM))**2
        XMV2=XINV(PSUM,PH)
        XM1=1D0/XM1/((XMV1-AMZ**2)**2+GAMZ**2*AMZ**2)
        XM2=1D0/XM2/((XMV2-AMZ**2)**2+GAMZ**2*AMZ**2)
        IF (XM2 .GT. XM1) THEN
          DO K=1,4
           PSUM(K)=PSUM(K)+PH(K)
          ENDDO
        ENDIF
      ENDDO
      XMAS1=XMAS(PSUM)
      END

      FUNCTION XINV(qp,qm)
C invariant calculated from 2 four-momenta
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION qp(4),qm(4)
      XINV=(qp(4)+qm(4))**2-(qp(3)+qm(3))**2-(qp(2)+qm(2))**2-(qp(1)+qm(1))**2
      XINV=ABS(XINV)
      END

      FUNCTION XMASANG(cmsene,qp,qm)
C mass calculated from directions of qp qm assuming there is just 
c one extra collinear with beam ph.
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION qp(4),qm(4)
      c1=qp(3)/SQRT(qp(1)**2+qp(2)**2+qp(3)**2)
      c2=qm(3)/SQRT(qm(1)**2+qm(2)**2+qm(3)**2)
      s1=SQRT(1-c1**2)
      s2=SQRT(1-c2**2)
      t1=s1/c1
      t2=s2/c2
      tt=ABS(t1+t2)
      xk=2*tt/(tt+SQRT(1+1/t1**2)+SQRT(1+1/t2**2))
      R=ABS(c1 +c2*s1/s2)/(1+s1/s2)
      xk=2*r/(1+r)
      xmasang=cmsene*SQRT(1-xk)
      END
      FUNCTION XMASANG1(cmsene,qp,qm,PH)
C mass calculated from directions of qp qm assuming there is just one extra visible ph.
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION qp(4),qm(4),ph(4)
      c1=(qp(1)*PH(1)+qp(2)*PH(2)+qp(3)*PH(3))
     $   /SQRT(qp(1)**2+qp(2)**2+qp(3)**2)/SQRT(PH(1)**2+PH(2)**2+PH(3)**2)
      c2=(qm(1)*PH(1)+qm(2)*PH(2)+qm(3)*PH(3))
     $   /SQRT(qm(1)**2+qm(2)**2+qm(3)**2)/SQRT(PH(1)**2+PH(2)**2+PH(3)**2)
      s1=SQRT(1-c1**2)
      s2=SQRT(1-c2**2)
      t1=s1/c1
      t2=s2/c2
      tt=ABS(t1+t2)
      xk=2*tt/(tt+SQRT(1+1/t1**2)+SQRT(1+1/t2**2))
      R=ABS(c1 +c2*s1/s2)/(1+s1/s2)
      xk=2*r/(1+r)
      xmasang1=cmsene*SQRT(ABS(1-xk))
      END

      SUBROUTINE PMSUM(P,Q,PH)
c adds two 4-vectors
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PH(4),P(4),Q(4)
      DO K=1,4
       PH(k)=P(k)+Q(K)
      ENDDO
      END

      FUNCTION XACOL(X,Y,N)
C     ********************
c calculates acollinearity 
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8    X(*),Y(*)
      DIMENSION X1(4),Y1(4)
      DATA PI /3.1415926535897932D0/
      S=0.D0
      X2=0.D0
      Y2=0.D0
      DO 9  I=1,N
      X1(I)=X(I)
    9 Y1(I)=Y(I)
      DO 10 I=1,N
      S=S+X1(I)*Y1(I)
      X2=X2+X1(I)**2
   10 Y2=Y2+Y1(I)**2
      XACOL=ACOS(S/SQRT(X2*Y2))*180.D0/PI
      RETURN
      END
      SUBROUTINE xprop_view(mode,PP)
      IMPLICIT REAL*8 (A-H,O-Z)
c this routine is digging out of yfs3 non-obervable 4-momentum of Z.
c but can also get it from hepevt, in LL approximation
c event must have beams first fs fermions later and photons at END
c IF you DO not use guts of yfs subr. karlud comment out 1 line.
      PARAMETER (NMXHEP=2000)
      COMMON/d_HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      REAL*8 phep,vhep
      REAL*8 pol,sig,err
      SAVE  /HEPEVT/
      REAL*8 PP(4),PS(4)
      SAVE ps
      DIMENSION qp(4),qm(4),psum(4),pm(4),ph(4)
      AMZ=91.17
      GAMZ=2.5


      IF (mode .EQ. 0) THEN
        DO K=1,4
           ps(k)=pp(k)
        ENDDO
        RETURN
      ENDIF

      DO K=1,4
       PP(k)=PHEP(K,1)
       PM(k)=PHEP(K,2)
       QP(k)=PHEP(K,3)
       QM(k)=PHEP(K,4)
       PSUM(K)=QP(k)+QM(K)
       PP(k)=PSUM(K)
      ENDDO
      DO L=5,NHEP
        DO K=1,4
          PH(k)=PHEP(K,L)
        ENDDO
        XM1=XINV(PP,PH)
        XM1=MIN(XM1,XINV(PM,PH))
        XM2=XINV(QP,PH)
        XM2=MIN(XM2,XINV(QM,PH))
        XMV1=(XMAS(PSUM))**2
        XMV2=XINV(PSUM,PH)
        XM1=1D0/XM1/((XMV1-AMZ**2)**2+GAMZ**2*AMZ**2)
        XM2=1D0/XM2/((XMV2-AMZ**2)**2+GAMZ**2*AMZ**2)
        IF (XM2 .GT. XM1) THEN
          DO K=1,4
           PSUM(K)=PSUM(K)+PH(K)
           PP(k)=PSUM(K)
          ENDDO
        ENDIF
      ENDDO

C case when it could be directly taken from host PROGRAM (in mode 0)
C overWRITEs calculated in LL pp. IF it is not available
C next 3 lines should be commented out.
      DO k=1,4
        pp(k)=ps(k)
      ENDDO

      END
