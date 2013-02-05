
      SUBROUTINE SAB_MakeRen
*////////////////////////////////////////////////////////////////////////////////
* Reading hist-file from MC
* Define entry KFf=7 corresponding to sum over all quarks
* Define AFB
*////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'PlotAll.h'
      INCLUDE '../RobAll.h'
      INTEGER           ikFf,iThe,iSel,iMod
      DOUBLE PRECISION  xTot,eTot,xFB, eFB,xAFB,eAFB,sum1,sum2
      DOUBLE PRECISION  CMSene
      INTEGER  Init
      DATA     Init/0/
      SAVE     Init
      IF(Init.NE.0) RETURN
      Init=1
*
      CMSene = m_xpar(1)
      m_Energy = '??????'
      IF( ABS(CMSene-189d0).LT.001) m_Energy = '189GeV'
      IF( ABS(CMSene-200d0).LT.001) m_Energy = '200GeV'
      IF( ABS(CMSene-206d0).LT.001) m_Energy = '206GeV'
      WRITE(*,*) 'SAB_MakeRen: m_Energy =',m_Energy
*
      CALL GLK_RenHst("NB  ", m_IdGen, m_kSAB, m_kSABren)
      CALL GLK_Delet(m_kSAB)
*---------- Define entry KFf=7 corresponding to sum over all quarks
      DO iThe =  m_iThe1, m_iThe2
         DO iSel =  m_iSel1, m_iSel2
            DO iMod = m_iMod1, m_iMod2
               sum1=0d0
               sum2=0d0
               DO iKFf =  1, 5
                  CALL  SAB_GetBin(m_kSABren,ikFf,iSel,iThe,iMod, xTot,eTot) !
                  sum1=sum1+xTot
                  sum2=sum2+eTot
               ENDDO
               CALL  SAB_SetBin(m_kSABren,7,iSel,iThe,iMod, sum1,sum2) ! KFf=7 contains the sum
            ENDDO
         ENDDO
      ENDDO
*---------- Next define AFB
      DO iKFf =  m_iKFf1, m_iKFf2
         DO iSel =  m_iSel1, m_iSel2
            DO iMod = m_iMod1, m_iMod2
               CALL  SAB_GetBin(m_kSABren,ikFf,iSel,m_iThe1,iMod, xTot,eTot) !
               CALL  SAB_GetBin(m_kSABren,ikFf,iSel,m_iThe2,iMod, xFB, eFB) !
               xAFB  =0d0
               eAFB  =0d0
               IF(xTOT.NE.0d0) THEN
                  xAFB  = xFB/xTot
                  eAFB  = SQRT( (xFB*eTot/xTot**2)**2 +(eFB/xTot)**2) !
               ENDIF
               CALL  SAB_SetBin(m_kSABren,ikFf,iSel,m_iThe2,iMod, xAFB,eAFB) !
            ENDDO
         ENDDO
      ENDDO
      END

      SUBROUTINE SAB_FilKKsem
*////////////////////////////////////////////////////////////////////////////////
* For the moment it fills xsections only, AFB will be also needed!
*////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'PlotAll.h'
      INCLUDE '../RobAll.h'
      INTEGER           i,j,k,ikFf,iThe,iSel,iMod,idum
      INTEGER           KeyFSR, KeyDis
      DOUBLE PRECISION  xec,afb,eRror,sum1,xer
      CHARACTER*5       chType
      INTEGER  Init
      DATA     Init/0/
      SAVE     Init
      IF(Init.NE.0) RETURN
      Init=1
*-----------------------
      CALL GLK_Book1(m_kSABsem,'Sigma and AFB semianalytical $', m_nBinSAB, 0d0,1d0*m_nBinSAB) !
*---------- Define semianalytical entry iMod = 10
      CALL KKsem_GetKeyFSR(KeyFSR)
      CALL GLK_Book1(i_Dummy,'dummy $', 100, 0d0,1d0) !
c[[[[[    test of IFI soft limit
      IF( KeyFSR .NE. 0) THEN
         chType = "VRHO "       ! test of IFI soft limit, not itegration over bin
         KeyDis =  90
         CALL KKsem_VVplot(KeyDis,chType, i_Misc1,'test1 of IFI soft limit $',i_Dummy) !
         KeyDis =  91
         CALL KKsem_VVplot(KeyDis,chType, i_Misc2,'test2 of IFI soft limit $',i_Dummy) !
      ENDIF
c]]]]]
      IF( KeyFSR .NE. 0) THEN
         KeyDis =  303302       ! O(alf2)LL
         KeyDis =  304302       ! ISR with O(alf3)LL exp(3/4*gam)
         KeyDis =  305302       ! The best ISR,  O(alf3)LL + O(alf2)NLL
         chType = "XCHI2"       ! ISR+FSR
      ELSE
         KeyDis=   301          ! 301 is 1-st order ISR
         KeyDis=   305          ! The best ISR,  O(alf3)LL + O(alf2)NLL
         chType = "VCHI2"       ! FSR only
      ENDIF
      eRror=0d0
      iMod = 10                 !<------- This is the Model entry for SemiAnalyt.
      DO iKFf =  m_iKFf1, m_iKFf2
         IF(iKFf.LE.5 .OR. (iKFf.GE.12 .AND. iKFf.LE.16 )) THEN ! Quarks, muon, tau, neutrina
            IF( NINT(m_xpar(400+iKFf)).NE.0) THEN           ! Only generated entries
               CALL BornV_SetKF(iKFf)
               WRITE(*,*) '*> SAB_FilKKsem: iKFf =',iKFf
c[[[[[ Simple Born of KKsem
***            CALL KKsem_SetKeyFoB(-10) ! KeyFoB=-10, Simple Born of KKsem
***            CALL BornV_SetKeyQCD(0)
c]]]]]
               IF( m_FlagSem .EQ. 1 ) THEN
c((((((((((((((( without AFB older version
***            CALL KKsem_SetKeyFoB( 0) ! back to normal
***            CALL KKsem_VVplot(KeyDis,chType, ivx_Semi,' O(alf3) ini+fin best Backward$',i_Dummy) !
c)))))))))))))))
                  CALL KKsem_SetKeyFoB(+1) ! Forvard
                  CALL KKsem_VVplot(KeyDis,chType, ivx_SemiF,' O(alf3) ini+fin best Forward$',i_Dummy) !
                  CALL KKsem_SetKeyFoB(-1) ! backward
                  CALL KKsem_VVplot(KeyDis,chType, ivx_SemiB,' O(alf3) ini+fin best Backward$',i_Dummy) !
                  CALL GLK_Operat(ivx_SemiF,  '+', ivx_SemiB ,  ivx_Semi,  1d0, 1d0) ! F+B
                  CALL GLK_Operat(ivx_SemiF,  '-', ivx_SemiB,   iva_Semi,  1d0, 1d0) !
                  CALL GLK_Operat(iva_Semi,   '/', ivx_Semi ,   iva_Semi,  1d0, 1d0) !
*******        CALL GLK_Operat(ivx_Semi,   '+', ivx_Semi,    ivx_Semi,  1d3, 0d0) ! --> picobarns
               ELSE
                  WRITE(*,*) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
                  WRITE(*,*) '$$$$$ KKsem disabled! $$$$$$$$$'
                  WRITE(*,*) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
                  DO i=1,30000000
                     k=k*i
                  ENDDO
                  CALL GLK_Clone1(i_Dummy,ivx_Semi,' KKsem disabled $')
                  CALL GLK_Clone1(i_Dummy,iva_Semi,' KKsem disabled $')
               ENDIF
               CALL Sab_AngTranslate(ivx_Semi,ix_Semi)
               CALL Sab_AngTranslate(iva_Semi,ia_Semi)
               DO iSel =  1, m_nVmax
                  iThe=1
                  CALL GLK_GetBin(ix_Semi,iSel,xec)
                  CALL SAB_SetBin(m_kSABsem, ikFf,iSel,iThe,iMod, xec,0d0) !
                  iThe=2
                  CALL GLK_GetBin(ia_Semi,iSel,afb)
                  CALL SAB_SetBin(m_kSABsem, ikFf,iSel,iThe,iMod, afb,0d0) !
               ENDDO !iSel
*(((      CALL GLK_SetNout(6)
               CALL GLK_Print(ivx_Semi)
               CALL GLK_Print( ix_Semi)
*)))      CALL GLK_SetNout(16)
               CALL GLK_Delet(ivx_SemiF)
               CALL GLK_Delet(ivx_SemiB)
               CALL GLK_Delet(ivx_Semi)
               CALL GLK_Delet(iva_Semi)
               CALL GLK_Delet(ix_Semi)
               CALL GLK_Delet(ia_Semi)
            ENDIF !iKFf
         ENDIF !iKFf
      ENDDO
*---------- Define entry KFf=7 corresponding to sum over all quarks
      CALL SAB_sumQks(m_kSABsem,iMod)
      CALL GLK_Delet(i_Dummy)
      END


      SUBROUTINE SAB_GetBin(iDB,ikFf,iSel,iThe,iMod, xbin,ebin)
*////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'PlotAll.h'
      INCLUDE '../RobAll.h'
      INTEGER           iDB,ikFf,iSel,iThe,iMod,iBin
      DOUBLE PRECISION  xbin,ebin
      iBin= (iKFf-m_iKFf1  )*(m_iSel2-m_iSel1+1)*(m_iThe2-m_iThe1+1)*(m_iMod2-m_iMod1+1)  !
     $                        +(iSel -m_iSel1  )*(m_iThe2-m_iThe1+1)*(m_iMod2-m_iMod1+1)  !
     $                                            +(iThe -m_iThe1  )*(m_iMod2-m_iMod1+1)  !
     $                                                                +(iMod -m_iMod1  )+1!
      IF( iBin .LT.1 .OR. iBin.GT.m_nBinSAB ) THEN
         WRITE(*,*) ' ----STOP in SAB_GetBin: iBin,m_nBinSAB=',iBin,m_nBinSAB !
         WRITE(*,*) ' ----ikFf,iSel,iThe,iMod=', ikFf,iSel,iThe,iMod !
         STOP
      ENDIF
      CALL  GLK_GetBin(iDB,iBin,xbin)
      CALL  GLK_GetErr(iDB,iBin,ebin)
      END

      SUBROUTINE SAB_SetBin(iDB,ikFf,iSel,iThe,iMod, xbin,ebin)
*////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'PlotAll.h'
      INCLUDE '../RobAll.h'
      INTEGER           iDB,ikFf,iSel,iThe,iMod,iBin
      DOUBLE PRECISION  xbin,ebin
      iBin= (iKFf-m_iKFf1  )*(m_iSel2-m_iSel1+1)*(m_iThe2-m_iThe1+1)*(m_iMod2-m_iMod1+1)  !
     $                        +(iSel -m_iSel1  )*(m_iThe2-m_iThe1+1)*(m_iMod2-m_iMod1+1)  !
     $                                            +(iThe -m_iThe1  )*(m_iMod2-m_iMod1+1)  !
     $                                                                +(iMod -m_iMod1  )+1!
      IF( iBin .LT.1 .OR. iBin.GT.m_nBinSAB ) THEN
         WRITE(*,*) ' ----STOP in SAB_SetBin: iBin,m_nBinSAB=',iBin,m_nBinSAB !
         WRITE(*,*) ' ----ikFf,iSel,iThe,iMod=', ikFf,iSel,iThe,iMod !
         STOP
      ENDIF
      CALL  GLK_SetBin(iDB,iBin,xbin)
      CALL  GLK_SetErr(iDB,iBin,ebin)
      END


      SUBROUTINE SAB_GetKFHist(iDB,iKF1,iKF2, iThe,iMod,iSel, Fact, id)
*////////////////////////////////////////////////////////////////////////////////
* Get histogram in flavour
      IMPLICIT NONE
      INCLUDE 'PlotAll.h'
      INCLUDE '../RobAll.h'
      INTEGER           iDB,ikFf,iThe,iMod,iSel1,iSel2,id
      INTEGER           isel,iBin,icont,nSel,nkFf,ikF1,ikF2
      DOUBLE PRECISION  Fact, xbin(500),ebin(500),xx,ee
      LOGICAL           GLK_Exist
*----------
      IF(GLK_Exist(id)) THEN
        WRITE(*,*) " ++++ KKsem_VVplot: warning deleted id= ",id
        CALL GLK_Delet(id)
      ENDIF
      nKFf = iKF2-iKF1+1
      CALL GLK_Book1(id,'KF codes $', nKFf, 0d0,1d0*nKFf) !
      icont=0
      DO iKFf = iKF1,iKF2
         icont=icont+1
         CALL SAB_GetBin(iDB,iKFf,iSel,iThe,iMod, xx, ee)
         xbin(icont) = xx*Fact
         ebin(icont) = ee*Fact
cc         write(*,*) 'SAB_GetHiSel-->  iThe,iSel, xbin(iSel), ebin(iSel)', iThe,iSel, xbin(iSel), ebin(iSel) !
      ENDDO
      CALL GLK_Pak(  id, xbin)
      CALL GLK_Pake( id, ebin)
      CALL GLK_idopt(id,'ERRO')
      END

      SUBROUTINE SAB_GetHistSel(iDB,ikFf,iThe,iMod, iSel1,iSel2, Fact, id) !
*////////////////////////////////////////////////////////////////////////////////
* Get histogram in selection
      IMPLICIT NONE
      INCLUDE 'PlotAll.h'
      INCLUDE '../RobAll.h'
      INTEGER           iDB,ikFf,iThe,iMod,iSel1,iSel2,id
      INTEGER           isel,iBin,icont,nSel
      DOUBLE PRECISION  Fact, xbin(500),ebin(500),xx,ee
      LOGICAL           GLK_Exist
*--------------------
      IF(GLK_Exist(id)) THEN
        WRITE(*,*) " ++++ KKsem_VVplot: warning deleted id= ",id
        CALL GLK_Delet(id)
      ENDIF      
      nSel= iSel2-iSel1+1
      CALL GLK_Book1(id,'Selections $', nSel, 0d0,1d0*nSel) !
      icont=0
      DO iSel=iSel1,iSel2
         icont=icont+1
         CALL SAB_GetBin(iDB,ikFf,iSel,iThe,iMod, xx, ee)
         xbin(icont) = xx*Fact
         ebin(icont) = ee*Fact
****     WRITE(*,*) 'SAB_GetHiSel--> iThe,iSel, xbin(iSel), ebin(iSel)', iThe,iSel, xbin(iSel), ebin(iSel) !
      ENDDO
      CALL GLK_Pak(  id, xbin)
      CALL GLK_Pake( id, ebin)
      CALL GLK_idopt(id,'ERRO')
      END

      SUBROUTINE SAB_GetHistV(iDB,ikFf,iThe,iMod, Fact, id) !
*////////////////////////////////////////////////////////////////////////////////
* Get histogram in selection
      IMPLICIT NONE
      INCLUDE 'PlotAll.h'
      INCLUDE '../RobAll.h'
      INTEGER           iDB,ikFf,iThe,iMod,id,iVmax,iSel
      DOUBLE PRECISION  Fact
      LOGICAL           GLK_Exist
      DOUBLE PRECISION  xbin(m_nbV),ebin(m_nbV),xx,ee,vmax
      INTEGER           i
*----------
*     WRITE(*,*) 'SAB_GetHistV: entered ============'
      IF(GLK_Exist(id)) THEN
        WRITE(*,*) " ++++ KKsem_VVplot: warning deleted id= ",id
        CALL GLK_Delet(id)
      ENDIF
      CALL GLK_Book1(id,'V-distribution $', m_nbV, m_loV,m_upV) !
      DO i=1, m_nbV
         xbin(i) =-10d0**(Mod(id,99))
         ebin(i) = 0d0
         vmax= m_loV +((m_upV-m_loV)/m_nbV)*i
         iSel= 0
         DO iVmax=1,m_nVmax
            IF( ABS(vmax-m_VmaxList(iVmax)) .LT.1d-4) iSel=iVmax
         ENDDO
         IF(iSel.NE.0) THEN
            CALL SAB_GetBin(iDB,ikFf,iSel,iThe,iMod, xx, ee)
****        WRITE(*,*) 'SAB_GetHistV:ikFf,iSel,iThe,iMod,xx, ee=',ikFf,iSel,iThe,iMod,xx, ee !
            xbin(i) = xx*fact
            ebin(i) = ee*fact
         ENDIF
      ENDDO
      CALL GLK_Pak(  id, xbin)
      CALL GLK_Pake( id, ebin)
      CALL GLK_idopt(id,'ERRO')
*     WRITE(*,*) 'SAB_GetHistV: exited ============'
      END

      SUBROUTINE Sab_AngTranslate(id1,id2)
*////////////////////////////////////////////////////////////////////////////////
* Translate from 100 bins to vMax bining
      IMPLICIT NONE
      INCLUDE 'PlotAll.h'
      INCLUDE '../RobAll.h'
      INTEGER           id1,id2,iVmax,iBin,j
      LOGICAL           GLK_Exist
      DOUBLE PRECISION  xbin(m_nVmax),ebin(m_nVmax)
*
      IF(GLK_Exist(id2)) THEN
        WRITE(*,*) " ++++ Sab_AngTranslate: warning deleted id= ",id2
        CALL GLK_Delet(id2)
      ENDIF      
      CALL GLK_Book1(id2,'From translation $', m_nVmax, 0d0,1d0*m_nVmax) !
      DO iVmax=1, m_nVmax
         j=NINT(100d0*m_VmaxList(iVmax))
cc         write(*,*) iVmax,j
         CALL GLK_GetBin(id1,     j, xbin(iVmax))
         CALL GLK_GetErr(id1,     j, ebin(iVmax))
      ENDDO
      CALL GLK_Pak(  id2, xbin)
      CALL GLK_Pake( id2, ebin)
      CALL GLK_idopt(id2,'ERRO')
      END

      SUBROUTINE   SAB_RDatZF1(DataFile1,iDB)
*/////////////////////////////////////////////////////////////////////////////
*// reading sigma and  AFB with error from diskfile, ZFITTER  results       //
*/////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'PlotAll.h'
      INCLUDE '../RobAll.h'
      CHARACTER*(*)     DataFile1
      INTEGER           iDB
      CHARACTER*4       cSAB
      CHARACTER*33      ch2
      INTEGER           line, indf, ikFf,iSel,iThe,iMod,iVmax
      DOUBLE PRECISION  ene,vmax,xs1,xs2,xs3,CMSene
      SAVE iSel
*-------------
      OPEN(20,File=DataFile1)
* read header
      WRITE(*,*) '========== SAB_RDatZF1 Enter ================'
      DO line=1,10000
         READ(20,'(a,a,I2)') cSAB
         IF(cSAB .EQ. '<be>' ) GOTO 200
      ENDDO
      WRITE(*,*) ' STOP in SAB_RDatZF1, end of data file',DataFile1
      STOP
 200  CONTINUE
* read lines
      DO line=1,10000
         READ(20,'(a,a,I2,f7.1,f8.4,3f10.4)') cSAB,ch2,indf,ene,vmax,xs1,xs2,xs3 !
         IF(cSAB .EQ. '<en>') GOTO 900
* Translation iVmax->iSel depends on table of Vmax !!!!
         iSel= 0
         DO iVmax=1,m_nVmax
            IF( ABS(vmax-m_VmaxList(iVmax)) .LT.1d-4) iSel=iVmax
         ENDDO
         iThe=1
         IF(cSAB .EQ. 'A_FB') iThe=2
         ikFf = m_KFfromINDF(indf)
         CMSene = m_xpar(1)
         IF( (ABS(CMSene-ene).LT.1d-3) .AND. (iSel.NE.0) .AND. (ikFf.NE.-1) ) THEN !
         WRITE(*,'(a,a,I2,f7.1,f8.4,3f10.4)') cSAB,ch2,indf,ene,vmax,xs1,xs2,xs3 !
***         WRITE(*,'(a,10i5)') 'indf,ikFf,iSel,iThe,iMod =',indf,ikFf,iSel,iThe,iMod !
            iMod = 2    ! IFIoff
            CALL SAB_SetBin(iDB, ikFf,iSel,iThe,iMod, xs3,0d0) !
            iMod = 1    ! IFI-exp
            CALL SAB_SetBin(iDB, ikFf,iSel,iThe,iMod, xs1,0d0) !
            iMod = 3    ! IFI-nexp
            CALL SAB_SetBin(iDB, ikFf,iSel,iThe,iMod, xs2,0d0) !
         ENDIF
      ENDDO
 900  CONTINUE
      CLOSE(20)
      WRITE(*,*) '========== SAB_RDatZF1 Exit  ================'
      END

      SUBROUTINE   SAB_RDatZF1a(DataFile1,iDB)
*/////////////////////////////////////////////////////////////////////////////
*// reading sigma and  AFB with error from diskfile, ZFITTER  results       //
*/////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'PlotAll.h'
      INCLUDE '../RobAll.h'
      CHARACTER*(*)     DataFile1
      INTEGER           iDB
      CHARACTER*4       cSAB
      CHARACTER*10      ch1
      CHARACTER*2       ch2
      CHARACTER*16      ch3
      INTEGER           line, indf, ikFf,iSel,iThe,iMod,iVmax
      DOUBLE PRECISION  ene,vmax,xs1,xs2,xs3,CMSene
      SAVE iSel
*-------------
      OPEN(20,File=DataFile1)
* read header
      WRITE(*,*) '========== SAB_RDatZF1a Enter ================'
      DO line=1,10000
         READ(20,'(a)') cSAB
         IF(cSAB .EQ. '<be>' ) GOTO 200
      ENDDO
      WRITE(*,*) ' STOP in SAB_RDatZF1, end of data file',DataFile1
      STOP
 200  CONTINUE
* read lines
      DO line=1,10000
         READ(20,'(a,a,a,a,I2,f7.1,f7.3,3f14.7)') cSAB,ch1,ch2,ch3,indf,ene,vmax,xs1,xs2,xs3 !
****     WRITE(*,'(a,a,a,a,I2,f7.1,f7.3,3f14.7)') cSAB,ch1,ch2,ch3,indf,ene,vmax,xs1,xs2,xs3 !
         IF(cSAB .EQ. '<en>') GOTO 900
* Translation iVmax->iSel depends on table of Vmax !!!!
         iSel= 0
         DO iVmax=1,m_nVmax
            IF( ABS(vmax-m_VmaxList(iVmax)) .LT.1d-4) iSel=iVmax
         ENDDO
         iThe=1
         IF(ch2 .EQ. 'AS') iThe=2
         ikFf = m_KFfromINDF(indf)
         CMSene = m_xpar(1)
         IF( (ABS(CMSene-ene).LT.1d-3) .AND. (iSel.NE.0) .AND. (ikFf.NE.-1) ) THEN !
****        WRITE(*,'(a,a,a,a,I2,f7.1,f7.3,3f14.7)') cSAB,'--------->',ch2,ch3,indf,ene,vmax,xs1,xs2,xs3 !
            iMod = 2    ! ZF IFIoff
            CALL SAB_SetBin(iDB, ikFf,iSel,iThe,iMod, xs3,0d0) ! off
            iMod = 1    ! ZF IFI-exp
            CALL SAB_SetBin(iDB, ikFf,iSel,iThe,iMod, xs1,0d0) ! exp
            iMod = 3    ! ZF IFI-nexp
            CALL SAB_SetBin(iDB, ikFf,iSel,iThe,iMod, xs2,0d0) ! on
         ENDIF
      ENDDO
 900  CONTINUE
      CLOSE(20)
      WRITE(*,*) '========== SAB_RDatZF1a Exit  ================'
      END

      SUBROUTINE   SAB_RDatZF1b(DataFile1,iDB)
*/////////////////////////////////////////////////////////////////////////////
*// reading sigma and  AFB with error from diskfile, ZFITTER 
*/////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'PlotAll.h'
      INCLUDE '../RobAll.h'
      CHARACTER*(*)     DataFile1
      INTEGER           iDB
      CHARACTER*4       cSAB
      CHARACTER*33      ch2
      INTEGER           line, indf, ikFf,iSel,iThe,iMod,iVmax
      DOUBLE PRECISION  ene,vmax,xs1,xs2,xs3,CMSene
      SAVE iSel
*-------------
      OPEN(20,File=DataFile1)
* read header
      WRITE(*,*) '========== SAB_RDatZF1b Enter ================'
      DO line=1,10000
         READ(20,'(a,a,I2)') cSAB
         IF(cSAB .EQ. '<be>' ) GOTO 200
      ENDDO
      WRITE(*,*) ' STOP in SAB_RDatZF1, end of data file',DataFile1
      STOP
 200  CONTINUE
* read lines
      DO line=1,10000
         READ(20,'(a,a,I2,f7.1,f8.4,3f12.4)') cSAB,ch2,indf,ene,vmax,xs1,xs2,xs3 !
         IF(cSAB .EQ. '<en>') GOTO 900
* Translation iVmax->iSel depends on table of Vmax !!!!
         iSel= 0
         DO iVmax=1,m_nVmax
            IF( ABS(vmax-m_VmaxList(iVmax)) .LT.1d-4) iSel=iVmax
         ENDDO
         iThe=1
         IF(cSAB .EQ. 'A_FB') iThe=2
         ikFf = m_KFfromINDF(indf)
         CMSene = m_xpar(1)
         IF( (ABS(CMSene-ene).LT.1d-3) .AND. (iSel.NE.0) .AND. (ikFf.NE.-1) ) THEN !
         WRITE(*,'(a,a,I2,f7.1,f8.4,3f12.4)') cSAB,ch2,indf,ene,vmax,xs1,xs2,xs3 !
***         WRITE(*,'(a,10i5)') 'indf,ikFf,iSel,iThe,iMod =',indf,ikFf,iSel,iThe,iMod !
            iMod = 2    ! IFIoff
            CALL SAB_SetBin(iDB, ikFf,iSel,iThe,iMod, xs3,0d0) !
            iMod = 1    ! IFI-exp
            CALL SAB_SetBin(iDB, ikFf,iSel,iThe,iMod, xs1,0d0) !
            iMod = 3    ! IFI-nexp
            CALL SAB_SetBin(iDB, ikFf,iSel,iThe,iMod, xs2,0d0) !
         ENDIF
      ENDDO
 900  CONTINUE
      CLOSE(20)
      WRITE(*,*) '========== SAB_RDatZF1b Exit  ================'
      END

      SUBROUTINE   SAB_RDatZF2(DataFile1,iDB,iMod)
*/////////////////////////////////////////////////////////////////////////////
*// reading sigma and  AFB with error from diskfile, KORALZ results         //
*/////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'PlotAll.h'
      INCLUDE '../RobAll.h'
      CHARACTER*(*)     DataFile1
      INTEGER           iDB,iMod
      CHARACTER*4       cSAB
      CHARACTER*17      ch2
      INTEGER           line, indf, ikFf,iSel,iThe,iVmax
      DOUBLE PRECISION  ene,vmax,xs1,xs2,xs3,as1,as2,as3,CMSene
      SAVE iSel
*-------------
      OPEN(20,File=DataFile1)
* read header
      WRITE(*,*) '========== SAB_RDatZF2 Enter ================'
      DO line=1,10000
         READ(20,'(a,a,I2)') cSAB
         IF(cSAB .EQ. '<be>' ) GOTO 200
      ENDDO
      WRITE(*,*) ' STOP in SAB_RDatZF2, end of data file',DataFile1
      STOP
 200  CONTINUE
* read lines
      DO line=1,10000
*********************************
         READ(20,'(   a,   a,   I2,  f7.1, f7.4, 9f14.4)') 
     $             cSAB, ch2, indf,   ene, vmax, xs3, as3 !
         WRITE(*,'(   a,   a,   I2,  f7.1, f7.4, 9f14.4)') 
     $             cSAB, ch2, indf,   ene, vmax, xs3, as3 !
*********************************
         IF(cSAB .EQ. '<en>') GOTO 900
* Translation iVmax->iSel depends on table of Vmax !!!!
         iSel= 0
         DO iVmax=1,m_nVmax
            IF( ABS(vmax-m_VmaxList(iVmax)) .LT.1d-4) iSel=iVmax
         ENDDO
         iThe=1
         ikFf = m_KFfromINDF(indf)
         CMSene = m_xpar(1)
         IF( (ABS(CMSene-ene).LT.1d-3) .AND. (iSel.NE.0) .AND. (ikFf.NE.-1) ) THEN !
            CALL SAB_SetBin(iDB, ikFf,iSel,iThe,iMod, xs3,0d0) !
         ENDIF
      ENDDO
 900  CONTINUE
      CLOSE(20)
      WRITE(*,*) '========== SAB_RDatZF2 Exit  ================'
      END

      SUBROUTINE   SAB_RDatZF2b(DataFile1,iDB,iMod,iBox)
*/////////////////////////////////////////////////////////////////////////////
*// reading sigma and  AFB with error from diskfile, KORALZ results         //
*/////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'PlotAll.h'
      INCLUDE '../RobAll.h'
      CHARACTER*(*)     DataFile1
      INTEGER           iDB,iMod,iBox
      CHARACTER*4       cSAB
      CHARACTER*23      ch1
      CHARACTER*2       ch2
      INTEGER           line, indf, ikFf,iSel,iThe,iVmax,iBx
      DOUBLE PRECISION  ene,vmax,xs1,xs2,xs3,as1,as2,as3,CMSene
      SAVE iSel
*-------------
      OPEN(20,File=DataFile1)
* read header
      WRITE(*,*) '========== SAB_RDatZF2b Enter ================'
      DO line=1,10000
         READ(20,'(a,a,I2)') cSAB
         IF(cSAB .EQ. '<be>' ) GOTO 200
      ENDDO
      WRITE(*,*) ' STOP in SAB_RDatZF2b, end of data file',DataFile1
      STOP
 200  CONTINUE
* read lines
      DO line=1,10000
*********************************
         READ(20,'(   a,   a,  I1,   a,   I2,  f7.1, f7.4, 9f14.4)') 
     $             cSAB, ch1, iBx, ch2, indf,   ene, vmax, xs3, as3 !
c         WRITE(*,'(   a,   a,  I1,   a,   I2,  f7.1, f7.4, 9f14.4)') 
c     $             cSAB, ch1, iBx, ch2, indf,   ene, vmax, xs3, as3 !
*********************************
         IF(cSAB .EQ. '<en>') GOTO 900
* Translation iVmax->iSel depends on table of Vmax !!!!
         iSel= 0
         DO iVmax=1,m_nVmax
            IF( ABS(vmax-m_VmaxList(iVmax)) .LT.1d-4) iSel=iVmax
         ENDDO
         ikFf = m_KFfromINDF(indf)
         CMSene = m_xpar(1)
         IF( (ABS(CMSene-ene).LT.1d-3) 
     $        .AND. (iSel.NE.0) 
     $        .AND. (iBx.EQ.iBox)
     $        .AND. (ikFf.NE.-9) ) THEN !
            iThe=1
            CALL SAB_SetBin(iDB, ikFf,iSel,iThe,iMod, xs3,0d0) !
            iThe=2
            CALL SAB_SetBin(iDB, ikFf,iSel,iThe,iMod, as3,0d0) !
         WRITE(*,'(   a,   a,   I1,   a,   I2,  f7.1, f7.4, 9f14.4)') 
     $             cSAB, ch1, iBox, ch2, indf,   ene, vmax, xs3, as3 !
         ENDIF
      ENDDO
 900  CONTINUE
      CALL  SAB_sumQks(iDB,iMod)
      CLOSE(20)
      WRITE(*,*) '========== SAB_RDatZF2b Exit  ================'
      END


      SUBROUTINE   SAB_RDatZF2c(DataFile1,iDB,iMod,iBox)
*/////////////////////////////////////////////////////////////////////////////
*// reading sigma and  AFB with error from diskfile, ZFITTER results        //
*/////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'PlotAll.h'
      INCLUDE '../RobAll.h'
      CHARACTER*(*)     DataFile1
      INTEGER           iDB,iMod,iBox
      CHARACTER*4       cSAB
      CHARACTER*17      ch1
      CHARACTER*2       ch2
      INTEGER           line, indf, ikFf,iSel,iThe,iVmax,iBx
      DOUBLE PRECISION  ene,vmax,xs1,xs2,xs3,as1,as2,as3,CMSene
      SAVE iSel
*-------------
      OPEN(20,FILE=DataFile1)
* read header
      WRITE(*,*) '========== SAB_RDatZF2c Enter ================'
      DO line=1,10000
         READ(20,'(a,a,I2)') cSAB
         IF(cSAB .EQ. '<be>' ) GOTO 200
      ENDDO
      WRITE(*,*) ' STOP in SAB_RDatZF2b, end of data file',DataFile1
      STOP
 200  CONTINUE
* read lines
      DO line=1,10000
*********************************
         READ(20,'(   a,   a,   I3,  f8.1, f7.4, 9f14.4)') 
     $             cSAB, ch1, indf,   ene, vmax, xs3, as3 !
c         WRITE(*,'(   a,   a,   I3,  f8.1, f7.4, 9f14.4)') 
c     $             cSAB, ch1, indf,   ene, vmax, xs3, as3 !
*********************************
         IF(cSAB .EQ. '<en>') GOTO 900
* Translation iVmax->iSel depends on table of Vmax !!!!
         iSel= 0
         DO iVmax=1,m_nVmax
            IF( ABS(vmax-m_VmaxList(iVmax)) .LT.1d-4) iSel=iVmax
         ENDDO
         ikFf = m_KFfromINDF(indf)
         CMSene = m_xpar(1)
         IF( (ABS(CMSene-ene).LT.1d-3) 
     $        .AND. (iSel.NE.0) 
cc     $        .AND. (iBx.EQ.iBox)
     $        .AND. (ikFf.NE.-9) ) THEN !
            iThe=1
            CALL SAB_SetBin(iDB, ikFf,iSel,iThe,iMod, xs3,0d0) !
            iThe=2
            CALL SAB_SetBin(iDB, ikFf,iSel,iThe,iMod, as3,0d0) !
         WRITE(*,'(   a,   a,   I3,  f8.1, f7.4, 9f14.8)') 
     $             cSAB, ch1, indf,   ene, vmax, xs3, as3 !
         ENDIF
      ENDDO
 900  CONTINUE
      CALL  SAB_sumQks(iDB,iMod)
      CLOSE(20)
      WRITE(*,*) '========== SAB_RDatZF2c Exit  ================'
      END


      SUBROUTINE   SAB_sumQks(iDB,iMod)
*/////////////////////////////////////////////////////////////////////////////
*// combines xsections and afb for quarks
*/////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'PlotAll.h'
      INCLUDE '../RobAll.h'
      INTEGER           iDB,iMod,iSel,iKFf,iThe
      DOUBLE PRECISION  xec,xer,afb,aer,sum1,sum2,zum1,zum2
*---------- Define entry KFf=7 corresponding to sum over all quarks
* for the moment xsection entry only!!!
      DO iSel =  1, m_nVmax
         sum1=0d0
         sum2=0d0
         zum1=0d0
         zum2=0d0
         DO iKFf =  1,5         ! Quarks only
            IF( NINT(m_xpar(400+iKFf)).NE.0) THEN           ! Only generated entries
               iThe=1
               CALL  SAB_GetBin(iDB,ikFf,iSel,iThe,iMod, xec,xer)
               sum1=sum1+xec
               zum1=zum1+xer**2
               iThe=2
               CALL  SAB_GetBin(iDB,ikFf,iSel,iThe,iMod, afb,aer)
****           WRITE(*,*) '=============',iDB,ikFf,iSel,iThe,iMod, xec,xer !
               sum2=sum2+xec*afb
               zum2=zum2+(xec*aer)**2+(xer*afb)**2
            ENDIF !iKFf
         ENDDO
         iThe=1
         xer = SQRT(zum1)
         CALL  SAB_SetBin(m_kSABsem,7,iSel,iThe,iMod, sum1,xer) ! KFf=7 contains the sum
         iThe=2
         afb =0d0
         aer =0d0
         IF(sum1.NE.0d0) afb = sum2/sum1
         IF(sum1.NE.0d0) aer = SQRT(zum2)/sum1 !<- this is slight underestimate
         CALL  SAB_SetBin(m_kSABsem,7,iSel,iThe,iMod,  afb,aer) ! KFf=7 contains the sum
         WRITE(*,*) '///////iMod,iSel,xec,xer; afb,aer=',imod,iSel,sum1,xer, afb,aer !
      ENDDO
      END



