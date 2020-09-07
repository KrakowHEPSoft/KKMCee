      SUBROUTINE BornV_StartEW(xpar_input)
*///////////////////////////////////////////////////////////////////
*//                                                               //
*//   Jump-start EW correction library                            //
*//   Here it is done by reading forfactors from disk file        //
*//                                                               //
*///////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
*------------------------------------------------------------------
      DOUBLE PRECISION  xpar_input(*)
      CALL BornV_ReadAll
      END                       !BornV_StartEW

      SUBROUTINE BornV_ReadAll
*///////////////////////////////////////////////////////////////////
*//       Reading from disk  all pretabulations                   //
*//       All quarks, muon and tau                                //
*///////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BornV.h'
      INTEGER    KFdown, KFup, KFstran, KFcharm, KFbotom, KFtop
      PARAMETER( KFdown  = 1,   KFup    = 2,
     $           KFstran = 3,   KFcharm = 4,
     $           KFbotom = 5,   KFtop   = 6)
      INTEGER    KFel,KFelnu,KFmu,KFmunu,KFtau,KFtaunu
      PARAMETER( KFel    = 11,  KFelnu  = 12,
     $           KFmu    = 13,  KFmunu  = 14,
     $           KFtau   = 15,  KFtaunu = 16)
      CHARACTER*40 TableFile
*----------------------------------------------------------------------
      TableFile= '../../dizet/table.down'
      CALL BornV_ReadFile(TableFile,KFdown)
      TableFile= '../../dizet/table.up'
      CALL BornV_ReadFile(TableFile,KFup)
      TableFile= '../../dizet/table.down'
      CALL BornV_ReadFile(TableFile,KFstran)
      TableFile= '../../dizet/table.up'
      CALL BornV_ReadFile(TableFile,KFcharm)
      TableFile= '../../dizet/table.botom'
      CALL BornV_ReadFile(TableFile,KFbotom)
      TableFile= '../../dizet/table.mu'
      CALL BornV_ReadFile(TableFile,KFmu)
      TableFile= '../../dizet/table.tau'
      CALL BornV_ReadFile(TableFile,KFtau)
      TableFile= '../../dizet/table.nue'
      CALL BornV_ReadFile(TableFile,KFelnu)
      TableFile= '../../dizet/table.numu'
      CALL BornV_ReadFile(TableFile,KFmunu)
      TableFile= '../../dizet/table.nutau'
      CALL BornV_ReadFile(TableFile,KFtaunu)

      END


      SUBROUTINE BornV_ReadFile(DiskFile,KFfin)
*///////////////////////////////////////////////////////////////////
*//            Reading from disk  single file                     //
*///////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'BXformat.h'
      INCLUDE 'BornV.h'
      CHARACTER*(*) DiskFile
      INTEGER KFfin, KFf
      CHARACTER*1 chr
      INTEGER i,j,k,ndisk,n,n1,n2
      DOUBLE PRECISION   amz, amh, amtop
      DOUBLE PRECISION   ww,cosi,gammw
*----------------------------------------------------------------------
      KFf = ABS(KFfin)
      IF( KFf.LT.1 .OR. KFf.GT.16 ) GOTO 900
      WRITE(6 ,*) 'Tables are READ from DiskFile  ',DiskFile
      ndisk=21
      OPEN(ndisk,FILE=DiskFile)
* header
      READ(ndisk,m_fmt0)  amz, amh, amtop, m_swsq, m_gammz, m_MW, m_GammW !
      WRITE(*,'(a)') 'amz,amh,amtop,swsq,gammz,amw,gammw= '
      WRITE(*,'(a,10f12.7)') '     =',amz, amh, amtop, m_swsq, m_gammz, m_MW, m_GammW
* basic range
      DO i=0, m_poin1
         READ(ndisk,m_fmt1) chr,n,ww
         READ(ndisk,m_fmt2) (m_cyy(i+1,k,KFf),k=1,m_poinG) ! EW
         READ(ndisk,m_fmt2) (m_syy(i+1,k,KFf),k=1,m_poinQ) ! QCD
      ENDDO
* Z pole range
      DO i=0,m_poin2
         DO  j=0,m_poTh2
            READ(ndisk,m_fmt1) chr,n,ww
            READ(ndisk,m_fmt2) (m_czz(i+1,j+1,k,KFf),k=1,m_poinG) ! EW
         ENDDO
         READ(ndisk,m_fmt2) (m_szz(i+1,k,KFf),k=1,m_poinQ) ! QCD
      ENDDO
* LEP2 range
      DO  i=0,m_poin3
         DO  j=0,m_poTh3
            READ(ndisk,m_fmt1) chr,n1,ww,n2,cosi
            READ(ndisk,m_fmt2) (m_ctt(i+1,j+1,k,KFf),k=1,m_poinG) ! EW
         ENDDO
         READ(ndisk,m_fmt2) (m_stt(i+1,k,KFf),k=1,m_poinQ) ! QCD
      ENDDO
* NLC range
      DO  i=0,m_poin4
         DO  j=0,m_poTh4
            READ(ndisk,m_fmt1) chr,n1,ww,n2,cosi
            READ(ndisk,m_fmt2) (m_clc(i+1,j+1,k,KFf),k=1,m_poinG) ! EW
         ENDDO
         READ(ndisk,m_fmt2) (m_slc(i+1,k,KFf),k=1,m_poinQ) ! QCD
      ENDDO
      CLOSE(ndisk)
      WRITE(m_out,bxope)
      WRITE(m_out,bxtxt) 'BornV  Reading from disk file:         '
      WRITE(m_out,bxtxt) DiskFile
      WRITE(m_out,bxl1f)   amz   ,   'Z mass             ','amz   ','a1'
      WRITE(m_out,bxl1f)   amh   ,   'Higgs mass         ','amh   ','a2'
      WRITE(m_out,bxl1f)   amtop ,   'Top mass           ','amtop ','a3'
      WRITE(m_out,bxl1f) m_swsq  ,   'sin**2(thetaW)     ','swsq  ','a3'
      WRITE(m_out,bxl1f) m_gammz ,   'Z width            ','gammz ','a3'
      WRITE(m_out,bxl1f) m_MW    ,   'W mass             ','amw   ','a3'
      WRITE(m_out,bxl1f) m_GammW ,   'W width            ','gammw ','a3'
      WRITE(m_out,bxclo)
* Check if tables on disk were produced with the same input params!!!
      IF( amz.EQ.0d0  .OR. amh.EQ.0d0  .OR. amtop.EQ.0d0  .OR.
     $     ABS(1-m_MZ/amz)      .GT. 1d-5  .OR.
     $     ABS(1-m_amh/amh)     .GT. 1d-5  .OR.
     $     ABS(1-m_amtop/amtop) .GT. 1d-5  ) THEN
         WRITE(    *,*) '+++ STOP in BornV_ReadFile: CHECK mz,mh,mtop'
         WRITE(m_out,*) '+++ STOP in BornV_ReadFile: CHECK mz,mh,mtop'
         WRITE(    *,*) 'm_MZ/amz, m_amh/amh, m_amtop/amtop =',m_MZ/amz,m_amh/amh,m_amtop/amtop!
         STOP
      ENDIF
      RETURN
 900  WRITE(m_out,*) '+++ BornV_ReadFile: wrong KFf= ', KFf
      WRITE(    6,*) '+++ BornV_ReadFile: wrong KFf= ', KFf
      END
