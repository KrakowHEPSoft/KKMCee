      PROGRAM MAIN
C     ************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common / cglib  / blan(20000)
      COMMON / INOUT  / NINP,NOUT
      COMMON / DNAME / DNAME
      CHARACTER*6 DNAME
                                       
      CALL GLIMIT(20000)                      
      NINP= 5 
      NOUT=16
      CALL GOUTPU(NOUT)
      OPEN( NOUT,file="output/dallas.output")
C ----restore histograms from BHLUMI
      NINPH=20
      OPEN(NINPH,file="input_a/ftn70")
      CALL GRFILE(NINPH,' ',' ')  
      CALL GRIN(   0,9999,0)
C ----restore second group of histograms from OLDBIS
!      NINPH=21
!      OPEN(NINPH,file="input_b/dump.hst")
!      CALL GRFILE(NINPH,' ',' ')  
!      CALL GRIN(   0,9999,0)                       
      CALL GPRINT(7030)
      CALL GPRINT(7031)
      CALL GPRINT(7032)
      CALL GPRINT(7033)
      CALL GPRINT(7034)
      CALL GPRINT(7035)
      fact=0
      do k=0,2
       CALL RENHIE(7030+2*k  ,FACT,50)
       CALL RENHIE(7030+2*k+1,FACT,4)
      enddo      
      CALL GOPERA(7030,'/',7032,7036,1.D0,1.D0)   
      CALL GOPERA(7034,'/',7030,7037,1.D0,1.D0)   
      CALL GOPERA(7034,'/',7032,7038,1.D0,1.D0)   
C ----initialize GPLOT
      CALL GPLINT(0)
      NOUFIG=11
      OPEN(NOUFIG,file="output/dallas.tex")
      rewind(noufig)
      CALL GPLCAP(-NOUFIG)
      N=10
C======================================
      CALL GPLTIT(' PLOTs of acopl ME3  $')
      YMAX=3.5d-2
      CALL GMAXIM(7030,YMAX)
      CALL GMAXIM(7032,YMAX)
      CALL GMAXIM(7034,YMAX)
      CALL GMINIM(7030,0D0)
      CALL GMINIM(7032,0D0)
      CALL GMINIM(7034,0D0)
      CALL GPLOT (7030,' ',' ',0)
      CALL GPLTIT(' PLOTs of acopl YFS  $')
      CALL GPLOT (7032,' ',' ',0)
      CALL GPLTIT(' PLOTs of acopl ME2  $')
      CALL GPLOT (7034,' ',' ',0)

      CALL GPLTIT('acopl ME2 vs ME3  $')
      CALL GPLOT (7030,' ',' ',0)
      CALL GPLOT (7034,'S',' ',0)

      CALL GPLTIT('acopl ME3 vs YFS  $')
      CALL GPLOT (7030,' ',' ',0)
      CALL GPLOT (7032,'S',' ',0)

      CALL GMAXIM(7036,1.5d0)
      CALL GMINIM(7036,0.5D0)
      CALL GMAXIM(7037,1.5d0)
      CALL GMINIM(7037,0.5D0)
      CALL GMAXIM(7038,1.5d0)
      CALL GMINIM(7038,0.5D0)

      CALL GPLTIT('acopl ME3 ovr YFS  $')
      CALL GPLOT (7036,' ',' ',0)
      CALL GPLTIT('acopl ME2 ovr YFS  $')
      CALL GPLOT (7038,' ',' ',0)
      CALL GPLTIT('acopl ME2 ovr ME3  $')
      CALL GPLOT (7037,' ',' ',0)
C
C======================================
      CALL GPLTIT(' PLOTs of acopl ME3  $')
      YMAX=3.5d-1
      CALL GMAXIM(7031,YMAX)
      CALL GMAXIM(7033,YMAX)
      CALL GMAXIM(7035,YMAX)
      CALL GMINIM(7031,0D0)
      CALL GMINIM(7033,0D0)
      CALL GMINIM(7035,0D0)
      CALL GPLOT (7031,' ',' ',0)
      CALL GPLTIT(' PLOTs of acopl YFS  $')
      CALL GPLOT (7033,' ',' ',0)
      CALL GPLTIT(' PLOTs of acopl ME2  $')
      CALL GPLOT (7035,' ',' ',0)
C
C======================================
C
      CALL GPLEND                          

      END
      SUBROUTINE RENHIE(ID,FACT,NB)
C     ****************************
C errors taken into account
C     INTRODUCES HISTOGRAM NORMALISATION
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION X(400),XE(400)
      FAC=FACT
      CALL GUNPAK(ID, X,'    ',1)
      CALL GUNPAK(ID,XE,'ERRO',1)
      CALL GRESET(ID,' ')
      SUM=0D0
      DO 10 I=1,NB
   10 SUM=SUM+X(I)
      IF(SUM.EQ.0D0) RETURN
      IF(FAC.EQ.0D0) FAC=1D0/SUM
      DO 20 I=1,NB
      X(I) = X(I)*FAC
      XE(I)=XE(I)*FAC
   20 CONTINUE
      CALL GPAK( ID,X )
      CALL GPAKE(ID,XE)
      
      END
