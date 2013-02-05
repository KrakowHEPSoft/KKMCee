      PROGRAM Z0_test
C     *****************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON  /    / BLAN(100000) 
      COMMON / INOUT  / NINP,NOUT,NOUT2
      COMMON / RANPAR / KEYRND
      DIMENSION XPAR(99),NPAR(99)

      CALL GLIMIT(100000)

      NOUT =46
      NOUT2=45

      NOUT3=99

c....input parameters>>>
C------option on running program :
C               KEYOPT=1  running PYTHIA56
C               KEYOPT=2  running simple MC for t-tbar 
C                       (to test compatibility with PYTHIA56)
C               KEYOPT=3  running option for Z-->l l +  1 photon 
C               KEYOPT=4  running option for Z-->l l +  2 photons 
C               KEYOPT=5  running simple MC for Z--> l l  + PHOTOS_20 final
      KEYOPT = 4
C final state multiplicity.
      KEYDIM  = 2+2
C----- random generator
      KEYRND = 1
C----- bremsstrahlung option
C          NPAR(1)=1 initial state
C          NPAR(1)=2 initial/final state
C          NPAR(1)=3 final state
      NPAR(1)=2
C----- Z0/gamma mixing switched on/off 
C          NPAR(2)=1 Z0/gamma
C          NPAR(2)=2 gamma only
C          NPAR(2)=3 Z0 only
C          NPAR(2)=4 no axial coupling
C          NPAR(2)=5 no vector coupling
C          NPAR(2)=6 only Z0 axial coupling
C          NPAR(2)=7 only Z0 vector coupling
      NPAR(2)=2
C------proton-proton CMS energy in GeV
      XPAR(1)= 176D0
C------initial quark  mass in GeV
      XPAR(2)= 0.511D-3  
C------final lepton mass in GeV
      XPAR(3)= 0.511D-3  
c.....trigerr parameters>>>low level trigger
C------minim photon energy in laboratory  frame
      XPAR(5)= 0.01D0*XPAR(1)/2D0
C------maximum photon energy in laboratory  frame
      XPAR(6)= 0.990D0*XPAR(1)/2D0 !0.99D0*XPAR(1)/2D0
c.....trigerr parameters>>>physical trigger
c------minimu photon transverse momenta in laboratory frame
      XPAR(7)=   10D0
c------minimu lepton transverse momenta in laboratory frame
      XPAR(8)=   2D0
c------minimum lepton-lepton mass
      XPAR(9)=   2D0
 
c------GSW parameters
       XPAR(50)  =   91.177D0
       XPAR(51)  =   2.4786D0
       XPAR(52)  =   0.232D0


C------number of requested events
      IF(KEYOPT.EQ.1) NPAR(10)=       100
      IF(KEYOPT.EQ.2) NPAR(10)=     10 
      IF(KEYOPT.EQ.3) NPAR(10)=   1 000
      IF(KEYOPT.EQ.4) NPAR(10)=   40 00!      1 000
      IF(KEYOPT.EQ.5) NPAR(10)=   400 000
      IF(KEYOPT.EQ.6) NPAR(10)=      1000
 
c-------number of repetation of calling PHOTOS from one gevent
      NPAR(11) = 100
C-----output ident.
      IF(KEYOPT.EQ.1) THEN
         NOUT =36
         NOUT2=35
         NOUTH=30
      ELSEIF(KEYOPT.EQ.2) THEN
         NOUT =56
         NOUT2=55
         NOUTH=50
      ELSEIF(KEYOPT.EQ.3) THEN
         NOUT =66
         NOUT2=65
         NOUTH=60
      ELSEIF(KEYOPT.EQ.4) THEN
         NOUT =76
         NOUT2=75
         NOUTH=70
      ELSEIF(KEYOPT.EQ.5) THEN
         NOUT =86
         NOUT2=85
         NOUTH=80
      ELSEIF(KEYOPT.EQ.6) THEN
         NOUT =96
         NOUT2=95
         NOUTH=90
       ENDIF
C------initialize histo output
      CALL GOUTPU(NOUT)
C....initialization of main routines
      IF(KEYOPT.EQ.1) THEN
c        CALL ROBOL0(-1)   
      ELSEIF(KEYOPT.EQ.2) THEN
C**        CALL Z0DEC0(-1,XPAR,NPAR)
C**        CALL BOKER1(-1,XPAR,NPAR)
c        CALL BOKER2(-1,XPAR,NPAR)
C**        CALL BOKER3(-1,XPAR,NPAR)
C**      ELSEIF(KEYOPT.EQ.3) THEN
C**        CALL Z0DEC1(-1,XPAR,NPAR)
C**        CALL BOK1PH(-1,XPAR,NPAR)
C**        CALL BOKER6(-1,XPAR,NPAR)
      ELSEIF(KEYOPT.EQ.4) THEN
        if (KEYDIM.eq.5) then
         CALL Z0DEC3(-1,XPAR,NPAR)
         CALL BOK3PH(-1,XPAR,NPAR)
        else
         CALL Z0DEC2(-1,XPAR,NPAR)
         CALL BOK2PH(-1,XPAR,NPAR)
        endif

c        CALL BOKER7(-1,XPAR,NPAR)
ccc        CALL TRZYNUE(-1,XPAR,NPAR)
C**      ELSEIF(KEYOPT.EQ.5) THEN
C**        CALL Z0DEC0(-1,XPAR,NPAR)
ccc        CALL BOKER2(-1,XPAR,NPAR)
C**        CALL Z0RAD0(-1,XPAR,NPAR)
C**      ELSEIF(KEYOPT.EQ.6) THEN
C**        CALL Z0DEC1(-1,XPAR,NPAR)
ccc        CALL BOKER2(-1,XPAR,NPAR)
C**        CALL Z0RAD1(-1,XPAR,NPAR)
      ENDIF
c.... initialization of PYTHIA56
c      CALL PREPYT(XPAR,NPAR)
C
C....generating mode
      DO 10 IEV=1,NPAR(10)
      IF(KEYOPT.EQ.1) THEN
      if(mod(iev,500).eq.1) write(6,*)    'event no=',iev
c         CALL PYEVNT 
c         CALL LUHEPC(1)
c         CALL HEPLUJ
c         CALL ROBOL0(0)
      ELSEIF(KEYOPT.EQ.2) THEN
      if(mod(iev,50 000).eq.1) write(6,*)    'event no=',iev
C**         CALL Z0DEC0(0,XPAR,NPAR)      
C**         CALL BOKER1(0,XPAR,NPAR)
c         CALL BOKER2(0,XPAR,NPAR)
C**         CALL BOKER3(0,XPAR,NPAR)
      ELSEIF(KEYOPT.EQ.3) THEN
C**      if(mod(iev,500).eq.1) write(6,*)    'event no=',iev
C**         CALL Z0DEC1(0,XPAR,NPAR)
C**         CALL BOK1PH(0,XPAR,NPAR)      
C**         CALL BOKER6(0,XPAR,NPAR)
      ELSEIF(KEYOPT.EQ.4) THEN
      if(mod(iev, 100 ).eq.1) write(6,*)    'event no=',iev
        if (KEYDIM.eq.5) then
         CALL Z0DEC3(0,XPAR,NPAR)
         CALL BOK3PH(0,XPAR,NPAR)
        else
         CALL Z0DEC2(0,XPAR,NPAR)
         CALL BOK2PH(0,XPAR,NPAR)
        endif

!        CALL BOKER7(0,XPAR,NPAR)
ccc        CALL TRZYNUE(0,XPAR,NPAR)
      ELSEIF(KEYOPT.EQ.5) THEN
C**      if(mod(iev, 100 ).eq.1) write(6,*)    'event no=',iev
C**        CALL Z0DEC0( 0,XPAR,NPAR)
ccc        CALL BOKER2( 0,XPAR,NPAR)
C**        CALL Z0RAD0( 0,XPAR,NPAR)
      ELSEIF(KEYOPT.EQ.6) THEN
C**      if(mod(iev, 100 ).eq.1) write(6,*)    'event no=',iev
C**        CALL Z0DEC1( 0,XPAR,NPAR)
ccc        CALL BOKER2( 0,XPAR,NPAR)
C**        CALL Z0RAD1( 0,XPAR,NPAR)
      ENDIF
   10 CONTINUE
C
C....post generation mode
      IF(KEYOPT.EQ.1) THEN
c         CALL PYSTAT(1)
c         CALL ROBOL0(1)
      ELSEIF(KEYOPT.EQ.2) THEN  
C**         CALL Z0DEC0( 1,XPAR,NPAR)
C**         CALL BOKER1( 1,XPAR,NPAR)
c         CALL BOKER2( 1,XPAR,NPAR)
C**         CALL BOKER3( 1,XPAR,NPAR)
      ELSEIF(KEYOPT.EQ.3) THEN  
C**         CALL Z0DEC1( 1,XPAR,NPAR)
C**         CALL BOK1PH( 1,XPAR,NPAR)
C**         CALL BOKER6( 1,XPAR,NPAR)
      ELSEIF(KEYOPT.EQ.4) THEN  
        if (KEYDIM.eq.5) then
         CALL Z0DEC3(1,XPAR,NPAR)
         CALL BOK3PH(1,XPAR,NPAR)
        else
         CALL Z0DEC2(1,XPAR,NPAR)
         CALL BOK2PH(1,XPAR,NPAR)
        endif
!         CALL BOKER7( 1,XPAR,NPAR)
ccc         CALL TRZYNUE( 1,XPAR,NPAR)
      ELSEIF(KEYOPT.EQ.5) THEN
C**        CALL Z0DEC0( 1,XPAR,NPAR)
ccc        CALL BOKER2( 1,XPAR,NPAR)
C**        CALL Z0RAD0( 1,XPAR,NPAR)
      ELSEIF(KEYOPT.EQ.6) THEN
C**        CALL Z0DEC1( 1,XPAR,NPAR)
ccc        CALL BOKER2( 1,XPAR,NPAR)
C**       CALL Z0RAD1( 1,XPAR,NPAR)
      ENDIF
C ------------WRITING HISTOS ON THE DISK ------------------------
      CALL GRFILE(NOUTH,' ','N')
      CALL GROUT( 0,ICY,' ')
      CALL GREND(DNAME)
C ------------THE END OF HISTO WRITING -------------------------
 
CC>>>>>>>>>>>>>>


      END



