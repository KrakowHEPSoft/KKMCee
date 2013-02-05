*//////////////////////////////////////////////////////////////////////////////
*// The package for plotting histograms in the LOGarithmic vertical scale    //
*// for the histogramming/plotting library GLK version 1.30 of S. Jadach.    //
*// Based on the GLK_Plot package from this library.                         //
*//////////////////////////////////////////////////////////////////////////////
*// Written by: Wieslaw Placzek                          CERN, November 1999 //
*// Last correction: 26.11.1999           by: WP                             //
*//////////////////////////////////////////////////////////////////////////////
*//////////////////////////////////////////////////////////////////////////////
*//     List of procedures,  non-user subprograms in brackets                //
*//////////////////////////////////////////////////////////////////////////////
*      SUBR/FUNC     1 PAR.  2 PAR.  3 PAR. 4 PAR. 5 PAR. 6 PAR.
*  ====================================================================
*      GLK_PlotLs      INT    CHR*1   CHR*1   INT   ----   ----
*     (GLK_PlFram1Ls)  INT      INT     DBL   DBL    DBL    DBL
*     (GLK_SAxisYLs)   INT      INT     DBL   DBL   ----   ----
*     (GLK_PlHistLs)   INT      INT     INT   DBL    DBL    DBL  INT  DBL  DBL
*     (GLK_PlHis2Ls)   INT      INT     INT   DBL    DBL    DBL  INT  DBL  DBL
*      GLK_Plot2Ls     INT    CHR*1   CHR*1 CHR*(*) CHR*16 CHR*16
*     (GLK_PlFrameLs)  INT      INT     DBL   DBL    DBL    DBL CHR*1
*     (GLK_PlMarkLs)   INT      INT     INT   DBL    DBL    DBL  ...
*//////////////////////////////////////////////////////////////////////////////

      SUBROUTINE GLK_PlotLs(ID,CH1,CH2,KDUM)
*     **************************************
      IMPLICIT NONE
      INCLUDE 'GLK.h'
      CHARACTER CH1,CH2,CHR
      CHARACTER*80 TITLE
      INTEGER          ID,KDUM
      DOUBLE PRECISION YY(m_MaxNb),YER(m_MaxNb)
      DOUBLE PRECISION YEL(m_MaxNb),YEU(m_MaxNb)
      LOGICAL GLK_EXIST
      INTEGER          IDUM,kax,kay,ioplog,iopsla,ioperb,iopsc1,iopsc2
      INTEGER          ker,NCHX,i
      DOUBLE PRECISION XL,XU,DXL,DXU,YL,YU,DltY,Yval,Yerr,Yerl
*-------------------------------------------------------------
      DATA CHR /' '/
* RETURN if histo non-existing
      IF(.NOT.GLK_EXIST(ID)) GOTO 900
* Unpack histogram
      CALL GLK_UNPAK(ID,YY ,'    ',IDUM)
      CALL GLK_UNPAK(ID,YER,'ERRO',IDUM)
      CALL GLK_HINBO1(ID,TITLE,NCHX,DXL,DXU)
      XL = DXL
      XU = DXU
      CALL GLK_RANGE1(ID,YL,YU)
* Min and max values of y log scale
      IF (YU.GT.0d0) THEN
         YU = LOG10(YU)
      ELSE
         WRITE(6,*)
     &        '+++ GLK_PlotLs: LOG scale impossible for Ymax < 0: Ymax = ',YU
         STOP
      ENDIF
      IF (YL.GT.0d0) THEN
         YL = LOG10(YL)
      ELSE
         YL = YU - 10d0
      ENDIF
      DltY = YU - YL
      IF (DltY.LE.0d0) THEN
         WRITE(6,*)' Wrong values of: Ymin, Ymax = ',YL,YU
         STOP
      ELSEIF (DltY.GT.20d0) THEN
         YL = YU - 20d0
      ENDIF
* Transform y values and their errors into log scale 
      DO i = 1,NCHX
         Yval = YY(i)
         IF (Yval.GT.0d0) THEN
            YY(i) = MAX( LOG10(Yval), YL )
            Yerr = YER(i)
            Yerl = Yval - Yerr 
            IF (Yerl.GT.0d0) THEN
               YEL(i) = MAX( LOG10(Yerl), YL )
            ELSE
               YEL(i) = YL
            ENDIF
            Yerl = Yval + Yerr
            YEU(i) = MIN( MAX(LOG10(Yerl),YL), YU )
         ELSE
            YY(i)  = YL 
            YEL(i) = YL
            YEU(i) = YL
         ENDIF
      ENDDO
      kax = 1200
      kay = 1200
      IF (ch1 .EQ. 'S') THEN
* Superimpose plot
        BACKSPACE(m_ltx)
        BACKSPACE(m_ltx)
      ELSE
* New frame only
        CHR=CH1
        CALL GLK_PlFram1Ls(kax,kay,XL,XU,YL,YU)
      ENDIF
      WRITE(m_ltx,'(A)')     '%========== next plot (line) =========='
      WRITE(m_ltx,'(A,I10)') '%==== HISTOGRAM ID=',ID
      WRITE(m_ltx,'(A,A70 )')'% ',TITLE
* Continuous line for functions
      CALL GLK_OptOut(id,ioplog,iopsla,ioperb,iopsc1,iopsc2)
      ker = ioperb - 1
      IF (iopsla .EQ. 2)  CHR='C'
* Suppress GLK_PlotLs assignments
      IF (CH2 .EQ. 'B')   CHR=' '
      IF (CH2 .EQ. '*')   CHR='*'
      IF (CH2 .EQ. 'C')   CHR='C'
* Various types of lines
      IF     (CHR .EQ. ' ') THEN
* Contour line used for histogram
          CALL GLK_PlHistLs(kax,kay,NCHX,YL,YU,YY,ker,YEL,YEU)
      ELSE IF(CHR .EQ. '*') THEN
* Marks in the midle of the bin
          CALL GLK_PlHis2Ls(kax,kay,NCHX,YL,YU,YY,ker,YEL,YEU)
      ELSE IF(CHR .EQ. 'C') THEN
* Slanted (dotted) line in plotting non-MC functions
          CALL GLK_PlCirc(kax,kay,NCHX,YL,YU,YY)
      ENDIF
!------------------------------------------------------------------
*------------------------------!
* Ending
*------------------------------!
      WRITE(m_ltx,'(2A)') m_BS,'end{picture} % close entire picture '
      WRITE(m_ltx,'(2A)') m_BS,'end{figure}'

      RETURN
  900 CALL GLK_Retu1('+++GLK_PlotLs: Nonexistig histo, skipped, id=' ,ID)
      END

      SUBROUTINE GLK_PlFram1Ls(kax,kay,XL,XU,YL,YU)
*     *********************************************
      IMPLICIT NONE
      INCLUDE 'GLK.h'
      INTEGER           kax,kay
      CHARACTER*80 TITLE
      DOUBLE PRECISION   TIPSX(20)
      DOUBLE PRECISION   XL,XU
      INTEGER            ntipy,ntipx,nchx,icont
      DOUBLE PRECISION   YU,YL
      DATA ICONT/0/
*----------------
      ICONT=ICONT+1
      IF(ICONT .GT. 1) WRITE(m_ltx,'(2A)') m_BS,'newpage'
*------------------------------!
*           Header
*------------------------------!
      WRITE(m_ltx,'(A)') ' '
      WRITE(m_ltx,'(A)') ' '
      WRITE(m_ltx,'(A)') '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     $%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      WRITE(m_ltx,'(A)') '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     $%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      WRITE(m_ltx,'(2A)') m_BS,'begin{figure}[!ht]'
      WRITE(m_ltx,'(2A)') m_BS,'centering'
*------------------------------!
* General Caption
*------------------------------!
      WRITE(m_ltx,'(4A)') m_BS,'caption{',m_BS,'small'
      IF(M_KEYTIT.EQ.0) THEN
        WRITE(m_ltx,'(A)')     TITLE
      ELSE
        WRITE(m_ltx,'(A)')     m_titch(1)
      ENDIF
      WRITE(m_ltx,'(A)') '}'
*------------------------------!
* Frames and labels
*------------------------------!
      WRITE(m_ltx,'(A)') '% =========== big frame, title etc. ======='
      WRITE(m_ltx,'(4A)') m_BS,'setlength{',m_BS,'unitlength}{0.1mm}'
      WRITE(m_ltx,'(2A)') m_BS,'begin{picture}(1600,1500)'
      WRITE(m_ltx,'(4A)')
     $     m_BS,'put(0,0){',m_BS,'framebox(1600,1500){ }}'
      WRITE(m_ltx,'(A)') '% =========== small frame, labeled axis ==='
      WRITE(m_ltx,'(4A,I4,A,I4,A)')
     $    m_BS,'put(300,250){',m_BS,'begin{picture}( ',kax,',',kay,')'
      WRITE(m_ltx,'(4A,I4,A,I4,A)')
     $    m_BS,'put(0,0){',m_BS,'framebox( ',kax,',',kay,'){ }}'
      WRITE(m_ltx,'(A)') '% =========== x and y axis ================'
      CALL GLK_SAxisX(kax,XL,XU,NTIPX,TIPSX)
      CALL GLK_SAxisYLs(kax,kay,YL,YU)
      WRITE(m_ltx,'(3A)') m_BS,'end{picture}}'
     $                ,'% end of plotting labeled axis'
      END

      SUBROUTINE GLK_SAxisYLs(kax,kay,YL,YU)
*     **************************************
* Plotting y-axis with long and short tips in LOG scale
      IMPLICIT NONE
      INCLUDE 'GLK.h'
      INTEGER          kax,kay
      DOUBLE PRECISION YL,YU
*
      DOUBLE PRECISION pll,plu,pdl,psl,pds,psb,pst
      INTEGER          nll,nlu,nlt,nst,its,n,k
*---------------------------------------------------
* Find number of long ticks
      nll = NINT(yl + 0.49999999d0)
      nlu = NINT(yu - 0.49999999d0)
      nlt = nlu - nll + 1
      WRITE(m_ltx,'(A)') '% .......GLK_SAxisYLs........ '
      WRITE(m_ltx,'(A,I4)') '%  NLT= ',nlt
* Positions of first and last long ticks
      pll = kay*(nll - YL)/(YU - YL)
      plu = kay*(nlu - YL)/(YU - YL)
* Distance between two long ticks
      pdl = (plu - pll)/MAX(nlt-1,1)
* Plotting long tics on vertical axis
      WRITE(m_ltx,1000)
     $     m_BS,'multiput(0,'  ,pll,  ')(0,'  ,pdl  ,'){'  ,nlt,  '}{',
     $     m_BS,'line(1,0){25}}'
      WRITE(m_ltx,1001)
     $     m_BS,'multiput('  ,kax,  ','  ,pll,  ')(0,'  ,pdl,  '){'  ,nlt,
     $     '}{',m_BS,'line(-1,0){25}}'
* Plotting short tics on vertical axis (a bit more complicated)
      DO its = 2,9
         nst = nlt - 1
         psl = pll + pdl*LOG10(DBLE(its))
         pds = pdl
         psb = psl - pds
         IF (psb.GE.0) THEN
            psl = psb
            nst = nst + 1
         ENDIF
         pst = psl + nst*pds
         IF (pst.LE.kay) nst = nst + 1
         WRITE(m_ltx,1000)
     $        m_BS,'multiput(0,'  ,psl,  ')(0,'  ,pds,  '){'  ,nst,  '}{',
     $        m_BS,'line(1,0){10}}'
         WRITE(m_ltx,1001)
     $        m_BS,'multiput('  ,kax,  ','  ,psl,  ')(0,'  ,pds,  '){'  ,
     $        nst,'}{',m_BS,'line(-1,0){10}}'
      ENDDO
 1000 FORMAT(2A,F8.2,A,F8.2,A,I4,3A)
 1001 FORMAT(2A,I4,A,F8.2,A,F8.2,A,I4,3A)
* Labeling of axis
      DO n = nll,nlu
         k = NINT(pll + (n-nll)*pdl)
         WRITE(m_ltx,'(2A,I4,5A,I4,A)')
     $        m_BS,'put(-25,'  ,k,  '){',m_BS,'makebox(0,0)[r]{',
     $        m_BS,'large $ 10^{'  ,n,  '} $}}'
      ENDDO
      END

      SUBROUTINE GLK_PlHistLs(kax,kay,NCHX,YL,YU,YY,ker,YEL,YEU)
*     **********************************************************
* Plotting contour line for histogram in LOG scale
*     ***********************
      IMPLICIT NONE
      INCLUDE 'GLK.h'
      INTEGER           kax,kay,NCHX,ker
      DOUBLE PRECISION  YL,YU,YY(*),YEL(*),YEU(*)
      CHARACTER*80 FMT1
*
      INTEGER           ix0,ix1,ix2,idx,ie,ierr,idy,ib,iy0,iy1
*-------------------------------------------------------------
      WRITE(m_ltx,'(4A,I4,A,I4,A)')
     $     m_BS,'put(300,250){',m_BS,'begin{picture}( ',kax,',',kay,')'
      WRITE(m_ltx,'(A)') '% ========== plotting primitives =========='
* Various types of line
      IF(m_tline .EQ. 1) THEN
         WRITE(m_ltx,'(2A)') m_BS,'thicklines '
      ELSE
         WRITE(m_ltx,'(2A)') m_BS,'thinlines '
      ENDIF
* Short macros for vertical/horizontal straight lines
      WRITE(m_ltx,'(8A)')
     $     m_BS,'newcommand{',m_BS,'x}[3]{',m_BS,'put(#1,#2){',
     $     m_BS,'line(1,0){#3}}}'
      WRITE(m_ltx,'(8A)')
     $     m_BS,'newcommand{',m_BS,'y}[3]{',m_BS,'put(#1,#2){',
     $     m_BS,'line(0,1){#3}}}'
      WRITE(m_ltx,'(8A)')
     $     m_BS,'newcommand{',m_BS,'z}[3]{',m_BS,'put(#1,#2){',
     $     m_BS,'line(0,-1){#3}}}'
* Error bars (slightly different for LOG scale)
      WRITE(m_ltx,'(8A)')
     $   m_BS,'newcommand{',m_BS,'e}[3]{',
     $   m_BS,'put(#1,#2){',m_BS,'line(0,1){#3}}}'
      ix0 = 0
      iy0 = 0
      DO ib = 1,NCHX
         ix1 = NINT( kax*(ib-0.00001d0)/NCHX )
         iy1 = NINT( kay*(YY(ib)-YL)/(YU-YL) ) 
         idy = iy1-iy0
         idx = ix1-ix0
         FMT1 = '(2(2A,I4,A,I4,A,I4,A))'
         IF (idy .GE. 0) THEN  
            IF (iy1.GE.0 .AND. iy1.LE.kay)
     $           WRITE(m_ltx,FMT1) m_BS,'y{',ix0,'}{',iy0,'}{',idy,'}',
     $           m_BS,'x{',ix0,'}{',iy1,'}{',idx,'}'
         ELSE
            IF (iy1.GE.0 .AND. iy1.LE.kay)
     $           WRITE(m_ltx,FMT1) m_BS,'z{',ix0,'}{',iy0,'}{',-idy,'}',
     $           m_BS,'x{',ix0,'}{',iy1,'}{',idx,'}'
         ENDIF
         ix0 = ix1
         iy0 = iy1
         IF (ker.EQ.1) THEN
            ix2  = NINT(kax*(ib-0.5d0)/NCHX)
            ierr = NINT( kay*(YEL(ib) - YL)/(YU-YL) )
            ie   = NINT( kay*(YEU(ib) - YEL(ib))/(YU-YL) )
            IF(iy1.GE.0 .AND. iy1.LE.kay .AND. ABS(ierr).LE.9999
     $           .AND. ie.LE.9999) WRITE(m_ltx,8000) m_BS,ix2,ierr,ie  
         ENDIF
      ENDDO   
 8000 FORMAT(4(A1,2He{,I4,2H}{,I5,2H}{,I4,1H}:1X ))
      WRITE(m_ltx,'(3A)') m_BS,'end{picture}}',
     $       ' % end of plotting histogram'
* Change line-style
      m_tline= m_tline+1
      IF(m_tline .GT. 2) m_tline=1
      END

      SUBROUTINE GLK_PlHis2Ls(kax,kay,NCHX,YL,YU,YY,ker,YEL,YEU)
*     **********************************************************
* Marks in the midle of the bin in LOG scale
*     **********************************
      IMPLICIT NONE
      INCLUDE 'GLK.h'
      DOUBLE PRECISION  YL,YU,YY(*),YEL(*),YEU(*)
      INTEGER           kax,kay,NCHX,ker
*
      INTEGER           iy1,ierr,ie,ix1,irad1,irad2,ib
*-----------------------------------------------------
      WRITE(m_ltx,'(4A,I4,A,I4,A)')
     $ m_BS,'put(300,250){',m_BS,'begin{picture}( ',kax,',',kay,')'
      WRITE(m_ltx,'(A)') '% ========== plotting primitives =========='
* Various types of mark
      irad1 = 6
      irad2 =10
      IF(m_tline .EQ. 1) THEN
*- Small filled circle
         WRITE(m_ltx,'(8A,I3,A)')
     $        m_BS,'newcommand{',m_BS,'R}[2]{',
     $        m_BS,'put(#1,#2){',m_BS,'circle*{',IRAD1,'}}}'
      ELSEIF(m_tline .EQ. 2) THEN
*- Small open circle
         WRITE(m_ltx,'(8A,I3,A)')
     $        m_BS,'newcommand{',m_BS,'R}[2]{',
     $        m_BS,'put(#1,#2){',m_BS,'circle{',IRAD1,'}}}'
      ELSEIF(m_tline .EQ. 3) THEN
*- Big filled circle
         WRITE(m_ltx,'(8A,I3,A)')
     $        m_BS,'newcommand{',m_BS,'R}[2]{',
     $        m_BS,'put(#1,#2){',m_BS,'circle*{',IRAD2,'}}}'
      ELSEIF(m_tline .EQ. 4) THEN
*- Big open circle
         WRITE(m_ltx,'(8A,I3,A)')
     $        m_BS,'newcommand{',m_BS,'R}[2]{',
     $        m_BS,'put(#1,#2){',m_BS,'circle{',IRAD2,'}}}'
*- Other symbols
      ELSEIF(m_tline .EQ. 5) THEN
         WRITE(m_ltx,'(10A)')
     $        m_BS,'newcommand{',m_BS,'R}[2]{',
     $        m_BS,'put(#1,#2){',m_BS,'makebox(0,0){$',m_BS,'diamond$}}}'
      ELSE
         WRITE(m_ltx,'(10A)')
     $        m_BS,'newcommand{',m_BS,'R}[2]{',
     $        m_BS,'put(#1,#2){',m_BS,'makebox(0,0){$',m_BS,'star$}}}'
      ENDIF
*- Error bars
      WRITE(m_ltx,'(8A)')
     $     m_BS,'newcommand{',m_BS,'E}[3]{',
     $     m_BS,'put(#1,#2){',m_BS,'line(0,1){#3}}}'
      DO ib = 1,NCHX
         ix1 = NINT(kax*(ib-0.5d0)/NCHX)
         iy1 = NINT(kay*(YY(ib)-YL)/(YU-YL))
         IF(iy1.GT.0 .AND. iy1.LE.kay) WRITE(m_ltx,7000) m_BS,ix1,iy1
         IF(ker.EQ.1) THEN
            ierr = NINT( kay*(YEL(ib) - YL)/(YU-YL) )
            ie   = NINT( kay*(YEU(ib) - YEL(ib))/(YU-YL) )
            IF(iy1.GT.0 .AND. iy1.LE.kay .AND. ABS(ierr).LE.9999
     $           .AND. ie.LE.9999) WRITE(m_ltx,8000) m_BS,ix1,ierr,ie   
         ENDIF
      ENDDO
 7000 FORMAT(4(A1,2HR{,I4,2H}{,I4,1H}:1X ))
 8000 FORMAT(4(A1,2HE{,I4,2H}{,I5,2H}{,I4,1H}:1X ))
      WRITE(m_ltx,'(3A)') m_BS,'end{picture}}',
     $    ' % end of plotting histogram'
* Change line-style
      m_tline= m_tline+1
      IF(m_tline .GT. 6) m_tline=1
      END

      SUBROUTINE GLK_Plot2Ls(id,CH1,CH2,chmark,chxfmt,chdumm)
*     *******************************************************
* New version, more user-friendly, of GLK_PlotLs
* INPUT:
*    ID          histogram identifier
*    ch1 = ' '   normal new plot
*        = 'S'   impose new plot on previous one
*    ch2 = ' '   ploting line default, contour
*        = '*'   error bars in midle of the bin
*        = 'R'   error bars at Right edge of the bin
*        = 'L'   error bars at Left  edge of the bin
*        = 'C'   slanted continuous smooth line
*    chmark =    TeX symbol for ploting points
*    chxfmt =    format (string) for labeling x-axis
*    chdumm =    dummy character (string) parameter 
*                (kept only for compatibility with GLK_Plot2) 
*----------------------------------------------------------------------
* NOTE: For LOG scale the format for y-axis label is generated      
*       automatically and does not have to be specified by the user. 
*       This is why the parameter chyfmt of GLK_Plot2 is not needed here. 
*----------------------------------------------------------------------
* Furthermore:
* Captions are defined by means of 
*    CALL gplcapt(capt) before CALL gplot2
*    where CHARACTER*64 capt(50) is content of 
*    caption, line by line, see also comments in gplcapt routine.
* Additional text as a TeX source text can be appended by means of
*    CALL gplabel(lines) after CALL gplot2
*    where CHARACTER*64 lines(50) is the TeX add-on.
*    this is used to decorate plot with
*    any kind marks, special labels and text on the plot.
***********************************************************************
      IMPLICIT NONE
      INTEGER id
      CHARACTER CH1,CH2,chmark*(*)
      CHARACTER*(*) chxfmt,chdumm
      INCLUDE 'GLK.h'
      SAVE
      DOUBLE PRECISION  YY(m_MaxNb),YER(m_MaxNb),YEL(m_MaxNb),YEU(m_MaxNb)
      CHARACTER*80 TITLE
*---------------------------------------------------------------------
      LOGICAL GLK_Exist
      INTEGER kax,kay,incr,ker,NCHX
      INTEGER iopsla,ioplog,ioperb,iopsc1,iopsc2,idum,i
      DOUBLE PRECISION   dxl,dxu,XU,XL,YU,YL,Yval,DltY,Yerr,Yerl
      CHARACTER CHR
      DATA CHR /' '/
* TeX Names of the error-bar command and of the point-mark command
      CHARACTER*1 chre, chrp1
      PARAMETER ( chre = 'E', chrp1= 'R' )
      CHARACTER*2 chrp
* TeX Name of the point-mark command
      CHARACTER*1 chrx(12)
      DATA  chrx /'a','b','c','d','f','g','h','i','j','k','l','m'/
*---------------------------------------------------------------------
* Return if histo non-existing
      IF(.NOT.GLK_Exist(id)) GOTO 900
* Unpack histogram
      CALL GLK_UnPak(id,yy ,'    ',idum)
      CALL GLK_UnPak(id,yer,'ERRO',idum)
      CALL GLK_hinbo1(id,TITLE,NCHX,dxl,dxu)
* The X-range as in histo
      XL = dxl
      XU = dxu
* Header
      kax = 1200
      kay = 1200
      IF(CH1 .EQ. 'S') THEN
* Superimpose plot
        incr = incr + 1
        BACKSPACE(m_ltx)
        BACKSPACE(m_ltx)
      ELSE
* New frame only
        incr = 1
        CHR = CH1
* The Y-range from first plot is preserved
        CALL GLK_Range1(id,yl,yu)
* Min and max values of y log scale
        IF (YU.GT.0d0) THEN
           YU = LOG10(YU)
        ELSE
           WRITE(6,*)'+++ GLK_Plot2Ls: LOG scale impossible for Ymax < 0: Ymax = ',Yu
           STOP
        ENDIF
        IF (YL.GT.0d0) THEN
           YL = LOG10(YL)
        ELSE
           YL = YU - 10d0
        ENDIF
        DltY = YU - YL
        IF (DltY.LE.0d0) THEN
           WRITE(6,*)'+++ GLK_Plot2Ls: Wrong values of: Ymin, Ymax = ',YL,YU
           STOP
        ELSEIF (DltY.GT.20d0) THEN
           YL = YU - 20d0
        ENDIF
* Transform y values and their errors into log scale 
        DO i = 1,NCHX
           Yval = YY(i)
           IF (Yval.GT.0d0) THEN
              YY(i) = MAX( LOG10(yval), YL )
              Yerr = Yer(i)
              Yerl = Yval - Yerr 
              IF (Yerl.GT.0d0) THEN
                 Yel(i) = MAX( LOG10(Yerl), YL )
              ELSE
                 Yel(i) = YL
              ENDIF
              Yerl = Yval + Yerr
              Yeu(i) = MIN( MAX(LOG10(Yerl),YL), YU )
           ELSE
              YY(i)  = YL 
              Yel(i) = YL
              Yeu(i) = YL
           ENDIF 
        ENDDO
        CALL GLK_LFrameLs(kax,kay,XL,XU,YL,YU,chxfmt)
      ENDIF
      chrp= chrp1//chrx(incr)
      WRITE(m_ltx,'(A)')    '%=GLK_Plot2Ls:  next plot (line) =========='
      WRITE(m_ltx,'(A,I10)')'%====HISTOGRAM ID=',ID
      WRITE(m_ltx,'(A,A70 )') '% ',TITLE
      CALL GLK_OptOut(id,ioplog,iopsla,ioperb,iopsc1,iopsc2)
      ker = ioperb - 1
* Default line type
      IF (iopsla .EQ. 2) THEN 
         CHR = 'C'
      ELSE
         CHR = ' '
      ENDIF
* User defined line-type
      IF (CH2 .EQ. 'B')   CHR=' '
* Marks in the midle of the bin
      IF (CH2 .EQ. '*')   CHR='*'
* Marks on the right edge of the bin
      IF (CH2 .EQ. 'R')   CHR='R'
* Marks on the left edge of the bin
      IF (CH2 .EQ. 'L')   CHR='L'
      IF (CH2 .EQ. 'C')   CHR='C'
* Various types of lines
      IF (CHR.EQ.' ') THEN
* Contour line used for histogram
          CALL GLK_PlHistLs(kax,kay,NCHX,YL,YU,YY,ker,YEL,YEU)
      ELSE IF(CHR.EQ.'*' .OR. CHR.EQ.'R' .OR. CHR.EQ.'L') THEN
* Marks on the right/left/midle of the bin
         CALL GLK_PlMarkLs(kax,kay,NCHX,YL,YU,YY,ker,YEL,YEU,
     &        chmark,chr,chrp,chre)
      ELSE IF(chr .EQ. 'C') THEN
* Slanted (dotted) line in plotting non-MC functions
          CALL GLK_PlCirc(kax,kay,nchx,yl,yu,yy)
       ELSE
          WRITE(6,*)' ++++ GPlot2Ls: Wrong mark position label: ',ch2
       ENDIF
      WRITE(m_ltx,'(2A)') m_BS,'end{picture} % close entire picture '
      IF(ABS(m_lint) .EQ. 2) THEN
         WRITE(m_ltx,'(A)') '%== GLK_Plot2Ls:  end of plot  =========='
      ELSE
         WRITE(m_ltx,'(2A)') m_BS,'end{figure}'
      ENDIF
      RETURN
  900 CALL GLK_Stop1('+++GLK_Plot2Ls: Nonexistig histo, skipped, id= ',ID)
      END

      SUBROUTINE GLK_LFrameLs(kax,kay,XL,XU,YL,YU,chxfmt)
*     ***************************************************
      IMPLICIT NONE
      INTEGER kax,kay
      CHARACTER chxfmt*(*)
      INCLUDE 'GLK.h'
      SAVE
*---------------------------------------------------
      CHARACTER*80 TITLE
      DOUBLE PRECISION    XL,XU,YL,YU
      INTEGER  icont,i
      DATA icont/0/
*---------------------------------------------------
      icont = icont + 1
*
      IF(icont .GT. 1) WRITE(m_ltx,'(2A)') m_BS,'newpage'
*------------------------------!
*           Header
*------------------------------!
      WRITE(m_ltx,'(A)') ' '
      WRITE(m_ltx,'(A)') ' '
      WRITE(m_ltx,'(A)')
     $'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      WRITE(m_ltx,'(A)')
     $'%%%%%%%%%%%%%%%%%%%%%%GLK_PlFrameLs%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      IF(ABS(m_lint) .EQ. 2) THEN
         WRITE(m_ltx,'(2A)') m_BS,'noindent'
      ELSE
         WRITE(m_ltx,'(2A)') m_BS,'begin{figure}[!ht]'
         WRITE(m_ltx,'(2A)') m_BS,'centering'
         WRITE(m_ltx,'(2A)') m_BS,'htmlimage{scale=1.4}'
      ENDIF
*------------------------------!
* General Caption
*------------------------------!
      IF(ABS(m_lint) .NE. 2) THEN
         WRITE(m_ltx,'(6A)')
     $        m_BS,'caption{',m_BS,'footnotesize',m_BS,'sf'
         DO i=1,m_KeyTit
            WRITE(m_ltx,'(A)')     m_titch(i)
         ENDDO
         WRITE(m_ltx,'(A)') '}'
      ENDIF
*------------------------------!
* Frames and labels
*------------------------------!
      WRITE(m_ltx,'(A)') '% =========== big frame, title etc. ======='
      WRITE(m_ltx,'(4A)') m_BS,'setlength{',m_BS,'unitlength}{0.1mm}'
      WRITE(m_ltx,'(2A)') m_BS,'begin{picture}(1600,1500)'
      IF( m_lint .LT. 0) THEN
* Big frame usefull for debuging
         WRITE(m_ltx,'(4A)')
     $        m_BS,'put(0,0){',m_BS,'framebox(1600,1500){ }}'
      ENDIF
      WRITE(m_ltx,'(A)') '% =========== small frame, labeled axis ==='
      WRITE(m_ltx,'(4A,I4,A,I4,A)')
     $    m_BS,'put(300,250){',m_BS,'begin{picture}( ',kax,',',kay,')'
      WRITE(m_ltx,'(4A,I4,A,I4,A)')
     $    m_BS,'put(0,0){',m_BS,'framebox( ',kax,',',kay,'){ }}'
      WRITE(m_ltx,'(A)') '% =========== x and y axis ================'
      CALL GLK_AxisX(kax,XL,XU,chxfmt)
      CALL GLK_SAxisYLs(kax,kay,YL,YU)
      WRITE(m_ltx,'(3A)') m_BS,'end{picture}}'
     $                ,'% end of plotting labeled axis'
      END

      SUBROUTINE GLK_PlMarkLs(kax,kay,NCHX,YL,YU,YY,ker,YEL,YEU,
     &                        chmark,chr,chr2,chr3)
*     **********************************************************
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//                       marks in the midle of the bin                             //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INTEGER     kax,kay,NCHX,ker
      DOUBLE PRECISION       YL,YU, YY(*),YEL(*),YEU(*)
      CHARACTER*1 chr
      CHARACTER   chmark*(*),chr2*(*),chr3*(*)
*---------------------------------------------------
      INCLUDE 'GLK.h'
      SAVE
      INTEGER    ib,ix1,iy1,ierr,ie
*---------------------------------------------------
      WRITE(m_ltx,'(4A,I4,A,I4,A)') m_BS,'put(300,250){',m_BS,'begin{picture}( ',kax,',',kay,')'
      WRITE(m_ltx,'(A)') '% ===GLK_PlMarkLs: plotting primitives ======'
* Color string, optionaly
      IF(m_KeyCol .EQ. 1) THEN
         WRITE(m_ltx,'(A)') m_Color
         m_KeyCol = 0
      ENDIF
* Plotting symbol
      WRITE(m_ltx,'(10A)') m_BS,'newcommand{',m_BS,chr2  , '}[2]{', m_BS,'put(#1,#2){',chmark,'}}'
* Error bar symbol
      WRITE(m_ltx,'(10A)')
     $   m_BS,'newcommand{',m_BS,chr3  , '}[3]{', m_BS,'put(#1,#2){',m_BS,'line(0,1){#3}}}'
* Print marks and values
      DO ib = 1,NCHX
         IF (chr .EQ. '*') THEN
            ix1 = NINT(kax*(ib-0.5d0)/nchx)    ! Midle of bin
         ELSEIF (chr .EQ. 'R') THEN
            ix1 = NINT(kax*(ib*1d0)/nchx)      ! Right edge of bin
         ELSEIF (chr .EQ. 'L') THEN
            ix1 = NINT(kax*(ib-1d0)/nchx)      ! Left edge of bin
         ELSE
            WRITE(6,*) '+++++ GLK_PlMarkLs: Wrong line type: ',chr
            RETURN
         ENDIF
         iy1 = NINT(kay*(YY(ib)-YL)/(YU-YL))
         IF(iy1.GT.0 .AND. iy1.LE.kay)
     $   WRITE(m_ltx,'(A,A,A,I4,A,I4,A)') 
     $               m_BS,chr2, '{' ,ix1, '}{' ,iy1, '}'
         IF(ker.EQ.1) THEN
            ierr = NINT( kay*(YEL(ib) - YL)/(YU-YL) )       ! bottom of error bar
            ie   = NINT( kay*(YEU(ib) - YEL(ib))/(YU-YL) )  ! total length of error bar
            IF(iy1.GT.0 .AND. iy1.LE.KAY .AND. ABS(ierr).LE.9999
     $         .AND. ie.LE.9999) WRITE(m_ltx,'(A,A,A,I4,A,I5,A,I4,A)') 
     $           m_BS, chr3,  '{'  ,ix1, '}{'  ,ierr, '}{'  ,ie,   '}'
         ENDIF
      ENDDO
      WRITE(m_ltx,'(3A)') m_BS,'end{picture}}',
     $    ' % end of plotting histogram'
      END
