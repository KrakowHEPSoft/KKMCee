
      SUBROUTINE ReaDataY(DiskFile,iReset,iTalk,imax,xpar)
*/////////////////////////////////////////////////////////////////////////////////////
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
ccc      iTalk = 1
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
         WRITE(  *,*) '*    KK2f_ReaDataX Starts  *'
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
         WRITE(  *,*) '*   KK2f_ReaDataX Ends   *'
         WRITE(  *,*) '**************************'
      ENDIF
      CLOSE(ninp)
      RETURN
*-----------
 990  WRITE(    *,*) '+++ ReaDataY: wrong index= ',index
      STOP
      RETURN
 991  WRITE(    *,*) '+++ ReaDataY: wrong end of data '
      STOP
      END


      SUBROUTINE fort_open(nout,fname)
!     **********************************
* -------------------------------------------------------
*  Interface used by c++ programs
* -------------------------------------------------------
      CHARACTER fname*(*)

      nout2 = nout
      OPEN(nout2,file=fname)
**      WRITE(6,'(A,A20,A)')    '======>',filename,'<========='
**      WRITE(nout,'(A,A20,A)') '======>',filename,'<========='
      END

      SUBROUTINE fort_close(nout)
!     **************************
* -------------------------------------------------------
*  Interface used by c++ programs
* -------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 

      nout2 = nout
      CLOSE(nout2)

      END

      SUBROUTINE ReaDataX(ninp,xpar,imax)
*     ***********************************
*
*  Single data card is:    (a1,i4,d15.0,a60)
*  First character * defines comment card!
*
*     ***********************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      DIMENSION xpar(*)
      character*6 beg6
      character*4 end4
      character*1 mark1
      character*60 comm60
      character*80 comm80

*      WRITE(  *,*) '***************************'
*      WRITE(  *,*) '*  Parser ReaDataX starts *'
*      WRITE(  *,*) '***************************'
* Clear content
      DO i=1,imax
         xpar(i)=0d0
      ENDDO

* Search for 'BeginX'
      DO line =1,10000
         READ(ninp,'(a6,a)') beg6,comm60
         IF(beg6 .EQ. 'BeginX') THEN
            WRITE(  *,'(a6,a)') beg6,comm60
            GOTO 200
         ENDIF
      ENDDO
 200  CONTINUE

* Read data, 'EndX' terminates data, '*' marks comment
      DO line =1,1000
         READ(ninp,'(a)') mark1
         IF(mark1 .EQ. ' ') THEN
            BACKSPACE(ninp)
            READ(ninp,'(a1,i4,d15.0,a60)') mark1,index,value,comm60
            WRITE(  *,'(a1,i4,g15.6,a60)') mark1,index,value,comm60
            IF( (index .LE. 0) .OR. (index .GE. imax)) GOTO 990
            xpar(index) = value
         ELSEIF(mark1 .EQ. 'E') THEN
            BACKSPACE(ninp)
            READ(ninp,'(a4,a)') end4,comm60
            WRITE(  *,'(a4,a)') end4,comm60
            IF(end4 .EQ. 'EndX') GOTO 300
            GOTO 991
         ELSEIF(mark1 .EQ. '*') THEN
            BACKSPACE(ninp)
            READ(ninp,'(a)') comm80
            WRITE(  *,'(a)') comm80
         ENDIF
      ENDDO
 300  CONTINUE

*      WRITE(  *,*) '***************************'
*      WRITE(  *,*) '* Parser ReaDataX ends    *'
*      WRITE(  *,*) '***************************'
      RETURN
 990  WRITE(*,*) '+++ ReaDataX: wrong index= ',index
      STOP
      RETURN
 991  WRITE(*,*) '+++ ReaDataX: wrong end of data '
      STOP
      END

      SUBROUTINE ReaDataN(ninp,npar,imax)
*     ***********************************
*
*  Single data card is:    (a1,i4,i15,a60)
*  First character * defines comment card!
*
*     ***********************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      DIMENSION npar(*)
      character*6 beg6
      character*4 end4
      character*1 mark1
      character*60 comm60
      character*80 comm80

*      WRITE(  *,*) '***************************'
*      WRITE(  *,*) '* Parser  ReaDataN starts *'
*      WRITE(  *,*) '***************************'
* Clear content
      DO i=1,imax
         npar(i)=0
      ENDDO

* Search for 'BeginN'
      DO line =1,10000
         READ(ninp,'(a6,a)') beg6,comm60
         IF(beg6 .EQ. 'BeginN') THEN
            WRITE(  *,'(a6,a)') beg6,comm60
            GOTO 200
         ENDIF
      ENDDO
 200  CONTINUE

* Read data, 'EndN' terminates data, '*' marks comment
      DO line =1,1000
         READ(ninp,'(a)') mark1
         IF(mark1 .EQ. ' ') THEN
            BACKSPACE(ninp)
            READ(ninp,'(a1,i4,i15,a60)') mark1,index,nvalue,comm60
            WRITE(  *,'(a1,i4,i15,a60)') mark1,index,nvalue,comm60
            IF( (index .LE. 0) .OR. (index .GE. imax)) GOTO 990
            npar(index) = nvalue
         ELSEIF(mark1 .EQ. 'E') THEN
            BACKSPACE(ninp)
            READ(ninp,'(a4,a)') end4,comm60
            WRITE(  *,'(a4,a)') end4,comm60
            IF(end4 .EQ. 'EndN') GOTO 300
            GOTO 991
         ELSEIF(mark1 .EQ. '*') THEN
            BACKSPACE(ninp)
            READ(ninp,'(a)') comm80
            WRITE(  *,'(a)') comm80
         ENDIF
      ENDDO
 300  CONTINUE

*      WRITE(  *,*) '***************************'
*      WRITE(  *,*) '* Parser ReaDataN ends    *'
*      WRITE(  *,*) '***************************'
      RETURN
 990  WRITE(*,*) '+++ ReaDataN: wrong index= ',index
      STOP
      RETURN
 991  WRITE(*,*) '+++ ReaDataN: wrong end of data '
      STOP
      END

