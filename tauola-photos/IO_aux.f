
      SUBROUTINE fort_open(nout,fname)
*/////////////////////////////////////////////////////////////////////////////////////
*//   Interface used by c++ programs                                                //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      CHARACTER fname*(*)
      INTEGER nout,nout2
      nout2 = nout
      WRITE( *,*) 'KK2f_fort_open: nout = ',nout,'   fname= ',fname
      OPEN(nout2,file=fname)
      END

      SUBROUTINE fort_close(nout)
*/////////////////////////////////////////////////////////////////////////////////////
*//   Interface used by c++ programs                                                //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INTEGER nout,nout2
      WRITE( *,*) 'KK2f_fort_close: nout = ',nout
      CLOSE(nout2)
      END

      SUBROUTINE ReaDataZ(iReset,imax,xpar)
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
cc      CHARACTER*(*)     DiskFile
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

      IF( iReset .EQ. 1) THEN
         OPEN(ninp,file='./KKMChh_defaults')
      ELSE
         OPEN(ninp,file='./pro.input')
      ENDIF

      IF(iTalk .EQ. 1) THEN
         WRITE(  *,*) '****************************'
         WRITE(  *,*) '*  DZface_ReaDataX Starts  *'
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
         WRITE(  *,*) '*   DZface_ReaDataX Ends   *'
         WRITE(  *,*) '**************************'
      ENDIF
      CLOSE(ninp)
      RETURN
*-----------
 990  WRITE(    *,*) '+++ DZface_ReaDataX: wrong index= ',index
      STOP
      RETURN
 991  WRITE(    *,*) '+++ DZface_ReaDataX: wrong end of data '
      STOP
      END
