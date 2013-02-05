      SUBROUTINE ReaDat(Dname,imax,xpar)
*     **********************************
*///////////////////////////////////////////////////////////////////////
*//                                                                   //
*//   Reads default xpar and actual xpar from disk file Dname         //
*//                                                                   //
*///////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      CHARACTER*60 Dname
*---
      INTEGER imax
      REAL*8  xpar(*)
*---
      CHARACTER*80 DiskFile
      INTEGER i,ninpy,ninpd
*-----------------------------------------------------------------------------
* Read DEFAULT input data for MC generator into xpar (no other action!)
* Name with full path has to be provided by the user
      DiskFile = '../../.KK2f_defaults'
      ninpy = 3
      OPEN( ninpy, file=DiskFile)
*     ***************************
      DO i=1,imax
         xpar(i)=0d0
      ENDDO
      CALL ReaDataX(ninpy,xpar,imax)
*     ******************************
      CLOSE(ninpy)
*---
      WRITE(6,*) '|ReaDat||',Dname,'|||'
      ninpd =4
      OPEN( ninpd, file=Dname)
*---
* Read actual data for MC generator, the defaults may now be overwritten!
*     *****************************
      CALL ReaDataX(ninpd,xpar,imax)
*     *****************************
*
      CLOSE(ninpd)
*
      END     

      SUBROUTINE ReaDataX(ninp,xpar,imax)
*     ***********************************
*///////////////////////////////////////////////////////////////////////
*//                                                                   //
*//  Single data card is:    (a1,i4,d15.0,a60)                        //
*//  First character * defines comment card!                          //
*//                                                                   //
*//  Note that this program does not clear xpar!!!                    //
*//  one has to do it before calling it, if necessary!!!              //
*//                                                                   //
*///////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INTEGER ninp,imax
      REAL*8  xpar(*)
*---
      CHARACTER*6 beg6
      CHARACTER*4 end4
      CHARACTER*1 mark1
      CHARACTER*60 comm60
      CHARACTER*80 comm80
*---
      INTEGER line,index
      REAL*8  value
*-------------------------------------------------------------------
*      WRITE(  *,*) '***************************'
*      WRITE(  *,*) '*  Parser ReaDataX starts *'
*      WRITE(  *,*) '***************************'

* Search for 'BeginX'
      DO line =1,100000
         READ(ninp,'(a6,a)') beg6,comm60
         IF(beg6 .EQ. 'BeginX') THEN
            WRITE(  *,'(a6,a)') beg6,comm60
            GOTO 200
         ENDIF
      ENDDO
 200  CONTINUE

* Read data, 'EndX' terminates data, '*' marks comment
      DO line =1,10000
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

