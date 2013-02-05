      SUBROUTINE ReaDat(Dname)
*     **********************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*60 Dname
*==========================================
      INTEGER imax
      PARAMETER (imax=10000)
      COMMON  /xpar/  xpar(imax)
*==========================================
      SAVE
      COMMON / inout  / ninp,nout
      CHARACTER*80 DiskFile
*-----------------------------------------------------------------------------
* Read DEFAULT input data for MC generator into xpar (no other action!)
* Name with full path has to be provided by the user
      DiskFile = '../../.KK2f_defaults'
      ninpX = 3
      OPEN( ninpX, file=DiskFile)
*     ***************************
      DO i=1,imax
         xpar(i)=0d0
      ENDDO
      CALL ReaDataX(ninpX,nout,xpar,imax)
*     ******************************
      CLOSE(ninpX)
*=============================================================================
      ninpD =4
      WRITE(6,*) '|ReaDat||',Dname,'|||'
      OPEN( ninpD, file=Dname)
* Read actual data for MC generator, the defaults may now be overwritten!
*     *****************************
      CALL ReaDataX(ninpD,6,xpar,imax)
*     *****************************
      CLOSE(ninpD)
*=============================================================================
* Initialize BornV class using the same input as in MC run
* This will provide data base on particle properties
* and  Born x-sections.
      CALL BornV_Initialize(xpar)
      END     


      SUBROUTINE ReaDataX(ninp,nout,xpar,imax)
*     ****************************************
*
*  Single data card is:    (a1,i4,d15.0,a60)
*  First character * defines comment card!
*
*  Note that this program does not clear xpar!!!
*  one has to do it before calling it, if necessary!!!
*
*     ***********************************
      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      DIMENSION xpar(*)
      character*6 beg6
      character*4 end4
      character*1 mark1
      character*60 comm60
      character*80 comm80

*      WRITE(nout,*) '***************************'
*      WRITE(nout,*) '*  Parser ReaDataX starts *'
*      WRITE(nout,*) '***************************'

* Search for 'BeginX'
      DO line =1,10000
         READ(ninp,'(a6,a)') beg6,comm60
         IF(beg6 .EQ. 'BeginX') THEN
            WRITE(nout,'(a6,a)') beg6,comm60
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
            WRITE(nout,'(a1,i4,g15.6,a60)') mark1,index,value,comm60
            IF( (index .LE. 0) .OR. (index .GE. imax)) GOTO 990
            xpar(index) = value
         ELSEIF(mark1 .EQ. 'E') THEN
            BACKSPACE(ninp)
            READ(ninp,'(a4,a)') end4,comm60
            WRITE(nout,'(a4,a)') end4,comm60
            IF(end4 .EQ. 'EndX') GOTO 300
            GOTO 991
         ELSEIF(mark1 .EQ. '*') THEN
            BACKSPACE(ninp)
            READ(ninp,'(a)') comm80
            WRITE(nout,'(a)') comm80
         ENDIF
      ENDDO
 300  CONTINUE

*      WRITE(nout,*) '***************************'
*      WRITE(nout,*) '* Parser ReaDataX ends    *'
*      WRITE(nout,*) '***************************'
      RETURN
 990  WRITE(nout,*) '+++ ReaDataX: wrong index= ',index
      WRITE(   *,*) '+++ ReaDataX: wrong index= ',index
      STOP
      RETURN
 991  WRITE(nout,*) '+++ ReaDataX: wrong end of data '
      WRITE(   *,*) '+++ ReaDataX: wrong end of data '
      STOP
      END

