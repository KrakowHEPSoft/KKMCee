*/////////////////////////////////////////////////////////////////////////////
* This program makes tables of E-Weak formfactors and stores them on disk.
*
* Execution>    make tables
* Execution>    make table.mu
* Execution>    make table.tau
* Execution>    make table.down
* Execution>    make table.up
* Execution>    make table.botom
* Execution>    make table.down-debug
*
*/////////////////////////////////////////////////////////////////////////////
* Notes:
* Re-initialization of EW package for new flavour is uncertain,
* therefore makefile invokes TabMain.exe several times, 
* hence input comes through standard input. 
* Input files are separate for each flavor, care must be
* taken to update top and Higgs mass everywhere.
* There is no danger, however, because program reading files checks if
* the input masses amz, amh, amtop were correctly aligned.
* On the other hand, since EW package may use standard output, the output
* table is writen in non-standard output. 
*/////////////////////////////////////////////////////////////////////////////
      PROGRAM tabmain
      IMPLICIT NONE
      INTEGER KFdown, KFup, KFstran, KFcharm, KFbotom, KFtop
      PARAMETER( 
     $     KFdown  = 1,   KFup    = 2,
     $     KFstran = 3,   KFcharm = 4,
     $     KFbotom = 5,   KFtop   = 6)
      INTEGER KFel,KFelnu,KFmu,KFmunu,KFtau,KFtaunu
      PARAMETER(
     $     KFel    = 11,  KFelnu  = 12,
     $     KFmu    = 13,  KFmunu  = 14,
     $     KFtau   = 15,  KFtaunu = 16)
      INTEGER KFfin
      INTEGER    iout
      INTEGER    iReset,imax,i
      PARAMETER (imax = 10000)
      DOUBLE PRECISION  xpar(imax)
      CHARACTER*40 DiskFile
*-------------------------------------------------------------------------------
      iout = 16
      OPEN(iout,FILE='./TabMain.output')
      CALL DZface_ReaDataX('./.KK2f_defaults', 1,imax,xpar)  ! reading general defaults
      CALL DZface_ReaDataX(    './input.data', 0,imax,xpar)  ! reading actual user input
*     Find active chanels

      DO i=401,416
         IF( xpar(i) .EQ. 1d0 ) THEN
            KFfin= i-400
* define disk file, check if KFfin is acceptable
            IF(    KFfin .EQ . 1) THEN 
               DiskFile= './table.down'
            ELSEIF(KFfin .EQ . 2)  THEN 
               DiskFile= './table.up'
            ELSEIF(KFfin .EQ . 3)  THEN 
               DiskFile= './table.stran'
            ELSEIF(KFfin .EQ . 4)  THEN 
               DiskFile= './table.charm'
            ELSEIF(KFfin .EQ . 5)  THEN 
               DiskFile= './table.botom'
            ELSEIF(KFfin .EQ .12)  THEN 
               DiskFile= './table.nue'
            ELSEIF(KFfin .EQ .13)  THEN 
               DiskFile= './table.mu'
            ELSEIF(KFfin .EQ .14)  THEN 
               DiskFile= './table.numu'
            ELSEIF(KFfin .EQ .15)  THEN 
               DiskFile= './table.tau'
            ELSEIF(KFfin .EQ .16)  THEN 
               DiskFile= './table.nutau'
            ELSE
               WRITE(*,   *) '#### STOP in TabMain, wrong KFfin=',KFfin
               WRITE(iout,*) '#### STOP in TabMain, wrong KFfin=',KFfin
               STOP
            ENDIF
            WRITE( *,'(a,a40)') 'Tabmain: DiskFile= ',DiskFile
 
            CALL DZface_Initialize( KFfin, xpar)  ! Set EW params and run Dizet
 
            CALL DZface_Tabluj  ! Calculate formfactor and store in tables
            CALL DZface_WriteFile(  DiskFile)     ! Write tables into DiskFile
         ENDIF
      ENDDO
      END
