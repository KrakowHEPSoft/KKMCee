*/////////////////////////////////////
*//  common blocks are absent here  //
*/////////////////////////////////////


      SUBROUTINE plotANG2(chak,title,id1,id2,ymin,ymax)
*//////////////////////////////////////////////////////////////////////////////
*// Purpose:    plot single angular distribution                             //
*//  id1,2      histograms MC (not cumulative)                               //
*//  ymin,ymax  uper lower limit in vertical scale of the plot               //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      CHARACTER*4 CHAK
      CHARACTER*80 TITLE
      INTEGER id1,id2
      REAL*8  ymin,ymax
*------------------------------------------------------------------------------
      INTEGER    IdGenYFS3
      PARAMETER (IdGenYFS3 = 6)
*------------------------------------------------------------------------------
*     Mark plots for plots
      CHARACTER*32 circle,disc,dot
      PARAMETER (circle  ='\\circle{6}')
      PARAMETER (disc    ='\\circle*{10}')
      PARAMETER (dot     ='\\circle*{6}')
      CHARACTER*16  fmtx,fmty
*------------------------------------------------------------------------------
      INTEGER idang1,idang2
*------------------------------------------------------------------------------
      CALL GLK_PlTitle(title)
      idang1  = id1 +500
      idang2  = id2 +500
      CALL GLK_RenHst(chak,IdGenYFS3,id1,idang1)
      CALL GLK_RenHst(chak,IdGenYFS3,id2,idang2)
*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      CALL GLK_Ymimax(idang1,  ymin,ymax)
      CALL GLK_Ymimax(idang2,  ymin,ymax)
      CALL GLK_idopt( idang1,   'ERRO')
      CALL GLK_idopt( idang2,   'ERRO')
      CALL GLK_Print( idang1)
      CALL GLK_Print( idang2)
*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      fmtx='f10.2'
      fmty='f10.3'
      CALL GLK_Plot2(idang1,' ',' ',circle  ,fmtx,fmty)
      CALL GLK_Plot2(idang2,'S',' ',disc  ,fmtx,fmty)
* purge local histos
      CALL GLK_Delet(idang1)
      CALL GLK_Delet(idang2)
      END

      SUBROUTINE plotANG1(chak,title,id1,ymin,ymax)
*//////////////////////////////////////////////////////////////////////////////
*// Purpose:    plot single angular distribution                             //
*//  id1,2      histograms MC (not cumulative)                               //
*//  ymin,ymax  uper lower limit in vertical scale of the plot               //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      CHARACTER*4 CHAK
      CHARACTER*80 TITLE
      INTEGER id1,id2
      REAL*8  ymin,ymax
*------------------------------------------------------------------------------
      INTEGER    IdGenYFS3
      PARAMETER (IdGenYFS3 = 6)
*------------------------------------------------------------------------------
*     Mark plots for plots
      CHARACTER*32 circle,disc,dot
      PARAMETER (circle  ='\\circle{6}')
      PARAMETER (disc    ='\\circle*{10}')
      PARAMETER (dot     ='\\circle*{6}')
      CHARACTER*16  fmtx,fmty
*------------------------------------------------------------------------------
      INTEGER idang1,idang2
*------------------------------------------------------------------------------
      CALL GLK_PlTitle(title)
      idang1  = id1 +500
      CALL GLK_RenHst(chak,IdGenYFS3,id1,idang1)
*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      CALL GLK_Ymimax(idang1,  ymin,ymax)
      CALL GLK_idopt( idang1,   'ERRO')
      CALL GLK_Print( idang1)
*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      fmtx='f10.2'
      fmty='f10.3'
      CALL GLK_Plot2(idang1,' ',' ',circle  ,fmtx,fmty)
* purge local histos
      CALL GLK_Delet(idang1)
      END


      SUBROUTINE plotAFB(chak,title,id1,id2,ymin,ymax)
*//////////////////////////////////////////////////////////////////////////////
*// Purpose:    calculation of AFB from 2 histos                             //
*//  chak       control of normalization                                     //
*//  title      title of the plot                                            //
*//  id1,2      histograms MC (not cumulative)                               //
*//  ymin,ymax  uper lower limit in vertical scale of the plot               //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      CHARACTER*4 CHAK
      CHARACTER*80 TITLE
      INTEGER id1,id2
      REAL*8  ymin,ymax
*------------------------------------------------------------------------------
      INTEGER    IdGenYFS3
      PARAMETER (IdGenYFS3 = 6)
*------------------------------------------------------------------------------
*     Mark plots for plots
      CHARACTER*32 circle,disc,dot
      PARAMETER (circle  ='\\circle{6}')
      PARAMETER (disc    ='\\circle*{10}')
      PARAMETER (dot     ='\\circle*{6}')
      CHARACTER*16  fmtx,fmty
*------------------------------------------------------------------------------
      INTEGER idcum1,idcum2
*------------------------------------------------------------------------------
      CALL GLK_PlTitle(title)
*! turn into cumulative normalized distribution or renormalize
      idcum1  = id1 +500
      idcum2  = id2 +500
      IF(chak .EQ. "CUMU") THEN
         CALL GLK_CumHis(     IdGenYFS3,id1,idcum1)
         CALL GLK_CumHis(     IdGenYFS3,id2,idcum2)
      ELSE
         CALL GLK_RenHst(chak,IdGenYFS3,id1,idcum1)
         CALL GLK_RenHst(chak,IdGenYFS3,id2,idcum2)
      ENDIF
      CALL GLK_Operat(idcum2,'/',idcum1,idcum2, 1d0, 1d0)
*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      CALL GLK_Ymimax(idcum2,  ymin,ymax)
      CALL GLK_idopt( idcum2,   'ERRO')
      CALL GLK_Print( idcum2)
*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      fmtx='f10.2'
      fmty='f10.3'
      CALL GLK_Plot2(idcum2,' ',' ',circle  ,fmtx,fmty)
* purge local histos
      CALL GLK_Delet(idcum1)
      CALL GLK_Delet(idcum2)
*------------
      END


      SUBROUTINE pldiff(chak,title,id,idan,ymin,ymax,xmag,ymag,born1)
*//////////////////////////////////////////////////////////////////////////////
*// Purpose:                                                                 //
*//          Plots difference of MC and semianalytical results               //
*//  chak       control of normalization                                     //
*//  title      title of the plot                                            //
*//  id         histogram MC (not cumulative)                                //
*//  idan       histogram analytical non-cumulative or cumulative            //
*//  ymin,ymax  uper lower limit in vertical scale of the plot               //
*//  xmag       magnification for difference MC-Analit.                      //
*//  ymag       magnification for ID and IDAN.                               //
*//  born1      normalization reference x-section                            //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      CHARACTER*4 CHAK
      CHARACTER*80 TITLE
      INTEGER id,idan
      REAL*8  ymin,ymax,xmag,ymag,born1
*------------------------------------------------------------------------------
      INTEGER    IdGenYFS3
      PARAMETER (IdGenYFS3 = 6)
      REAL*8  born
      INTEGER idcum,idan2,idiff
*------------------------------------------------------------------------------
*     Mark plots for plots
      CHARACTER*32 circle,disc,dot
      PARAMETER (circle  ='\\circle{6}')
      PARAMETER (disc    ='\\circle*{10}')
      PARAMETER (dot     ='\\circle*{6}')
      CHARACTER*16  fmtx,fmty
*------------------------------------------------------------------------------
      born=abs(born1)
      CALL GLK_PlTitle(title)
*! turn into cumulative normalized distribution or renormalize
      idcum  = id+500
      IF(chak .EQ. "CUMU") THEN
         CALL GLK_CumHis(     IdGenYFS3,id,idcum)
      ELSE
         CALL GLK_RenHst(chak,IdGenYFS3,id,idcum)
      ENDIF
*! divide by Born
      idan2 = idan+500
      IF( born1 .EQ. 0D0) GOTO 901
      CALL GLK_Operat(idcum,'+',idcum,idcum, 0d0, 1/born)
      CALL GLK_Operat(idan ,'+',idan ,idan2, 0d0, 1/born)
      IF(born1 .LT. 0d0) THEN
*        ! optionaly subtract one, obsolete!!!
         WRITE(*,*) 'WARNIG pldiff: negative BORN=',born
      ENDIF
*! the difference MC-Analit. times XMAG
      idiff =  id+502
      CALL GLK_Operat(idcum,'-',idan2,idiff, xmag, xmag)
*! shrink/magnify original MC and Analitical curves
      CALL GLK_Operat(idcum,'+',idcum,idcum, 0d0, ymag)
      CALL GLK_Operat(idan2,'+',idan2,idan2, 0d0, ymag)
* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      CALL GLK_Ymimax(idcum,  ymin,ymax)
      CALL GLK_Ymimax(idan2,  ymin,ymax)
      CALL GLK_Ymimax(idiff,  ymin,ymax)
      CALL GLK_idopt(idcum,   'ERRO')
      CALL GLK_idopt(idiff,   'ERRO')
      CALL GLK_Print(idcum)
      CALL GLK_Print(idan2)
      CALL GLK_Print(idiff)
* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      fmtx='f10.2'
      fmty='f10.3'
      CALL GLK_Plot2(idcum,' ',' ',dot    ,fmtx,fmty)
      CALL GLK_Plot2(idan2,'S','*',circle ,fmtx,fmty)
      CALL GLK_Plot2(idiff,'S','*',dot    ,fmtx,fmty)
* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
* purge local histos
      CALL GLK_Delet(idcum)
      CALL GLK_Delet(idan2)
      CALL GLK_Delet(idiff)
      RETURN
 901  WRITE(6,*) ' PLDIFF: wrong born!!!!'
      END

      SUBROUTINE pld2mc(chak,title,id1,id2,ymin,ymax,xmag,ymag,born1)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//  TITLE      title of the plot                                            //
*//  ID1         histogram MC (not cumulative)                               //
*//  ID2         histogram MC (not cumulative)                               //
*//  YMIN,YMAX  uper lower limit in vertical scale of the plot               //
*//  XMAG       magnification for difference MC-Analit.                      //
*//  YMAG       magnification for ID and IDAN.                               //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INTEGER id1,id2
      REAL*8  ymin,ymax,xmag,ymag,born1
      CHARACTER*4 chak
      CHARACTER*80 title
*------------------------------------------------------------------------------
      INTEGER    IdGenYFS3
      PARAMETER (IdGenYFS3 = 6)
      INTEGER    idan,idan2,idiff,idcum1,idcum2
      REAL*8     born
*------------------------------------------------------------------------------
*     Mark plots for plots
      CHARACTER*32 circle,disc,dot
      PARAMETER (circle  ='\\circle{6}')
      PARAMETER (disc    ='\\circle*{10}')
      PARAMETER (dot     ='\\circle*{6}')
      CHARACTER*16  fmtx,fmty
*------------------------------------------------------------------------------
      born=abs(born1)
      CALL GLK_PlTitle(title)
*! turn into cumulative normalized distribution or renormalize
      idcum1  = id1+500
      idcum2  = id2+500
      IF(chak .EQ. "CUMU") THEN
         CALL GLK_CumHis(IdGenYFS3,id1,idcum1)
         CALL GLK_CumHis(IdGenYFS3,id2,idcum2)
      ELSE
         CALL GLK_RenHst(chak,IdGenYFS3,id1,idcum1)
         CALL GLK_RenHst(chak,IdGenYFS3,id2,idcum2)
      ENDIF
*! divide by Born
      idan2 = idan+500
      IF( born1 .EQ. 0D0) GOTO 901
      CALL GLK_Operat(idcum1,'+',idcum1,idcum1, 0d0, 1/born)
      CALL GLK_Operat(idcum2,'+',idcum2,idcum2, 0d0, 1/born)
*! optionaly subtract one
      IF(born1 .LT. 0d0) THEN
         WRITE(*,*) 'WARNIG pldiff: negative BORN=',born
      ENDIF
*! the difference MC-Analit. times XMAG
      idiff =  id1+502
      CALL GLK_Operat(idcum1,'-',idcum2,idiff, xmag, xmag)
*! shrink/magnify original MC and Analitical curves
      CALL GLK_Operat(idcum1,'+',idcum1,idcum1, 0d0, ymag)
      CALL GLK_Operat(idcum2,'+',idcum2,idcum2, 0d0, ymag)
* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      CALL GLK_Ymimax(idcum1,    ymin,ymax)
      CALL GLK_Ymimax(idcum2,    ymin,ymax)
      CALL GLK_Ymimax( idiff,    ymin,ymax)
      CALL GLK_idopt(idcum1,   'ERRO')
      CALL GLK_idopt(idcum2,   'ERRO')
      CALL GLK_idopt( idiff,   'ERRO')
      CALL GLK_Print(idcum1)
      CALL GLK_Print(idcum2)
      CALL GLK_Print( idiff)
* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      fmtx='f10.2'
      fmty='f10.3'
      CALL GLK_Plot2(idcum1,' ','B',dot    ,fmtx,fmty)
      CALL GLK_Plot2(idcum2,'S','B',circle ,fmtx,fmty)
      CALL GLK_Plot2(idiff ,'S','*',disc   ,fmtx,fmty)
* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c* purge local histos
      CALL GLK_Delet(idcum1)
      CALL GLK_Delet(idcum2)
      CALL GLK_Delet( idiff)
      RETURN
 901  WRITE(6,*) ' pld2mc: wrong born!!!!'
      END



      SUBROUTINE plt1mc(chak,title,id1,ymin,ymax,born1)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*// Absolute normalization                                                   //
*//  TITLE      title of the plot                                            //
*//  ID1         histogram MC (not cumulative)                               //
*//  YMIN,YMAX  uper lower limit in vertical scale of the plot               //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INTEGER id1
      REAL*8  ymin,ymax,born1
      CHARACTER*4 CHAK
      CHARACTER*80 TITLE
*------------------------------------------------------------------------------
      INTEGER    IdGenYFS3
      PARAMETER (IdGenYFS3 = 6)
      REAL*8  born
      INTEGER idcum1
*------------------------------------------------------------------------------
*     Mark plots for plots
      CHARACTER*32 circle,disc,dot
      PARAMETER (circle  ='\\circle{6}')
      PARAMETER (disc    ='\\circle*{10}')
      PARAMETER (dot     ='\\circle*{6}')
      CHARACTER*16  fmtx,fmty
*------------------------------------------------------------------------------
      born=abs(born1)
      CALL GLK_PlTitle(title)
*! turn into cumulative normalized distribution or renormalize
      idcum1  = id1+500
      IF(chak .EQ. "CUMU") THEN
         CALL GLK_CumHis(IdGenYFS3,id1,idcum1)
      ELSE
         CALL GLK_RenHst(chak,IdGenYFS3,id1,idcum1)
      ENDIF
*! divide by Born
      IF( born1 .EQ. 0D0) GOTO 901
      CALL GLK_Operat(idcum1,'+',idcum1,idcum1, 0d0, 1/born)
*! optionaly subtract one
      IF(born1 .LT. 0d0) THEN
         WRITE(*,*) 'WARNIG pldiff: negative BORN=',born
      ENDIF
* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      CALL GLK_Ymimax(idcum1,    ymin,ymax)
      CALL GLK_idopt(idcum1,   'ERRO')
      CALL GLK_Print(idcum1)
* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      fmtx='f10.2'
      fmty='f10.3'
      CALL GLK_Plot2(idcum1,' ','B',dot    ,fmtx,fmty)
* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
* purge local histos
      CALL GLK_Delet(idcum1)
      RETURN
 901  WRITE(6,*) ' plt1mc: wrong born!!!!'
      END

      SUBROUTINE plt1an(title,id1,id2,ymin,ymax,born1)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//  TITLE      title of the plot                                            //
*//  ID1         histogram MC (not cumulative)                               //
*//  YMIN,YMAX  uper lower limit in vertical scale of the plot               //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INTEGER   id1,id2
      REAL*8    ymin,ymax,born1
      CHARACTER*4 CHAK
      CHARACTER*80 TITLE
      REAL*8    Born
*------------------------------------------------------------------------------
*     Mark plots for plots
      CHARACTER*32 circle,disc,dot
      PARAMETER (circle  ='\\circle{6}')
      PARAMETER (disc    ='\\circle*{10}')
      PARAMETER (dot     ='\\circle*{6}')
      CHARACTER*16  fmtx,fmty
*------------------------------------------------------------------------------
      born=abs(born1)
      CALL GLK_PlTitle(title)
*! divide by Born
      IF( born1 .EQ. 0D0) GOTO 901
      CALL GLK_Operat(id1,'+',id1,id2, 0d0, 1/born)
* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      CALL GLK_Ymimax(id2,    ymin,ymax)
      CALL GLK_Print( id2)
* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      fmtx='f10.2'
      fmty='f10.3'
      CALL GLK_Plot2(id2,' ','B',dot    ,fmtx,fmty)
* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      RETURN
 901  WRITE(6,*) ' plt1an: wrong born!!!!'
      END

      SUBROUTINE pld2an(title,id1,id2,ymin,ymax,xmag,ymag,born1)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*// Absolute normalization                                                   //
*//  TITLE      title of the plot                                            //
*//  id1        histogram, analytical                                        //
*//  id2        histogram, analytical                                        //
*//  ymin,ymax  uper lower limit in vertical scale of the plot               //
*//  xmag       magnification for difference Analit-Analit.                  //
*//  ymag       magnification for ID and id1,id2                             //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INTEGER id1,id2
      REAL*8  ymin,ymax,xmag,ymag,born1
      CHARACTER*4 chak
      CHARACTER*80 title
*------------------------------------------------------------------------------
      INTEGER    IdGenYFS3
      PARAMETER (IdGenYFS3 = 6)
      INTEGER    idan,idan2,idiff,idcum1,idcum2
      REAL*8     born
*------------------------------------------------------------------------------
*     Mark plots for plots
      CHARACTER*32 circle,disc,dot
      PARAMETER (circle  ='\\circle{6}')
      PARAMETER (disc    ='\\circle*{10}')
      PARAMETER (dot     ='\\circle*{6}')
      CHARACTER*16  fmtx,fmty
*------------------------------------------------------------------------------
      born=abs(born1)
      CALL GLK_PlTitle(title)
*! turn into cumulative normalized distribution or renormalize
      idcum1  = id1+500
      idcum2  = id2+500
*! divide by Born
      idan2 = idan+500
      IF( born1 .LE. 0D0) GOTO 901
      CALL GLK_Operat(id1,'+',id1,idcum1, 0d0, 1/born)
      CALL GLK_Operat(id2,'+',id2,idcum2, 0d0, 1/born)
*! the difference Analit-Analit. times XMAG
      idiff =  id1+600
      CALL GLK_Operat(idcum1,'-',idcum2,idiff, xmag, xmag)
*! shrink/magnify original MC and Analitical curves
      CALL GLK_Operat(idcum1,'+',idcum1,idcum1, 0d0, ymag)
      CALL GLK_Operat(idcum2,'+',idcum2,idcum2, 0d0, ymag)
* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      CALL GLK_Ymimax(idcum1,    ymin,ymax)
      CALL GLK_Ymimax(idcum2,    ymin,ymax)
      CALL GLK_Ymimax( idiff,    ymin,ymax)
      CALL GLK_idopt(idcum1,   'ERRO')
      CALL GLK_idopt(idcum2,   'ERRO')
      CALL GLK_idopt( idiff,   'ERRO')
      CALL GLK_Print(idcum1)
      CALL GLK_Print(idcum2)
      CALL GLK_Print( idiff)
* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      fmtx='f10.2'
      fmty='f10.3'
      CALL GLK_Plot2(idcum1,' ','B',dot    ,fmtx,fmty)
      CALL GLK_Plot2(idcum2,'S','B',circle ,fmtx,fmty)
      CALL GLK_Plot2(idiff ,'S','*',disc   ,fmtx,fmty)
* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c* purge local histos
      CALL GLK_Delet(idcum1)
      CALL GLK_Delet(idcum2)
      CALL GLK_Delet( idiff)
      RETURN
 901  WRITE(6,*) ' pld2mc: wrong born!!!!',born
      END





      SUBROUTINE pltech(chak,title,idmc2,idmc21,idan2,idan21,ymin,ymax,xmag,ymag,born)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//  It prints: MC O(alf2), Analytical O(alf2) and difference of two,        //
*//             MC O(alf2-alf1), Analytical O(alf2-alf1)                     //
*//   Called in susini and karfin                                            //
*//  CHAK       if CHAK="CUMU" then IDMC2 IDMC21 are cumulated               //
*//             and normalized in nanobarns using CUMHIS                     //
*//             else RENHST is invoked                                       //
*//  TITLE      title of the plot                                            //
*//  IDMC2      histogram MC O(alf2)                                         //
*//  IDMC21     histogram MC O(alf2-alf1)                                    //
*//  IDAN2      histogram Analytical O(alf2)                                 //
*//  IDAN21     histogram Analytical O(alf2-alf1)                            //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      CHARACTER*4 CHAK
      CHARACTER*80 TITLE
      REAL*8    ymin,ymax,xmag,ymag,born
      INTEGER   idmc2,idmc21,idan2,idan21
*---
      INTEGER    IdGenYFS3
      PARAMETER (IdGenYFS3 = 6)
*---
      INTEGER    jmonte,jphysm,jtechn,janalt,jphysa
*------------------------------------------------------------------------------
      CALL GLK_PlTitle(title)
*! turn into cumulative normalized distribution or renormalize
      jmonte   = idmc2  +200
      jphysm   = idmc21 +200
      IF(CHAK.eq."CUMU") THEN
         CALL GLK_CumHis(IdGenYFS3,idmc2, jmonte)   ! MC O(alf2)
         CALL GLK_CumHis(IdGenYFS3,idmc21,jphysm)   ! MC O(alf2-alf1)
      ELSE
         CALL GLK_RenHst(chak,IdGenYFS3,idmc2, jmonte)
         CALL GLK_RenHst(chak,IdGenYFS3,idmc21,jphysm)
      ENDIF
*! Technical precision MC-Analit. O(alf2) times XMAG/BORN
      jtechn =  idmc2+500
      CALL GLK_Operat(jmonte,'-',idan2,   jtechn, xmag/born, xmag/born)
*! Physical precision  O(alf2-alf1) times XMAG
      jphysa =   idan21+500
      CALL GLK_Operat(jphysm,'+',jphysm,  jphysm,  0d0, xmag/born)
      CALL GLK_Operat(idan21,'+',idan21,  jphysa,  0d0, xmag/born)

*! Originals divided by Born and (de)magnify by YMAG
      janalt =idan2+500
      CALL GLK_Operat(jmonte, '+',jmonte, jmonte, 0d0, ymag/born)
      CALL GLK_Operat(idan2,  '+', idan2, janalt, 0d0, ymag/born)
* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      CALL GLK_Ymimax(jmonte,    ymin,ymax)
      CALL GLK_Ymimax(janalt,    ymin,ymax)
      CALL GLK_Ymimax(jtechn,    ymin,ymax)
      CALL GLK_Ymimax(jphysm,    ymin,ymax)
      CALL GLK_Ymimax(jphysa,    ymin,ymax)
      CALL GLK_idopt(jmonte,   'ERRO')
      CALL GLK_idopt(jtechn,   'ERRO')
      CALL GLK_idopt(jphysm,   'ERRO')
      CALL GLK_Print(jmonte)
      CALL GLK_Print(janalt)
      CALL GLK_Print(jtechn)
      CALL GLK_Print(jphysm)
      call GLK_Print(jphysa)
* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
* technical
      CALL GLK_PlSet('DMOD',2d0)
      CALL GLK_Plot(jmonte,   ' ','B',0)
      CALL GLK_PlSet('DMOD',2d0)
      CALL GLK_Plot(janalt,   'S','*',0)
      CALL GLK_PlSet('DMOD',1d0)
      CALL GLK_Plot(jtechn  , 'S','*',0)
* physical
      CALL GLK_PlSet('DMOD',2d0)
      CALL GLK_Plot(jphysm,   ' ','B',0)
      CALL GLK_PlSet('DMOD',2D0)
      CALL GLK_Plot(jphysa,   'S','*',0)
* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      CALL GLK_Delet(jmonte)
      CALL GLK_Delet(janalt)
      CALL GLK_Delet(jphysm)
      CALL GLK_Delet(jtechn)
      CALL GLK_Delet(jphysa)
      END
