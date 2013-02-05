*/////////////////////////////////////////////////////////////////////////////////////////
*//                                                                                     //
*//                                                                                     //
*/////////////////////////////////////////////////////////////////////////////////////////
*
*  Class members:
      INTEGER       m_imax,         m_IdGen
      PARAMETER   ( m_imax = 10000, m_IdGen=6 )
      CHARACTER*132 m_GrandFormat
      PARAMETER (   m_GrandFormat ='(a1,a13,a1,a14,a1,a1,a1,f7.2,a1,f14.7,a1,f14.7,a1,f14.7,a1,a9,a1,a9,a1,a26,a1)') !
*-----------------------------------------------------------
*     $       '(a1,  ! Separator        position    1
*     $        a13,  ! observable label position    2-14
*     $        a1,   ! Separator        position   15
*     $        a14,  ! program name     position   16-29
*     $        a1,   ! Separator        position   30
*     $        a1,   ! Observable type  S or A     31
*     $        a1,   ! Separator        position   32
*     $        f7.1, ! CMS energy       position   33-39
*     $        a1,   ! Separator        position   40
*     $        g14.4,! prediction       position   41-54
*     $        a1,   ! Separator        position   55
*     $        g14.4,! Statistical err. position   56-69
*     $        a1,   ! Separator        position   70
*     $        g14.4,! Systematic error position   71-84
*     $        a1,   ! Separator        position   85
*     $        a9,   ! Prediction units position   86-94
*     $        a1,   ! Separator        position   95
*     $        a9,   ! Error      units position   96-104
*     $        a1,   ! Separator        position  105
*     $        a26,  ! Comment          position  106-131
*     $        a1,   ! Separator        position  132
*-----------------------------------------------------------
*
*  Histogram pointers xsections
      INTEGER     ix_Best,ix_Ceex1,ix_NoInt,ix_Ceex1Noint,ix_EEX3,ix_EEX2,ix_Pair !
      PARAMETER(  ix_Best  = 10203,     ix_Ceex1     = 10202,
     $            ix_NoInt = 10253,     ix_Ceex1Noint= 10252,
     $            ix_EEX3  = 10074,     ix_EEX2      = 10073,
     $            ix_Pair  = 10263)
*
      INTEGER     ix_IFI,ix_Phpr,ix_Ord3,ix_Int
      PARAMETER(  ix_IFI   = 20253,     ix_Phpr      = 20202,
     $            ix_Ord3  = 20073,     ix_Int       = 20953 )
*  Histogram pointers asymetries
      INTEGER     ia_Best,ia_Ceex1,ia_NoInt,ia_Ceex1Noint,ia_EEX3,ia_EEX2 !
      PARAMETER(  ia_Best  = 30203,     ia_Ceex1     = 30202,
     $            ia_NoInt = 30253,     ia_Ceex1Noint= 30252,
     $            ia_EEX3  = 30074,     ia_EEX2      = 30073)
*
      INTEGER     ia_IFI,ia_Phpr,ia_Ord3,ia_Int,ia_Work
      PARAMETER(  ia_IFI   = 40253,     ia_Phpr      = 40202,
     $            ia_Ord3  = 40073,     ia_Int       = 40953,     ia_Work  = 40853  )!
*  Histogram pointers semianalytical plots
      INTEGER     ivx_O3bestF, ivx_O3bestB, ivx_O3best, iva_O3best
      PARAMETER(  ivx_O3bestF= 51305, ivx_O3bestB= 52305,
     $            ivx_O3best = 53305, iva_O3best = 54305)
      INTEGER     ix_O3best,          ia_O3best
      PARAMETER(  ix_O3best  = 55305, ia_O3best  = 56305)
* Histograms for filling database
      INTEGER     ivx_SemiF, ivx_SemiB, ivx_Semi, iva_Semi
      PARAMETER(  ivx_SemiF= 61305,   ivx_SemiB= 62305,
     $            ivx_Semi = 63305,   iva_Semi = 64305) !
      INTEGER     ix_Semi,            ia_Semi
      PARAMETER(  ix_Semi  = 65305,   ia_Semi  = 66305)
cc      INTEGER     ivx_Semi,           ix_Semi
cc      PARAMETER(  ivx_Semi = 61305,   ix_Semi  = 62305)
*
      INTEGER     ix_ZFter,        ix_ZFNoInt,        ix_ZFInt,        ix_ZFIntExp !
      PARAMETER(  ix_ZFter=70001,  ix_ZFNoInt=70002,  ix_ZFInt=70003,  ix_ZFIntExp=70004  ) !
      INTEGER     ia_ZFter,        ia_ZFNoInt,        ia_ZFInt,        ia_ZFIntExp !
      PARAMETER(  ia_ZFter=70101,  ia_ZFNoInt=70102,  ia_ZFInt=70103,  ia_ZFIntExp=70104  ) !

      INTEGER     i_Dummy,         ix_Diff1,         ix_Diff2,         ix_Diff3 !
      PARAMETER(  i_Dummy = 10009, ix_Diff1 = 10001, ix_Diff2 = 10002, ix_Diff3 = 10003) !
      INTEGER                      ia_Diff1,         ia_Diff2,         ia_Diff3 !
      PARAMETER(                   ia_Diff1 = 10101, ia_Diff2 = 10102, ia_Diff3 = 10103) !

      INTEGER     i_Misc1,       i_Misc2,        i_Misc3
      PARAMETER(  i_Misc1=10091, i_Misc2=10092, i_Misc3=10093)
* Translation matrix from INDF to KF (-9 stand for forbiden)
      INTEGER     m_KFfromINDF(-1:11)
*                        nue, nu,  e, mu, tau,  u,  d,  c,  s,  t,  b, q's,Bhabha 
      DATA m_KFfromINDF / 12, 14, -9, 13,  15,  2,  1,  4,  3, -9,  5,   7,    -9/!
      SAVE m_KFfromINDF

* For plotting v-distribution
      DOUBLE PRECISION   m_loV,m_upV
      INTEGER            m_nbV
      PARAMETER(         m_loV=0d0, m_upV=1d0, m_nbV=100)
*------------------------------------------------------------------------------
      CHARACTER*13   m_QuLabel(9)
      DATA  m_QuLabel /
     $               'IAleph1      ',
     $               'IAleph2      ',
     $               'IDelphi1     ',
     $               'IDelphi2     ',
     $               'ILT1         ',
     $               'ILT2         ',
     $               'ILT3         ',
     $               'IOpal1       ',
     $               'IOpal2       '/
*------------------------------------------------------------------------------
      CHARACTER*13   m_Mulabel(29)
      DATA  m_MuLabel /
     $               'IAleph5      ',
     $               'IAleph6      ',
     $               'IDelphi5     ',
     $               'IDelphi6     ',
     $               'ILT9         ',
     $               'ILT10        ',
     $               'ILT11        ',
     $               'IOpal6       ',
     $               'IOpal7       ',
     $               'IOpal8       ',
     $               'IOpal9       ',
     $               'Aleph5       ',      ! REALISTIC START HERE
     $               'Aleph6       ',
     $               'Delphi4      ',
     $               'Delphi5      ',
     $               'LT9          ',
     $               'LT10         ',
     $               'Opal6        ',
     $               'Opal7        ',
     $               'ALEPH-12     ',       ! PHOTONIC START HERE 
     $               'ALEPH-15     ',
     $               'DELPHI9      ',
     $               'DELPHI12     ',
     $               'LT15         ',
     $               'LT18         ',
     $               'Opal14       ',
     $               'Opal15       ',
     $               'Opal20       ',
     $               'Opal21       '/
*------------------------------------------------------------------------------
      CHARACTER*13   m_TauLabel(11)
      DATA  m_TauLabel /
     $               'IAleph7      ',
     $               'IAleph8      ',
     $               'IDelphi7     ',
     $               'IDelphi8     ',
     $               'ILT12        ',
     $               'ILT13        ',
     $               'ILT14        ',
     $               'IOpal10      ',
     $               'IOpal11      ',
     $               'IOpal12      ',
     $               'IOpal13      '/
*------------------------------------------------------------------------------
*     Mark for plots
      CHARACTER*32 m_star,m_diamond,m_circle,m_ring,m_times,m_disc,m_plus,m_box,m_dot !
      PARAMETER (m_diamond ='\\makebox(0,0){\\Large $\\diamond$}')
      PARAMETER (m_star    ='\\makebox(0,0){\\Large $\\star$}')
      PARAMETER (m_circle  ='\\circle{30}')
      PARAMETER (m_ring    ='\\circle{20}')
      PARAMETER (m_times   ='\\makebox(0,0){\\Large $\\times$}')
      PARAMETER (m_disc    ='\\circle*{20}')
      PARAMETER (m_plus    ='\\makebox(0,0){\\Large $+$}')
      PARAMETER (m_box     ='\\makebox(0,0){\\Large $\\Box$}') !!! does not work???
      PARAMETER (m_dot     ='\\circle*{10}')
*
*//////////////////////////////////////////////////////////////////////////////
      DOUBLE PRECISION   m_xpar
      INTEGER            m_out,      m_FlagSem
      CHARACTER*26       m_Comment
      CHARACTER*6        m_Energy
*
      COMMON /c_PlotAll/
     $  m_xpar(m_imax),                 ! MC input data
     $  m_out,                          ! output unit
     $  m_Comment,                      ! ????
     $  m_Energy,                       ! Energy label
     $  m_FlagSem                       ! flag for semi on/off
      SAVE   /c_PlotAll/ 
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//                      End of CLASS  Plotter                               //
*//////////////////////////////////////////////////////////////////////////////
