

 ***************************************************************************
 *   DZface_Initialize, Interface to Dizet 6.xx                            *
 *      91.18700000                 Z mass                   amz        a1 *
 *     100.00000000                 Higgs mass               amh        a2 *
 *     175.00000000                 Top mass                 amtop      a3 *
 *               11                 KF code of beam          KFini      a5 *
 *               13                 KF of final fermion      KFfin      a6 *
 *                2                 IV code of beam          IVini      a7 *
 *                2                 IV of final fermion      IVfin      a8 *
 *                1                 EW box switch            ibox       a9 *
 *     128.86674175                 QED alfa inv. at Z       alfinv     a1 *
 *        .12500000                 QCD alfa at Z mass       alfQCD     a2 *
 ***************************************************************************

 DIZET flags, see routine Dizet for explanation:
  Ihvp = 1 Iamt4 = 4
  Iqcd = 3 Imoms = 1
 Imass = 0 Iscre = 0
 Ialem = 3 Imask = 0
 Iscal = 0 Ibarb = 2
 IFtjr = 1 Ifacr = 0
 IFact = 0 Ihigs = 0
 Iafmt = 1 Iewlc = 1
 Iczak = 1 Ihig2 = 1
 Iale2 = 3 Igfer = 2
 Iddzz = 1
    
 Alpha-QED   (MZ)  = .007759
 Alfa strong (MZ)  =   .1250
 Alfa strong (Mt)  =   .1132
zpard(20): QCD corr.fact. to Z-width (no b)  =   1.0334403654
zpard(25): QCD corr.fact. to Z-width (into b)=   1.0420641784
    
 zpar-matrix: standard output of dizet:
    zpar( 1)=   .03542132
    zpar( 2)=   .01157806
    zpar( 3)=   .22302485
    zpar( 4)=  1.16637000
    zpar( 5)=   .23120067
    zpar( 6)=   .23158158
    zpar( 7)=   .23158158
    zpar( 8)=   .23158158
    zpar( 9)=   .23147509
    zpar(10)=   .23134812
    zpar(11)=   .23147509
    zpar(12)=   .23134812
    zpar(13)=   .00000000
    zpar(14)=   .23294897
    zpar(15)=   .12500000
    zpar(16)=  1.00000000
    zpar(17)=  1.04206385
    zpar(18)=  1.05007810
    zpar(19)=  1.04145461
    zpar(20)=  1.03344037
    zpar(21)=  1.04207905
    zpar(22)=  1.04992189
    zpar(23)=  1.04145461
    zpar(24)=  1.03344037
    zpar(25)=  1.04206418
    zpar(26)=  1.05007843
    zpar(27)=  1.04201871
    zpar(28)=  1.02759398
    zpar(29)=  -.00002603
    zpar(30)=   .18568077


 ***************************************************************************
 *                     DZface_Initializion ended                           *
 ***************************************************************************

 pretabulation, basic LEP1 range
 a: i,ww=  0 0.100000000000000002E-01
 a: i,ww=  10 0.218741946339696160E-01
 a: i,ww=  20 0.478480390884785112E-01
 a: i,ww=  30 0.104663731987516523
 a: i,ww=  40 0.228943484461256785
 a: i,ww=  50 0.500795433928473144
 a: i,ww=  60 1.09544967935546911
 a: i,ww=  70 2.39620794979411489
 a: i,ww=  80 5.24151190772617248
 a: i,ww=  90 11.4653851645871736
 a: i,ww=  100 25.0796066643607674
 a: i,ww=  110 54.8596197519628674
 a: i,ww=  120 120.001000000000005
 pretabulation, near Z0:
 b: i,ww=  0 86.1884912220051831
 b: i,ww=  10 91.1869999999999976
 b: i,ww=  20 96.1855087779948121
 LEP2 energy zone: pretabulation starts
 c: i,ww=  0 120.001000000000005
 c: i,ww=  10 130.001000000000005
 c: i,ww=  20 140.001000000000005
 c: i,ww=  30 150.001000000000005
 c: i,ww=  40 160.001000000000005
 c: i,ww=  50 170.001000000000005
 c: i,ww=  60 180.001000000000005
 c: i,ww=  70 190.001000000000005
 c: i,ww=  80 200.001000000000005
 c: i,ww=  90 210.001000000000005
 c: i,ww=  100 220.001000000000005
 c: i,ww=  110 230.001000000000005
 c: i,ww=  120 240.001000000000005
 NLC energy range: pretabulation starts
 d: i,ww=  0 240.001000000000005
 d: i,ww=  10 340.000999999999976
 d: i,ww=  20 440.000999999999976
 d: i,ww=  30 540.000999999999976
 d: i,ww=  40 640.000999999999976
 d: i,ww=  50 740.000999999999976
 d: i,ww=  60 840.000999999999976
 d: i,ww=  70 940.000999999999976
 d: i,ww=  80 1040.00099999999998
 ... pretabulatin finished  now !
