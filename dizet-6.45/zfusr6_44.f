*==========================================================================*
*                                                                          *
* This is ZFITTER version 6.44, 15 January 2013
*
* Update 15 Jan 2013:
* For the program, description, updates etc. see (January 2013):            *
* http://zfitter.com (responsible T. Riemann)
* and 
* http://zfitter.desy.de (S. Riemann, under responsibility of 
* DESY Board of Directors)
*
* For program licence and 'conditions of use' see: 
* http://cpc.cs.qub.ac.uk/licence/licence.html
* and the webpages; in case please contact the spokesperson
*----------------------------------------------------------------------*
* VERSION 6.44 (15 January 2013)
* 
* new flag IBAIKOV is used, is set fixed inside subroutine QCDCOF
* for the modifications see: "mod. 14 Jan 2013"
* The changes are made in order to implement the QDC corrections 
* as described in:
* P. Baikov, K. Chetyrkin, J. Kuehn, J. Rittinger
* "Complete QCD Corrections to Hadronic Z-Decays in Order alpha_s^4"
* following exactly: arXiv:1201.5804v3 (2 May 2012)
* published as Phys.Rev.Lett. 108 (2012) 222003
* contact: Tord Riemann tordriemann@gmail.com
*---------------------------------------------------------------------
*
* ZFITTER group since February 2013:
* A. Akhundov, A. Arbuzov, B. Bardin, P. Christova, L. Kalinovskaya, 
* A. Olshevsky, S. Riemann, T. Riemann                            
*
*==========================================================================*
*                                                                          *
* The authors of ZFITTER are:                                              *
*--------------------------------------------------------------------------*
* D. Bardin,               M. Bilenky (1987-1994), A. Chizhov (1987-1991), *
* P. Christova (since 1999),  O. Fedorenko (1990),                         *
* M. Jack (since 1999),    L. Kalinovskaya (since 1997),                   *
* A. Olshevsky,            S. Riemann,             T. Riemann,             *
* M. Sachwitz (1987-1991), A. Sazonov (1987-1991), Yu. Sedykh (1989-1991), *
* I. Sheer (1991-1992),    L. Vertogradov (1990)                           *
*                                                                          *
*==========================================================================*
*                                                                          *
* ZFITTER relies on the following articles:                                *
*--------------------------------------------------------------------------*
*  Update description: 
* ZFITTER: a semi-analytical program for fermion pair production in 
* e+ e- annihilation, from version 6.21 to version 6.42
* [hep-ph/0507146]
* Comput. Phys. Commun. 174 (2006) 728-758. 
*
* Program description:
* D. Bardin et al.:                                                        *
*  CERN-TH.6443/92 (March 1992)    [hep-ph/9412201]                        *
*  old description of package                                              *
* New descripton of the package is published in CPC                        *
* D. Bardin, M. Bilenky, M. Jack, P. Christova,                            *
* L. Kalinovskaya, A. Olshevsky, S. Riemann, T. Riemann:                   *
*  Comput.Phys.Commun. 133 (2001) 229-395,2001                             *
* See also,                                                                *
* D. Bardin et al.:                                                        *
*  Z. Physik C44 (1989) 493        [KEK library entry 8906215]             *
*  -- weak library                                                         *
*  Nucl. Phys. B351 (1991) 1       [hep-ph/9801208]                        *
*  -- angular distribution and general QED approach                        *
*  Phys. Letters B 255 (1991) 290  [hep-ph/9801209]                        *
*  -- QED with s'- and  acceptance cut                                     *
* P. Christova, M. Jack, T. Riemann:                                       *
*  DESY 99-015                     [hep-ph/9902408]                        *
*  -- QED corr's. with acollinearity and fermion energy cut                *
*     unpublished so far:                                                  *
*  -- QED with acollinearity and fermion energy and acceptance cut         *
*--------------------------------------------------------------------------*
* also:                                                                    *
* A. Akhundov, D. Bardin, T. Riemann, Nucl. Phys. B271 (1986) 1            *
* (Weak corrs. to Z decay)                                                 *
* Many further references may be found in the description                  *
* in preparation.                                                          *
*==========================================================================*
*                                                                          *
* The package ZFITTER calls the following packages:                        *
*--------------------------------------------------------------------------*
* DIZET  - one loop weak library -                                         *
* Authors:                                                                 *
* A. Akhundov (1985-1989), D. Bardin,           M. Bilenky (1987-1994),    *
* P. Christova,            L. Kalinovskaya (since 1997),                   *
* S. Riemann,              T. Riemann,          M. Sachwitz (1987-1991),   *
* H. Vogt (1989)                                                           *
*--------------------------------------------------------------------------*
* m2tcor - two loop weak library -                                         *
* Author:                                                                  *
* G. Degrassi                                                              *
*--------------------------------------------------------------------------*
* hadr5  - hadronic vacuum polarization -                                  *
* Author:                                                                  *
* F. Jegerlehner                                                           *
*--------------------------------------------------------------------------*
* BHANG  - a supplement for Bhabha scattering -                            *
* Author:                                                                  *
* M. Bilenky                                                               *
*--------------------------------------------------------------------------*
* bkqcdl - a QCD-library based on                                          *
*          functions V(r), A(r), F(x) provided by                          *
* B. Kniehl                                                                *
*--------------------------------------------------------------------------*
* pairho - higher order pair corrections to annihilation                   *
* Author:                                                                  *
* A.B. Arbuzov (hep-ph/9907500)                                            *
*==========================================================================*
*
* ZFITTER v. 6_43 (17-06-2008):                                            *
*------------------------------
* a bug corrected in dizet6_43.f (see there for details), no numerical 
* influence:  a negligible effect, 0.9 MeV (!) on the central MH value
*
* -----------------------------                                            *
*                                                                          *
* ZFITTER v. 6_42 (18-05-2005):                                            *
* -----------------------------                                            *
* Correction of a bug in the new implementation of universal 2-loop        *
* corrections for  e+e- -> bb  in version 6.41. There was an error in the  *
* treatment of the form factor rho [see (3.3.1) in manual for definiton].  *
* The removal of this error only affects the computation of the bb cross-  *
* section, whereas the b-asymmetries were produced correctly by version    *
* 6.41 and stay unaltered.                                                 *
*==========================================================================*
*                                                                          *
* ZFITTER v. 6_41 (15-10-2004):                                            *
* -----------------------------                                            *
* 1. Changes made to function ROKANC in order to accomodate proper         *
*    treatment of new 2-loop corrections for kappa (sw_eff).               *
*    Also the treatment of the bb final state has been altered:            *
*    Previously the corrections using ROKANC (ZUTHSM) were all calculated  *
*    up to 1-loop level for e+e- -> bb, since 2-loop corrections for the   *
*    Zbb vertex are not known. This means that also the (known) universal  *
*    corrections to the initial state Zee vertex were not included. This   *
*    results in an inconsistency with respect to the usage of ZUXSA, which *
*    uses effective couplings and ALWAYS includes higher order corrections *
*    to the initial state Zee vertex. In the new version, also ROKANC      *
*    includes higher order corrections to the initial state vertex for the *
*    bb final state. This change does affect all AMT4>3. For AMT4=6, the   *
*    new 2-loop corrections to kappa are included. The corrections to the  *
*    Zbb vertex and the boxes are calculated at 1-loop level as before.    *
* 2. A bug in the transfer of flags in version 6.40 has been corrected     * 
*==========================================================================*
*                                                                          *
* ZFITTER v. 6_40 (29-04-2004):                                            *
* -----------------------------                                            *
* 1. Last results by M. Awramik, M. Czakon, A. Freitas and G. Weiglein     *
*    [hep-ph/0311148] on M_w are implemented at AMT4=6.                    *
*    I.e. if AMT4=4 then the old treatment of the leading and next-to-lead.*
*    two loop EWRC by Degrassi et al. is used (as before) and              *
*         if AMT4=5 then only for the calculation of M_w the result        *
*         including complete 2-loop fermionic corrections is used          *
*    For AMT4=6 M_w is calculated using the new result including complete  *
*     fermionic and bosonic 2-loop corrections and leading 3-loop corr.    *
*     by Faisst et al. [hep-ph/0302275]                                    *
*     AMT4=6 also activates the new sw_eff result (see below)              *
*    For the rho formfactor of Zff interactions still the package of       *
*    Degrassi et al. is used                                               *
* 2. New result by M. Awramik, M. Czakon, A. Freitas and G. Weiglein       *
*    [hep-ph/0407317] for kappa (sw_eff) implemented at AMT4=6             *
*    This result incorporates complete fermionic 2-loop corrections and    *
*    leading electroweak 3-loop corrections by Faisst et al.               *
*    I.e. for AMT4=6 the newest results for M_w and kappa(sw_eff) are used *
* 3. Flag DMWW extended to simulate some theoretical uncertainties in      *
*    new result for MW                                                     *
*    DMWW=0 ! +-1 external flag with values 0,+-1 [sigma]                  *
* 4. New flag DSWW introduced to simulate theoretical uncertainties in     *
*    new 2-loop result for kappa(sw_eff)                                   *
*    DSWW=0 ! +-1 external flag with values 0,+-1 [sigma]                  *
*==========================================================================*
*                                                                          *
* ZFITTER v. 6_36 (21-06-2001):                                            *
* -----------------------------                                            *
* 1. Some improvements for INDF=-1 following common work with SJ and ZW    *
*                                  at CERN in March 2001.                  *
*    Adding WEAK=-1 user option, see description in ZUFLAG.                *
* 2. Some improvements for AMT4= 5 following common work with GW in March. *
*    AMT4= 5 implements complete for Mw (and NON-complete for \sin^2_eff)  *
*    calculations by A.Freitas, W.Hollik, W.Walter and G.Weiglein          *
*                    CERN-TH/2000-194; hep-ph/0007091                      *
* 3. Flag AFMT is extended:                                                *
*     AFMT       0  WITHOUT AFMT CORRECTION                                *
*                1   WITH   AFMT CORRECTION O(Gf*m^2_t*als^2)              *
*                2   WITH   ChKS CORRECTION + O(Gf*M^2_z*als^2+log(m^2_t)  *
* (Now default)  3   WITH   ChKS CORRECTION + O(Gf*M^2_z/m^2_t*als^2)      *
* ChKS=Chetyrkin, Kuehn, Steinhouser (AFMT-extended) from hep-ph/9504413   *
* 4. Martin's flag to simulate some theoretical uncertainties in MW        *
*    DMWW=0 ! +-1 new external flag with values 0,+-1 [sigma]              *
*==========================================================================*
*                                                                          *
* ZFITTER v. 6_35 (18-02-2001):                                            *
* -----------------------------                                            *
* 1. Some improvements in ZU_APV branch (to be described in a paper)       *
* 2. Martin Gruenewald fixes for AMT4=5                                    *
* 3. Flag TUPV to simulate EW TUs for APV                                  *
*    TUPV=1 default                                                        *
*    TUPV=2/3 to vary (actually, 1 and 3 are sufficient)                   *
*==========================================================================*
*                                                                          *
* ZFITTER v. 6_34 (26-01-2001):                                            *
* -----------------------------                                            *
*                                                                          *
* 1. Consistent treatment of APV in OMS is added.                          *
*    Accessible via user SUBROUTINE ZU_APV.                                *
*    (An OMS-consistent re-calculation of APV is realized                  *
*     by D.Bardin, P.Christova and L.Kalinovskaya, January 2001)           *
*                                                                          *
* 2. IBA cross-section for the process ee-->\nu_e\nu_e is added.           *
*    (By D.Bardin and T.Riemann within a work together                     *
*                               with S.Jadach and Z.Was, to be published.) *
*    Accessible via COSCUT with the aid of zfEENN_34.f                     *
*    Governed by flag ENUE:                                                *
*    ENUE=-1 s-channel only                                                *
*    ENUE= 0 s- and t-channels                                             *
*    ENUE=+1 s+t complete with s-t interference                            *
*                                                                          *
*==========================================================================*
*                                                                          *
* ZFITTER v. 6_33 (12-12-2000) this version was not officially released    *
* ---------------------------------------------------------------------    *
* 1. Several tiny bugs, mostly in comments are fixed.                      *
* 2. Following Martin Gruenewald suggestion the DIZET v.6.30 allowing for  *
*    a fit of V_tb is joined with ZFITTER v.6.30.                          *
*    This had, some consequences, since some CALL's had to be changed.     * 
*    In some cases I followed the principle to replace in the argument     *
*    list of a SUBROUTINE (i.g. ZDIZET, etc...):                           *
*    *,DAL5H,ALFAS,* to *,DAL5H,V_TB,ALFAS,*                               *
*    in order to accomodate V_tb on the same footing as all the other      *
*    fitting parameters.                                                   *
*    Note, however, that DIZET argument list was changed completely.       *
*    For the user INTERFACES (ZUWEAK, ZUTHSM, ZUATSM, ZUTPSM, ZULRSM)      *
*    wrapper routines were added as Martin Gruenewald suggested:           *
*                                                                          *
*       SUBROUTINE ZUWEAK(ZMASS,TMASS,HMASS,DAL5H,ALFAS)                   *
* *     ========== =====================================                   *
* *     ZUWEAK-WRAPPER                                                     *
* *     ==============                                                     *
*       IMPLICIT REAL*8(A-Z)                                               *
*       V_TB=1D0                                                           *
*       CALL  ZVWEAK(ZMASS,TMASS,HMASS,DAL5H,V_TB,ALFAS)                   *
*       RETURN                                                             *
*       END                                                                *
*                                                                          *
* 3. Last results by A. Freitas, W. Holik, W. Walter, and G. Weiglein      *
*    [hep-ph/0007091] on M_w are implemented at AMT4=5.                    *
*    I.e. if AMT4=4 then the old treatment of the leading and next-to-lead.*
*    two loop EWRC is used (as before) and                                 *
*         if AMT4=5 then only for the calculation of M_w new results of    *
*         above quoted gentlemen are used.                                 *
*    For Z-decay and EWFF of the process ee->ff one used Degrassi package. *
*    For this reason many IFs *AMT4.EQ.4* became *AMT4.GE.4*               *
*                                                                          *
* 4. About AZ1 difference in ZF6.23 and lower compared to ZF6.30 and higher*
*    (DB forgot to place a comment in the description of v.6.30 below.)    *
*    This is long story showing up as a tiny difference for f-b asymmetry. *
*    A flag ASCR is added. If ASCR=0 then backcompatibulity with ZF6.23    *
*    and below (not all below, maybe 5-th and 6-th families).              *
*    If ASCR=1 then new treatment (the same as very old), this is default. *
*                                                                          *
*==========================================================================*
*                                                                          *
* ZFITTER v. 6_30 (23-06-2000):                                            *
* -----------------------------                                            *
* 1. A new option governed by new flag FUNA is implemented.                *
*    FUNA=0 old treatment                                                  *
*    FUNA=1 new treatment                                                  *
*      This is a new treatment of the second order ISR QED corrections,    *
*      in the presence of angular acceptance cut ANG0,ANG1.                *
*      It is based on a new calculation by A. Arbuzov (to appear as hep-ph)*
*      It is compatible with the use of ICUT=1,2,3.                        *
*                                                                          *
* 2. Meaning of INTF flag is extended in order to accomodate new           *
*    implementation of an exponentiation of IFI QED corrections            *
*    realized by A. Arbuzov (to appear as hep-ph).                         *
*    INTF=0,1 old options                                                  *
*    INTF=2   exponentiated IFI                                            *
*                                                                          *
* 3. Final state pair production corrections are implemented (A. Arbuzov); *
*    the option is governed by new flag FSPP:                              *
*    FSPP=0 without FSR pairs;                                             *
*    FSPP=1  with   FSR pairs (additive treatment);                        *
*    FSPP=2  with   FSR pairs (multiplicative treatment).                  *
*                                                                          *
* 4. For the FSPP the cut on the invariant mass of the secondary pair      * 
*    is accessible.                                                        *
*    In order to accomodate this cut value, the variable SIPP of the       *
*    SUBROUTINE ZUCUTS(INDF,ICUT,ACOL,EMIN,S_PR,ANG0,ANG1,SIPP)            *
*    is used.                                                              *
*    Therefore, the meaning of the variable SIPP has been changed.         *
*    It has nothing to do with cutting of ISPP; there is no possibility to *
*    cut secondary pairs for ISPP, whre the primary pair cut should        *
*    be equal to S_PR                                                      *
*                                                                          *
* 5. ZFITTER v. 6_30 should be used together with DIZET v. 6_23            *
*    Changes in DIZET v. 6_30 compared to v. 6_21 are due to bug fixes:    *
*    1) A bug in calculation of running \alpha_em (D. Bardin)              *
*    2) A bug in calculation of \Gamma_W is fixed (L. Kalinovskaya)        *
*                                                                          *
*    An option to fit V_TB is available too.                               *
*    It might be accessed with DIZET v. 6_30 together with the use of      *
*    zwidthtb6_30.f (D. Bardin, L. Kalinovskaya, A. Olshevsky, March 2000) *
*    However, since for this case                                          *
*    the argument list of DIZET is changed to accomodate this possibility. *
*    On the other hand, the fit of V_tb is not used too widely and its use *
*    is limited by LEP1 energies. Given two above facts, DIZET v. 6_30     *
*    is not interfaced with ZFITTER v. 6_30. This might be done upon       *
*    an explicit request of LEP community.                                 *
*                                                                          *
* 6. A new option governed by the flag IPTO is implemented.                *
*    IPTO=0,1,2,3 old treatment (see below)                                *
*    IPTO=-1 new treatment                                                 *
*         allows to calculate the pure virtual pair contribution separately*
*==========================================================================*
*                                                                          *
* ZFITTER v. 6_23 (13-12-99):                                              *
* ---------------------------                                              *
* A bug of version 6_22 pointed out by Martin Gruenewald is fixed          *
*                                   by Andrej Arbuzov                      *
*==========================================================================*
*                                                                          *
* ZFITTER v. 6_22 (15-10-99):                                              *
* ---------------------------                                              *
* Second order corrections to angular distributions and A_FB:              *
* photonic RC are improved, O(alpha^2) LLA pairs are added                 *
* Flag FBHO = 0 (old version), FBHO = 1 (new treatment)                    *
* (to be described in an extended version of hep-ph/9907500)               *
*==========================================================================*
*                                                                          *
* ZFITTER v. 6_21 (13-10-99):                                              *
* ---------------------------                                              *
* Differs from the versions from 26-07-99 and 01-10-99 only by bug fixes:  *
* Wynhoff.bug, a bug affecting angular distribution with IFI=1, etc.       *
*==========================================================================*
*                                                                          *
* ZFITTER v. 6_21 (26-07-99):                                              *
* ---------------------------                                              *
*  Treatment of pairs is extended a lot                                    *
*  (see A.B. Arbuzov, hep-ph/9907500)                                      *
*                                                                          *
*  ISPP = 2 allows to include singlet channel and higher orders            *
*           according to A. Arbuzov                                        *
*  ISPP = 3 coded according to JMS paper                                   *
*  ISPP = 4 JMS + extended treatment of hadronic pairs                     *
*                                                                          *
*  IPFC = 1 only electron   pairs                                          *
*  IPFC = 2 only muon       pairs                                          *
*  IPFC = 3 only tau-lepton pairs                                          *
*  IPFC = 4 only hadron     pairs                                          *
*  IPFC = 5 all summed                                                     *
*  IPFC = 6 leptonic pairs (without hadrons)                               *
*                                                                          *
*  flag IPSC =0,1,2,3 governs the siglet pair contributions                *
*                     (see description in ZUFLAG).                         *
*                                                                          *
* Attention! Present default is IPSC=0                                     *
* The proper treatment of singlet pairs is presently under discussion!!!   *
*                                                                          *
*  flag IPTO =0,1,2,3 governs the higher order pair contributions          *
*                     (see description in ZUFLAG).                         *
*                                                                          *
*==========================================================================*
*                                                                          *
* ZFITTER v. 6_20 (30-06-99):                                              *
* ---------------------------                                              *
*                                                                          *
* Angular distributions of cross sections in c=\cos\theta (= CSA)          *
* Update of branch with acollinearity/minimal energies cut to the final    *
* state fermions: corrected hard initial state and initial-final state     *
* interference flux functions                                              *
*                                                                          *
* Modified subroutines/functions in v. 6_1:                                *
*                                                                          *
* ZANCUT <-- COSCUT <-- HARD                                               *
*                                                                          *
* New functions called from updated package acol6_1.f in subroutine HARD:  *
*                                                                          *
*   HCUTACOL and HCUTC (replace original code HCUT,HINIM,HINTFM)           *
*                                                                          *
*   SAVAR (update!), DLOGCP,DLOGCM      (logarithms,polynomials)           *
*   DRADTINI  <--- DRTFUN0,DRTFUN1      (sym. flux functions in c)         *
*   DRADFBINI <--- DRFBFUN0,DRFBFUN1    (anti-sym. flux functions in c)    *
*                                                                          *
* New definition of flag ICUT/IRCUT in analogy to calculation of total     *
* cross sections/asymmetries for the angular distributions:                *
*                                                                          *
*     IRCUT (int/read) =  0      ACOL+EMIN cuts         (original code)    *
*                      =  1      S_PR cut               (original code)    *
*                      =  2 or 3 ACOL+EMIN+cuts         (new code)         *
*     COSCUT uses IRCUT=  3                                                *
*                                                                          *
*==========================================================================*
*                                                                          *
* ZFITTER v. 6_11 (22-06-99):                                              *
* 1) Little bug fixes, indicated by J.Mhich with a fix in SETCUT as        *
*                                           recommended by M.Gruenewald    *
* 2) Polishing of FORTRAN while writing description, no change of numbers  *
* ---------------------------                                              *
*                                                                          *
* ZFITTER v. 6_10 (27-05-99):                                              *
* ---------------------------                                              *
*                                                                          *
* ATTENTION!!! ISPP at LEP2 with loose cuts, allowing radiative return, is *
*          still under investigation. These numbers may not be reliable!!! *
*          [Courtesy Michael Kobel and Pat Ward (OPAL)]                    *
*                                                                          *
* Use of BOXD=2 is forbidden, will be allowed from 6.20                    *
*                                                                          *
* MISD flag is implemented (Model Independent, S-Dependence)               *
* MISD=0 S=M^2_z in EWRC, old treatment                                    *
* MISD=1 default, new option,                                              *
*        this ensures equal numbers from all interfaces                    *
*        for all \sqrt{S} up to 100 GeV                                    *
*                                                                          *
*==========================================================================*
*                                                                          *
* ZFITTER v. 6_08 (19-05-99):                                              *
*                                                                          *
* ZUCUTS is replaced by Tatsuo's version, should have NO EFFECT on numbers *
*                                                                          *
* ZFITTER v. 6_07 (18-05-99):                                              *
*                                                                          *
* In ZUXSEC bug 2*u+3*d ---> 2*u+2*d+b  corrected                          *
* Some cosmetic bugs corrected, ZUINFO prints ICUT parameter as expected   *
*                                                                          *
* ZFITTER v. 6_06 (17-05-99):                                              *
*                                                                          *
* 0. BUG fixes: Two bugs indicated by Tatsuo Kamamoto 12/05/99 fixed       *
*               Bugs indicated by Guenter Quast 13/05/99 fixed,            *
*               because of that the date of release is changed to 13/05/99 *
*                                                                          *
* Main changes in comparison to older versions 6_04 and 6_03:              *
*                                                                          *
* 1. MISC-flag is implemented                                              *
*    MISC=0 (default) non-scaled \rho's are used, AROTFZ-array             *
*    MISC=1               scaled \rho's are used, ARROFZ-array             *
*                                                                          *
* 2. Options to simulate theoretical precision of Pseudo-Observables       *
*    are updated as compared to 5_20.                                      *
*    The options:                                                          *
*    SCAL(ISCAL)=0,4,4  |                                                  *
*    HIGS(IHIGS)=0,1,1  |                                                  *
*    SCRE(ISCRE)=0,2,1  \                                                  *
*    EXPR(IFACR)=0,2,1    -  (Please, do not use the other flag values!)   *
*    EXPF(IFACT)=0,2,1  /                                                  *
*    HIG2(IHIG2)=0,1,1  |                                                  *
*    produce following uncertainty intervals                               *
*    for M_z=91.1867, M_t=173.8, M_h=100, \als=0.119:                      *
*                                                                          *
*                                          5.12                            *
*                                min      central      max      interval   *
* 1) M_w                     = 80.3669    80.3738    80.3742    0.73  MeV  *
* 2) \sigma^{0}_{had}        = 41.4777    41.4777    41.4785    0.8   pb   *
* 4) \Gamma_{e}              = 83.983     83.995     83.999     0.016 MeV  *
*10) \Gamma_{b}              = 375.886    375.993    375.995    0.109 MeV  *
*13) \Gamma_{Z}              = 2.49527    2.49573    2.49578    0.51  MeV  *
*14) R_e                     = 20.7406    20.7420    20.7421    0.0015     *
*15) R_b                     = 0.215786   0.215811   0.215813   0.000027   *
*17) sin^2\theta^{lept}_{eff}= 0.231594   0.231601   0.231653   0.000059   *
*18) sin^2\theta^{b   }_{eff}= 0.232941   0.232950   0.233004   0.000063   *
*20) A^{0,l}_{FB}            = 0.015984   0.016074   0.016086   0.000102   *
*21) A^{0,b}_{FB}            = 0.102327   0.102617   0.102657   0.000330   *
*23) {\cal{A}}_{e}           = 0.145988   0.146396   0.146452   0.000464   *
*24) {\cal{A}}_{b}           = 0.934573   0.934607   0.934613   0.000040   *
*                                                                          *
* The complete Table of POs to appear in ZFITTER description               *
*                                                                          *
* The pseudo observables as defined in the PCP report are                  *
* stored in the following COMMON blocks                                    *
*                                                                          *
* effective sinuses      : array ARSEFZ in COMMON /CDZRKZ/                 *
* rho_f                  : array AROTFZ in COMMON /CDZRKZ/                 *
* partial Z decay widths : array WIDTHS in COMMON /ZUPARS/                 *
*                                                                          *
* 3. Simplification of expressions in acol6_05.f:                          *
*   much shorter expressions for O(alpha) hard interference                *
*   corrections for branch with acollinearity and acceptance cut (ICUT=3); *
*                                                                          *
*   former subr. HARDINTACOL in acol6_05.f removed and hard intf.          *
*   expressions with cuts now treated parallel to initial-state terms      *
*   in subr. RADTIN and RADFBIN and functions called there;                *
*   hard final-state corrections with cuts (not used yet) treated          *
*   separately in subr. RADTFIN and RADFBFIN                               *
*                                                                          *
* - Removing of some Fortran bugs in acoll. branch (2 undefined logarithms *
*   [1 bug for ICUT=2 had already been present in acol6_04.f],             *
*   one missing exponent in acol6_05.f)                                    *
*                                                                          *
* - Just some optical 'improvements' (change of the order of the           *
*   different subroutines in acol6_05.f; subr. which are not used          *
*   at the end of acol6_05.f; shorter names for some subr. in acol6_05.f,  *
*   restructuring some (not used) parts in SHARD, AHARD in zfitr6_05.f)    *
*                                                                          *
*==========================================================================*
*                                                                          *
* ZFITTER v. 6_04 (19-04-99):                                              *
* ZFITTER v. 6_03 (01-04-99):                                              *
* This version replaces 5.xx series!                                       *
*                                                                          *
* The following changes are activated in zfitr6_03.f:                      *
* SUBROUTINE ZUCUTS                                                        *
* - a print added                                                          *
* SUBROUTINE SCUT                                                          *
* - comment added on ICUT                                                  *
* FUNCTIONs SHARD, AHARD                                                   *
* - new chains with IRCUT=2,3 activated; calls to package acol6_03.f       *
* FUNCTION FAL2(Z)                                                         *
* - debugged; see there for details                                        *
* SUBROUTINE FUNFIN(S,AMF2,QF2,R,SFIN,AFIN)                                *
* - debugged; see there for details                                        *
* further comments see inside package acol6_03.f                           *
*                                                                          *
* List of some important changes for LEP1 use (compared to ZFITTER_5_21):  *
*                                                                          *
* Bug fixes                                                                *
* ---------                                                                *
* Several bugs affected calculation of RO's with ALEM=2                    *
* (pointed out by Martin Gruenewald) are fixed.                            *
* In particular a bug affecting A_fb for light (u,d,s) quarks              *
*==========================================================================*
* ZFITTER v. 5.21                                                          *
* An intermediate version between 5.xx series intended for LEP1 use and    *
* fortcoming version 6.00 suitable both at LEP1 and LEP2 energies.         *
*                                                                          *
* List of some important changes for LEP1 use (compared to ZFITTER_5_20):  *
*                                                                          *
* 0. Bug in IFI for hadrons (color factor double-counting) is fixed.       *
*    ----------------------                                                *
*                                                                          *
* 1. New treatment of ZMIN-cut                                             *
*    -------------------------                                             *
*    BEWARE!!! ZUCUTS argument list is modified!!!                         *
*                                                                          *
* Old: SUBROUTINE ZUCUTS(INDF,ICUT,ACOL,EMIN,S_PR,ANG0,ANG1)               *
* New: SUBROUTINE ZUCUTS(INDF,ICUT,ACOL,EMIN,S_PR,ANG0,ANG1,SIPP)          *
*                                                           SIPP=ZMIN*S    *
* Last argument (SIPP) provides user's access to the parameter ZMIN which  *
*                      governs cut on the secondary pairs in ISPP and FSPP *
*                                                                          *
* 2. New FLAG FSRS is added (Final State Radiation Scale)                  *
*    ----------------------------------------------------                  *
*    FSRS, FSRS=1 default                                                  *
*    ----                                                                  *
*    if FSRS= 0, FSR corrections are calculated with the scale \alpha(0)   *
*    if FSRS= 1, FSR correstions are calculated with the scale \alpha(S)   *
*                                                                          *
* 3. Use at LEP2 energies.                                                 *
*    ---------------------                                                 *
* For use at LEP2 energies a special attention should be paid on some      *
* flags:                                                                   *
*                                                                          *
* 1) FINR should be reset to 0, since such a value seems to be more        *
*    adequate to LEP2 analysis.                                            *
* 2) Use default ISPP=1, but be careful with setting of SIPP variable in   *
*    the SUBROUTINE ZUCUTS. It looks very reasonable to use SIPP=S_PR,     *
*    as was suggested by Michael Kobel (OPAL).                             *
*    It should be understood, however, that authors of ZFITTER do not take *
*    any responsibility for the concrete choice of flags and parameters.   *
*    This should be responsibility of the users.                           *
*                                                                          *
*==========================================================================*
*                                                                          *
* ZFITTER v. 5.20                                                          *
* List of some important changes for LEP1 use (compared to ZFITTER_5_15):  *
*                                                                          *
* Look inside ZUFLAG and DIZET for updated descriptions of flags.          *
* ---------------------------------------------------------------          *
*                                                                          *
* 1. Meaning of s' - cut was critically revisited.                         *
*    ---------------------------------------------                         *
*    It is governed now by flag FINR with the following meaning:           *
*                                                                          *
*    FINR, FINR=1 default                                                  *
*    ----                                                                  *
*    if FINR=-1, NO FSR and s'=M^2, where the M^2 is the final state       *
*                                   fermions invariant mass                *
*    if FINR= 0, FSR is treated inclusively, however consistently with     *
*                    imposing an s'-cut which in this case means           *
*                    the invariant mass of the virtual particle (Z or G)   *
*    if FINR= 1, FSR is taken into account and the M^2 cut on the final    *
*                    state fermions invariant mass may be imposed          *
*                                                                          *
* Note: only in latter case the Initial-Final Interference is meaningful.  *
*                                                                          *
* IMPORTANT:                                                               *
* If FINR=0,1 then SBRN, ABRN (SIGBRN,AFBBRN) include FSR corrections.     *
* In other words:  SQED, AQED (SIGQED,AFBQED) differ from corresponding    *
* BORN values only in the presence of ISR QED radiation                    *
*                                                                          *
* 2. SFAST interface is enabled to calculate realistic observables         *
*    ---------------    with BOTH s' and M^2 cuts and with all other       *
*                       options but ACOL/EMIN cuts.                        *
*    As before, SFASTis accessible at ICUT=-1.                             *
*                                                                          *
* 3. The meaning of CONV flag is revised.                                  *
*    ------------------------------------                                  *
*    CONV, CONV=1 default                                                  *
*    ----                                                                  *
*    if CONV=0, no convolution of running couplings as before              *
*    if CONV=1, \alpha_{em}(s') and convoluted                             *
*    if CONV=2, everything, including EWRC are calculated at s'            *
*               and convoluted (the latter is very time consuming)         *
*                                                                          *
* IMPORTANT:                                                               *
* It was checked that for LEP1 use it is absolutely not needed to run      *
* with CONV=2, for instance at M_Z+3 its effect is 1.5e-4.                 *
*                                                                          *
* 4. FTJR and CZAK flags are revised.                                      *
*    --------------------------------                                      *
*                                                                          *
*    FTJR, FTJR=1 default.                                                 *
*    CZAK, CZAK=1 default.                                                 *
*    ----                                                                  *
*    if FTJR/CZAK=0 corrections OFF                                        *
*    if FTJR/CZAK=1 corrections ON                                         *
*    if FTJR/CZAK=2 corrections ON for PO's and OFF for RO's               *
*                   (this is for comparison with TOPAZ0 only               *
*                    and shouldn't be used in "users' runs")               *
*    Nothing is changed in the phylosophy of using CZAK flag:              *
*                   OFF for MI fits                                        *
*                   ON  for SM fits                                        *
*    However, by now it is open what to do with FTJR correction.           *
*    Physically, it is of the same origin: mixed O(\alpha\als) correction. *
*    Technically, it is implemented differently. A part of FTJR correction *
*                 is 'factorized' with FSR QCD correction and is added to  *
*                 the EW form factors of e^+e^- --> b\bar{b} amplitude.    *
*                                                                          *
* 5. FOT2 flag is COMPLETELY revised.                                      *
*    --------------------------------                                      *
*                                                                          *
*    FOT2, FOT2=3 default (recommended)                                    *
*    ----                                                                  *
*    if FOT2=-1 then NO ISR QED radiation at all                           *
*    if FOT2=0  then the complete additive O(\alpha) radiator              *
*            1  then the complete additive O(\alpha^2) radiator            *
*               with 'some' LLA corrections                                *
*               (supposed to be critically revisited)                      *
*    if FOT2=2  then the complete additive O(\alpha^2) radiator            *
*    if FOT2=3  then the complete additive O(\alpha^2) + LLA O(\alpha^3)   *
*    if FOT2=4  then the QED-E radiator, Eq. (3.31)-(3.32) of YR (1989)    *
*    by F.A. Berends and W.L. van Neerven (unpublished) is implemented     *
*    if FOT2=5 then LLA 'pragmatic' factorized radiator, a'la M. Skrzypek  *
*    IMPORTANT!!!                                                          *
*       FOT2=5 is compatible only with s' and M^2 cuts.                    *
*    IMPORTANT!!! As seen, FOT2 -flag is completely decoupled with initial *
*                 state pair production corrction, which is governed       *
*                 by flag ISPP with the following meaning:                 *
*                                                                          *
* 6. ISPP - Initial State Pair Production.                                 *
*    -------------------------------------                                 *
*                                                                          *
*    ISPP, ISPP=1 default (recommended)                                    *
*    ----                                                                  *
*    if ISPP=-1 then the treatment of ISPP as in all old ZFITTERs          *
*    if ISPP=0  no initial state pair corrections                          *
*    if ISPP=1  initial state pair corrections are implemented as in KKKS  *
*               without any artificial adjusting parameters as it was in   *
*               older versions                                             *
*                                                                          *
* 7. INCL flag is eliminated as suggested by Olshevsky.                    *
*    -----------------------                                               *
*    In its place a new flag appears.                                      *
*                                                                          *
* 8. New flag DIAG.                                                        *
*    --------------                                                        *
*                                                                          *
*    DIAG, DIAG=1 default                                                  *
*    ----                                                                  *
*    if DIAG=-1 then only Z-exchange diagram is take into account          *
*    if DIAG= 0 then G and Z squares only, i.e neglecting GZ-interference  *
*    if DIAG= 1 complete treatment - | G + Z |^2                           *
*                                                                          *
*==========================================================================*
*                                                                          *
* ZFITTER v. 5.15                                                          *
* bug fixes (courtesy Tatsuo Kawamoto)                                     *
*                                                                          *
*==========================================================================*
* ZFITTER v. 5.14/5.13                                                     *
* List of some important changes for LEP1 use (compared to ZFITTER_5_12):  *
*--------------------------------------------------------------------------*
*                                                                          *
* 0. Lepton masses are updated as in PDG'98                                *
*------------------------------------------                                *
*                                                                          *
* 1. Running \alpha_{em}                                                   *
*-----------------------                                                   *
*                                                                          *
* In version 5.13 we realised the idea to release for fit                  *
* the hadronic vacuum polarization for five flavours                       *
*         \Delta\alpha^{(5)}_{had}(M_Z)=DAL5H                              *
* instead of                                                               *
*                      \alpha_{em}(M_Z)=ALQED                              *
* is it was in ALL earlier versions.                                       *
* In all user interfaces, ZU****, the number and sequence of arguments     *
* is not changed, however, instead of ALQED one should use now DAL5H.      *
* Furthermore, recent calculations of M. Steinhauser (hep-ph/9803313) of   *
* the three-loop corrections to the leptonic contribution to \Delta\alpha  *
* are implemented.                                                         *
* Because of these, the meaning of ALE2 flag has been revised:             *
*                                                                          *
*    ALE2, ALE2=3 default                                                  *
*    ----                                                                  *
*       ALE2=0 retained for back compatibility with version 5.12,          *
*              its meaning being the same as described in 5.12.            *
*   For ALE2>0:                                                            *
*              The leptonic \Delta\alpha is calculated with:               *
*    if ALE2=1 one  -loop \                                                *
*       ALE2=2 two  -loop  - corrections                                   *
*       ALE2=3 three-loop /                                                *
*                                                                          *
* Since the three-loop corrections became available, ALE2 - flag should not*
* be treated as a "theoretical option" anymore and it is recommended to use*
* its default setting. Using the other values is recommended only in order *
* to see the influence of higher orders.                                   *
*                                                                          *
* The influence of the three-loop corrections is very tiny.                *
* For the vast majority of numbers printed by ZFTEST it is even unseen     *
* within digits printed.                                                   *
*                                                                          *
* 2. Fermi constant                                                        *
*------------------                                                        *
*                                                                          *
* Recent results of T. van Ritbergen and R. Stuart  (hep-ph/9808283) on    *
* the two-loop QED corrections for \mu-decay are taken into account with   *
* the aid of new flag:                                                     *
*                                                                          *
*    GFER, GFER=2 default                                                  *
*    ----                                                                  *
*       GFER=0 retained for back compatibility with version 5.12.          *
*   For GFER>0:                                                            *
*              The Fermi constant corresponds to QED corrections           *
*              taken into account in the:                                  *
*       GFER=1 one-loop \                                                  *
*                        - approximation                                   *
*       GFER=2 two-loop /                                                  *
*                                                                          *
* The influence of the two-loop corrections is tiny,                       *
* although seen within digits printed by ZFTEST.                           *
*                                                                          *
* 3. \rho_f                                                                *
*----------                                                                *
*                                                                          *
* Flavour dependent \rho^f - factors for partial decay widths              *
* are printed together with effective sine^2\Theta^f_{eff}.                *
*                                                                          *
*--------------------------------------------------------------------------*
* ZFITTER v. 5.12                                                          *
* List of some important changes for LEP1 use (compared to ZFITTER_5_11):  *
*--------------------------------------------------------------------------*
*                                                                          *
* 0. Bug fixes for CZAK=2, pointed by Martin.                              *
*    Bug fixes IF(AMT4...) --> IF(IAMT4...), pointed by Guenter.           *
*                                                                          *
* 1. One flag is extended:                                                 *
*                                                                          *
*    FOT2, FOT2=4 default (recommended)                                    *
*    ----                                                                  *
*    if FOT2=5 then the QED-E radiator, Eq. (3.31)-(3.32) of YR (1989)     *
*    by F.A. Berends and W.L. van Neerven (unpublished) is implemented.    *
*                                                                          *
* 2. New flags:                                                            *
*                                                                          *
*    HIG2, HIG2=1 default                                                  *
*    ----                                                                  *
*    if HIG2=0 two loop Higgs correction for \Delta r is not included      *
*    if HIG2=1                                        is included          *
*                                                                          *
*    ALE2, ALE2=1 default                                                  *
*    ----                                                                  *
*    if ALE2=0 two loop constant correction for \Delta\alpha               *
*                                                     is not included      *
*    if ALE2=1                                        is included          *
*                                                                          *
*    Comment: Both effects were de-activated in 5.10 by convention.        *
*             In 5.12 they are reactivated. The effect of ALE2=0/1         *
*             is totally negligible. The effect of HIG2=0/1 is totally     *
*             negligible for M_h < 300. It grows then with growing of M_h. *
*             At M_h=1000 it results in -2 MeV for M_w.                    *
*                                                                          *
* 3. Options to simulate theoretical precision of Pseudo-Observables.      *
*    All options of version 4.9 used in 1994 LEP1 Workshop on              *
*    "Precision calculation for Z-resonance" are made compatible with      *
*    AMT4=4 which gives, in principle, an opportunity to produce           *
*    the "blue band". However, we are yet discussing about options.        *
*                                                                          *
*    The options:                                                          *
*    SCAL(ISCAL)=0,4,4  |                                                  *
*    HIGS(IHIGS)=0,1,1  |                                                  *
*    SCRE(ISCRE)=0,2,1  \                                                  *
*    EXPR(IFACR)=0,2,1    -  (Please, do not use the other flag values!)   *
*    EXPF(IFACT)=0,2,1  /                                                  *
*    HIG2(IHIG2)=0,1,1  |                                                  *
*    ALE2(IALE2)=0,1,1  |                                                  *
*    produce following uncertainty intervals                               *
*    for M_z=91.1888, M_t=175, M_h=300, \als=.125 (to compare with 95-03)  *
*                                                                          *
*                                      5.12                   |    4.9     *
*                                                             | (YR 1994)  *
*    M_w                     = 80.2998 -:- 80.3089   9.1  MeV |  14  MeV   *
*    sin^2\theta^{lept}_{eff}= 0.232090-:- 0.232180  0.000090 |  0.00020   *
*    \Gamma_{Z}              = 2.49567 -:- 2.49631   0.64 MeV |  1.1 MeV   *
*    R_e                     = 20.7708 -:- 20.7729   0.0021   |  0.0070    *
*                                                                          *
*    and much more narrow intervals for M_h=100 GeV which are shown now    *
*    for M_z=91.1867, M_t=175.6, M_h=100, \als=0.120                       *
*                                                                          *
*                                          5.12                            *
*                                min      central      max      interval   *
*    M_w                     = 80.3786    80.3860    80.3862    1.6  MeV   *
*    sin^2\theta^{lept}_{eff}= 0.231509   0.231516   0.231577   0.000068   *
*    \Gamma_{Z}              = 2.49633    2.49681    2.49687    0.54 MeV   *
*    R_e                     = 20.7471    20.7486    20.7486    0.0015     *
*                                                                          *
*--------------------------------------------------------------------------*
* ZFITTER v. 5.11                                                          *
* List of some important changes for LEP1 use (compared to ZFITTER_5_10):  *
*--------------------------------------------------------------------------*
*                                                                          *
* 1. Several tiny problems pointed by users of 5.10 are fixed              *
*    (none with any influence on numbers)                                  *
*                                                                          *
* 2. One more flag is extended:                                            *
*                                                                          *
*    FOT2, FOT2=4 default                                                  *
*    ----                                                                  *
*    if FOT2=4 then the \alpha^3 corrections                               *
*    by G. Montagna, O. Nicrosini and F. Piccinini, hep-ph/9611463,        *
*    are implemented.                                                      *
*                                                                          *
* 3. The implementation of Czarnecki/Kuehn corrections, CZAKUE, and        *
*    meaning of the flag CZAK is revised,                                  *
*                                                                          *
*    CZAK, CZAK=1 default (recommended)                                    *
*    ----                                                                  *
*    if CZAK=0 then CZAKUE corrections are not implemented;                *
*    if CZAK=1 then CZAKUE corrections are implemented                     *
*              as suggested by db, i.e. with additional suppression        *
*              by BW and therefore with a convolution. Moreover,           *
*              as suggested by Martin Grunewald, they are done             *
*              channel-dependent;                                          *
*    if CZAK=2 then Czarnecki/Kuehn corrections are implemented            *
*              as suggested by Martin Grunewald, i.e without additional    *
*              suppression by BW as constant (no convolution),             *
*              channel-dependent factors.                                  *
*                                                                          *
* 4. New flag:                                                             *
*                                                                          *
*    PREC, PREC=1 default                                                  *
*    ----                                                                  *
*    1 < PREC < 99 is an integer number which ANY precision govering ANY   *
*    numerical integration is devided to, increasing thereby the numerical *
*    precision of computation. In some cases when some numerical instability
*    while running 5_10 was registered, it was sufficient to use PREC=3    *
*    in some other cases (Ptau) only PREC=30 solved the instability.       *
*                                                                          *
* 5. The options to simulate theoretical uncertainties are still           *
*    not compatible with AMT4=4. Hopefully, a set of plausible options     *
*    for AMT4=4 will be developped during next week: 11/06-18/06 1998.     *
*                                                                          *
* 6. The interval of allowed \sqrt(s) is:                                  *
*        9.5 < \sqrt(s) < 350 GeV,                                         *
*    i.e. the code is not intended to run                                  *
*    below bb and above the tt-thresholds                                  *
*                                                                          *
*--------------------------------------------------------------------------*
* ZFITTER v. 5.10                                                          *
* List of some important changes for LEP1 use (compared to ZFITTER_beta):  *
*--------------------------------------------------------------------------*
*                                                                          *
* 1. Two flags are extended:                                               *
*                                                                          *
*    AMT4, AMT4=4 default                                                  *
*    ----                                                                  *
*    if AMT4=4 then the two loop electroweak corrections                   *
*              by Gambino and Degrassi are ON                              *
*                                                                          *
*    QCDC, QCDC=3 default                                                  *
*    ----                                                                  *
*    if QCDC=3 then mixed corrections O(\alpha_em\alpha_st)                *
*              are calculated with B. Kniehl library                       *
*                                                                          *
* 2. New flag:                                                             *
*                                                                          *
*    CZAK, CZAK=1 default                                                  *
*    ----                                                                  *
*    if CZAK=0 then Czarnecki/Kuehn corrections are not implemented        *
*    if CZAK=1 then Czarnecki/Kuehn corrections are implemented            *
*                                                                          *
* 3. The options to simulate theoretical uncertainties are presently       *
*    not compatible with AMT4=4. Hopefully, a set of plausible options     *
*    for AMT4=4 will be developed after Moriond'98 well before Canadian    *
*                                                             Rochester    *
*                                                                          *
* 4. The interval of allowed \sqrt(s) is presently limited by              *
*    9 < \sqrt(s) < 151 GeV. This has to be improved soon.                 *
*                                                                          *
*==========================================================================*
*                                                                          *
      SUBROUTINE ZFTEST(IMISC)                                             
*     ========== =============                                             *
****************************************************************************
*                                                                          *
*     SUBR. ZFTEST                                                         *
*                                                                          *
*     Example program to demonstrate the use of the ZFITTER package.       *
*                                                                          *
****************************************************************************
*
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 XALLCH,XFOTF
*
      DIMENSION XS(0:11,6),AFB(0:11,5),TAUPOL(2),TAUAFB(2),ALRI(0:11,2)
*
* constants
*
      PARAMETER(ALFAI=137.0359895D0,ALFA=1.D0/ALFAI,CONS=1.D0)
      PARAMETER(ZMASS=91.1867D0,TMASS=173.80D0,HMASS=120.D0)
      PARAMETER(ALFAS=.1184D0)
*     PARAMETER(ZMASS=91.1867D0,TMASS=173.8D0,HMASS=100.D0)
*     PARAMETER(ALFAS=.119D0)
      PARAMETER(RSMN=87.D0,DRS=1.D0,NRS=9)
      PARAMETER(ANG0=35D0,ANG1=145D0)
      PARAMETER(QE=-1.D0,AE=-.5D0,QU= 2.D0/3.D0,AU= .5D0,
     +                            QD=-1.D0/3.D0,AD=-.5D0)
*
* ZFITTER common blocks
*
      COMMON /ZUPARS/QDF,QCDCOR(0:14),ALPHST,SIN2TW,S2TEFF(0:11),
     & WIDTHS(0:11)
*
      COMMON /CDZRKZ/ARROFZ(0:10),ARKAFZ(0:10),ARVEFZ(0:10),ARSEFZ(0:10)
     &              ,AROTFZ(0:10),AIROFZ(0:10),AIKAFZ(0:10),AIVEFZ(0:10)
      COMMON /CDZAUX/PARTZA(0:10),PARTZI(0:10),RENFAC(0:10),SRENFC(0:10)
*
      COMMON /EWFORM/ XALLCH(5,4),XFOTF
*
*-----------------------------------------------------------------------
*
* initialize
*
      CALL ZUINIT(1)
*
* set ZFITTER flags and print flag values
*
      CALL ZUFLAG('PRNT',1)
      CALL ZUFLAG('MISC',0)
      CALL ZUFLAG('MISD',1)
      CALL ZUFLAG('ISPP',2)
      CALL ZUFLAG('AMT4',6)
*      CALL ZUFLAG('ALEM',2)
      CALL ZUFLAG('ASCR',1)
      CALL ZUFLAG('DMWW',0)
      CALL ZUFLAG('DSWW',0)
*
      ICUTC=-1
* 
* do weak sector calculations
*
      V_TB =1D0
*      DAL5H=2.7572D-02 
c       DAL5H=2.750D-02 
       DAL5H=2.7572D-02 
*     DAL5H=2.758D-02!+0.00035d0
      CALL ZVWEAK(ZMASS,TMASS,HMASS,DAL5H,V_TB,ALFAS)
c      CALL ZUFLAG('PRNT',0)
*
* Test CALL ZU_APV 
*
      UMASS=.1D0
      DMASS=.1D0
      CALL ZU_APV(ZMASS,TMASS,HMASS,SIN2TW,UMASS,DMASS,C1U,C1D,C2U,C2D)
*
      PRINT *       
      PRINT 1,C1U,C1D,C2U,C2D
 1    FORMAT(1X,'C1U(D)2U(D)=',4(F12.10,2X))
      PRINT *       
*
* define cuts for fermion channels and print cut values
*
c      CALL ZUCUTS(11,0,15.D0,10.D0,0.D0,ANG0,ANG1,0.D0)
      CALL ZUINFO(0)
c      print*,s2teff
*
* make table of cross sections and asymmetries
*
      PI   = DACOS(-1.D0)
*
* DO loop over S
*
      DO I = 6,6
*
* array of RS=SQRT(S)
*
        IF(I.EQ.1) RS=35D0
        IF(I.EQ.2) RS=65D0
        IF(I.EQ.3) RS=ZMASS-2D0
        IF(I.EQ.4) RS=ZMASS
        IF(I.EQ.5) RS=ZMASS+2D0
        IF(I.EQ.6) RS=100D0
        IF(I.EQ.7) RS=140D0
        IF(I.EQ.8) RS=175D0
*
* Changing of default, for instance:
* 1) NO PAIRS at PETRA and TRISTAN  
*
        IF(I.LT.3.) CALL ZUFLAG('ISPP',0)
        IF(I.GE.3.) CALL ZUFLAG('ISPP',1)
*
* The OUTPUT is adjusted for hp WS at Zeuthen
*
      IF(I.EQ.5.OR.I.EQ.8) THEN
        PRINT *
        PRINT *
        PRINT *
          ELSE
        PRINT *
      ENDIF
* table header
*        PRINT *,' SQRT(S) = ',REAL(RS)
*        PRINT *
*        PRINT 1000
 1000   FORMAT(1X,'     '
     +  ,'<-------------- Cross Section -------------->'
     +  ,'  <------- Asymmetry ------->','  <--Tau_Pol-->',
     +   '  <----A_LR--->')
*        PRINT 1001
 1001   FORMAT(1X,'INDF   ZUTHSM   ZUXSEC    ZUXSA   ZUXSA2   ZUXAFB'
     +,'   ZUTHSM  ZUXSA ZUXSA2 ZUXAFB ZUTPSM  ZUTAU  ZULRSM  ZUALR')
*
* loop over fermion indicies
*
*ls      DO INDF = 0,11
        DO INDF = 4,4
         S=RS**2
         IF(INDF.NE.11) 
     +   CALL ZUCUTS(INDF,ICUTC,0D0,0D0,1D-2*S,0D0,180D0,.25D0*S)
* standard model interf. (INTRF=1)
         CALL ZVTHSM(INDF,RS,ZMASS,TMASS,HMASS,DAL5H,V_TB,ALFAS,
     +     XS(INDF,1),AFB(INDF,1))
         IF(INDF.EQ.3) CALL ZVTPSM(RS,ZMASS,TMASS,HMASS,DAL5H,V_TB,
     +     ALFAS,TAUPOL(1),TAUAFB(1))
* cross section interf. (INTRF=2)
           GAME = WIDTHS(   1)/1000.
           GAMF = GAME          
           GAMZ = WIDTHS(11)/1000.
           IF(INDF.NE.11) GAMF = WIDTHS(INDF)/1000.          
           CALL ZUXSEC(INDF,RS,ZMASS,GAMZ,GAME,GAMF,XS(INDF,2))

* cross section & forward-backward asymmetry interf. (INTRF=3)
         IF(INDF.NE.0 .AND. INDF.NE.10) THEN
           IF(IMISC.EQ.0) THEN
             ROEE= AROTFZ(1)
           ELSE
             ROEE= ARROFZ(1)
           ENDIF
          IF(IMISC.EQ.0) THEN
            GAE= SQRT(AROTFZ(1))/2
          ELSE
            GAE= SQRT(ARROFZ(1))/2
          ENDIF
          GVE  = ARVEFZ(1)*GAE
*
           IF(INDF.LT.11) THEN
             IF(IMISC.EQ.0) THEN
               ROFI= AROTFZ(INDF)
               GAF = SQRT(AROTFZ(INDF))/2
             ELSE
               ROFI= ARROFZ(INDF)
               GAF = SQRT(ARROFZ(INDF))/2
             ENDIF
             GVF = ARVEFZ(INDF)*GAF
           ELSEIF(INDF.EQ.11) THEN
c             IF(IMISC.EQ.0) THEN
c               GAF = SQRT(AROTFZ(1))/2
c             ELSE
c               GAF = SQRT(ARROFZ(1))/2
c             ENDIF
             GAF = GAE
             GVF = ARVEFZ(1)*GAF
           ENDIF
           CALL ZUXSA(INDF,RS,ZMASS,GAMZ,0,GVE,GAE,GVF,GAF,
CB: MODE=1 CALL ZUXSA(INDF,RS,ZMASS,GAMZ,1,GVE/GAE,ROEE,GVF/GAF,ROFI,
     +     XS(INDF,3),AFB(INDF,3))
         ENDIF
c           if(indf.eq.1) print*,'Ae,gae,gva', 
c     &                           2*GVE*GAE/(GVE**2+GAE**2),gae,gve
c           if(indf.eq.6) print*,'Ac,gac,gvc', 
c     &                           2*GVF*GAF/(GVF**2+GAF**2),gaf,gvf
c           if(indf.eq.9) print*,'Ab,gab,gvb', 
c     &                           2*GVF*GAF/(GVF**2+GAF**2),gaf,gvf
          
* tau polarization interf. (INTRF=3)
         IF(INDF.EQ.3) CALL ZUTAU(RS,ZMASS,GAMZ,0,GVE,GAE,GVF,GAF,
     +     TAUPOL(2),TAUAFB(2))
* left-right polarization Asymmetry (IBRA=4)
         IF((INDF.GE.4.AND.INDF.LE.7).OR.INDF.EQ.9) THEN
          POL=.63D0
          CALL
     &    ZVLRSM(INDF,RS,ZMASS,TMASS,HMASS,DAL5H,V_TB,ALFAS,POL,
     &           XSPL,XSMI)
          ALRI(INDF,1)=(XSMI-XSPL)/(XSMI+XSPL)/POL
         ENDIF
         IF(INDF.EQ.10) THEN
          POL=.65D0
          XSMINS=0D0
          XSPLUS=0D0
          DO 2004 IALR=4,9
           IF(IALR.EQ.8) GO TO 2004
           CALL
     &     ZVLRSM(IALR,RS,ZMASS,TMASS,HMASS,DAL5H,V_TB,ALFAS,POL,
     &            XSPL,XSMI)
           XSMINS=XSMINS+XSMI
           XSPLUS=XSPLUS+XSPL
2004      CONTINUE
          ALRI(INDF,1)=(XSMINS-XSPLUS)/(XSMINS+XSPLUS)/POL
         ENDIF
*
* cross section & forward-backward asymmetry interf. for gv**2 and ga**2(IBRA=4)
*
         IF(INDF.GE.1 .AND. INDF.LE.3 .OR. INDF.EQ.11) THEN
           IF(INDF.LT.11) THEN
             IF(IMISC.EQ.0) THEN
               ROFI= AROTFZ(INDF)
               GAF = SQRT(AROTFZ(INDF))/2
             ELSE
               ROFI= ARROFZ(INDF)
               GAF = SQRT(ARROFZ(INDF))/2
             ENDIF
             GVF = ARVEFZ(INDF)*GAF
           ELSEIF(INDF.EQ.11) THEN
             IF(IMISC.EQ.0) THEN
               ROFI= AROTFZ(1)
               GAF = SQRT(AROTFZ(1))/2
             ELSE
               ROFI= ARROFZ(1)
               GAF = SQRT(ARROFZ(1))/2
             ENDIF
             GVF = ARVEFZ(1)*GAF
           ENDIF
           GVF2= GVF**2
           GAF2= GAF**2
           CALL ZUXSA2(INDF,RS,ZMASS,GAMZ,0,GVF2,GAF2,
CB: MODE=1 CALL ZUXSA2(INDF,RS,ZMASS,GAMZ,1,GVF2/GAF2,ROEE**2,
     +     XS(INDF,4),AFB(INDF,4))
*
* cross section & forward-backward asymmetry interf. for gve*gae*gvf*gaf, 
* gve**2+gae**2 and gvf**2+gaf**2 (IBRA=6)
*
           PFOUR=GVE*GVF*GAE*GAF
           PVAE2=GVE**2+GAE**2
           PVAF2=GVF**2+GAF**2
           IF(INDF.NE.11) THEN
            CALL ZUXAFB(INDF,RS,ZMASS,GAMZ,PFOUR,PVAE2,PVAF2,
     +      XS(INDF,5),AFB(INDF,5))
           ENDIF
         ENDIF
* S-matrix interf. (ISMA=1, via SMATASY)
* results
         IF(INDF.EQ.0) THEN
*           PRINT 9000,INDF,(XS(INDF,J),J=1,2)
         ELSEIF(INDF.EQ.1) THEN
*           PRINT 9010,INDF,(XS(INDF,J),J=1,5),AFB(INDF,1),
*     +      (AFB(INDF,J),J=3,5)
         ELSEIF(INDF.EQ.2) THEN
*           PRINT 9010,INDF,(XS(INDF,J),J=1,5),AFB(INDF,1),
*     +      (AFB(INDF,J),J=3,5)
         ELSEIF(INDF.EQ.3) THEN
*           PRINT 9015,INDF,(XS(INDF,J),J=1,5),AFB(INDF,1),
*     +      (AFB(INDF,J),J=3,5),(TAUPOL(J),J=1,2)
         ELSEIF(INDF.EQ.10) THEN
*           PRINT 9025,INDF,(XS(INDF,J),J=1,2),(ALRI(INDF,J),J=1,2)
         ELSEIF(INDF.EQ.11) THEN
*           PRINT 9011,INDF,(XS(INDF,J),J=1,4),AFB(INDF,1),
*     +      (AFB(INDF,J),J=3,4)
         ELSE
*           PRINT 9020,INDF,(XS(INDF,J),J=1,3),AFB(INDF,1),AFB(INDF,3)
*     +                  ,(ALRI(INDF,J),J=1,2)
          ENDIF
        ENDDO
        PRINT *
      ENDDO
      RETURN
 9000 FORMAT(1X,I4, 2F9.5)
 9010 FORMAT(1X,I4, 5F9.5, 2X,4F7.4)
 9011 FORMAT(1X,I4, 4F9.5,11X,4F7.4)
 9015 FORMAT(1X,I4, 5F9.5, 2X,6F7.4)
 9020 FORMAT(1X,I4, 3F9.5,20X,2F7.4,28X,2F7.4)
 9025 FORMAT(1X,I4, 2F9.5,71X,2F7.4)
*                                                             END ZFTEST
      END
 
      SUBROUTINE ZUINIT(IPRINT)
*     ========== ==============
************************************************************************
*
*     SUBR. ZUINIT
*
*     Initializes ZFITTER common blocks
*
************************************************************************
*
      IMPLICIT REAL*8 (A-H,O-Z)
*
* program version
*
      COMMON /PROVER/ IZFVER,IZFDAT
*
* flags
*
      PARAMETER(NFLGMX=46)
      COMMON /ZUFLGS/ IFLAGS(NFLGMX),CFLAGS(NFLGMX)
      CHARACTER CFLAGS*4
      PARAMETER (IFAFBC= 1,IFSCAL= 2,IFSCRE= 3,IFAMT4= 4,IFBORN= 5,
     & IFBOXD= 6,IFCONV= 7,IFFINR= 8,IFFOT2= 9,IFGAMS=10,IFDIAG=11,
     & IFINTF=12,IFBARB=13,IFPART=14,IFPOWR=15,IFPRNT=16,IFALEM=17,
     & IFQCDC=18,IFVPOL=19,IFWEAK=20,IFFTJR=21,IFEXPR=22,IFEXPF=23,
     & IFHIGS=24,IFAFMT=25,IFCZAK=26,IFPREC=27,IFHIG2=28,IFALE2=29,
     & IFGFER=30,IFISPP=31,IFFSRS=32,IFMISC=33,IFMISD=34,IFIPFC=35,
     & IFIPSC=36,IFIPTO=37,IFFBHO=38,IFFSPP=39,IFFUNA=40,IFASCR=41,
     & IFSFSR=42,IFENUE=43,IFTUPV=44,IFDMWW=45,IFDSWW=46)
*
      COMMON/PARTZW/PARTZ(0:11),PARTW(3)
      COMMON /EXPERT/ IMASK,IMASS
*
* cuts
*
      COMMON /ZUDATA/ ACOLIN(0:11),EF_MIN(0:11),SPRIME(0:11),
     + ANGMIN(0:11),ANGMAX(0:11),SPRIPP(0:11),IRCUTS(0:11),IRFAST(0:11)
*
*-----------------------------------------------------------------------
*
* program version and date
* note: IZFVER = (x.y)*10; IZFDAT = yymmdd
*
      IZFVER = 643
      IZFDAT = 080617
*
* flag names
*
      CFLAGS(IFAFBC) = 'AFBC'
      CFLAGS(IFSCAL) = 'SCAL'
      CFLAGS(IFSCRE) = 'SCRE'
      CFLAGS(IFAMT4) = 'AMT4'
      CFLAGS(IFBORN) = 'BORN'
      CFLAGS(IFBOXD) = 'BOXD'
      CFLAGS(IFCONV) = 'CONV'
      CFLAGS(IFFINR) = 'FINR'
      CFLAGS(IFFOT2) = 'FOT2'
      CFLAGS(IFGAMS) = 'GAMS'
      CFLAGS(IFDIAG) = 'DIAG'
      CFLAGS(IFINTF) = 'INTF'
      CFLAGS(IFBARB) = 'BARB'
      CFLAGS(IFPART) = 'PART'
      CFLAGS(IFPOWR) = 'POWR'
      CFLAGS(IFPRNT) = 'PRNT'
      CFLAGS(IFALEM) = 'ALEM'
      CFLAGS(IFQCDC) = 'QCDC'
      CFLAGS(IFVPOL) = 'VPOL'
      CFLAGS(IFWEAK) = 'WEAK'
      CFLAGS(IFFTJR) = 'FTJR'
      CFLAGS(IFEXPR) = 'EXPR'
      CFLAGS(IFEXPF) = 'EXPF'
      CFLAGS(IFHIGS) = 'HIGS'
      CFLAGS(IFAFMT) = 'AFMT'
      CFLAGS(IFCZAK) = 'CZAK'
      CFLAGS(IFPREC) = 'PREC'
      CFLAGS(IFHIG2) = 'HIG2'
      CFLAGS(IFALE2) = 'ALE2'
      CFLAGS(IFGFER) = 'GFER'
      CFLAGS(IFISPP) = 'ISPP'
      CFLAGS(IFFSRS) = 'FSRS'
      CFLAGS(IFMISC) = 'MISC'
      CFLAGS(IFMISD) = 'MISD'
      CFLAGS(IFIPFC) = 'IPFC'
      CFLAGS(IFIPSC) = 'IPSC'
      CFLAGS(IFIPTO) = 'IPTO'
      CFLAGS(IFFBHO) = 'FBHO'
      CFLAGS(IFFSPP) = 'FSPP'
      CFLAGS(IFFUNA) = 'FUNA'
      CFLAGS(IFASCR) = 'ASCR'
      CFLAGS(IFSFSR) = 'SFSR'
      CFLAGS(IFENUE) = 'ENUE'
      CFLAGS(IFTUPV) = 'TUPV'
      CFLAGS(IFDMWW) = 'DMWW'
      CFLAGS(IFDSWW) = 'DSWW'
*
* default values
*
      IFLAGS(IFAFBC) = 1
      IFLAGS(IFSCAL) = 0
      IFLAGS(IFSCRE) = 0
      IFLAGS(IFAMT4) = 0
      IFLAGS(IFBORN) = 0
      IFLAGS(IFBOXD) = 2
      IFLAGS(IFCONV) = -1
      IFLAGS(IFFINR) = -1
      IFLAGS(IFFOT2) = 3
      IFLAGS(IFGAMS) = 1
      IFLAGS(IFDIAG) = 1
      IFLAGS(IFINTF) = 0
      IFLAGS(IFBARB) = 2
      IFLAGS(IFPART) = 0
      IFLAGS(IFPOWR) = 0
      IFLAGS(IFPRNT) = 1
      IFLAGS(IFALEM) = 0
      IFLAGS(IFQCDC) = 0
      IFLAGS(IFVPOL) = 1
      IFLAGS(IFWEAK) = 1
      IFLAGS(IFFTJR) = 0
      IFLAGS(IFEXPR) = 0
      IFLAGS(IFEXPF) = 0
      IFLAGS(IFHIGS) = 0
      IFLAGS(IFAFMT) = 1
      IFLAGS(IFCZAK) = 2
      IFLAGS(IFPREC) = 10
      IFLAGS(IFHIG2) = 0
      IFLAGS(IFALE2) = 3
      IFLAGS(IFGFER) = 2
      IFLAGS(IFISPP) = 2
      IFLAGS(IFFSRS) = 1
      IFLAGS(IFMISC) = 0
      IFLAGS(IFMISD) = 1
      IFLAGS(IFIPFC) = 5
      IFLAGS(IFIPSC) = 0
      IFLAGS(IFIPTO) = 3
      IFLAGS(IFFBHO) = 0
      IFLAGS(IFFSPP) = 0
      IFLAGS(IFFUNA) = 0
      IFLAGS(IFASCR) = 1
      IFLAGS(IFSFSR) = 1
      IFLAGS(IFENUE) = 1
      IFLAGS(IFTUPV) = 1
      IFLAGS(IFDMWW) = 0
      IFLAGS(IFDSWW) = 0
*
* experts only...
*
      IMASK = 0
      IMASS = 1
*
* cuts
*
      DO I = 0,11
*cbardin!!! default changed
        IRCUTS(I) =   1
        IRFAST(I) =   1
        ACOLIN(I) =   0D0
        EF_MIN(I) =   0D0
        SPRIME(I) =   0D0
        ANGMIN(I) =   0D0
        ANGMAX(I) = 180D0
        SPRIPP(I) =   0D0
      ENDDO
*
* print header
*
      IF(IPRINT.NE.0) THEN
      PRINT *,'******************************************************'
      PRINT *,'******************************************************'
      IVER = IZFVER/100
      ISEQ = MOD(IZFVER,100)
      PRINT '('' **           This is ZFITTER version '',I1,''.'',
     & I2.2,''           **'')',IVER,ISEQ
      IYY = IZFDAT/10000
      IMM = MOD(IZFDAT,10000)/100
      IDD = MOD(IZFDAT,100)
      PRINT '('' **                   '',I2.2,''/'',I2.2,''/'',I2.2,
     & ''                       **'')',IYY,IMM,IDD
      PRINT *,'******************************************************'
      PRINT *,'** http://www.ifh.de/theory/publist.html            **'
      PRINT *,'******************************************************'
      PRINT *
*
* print defaults
*
      PRINT *,'ZUINIT> ZFITTER defaults:'
      PRINT *
      CALL ZUINFO(0)
      CALL ZUINFO(1)
      ENDIF
*                                                             END ZUINIT
      END
 
      SUBROUTINE ZUINFO(MODE)
*     ========== ============
************************************************************************
*
*     SUBR. ZUINFO
*
*     Routine prints ZFITTER information.
*
*     MODE (int/read) = 0 print flag info
*                       1 print cut info
*
*     Called by USER,ZUINIT
*
************************************************************************
*
      IMPLICIT REAL*8(A-H,O-Z)
*
* flags
*
      PARAMETER(NFLGMX=46)
      COMMON /ZUFLGS/ IFLAGS(NFLGMX),CFLAGS(NFLGMX)
      CHARACTER CFLAGS*4
      PARAMETER (IFAFBC= 1,IFSCAL= 2,IFSCRE= 3,IFAMT4= 4,IFBORN= 5,
     & IFBOXD= 6,IFCONV= 7,IFFINR= 8,IFFOT2= 9,IFGAMS=10,IFDIAG=11,
     & IFINTF=12,IFBARB=13,IFPART=14,IFPOWR=15,IFPRNT=16,IFALEM=17,
     & IFQCDC=18,IFVPOL=19,IFWEAK=20,IFFTJR=21,IFEXPR=22,IFEXPF=23,
     & IFHIGS=24,IFAFMT=25,IFCZAK=26,IFPREC=27,IFHIG2=28,IFALE2=29,
     & IFGFER=30,IFISPP=31,IFFSRS=32,IFMISC=33,IFMISD=34,IFIPFC=35,
     & IFIPSC=36,IFIPTO=37,IFFBHO=38,IFFSPP=39,IFFUNA=40,IFASCR=41,
     & IFSFSR=42,IFENUE=43,IFTUPV=44,IFDMWW=45,IFDSWW=46)
*
* cuts
*
      COMMON /ZUDATA/ ACOLIN(0:11),EF_MIN(0:11),SPRIME(0:11),
     + ANGMIN(0:11),ANGMAX(0:11),SPRIPP(0:11),IRCUTS(0:11),IRFAST(0:11)
*
*-----------------------------------------------------------------------
*
* flags
*
      IF(MODE.EQ.0) THEN
        PRINT *,'ZFITTER flag values:'
        PRINT 9000,(CFLAGS(I),IFLAGS(I),I=1,NFLGMX)
        PRINT *
      ENDIF
*
* cuts
*
      IF(MODE.EQ.1) THEN
        PRINT *
        PRINT *,'ZFITTER cut values:'
        PRINT *,'  INDF  ICUT    ACOL    EMIN    S_PR',
     +   '    ANG0    ANG1     SPP'
        DO INDF = 4,4
        IF(IRCUTS(INDF).EQ.1.AND.IRFAST(INDF).EQ.1) THEN
          ICUTPR = -1
        ELSE
          ICUTPR = IRCUTS(INDF)
        ENDIF
          PRINT 9100,INDF,ICUTPR,ACOLIN(INDF),EF_MIN(INDF),SPRIME(INDF),
     +     ANGMIN(INDF),ANGMAX(INDF),SPRIPP(INDF)
        ENDDO
        PRINT *
        PRINT *
      ENDIF
      RETURN
 9000 FORMAT((5(1X,A4,':',I2)))
 9100 FORMAT(1X,2I6,6F8.2)
*                                                             END ZUINFO
      END
 
      SUBROUTINE ZUFLAG(CHFLAG,IVALUE)
*     ========== =====================
***************************************************************************
*                                                                         *
*     SUBR. ZUFLAG(CHFLAG,IVALUE)                                         *
*                                                                         *
*     Routine sets ZFITTER flags.  This routine must be called before     *
*     calling ZUWEAK.                                                     *
*                                                                         *
*     CHFLAG (char/read) = character identifer of ZFITTER flag            *
*     IVALUE (int/read)  = flag value                                     *
*                                                                         *
***************************************************************************
*                                                                         *
*     Flag   Value  Meaning                                               *
*                                                                         *
***************************************************************************
*     AFBC       0  NO AFB CALCULATION                                    *
*                1  BOTH SIGMA AND AFB                                    *
***************************************************************************
*     AFMT       0  WITHOUT AFMT CORRECTION                               *
*                1   WITH   AFMT CORRECTION O(Gf*m^2_t*als^2)             *
*                2   WITH   AFMT CORRECTION + O(Gf*M^2_z*als^2+log(m^2_t) *
*                3   WITH   AFMT CORRECTION + O(Gf*M^2_z/m^2_t*als^2)     *
***************************************************************************
*     ALEM       0 OR 2: DALH5(AMZ) MUST BE SUPPLIED BY THE USER AS       *
*                ------             INPUT TO THE DIZET PACKAGE            *
*                        USING THIS INPUT THE CODE CALCULATES ALPHA_EM    *
*                ---------------------------------------------------------*
*                1 OR 3: DALH5(AMZ) IS CALCULATED BY THE PROGRAM          *
*                ------             USING A PARAMETRIZATION (IHVP)        *
*                                   THEN ALPHA_EM IS CALCULATED           *
*                ---------------------------------------------------------*
*                The scale of ALPHA_EM(scale) is governed in addition     *
*                by the flag CONV:                                        *
*                Values ALEM=0,1 are accessible only at CONV=0.           * 
*                Then ALPHA_EM is calculated at M_Z for ALEM=0,1          *
*                              and at S for ALEM=2,3                      *
*                Values ALEM=2,3 are accessible for CONV=0,1,2.           *
*                Then ALPHA_EM is calculated at S for CONV=0              *
*                              and at S' for CONV=1,2                     *
*                ---------------------------------------------------------*
*                Recommended values: ALEM=2,3                             *
***************************************************************************
*     ALE2       0  WITHOUT TWO-LOOP CONSTANT CORRECTIONS IN DELTA_ALPHA  *
*                   THIS IS FOR A BACK COMPATIBILITY ONLY                 *
*                1  -- WITH   ONE-LOOP CORRECTIONS\                       *
*                2  -- WITH   TWO-LOOP CORRECTIONS | FOR LEPTONIC DALPHA  *
*                3  -- WITH THREE-LOOP CORRECTIONS/                       *
***************************************************************************
*     AMT4      -1  BACKWARDS COMPATIBILITY                               *
*                0  NO   MT**4 CORRECTIONS,                               *
*                1  WITH MT**4 CORRECTIONS RESUMED, SEE MPI-PREPRINT      *
*                2  WITH MT**4 CORRECTIONS AS IN HALZEN & KNIEHL (1990)   *
*                3  WITH MT**4 CORRECTIONS AS IN A.SIRLIN ET AL. (1990)   *
*                4  WITH 2EWRC CORRECTIONS AS IN DEGRASSI ET AL. (1998)   *
*                5  WITH 2EWRC CORRECTIONS AS IN WEIGLEIN ET AL. (2000)   *
*                6  WITH 2EWRC CORRECTIONS AS IN AWRAMIK ET AL. (2004)    *
***************************************************************************
*     BARB       0  19-2*PI2, MT > MH LIMIT                               *
*                1  R. BARBIERI APPROXIMATE EXPRESSION, MH > MT           *
*                2  R. BARBIERI EXACT EXPRESSION (FIT TO EXACT FORMULAE)  *
***************************************************************************
*     BORN       0  QED- CONVOLUTED OBSERVABLES                           *
*                1  EW- CORRECTED 'IMPROVED BORN'                         *
***************************************************************************
*     BOXD       0  WW AND ZZ-BOX CORRS.(WHICH ARE SMALL AND IN CONTRAST  *
*                   TO ELECTROWEAK FORM FACTORS DEPENDENT ON COS THETA)   *
*                   ARE NOT INCLUDED IN THE CROSS-SECTION                 *
*                1  THESE ARE INCLUDED VIA DZEWBX                         *
*                2  THESE ARE INCLUDED VIA ROKANC                         *
***************************************************************************
*     CONV      -1  ALPHA_EM(0)                                           *
*                0  ALPHA_EM(S)                                           *
*                1  ALPHA_EM(S_PRIME), CONVOLUTED                         *
*                2  EWRC AND ALSTR ARE CONVOLUTED                         *
***************************************************************************
*     CZAK       0  NONFACTORIZABLE O(\alpha\als) CORRECTIONS - OFF       *
*                1                                            - ON        *
*                2                          ON for PO's and OFF for RO's  *
***************************************************************************
*     DIAG      -1  ONLY  Z -EXCHANGE DIAGRAMS ARE TAKEN INTO ACCOUNT     *
*                0  Z- AND G-EXCHANGE DIAGRAMS ARE TAKEN INTO ACCOUNT     *
*                1  Z-, G-EXCHANGE AND ZG-INTERFERENCE ARE INCLUDED       *
***************************************************************************
*     EXPR       0  NON   -EXPANDED DR                                    *
*                1  PARTLY-EXPANDED DR                                    *
*                2  FULLY -EXPANDED DR, AN ANALOG OF OMS-II SCHEME        *
***************************************************************************
*     EXPF       0  NON   -EXPANDED RHO AND KAPPA analogs of   ff level   *
*                1  PARTLY-EXPANDED RHO AND KAPPA    those     ff level   *
*                2  FULLY -EXPANDED RHO AND KAPPA  in  IFACR   ff level   *
***************************************************************************
*     FBHO       0  PHOTONIC 2ND ORDER RC IN A_FB AS IN v.6_21            *
*                1  MODIFIED TREATMENT WITH LLA PAIRS (OPTIONALLY)        *
***************************************************************************
*     FINR      -1  FINAL STATE QED AND QCD CORRECTIONS ARE NOT APPLIED   *
*                0  FINAL STATE QED CORRECTION AS FACTOR =1+3/4/PI*QF**2  *
*                1  WITH FINAL STATE QED CORRECTIONS+EXPONENTIATION       *
***************************************************************************
*     FOT2      -1  NO ISR QED CONVOLUTION AT ALL                         *
*                0  COMPLETE ALPHA   ADDITIVE RADIATOR                    *
*                1  WITH LOGARITHMIC HARD CORRECTIONS                     *
*                   (supposed to be critically revisited)                 *
*                2  COMPLETE ALPHA^2 ADDITIVE RADIATOR                    *
*                3  COMPLETE ALPHA^3 ADDITIVE RADIATOR                    *
*                4  OPTIONAL ALPHA^3 ADDITIVE RADIATOR FOR TEOR. ERRORS   *
*                5  LLA YFS  ALPHA^3 FACTORIZED RADIATOR A'LA JADACH ET AL*
***************************************************************************
*     FSPP       0  FINAL STATE PAIRS ARE NOT ACTIVE                      *
*                1  FINAL STATE PAIRS ARE ACTIVE (ADDITIVE COR.)          *
*                2  FINAL STATE PAIRS ARE ACTIVE (MULTIPLICATIVE COR.)    *
***************************************************************************
*     FSRS       0  FINAL STATE RADIATION SCALE \alpha(0)                 *
*                1  FINAL STATE RADIATION SCALE \alpha(S)                 * 
***************************************************************************
*     FTJR       0  FLEISHER TARASOV JEGERLEHLER RACZKA - OFF             *
*                1  FTJR ON                                               *
*                2  FTJR ON for PO's and OFF for RO's                     *
***************************************************************************
*     FUNA       0  ORIGINAL, NON-CORRECTED ANGULAR DISTRIBUTION          *
*                1  NEWLY CORRECTED ANGULAR DISTRIBUTION (BY A. ARBUZOV)  *
*                   FOR HIGHER ORDER PHOTONIC LLA CORRECTIONS             *
***************************************************************************
*     GAMS       0  CONSTANT TOTAL Z0 WIDTH                               *
*                1  S DEPENDENT TOTAL Z0 WIDTH                            *
***************************************************************************
*     GFER       0  FOR BACK COMPATIBILITY WITH 5.12                      *
*                1  ONE-LOOP QED CORRECTIONS FOR FERMI CONSTANT           *
*                2  TW0-LOOP QED CORRECTIONS FOR FERMI CONSTANT           *
***************************************************************************
*     HIGS       0  LEADING HIGGS CONTRIBUTION IS NOT RESUMMED            *
*      I         1  LEADING HIGGS CONTRIBUTION IS     RESUMMED            *
***************************************************************************
*     HIG2       0  WITHOUT TWO-LOOP HIGGS  CORRECTIONS                   * 
*                1  WITH   --/--/--/--/--/--/--/--/--/                    * 
***************************************************************************
*     MISC       0  (default) non-scaled \rho's are used, AROTFZ-array    *
*     MISC       1                scaled \rho's are used, ARROFZ-array    *
***************************************************************************
*     MISD          Model Independent, S-Dependence                       *
*     MISD       0  default, S=M^2_z in EWRC, old treatment               *
*     MISD       1  ensures equal numbers from all interfaces             *
*                   for all \sqrt{S} up to 100 GeV                        *
***************************************************************************
*     INTF       0  NO INITIAL-FINAL INTERFERENCE                         *
*                1  WITH I-F INT. IN THE LOWEST ORDER IN ALPHA            *
*                2  EXPONENTIATED  I-F INT.                               *
***************************************************************************
*     IPFC          Pair Flavour Content                                  *
*                1  only electron   pairs                                 *
*                2  only muon       pairs                                 *
*                3  only tau-lepton pairs                                 *
*                4  only hadron     pairs                                 *
*                5  all summed                                            *
*                6  leptonic pairs (without hadrons)                      *
***************************************************************************
*     IPSC          Pair Singlet-channel Contributions                    *
*                   works with ISPP = 2                                   *
*                0  only non-singlet pairs                                *
*                1  LLA singlet pairs according to BNB                    *
*                2  complete O(alpha^2) singlet pairs (BNB)               *
*                   (BNB=Berends, van Neerven, Burgers)                   *
*                3  Singlet pairs up to (alpha*L)**3                      *
***************************************************************************
*     IPTO          Pair Third (and higher) Order contributions           *
*                   works with ISPP = 2                                   *
*                   (see A.B. Arbuzov, hep-ph/9907500)                    *
*               -1  only O(alpha^2) virtual pairs (new)                   *
*                0  only O(alpha^2) contributions                         *
*                1  O(alpha^3) pairs                                      *
*                2  some "non-standard" O(alpha^3) LLA pairs added        *
*                3  O(alpha^4) LLA electron pairs added                   *
***************************************************************************
*     ISPP      -1  IS PAIRS ARE TREATED AS VERSIONS UP TO 5.14           *
*                0  WITHOUT ISR PAIRS                                     *
*                1  WITH    ISR PAIRS, KKKS with re-weighting, PCP version*
*                2  WITH    ISP PAIRS ACCORDING TO A.B. Arbuzov           *
*                3  WITH    ISR PAIRS, JMS                                *
*                4  WITH    ISR PAIRS, JMS + extended treatment of hadr.  *
***************************************************************************
*     PART       0  (S+T) CROSS-SECTION AND ASYMMETRY WITH BHANG INDF=11) *
*                1   S CHANNEL ONLY                                       *
*                2   T CHANNEL ONLY                                       *
*                3   S-T INTERFERENCE ONLY                                *
***************************************************************************
*     POWR       0  FINAL FERMION MASSES IN KINEMATICAL FACTORS SET TO 0  *
*                   (FOR TESTS ONLY)                                      *
*                1  FINAL RERMION MASSES IN KINEMATICAL FACTORS RETAINED  *
***************************************************************************
*     PREC          PREC=1 default                                        *
*  1 < PREC < 99    is an integer number which ANY precision govering ANY *
* numerical integration is devided to, increasing thereby the numerical   *
* precision of computation. In some cases when some numerical instability *
* while running 5_10 was registered, it was sufficient to use PREC=3      *
* in some other cases (Ptau) only PREC=30 solved the instability.         *
***************************************************************************
*     PRNT       0  SUPRESSES ZUWEAK PRINTING                             *
*                1  ZUWEAK PRINTS VARIOUS QUANTITIES ON EACH CALL         *
***************************************************************************
*     QCDC       0  NO QCD CORRS. TO VECTOR BOSON SELF ENERGIES           *
*                   IN DELTA-R, WIDTHS, CROSS SECTION                     *
*                1  APPROX. FAST QCD CORRS. (NOT REALISED FOR  W-WIDTH)   *
*                   IMPORTANT NOTICE: THESE ARE RELIABLE ONLY FOR LEP-I   *
*                2  EXACT FORMULAE (SEE REF. IN DIZET)                    *
*                3  KNIEH'S QCD-LIBRARY                                   *
***************************************************************************
*     SCRE       0  SCALE OF THE REMAINDERS = 1                           *
*                1  ----------------------- RENORD                        *
*                2  ----------------------- RENORM (not recommended, n.r.)*
***************************************************************************
*     SCAL       0  alpha_s(m_t)                                          *
*                1  alpha_s(xi*m_t)+delta  |                              *
*                2  alpha_s(xi*m_t)        | Kniehl's                     *
*                3  alpha_s(xi*m_t)-delta  |                              *
*                4  alpha_s(.35*m_t)                                      *
***************************************************************************
*     VPOL       1  HADR. VAC. POLARIZATION OF JEGERLEHNER(1991)          *
*                2                          OF JEGERLEHNER(1988)          *
*                3                          OF BURKHARDT  (1989)          *
***************************************************************************
*     WEAK      -1  VALID ONLY FOR INDF=-1, THE SAME AS 1 AND 2 BUT       *
*                   WITHOUT EW CORRERTIONS FOR W t-CHANNEL EXCHANGE       *
*                0  WEAK LOOP CORRECTIONS OFF                             *
*                1  WEAK LOOPS ON                                         *
*                2  SOME HIGHER ORDER EWRC THAT DO NOT PROPAGATE VIA      *
*                   DIZET ARE DISABLED (ADDIME, ADDIMF). THIS IS DONE     *
*                   BECAUSE THEY ARE SMALL AT LEP2 AND BECAUSE THEY CAN'T *
*                   BE USED BY THE OTHER CODES WHICH USE ONLY DIZET       *
***************************************************************************
* Hidden flaggs:                                                          *
* Asymmetry crisis:                                                       *
*     ASCR       1  Very old treatment: ZF4... to ZF5...                  *
*                0  IZA treatment: ZF5... to ZF6.23                       *
*                1  NTA treatment: ZF6.30 and later                       *
***************************************************************************
* Jadach's switch:                                                        *
*     SFSR       1  all included                                          *
*                0  3/4*... QED is excluded                               *
*               -1  both QED and mixed QED x QCD are excluded             *
***************************************************************************
*                                                                         *
*     Called by USER                                                      *
*                                                                         *
***************************************************************************
*
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER CHFLAG*(*)
*
* flags
*
      PARAMETER(NFLGMX=46)
      COMMON /ZUFLGS/ IFLAGS(NFLGMX),CFLAGS(NFLGMX)
      CHARACTER CFLAGS*4
      PARAMETER (IFAFBC= 1,IFSCAL= 2,IFSCRE= 3,IFAMT4= 4,IFBORN= 5,
     & IFBOXD= 6,IFCONV= 7,IFFINR= 8,IFFOT2= 9,IFGAMS=10,IFDIAG=11,
     & IFINTF=12,IFBARB=13,IFPART=14,IFPOWR=15,IFPRNT=16,IFALEM=17,
     & IFQCDC=18,IFVPOL=19,IFWEAK=20,IFFTJR=21,IFEXPR=22,IFEXPF=23,
     & IFHIGS=24,IFAFMT=25,IFCZAK=26,IFPREC=27,IFHIG2=28,IFALE2=29,
     & IFGFER=30,IFISPP=31,IFFSRS=32,IFMISC=33,IFMISD=34,IFIPFC=35,
     & IFIPSC=36,IFIPTO=37,IFFBHO=38,IFFSPP=39,IFFUNA=40,IFASCR=41,
     & IFSFSR=42,IFENUE=43,IFTUPV=44,IFDMWW=45,IFDSWW=46)
*
      INTEGER*4 IFLGMN(NFLGMX),IFLGMX(NFLGMX)
      DATA IFLGMN(IFAFBC)/0/  ,IFLGMX(IFAFBC)/1/
      DATA IFLGMN(IFSCAL)/0/  ,IFLGMX(IFSCAL)/4/
      DATA IFLGMN(IFSCRE)/0/  ,IFLGMX(IFSCRE)/2/
      DATA IFLGMN(IFAMT4)/-1/ ,IFLGMX(IFAMT4)/6/
      DATA IFLGMN(IFBORN)/0/  ,IFLGMX(IFBORN)/1/
      DATA IFLGMN(IFBOXD)/0/  ,IFLGMX(IFBOXD)/2/
      DATA IFLGMN(IFCONV)/-1/ ,IFLGMX(IFCONV)/2/
      DATA IFLGMN(IFFINR)/-1/ ,IFLGMX(IFFINR)/2/
      DATA IFLGMN(IFFOT2)/-1/ ,IFLGMX(IFFOT2)/5/
      DATA IFLGMN(IFGAMS)/0/  ,IFLGMX(IFGAMS)/1/
      DATA IFLGMN(IFDIAG)/-1/ ,IFLGMX(IFDIAG)/1/
      DATA IFLGMN(IFINTF)/0/  ,IFLGMX(IFINTF)/2/
      DATA IFLGMN(IFBARB)/0/  ,IFLGMX(IFBARB)/2/
      DATA IFLGMN(IFPART)/0/  ,IFLGMX(IFPART)/3/
      DATA IFLGMN(IFPOWR)/0/  ,IFLGMX(IFPOWR)/1/
      DATA IFLGMN(IFPRNT)/0/  ,IFLGMX(IFPRNT)/1/
      DATA IFLGMN(IFALEM)/0/  ,IFLGMX(IFALEM)/3/
      DATA IFLGMN(IFQCDC)/0/  ,IFLGMX(IFQCDC)/3/
      DATA IFLGMN(IFVPOL)/1/  ,IFLGMX(IFVPOL)/3/
      DATA IFLGMN(IFWEAK)/-1/ ,IFLGMX(IFWEAK)/2/
      DATA IFLGMN(IFFTJR)/0/  ,IFLGMX(IFFTJR)/2/
      DATA IFLGMN(IFEXPR)/0/  ,IFLGMX(IFEXPR)/2/
      DATA IFLGMN(IFEXPF)/0/  ,IFLGMX(IFEXPF)/2/
      DATA IFLGMN(IFHIGS)/0/  ,IFLGMX(IFHIGS)/1/
      DATA IFLGMN(IFAFMT)/0/  ,IFLGMX(IFAFMT)/3/
      DATA IFLGMN(IFCZAK)/0/  ,IFLGMX(IFCZAK)/2/
      DATA IFLGMN(IFPREC)/1/  ,IFLGMX(IFPREC)/99/
      DATA IFLGMN(IFHIG2)/0/  ,IFLGMX(IFHIG2)/1/
      DATA IFLGMN(IFALE2)/0/  ,IFLGMX(IFALE2)/3/
      DATA IFLGMN(IFGFER)/0/  ,IFLGMX(IFGFER)/2/
      DATA IFLGMN(IFISPP)/-1/ ,IFLGMX(IFISPP)/4/
      DATA IFLGMN(IFFSRS)/0/  ,IFLGMX(IFFSRS)/1/
      DATA IFLGMN(IFMISC)/0/  ,IFLGMX(IFMISC)/1/
      DATA IFLGMN(IFMISD)/0/  ,IFLGMX(IFMISD)/1/
      DATA IFLGMN(IFIPFC)/1/  ,IFLGMX(IFIPFC)/6/
      DATA IFLGMN(IFIPSC)/0/  ,IFLGMX(IFIPSC)/3/
      DATA IFLGMN(IFIPTO)/-1/ ,IFLGMX(IFIPTO)/3/
      DATA IFLGMN(IFFBHO)/0/  ,IFLGMX(IFFBHO)/1/
      DATA IFLGMN(IFFSPP)/0/  ,IFLGMX(IFFSPP)/2/
      DATA IFLGMN(IFFUNA)/0/  ,IFLGMX(IFFUNA)/1/
      DATA IFLGMN(IFASCR)/0/  ,IFLGMX(IFASCR)/1/
      DATA IFLGMN(IFSFSR)/-1/ ,IFLGMX(IFSFSR)/1/
      DATA IFLGMN(IFENUE)/-1/ ,IFLGMX(IFENUE)/1/
      DATA IFLGMN(IFTUPV)/ 1/ ,IFLGMX(IFTUPV)/3/
      DATA IFLGMN(IFDMWW)/-1/ ,IFLGMX(IFDMWW)/1/
      DATA IFLGMN(IFDSWW)/-1/ ,IFLGMX(IFDSWW)/1/
*
*-----------------------------------------------------------------------
*
      DO 10 I = 1,NFLGMX
        IF(CHFLAG.NE.CFLAGS(I)) GOTO 10
        IF(IVALUE.LT.IFLGMN(I) .OR. IVALUE.GT.IFLGMX(I)) THEN
          PRINT *,'ZUFLAG> IVALUE out of range; CHFLAG/IVALUE = ',
     +     CHFLAG,' ',IVALUE
          RETURN
        ENDIF
        IFLAGS(I) = IVALUE
        RETURN
   10 ENDDO
      PRINT*,'ZUFLAG> CHFLAG not found; CHFLAG = ',CHFLAG
      RETURN
*                                                             END ZUFLAG
      END
 
      SUBROUTINE ZUWEAK(ZMASS,TMASS,HMASS,DAL5H,ALFAS)
*     ========== =====================================
*     ZUWEAK-WRAPPER
*     ==============
      IMPLICIT REAL*8(A-Z)
      V_TB=1D0
      CALL  ZVWEAK(ZMASS,TMASS,HMASS,DAL5H,V_TB,ALFAS)
      RETURN
      END

      SUBROUTINE ZVWEAK(ZMASS,TMASS,HMASS,DAL5H,V_TB,ALFAS)
*     ========== ==========================================
************************************************************************
*                                                                      *
*     SUBR. ZVWEAK(ZMASS,TMASS,HMASS,DAL5H,V_TB,ALFAS)                 *
*     SUBR. ZUWEAK(ZMASS,TMASS,HMASS,DAL5H,     ALFAS)                 *
*                                                                      *
*     Routine initializes weak part of ZFITTER.  (Replaces ZINITF in   *
*     pre 3.05 ZFITTER versions)                                       *
*                                                                      *
*     ZMASS (real/read) = Z0 mass (GeV)                                *
*     TMASS (real/read) = top mass (GeV)                               *
*     HMASS (real/read) = Higgs mass (GeV)                             *
*     V_TB  (real/read) = the Vtb mixing matrix element                *
*     ALFAS (real/read) = the strong coupling constant (Q**2 = M(Z)**2)*
*                                                                      *
*     Called by USER                                                   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8(A-H,O-W,Y-Z)
      IMPLICIT COMPLEX*16(X)
      PARAMETER(ALFAI=137.0359895D0)
      COMMON /FRINIT/ NPAR(30),ZPAR(31)
      COMMON /PSCONS/ SW2,AMZ,GAMZ
      DIMENSION XFZ(4)
*
      COMMON /PARTZW/ PARTZ(0:11),PARTW(3)
      COMMON /EWFORM/ XALLCH(5,4),XFOTF
      COMMON /ZFCHMS/ ALLCH(0:11),ALLMS(0:11)
*
      COMMON /CDZRKZ/ARROFZ(0:10),ARKAFZ(0:10),ARVEFZ(0:10),ARSEFZ(0:10)
     &              ,AROTFZ(0:10),AIROFZ(0:10),AIKAFZ(0:10),AIVEFZ(0:10)
*
      COMMON /CDZFER/CLM(8),AML(8),CQM(8),AMQ(8),VB,VT,VB2,VB2T,VT2,VT2T
      COMMON /CDZTHR/ AMTH(6)
      COMMON /CDZRUN/ CMQRUN(8)

*  ls
*     COMMON/ZPARD_/ZPARD

      COMMON/DR1/DR1
      DIMENSION ACOST(3),ASQS(3)
      DATA ACOST/-0.9D0,0,0.9D0/,ASQS/100D0,200D0,300D0/
* ls

*
* flags
*
      PARAMETER(NFLGMX=46)
      COMMON /ZUFLGS/ IFLAGS(NFLGMX),CFLAGS(NFLGMX)
      CHARACTER CFLAGS*4
      PARAMETER (IFAFBC= 1,IFSCAL= 2,IFSCRE= 3,IFAMT4= 4,IFBORN= 5,
     & IFBOXD= 6,IFCONV= 7,IFFINR= 8,IFFOT2= 9,IFGAMS=10,IFDIAG=11,
     & IFINTF=12,IFBARB=13,IFPART=14,IFPOWR=15,IFPRNT=16,IFALEM=17,
     & IFQCDC=18,IFVPOL=19,IFWEAK=20,IFFTJR=21,IFEXPR=22,IFEXPF=23,
     & IFHIGS=24,IFAFMT=25,IFCZAK=26,IFPREC=27,IFHIG2=28,IFALE2=29,
     & IFGFER=30,IFISPP=31,IFFSRS=32,IFMISC=33,IFMISD=34,IFIPFC=35,
     & IFIPSC=36,IFIPTO=37,IFFBHO=38,IFFSPP=39,IFFUNA=40,IFASCR=41,
     & IFSFSR=42,IFENUE=43,IFTUPV=44,IFDMWW=45,IFDSWW=46)
*
      COMMON /EXPERT/ IMASK,IMASS
      COMMON /PRECIS/ NPREC
*
* parameters
*
      COMMON /ZUPARS/QDF,QCDCOR(0:14),ALPHST,SIN2TW,S2TEFF(0:11),
     & WIDTHS(0:11)
*
* cuts
*
      COMMON /ZUDATA/ ACOLIN(0:11),EF_MIN(0:11),SPRIME(0:11),
     + ANGMIN(0:11),ANGMAX(0:11),SPRIPP(0:11),IRCUTS(0:11),IRFAST(0:11)
*
      CHARACTER FCHAN(0:11)*10,FCHANW(1:3)*10
      DATA FCHAN/'nu,nubar','e+,e-','mu+,mu-','tau+,tau-','u,ubar',
     + 'd,dbar','c,cbar','s,sbar','t,tbar','b,bbar','hadron','total'/
      DATA FCHANW/'lept,nubar','down,ubar','total'/
*
*-----------------------------------------------------------------------
*
* filling of COMMON /PRECIS/
*
      NPREC=IFLAGS(IFPREC)
*
* call DIZET
*
      CALL ZDIZET(ZMASS,TMASS,HMASS,DAL5H,V_TB,ALFAS,ALQED,ALFAT)
*
* FILLING OF FLAG VECTOR NPAR FOR ZFITTER
*
      NPAR(1)  = IFLAGS(IFWEAK)
      NPAR(2)  = IFLAGS(IFVPOL)
      NPAR(3)  = IFLAGS(IFQCDC)
      NPAR(4)  = IFLAGS(IFBOXD)
      NPAR(5)  = IFLAGS(IFGAMS)
      NPAR(7)  = IFLAGS(IFDIAG)
      NPAR(8)  = IFLAGS(IFINTF)
      NPAR(9)  = IFLAGS(IFFINR)
      NPAR(10) = IFLAGS(IFFOT2)
      NPAR(13) = IFLAGS(IFAFBC)
      NPAR(14) = IFLAGS(IFBORN)
      NPAR(15) = IFLAGS(IFSCAL)
      NPAR(16) = IFLAGS(IFAMT4)
      NPAR(17) = IFLAGS(IFPART)
      NPAR(18) = IMASK
      NPAR(19) = IFLAGS(IFEXPF)
      NPAR(20) = IFLAGS(IFALEM)
      NPAR(21) = IFLAGS(IFALE2)
*
* PRE-FILLING OF PARAMETER VECTOR ZPAR FOR ZFITTER (SOME OF THEM WILL BE
*                                                   REINITIALISED LATER)
* DESCRIPTION OF PARAMETER VECTOR ZPAR
*
* ZPAR   I
*  1     I QE ...... ELECTRON BEAM CHARGE = -1.
*  2     I QF .......FINAL FERMION CHARGE
*  3     I COLOR ... FERMION COLOR FACTOR (1 OR 3)
*  4     I AMF...... FERMION MASS
*  5     I AMZ...... Z-BOSON MASS
*  6     I ALFAS.... STRONG COUPLING CONSTANT (.12)
*7-21,31 I QCDCOR... ALST/PI (1.04)
* 22     I ALAMP.... AVERAGE POSITRON BEAM POLARIZATION (BETWEEN -1.D0
* 23     I ALAME....         ELECTRON                    AND +1.D0 )
* 24     I HELP..... HELICITY OF POSITIVE OUTGOING FERMION
* 25     I HELM.....             NEGATIVE
*                    HELICITY -1.D0 OR +1.D0 OR, IF AVEREGED, EQUAL 0.D0
* 26     I DELTA.... PHOTON CUT PARAMETER
*                    =2*K(ZERO)MAX/SQRT(S)
*                    =1. - S(PRIME)/S
* 27     I ACOL..... ACOLLINARITY CUT IN DEGREES
* 28     I EMIN..... MUON ENERGIES CUT IN GEV
* 29     I ANG1..... ANG1,ANG2 ARE LIMITS ON ANGLE-INTEGRATION I.E.
* 30     I ANG2..... GEOMETRICAL ACCEPTANCE CUTS. THEY MAY BE EVEN NOT
*        I ......... EQUAL AND LIE IN ONE HEMISPHERE. ZCUT WILL RETURN
*        I ......... SIGMA=INT(FROM C1 TO C2) OF DSIGMA/DCOSTH
* 31     I SPRIPP... A CUT FOR PAIR PRODUCTION
*
      ZPAR(1)  =-1.D0
      ZPAR(2)  = 0.D0
      ZPAR(3)  = 1.D0
      ZPAR(4)  = 0.D0
      ZPAR(5)  = ZMASS
      DO IQCDC=0,14
       ZPAR(7+IQCDC) = QCDCOR(IQCDC)
      ENDDO
      ZPAR(22) = 0.D0
      ZPAR(23) = 0.D0
      ZPAR(24) = 0.D0
      ZPAR(25) = 0.D0
      ZPAR(26) = 0.D0
      ZPAR(27) = 1.D1
      ZPAR(28) = 1.D1
      ZPAR(29) = 180.D0
      ZPAR(30) = 0.D0
      ZPAR(31) = 0.D0
*
* MORE FILLING
*
      AMZ =ZMASS
      GAMZ=PARTZ(11)/1D3
      SW2 =SIN2TW
*
      DO 1 I=0,11
        WIDTHS(I)=PARTZ(I)
 1    CONTINUE
      IF(IFLAGS(IFPRNT).EQ.1) THEN
        PRINT *
        PRINT *,'ZFITTER input parameters:'
        PRINT *,'DAL5H =',DAL5H
        PRINT *,'ALQED5=',1D0/ALQED
        PRINT 9000,'ZMASS',ZMASS,'TMASS',TMASS
        PRINT 9010,'HMASS',HMASS
        PRINT 9000,'DAL5H',DAL5H,'ALQED5',1D0/ALQED
        PRINT 9000,'ALFAS',ALFAS,'ALFAT',ALFAT
        PRINT *
        PRINT *,'ZFITTER intermediate results:'
        WMASS=ZMASS*SQRT(1D0-SIN2TW)
        PRINT 9000,'WMASS',WMASS,'SIN2TW',SIN2TW
        PRINT *
        PRINT 9000,'ALPHST',ALPHST
        PRINT 9100,'QCDCOR',QCDCOR
        PRINT *
        PRINT *,'CHANNEL         WIDTH         RHO_F_R        RHO_F_T '
     + ,'       SIN2_EFF'
        PRINT *,'-------        -------       --------       --------'
     + ,'       --------'
        DO 2 I=0,9
          PRINT '(1X,A10,2X,F10.3,7X,F8.6,7X,F8.6,7X,F8.6)',
     +    FCHAN(I),WIDTHS(I),ARROFZ(I),AROTFZ(I),ARSEFZ(I)
 2      CONTINUE
        DO 21 I=10,11
          PRINT '(1X,A10,2X,F10.3)',FCHAN(I),WIDTHS(I)
 21     CONTINUE
        PRINT *
        PRINT *,'W-widths'
        DO 20 I=1,3
          PRINT '(1X,A10,2X,F10.3)',FCHANW(I),PARTW(I)
 20     CONTINUE

      ENDIF
*
* INITIALISATION OF PARTIAL CHANNELS
*
      IBOXF=0
*
      IBFLA = 0
      QE    = -1.D0
*
      S     = ZMASS**2
      Q2    = S/2D0
      U2    = Q2-S
*
c$$$* NEUTRINO CHANNELS INITIALIZATION
c$$$*
c$$$      QF       = 0.D0
c$$$      ALLCH(0) = QF
c$$$      ALLMS(0) = 0.D0
c$$$      CALL ROKANC(IBOXF,IBFLA,U2,-S,-Q2,QE,QF,XFZ,XFZT,XFZTA)
c$$$      DO 3 I = 1,4
c$$$ 3    XALLCH(1,I) = XFZ(I)
c$$$      S2TEFF(0) = SIN2TW*DREAL(XFZ(3))
c$$$*
c$$$* LEPTON   CHANNELS INITIALIZATION
c$$$*
c$$$      QF=-1D0
c$$$      DO 4 I = 1,3
c$$$ 4    ALLCH(I) = QF
c$$$      ALLCH(11)= QF
c$$$      ALLMS(1) = AML(2)
c$$$      ALLMS(11)= AML(2)
c$$$      ALLMS(2) = AML(4)
c$$$      ALLMS(3) = AML(6)
c$$$      CALL ROKANC(IBOXF,IBFLA,U2,-S,-Q2,QE,QF,XFZ,XFZT,XFZTA)
c$$$      DO 5 I = 1,4
c$$$ 5    XALLCH(2,I) = XFZ(I)
c$$$      S2TEFF(1) = SIN2TW*DREAL(XFZ(3))
c$$$      S2TEFF(2) = SIN2TW*DREAL(XFZ(3))
c$$$      S2TEFF(3) = SIN2TW*DREAL(XFZ(3))
c$$$*
* U-QUARK CHANNELS INITIALIZATION
*
      QE=-1
      QF       = 2.D0/3D0
      ALLCH(4) = QF
      ALLCH(6) = QF
      ALLCH(8) = QF
      IF(IMASK.EQ.0) THEN
        ALLMS(4) = AMQ (1)
        ALLMS(6) = CMQRUN(3)
          ELSE
        ALLMS(4) = AMTH(1)
        ALLMS(6) = AMTH(3)
      ENDIF
      ALLMS(8) = TMASS
      S=10000

      DO IS = 1,3
      DO IC = 1,3
         COCOST=ACOST(IC)
         S=ASQS(IS)**2
       print *,S,COCOST
       Q2=S/2D0*(1D0-COCOST)
       U2=Q2-S
       IBOXF=2
       IBFLA=0
      
      CALL ROKANC(IBOXF,IBFLA,U2,-S,-Q2,QE,QF,XFZ,XFZT,XFZTA)
      DO 7 I = 1,4
 7    XALLCH(3,I) = XFZ(I)
      S2TEFF(4) = SIN2TW*DREAL(XFZ(3))
      S2TEFF(6) = SIN2TW*DREAL(XFZ(3))
      S2TEFF(8) = 0.D0
*      print*,"1-IAMT4=",IAMT4
*      print*,"eett"
*      print *,'after ROKANC'
      PI=ATAN(1D0)*4D0
      ALPHA=1/137.035 989 5 D0
      AL4PI=ALPHA/4/PI

      DR=DR1
      DR=0d0
*       DR=ZPARD(1)! *AL4PI
       
       FLL=(XFZ(1))*SW2/AL4PI
*      FLL=(XFZ(1)-1D0+DR)*SW2/AL4PI
       FQL=FLL+(XFZ(2)-1D0)*SW2/AL4PI
       FLQ=FLL+(XFZ(3)-1D0)*SW2/AL4PI
       FQQ=FLL+(XFZ(4)-1D0)*SW2/AL4PI
c$$$       print *,'AL4PI =',AL4PI,ALPHA,PI       
c$$$       print *,'RS =',RS
c$$$       print *,'DR =',DR
c$$$       print *,'SW2=',SW2
c$$$       print *,'AMZ=',ZMASS
c$$$       print *,'AMW=',WMASS      
c$$$       print *,'AMH=',HMASS       
c$$$       print *,'AMT=',TMASS       
c$$$       print *,'AL4PI=',AL4PI    
       print *,'FLL=',FLL
*       print *,'FQL=',FQL
*       print *,'FLQ=',FLQ
*       print *,'FQQ=',FQQ
*       print *,'XFZ=',XFZ
*      stop
      ENDDO
      ENDDO   
      stop
*     
* D-QUARK CHANNEL INITIALIZATION
*
      QF       = -1.D0/3D0
      ALLCH(5) = QF
      ALLCH(7) = QF
      IF(IMASK.EQ.0) THEN
        ALLMS(5) = AMQ (2)
        ALLMS(7) = AMQ (4)
      ELSE
        ALLMS(5) = AMTH(2)
        ALLMS(7) = AMTH(4)
      ENDIF
      CALL ROKANC(IBOXF,IBFLA,U2,-S,-Q2,QE,QF,XFZ,XFZT,XFZTA)
      DO 9 I = 1,4
 9    XALLCH(4,I) = XFZ(I)
      S2TEFF(5) = SIN2TW*DREAL(XFZ(3))
      S2TEFF(7) = SIN2TW*DREAL(XFZ(3))
*
* T-QUARK CHANNEL INITIALIZATION
*
      ALLMS(8)=AMQ(5)
*
* B-QUARK CHANNEL INITIALIZATION
*
      IBFLA=1
      QF       = -1.D0/3D0
      ALLCH(9) = QF
      IF(IMASK.EQ.0) THEN
        ALLMS(9) = CMQRUN(6)
      ELSE
        ALLMS(9) = AMTH(6)
      ENDIF
      CALL ROKANC(IBOXF,IBFLA,U2,-S,-Q2,QE,QF,XFZ,XFZT,XFZTA)
      DO 11 I = 1,4
 11   XALLCH(5,I) = XFZ(I)
      S2TEFF(9) = SIN2TW*DREAL(XFZ(3))
*
* Here XFOTF is always at M^2_Z
*
      XFOTF      = XFZT
*
* Attention, ALLCH for INDF=10 is set to 1 in order to account for IFI 
*
      IF(IFLAGS(IFMISD).EQ.1) THEN
        ALLCH(10)  = 1D0
      ELSE
        ALLCH(10)  = 1D0
      ENDIF
      ALLMS(10)  = 0D0
      S2TEFF(10) = 0D0
      S2TEFF(11) = 0D0
*
      RETURN
 9000 FORMAT(1X,A6,' = ',F10.5,';',1X,A6,' = ',F10.5)
 9010 FORMAT(1X,A6,' = ',F10.5)
 9100 FORMAT(1X,A6,' = ',16(F8.5))
*                                                             END ZUWEAK
      END

      SUBROUTINE ZUCUTS(INDF,ICUT,ACOL,EMIN,S_PR,ANG0,ANG1,SIPP)
*=============== ===============================================
* Modified by T. Riemann on 04/03/99 04:00am
* Modified by T.K    18.05.1999
************************************************************************
*                                                                      *
*     SUBR. ZUCUTS(INDF,ICUT,ACOL,EMIN,S_PR,ANG0,ANG1,SIPP)            *
*                                                                      *
*     Routine defines ZFITTER cuts for fermion channel.                *
*                                                                      *
*     INDF (int/read)  = fermion index                                 *
*     ICUT (int/read)  = -1 no cuts to be used (original code)         *
*                      =  0 ACOL+EMIN cuts     (original code)         *
*                      =  1 S_PR cut           (original code)         *
*          =  2 ACOL+EMIN cuts, no ACCEPTANCE CUT (new code)           *
*          =  3 ACOL+EMIN+ACCEPTANCE cuts         (new code)           *
*                                                                      *
*     ACOL (real/read) = acolinearity cut                              *
*     EMIN (real/read) = minimum fermion energy                        *
*     S_PR (real/read) = s' (variable related to maximum photon energy)*
*     ANG0 (real/read) = minimum polar angle theta (deg)               *
*     ANG1 (real/read) = maximum polar angle theta (deg)               *
*     SIPP (real/read) = s cut for for the secondary pairs, ISPP+FSPP  *
*                                                                      *
*     Called by USER                                                   *
*                                                                      *
************************************************************************
*  ZUCUTS:              parameter | internal ZFITTER variables         *
*----------------------------------------------------------------------*
*                            ICUT | IRCUTS/IRCUT IRFAST/IFAST IFUNFIN  *
* s'(M^2)/no-theta cuts   (old)-1 |       1            1        1      *
* acol/emin/theta cuts    (old) 0 |       0            0        0      *
* s'(M^2)/theta cuts      (old) 1 |       1            0        1      *
* acol/emin/no-theta cuts (new) 2 |       2            1        1      *
* acol/emin/theta cuts    (new) 3 |       3            0        1      *
*----------------------------------------------------------------------*
*
      IMPLICIT REAL*8(A-H,O-Z)
*
* cuts
*
      COMMON /ZUDATA/ ACOLIN(0:11),EF_MIN(0:11),SPRIME(0:11),
     + ANGMIN(0:11),ANGMAX(0:11),SPRIPP(0:11),IRCUTS(0:11),IRFAST(0:11)
*
      COMMON /INDFIN/ IFUNFIN
*
* Sets ifunfin = 0 if one wants to reproduce 5.xx numbers
*
      DATA NNN0,NNN2,NNN3/1,1,1/
*
*-----------------------------------------------------------------------
*
      IF(INDF.LT.-1.OR.INDF.GT.11) THEN
        PRINT *,'ZUCUTS> fermion index out of range:  INDF = ',INDF
        RETURN
      ENDIF
*
      IF(ICUT.LT.-1.OR.ICUT.GT.3) THEN
        PRINT *,'ZUCUTS> cut number out of range:  ICUT = ',ICUT      
        RETURN
      ENDIF
*
      IFUNFIN=1
*
      IRCUTS(INDF) = ABS(ICUT)
*
      if(ICUT.eq.-1) THEN
        IFUNFIN=1
        IRFAST(INDF) = 1
      endif
*
cbardin
cb    if(ICUT.eq.0.and.INDF.ne.11) THEN
      if(ICUT.eq.0) THEN
        IFUNFIN = 0
        IRFAST(INDF) = 0
        if(NNN0.EQ.1)  then
          PRINT*,'USE OF AN OBSOLETED OPTION, ICUT=0;'
          PRINT*,'Might be useful for backcompatibility with 5.xx;'
          PRINT*,'Presently, ICUT=2,3 recommended for realistic cuts.'
          NNN0=NNN0+1
        endif
      endif
*
      if(ICUT.eq.1) THEN
        IFUNFIN=1
        IRFAST(INDF)=0
      endif
*
      if(ICUT.eq.2) THEN
        IFUNFIN = 1
        IRFAST(INDF)=1
        IF(NNN2.EQ.1) then
          PRINT *,'USE OF PACKAGE acol.f'
          PRINT *,'by P. CHRISTOVA, M. JACK, T. RIEMANN,' 
          PRINT *,'DESY 99-015, HEP-PH/9902408'
          NNN2=NNN2+1
        ENDIF
      endif
*
      if(ICUT.eq.3) THEN
        IFUNFIN = 1
        IRFAST(INDF)=0
        IF(NNN3.EQ.1) THEN
          PRINT *,'USE OF PACKAGE acol.f'
          PRINT *,'by P. CHRISTOVA, M. JACK, T. RIEMANN,' 
          PRINT *,
     +    'DESY 99-015, hep-ph/9902408 AND PUBLICATION IN PREPARATION'
          NNN3=NNN3+1
        ENDIF
      endif
* 
      ACOLIN(INDF) = ACOL
      EF_MIN(INDF) = EMIN
      SPRIME(INDF) = S_PR
      ANGMIN(INDF) = ANG0
      ANGMAX(INDF) = ANG1
      SPRIPP(INDF) = SIPP
      RETURN
*                                                            END ZUCUTS
      END
      
      SUBROUTINE ZUTHSM(INDF,SQRS,ZMASS,TMASS,HMASS,DAL5H,ALFAS,XS,AFB)
*     ========== ======================================================
*     ZUTHSM-WRAPPER
*     ==============
      IMPLICIT REAL*8(A-H,O-Z)
      V_TB=1D0
      CALL  ZVTHSM(INDF,SQRS,ZMASS,TMASS,HMASS,DAL5H,V_TB,ALFAS,XS,AFB)
      RETURN
      END

      SUBROUTINE 
     &      ZVTHSM(INDF,SQRS,ZMASS,TMASS,HMASS,DAL5H,V_TB,ALFAS,XS,AFB)
*     ==================================================================
************************************************************************
*                                                                      *
*     ZVTHSM(INDF,SQRS,ZMASS,TMASS,HMASS,DAL5H,V_TB,ALFAS,XS*,AFB*)    *
*     ZUTHSM(INDF,SQRS,ZMASS,TMASS,HMASS,DAL5H,     ALFAS,XS*,AFB*)    *
*                                                                      *
*     Routine returns theoretical standard model cross section and     *
*     forward-backward asymmetry.                                      *
*                                                                      *
*     INDF  (int/read)   = fermion index                               *
*     SQRS  (real/read)  = sqrt(s)                                     *
*     ZMASS (real/read)  = Z0 mass (GeV)                               *
*     TMASS (real/read)  = top mass (GeV)                              *
*     HMASS (real/read)  = Higgs mass (GeV)                            *
*     DAL5H (real/read)  = hadronic 5-flavour contribution to DALPHA   *
*     V_TB  (real/read) = the Vtb mixing matrix element                *
*     ALFAS (real/read)  = \alpha_s(s)                                 *
*     XS    (real/write) = theoretical cross section                   *
*     AFB   (real/write) = theoretical asymmetry                       *
*                                                                      *
*     Called by USER                                                   *
*                                                                      *
*     NOTE:  If the ZFITTER flag 'BORN' is set to 1 XS, AFB return     *
*            improved Born values.  Forward-backward asymmetry is      *
*            only calculated if the AFBC flag has been set to 1.       *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8(A-H,O-Z)
*
* flags
*
      PARAMETER(NFLGMX=46)
      COMMON /ZUFLGS/ IFLAGS(NFLGMX),CFLAGS(NFLGMX)
      CHARACTER CFLAGS*4
      PARAMETER (IFAFBC= 1,IFSCAL= 2,IFSCRE= 3,IFAMT4= 4,IFBORN= 5,
     & IFBOXD= 6,IFCONV= 7,IFFINR= 8,IFFOT2= 9,IFGAMS=10,IFDIAG=11,
     & IFINTF=12,IFBARB=13,IFPART=14,IFPOWR=15,IFPRNT=16,IFALEM=17,
     & IFQCDC=18,IFVPOL=19,IFWEAK=20,IFFTJR=21,IFEXPR=22,IFEXPF=23,
     & IFHIGS=24,IFAFMT=25,IFCZAK=26,IFPREC=27,IFHIG2=28,IFALE2=29,
     & IFGFER=30,IFISPP=31,IFFSRS=32,IFMISC=33,IFMISD=34,IFIPFC=35,
     & IFIPSC=36,IFIPTO=37,IFFBHO=38,IFFSPP=39,IFFUNA=40,IFASCR=41,
     & IFSFSR=42,IFENUE=43,IFTUPV=44,IFDMWW=45,IFDSWW=46)
*
* parameters
*
      COMMON /ZUPARS/QDF,QCDCOR(0:14),ALPHST,SIN2TW,S2TEFF(0:11),
     & WIDTHS(0:11)
*
* cuts
*
      COMMON /ZUDATA/ ACOLIN(0:11),EF_MIN(0:11),SPRIME(0:11),
     + ANGMIN(0:11),ANGMAX(0:11),SPRIPP(0:11),IRCUTS(0:11),IRFAST(0:11)
*
* commons for ZFITTER internal use
*
      COMMON/ZFCHMS/ALLCH(0:11),ALLMS(0:11)
      COMMON/FRINIT/NPAR(30),ZPAR(31)
      COMMON/PARTZW/PARTZ(0:11),PARTW(3)
      DIMENSION QCDCCR(0:14)
      DIMENSION INDQ(10)
      DATA INDQ /0,0,0,0,1,2,3,4,5,6/
*
*-----------------------------------------------------------------------
*
      IF(INDF.LT.-1 .OR. INDF.GT.11) THEN
        PRINT *,'ZUTHSM> fermion index out of range:  INDF = ',INDF
        RETURN
      ENDIF
      XS     = 0.D0
      AFB    = 0.D0
      IF(INDF.EQ.8) RETURN
*
*-----------------------------------------------------------------------
*
      IF(SQRS.LT.9D0.OR.SQRS.GT.350D0) THEN
        PRINT*,'REQUIRED SQRS  ',SQRS,'  IS OUT OF ALLOWED RANGE'
        STOP
      ENDIF
*
* INTERF. 1
*
      INTRF = 1
*
* call DIZET
*
      CALL ZDIZET(ZMASS,TMASS,HMASS,DAL5H,V_TB,ALFAS,ALQED,ALFAT)
*
* CALCULATION OF QCD - FACTORS AS FUNCTIONS OF S
*

      S = SQRS*SQRS
      ALFAST=ALFAS
*
* call QCDCOF (ALQED=ALQEDZ)
*
      CALL ZQCDCO(SQRS,TMASS,SIN2TW,ALQED,ALFAST,ALFATT,ALPHXI,QCDCCR)
*
      GAMZ     = PARTZ(11)/1D3
      IFAST    = IRFAST(INDF)
      NPAR(11) = IRCUTS(INDF)
      NPAR(13) = 1
      DO IQCDC=0,14
        ZPAR(7+IQCDC) = QCDCCR(IQCDC)
      ENDDO
      ZPAR(26) = 1.D0-SPRIME(INDF)/S
      ZPAR(27) = ACOLIN(INDF)
      ZPAR(28) = EF_MIN(INDF)
      ZPAR(29) = ANGMAX(INDF)
      ZPAR(30) = ANGMIN(INDF)
      ZPAR(31) = SPRIPP(INDF)
*
      IMN = INDF
      IMX = INDF
      IF(INDF.EQ.10) THEN
        IMN = 4
        IMX = 9
      ENDIF
*
      DO 10 I = IMN,IMX
        IF(I.EQ.8) GOTO 10
        ZPAR(2) = ALLCH(I)
        ZPAR(4) = ALLMS(I)
        IF(INDF.NE.11) THEN
          CALL ZCUT(INTRF,IFAST,I,S,ZMASS,GAMZ,WIDTHS,SIN2TW,NPAR,ZPAR,
     &              SIGBRN,SIGQED,AFBBRN,AFBQED)
        ELSE
          CALL BHANG(INTRF,       S,ZMASS,GAMZ,WIDTHS,SIN2TW,NPAR,ZPAR,
     &               SIGBRN,SIGQED,AFBBRN,AFBQED)
        ENDIF
        SBORN = SIGBRN
        STOT  = SIGQED
        ABORN = AFBBRN
        ATOT  = AFBQED
*
        IF(IFLAGS(IFBORN).EQ.1) THEN
          XS  = XS+SBORN
          AFB = AFB+ABORN
        ELSE
          XS  = XS+STOT
          AFB = AFB+ATOT
        ENDIF
   10 ENDDO
      RETURN
*                                                             END ZUTHSM
      END

      SUBROUTINE ZUATSM(INDF,SQRS,ZMASS,TMASS,HMASS,DAL5H,ALFAS,CSA,DXS) 
*     ========== =======================================================
*     ZUATSM-WRAPPER
*     ==============
      IMPLICIT REAL*8(A-H,O-Z)
      V_TB=1D0
      CALL  ZVATSM(INDF,SQRS,ZMASS,TMASS,HMASS,DAL5H,V_TB,ALFAS,CSA,DXS)
      RETURN
      END

      SUBROUTINE 
     &      ZVATSM(INDF,SQRS,ZMASS,TMASS,HMASS,DAL5H,V_TB,ALFAS,CSA,DXS)
*     ==================================================================
************************************************************************
*                                                                      *
*     ZVATSM(INDF,SQRS,ZMASS,TMASS,HMASS,DAL5H,V_TB,ALFAS,CSA,DXS*)    *
*     ZUATSM(INDF,SQRS,ZMASS,TMASS,HMASS,DAL5H,     ALFAS,CSA,DXS*)    *
*                                                                      *
*     Routine returns theoretical standard model cross section         *
*     differential in the Cosine of the Scattering Angle, CSA.         *
*                                                                      *
*     INDF  (int/read)   = fermion index                               *
*     SQRS  (real/read)  = sqrt(s)                                     *
*     ZMASS (real/read)  = Z0 mass (GeV)                               *
*     TMASS (real/read)  = top mass (GeV)                              *
*     HMASS (real/read)  = Higgs mass (GeV)                            *
*     DAL5H (real/read)  = hadronic 5-flavour contribution to DALPHA   *
*     V_TB  (real/read) = the Vtb mixing matrix element                *
*     ALFAS (real/read)  = \alpha_s(s)                                 *
*     CSA   (real/read)  = cosine of scattering angle                  *    
*     DXS   (real/write) = theoretical differential cross section      *
*                                                                      *
*     Called by USER                                                   *
*                                                                      *
*     NOTE:  If the ZFITTER flag 'BORN' is set to 1 XS returns         *
*            improved Born value.                                      *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8(A-H,O-Z)
*
* flags
*
      PARAMETER(NFLGMX=46)
      COMMON /ZUFLGS/ IFLAGS(NFLGMX),CFLAGS(NFLGMX)
      CHARACTER CFLAGS*4
      PARAMETER (IFAFBC= 1,IFSCAL= 2,IFSCRE= 3,IFAMT4= 4,IFBORN= 5,
     & IFBOXD= 6,IFCONV= 7,IFFINR= 8,IFFOT2= 9,IFGAMS=10,IFDIAG=11,
     & IFINTF=12,IFBARB=13,IFPART=14,IFPOWR=15,IFPRNT=16,IFALEM=17,
     & IFQCDC=18,IFVPOL=19,IFWEAK=20,IFFTJR=21,IFEXPR=22,IFEXPF=23,
     & IFHIGS=24,IFAFMT=25,IFCZAK=26,IFPREC=27,IFHIG2=28,IFALE2=29,
     & IFGFER=30,IFISPP=31,IFFSRS=32,IFMISC=33,IFMISD=34,IFIPFC=35,
     & IFIPSC=36,IFIPTO=37,IFFBHO=38,IFFSPP=39,IFFUNA=40,IFASCR=41,
     & IFSFSR=42,IFENUE=43,IFTUPV=44,IFDMWW=45,IFDSWW=46)
*
* parameters
*
      COMMON /ZUPARS/QDF,QCDCOR(0:14),ALPHST,SIN2TW,S2TEFF(0:11),
     & WIDTHS(0:11)
*
* cuts
*
      COMMON /ZUDATA/ ACOLIN(0:11),EF_MIN(0:11),SPRIME(0:11),
     + ANGMIN(0:11),ANGMAX(0:11),SPRIPP(0:11),IRCUTS(0:11),IRFAST(0:11)
*
* commons for ZFITTER internal use
*
      COMMON/ZFCHMS/ALLCH(0:11),ALLMS(0:11)
      COMMON/FRINIT/NPAR(30),ZPAR(31)
      COMMON/PARTZW/PARTZ(0:11),PARTW(3)
      DIMENSION QCDCCR(0:14)
      DIMENSION INDQ(10)
      DATA INDQ /0,0,0,0,1,2,3,4,5,6/
*
*-----------------------------------------------------------------------
*
      IF(INDF.LT.-1 .OR. INDF.GT.11) THEN
        PRINT *,'ZUATSM> fermion index out of range:  INDF = ',INDF
        RETURN
      ENDIF
*
***   IF(IFLAGS(IFFOT2).EQ.5) THEN
***     PRINT *,'ZUATSM> can not work in CF3 mode! Program STOPPED'
***     STOP
***   ENDIF
*
      IF(IFLAGS(IFISPP).EQ.-1.OR.IFLAGS(IFISPP).EQ.1) THEN
        PRINT *,'ZUATSM> can not work with ISPP=-1,+1! Program STOPPED'
        STOP
      ENDIF
*
      IF(IRCUTS(INDF).EQ.0) THEN
        PRINT *,'ZUATSM> can not work with ICUT=0! Program STOPPED'
        STOP
      ENDIF
*
      DXS = 0.D0
      IF(INDF.EQ.8) RETURN
*
*-----------------------------------------------------------------------
*
      IF(SQRS.LT.9.5D0.OR.SQRS.GT.350D0) THEN
        PRINT*,'REQUIRED SQRS  ',SQRS,'  IS OUT OF ALLOWED RANGE'
        STOP
      ENDIF
*
* INTERF. 1
*
      INTRF = 1
*
* call DIZET
*
      CALL ZDIZET(ZMASS,TMASS,HMASS,DAL5H,V_TB,ALFAS,ALQED,ALFAT)
*
* CALCULATION OF QCD - FACTORS AS FUNCTIONS OF S
*
      S = SQRS*SQRS
      ALFAST=ALFAS
*
* call QCDCOF (ALQED=ALQEDZ)
*
      CALL ZQCDCO(SQRS,TMASS,SIN2TW,ALQED,ALFAST,ALFATT,ALPHXI,QCDCCR)
*
      GAMZ     = PARTZ(11)/1D3
      IFAST    = IRFAST(INDF)
      NPAR(11) = IRCUTS(INDF)
      NPAR(13) = 1
      DO IQCDC=0,14
        ZPAR(7+IQCDC) = QCDCCR(IQCDC)
      ENDDO
      ZPAR(26) = 1.D0-SPRIME(INDF)/S
      ZPAR(27) = ACOLIN(INDF)
      ZPAR(28) = EF_MIN(INDF)
      ZPAR(29) = ANGMAX(INDF)
      ZPAR(30) = ANGMIN(INDF)
      ZPAR(31) = SPRIPP(INDF)
*
      IMN = INDF
      IMX = INDF
      IF(INDF.EQ.10) THEN
        IMN = 4
        IMX = 9
      ENDIF
      DO 10 I = IMN,IMX
        IF(I.EQ.8) GOTO 10
        ZPAR(2) = ALLCH(I)
        ZPAR(4) = ALLMS(I)
        IF(INDF.NE.11) THEN
          CALL ZANCUT(INTRF,IFAST,I,CSA,
     +              S,ZMASS,GAMZ,WIDTHS,SIN2TW,NPAR,ZPAR,SIGBRN,SIGQED)
        ELSE
* CALL to BHANG is not correct
          PRINT *,'INDF=11 does not exist in this branch'
cbard     CALL BHANG(INTRF,       S,ZMASS,GAMZ,WIDTHS,SIN2TW,NPAR,ZPAR,
cbard+     SIGBRN,SIGQED,AFBBRN,AFBQED)
          stop
        ENDIF
        SBORN = SIGBRN
        STOT  = SIGQED
*
        IF(IFLAGS(IFBORN).EQ.1) THEN
         DXS = DXS+SBORN
          ELSE
         DXS = DXS+STOT
        ENDIF
   10 ENDDO
      RETURN
*                                                             END ZUATSM
      END
     
      SUBROUTINE 
     &          ZUTPSM(SQRS,ZMASS,TMASS,HMASS,DAL5H,ALFAS,TAUPOL,TAUAFB)
*     ==================================================================
*     ZUTPSM-WRAPPER
*     ==============
      IMPLICIT REAL*8(A-H,O-Z)
      V_TB=1D0
      CALL ZVTPSM(SQRS,ZMASS,TMASS,HMASS,DAL5H,V_TB,ALFAS,TAUPOL,TAUAFB)
      RETURN
      END

      SUBROUTINE
     &     ZVTPSM(SQRS,ZMASS,TMASS,HMASS,DAL5H,V_TB,ALFAS,TAUPOL,TAUAFB)
*     ==================================================================
************************************************************************
*                                                                      *
*     ZVTPSM(SQRS,ZMASS,TMASS,HMASS,DAL5H,V_TB,ALFAS,TAUPOL*,TAUAFB*)  *
*     ZUTPSM(SQRS,ZMASS,TMASS,HMASS,DAL5H,     ALFAS,TAUPOL*,TAUAFB*)  *
*                                                                      *
*     Routine returns theoretical standard model tau polarization and  *
*     tau polarization asymmetry.                                      *
*                                                                      *
*     SQRS  (real/read)  = sqrt(s)                                     *
*     ZMASS (real/read)  = Z0 mass (GeV)                               *
*     TMASS (real/read)  = top mass (GeV)                              *
*     HMASS (real/read)  = Higgs mass (GeV)                            *
*     DAL5H (real/read)  = hadronic 5-flavour contribution to DALPHA   *
*     V_TB  (real/read) = the Vtb mixing matrix element                *
*     ALFAS (real/read)  = \alpha_s(s)                                 *
*     TAUPOL(real/write) = theoretical tau polarization                *
*     TAUAFB(real/write) = theoretical tau polarization asymmetry      *
*                                                                      *
*     Called by USER                                                   *
*                                                                      *
*     NOTE:  If the ZFITTER flag 'BORN' is set to 1 TAUPOL, TAUAFB     *
*            return improved Born values.                              *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8(A-H,O-Z)
*
* flags
*
      PARAMETER(NFLGMX=46)
      COMMON /ZUFLGS/ IFLAGS(NFLGMX),CFLAGS(NFLGMX)
      CHARACTER CFLAGS*4
      PARAMETER (IFAFBC= 1,IFSCAL= 2,IFSCRE= 3,IFAMT4= 4,IFBORN= 5,
     & IFBOXD= 6,IFCONV= 7,IFFINR= 8,IFFOT2= 9,IFGAMS=10,IFDIAG=11,
     & IFINTF=12,IFBARB=13,IFPART=14,IFPOWR=15,IFPRNT=16,IFALEM=17,
     & IFQCDC=18,IFVPOL=19,IFWEAK=20,IFFTJR=21,IFEXPR=22,IFEXPF=23,
     & IFHIGS=24,IFAFMT=25,IFCZAK=26,IFPREC=27,IFHIG2=28,IFALE2=29,
     & IFGFER=30,IFISPP=31,IFFSRS=32,IFMISC=33,IFMISD=34,IFIPFC=35,
     & IFIPSC=36,IFIPTO=37,IFFBHO=38,IFFSPP=39,IFFUNA=40,IFASCR=41,
     & IFSFSR=42,IFENUE=43,IFTUPV=44,IFDMWW=45,IFDSWW=46)
*
* parameters
*
      COMMON /ZUPARS/QDF,QCDCOR(0:14),ALPHST,SIN2TW,S2TEFF(0:11),
     & WIDTHS(0:11)
*
* cuts
*
      COMMON /ZUDATA/ ACOLIN(0:11),EF_MIN(0:11),SPRIME(0:11),
     + ANGMIN(0:11),ANGMAX(0:11),SPRIPP(0:11),IRCUTS(0:11),IRFAST(0:11)
*
* commons for ZFITTER internal use
*
      COMMON/ZFCHMS/ALLCH(0:11),ALLMS(0:11)
      COMMON/FRINIT/NPAR(30),ZPAR(31)
      COMMON/PARTZW/PARTZ(0:11),PARTW(3)
      DIMENSION QCDCCR(0:14)
      DIMENSION INDQ(10)
      DATA INDQ /0,0,0,0,1,2,3,4,5,6/
*
*-----------------------------------------------------------------------
*
* INTERF. 1
*
      INTRF = 1
*
* call DIZET
*
      CALL ZDIZET(ZMASS,TMASS,HMASS,DAL5H,V_TB,ALFAS,ALQED,ALFAT)
cdl
       stop
*
* CALCULATION OF QCD - FACTORS AS FUNCTIONS OF S
*
      S = SQRS*SQRS
      ALFAST=ALFAS
*
* call QCDCOF (ALQED=ALQEDZ)
*
      CALL ZQCDCO(SQRS,TMASS,SIN2TW,ALQED,ALFAST,ALFATT,ALPHXI,QCDCCR)
*
      GAMZ     = PARTZ(11)/1D3
      NPAR(11) = IRCUTS(3)
      ICURR    = NPAR(13)
      NPAR(13) = 1
      ZPAR(2)  = ALLCH(3)
      ZPAR(4)  = ALLMS(3)
      DO IQCDC=0,14
        ZPAR(7+IQCDC) = QCDCCR(IQCDC)
      ENDDO
      ZPAR(26) = 1.D0-SPRIME(3)/S
      ZPAR(27) = ACOLIN(3)
      ZPAR(28) = EF_MIN(3)
      ZPAR(29) = ANGMAX(3)
      ZPAR(30) = ANGMIN(3)
      ZPAR(31) = SPRIPP(3)
*
      ZPAR(24) = 0.D0
      ZPAR(25) = 1.D0
      CALL ZCUT(INTRF,IRFAST(3),3,S,ZMASS,GAMZ,WIDTHS,SIN2TW,NPAR,ZPAR,
     + SBORPL,SQEDPL,AFBPL,AFQPL)
      ZPAR(25) =-1.D0
      CALL ZCUT(INTRF,IRFAST(3),3,S,ZMASS,GAMZ,WIDTHS,SIN2TW,NPAR,ZPAR,
     + SBORMI,SQEDMI,AFBMI,AFQMI)
*
      IF(IFLAGS(IFBORN).EQ.1) THEN
        TAUPOL = (SBORPL-SBORMI)/(SBORPL+SBORMI)
        TAUAFB = (AFBPL*SBORPL - AFBMI*SBORMI)/(SBORPL+SBORMI)
      ELSE
        TAUPOL = (SQEDPL-SQEDMI)/(SQEDPL+SQEDMI)
        TAUAFB = (AFQPL*SQEDPL - AFQMI*SQEDMI)/(SQEDPL+SQEDMI)
      ENDIF
*
      ZPAR(25) = 0.D0
      NPAR(13) = ICURR
      RETURN
*                                                             END ZUTPSM
      END
 
      SUBROUTINE
     &     ZULRSM(INDF,SQRS,ZMASS,TMASS,HMASS,DAL5H,ALFAS,POL,XSPL,XSMI)
*     ==================================================================
*     ZULRSM-WRAPPER
*     ==============
      IMPLICIT REAL*8(A-H,O-Z)
      V_TB=1D0
      CALL 
     &ZVLRSM(INDF,SQRS,ZMASS,TMASS,HMASS,DAL5H,V_TB,ALFAS,POL,XSPL,XSMI)
      RETURN
      END
      
      SUBROUTINE
     &ZVLRSM(INDF,SQRS,ZMASS,TMASS,HMASS,DAL5H,V_TB,ALFAS,POL,XSPL,XSMI)
*     ==================================================================
************************************************************************
*                                                                      *
*   ZVLRSM(INDF,SQRS,ZMASS,TMASS,HMASS,DAL5H,V_TB,ALFAS,POL,XSPL,XSMI) *
*   ZULRSM(INDF,SQRS,ZMASS,TMASS,HMASS,DAL5H,     ALFAS,POL,XSPL,XSMI) *
*                                                                      *
*     Routine returns theoretical standard model left-right            *
*     polarization asymmetry.                                          *
*                                                                      *
*     SQRS  (real/read)  = sqrt(s)                                     *
*     ZMASS (real/read)  = Z0 mass (GeV)                               *
*     TMASS (real/read)  = top mass (GeV)                              *
*     HMASS (real/read)  = Higgs mass (GeV)                            *
*     DAL5H (real/read)  = hadronic 5-flavour contribution to DALPHA   *
*     V_TB  (real/read) = the Vtb mixing matrix element                *
*     ALFAS (real/read)  = \alpha_s(s)                                 *
*     POL   (real/read)  = degree of the longitudinal e^- polarization *
*     XSPL  (real/write) = theoretical cross section for a given +POL  *
*     XSMI  (real/write) = theoretical cross section for a given -POL  *
*     Called by USER                                                   *
*                                                                      *
*     NOTE:  If the ZFITTER flag 'BORN' is set to 1 ALR                *
*            return improved Born values.                              *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8(A-H,O-Z)
*
* flags
*
      PARAMETER(NFLGMX=46)
      COMMON /ZUFLGS/ IFLAGS(NFLGMX),CFLAGS(NFLGMX)
      CHARACTER CFLAGS*4
      PARAMETER (IFAFBC= 1,IFSCAL= 2,IFSCRE= 3,IFAMT4= 4,IFBORN= 5,
     & IFBOXD= 6,IFCONV= 7,IFFINR= 8,IFFOT2= 9,IFGAMS=10,IFDIAG=11,
     & IFINTF=12,IFBARB=13,IFPART=14,IFPOWR=15,IFPRNT=16,IFALEM=17,
     & IFQCDC=18,IFVPOL=19,IFWEAK=20,IFFTJR=21,IFEXPR=22,IFEXPF=23,
     & IFHIGS=24,IFAFMT=25,IFCZAK=26,IFPREC=27,IFHIG2=28,IFALE2=29,
     & IFGFER=30,IFISPP=31,IFFSRS=32,IFMISC=33,IFMISD=34,IFIPFC=35,
     & IFIPSC=36,IFIPTO=37,IFFBHO=38,IFFSPP=39,IFFUNA=40,IFASCR=41,
     & IFSFSR=42,IFENUE=43,IFTUPV=44,IFDMWW=45,IFDSWW=46)
*
* parameters
*
      COMMON /ZUPARS/QDF,QCDCOR(0:14),ALPHST,SIN2TW,S2TEFF(0:11),
     & WIDTHS(0:11)
*
* cuts
*
      COMMON /ZUDATA/ ACOLIN(0:11),EF_MIN(0:11),SPRIME(0:11),
     + ANGMIN(0:11),ANGMAX(0:11),SPRIPP(0:11),IRCUTS(0:11),IRFAST(0:11)
*
* commons for ZFITTER internal use
*
      COMMON/ZFCHMS/ALLCH(0:11),ALLMS(0:11)
      COMMON/FRINIT/NPAR(30),ZPAR(31)
      COMMON/PARTZW/PARTZ(0:11),PARTW(3)
      DIMENSION QCDCCR(0:14)
      DIMENSION INDQ(10)
      DATA INDQ /0,0,0,0,1,2,3,4,5,6/
*
*-----------------------------------------------------------------------
*
      IF(INDF.LT.-1 .OR. INDF.GT.11) THEN
        PRINT *,'ZULRSM> fermion index out of range:  INDF = ',INDF
        RETURN
      ENDIF
      XS     = 0.D0
      AFB    = 0.D0
      IF(INDF.EQ.8) RETURN
*
*-----------------------------------------------------------------------
*
      IF(SQRS.LT.9D0.OR.SQRS.GT.350D0) THEN
        PRINT*,'REQUIRED SQRS  ',SQRS,'  IS OUT OF ALLOWED RANGE'
        STOP
      ENDIF
*
* INTERF. 1
*
      INTRF = 1
*
* call DIZET
*
      CALL ZDIZET(ZMASS,TMASS,HMASS,DAL5H,V_TB,ALFAS,ALQED,ALFAT)
*
* CALCULATION OF QCD - FACTORS AS FUNCTIONS OF S
*
      S = SQRS*SQRS
      ALFAST=ALFAS
*
* call QCDCOF (ALQED=ALQEDZ)
*
      CALL ZQCDCO(SQRS,TMASS,SIN2TW,ALQED,ALFAST,ALFATT,ALPHXI,QCDCCR)
*
      GAMZ     = PARTZ(11)/1D3
      IFAST    = IRFAST(INDF)
      NPAR(11) = IRCUTS(INDF)
      ICURR    = NPAR(13)
      NPAR(13) = 1
      DO IQCDC=0,14
        ZPAR(7+IQCDC) = QCDCCR(IQCDC)
      ENDDO
      ZPAR(26) = 1.D0-SPRIME(INDF)/S
      ZPAR(27) = ACOLIN(INDF)
      ZPAR(28) = EF_MIN(INDF)
      ZPAR(29) = ANGMAX(INDF)
      ZPAR(30) = ANGMIN(INDF)
      ZPAR(31) = SPRIPP(INDF)
*
      ZPAR(22) = 0D0
      ZPAR(23) = POL
      I=INDF
      ZPAR(2) = ALLCH(I)
      ZPAR(4) = ALLMS(I)
      CALL ZCUT(INTRF,IFAST,I,S,ZMASS,GAMZ,WIDTHS,SIN2TW,NPAR,ZPAR,
     + SBORPL,SQEDPL,AFBPL,AFQPL)
      ZPAR(23) =-POL
      CALL ZCUT(INTRF,IFAST,I,S,ZMASS,GAMZ,WIDTHS,SIN2TW,NPAR,ZPAR,
     + SBORMI,SQEDMI,AFBMI,AFQMI)
*
      IF(IFLAGS(IFBORN).EQ.1) THEN
        XSPL= SBORPL
        XSMI= SBORMI
      ELSE
        XSPL= SQEDPL
        XSMI= SQEDMI
      ENDIF
*
      ZPAR(23) = 0D0
      NPAR(13) = ICURR
      RETURN
*                                                             END ZULRSM
      END

      SUBROUTINE 
     &      ZU_APV(ZMASS,TMASS,HMASS,SIN2TW,UMASS,DMASS,C1U,C1D,C2U,C2D)
*     ==================================================================
************************************************************************
*                                                                      *
*     ZU_APV(ZMASS,TMASS,HMASS,SIN2TW,UMASS,DMASS,C1U,C1D,C2U,C2D)     *
*                                                                      *
*     Routine returns theoretical standard model predictions for APV   *
*     parameters: C1U,C1D,C2U,C2D (signs convention as in PDG)         *
*     OMS-consistent calculation by                                    *
*     D.Bardin, P.Christova and L.Kalinovskaya, January 2001           *
*     following W.J.Marciano & A.Sirlin, PRD 27 (1983) 552-556         *
*                                                                      *
*     ZMASS (real/read)  = Z0    mass (GeV)                            *
*     TMASS (real/read)  = top   mass (GeV)                            *
*     HMASS (real/read)  = Higgs mass (GeV)                            *
*     SIN2TW(real/read)  = \sin of weak mixing angle                   *
*     UMASS (real/read)  = u-quark mass (constituent)                  *
*     DMASS (real/read)  = d-quark mass (constituent)                  *
*                                                                      *
*     C1U(D)(real/write) = theoretical APV parameters for QW           *
*     C2U(D)(real/write) = theoretical APV C2's parameters             *
*                                                                      *
*     Called by USER                                                   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8(A-H,O-W,Y,Z)
      IMPLICIT COMPLEX*16(X)
      INTEGER*4 EWFFTR,GAMZTR,FERMTR
      COMPLEX*16 CV1,CA1,CF1
*
      PARAMETER(ALFAI=137.0359895D0,EMASS=.51099907D-3)!,GMU=1.16639D-5)
*
* flags
*
      PARAMETER(NFLGMX=46)
      COMMON /ZUFLGS/ IFLAGS(NFLGMX),CFLAGS(NFLGMX)
      COMMON /CDZZWG/ AMZ,AMH,GMU,A0,GAMZ,GAMW,ALSZ,ALST,CALXI,CALQED
      CHARACTER CFLAGS*4
      PARAMETER (IFAFBC= 1,IFSCAL= 2,IFSCRE= 3,IFAMT4= 4,IFBORN= 5,
     & IFBOXD= 6,IFCONV= 7,IFFINR= 8,IFFOT2= 9,IFGAMS=10,IFDIAG=11,
     & IFINTF=12,IFBARB=13,IFPART=14,IFPOWR=15,IFPRNT=16,IFALEM=17,
     & IFQCDC=18,IFVPOL=19,IFWEAK=20,IFFTJR=21,IFEXPR=22,IFEXPF=23,
     & IFHIGS=24,IFAFMT=25,IFCZAK=26,IFPREC=27,IFHIG2=28,IFALE2=29,
     & IFGFER=30,IFISPP=31,IFFSRS=32,IFMISC=33,IFMISD=34,IFIPFC=35,
     & IFIPSC=36,IFIPTO=37,IFFBHO=38,IFFSPP=39,IFFUNA=40,IFASCR=41,
     & IFSFSR=42,IFENUE=43,IFTUPV=44,IFDMWW=45,IFDSWW=46)
*
      COMMON/TREATM/EWFFTR,GAMZTR,FERMTR
      COMMON/CDZDEG/DROBAR,DROBLO,DRREMD
      COMMON/CDZAPV/DRHOD,DKDREM,DRDREM
*
* APV TUs
*
      ITUPV=IFLAGS(IFTUPV)
      IF(ITUPV.EQ.1) THEN
        TUAPV1=0D0
        TUAPV2=0D0
      ELSEIF(ITUPV.EQ.2) THEN
        TUAPV1=1D0
        TUAPV2=1D0
      ELSEIF(ITUPV.EQ.3) THEN
        TUAPV1=1D0
        TUAPV2=0D0
      ENDIF
*
* APV constants
*

      FERMTR=3
      PI=4D0*ATAN(1D0)
      QE=-1D0
      QU=+2D0/3
      QD=-1D0/3
*
      COS2TW=1D0-SIN2TW
      ZM2=ZMASS**2
      TM2=TMASS**2
      HM2=HMASS**2
      WM2=ZM2*COS2TW
      WMASS=SQRT(WM2)
      HS2=WM2
*
      CALL INEETT(HS2,TMASS,WMASS,ZMASS,HMASS)
      CALL ROKAP0(HS2,ZM2,TM2,HM2,WM2,DRHO0,DKAP0,DR0LEA,DROLEA,DKALEA)
*
      DROREM=DRHO0-DR0LEA
      DKAREM=DKAP0-DKALEA
*
      CW2=COS2TW
      SW2=SIN2TW
      QTM=2D0/3
      QBM=1D0/3
      XCR=DCMPLX(ZM2/(4D0*TM2))
      XCX=DCMPLX(WM2/TM2)
      XROQCD=1D0/SW2*TM2/WM2*(
     &      -1D0/4D0*((1D0-4D0*QTM*SW2)**2*CV1(XCR)+CA1(XCR))+CF1(XCX)
     &      +1D0/8D0*ZM2/TM2*(1D0-2D0*QBM*SW2)*LOG(ZM2/TM2))
*
      DROQCD=1D0/ALFAI/PI*ALST/PI*(DREAL(XROQCD)-1D0/2D0/SW2*LOG(CW2))
*
      CONVF=SQRT(2D0)*GMU*ZM2*SIN2TW*COS2TW/PI*ALFAI
      SINMS=(1D0+(-1D0/4/ALFAI/PI*CW2/SW2**2*DROLEA*CONVF+CW2/SW2*DRHOD)
     &               *(1D0-DRREMD)-CW2/SW2*(DROQCD+TBQCD3)*CONVF)*SW2
*     print *,'SINMSB=',SINMS
*
      TBQCD0=AFMT3(ALST,TM2,ZM2,SW2)
* Below is pure AFMT-correction to \Delta\rho, but not for Dktb yet!!!
      TBQCD3=-.75D0/4/ALFAI/PI/SW2*TM2/WM2
     &          *(TBQCD0-(-ALST/PI*2D0/3D0*(1D0+PI**2/3D0)))     
*
      RHOPV= 1D0+1D0/4/ALFAI/PI/SW2*(DR0LEA+DROREM)+DROQCD+TBQCD3
     &      -DRHOD*TUAPV1
      AKAPV=(1D0+1D0/4/ALFAI/PI/SW2*DKAREM+DKDREM*TUAPV2)
     &     *(1D0+(-1D0/4/ALFAI/PI*CW2/SW2**2*DROLEA*CONVF+CW2/SW2*DRHOD)
     &               *(1D0-DRREMD)-CW2/SW2*(DROQCD+TBQCD3)*CONVF)
*
      VE=1D0-4D0*SIN2TW
      VU=1D0-8D0/3*SIN2TW
      VD=1D0-4D0/3*SIN2TW
*
      C1U=-(+.5D0*RHOPV*(1D0-8D0/3*AKAPV*SIN2TW)
     &   +1D0/ALFAI/PI*(-1D0/4*VU
     &       -1D0/9*VE*(LOG(ZMASS**2/EMASS**2)+1D0/6)
     &       +1D0/2/SIN2TW
     &       +1D0/2*VE*(LOG(ZMASS**2/UMASS**2)+3D0/2)
     &       +3D0/64/SIN2TW/COS2TW*VU*(1D0+VE**2)
     &                 )
     &     )
*
      C1D=-(-.5D0*RHOPV*(1D0-4D0/3*AKAPV*SIN2TW)
     &   +1D0/ALFAI/PI*(+1D0/4*VD
     &       +1D0/18*VE*(LOG(ZMASS**2/EMASS**2)+1D0/6)
     &       -1D0/8/SIN2TW
     &       +1D0/4*VE*(LOG(ZMASS**2/DMASS**2)+3D0/2)
     &       +3D0/64/SIN2TW/COS2TW*VD*(1D0+VE**2)
     &                 )
     &     )
*
      C2U=-(+.5D0*RHOPV*(1D0-4D0*AKAPV*SIN2TW)
     &   +1D0/ALFAI/PI*(-1D0/9*VE
     &       +1D0/9*(LOG(WM2/UMASS**2)+1D0/6)
     &       -1D0/9*VU*(LOG(ZMASS**2/UMASS**2)+1D0/6)
     &       +1D0/2/SIN2TW
     &       +1D0/2*VU*(LOG(ZMASS**2/UMASS**2)+3D0/2)
     &       +3D0/64/SIN2TW/COS2TW*VE*(1D0+VU**2)
     &                 )
     &     )
*
      C2D=-(-.5D0*RHOPV*(1D0-4D0*AKAPV*SIN2TW)
     &   +1D0/ALFAI/PI*(+1D0/36*VE
     &       -2D0/9*(LOG(WM2/DMASS**2)+1D0/6)
     &       +1D0/18*VD*(LOG(ZMASS**2/DMASS**2)+1D0/6)
     &       -1D0/8/SIN2TW
     &       +1D0/4*VD*(LOG(ZMASS**2/DMASS**2)+3D0/2)
     &       +3D0/64/SIN2TW/COS2TW*VE*(1D0+VD**2)
     &                 )
     &     )
*
      RETURN
*                                                             END ZU_APV
      END

      SUBROUTINE ZDIZET(ZMASS,TMASS,HMASS,DAL5H,V_TB,ALFAS,ALQED,ALFAT)
*     =================================================================
************************************************************************
*                                                                      *
*     SUBR. ZDIZET(ZMASS,TMASS,HMASS,DAL5H,V_TB,ALFAS,ALQED,ALFAT)     *
*                                                                      *
*     CPU saving interface to DIZET.                                   *
*                                                                      *
*     ZMASS (real/read) = Z0 mass (GeV)                                *
*     TMASS (real/read) = top mass (GeV)                               *
*     HMASS (real/read) = Higgs mass (GeV)                             *
*     DAL5H (real/read) = Hadronic vacuum polarization at M(Z)**2      *
*     ALFAS (real/read) = the strong coupling constant (Q**2 = M(Z)**2)*
*     V_TB  (real/read) = the V_TB mixing matrix element               *
*     ALQED (real/outp)=  running \alpha_em(M(Z))                      *
*     ALFAT (real/outp) = the strong coupling constant (Q**2 = M(T)**2)*
*                                                                      *
*     Called by ZUTHSM,ZUTPSM,ZUATSM                                   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(ALFAI=137.0359895D0)
*
* flags
*
      PARAMETER(NFLGMX=46)
      COMMON /ZUFLGS/ IFLAGS(NFLGMX),CFLAGS(NFLGMX)
      CHARACTER CFLAGS*4
      PARAMETER (IFAFBC= 1,IFSCAL= 2,IFSCRE= 3,IFAMT4= 4,IFBORN= 5,
     & IFBOXD= 6,IFCONV= 7,IFFINR= 8,IFFOT2= 9,IFGAMS=10,IFDIAG=11,
     & IFINTF=12,IFBARB=13,IFPART=14,IFPOWR=15,IFPRNT=16,IFALEM=17,
     & IFQCDC=18,IFVPOL=19,IFWEAK=20,IFFTJR=21,IFEXPR=22,IFEXPF=23,
     & IFHIGS=24,IFAFMT=25,IFCZAK=26,IFPREC=27,IFHIG2=28,IFALE2=29,
     & IFGFER=30,IFISPP=31,IFFSRS=32,IFMISC=33,IFMISD=34,IFIPFC=35,
     & IFIPSC=36,IFIPTO=37,IFFBHO=38,IFFSPP=39,IFFUNA=40,IFASCR=41,
     & IFSFSR=42,IFENUE=43,IFTUPV=44,IFDMWW=45,IFDSWW=46)
*
      COMMON /EXPERT/ IMASK,IMASS
*
* parameters
*
      COMMON /ZUPARS/QDF,QCDCOR(0:14),ALPHST,SIN2TW,S2TEFF(0:11),
     & WIDTHS(0:11)
*
* commons for ZFITTER internal use
*
      COMMON/ZFCHMS/ALLCH(0:11),ALLMS(0:11)
      COMMON/FRINIT/NPAR(30),ZPAR(31)
      COMMON/PARTZW/PARTZ(0:11),PARTW(3)
      COMMON/CALQED/ALQEDZ,ALQEDS
      COMMON/CDAL5H/CDAL5H
*      COMMON/ZPARD_/ZPARD
*     
      DIMENSION NPARD(25),NPARD0(25),ZPARD(30)
*
      SAVE NPARD0,ZMASS0,TMASS0,HMASS0,DAL5H0,ALFAS0,V_TB0
*
      LOGICAL*4 FIRST
      SAVE FIRST
      DATA FIRST/.TRUE./
*
*-----------------------------------------------------------------------
*
      NPARD(1) = IFLAGS(IFVPOL)
cdl      NPARD(2) = IFLAGS(IFAMT4)
      NPARD(2) = 0!5!6
      NPARD(3) = 3!IFLAGS(IFQCDC)
      NPARD(4) = 3!2!3
      NPARD(5) = IMASS
      NPARD(6) = IFLAGS(IFSCRE)
      NPARD(7) = IFLAGS(IFALEM)
      NPARD(8) = IMASK
      NPARD(9) = IFLAGS(IFSCAL)
      NPARD(10)= IFLAGS(IFBARB)
      NPARD(11)= IFLAGS(IFFTJR)
      NPARD(12)= IFLAGS(IFEXPR)
      NPARD(13)= IFLAGS(IFEXPF)
      NPARD(14)= IFLAGS(IFHIGS)
      NPARD(15)= IFLAGS(IFAFMT)
      NPARD(16)= 1
      NPARD(17)= IFLAGS(IFCZAK)
      NPARD(18)= IFLAGS(IFHIG2)
      NPARD(19)= IFLAGS(IFALE2)
      NPARD(20)= IFLAGS(IFGFER)
      NPARD(21)= 1
      NPARD(22)= 0
      NPARD(23)= IFLAGS(IFSFSR)
      NPARD(24)= IFLAGS(IFDMWW)
      NPARD(25)= IFLAGS(IFDSWW)
*
* return if DIZET input parameters have not changed since last call
*
      IF(FIRST) THEN
        FIRST = .FALSE.
        GOTO 10
      ENDIF
*cbardin!!!!
      DO I = 1,25
        IF(NPARD(I).NE.NPARD0(I)) GOTO 10
      ENDDO
      IF( ZMASS.NE.ZMASS0.OR.TMASS.NE.TMASS0.OR.HMASS.NE.HMASS0
     &.OR.ALFAS.NE.ALFAS0.OR.DAL5H.NE.DAL5H0.OR.V_TB .NE.V_TB0)
     & GOTO 10
      RETURN
*
   10 CONTINUE
*
* save DIZET input parameters
*
      DO I = 1,25
        NPARD0(I) = NPARD(I)
      ENDDO
      ZMASS0 = ZMASS
      TMASS0 = TMASS
      HMASS0 = HMASS
      DAL5H0 = DAL5H
      ALFAS0 = ALFAS
      V_TB0  = V_TB
*
* call DIZET
*
      PI    = DACOS(-1.D0)
      W_MAS = 80D0
      W_MASS = 80.385d0
      W_MASS = 80.4514958d0
      Z_MAS = ZMASS
*
* protection at HMASS=2*Z_MAS
*
      AHMASS=HMASS
      if(ABS(HMASS-2d0*Z_MAS).LT.5d-5) AHMASS=2d0*Z_MAS-1d-4
*
      CALL DIZET (NPARD,W_MASS,Z_MAS,TMASS,AHMASS,DAL5H,V_TB
     &           ,ALFAS,ALQED,ALFAT,ZPARD,PARTZ,PARTW)
*
* Renovation of WIDTHS BEWARE!!!!!!!!!!!!!
*
      DO 1 I=0,11
        WIDTHS(I)=PARTZ(I)
 1    CONTINUE
*
      SIN2TW = ZPARD(3)
      ALPHST = ZPARD(15)
      DO IQCDC=0,14
      QCDCOR(IQCDC) = ZPARD(16+IQCDC)
      ENDDO
*
* ALQEDZ=ALQED(M^2_Z), ALQEDS=ALQED(S), here yet equal
*
      ALQEDZ=ALQED
      ALQEDS=ALQED
      CDAL5H=DAL5H
*
      RETURN
*                                                             END ZDIZET
      END
 
      SUBROUTINE
     &       ZQCDCO(SQRS,TMASS,SINEFF,ALQED,ALFAST,ALFATT,ALPHXI,QCDCCR)
*     ==================================================================
************************************************************************
*                                                                      *
*  SUBR.ZQCDCO(SQRS,TMASS,SIN2TW,ALQED,ALFAST*,ALFATT*,ALPHXI*,QCDCCR*)*
*                                                                      *
*     CPU saving interface to QCDCOF                                   *
*                                                                      *
*     SQRS  (real/read) = SQRT(S)  (GeV)                               *
*     TMASS (real/read) = top mass (GeV)                               *
*     SIN2TW(real/read) = a weak mixing angle                          *
*     ALFAST(real/read) = the strong coupling constant (Q**2 = M(Z)**2)*
*                                                                      *
*     Called by ZUATSM                                                 *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(ALFAI=137.0359895D0)
*
* flags
*
      PARAMETER(NFLGMX=46)
      COMMON /ZUFLGS/ IFLAGS(NFLGMX),CFLAGS(NFLGMX)
      CHARACTER CFLAGS*4
      PARAMETER (IFAFBC= 1,IFSCAL= 2,IFSCRE= 3,IFAMT4= 4,IFBORN= 5,
     & IFBOXD= 6,IFCONV= 7,IFFINR= 8,IFFOT2= 9,IFGAMS=10,IFDIAG=11,
     & IFINTF=12,IFBARB=13,IFPART=14,IFPOWR=15,IFPRNT=16,IFALEM=17,
     & IFQCDC=18,IFVPOL=19,IFWEAK=20,IFFTJR=21,IFEXPR=22,IFEXPF=23,
     & IFHIGS=24,IFAFMT=25,IFCZAK=26,IFPREC=27,IFHIG2=28,IFALE2=29,
     & IFGFER=30,IFISPP=31,IFFSRS=32,IFMISC=33,IFMISD=34,IFIPFC=35,
     & IFIPSC=36,IFIPTO=37,IFFBHO=38,IFFSPP=39,IFFUNA=40,IFASCR=41,
     & IFSFSR=42,IFENUE=43,IFTUPV=44,IFDMWW=45,IFDSWW=46)
*
      COMMON /EXPERT/ IMASK,IMASS
*
* parameters
*
      COMMON /ZUPARS/QDF,QCDCOR(0:14),ALPHST,SIN2TW,S2TEFF(0:11),
     & WIDTHS(0:11)
*
* commons for ZFITTER internal use
*
      COMMON/ZFCHMS/ALLCH(0:11),ALLMS(0:11)
      COMMON/FRINIT/NPAR(30),ZPAR(31)
      COMMON/PARTZW/PARTZ(0:11),PARTW(3)
      DIMENSION QCDCCR(0:14),QCDCC0(0:14)
      SAVE SQRS0,TMASS0,SINEF0,ALQED0,ALFAS0,QCDCC0
*
      LOGICAL*4 FIRST
      SAVE FIRST
      DATA FIRST/.TRUE./
*
* return if QCDCOF input parameters have not changed since last call
*
      IF(FIRST) THEN
        FIRST = .FALSE.
        GOTO 10
      ENDIF
      IF(SQRS.NE.SQRS0 .OR. TMASS.NE.TMASS0 .OR. SINEFF.NE.SINEF0
     + .OR. ALQED.NE.ALQED0 .OR. ALFAST.NE.ALFAS0) GOTO 10
        DO ICR=0,14
          QCDCCR(ICR)=QCDCC0(ICR)
        ENDDO
      RETURN
*
* save QCDCOF input parameters
*
   10 CONTINUE
      SQRS0  = SQRS
      TMASS0 = TMASS
      SINEF0 = SINEFF
      ALQED0 = ALQED
      ALFAS0 = ALFAST
*
* call QCDCOF
*
      CALL QCDCOF(SQRS,TMASS,SIN2TW,ALQED,ALFAST,ALFATT,ALPHXI,QCDCCR)
*
* save QCDCOF output parameters
*
      DO ICR=0,14
        QCDCC0(ICR)=QCDCCR(ICR)
      ENDDO
*
      RETURN
*                                                             END ZQCDCO
      END
 
      SUBROUTINE ZUXSEC(INDF,SQRS,ZMASS,GAMZ0,GAMEE,GAMFF,XS)
*     ========== ============================================
************************************************************************
*                                                                      *
*     This is ZFITTER `interface'                                      *
*                                                                      *
*     SUBR. ZUXSEC(INDF,SQRS,ZMASS,GAMZ0,GAMEE,GAMFF,XS*)              *
*                                                                      *
*     Routine calculates cross section.                                *
*                                                                      *
*     INDF  (int/read)   = fermion index                               *
*     SQRS  (real/read)  = sqrt(s)                                     *
*     ZMASS (real/read)  = Z0 mass                                     *
*     GAMZ0 (real/read)  = Z0 width                                    *
*     GAMEE (real/read)  = e+e- partial width                          *
*     GAMFF (real/read)  = fermion-fermion partial width               *
*     XS    (real/write) = theoretical cross section                   *
*                                                                      *
*     Called by USER                                                   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8(A-H,O-Z)
*
* flags
*
      PARAMETER(NFLGMX=46)
      COMMON /ZUFLGS/ IFLAGS(NFLGMX),CFLAGS(NFLGMX)
      CHARACTER CFLAGS*4
      PARAMETER (IFAFBC= 1,IFSCAL= 2,IFSCRE= 3,IFAMT4= 4,IFBORN= 5,
     & IFBOXD= 6,IFCONV= 7,IFFINR= 8,IFFOT2= 9,IFGAMS=10,IFDIAG=11,
     & IFINTF=12,IFBARB=13,IFPART=14,IFPOWR=15,IFPRNT=16,IFALEM=17,
     & IFQCDC=18,IFVPOL=19,IFWEAK=20,IFFTJR=21,IFEXPR=22,IFEXPF=23,
     & IFHIGS=24,IFAFMT=25,IFCZAK=26,IFPREC=27,IFHIG2=28,IFALE2=29,
     & IFGFER=30,IFISPP=31,IFFSRS=32,IFMISC=33,IFMISD=34,IFIPFC=35,
     & IFIPSC=36,IFIPTO=37,IFFBHO=38,IFFSPP=39,IFFUNA=40,IFASCR=41,
     & IFSFSR=42,IFENUE=43,IFTUPV=44,IFDMWW=45,IFDSWW=46)
*
* parameters
*
      COMMON /ZUPARS/QDF,QCDCOR(0:14),ALPHST,SIN2TW,S2TEFF(0:11),
     & WIDTHS(0:11)
*
* cuts
*
      COMMON /ZUDATA/ ACOLIN(0:11),EF_MIN(0:11),SPRIME(0:11),
     + ANGMIN(0:11),ANGMAX(0:11),SPRIPP(0:11),IRCUTS(0:11),IRFAST(0:11)
*
* commons for ZFITTER internal use
*
      COMMON/ZFCHMS/ALLCH(0:11),ALLMS(0:11)
      COMMON/GAMFIT/GAMEEC,GAMFIC
      COMMON/FRINIT/NPAR(30),ZPAR(31)
*
*-----------------------------------------------------------------------
*
      IF(INDF.LT.-1 .OR. INDF.GT.11) THEN
        PRINT *,'ZUXSEC> fermion index out of range:  INDF = ',INDF
        RETURN
      ENDIF
      XS = 0.D0
      IF(INDF.EQ.8) RETURN
*
*-----------------------------------------------------------------------
*
      IF(SQRS.LT.9.5D0.OR.SQRS.GT.350D0) THEN
        PRINT*,'REQUIRED S IS OUT OF ALLOWED RANGE'
        STOP
      ENDIF
*
* INTERF. 2
*
      INTRF = 2
*
* This forces total hadronic cross section to be calculated with s'-cut
*cbardin
cb      NPAR(9)=0
*
      S        = SQRS*SQRS
      NPAR(11) = IRCUTS(INDF)
      NPAR(13) = 0
      ZPAR(2)  = ALLCH(INDF)
      ZPAR(4)  = ALLMS(INDF)
cbardin

      DO IQCDC=0,14
      ZPAR(7+IQCDC) = QCDCOR(IQCDC)
      ENDDO
      ZPAR(26) = 1.D0-SPRIME(INDF)/S
      ZPAR(27) = ACOLIN(INDF)
      ZPAR(28) = EF_MIN(INDF)
      ZPAR(29) = ANGMAX(INDF)
      ZPAR(30) = ANGMIN(INDF)
      ZPAR(31) = SPRIPP(INDF)
      GAMEEC   = GAMEE
      GAMFIC   = GAMFF
*
      IF(INDF.NE.11) THEN
        CALL ZCUT(INTRF,IRFAST(INDF),INDF,S,ZMASS,GAMZ0,WIDTHS,SIN2TW,
     +   NPAR,ZPAR,SBORNB,STOT,AFBDUM,AFBDUM)
      ELSE
        CALL BHANG(INTRF,                 S,ZMASS,GAMZ0,WIDTHS,SIN2TW,
     +   NPAR,ZPAR,SBORNB,STOT,AFBDUM,AFBDUM)
      ENDIF
*
      IF(IFLAGS(IFBORN).EQ.1) THEN
        XS = SBORNB
      ELSE
        XS = STOT
      ENDIF
      RETURN
*                                                             END ZUXSEC
      END
 
      SUBROUTINE ZUXSA(INDF,SQRS,ZMASS,GAMZ0,MODE,GVE,XE,GVF,XF,XS,AFB)
*     ========== ======================================================
************************************************************************
*                                                                      *
*     This is ZFITTER `interface'                                      *
*                                                                      *
*     SUBR. ZUXSA(INDF,SQRS,ZMASS,GAMZ0,MODE,GVE,XE,GVF,XF,XS*,AFB*)   *
*                                                                      *
*     Routine calculates cross section and asymmetry.                  *
*                                                                      *
*     INDF (int/read)   = fermion index                                *
*     SQRS (real/read)  = sqrt(s)                                      *
*     ZMASS(real/read)  = Z0 mass                                      *
*     GAMZ0(real/read)  = Z0 width                                     *
*     MODE (int/read)   = meaning of XE/XF (electron/fermion)          *
*                       = 0 axial vector couplings                     *
*                       = 1 rho                                        *
*     GVE  (real/read)  = effective vector coupling (electron)         *
*     XE   (real/read)  = axial vector coupling or rho (electron)      *
*     GVF  (real/read)  = effective vector coupling (fermion)          *
*     XF   (real/read)  = axial vector coupling or rho (fermion)       *
*     XS   (real/write) = theoretical cross section                    *
*     AFB  (real/write) = theoretical charge asymmetry                 *
*                                                                      *
*     Called by USER                                                   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8(A-H,O-Z)
*
* flags
*
      PARAMETER(NFLGMX=46)
      COMMON /ZUFLGS/ IFLAGS(NFLGMX),CFLAGS(NFLGMX)
      CHARACTER CFLAGS*4
      PARAMETER (IFAFBC= 1,IFSCAL= 2,IFSCRE= 3,IFAMT4= 4,IFBORN= 5,
     & IFBOXD= 6,IFCONV= 7,IFFINR= 8,IFFOT2= 9,IFGAMS=10,IFDIAG=11,
     & IFINTF=12,IFBARB=13,IFPART=14,IFPOWR=15,IFPRNT=16,IFALEM=17,
     & IFQCDC=18,IFVPOL=19,IFWEAK=20,IFFTJR=21,IFEXPR=22,IFEXPF=23,
     & IFHIGS=24,IFAFMT=25,IFCZAK=26,IFPREC=27,IFHIG2=28,IFALE2=29,
     & IFGFER=30,IFISPP=31,IFFSRS=32,IFMISC=33,IFMISD=34,IFIPFC=35,
     & IFIPSC=36,IFIPTO=37,IFFBHO=38,IFFSPP=39,IFFUNA=40,IFASCR=41,
     & IFSFSR=42,IFENUE=43,IFTUPV=44,IFDMWW=45,IFDSWW=46)
*
* parameters
*
      COMMON /ZUPARS/QDF,QCDCOR(0:14),ALPHST,SIN2TW,S2TEFF(0:11),
     & WIDTHS(0:11)
*
* cuts
*
      COMMON /ZUDATA/ ACOLIN(0:11),EF_MIN(0:11),SPRIME(0:11),
     + ANGMIN(0:11),ANGMAX(0:11),SPRIPP(0:11),IRCUTS(0:11),IRFAST(0:11)
*
* commons for ZFITTER internal use
*
      COMMON/ZFCHMS/ALLCH(0:11),ALLMS(0:11)
      COMMON/ASSFIT/ROEC,GVEC,GAEC,ROFC,GVFC,GAFC
      COMMON/FRINIT/NPAR(30),ZPAR(31)
      COMMON/CDZAUX/PARTZA(0:10),PARTZI(0:10),RENFAC(0:10),SRENFC(0:10)
*
      DIMENSION INDQ(10)
      DATA INDQ /0,0,0,0,1,2,3,4,5,6/
*
*-----------------------------------------------------------------------
*
      IF(INDF.EQ.-1 .OR. INDF.EQ.10) THEN
        PRINT *,'ZUXSA>  fermion index out of range:  INDF = ',INDF
        RETURN
      ENDIF
      XS  = 0.D0
      AFB = 0.D0
      IF(INDF.EQ.8) RETURN
*
*-----------------------------------------------------------------------
*
      IF(SQRS.LT.9.5D0.OR.SQRS.GT.350D0) THEN
        PRINT*,'REQUIRED S IS OUT OF ALLOWED RANGE'
        STOP
      ENDIF
*
* INTERF. 3
*
      INTRF = 3
*
      IF(MODE.EQ.1) THEN
        DRE= XE
        VE = GVE
        IF(IFLAGS(IFMISC).EQ.0) THEN
          AE = 1D0
        ELSE
          AE = SRENFC(1)
        ENDIF
        DRF= XF
        VF = GVF
        AF = AE
        IF(INDF.LE.11) THEN 
         IF(IFLAGS(IFMISC).EQ.0) THEN
          AF = 1D0
         ELSE
          AF = SRENFC(INDF)
         ENDIF
        ENDIF
      ELSE
        IF(IFLAGS(IFMISC).EQ.0) THEN
          DRE = 1D0
        ELSE
          DRE = RENFAC(1)
        ENDIF
        VE = 2D0*GVE
        AE = 2D0*XE
        DRF= DRE
        IF(INDF.LT.11) THEN
         IF(IFLAGS(IFMISC).EQ.0) THEN
          DRF = 1D0
         ELSE
          DRF = RENFAC(INDF)
         ENDIF
        ENDIF
        VF = 2D0*GVF
        AF = 2D0*XF
      ENDIF
*
      S        = SQRS*SQRS
      NPAR(11) = IRCUTS(INDF)
      NPAR(13) = 1
      ZPAR(2)  = ALLCH(INDF)
      ZPAR(4)  = ALLMS(INDF)
      DO IQCDC=0,14
      ZPAR(7+IQCDC) = QCDCOR(IQCDC)
      ENDDO
      ZPAR(26) = 1.D0-SPRIME(INDF)/S
      ZPAR(27) = ACOLIN(INDF)
      ZPAR(28) = EF_MIN(INDF)
      ZPAR(29) = ANGMAX(INDF)
      ZPAR(30) = ANGMIN(INDF)
      ZPAR(31) = SPRIPP(INDF)
      ROEC     = DRE
      GVEC     = VE
      GAEC     = AE
      ROFC     = DRF
      GVFC     = VF
      GAFC     = AF
*
      IF(INDF.NE.11) THEN
        CALL ZCUT(INTRF,IRFAST(INDF),INDF,S,ZMASS,GAMZ0,WIDTHS,SIN2TW,
     +   NPAR,ZPAR,SBORNB,STOT,ABORN,ATOT)
      ELSE
        CALL BHANG(INTRF,                 S,ZMASS,GAMZ0,WIDTHS,SIN2TW,
     +   NPAR,ZPAR,SBORNB,STOT,ABORN,ATOT)
      ENDIF
*
      IF(IFLAGS(IFBORN).EQ.1) THEN
        XS  = SBORNB
        AFB = ABORN
      ELSE
        XS  = STOT
        AFB = ATOT
      ENDIF
      RETURN
*                                                             END ZUXSA
      END
 
      SUBROUTINE ZUXSA2(INDF,SQRS,ZMASS,GAMZ0,MODE,GV2,X2,XS,AFB)
*     ========== ================================================
************************************************************************
*                                                                      *
*     This is ZFITTER `interface'                                      *
*                                                                      *
*     SUBR. ZUXSA2(INDF,SQRS,ZMASS,GAMZ0,MODE,GV2,X2,XS*,AFB*)         *
*                                                                      *
*     Routine calculates cross section and asymmetry.                  *
*                                                                      *
*     INDF (int/read)   = fermion index                                *
*     SQRS (real/read)  = sqrt(s)                                      *
*     ZMASS(real/read)  = Z0 mass                                      *
*     GAMZ0(real/read)  = Z0 width                                     *
*     MODE (int/read)   = meaning of X2                                *
*                       = 0 axial vector coupling squared              *
*                       = 1 rho squared                                *
*     GV2  (real/read)  = effective vector coupling squared            *
*     X2   (real/read)  = eff. axial vector coupling or rho, squared   *
*     XS   (real/write) = theoretical cross section                    *
*     AFB  (real/write) = theoretical charge asymmetry                 *
*                                                                      *
*     Called by USER                                                   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8(A-H,O-Z)
*
* flags
*
      PARAMETER(NFLGMX=46)
      COMMON /ZUFLGS/ IFLAGS(NFLGMX),CFLAGS(NFLGMX)
      CHARACTER CFLAGS*4
      PARAMETER (IFAFBC= 1,IFSCAL= 2,IFSCRE= 3,IFAMT4= 4,IFBORN= 5,
     & IFBOXD= 6,IFCONV= 7,IFFINR= 8,IFFOT2= 9,IFGAMS=10,IFDIAG=11,
     & IFINTF=12,IFBARB=13,IFPART=14,IFPOWR=15,IFPRNT=16,IFALEM=17,
     & IFQCDC=18,IFVPOL=19,IFWEAK=20,IFFTJR=21,IFEXPR=22,IFEXPF=23,
     & IFHIGS=24,IFAFMT=25,IFCZAK=26,IFPREC=27,IFHIG2=28,IFALE2=29,
     & IFGFER=30,IFISPP=31,IFFSRS=32,IFMISC=33,IFMISD=34,IFIPFC=35,
     & IFIPSC=36,IFIPTO=37,IFFBHO=38,IFFSPP=39,IFFUNA=40,IFASCR=41,
     & IFSFSR=42,IFENUE=43,IFTUPV=44,IFDMWW=45,IFDSWW=46)
*
* parameters
*
      COMMON /ZUPARS/QDF,QCDCOR(0:14),ALPHST,SIN2TW,S2TEFF(0:11),
     & WIDTHS(0:11)
*
* cuts
*
      COMMON /ZUDATA/ ACOLIN(0:11),EF_MIN(0:11),SPRIME(0:11),
     + ANGMIN(0:11),ANGMAX(0:11),SPRIPP(0:11),IRCUTS(0:11),IRFAST(0:11)
*
* commons for ZFITTER internal use
*
      COMMON/ZFCHMS/ALLCH(0:11),ALLMS(0:11)
      COMMON/AS2FIT/ROL2C,GVL2C,GAL2C
      COMMON/FRINIT/NPAR(30),ZPAR(31)
      COMMON/CDZAUX/PARTZA(0:10),PARTZI(0:10),RENFAC(0:10),SRENFC(0:10)
*
      DIMENSION INDQ(10)
      DATA INDQ /0,0,0,0,1,2,3,4,5,6/
*
*-----------------------------------------------------------------------
*
      IF(INDF.EQ.-1 .OR. INDF.EQ.0 .OR. INDF.GE.4.AND.INDF.LE.10) THEN
        PRINT *,'ZUXSA2> fermion index out of range:  INDF = ',INDF
        RETURN
      ENDIF
*
*-----------------------------------------------------------------------
*
      IF(SQRS.LT.9.5D0.OR.SQRS.GT.350D0) THEN
        PRINT*,'REQUIRED S IS OUT OF ALLOWED RANGE'
        STOP
      ENDIF
*
* INTERF. 4
*
      INTRF = 4
*
      IF(MODE.EQ.1) THEN
        DR2 = X2
        DV2 = GV2
        IF(IFLAGS(IFMISC).EQ.0) THEN
          DA2 = 1D0
        ELSE        
          DA2 = RENFAC(1)
        ENDIF
      ELSE
        IF(IFLAGS(IFMISC).EQ.0) THEN
          DR2 = 1D0
        ELSE        
          DR2 = RENFAC(1)**2
        ENDIF
        DV2 = 4D0*GV2
        DA2 = 4D0*X2
      ENDIF
*
      S        = SQRS*SQRS
      NPAR(11) = IRCUTS(INDF)
      NPAR(13) = 1
      ZPAR(2)  = ALLCH(INDF)
      ZPAR(4)  = ALLMS(INDF)
      DO IQCDC=0,14
      ZPAR(7+IQCDC) = QCDCOR(IQCDC)
      ENDDO
      ZPAR(26) = 1.D0-SPRIME(INDF)/S
      ZPAR(27) = ACOLIN(INDF)
      ZPAR(28) = EF_MIN(INDF)
      ZPAR(29) = ANGMAX(INDF)
      ZPAR(30) = ANGMIN(INDF)
      ZPAR(31) = SPRIPP(INDF)
      ROL2C    = DR2
      GVL2C    = DV2
      GAL2C    = DA2
*
      IF(INDF.NE.11) THEN
        CALL ZCUT(INTRF,IRFAST(INDF),INDF,S,ZMASS,GAMZ0,WIDTHS,SIN2TW,
     +   NPAR,ZPAR,SBORNB,STOT,ABORN,ATOT)
      ELSE
        CALL BHANG(INTRF,                 S,ZMASS,GAMZ0,WIDTHS,SIN2TW,
     +   NPAR,ZPAR,SBORNB,STOT,ABORN,ATOT)
      ENDIF
*
      IF(IFLAGS(IFBORN).EQ.1) THEN
        XS  = SBORNB
        AFB = ABORN
      ELSE
        XS  = STOT
        AFB = ATOT
      ENDIF
      RETURN
*                                                            END ZUXSA2
      END
 
      SUBROUTINE ZUXAFB(INDF,SQRS,ZMASS,GAMZ0,PFOUR,PVAE2,PVAF2,XS,AFB)
*     ========== ======================================================
************************************************************************
*                                                                      *
*     This is ZFITTER `interface'                                      *
*                                                                      *
*     SUBR. ZUXAFB(INDF,SQRS,ZMASS,GAMZ0,PFOUR,PVAE2,PVAF2,XS*,AFB*)   *
*                                                                      *
*     Routine calculates cross section and asymmetry.                  *
*                                                                      *
*     INDF (int/read)   = fermion index                                *
*     SQRS (real/read)  = sqrt(s)                                      *
*     ZMASS(real/read)  = Z0 mass                                      *
*     GAMZ0(real/read)  = Z0 width                                     *
*     MODE (int/read)   = 0 axial vector coupling squared              *
*                       = 1 disabled here                              *
*     PFOUR(real/read)  = VE*AE*VF*AF                                  *
*     PVAE2(real/read)  = VE**2+AE**2                                  *
*     PVAF2(real/read)  = VF**2+AF**2                                  *
*     XS   (real/write) = theoretical cross section                    *
*     AFB  (real/write) = theoretical charge asymmetry                 *
*                                                                      *
*     Called by USER                                                   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8(A-H,O-Z)
*
* flags
*
      PARAMETER(NFLGMX=46)
      COMMON /ZUFLGS/ IFLAGS(NFLGMX),CFLAGS(NFLGMX)
      CHARACTER CFLAGS*4
      PARAMETER (IFAFBC= 1,IFSCAL= 2,IFSCRE= 3,IFAMT4= 4,IFBORN= 5,
     & IFBOXD= 6,IFCONV= 7,IFFINR= 8,IFFOT2= 9,IFGAMS=10,IFDIAG=11,
     & IFINTF=12,IFBARB=13,IFPART=14,IFPOWR=15,IFPRNT=16,IFALEM=17,
     & IFQCDC=18,IFVPOL=19,IFWEAK=20,IFFTJR=21,IFEXPR=22,IFEXPF=23,
     & IFHIGS=24,IFAFMT=25,IFCZAK=26,IFPREC=27,IFHIG2=28,IFALE2=29,
     & IFGFER=30,IFISPP=31,IFFSRS=32,IFMISC=33,IFMISD=34,IFIPFC=35,
     & IFIPSC=36,IFIPTO=37,IFFBHO=38,IFFSPP=39,IFFUNA=40,IFASCR=41,
     & IFSFSR=42,IFENUE=43,IFTUPV=44,IFDMWW=45,IFDSWW=46)
*
* parameters
*
      COMMON /ZUPARS/QDF,QCDCOR(0:14),ALPHST,SIN2TW,S2TEFF(0:11),
     & WIDTHS(0:11)
*
* cuts
*
      COMMON /ZUDATA/ ACOLIN(0:11),EF_MIN(0:11),SPRIME(0:11),
     + ANGMIN(0:11),ANGMAX(0:11),SPRIPP(0:11),IRCUTS(0:11),IRFAST(0:11)
*
* commons for ZFITTER internal use
*
       COMMON/ZFCHMS/ALLCH(0:11),ALLMS(0:11)
       COMMON/AFBFIT/FOUR,VAE2,VAF2
       COMMON/FRINIT/NPAR(30),ZPAR(31)
*
      DIMENSION INDQ(10)
      DATA INDQ /0,0,0,0,1,2,3,4,5,6/
*
*----------------------------------------------------------------------
*
      IF(INDF.EQ.-1 .OR. INDF.EQ.0 .OR. INDF.GE.4.AND.INDF.LE.11) THEN
        PRINT *,'ZUXAFB> fermion index out of range:  INDF = ',INDF
        RETURN
      ENDIF
*
*-----------------------------------------------------------------------
*
      IF(SQRS.LT.9.5D0.OR.SQRS.GT.350D0) THEN
        PRINT*,'REQUIRED S IS OUT OF ALLOWED RANGE'
        STOP
      ENDIF
*
      FOUR=16*PFOUR
      VAE2=4*PVAE2
      VAF2=4*PVAF2
*
* INTERF. 5
*
      INTRF = 5
*
      S        = SQRS*SQRS
      NPAR(11) = IRCUTS(INDF)
      NPAR(13) = 1
      ZPAR(2)  = ALLCH(INDF)
      ZPAR(4)  = ALLMS(INDF)
      DO IQCDC=0,14
      ZPAR(7+IQCDC) = QCDCOR(IQCDC)
      ENDDO
      ZPAR(26) = 1.D0-SPRIME(INDF)/S
      ZPAR(27) = ACOLIN(INDF)
      ZPAR(28) = EF_MIN(INDF)
      ZPAR(29) = ANGMAX(INDF)
      ZPAR(30) = ANGMIN(INDF)
      ZPAR(31) = SPRIPP(INDF)
*
      IF(INDF.NE.11) THEN
        CALL ZCUT(INTRF,IRFAST(INDF),INDF,S,ZMASS,GAMZ0,WIDTHS,SIN2TW,
     +   NPAR,ZPAR,SBORNB,STOT,ABORN,ATOT)
      ELSE
        CALL BHANG(INTRF,                 S,ZMASS,GAMZ0,WIDTHS,SIN2TW,
     +   NPAR,ZPAR,SBORNB,STOT,ABORN,ATOT)
      ENDIF
*
      IF(IFLAGS(IFBORN).EQ.1) THEN
        XS  = SBORNB
        AFB = ABORN
      ELSE
        XS  = STOT
        AFB = ATOT
      ENDIF
      RETURN
*                                                            END ZUXAFB
      END
 
      SUBROUTINE 
     +          ZUTAU(SQRS,ZMASS,GAMZ0,MODE,GVE,XE,GVF,XF,TAUPOL,TAUAFB)
*     ==================================================================
************************************************************************
*                                                                      *
*     This is ZFITTER `interface'                                      *
*                                                                      *
*     SUBR. ZUTAU(SQRS,ZMA,GAMZ,MODE,GVE,XE,GVF,XF,TAUPOL*,TAUAFB*)    *
*                                                                      *
*     Routine calculates tau polarisation and polarization asymmetry.  *
*                                                                      *
*     SQRS   (real/read)  = sqrt(s)                                    *
*     ZMASS  (real/read)  = Z0 mass                                    *
*     GAMZ0  (real/read)  = Z0 width                                   *
*     MODE   (int/read)   = meaning of XE/XF (electron/fermion)        *
*                         = 0 axial vector couplings                   *
*                         = 1 rho                                      *
*     GVE    (real/read)  = effective vector coupling (electron)       *
*     XE     (real/read)  = axial vector coupling or rho (electron)    *
*     GVF    (real/read)  = effective vector coupling (fermion)        *
*     XF     (real/read)  = axial vector coupling or rho (fermion)     *
*     TAUPOL (real/write) = tau polarization                           *
*     TAUAFB (real/write) = tau polarization asymmetry                 *
*                                                                      *
*     Called by USER                                                   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8(A-H,O-Z)
*
* flags
*
      PARAMETER(NFLGMX=46)
      COMMON /ZUFLGS/ IFLAGS(NFLGMX),CFLAGS(NFLGMX)
      CHARACTER CFLAGS*4
      PARAMETER (IFAFBC= 1,IFSCAL= 2,IFSCRE= 3,IFAMT4= 4,IFBORN= 5,
     & IFBOXD= 6,IFCONV= 7,IFFINR= 8,IFFOT2= 9,IFGAMS=10,IFDIAG=11,
     & IFINTF=12,IFBARB=13,IFPART=14,IFPOWR=15,IFPRNT=16,IFALEM=17,
     & IFQCDC=18,IFVPOL=19,IFWEAK=20,IFFTJR=21,IFEXPR=22,IFEXPF=23,
     & IFHIGS=24,IFAFMT=25,IFCZAK=26,IFPREC=27,IFHIG2=28,IFALE2=29,
     & IFGFER=30,IFISPP=31,IFFSRS=32,IFMISC=33,IFMISD=34,IFIPFC=35,
     & IFIPSC=36,IFIPTO=37,IFFBHO=38,IFFSPP=39,IFFUNA=40,IFASCR=41,
     & IFSFSR=42,IFENUE=43,IFTUPV=44,IFDMWW=45,IFDSWW=46)
*
* parameters
*
      COMMON /ZUPARS/QDF,QCDCOR(0:14),ALPHST,SIN2TW,S2TEFF(0:11),
     & WIDTHS(0:11)
*
* cuts
*
      COMMON /ZUDATA/ ACOLIN(0:11),EF_MIN(0:11),SPRIME(0:11),
     + ANGMIN(0:11),ANGMAX(0:11),SPRIPP(0:11),IRCUTS(0:11),IRFAST(0:11)
*
* commons for ZFITTER internal use
*
      COMMON/ZFCHMS/ALLCH(0:11),ALLMS(0:11)
      COMMON/ASSFIT/ROEC,GVEC,GAEC,ROFC,GVFC,GAFC
      COMMON/FRINIT/ NPAR(30),ZPAR(31)
*
*-----------------------------------------------------------------------
*
      IF(SQRS.LT.9.5D0.OR.SQRS.GT.350D0) THEN
        PRINT*,'REQUIRED S IS OUT OF ALLOWED RANGE'
        STOP
      ENDIF
*
* INTERF. 3
*
      INTRF = 3
*
      VE   = 2.D0*GVE
      VF   = 2.D0*GVF
      IF(MODE.EQ.1) THEN
        DRE = XE
        AE  = 1.D0
        DRF = XF
        AF  = 1.D0
      ELSE
        DRE = 1.D0
        AE  = 2.D0*XE
        DRF = 1.D0
        AF  = 2.D0*XF
      ENDIF
*
      S        = SQRS*SQRS
      NPAR(11) = IRCUTS(3)
      ICURR    = NPAR(13)
      NPAR(13) = 1
      ZPAR(2)  = ALLCH(3)
      ZPAR(4)  = ALLMS(3)
      DO IQCDC=0,14
      ZPAR(7+IQCDC) = QCDCOR(IQCDC)
      ENDDO
      ZPAR(26) = 1.D0-SPRIME(3)/S
      ZPAR(27) = ACOLIN(3)
      ZPAR(28) = EF_MIN(3)
      ZPAR(29) = ANGMAX(3)
      ZPAR(30) = ANGMIN(3)
      ZPAR(31) = SPRIPP(3)
      ROEC     = DRE
      GVEC     = VE
      GAEC     = AE
      ROFC     = DRF
      GVFC     = VF
      GAFC     = AF
*
      ZPAR(24) = 0.D0
      ZPAR(25) = 1.D0
      CALL ZCUT(INTRF,IRFAST(3),3,S,ZMASS,GAMZ0,WIDTHS,SIN2TW,
     + NPAR,ZPAR,SBORPL,SQEDPL,ABPL,AQPL)
      ZPAR(25) =-1.D0
      CALL ZCUT(INTRF,IRFAST(3),3,S,ZMASS,GAMZ0,WIDTHS,SIN2TW,
     + NPAR,ZPAR,SBORMI,SQEDMI,ABMI,AQMI)
*
      IF(IFLAGS(IFBORN).EQ.1) THEN
        TAUPOL = (SBORPL-SBORMI)/(SBORPL+SBORMI)
        TAUAFB = (ABPL*SBORPL - ABMI*SBORMI)/(SBORPL+SBORMI)
      ELSE
        TAUPOL = (SQEDPL-SQEDMI)/(SQEDPL+SQEDMI)
        TAUAFB = (AQPL*SQEDPL - AQMI*SQEDMI)/(SQEDPL+SQEDMI)
      ENDIF
*
      ZPAR(25) = 0.D0
      NPAR(13) = ICURR
      RETURN
*                                                              END ZUTAU
      END

      SUBROUTINE 
     +          ZUALR(SQRS,ZMASS,GAMZ0,MODE,GVE,XE,GVF,XF,TAUPOL,TAUAFB)
*     ==================================================================
************************************************************************
*                                                                      *
*     This is ZFITTER `interface'                                      *
*                                                                      *
*     SUBR. ZUALR(SQRS,ZMA,GAMZ,MODE,GVE,XE,GVF,XF,???????????????)    *
*                                                                      *
*     Routine calculates tau polarisation and polarization asymmetry.  *
*                                                                      *
*     SQRS   (real/read)  = sqrt(s)                                    *
*     ZMASS  (real/read)  = Z0 mass                                    *
*     GAMZ0  (real/read)  = Z0 width                                   *
*     MODE   (int/read)   = meaning of XE/XF (electron/fermion)        *
*                         = 0 axial vector couplings                   *
*                         = 1 rho                                      *
*     GVE    (real/read)  = effective vector coupling (electron)       *
*     XE     (real/read)  = axial vector coupling or rho (electron)    *
*     GVF    (real/read)  = effective vector coupling (fermion)        *
*     XF     (real/read)  = axial vector coupling or rho (fermion)     *
*     TAUPOL (real/write) = tau polarization                           *
*     TAUAFB (real/write) = tau polarization asymmetry                 *
*                                                                      *
*     Called by USER                                                   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8(A-H,O-Z)
*
* flags
*
      PARAMETER(NFLGMX=46)
      COMMON /ZUFLGS/ IFLAGS(NFLGMX),CFLAGS(NFLGMX)
      CHARACTER CFLAGS*4
      PARAMETER (IFAFBC= 1,IFSCAL= 2,IFSCRE= 3,IFAMT4= 4,IFBORN= 5,
     & IFBOXD= 6,IFCONV= 7,IFFINR= 8,IFFOT2= 9,IFGAMS=10,IFDIAG=11,
     & IFINTF=12,IFBARB=13,IFPART=14,IFPOWR=15,IFPRNT=16,IFALEM=17,
     & IFQCDC=18,IFVPOL=19,IFWEAK=20,IFFTJR=21,IFEXPR=22,IFEXPF=23,
     & IFHIGS=24,IFAFMT=25,IFCZAK=26,IFPREC=27,IFHIG2=28,IFALE2=29,
     & IFGFER=30,IFISPP=31,IFFSRS=32,IFMISC=33,IFMISD=34,IFIPFC=35,
     & IFIPSC=36,IFIPTO=37,IFFBHO=38,IFFSPP=39,IFFUNA=40,IFASCR=41,
     & IFSFSR=42,IFENUE=43,IFTUPV=44,IFDMWW=45,IFDSWW=46)
*
* parameters
*
      COMMON /ZUPARS/QDF,QCDCOR(0:14),ALPHST,SIN2TW,S2TEFF(0:11),
     & WIDTHS(0:11)
*
* cuts
*
      COMMON /ZUDATA/ ACOLIN(0:11),EF_MIN(0:11),SPRIME(0:11),
     + ANGMIN(0:11),ANGMAX(0:11),SPRIPP(0:11),IRCUTS(0:11),IRFAST(0:11)
*
* commons for ZFITTER internal use
*
      COMMON/ZFCHMS/ALLCH(0:11),ALLMS(0:11)
      COMMON/ASSFIT/ROEC,GVEC,GAEC,ROFC,GVFC,GAFC
      COMMON/FRINIT/ NPAR(30),ZPAR(31)
*
*-----------------------------------------------------------------------
*
      IF(SQRS.LT.9.5D0.OR.SQRS.GT.350D0) THEN
        PRINT*,'REQUIRED S IS OUT OF ALLOWED RANGE'
        STOP
      ENDIF
*
* INTERF. 6
*
*     INTRF = 6
*
* EMPTY!!!!!!!!!!!
*                                                              END ZUALR
      END
