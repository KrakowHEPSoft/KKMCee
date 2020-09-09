AM_CPPFLAGS = -I$(top_srcdir)/bornv
METASOURCES = AUTO
lib_LIBRARIES = libBornvTabs.a libBornvDizet.a
libBornvDizet_a_SOURCES = BornV.f BStra.f FoamC.f IRC.f MBrA.f\
	MBrB.f ranlux.f Vesk1.f
libBornvTabs_a_SOURCES = BornV.f BornV_StartEW.f BStra.f FoamC.f\
	IRC.f MBrA.f MBrB.f ranlux.f Vesk1.f

##########################################################################
mbr:
	(sed -e 's#MBrA#MBrB#'    MBrA.h > MBrB.h)
	(sed -e 's#MBrA#MBrB#'    MBrA.f > MBrB.f)
foam:
	(sed -e 's#FoamA#FoamB#'g    FoamA.h > FoamB.h)
	(sed -e 's#FoamA#FoamB#'g    FoamA.f > FoamB.f)
	(sed -e 's#FoamA#FoamC#'g    FoamA.h > FoamC.h)
	(sed -e 's#FoamA#FoamC#'g    FoamA.f > FoamC.f)
##############################################################################

AM_FFLAGS = @MY_F77_FLAGS@
