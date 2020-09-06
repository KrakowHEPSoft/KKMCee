INCLUDES = -I$(top_srcdir)/Rhadr $(all_includes) -I$(ROOTINCDIR)
METASOURCES = AUTO
bin_PROGRAMS = MainP
MainP_SOURCES = MainP.cxx Semaph.cxx ROBOL.cxx ROBOL2.cxx KKMC.cxx PartLund.cxx
CXXFLAGS = -g -O2 -std=c++11 -lstdc++

############################################
#  Tables (re)generated in flight.
#  To be xchecked if it works correctly!!!
############################################
KK2f_Libs2= \
        $(top_builddir)/KK2f/libKK2f.a     $(top_builddir)/glibk/libGlibk.a \
        $(top_builddir)/jetset/libJetset.a $(top_builddir)/bornv/libBornvDizet.a \
        $(top_builddir)/tauola/libTauola.a $(top_builddir)/photos/libPhotos.a \
        $(top_builddir)/dizet/libDizet.a \
        $(top_builddir)/KKsem/libKKsem.a  -lstdc++



MainP_LDADD = $(KK2f_Libs2) $(ROOTAUXLIBS) $(ROOTGLIBS) -lgfortran
MainP_LDFLAGS = -L$(ROOTLIBDIR)

SUBDIRS = Plots

########################################################################
###                   USER PART
###  make makprod       <=== mandatory
###------------------------------------
###  make demo-start
###  make Inclusive-start
###  make figinc-pdf
###  make Mu-start
# Note on linking EW library:
#  (a) EW formfactors read from tables 
#      link /bornv/libBornvTabs.a
#  (b) EW formfactors created in flight, now it is DEFAULT!!!
#      link /bornv/libBornvDizet.a   and /dizet/libDizet.a
########################################################################
check_all:
	(cd ../KK2f;   make install)
	(cd ../bornv;  make install)
	(cd ../glibk;  make install)
	(cd ../dizet;  make install)
	(cd ../tauola; make install)
	(cd ../photos; make install)
	(cd ../jetset; make install)
	(make install)
#-----------------------------------------------------------------------
demoC-start: check_all
	(cd ./demoC;  make start )
demoC-stop:
	(cd ./demoC;  make stop )
demoC-cont:
	(cd ./demoC;  make cont )
#-----------------------------------------------------------------------
#                  update production makefiles
makprod:
	(sed -e 's#<main>#MainP#'  -e 's#<dset>#demoC#' -e 's#<class>#8nh#' ./makefile.production > ./demoC/makefile)
#######################################################################
AM_FFLAGS = @MY_F77_FLAGS@
