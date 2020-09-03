
METASOURCES = AUTO
bin_PROGRAMS = DEMO ProdMC

DEMO_SOURCES = demo.f
ProdMC_SOURCES = ProdMC.f
DEMO_LDADD = $(top_builddir)/KK2f/libKK2f.a \
	$(top_builddir)/bornv/libBornvDizet.a $(top_builddir)/dizet/libDizet.a $(top_builddir)/photos/libPhotos.a \
	$(top_builddir)/tauola/libTauola.a $(top_builddir)/jetset/libJetset.a $(top_builddir)/glibk/libGlibk.a \
        $(top_builddir)/tauola/glib.a -lstdc++

ProdMC_LDADD = $(top_builddir)/KK2f/libKK2f.a \
	$(top_builddir)/bornv/libBornvDizet.a $(top_builddir)/dizet/libDizet.a $(top_builddir)/photos/libPhotos.a \
	$(top_builddir)/tauola/libTauola.a $(top_builddir)/jetset/libJetset.a $(top_builddir)/glibk/libGlibk.a \
	$(top_builddir)/tauola/glib.a -lstdc++

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
demo-start: check_all
	(cd ./demo;  make start )
#-----------------------------------------------------------------------
demo_R-start:  check_all
	(cd ./demo_R;  make start )
#-----------------------------------------------------------------------
Inclusive-start:   check_all
	(cd ./Inclusive;  make start )
Inclusive-stop:
	(cd ./Inclusive;  make stop )
#-----------------------------------------------------------------------
figinc-pdf:
	(cd Plots; make figinc-pdf)
#-----------------------------------------------------------------------
Test-start: check_all
	(cd ./Test;  make start )
Test-stop:
	(cd ./Test;  make stop )
Test-cont:
	(cd ./Test;  make cont )
#-----------------------------------------------------------------------
Mu-start: check_all
	(cd ./Mu;  make start )
Mu-stop:
	(cd ./Mu;  make stop )
Mu-cont:
	(cd ./Mu;  make cont )
Mu-plot:
	(cd ./Plots; make figmix)
#-----------------------------------------------------------------------
Tau-start: check_all
	(cd ./Tau;  make start )
Tau-stop:
	(cd ./Tau;  make stop )
#-----------------------------------------------------------------------
Down-start: check_all
	(cd ./Down;  make start )
#-----------------------------------------------------------------------
Up-start: check_all
	(cd ./Up;  make start )
#-----------------------------------------------------------------------
Botom-start: check_all
	(cd ./Botom;  make start )
#-----------------------------------------------------------------------
Beast-start:
	(cd ./Beast;  make start )
#-----------------------------------------------------------------------
#                  update production makefiles
#-----------------------------------------------------------------------
makprod:
	(sed -e 's#<main>#DEMO#'   -e 's#<dset>#demo#'  -e 's#<class>#long#' ./makefile.production > ./demo/makefile)
	(sed -e 's#<main>#ProdMC#' -e 's#<dset>#Test#'  -e 's#<class>#long#' ./makefile.production > ./Test/makefile)
	(sed -e 's#<main>#ProdMC#' -e 's#<dset>#Mu#'    -e 's#<class>#long#' ./makefile.production > ./Mu/makefile)
	(sed -e 's#<main>#ProdMC#' -e 's#<dset>#Tau#'   -e 's#<class>#long#' ./makefile.production > ./Tau/makefile)
	(sed -e 's#<main>#ProdMC#' -e 's#<dset>#Down#'  -e 's#<class>#long#' ./makefile.production > ./Down/makefile)
	(sed -e 's#<main>#ProdMC#' -e 's#<dset>#Up#'    -e 's#<class>#long#' ./makefile.production > ./Up/makefile)
	(sed -e 's#<main>#ProdMC#' -e 's#<dset>#Botom#' -e 's#<class>#long#' ./makefile.production > ./Botom/makefile)
	(sed -e 's#<main>#ProdMC#' -e 's#<dset>#Inclusive#' -e 's#<class>#long#' ./makefile.production > ./Inclusive/makefile)
	(sed -e 's#<main>#ProdMC#' -e 's#<dset>#Beast#' -e 's#<class>#long#' ./makefile.production > ./Beast/makefile)
	(sed -e 's#<main>#ProdMC#' -e 's#<dset>#demo_R#' -e 's#<class>#long#' ./makefile.production > ./demo_R/makefile)
#######################################################################
AM_FFLAGS = @MY_F77_FLAGS@
