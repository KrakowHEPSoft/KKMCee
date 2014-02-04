INCLUDES = -I$(top_srcdir)/MaMar $(all_includes) -I$(ROOTINCDIR)
METASOURCES = AUTO

SUBDIRS = test0 test2 Plots

bin_PROGRAMS = Test Main2

Test_SOURCES = Test.cxx KKMC.cxx PartLund.cxx ROBOL.cxx Semaph.cxx
Test_LDADD = $(top_builddir)/KK2f/libKK2f.a $(top_builddir)/glibk/libGlibk.a \
	$(top_builddir)/jetset/libJetset.a $(top_builddir)/bornv/libBornvDizet.a \
	$(top_builddir)/tauola/libTauola.a $(top_builddir)/photos/libPhotos.a $(top_builddir)/dizet/libDizet.a \
	$(ROOTAUXLIBS) $(ROOTGLIBS)
Test_LDFLAGS = -L$(ROOTLIBDIR) -lgfortran
SUBDIRS = test0 test2 Plots

Main2_SOURCES = Main2.cxx KKMC.cxx PartLund.cxx ROBOL2.cxx Semaph.cxx
Main2_LDADD = $(top_builddir)/KK2f/libKK2f.a $(top_builddir)/glibk/libGlibk.a \
	$(top_builddir)/jetset/libJetset.a $(top_builddir)/bornv/libBornvDizet.a \
	$(top_builddir)/tauola/libTauola.a $(top_builddir)/photos/libPhotos.a $(top_builddir)/dizet/libDizet.a \
	$(ROOTAUXLIBS) $(ROOTGLIBS)
Main2_LDFLAGS = -L$(ROOTLIBDIR) -lgfortran

AM_FFLAGS = @MY_F77_FLAGS@

########################################################################
###                   USER PART
###  make test0-start
###  make Plot1
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
test0-start:	check_all
	(cd ./test0;  make start )
test0-stop:
	(cd ./test0;  make stop )
Plot1:		check_all
	(cd Plots; make Plot1; ./Plot1)
#-----------------------------------------------------------------------
