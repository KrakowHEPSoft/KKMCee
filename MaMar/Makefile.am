INCLUDES = -I$(top_srcdir)/MaMar $(all_includes) -I$(ROOTINCDIR)
METASOURCES = AUTO
bin_PROGRAMS = Test
Test_SOURCES = Test.cxx KKMC.cxx PartLund.cxx ROBOL.cxx Semaph.cxx

Test_LDADD = $(top_builddir)/KK2f/libKK2f.a $(top_builddir)/glibk/libGlibk.a \
	$(top_builddir)/jetset/libJetset.a $(top_builddir)/bornv/libBornvDizet.a \
	$(top_builddir)/tauola/libTauola.a $(top_builddir)/photos/libPhotos.a $(top_builddir)/dizet/libDizet.a \
	$(ROOTAUXLIBS) $(ROOTGLIBS)

Test_LDFLAGS = -L$(ROOTLIBDIR) -lgfortran
SUBDIRS = test0 Plots

#-----------------------------------------------------------------------
test0-start:	Test
	(cd ./test0;  make start )
test0-stop:
	(cd ./test0;  make stop )
#-----------------------------------------------------------------------
