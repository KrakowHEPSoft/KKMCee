INCLUDES = -I$(top_srcdir)/KKsem -I$(top_srcdir)/MaMar $(all_includes) \
	-I$(ROOTINCDIR)
METASOURCES = AUTO

bin_PROGRAMS = Plot1 Plot2 PlotFSR

########### KK2f MC libraries ###########
$(top_builddir)/KK2f/libKK2f.a:
	(cd $(top_builddir)/KK2f/; make)
$(top_builddir)/glibk/libGlibk.a:
	(cd $(top_builddir)/glibk/; make)
$(top_builddir)/jetset/libJetset.a:
	(cd $(top_builddir)/jetset/; make)
$(top_builddir)/bornv/libBornvDizet.a:
	(cd $(top_builddir)/bornv/; make)
$(top_builddir)/tauola/libTauola.a:
	(cd $(top_builddir)/tauola/; make)
$(top_builddir)/photos/libPhotos.a:
	(cd $(top_builddir)/photos/; make)
$(top_builddir)/dizet/libDizet.a:
	(cd $(top_builddir)/dizet/; make)
############### extras ##################
$(top_builddir)/KKsem/libKKsem.a:
	(cd $(top_builddir)/KKsem/; make)
$(top_builddir)/bornv/libBornvTabs.a:
	(cd $(top_builddir)/bornv/; make)
#########################################
	
KK2f_Libs2= \
	$(top_builddir)/KK2f/libKK2f.a $(top_builddir)/glibk/libGlibk.a \
	$(top_builddir)/jetset/libJetset.a $(top_builddir)/bornv/libBornvTabs.a \
	$(top_builddir)/tauola/libTauola.a $(top_builddir)/photos/libPhotos.a \
	$(top_builddir)/dizet/libDizet.a \
	$(top_builddir)/KKsem/libKKsem.a

Plot1_SOURCES = KKsem.cxx Plot1.cxx
Plot1_LDADD   = $(KK2f_Libs2) $(ROOTAUXLIBS) $(ROOTGLIBS)
Plot1_LDFLAGS = -L$(ROOTLIBDIR) -lgfortran

Plot2_SOURCES = Plot2.cxx
Plot2_LDADD   = $(KK2f_Libs2) $(ROOTAUXLIBS) $(ROOTGLIBS)
Plot2_LDFLAGS = -L$(ROOTLIBDIR) -lgfortran

PlotFSR_SOURCES = KKsem.cxx PlotFSR.cxx
PlotFSR_LDADD   = $(KK2f_Libs2) $(ROOTAUXLIBS) $(ROOTGLIBS)
PlotFSR_LDFLAGS = -L$(ROOTLIBDIR) -lgfortran

######################### private #############################
plot1: Plot1
	./Plot1
###############################################################