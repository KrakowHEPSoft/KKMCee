INCLUDES = -I$(top_srcdir)/KKsem -I$(top_srcdir)/MaMar $(all_includes) \
	-I$(ROOTINCDIR)
METASOURCES = AUTO

CXXFLAGS = -g -O2 -std=c++11 -lstdc++  

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
	$(top_builddir)/KKsem/libKKsem.a -lstdc++

Plot1_SOURCES = KKsem.cxx Plot1.cxx
Plot1_LDADD   = $(KK2f_Libs2) $(ROOTAUXLIBS) $(ROOTGLIBS) -lgfortran 
Plot1_LDFLAGS = -L$(ROOTLIBDIR) 

Plot2_SOURCES = Plot2.cxx
Plot2_LDADD   = $(KK2f_Libs2) $(ROOTAUXLIBS) $(ROOTGLIBS) -lgfortran 
Plot2_LDFLAGS = -L$(ROOTLIBDIR) 

PlotFSR_SOURCES = KKsem.cxx PlotFSR.cxx
PlotFSR_LDADD   = $(KK2f_Libs2) $(ROOTAUXLIBS) $(ROOTGLIBS) -lFoam -lgfortran 
PlotFSR_LDFLAGS = -L$(ROOTLIBDIR) 

######################### private #############################
plot1: Plot1
	./Plot1
###############################################################
