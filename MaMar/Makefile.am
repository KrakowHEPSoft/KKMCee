# ****** This is input data file for Automake ****** 
# Automake creates Makefile using all information defined below
INCLUDES = -I$(top_srcdir)/MaMar $(all_includes) -I$(ROOTINCDIR) -I$(top_srcdir)/MCdev -I$(top_srcdir)/Foam2
METASOURCES = AUTO
SUBDIRS = test0 workAFB workZinv workFSR workFOAM workFOKK Plots
#---------------------------------------------------------------
#   ROOT Dictionary file definition
dictC = KKfmDict.cxx
dictH = KKfmDict.h
nodist_libKKfm_la_SOURCES = $(dictC)
noinst_HEADERS = LinkDef.h
CLEANFILES = $(dictC) $(dictH)
#   ROOT Directory generation:
ROOTCINTHDRS= $(include_HEADERS) $(noinst_HEADERS)
$(dictC) : $(ROOTCINTHDRS)
	$(ROOTCINT) -f $(dictC) -c $(INCLUDES) $(ROOTCINTHDRS)
#---------------------------------------------------------------
#  Local library
#Header files:
include_HEADERS = TMCgenFOAM.h TRobolFOAM.h
#Additional flags:
lib_LTLIBRARIES = libKKfm.la
libKKfm_la_SOURCES =  TMCgenFOAM.cxx TRobolFOAM.cxx
libKKfm_la_LIBADD =  $(top_builddir)/MCdev/libMCdev.la \
    $(ROOTAUXLIBS) $(ROOTGLIBS) -lFoam $(FLIBS)
libKKfm_la_LDFLAGS = -L$(ROOTLIBDIR)
# explicit target in ../
$(top_builddir)/MCdev/libMCdev.la: $(top_builddir)/MCdev/*.h $(top_builddir)/MCdev/*.cxx
	(cd $(top_builddir)/MCdev/; make)

#---------------------------------------------------------------

########### KK2f MC libraries #############
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
############################################
KK2f_Libs= \
	$(top_builddir)/KK2f/libKK2f.a     $(top_builddir)/glibk/libGlibk.a \
	$(top_builddir)/jetset/libJetset.a $(top_builddir)/bornv/libBornvDizet.a \
	$(top_builddir)/tauola/libTauola.a $(top_builddir)/photos/libPhotos.a \
	$(top_builddir)/dizet/libDizet.a

#---------------------------------------------------------------
# Local executables
bin_PROGRAMS = Prod0 ProdZinv ProdFSR ProdFoam MainFOAM
AM_FFLAGS = @MY_F77_FLAGS@

MainFOAM_SOURCES = MainFOAM.cxx 
MainFOAM_LDADD = $(top_builddir)/MCdev/.libs/libMCdev.a \
				 $(top_builddir)/Foam2/.libs/libFOAM.a  \
    			 $(ROOTAUXLIBS) $(ROOTGLIBS) -lFoam $(FLIBS) \
    			 $(top_builddir)/MaMar/libKKfm.la \
    			 $(KK2f_Libs)
MainFOAM_LDFLAGS = -L$(ROOTLIBDIR) -rpath $(ROOTRPATH) -g

Prod0_SOURCES = Prod0.cxx KKMC.cxx PartLund.cxx ROBOL.cxx Semaph.cxx
Prod0_LDADD   = $(KK2f_Libs) $(ROOTAUXLIBS) $(ROOTGLIBS)
Prod0_LDFLAGS = -L$(ROOTLIBDIR) -lgfortran

ProdZinv_SOURCES = ProdZinv.cxx KKMC.cxx PartLund.cxx ROBOL2.cxx Semaph.cxx
ProdZinv_LDADD   = $(KK2f_Libs) $(ROOTAUXLIBS) $(ROOTGLIBS)
ProdZinv_LDFLAGS = -L$(ROOTLIBDIR) -lgfortran

ProdFSR_SOURCES = ProdFSR.cxx KKMC.cxx PartLund.cxx RoboFSR.cxx Semaph.cxx
ProdFSR_LDADD   = $(KK2f_Libs) $(ROOTAUXLIBS) $(ROOTGLIBS)
ProdFSR_LDFLAGS = -L$(ROOTLIBDIR) -lgfortran

ProdFoam_SOURCES =  ProdFoam.cxx KKMC.cxx KKfoam.cxx  RoboFoam.cxx Semaph.cxx
ProdFoam_LDADD   = $(KK2f_Libs) $(ROOTAUXLIBS) $(ROOTGLIBS) -lFoam
ProdFoam_LDFLAGS = -L$(ROOTLIBDIR) -lgfortran

########################################################################
###                   USER PART
###  make test0-start
###  make Plot1
#-----------------------------------------------------------------------
test0-start:	all
	(cd ./test0;  make start )
test0-stop:
	(cd ./test0;  make stop )
workAFB-start:	all
	(cd ./workAFB;  make start )
workAFB-stop:
	(cd ./workAFB;  make stop )
Plot1:
	(cd Plots; make Plot1; ./Plot1)
PlotAFB:
	(cd Plots; make PlotAFB; ./PlotAFB)
PlotAFB2:
	(cd Plots; make PlotAFB2; ./PlotAFB2)
PlotIFI:
	(cd Plots; make PlotIFI; ./PlotIFI)
#-----------------------------------------------------------------------
###   make workZinv-start
workZinv-start:	all
	(cd ./workZinv;  make start )
workZinv-stop:
	(cd ./workZinv;  make stop )
Plot2:
	(cd Plots; make Plot2; ./Plot2)
#-----------------------------------------------------------------------
###   make workFSR-start
workFSR-start: all
	(cd ./workFSR; make start)
workFSR-stop:
	(cd ./workFSR; make stop)
PlotFSR:
	(cd Plots; make PlotFSR; ./PlotFSR)
#-----------------------------------------------------------------------
###   make workFOAM-start
workFOAM-start: all
	(cd ./workFOAM; make start)
workFOAM-stop:
	(cd ./workFOAM; make stop)
#-----------------------------------------------------------------------
###   make workFOKK-start
workFOKK-start: all
	(cd ./workFOKK; make start)
workFOKK-stop:
	(cd ./workFOKK; make stop)
#-----------------------------------------------------------------------
