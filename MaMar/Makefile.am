# ****** This is input data file for Automake ****** 
# Automake creates Makefile using all information defined below
INCLUDES = -I$(top_srcdir)/MaMar $(all_includes) -I$(ROOTINCDIR) -I$(top_srcdir)/MCdev -I$(top_srcdir)/Foam2
METASOURCES = AUTO
SUBDIRS = test0 workAFB workZinv workFSR workFOAM workFoam0 workKKMC Plots
#///////////////////////////////////////////////////////////////
lib_LTLIBRARIES = libKKfm.la libKKmc.la
CLEANFILES = $(dict1C) $(dict1H) $(dict2C) $(dict2H)
#///////////////////////////////////////////////////////////////
#   ROOT Dictionary 1 file definition
dict1C = KKfmDict.cxx
dict1H = KKfmDict.h
#  Header files:
include_HEADERS1 = TMCgenFOAM.h TRobolFOAM.h
noinst_HEADERS1 = LinkDef1.h
#  ROOT Directory generation:
ROOTCINTHDRS1= $(include_HEADERS1) $(noinst_HEADERS1)
$(dict1C) : $(ROOTCINTHDRS1)
	$(ROOTCINT) -f $(dict1C) -c $(INCLUDES) $(ROOTCINTHDRS1)
#---------------------------------------------------------------
#  Local library 1
#  Additional flags:
nodist_libKKfm_la_SOURCES = $(dict1C)
libKKfm_la_SOURCES =  TMCgenFOAM.cxx TRobolFOAM.cxx
libKKfm_la_LIBADD =  $(top_builddir)/MCdev/libMCdev.la \
    $(ROOTAUXLIBS) $(ROOTGLIBS) -lFoam $(FLIBS)
libKKfm_la_LDFLAGS = -L$(ROOTLIBDIR)
#///////////////////////////////////////////////////////////////
#   ROOT Dictionary 2 file definition
dict2C = KKmcDict.cxx
dict2H = KKmcDict.h
#  Header files:
include_HEADERS2 = TMCgenKKMC.h TRobolKKMC.h TPartLund.h
noinst_HEADERS2 = LinkDef2.h
#  ROOT Directory generation:
ROOTCINTHDRS2= $(include_HEADERS2) $(noinst_HEADERS2)
$(dict2C) : $(ROOTCINTHDRS2)
	$(ROOTCINT) -f $(dict2C) -c $(INCLUDES) $(ROOTCINTHDRS2)
#---------------------------------------------------------------
#  Local library 2
#  Additional flags:
nodist_libKKmc_la_SOURCES = $(dict2C)
libKKmc_la_SOURCES =  TMCgenKKMC.cxx TRobolKKMC.cxx TPartLund.cxx
libKKmc_la_LIBADD =  $(top_builddir)/MCdev/libMCdev.la \
    $(ROOTAUXLIBS) $(ROOTGLIBS) -lFoam $(FLIBS)
libKKmc_la_LDFLAGS = -L$(ROOTLIBDIR)
#///////////////////////////////////////////////////////////////
#  explicit target in ../
$(top_builddir)/MCdev/libMCdev.la: $(top_builddir)/MCdev/*.h $(top_builddir)/MCdev/*.cxx
	(cd $(top_builddir)/MCdev/; make)
#///////////////////////////////////////////////////////////////

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
############### extras ##################
$(top_builddir)/KKsem/libKKsem.a:
	(cd $(top_builddir)/KKsem/; make)
############################################
KK2f_Libs= \
	$(top_builddir)/KK2f/libKK2f.a     $(top_builddir)/glibk/libGlibk.a \
	$(top_builddir)/jetset/libJetset.a $(top_builddir)/bornv/libBornvDizet.a \
	$(top_builddir)/tauola/libTauola.a $(top_builddir)/photos/libPhotos.a \
	$(top_builddir)/dizet/libDizet.a \
	$(top_builddir)/KKsem/libKKsem.a
#/////////////////////////////////////////////////////////////////////////
# Local executables with ROOT persistency
bin_PROGRAMS = Prod0 ProdZinv ProdFSR ProdFoam0 MainFOAM MainKKMC
AM_FFLAGS = @MY_F77_FLAGS@
#-------------------------------------------------------------------------
MainFOAM_SOURCES = MainFOAM.cxx 
MainFOAM_LDADD = $(top_builddir)/MCdev/.libs/libMCdev.a \
				 $(top_builddir)/Foam2/.libs/libFOAM.a  \
    			 $(ROOTAUXLIBS) $(ROOTGLIBS) -lFoam $(FLIBS) \
    			 $(top_builddir)/MaMar/libKKfm.la \
    			 $(KK2f_Libs)
MainFOAM_LDFLAGS = -L$(ROOTLIBDIR) -rpath $(ROOTRPATH) -no-install -g
#-------------------------------------------------------------------------
MainKKMC_SOURCES = MainKKMC.cxx 
MainKKMC_LDADD = $(top_builddir)/MCdev/.libs/libMCdev.a \
				 $(top_builddir)/Foam2/.libs/libFOAM.a  \
    			 $(ROOTAUXLIBS) $(ROOTGLIBS) -lFoam $(FLIBS) \
    			 $(top_builddir)/MaMar/libKKmc.la \
    			 $(KK2f_Libs)
MainKKMC_LDFLAGS = -L$(ROOTLIBDIR) -rpath $(ROOTRPATH) -no-install -g
#-------------------------------------------------------------------------
#/////////////////////////////////////////////////////////////////////////
# Local executables witout ROOT persistency
Prod0_SOURCES = Prod0.cxx KKMC.cxx PartLund.cxx ROBOL.cxx Semaph.cxx
Prod0_LDADD   = $(KK2f_Libs) $(ROOTAUXLIBS) $(ROOTGLIBS)
Prod0_LDFLAGS = -L$(ROOTLIBDIR) -lgfortran

ProdZinv_SOURCES = ProdZinv.cxx KKMC.cxx PartLund.cxx ROBOL2.cxx Semaph.cxx
ProdZinv_LDADD   = $(KK2f_Libs) $(ROOTAUXLIBS) $(ROOTGLIBS)
ProdZinv_LDFLAGS = -L$(ROOTLIBDIR) -lgfortran

ProdFSR_SOURCES = ProdFSR.cxx KKMC.cxx PartLund.cxx RoboFSR.cxx Semaph.cxx
ProdFSR_LDADD   = $(KK2f_Libs) $(ROOTAUXLIBS) $(ROOTGLIBS)
ProdFSR_LDFLAGS = -L$(ROOTLIBDIR) -lgfortran

ProdFoam0_SOURCES =  ProdFoam0.cxx KKMC.cxx KKfoam.cxx  RoboFoam.cxx Semaph.cxx
ProdFoam0_LDADD   = $(KK2f_Libs) $(ROOTAUXLIBS) $(ROOTGLIBS) -lFoam
ProdFoam0_LDFLAGS = -L$(ROOTLIBDIR) -lgfortran

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
###   make workFoam0-start
workFoam0-start: all
	(cd ./workFoam0; make start)
workFoam0-stop:
	(cd ./workFoam0; make stop)
#-----------------------------------------------------------------------
###   make workFOAM-start
workFOAM-start: all
	(cd ./workFOAM; make start)
workFOAM-stop:
	(cd ./workFOAM; make stop)
#-----------------------------------------------------------------------
###   make workKKMC-start
workKKMC-start: all
	(cd ./workKKMC; make start)
workKKMC-stop:
	(cd ./workKKMC; make stop)
##-----------------------------------------------------------------------
