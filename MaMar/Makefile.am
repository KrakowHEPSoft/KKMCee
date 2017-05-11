INCLUDES = -I$(top_srcdir)/MaMar $(all_includes) -I$(ROOTINCDIR)
METASOURCES = AUTO
SUBDIRS = test0 workAFB workZinv workFSR workFOAM Plots

bin_PROGRAMS = Prod0 ProdZinv ProdFSR ProdFoam
AM_FFLAGS = @MY_F77_FLAGS@

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

Prod0_SOURCES = Prod0.cxx KKMC.cxx PartLund.cxx ROBOL.cxx Semaph.cxx
Prod0_LDADD   = $(KK2f_Libs) $(ROOTAUXLIBS) $(ROOTGLIBS)
Prod0_LDFLAGS = -L$(ROOTLIBDIR) -lgfortran

ProdZinv_SOURCES = ProdZinv.cxx KKMC.cxx PartLund.cxx ROBOL2.cxx Semaph.cxx
ProdZinv_LDADD   = $(KK2f_Libs) $(ROOTAUXLIBS) $(ROOTGLIBS)
ProdZinv_LDFLAGS = -L$(ROOTLIBDIR) -lgfortran

ProdFSR_SOURCES = ProdFSR.cxx KKMC.cxx PartLund.cxx RoboFSR.cxx Semaph.cxx
ProdFSR_LDADD   = $(KK2f_Libs) $(ROOTAUXLIBS) $(ROOTGLIBS)
ProdFSR_LDFLAGS = -L$(ROOTLIBDIR) -lgfortran

ProdFoam_SOURCES =  KKMC.cxx KKfoam.cxx  RoboFoam.cxx ProdFoam.cxx Semaph.cxx
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
