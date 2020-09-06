
METASOURCES = AUTO
bin_PROGRAMS = yfspro
yfspro_SOURCES = yfspro.f
yfspro_LDADD = $(top_builddir)/KK2f/libKK2f.a $(top_builddir)/glibk/libGlibk.a\
	$(top_builddir)/jetset/libJetset.a $(top_builddir)/bornv/libBornvDizet.a\
	$(top_builddir)/tauola/libTauola.a $(top_builddir)/photos/libPhotos.a\
	$(top_builddir)/dizet/libDizet.a -lstdc++
SUBDIRS = libfig pubfig

########################################################################
###                user part
###=====================================================================
#  make makprod          produces makefiles in work dirs from template
#  make E189GeV-start
#  make E189GeV-stop
########################################################################
MCOMMAND = yfspro
#-----------------------------------------------------------------------
E500GeV-start:	$(MCOMMAND)
	(cd ./E500GeV;  make start )
E500GeV-stop:
	(cd ./E500GeV;  make stop )
#-----------------------------------------------------------------------
E189GeV-start:	$(MCOMMAND)
	(cd ./E189GeV;  make start )
E189GeV-batch:	$(MCOMMAND)
	(cd ./E189GeV;  make batch )
E189GeV-stop:
	(cd ./E189GeV;  make stop )
E189GeV-cont:
	(cd ./E189GeV;  make cont )
E189GeV-cont-debug:	$(MCOMMAND)
	(cp E189GeV/pro.hst.core       E189GeV/pro.hst )
	(cp E189GeV/semaphore.core     E189GeV/semaphore )
	(cp E189GeV/pro.input.core     E189GeV/pro.input )
###	(cd ./E189GeV;  make cont )
	(cd ./E189GeV;  make cont-debug )
#-----------------------------------------------------------------------
E120GeV-start:	$(MCOMMAND)
	(cd ./E120GeV;  make start )
E120GeV-stop:
	(cd ./E120GeV;  make stop )
E120GeV-cont:
	(cd ./E120GeV;  make cont )
#-----------------------------------------------------------------------
E91GeV-start:	$(MCOMMAND)
	(cd ./E91GeV;  make start )
E91GeV-stop:
	(cd ./E91GeV;  make stop )
E91GeV-cont:
	(cd ./E91GeV;  make cont )
#-----------------------------------------------------------------------
ini140-start:	$(MCOMMAND)
	(cd ./ini140;  make start )
ini140-stop:
	(cd ./ini140;  make stop )
ini140-plot:
	(cd ./profig; make figini-dvi)
ini140-cont-debug:	$(MCOMMAND)
	(cp ini140/pro.hst.75396 ini140/pro.hst )
	(cp ini140/semaphore.75396 ini140/semaphore )
	(cd ./ini140;  make cont )
#-----------------------------------------------------------------------
ini200-start:	$(MCOMMAND)
	(cd ./ini200;  make start )
ini200-stop:
	(cd ./ini200;  make stop )
ini200-cont:
	(cd ./ini200;  make cont )
#-----------------------------------------------------------------------
fin140-start:	$(MCOMMAND)
	(cd ./fin140;  make start )
fin140-stop:
	(cd ./fin140;  make stop )
#-----------------------------------------------------------------------
mix10-start:	$(MCOMMAND)
	(cd ./mix10;  make start )
mix10-stop:
	(cd ./mix10;  make stop )
#-----------------------------------------------------------------------
mix89-start:	$(MCOMMAND)
	(cd ./mix89;  make start )
mix89-stop:
	(cd ./mix89;  make stop )
#-----------------------------------------------------------------------
mix91-start:	$(MCOMMAND)
	(cd ./mix91;  make start )
mix91-stop:
	(cd ./mix91;  make stop )
#-----------------------------------------------------------------------
mix93-start:	$(MCOMMAND)
	(cd ./mix93;  make start )
mix93-stop:
	(cd ./mix93;  make stop )
#-----------------------------------------------------------------------
mix200-start:	$(MCOMMAND)
	(cd ./mix200;  make start )
mix200-stop:
	(cd ./mix200;  make stop )
mix200-cont:
	(cd ./mix200;  make cont )
mix200-plot:
	(cd ./profig; rm -f figmix.dvi; make figmix-ps)
mix200-start-debug:	$(MCOMMAND)
	(cd ./mix200;  make start-debug )
mix200-cont-debug:	$(MCOMMAND)
	(cp mix200/pro.hst.382 mix200/pro.hst )
	(cp mix200/semaphore.382 mix200/semaphore )
	(cd ./mix200;  make cont )
###	(cd ./mix200;  make cont-debug )
#-----------------------------------------------------------------------
mix1000-start:	$(MCOMMAND)
	(cd ./mix1000;  make start )
mix1000-stop:
	(cd ./mix1000;  make stop )
#=======================================================================
# Old historical
ini40-start:	$(MCOMMAND)
	(cd ./ini40;  make start )
ini40-stop:
	(cd ./ini40;  make stop )
ini40-plot:
	(cd ./profig; make ini40)
#----------------
# Old historical
fin40-start:	$(MCOMMAND)
	(cd ./fin40;  make start )
fin40-stop:
	(cd ./fin40;  make stop )
#=======================================================================
figini:
	(cd ./profig; rm -f figini.dvi; make figini-ps)
figfin:
	(cd ./profig; rm -f figfin.dvi; make figfin-ps)
figbet:
	(cd ./profig; rm -f figbet.dvi; make figbet-ps)
figmix:
	(cd ./profig; rm -f figmix.dvi; make figmix-ps)
figmix-debug:
	(cd ./profig; rm -f figmix.dvi; make figmix-debug)
#-----------------------------------------------------------------------
#                  update production makefiles
#-----------------------------------------------------------------------
makprod:
	(sed -e 's#<main>#yfspro#' -e 's#<dset>#E500GeV#' -e 's#<class>#8nh#' ./makefile.production > ./E500GeV/makefile)
	(sed -e 's#<main>#yfspro#' -e 's#<dset>#E189GeV#' -e 's#<class>#8nh#' ./makefile.production > ./E189GeV/makefile)
	(sed -e 's#<main>#yfspro#' -e 's#<dset>#E120GeV#' -e 's#<class>#8nh#' ./makefile.production > ./E120GeV/makefile)
	(sed -e 's#<main>#yfspro#' -e 's#<dset>#E91GeV#'  -e 's#<class>#8nh#' ./makefile.production > ./E91GeV/makefile)
	(sed -e 's#<main>#yfspro#' -e 's#<dset>#ini200#'  -e 's#<class>#8nh#' ./makefile.production > ./ini200/makefile)
	(sed -e 's#<main>#yfspro#' -e 's#<dset>#ini140#'  -e 's#<class>#8nh#' ./makefile.production > ./ini140/makefile)
	(sed -e 's#<main>#yfspro#' -e 's#<dset>#fin140#'  -e 's#<class>#8nh#' ./makefile.production > ./fin140/makefile)
#----------------------------------- oldies
	(sed -e 's#<main>#yfspro#' -e 's#<dset>#ini40#'  -e 's#<class>#8nh#' ./makefile.production > ./ini40/makefile)
	(sed -e 's#<main>#yfspro#' -e 's#<dset>#fin40#'  -e 's#<class>#8nh#' ./makefile.production > ./fin40/makefile)
	(sed -e 's#<main>#yfspro#' -e 's#<dset>#mix10#'  -e 's#<class>#8nh#' ./makefile.production > ./mix10/makefile)
	(sed -e 's#<main>#yfspro#' -e 's#<dset>#mix91#'  -e 's#<class>#8nh#' ./makefile.production > ./mix91/makefile)
	(sed -e 's#<main>#yfspro#' -e 's#<dset>#mix200#' -e 's#<class>#8nh#' ./makefile.production > ./mix200/makefile)
#----------------------------------- end oldies ----------------------
AM_FFLAGS = @MY_F77_FLAGS@
