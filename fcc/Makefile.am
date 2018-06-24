# ****** This is input data file for Automake ****** 
# Automake creates Makefile using all information defined below
METASOURCES = AUTO
SUBDIRS = workDEV
############################################
KK2f_Libs= \
	$(top_builddir)/KK2f/libKK2f.a     $(top_builddir)/glibk/libGlibk.a \
	$(top_builddir)/jetset/libJetset.a $(top_builddir)/bornv/libBornvDizet.a \
	$(top_builddir)/tauola/libTauola.a $(top_builddir)/photos/libPhotos.a \
	$(top_builddir)/dizet/libDizet.a \
	$(top_builddir)/KKsem/libKKsem.a
#
INCLUDES = -I$(top_srcdir)/fcc $(all_includes) -I$(ROOTINCDIR) -I$(top_srcdir)/MCdev -I$(top_srcdir)/Foam2
#///////////////////////////////////////////////////////////////
lib_LTLIBRARIES = libDEV.la
CLEANFILES = $(dict1C) $(dict1H)
#///////////////////////////////////////////////////////////////
#   ROOT Dictionary 1 file definition
dict1C = DevDict.cxx
dict1H = DevDict.h
#  Header files:
include_HEADERS1 = TMCgenDEV.h TRobolDEV.h
noinst_HEADERS1 = LinkDef1.h
#  ROOT Directory generation:
ROOTCINTHDRS1= $(include_HEADERS1) $(noinst_HEADERS1)
$(dict1C) : $(ROOTCINTHDRS1)
	$(ROOTCINT) -f $(dict1C) -c $(INCLUDES) $(ROOTCINTHDRS1)
#---------------------------------------------------------------
#  Local library 1
#  Additional flags:
nodist_libDEV_la_SOURCES = $(dict1C)
libDEV_la_SOURCES = TMCgenDEV.cxx TRobolDEV.cxx
libDEV_la_LIBADD  = $(top_builddir)/MCdev/libMCdev.la \
    				$(ROOTAUXLIBS) $(ROOTGLIBS) $(FLIBS)
libDEV_la_LDFLAGS = -L$(ROOTLIBDIR)
#///////////////////////////////////////////////////////////////
# Local executables with ROOT persistency
bin_PROGRAMS = MainDEV
AM_FFLAGS = @MY_F77_FLAGS@
#---------------------------------------------------------------
MainDEV_SOURCES = MainDEV.cxx 
MainDEV_LDADD = $(top_builddir)/fcc/libDEV.la \
				$(ROOTAUXLIBS) $(ROOTGLIBS) $(FLIBS) \
				$(KK2f_Libs)
MainDEV_LDFLAGS = -L$(ROOTLIBDIR) -rpath $(ROOTRPATH) -no-install -g
