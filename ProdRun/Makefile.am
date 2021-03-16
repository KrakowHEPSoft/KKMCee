SUBDIRS = workMini work1 workFoam

all_includes :=  -I$(top_srcdir)/SRCee -I$(top_srcdir)/MCdev -I$(top_srcdir)/Foam2 $(all_includes)
AM_CPPFLAGS = $(all_includes) $(ROOTCFLAGS) # used while compiling
AM_CXXFLAGS = $(ROOTAUXCFLAGS) -g -O2       # used while compiling and linking

#///////////////////////////////////////////////////////////////
lib_LTLIBRARIES = libProdKK.la  libProdFoam.la

CLEANFILES = $(dict2C) $(dict2H)

#################################################################[[[
#   ROOT Dictionary for local library 1
dict1C = ProdFoamDict.cxx
dict1H = ProdFoamDict.h
CLEANFILES = $(dict1C) $(dict1H)
#  Header files:
include_HEADERS1 = TRobolFoam.h
noinst_HEADERS1 = LinkDef1.h
#  ROOT Directory generation:
ROOTCINTHDRS1= $(include_HEADERS1) $(noinst_HEADERS1)
$(dict1C) : $(ROOTCINTHDRS1)
	$(ROOTCINT) -f $(dict1C) -c $(all_includes)  $(ROOTCINTHDRS1)
	if [ -f "$(dict1H:.h=_rdict.pcm)" -a -d ".libs/" ]; then ln -f -s "../$(dict1H:.h=_rdict.pcm)" ".libs/"; fi # a brutal fix for ROOT 6 and libtool
##-----------------------Local library 1--------------------------
nodist_libProdFoam_la_SOURCES = $(dict1C)
libProdFoam_la_SOURCES =  TRobolFoam.cxx
libProdFoam_la_LIBADD =  \
    $(top_builddir)/MCdev/libMCdev.la \
    $(ROOTAUXLIBS) $(ROOTGLIBS) -lFoam $(FLIBS)
libProdFoam_la_LDFLAGS = -L$(ROOTLIBDIR)
#################################################################]]]

#///////////////////////////////////////////////////////////////
#   ROOT Dictionary 2 file definition, Local library 2
dict2C = KKmcDict.cxx
dict2H = KKmcDict.h
#  Header files:
include_HEADERS2 = TRobolKKMC.h
noinst_HEADERS2 = LinkDef2.h
#  ROOT Directory generation:
ROOTCINTHDRS2= $(include_HEADERS2) $(noinst_HEADERS2)
$(dict2C) : $(ROOTCINTHDRS2)
	$(ROOTCINT) -f $(dict2C) -c $(all_includes) $(ROOTCINTHDRS2)
	if [ -f "$(dict2H:.h=_rdict.pcm)" -a -d ".libs/" ]; then ln -f -s "../$(dict2H:.h=_rdict.pcm)" ".libs/"; fi # a brutal fix for ROOT 6 and libtool
#---------------------------------------------------------------
#  Additional flags:
nodist_libProdKK_la_SOURCES = $(dict2C)
libProdKK_la_SOURCES =  TRobolKKMC.cxx
libProdKK_la_LIBADD =  \
	$(top_builddir)/MCdev/libMCdev.la \
	$(top_builddir)/SRCee/libKKee.la \
    -L$(ROOTLIBDIR) $(ROOTGLIBS) -lFoam $(FLIBS) $(ROOTAUXLIBS)
#///////////////////////////////////////////////////////////////

#///////////////////////////////////////////////////////////////
#  explicit target in ../MCdev and ../SRCee
$(top_builddir)/MCdev/libMCdev.la: $(top_builddir)/MCdev/*.h $(top_builddir)/MCdev/*.cxx
	(cd $(top_builddir)/MCdev/; make)
$(top_builddir)/SRCee/libKKee.la:  $(top_builddir)/SRCee/*.h $(top_builddir)/SRCee/*.cxx
	(cd $(top_builddir)/SRCee/; make)

#/////////////////////////////////////////////////////////////
# Local executables with ROOT persistency
bin_PROGRAMS = MainMini MainKKMC  MainFoam
#/////////////////////////////////////////////////////////
MainMini_SOURCES = MainMini.cxx
MainMini_LDADD   = \
			$(top_builddir)/SRCee/libKKee.la \
			$(ROOTAUXLIBS) $(ROOTGLIBS) -lgfortran -lFoam $(FLIBS)
MainMini_LDFLAGS = -L$(ROOTLIBDIR) -rpath $(ROOTRPATH) -no-install -g
#/////////////////////////////////////////////////////////
#-------------------------------------------------------------------------
MainKKMC_SOURCES = MainKKMC.cxx 
MainKKMC_LDADD = \
				$(top_builddir)/MCdev/.libs/libMCdev.a \
				$(top_builddir)/Foam2/.libs/libFOAM.a  \
				$(ROOTAUXLIBS) $(ROOTGLIBS) -lFoam $(FLIBS) \
				$(top_builddir)/ProdRun/libProdKK.la
MainKKMC_LDFLAGS = -L$(ROOTLIBDIR) -rpath $(ROOTRPATH) -no-install -g
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
MainFoam_SOURCES = MainKKMC.cxx 
MainFoam_LDADD = \
				$(top_builddir)/MCdev/.libs/libMCdev.a \
				$(top_builddir)/Foam2/.libs/libFOAM.a  \
				$(ROOTAUXLIBS) $(ROOTGLIBS) -lFoam $(FLIBS) \
				$(top_builddir)/ProdRun/libProdFoam.la
MainFoam_LDFLAGS = -L$(ROOTLIBDIR) -rpath $(ROOTRPATH) -no-install -g
#-------------------------------------------------------------------------
