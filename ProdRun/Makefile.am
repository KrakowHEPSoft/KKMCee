SUBDIRS = workMini work1 workFoam workNU workTau

all_includes :=  -I$(top_srcdir)/SRCee -I$(top_srcdir)/MCdev -I$(top_srcdir)/Foam2 $(HEPMC_INCLUDES) $(PHOTOS_CPPFLAGS) $(all_includes)
AM_CPPFLAGS = $(all_includes) $(ROOTCFLAGS) # used while compiling
AM_CXXFLAGS = $(ROOTAUXCFLAGS) -g -O2       # used while compiling and linking

#/////////////////////////////////////////////////////////////////////////
# -lFoam is optional in case FOAM of ROOT is used instead of local Foam2
#LIBADD_RTF = -L$(ROOTLIBDIR) $(ROOTGLIBS) -lFoam $(ROOTAUXLIBS) $(FLIBS)
#/////////////////////////////////////////////////////////////////////////
lib_LTLIBRARIES = libProdKK.la  libProdFoam.la

CLEANFILES = $(dict1C) $(dict1H) $(dict2C) $(dict2H)
#################################################################[[[
#   ROOT Dictionary for local library 1
dict1C = ProdFoamDict.cxx
dict1H = ProdFoamDict.h
#  Header files:
include_HEADERS1 = TRobolFoam.h
noinst_HEADERS1 = LinkDef1.h
#  ROOT Directory generation:
ROOTCINTHDRS1= $(include_HEADERS1) $(noinst_HEADERS1)
$(dict1C) : $(ROOTCINTHDRS1)
	$(ROOTCINT) -f $(dict1C) -c $(all_includes)  $(ROOTCINTHDRS1)
	if [ -f "$(dict1H:.h=_rdict.pcm)" -a -d ".libs/" ]; then ln -f -s "../$(dict1H:.h=_rdict.pcm)" ".libs/"; fi # a brutal fix for ROOT 6 and libtool
##-----------------------Local library 1--------------------------
#  Additional flags:
nodist_libProdFoam_la_SOURCES = $(dict1C)
libProdFoam_la_SOURCES =  TRobolFoam.cxx
libProdFoam_la_LIBADD =  \
	$(top_builddir)/SRCee/libKKee.la \
	$(top_builddir)/SRCee/libKKfm.la \
	$(LDFLAGS_HEPMC)
libProdFoam_la_LDFLAGS = #-Wl,-z,defs
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
#------------------------Local library 2---------------------------
#  Additional flags:
nodist_libProdKK_la_SOURCES = $(dict2C)
libProdKK_la_SOURCES =  TRobolKKMC.cxx
libProdKK_la_LIBADD =  \
	$(top_builddir)/SRCee/libKKee.la \
	$(PHOTOS_LDFLAGS) $(LDFLAGS_HEPMC)
libProdKK_la_LDFLAGS = -R $(PHOTOS_LIBDIR) #-Wl,-z,defs
#///////////////////////////////////////////////////////////////

#/////////////////////////////////////////////////////////////
# Local executables with ROOT persistency
bin_PROGRAMS = MainKKMC MainFoam MainMini MainPers
#/////////////////////////////////////////////////////////
#-------------------------------------------------------------------------
# option -no-install needed for debuger
#-------------------------------------------------------------------------
MainKKMC_SOURCES = MainMCdev.cxx 
MainKKMC_LDADD = \
			libProdKK.la
MainKKMC_LDFLAGS = -R $(PHOTOS_LIBDIR) $(PHOTOS_LDFLAGS) $(LDFLAGS_HEPMC) -no-install -g
#-------------------------------------------------------------------------
# option -no-install needed for debuger
#-------------------------------------------------------------------------
MainFoam_SOURCES = MainMCdev.cxx 
MainFoam_LDADD = \
			libProdFoam.la
MainFoam_LDFLAGS = -R $(PHOTOS_LIBDIR) $(PHOTOS_LDFLAGS) $(LDFLAGS_HEPMC) -no-install -g
#/////////////////////////////////////////////////////////
MainMini_SOURCES = MainMini.cxx
MainMini_LDADD   = \
			$(top_builddir)/SRCee/libKKee.la
MainMini_LDFLAGS = -R $(PHOTOS_LIBDIR) $(PHOTOS_LDFLAGS) $(LDFLAGS_HEPMC) -no-install -g
#/////////////////////////////////////////////////////////
MainPers_SOURCES = MainPers.cxx
MainPers_LDADD   = \
			$(top_builddir)/SRCee/libKKee.la
MainPers_LDFLAGS = -R $(PHOTOS_LIBDIR) $(PHOTOS_LDFLAGS) $(LDFLAGS_HEPMC) -no-install -g
#-------------------------------------------------------------------------
