SUBDIRS = workMini work1

all_includes :=  -I$(top_srcdir)/SRCee -I$(top_srcdir)/MCdev -I$(top_srcdir)/Foam2 $(all_includes)
AM_CPPFLAGS = $(all_includes) $(ROOTCFLAGS) # used while compiling
AM_CXXFLAGS = $(ROOTAUXCFLAGS) -g -O2       # used while compiling and linking

#///////////////////////////////////////////////////////////////
lib_LTLIBRARIES = libKKmc.la
CLEANFILES = $(dict2C) $(dict2H)
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
nodist_libKKmc_la_SOURCES = $(dict2C)
libKKmc_la_SOURCES =  TRobolKKMC.cxx
libKKmc_la_LIBADD =  \
	$(top_builddir)/MCdev/libMCdev.la \
	$(top_builddir)/SRCee/libKKee.la \
    -L$(ROOTLIBDIR) $(ROOTGLIBS) -lFoam $(FLIBS) $(ROOTAUXLIBS)
#///////////////////////////////////////////////////////////////


#/////////////////////////////////////////////////////////////
# Local executables with ROOT persistency
#-------------------------------------------------------------------------
bin_PROGRAMS = MainMini MainKKMC
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
				$(top_builddir)/ProdRun/libKKmc.la
MainKKMC_LDFLAGS = -L$(ROOTLIBDIR) -rpath $(ROOTRPATH) -no-install -g
#-------------------------------------------------------------------------
