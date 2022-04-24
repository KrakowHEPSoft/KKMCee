#================================================================
#########INCLUDES = -I$(top_srcdir)/Foam2 -I$(ROOTINCDIR) $(all_includes)
AM_CPPFLAGS = $(all_includes) $(ROOTCFLAGS) # used while compiling
AM_CXXFLAGS = $(ROOTAUXCFLAGS) -g -O2       # used while compiling and linking
#########
METASOURCES = AUTO
#---------------------------------------------------------------
# ROOT specific part added by hand
#CXXFLAGS = -g -O2 -std=c++11
dictC = FOAMDict.cxx
dictH = FOAMDict.h
nodist_libFOAM_la_SOURCES =$(dictC)
noinst_HEADERS = LinkDef.h
CLEANFILES = $(dictC) $(dictH)
# ROOT Directory generation:
ROOTCINTHDRS= $(include_HEADERS) $(noinst_HEADERS)
$(dictC) : $(ROOTCINTHDRS)
	$(ROOTCINT) -f $(dictC) -c $(all_includes) $(ROOTCINTHDRS)
	if [ -f "$(dictH:.h=_rdict.pcm)" -a -d ".libs/" ]; then ln -f -s "../$(dictH:.h=_rdict.pcm)" ".libs/"; fi # a brutal fix for ROOT 6 and libtool
#---------------------------------------------------------------
# part interpreted by automake, edited by automake manager
include_HEADERS = TFMAXWT.h TFOAM_INTEGRAND.h TFOAM.h ROOT_DEF.h TFHST.h
lib_LTLIBRARIES = libFOAM.la
libFOAM_la_SOURCES = TFMAXWT.cxx TFOAM_INTEGRAND.cxx TFOAM.cxx
libFOAM_la_LIBADD = \
	-L$(ROOTLIBDIR) $(ROOTGLIBS) -lFoam $(ROOTAUXLIBS)
libFOAM_la_LDFLAGS = #-Wl,-z,defs
