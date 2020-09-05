#================================================================
# part interpreted by automake, edited by automake manager
INCLUDES = -I$(top_srcdir)/Foam2 -I$(ROOTINCDIR) $(all_includes)
METASOURCES = AUTO
#---------------------------------------------------------------
# ROOT specific part added by hand
CXXFLAGS = -g -O2 -std=c++11
dictC = FOAMDict.cxx
dictH = FOAMDict.h
nodist_libFOAM_la_SOURCES =$(dictC)
noinst_HEADERS = LinkDef.h
CLEANFILES = $(dictC) $(dictH)
# ROOT Directory generation:
ROOTCINTHDRS= $(include_HEADERS) $(noinst_HEADERS)
$(dictC) : $(ROOTCINTHDRS)
	$(ROOTCINT)  -f $(dictC) -c $(INCLUDES) $(ROOTCINTHDRS)
#---------------------------------------------------------------
# part interpreted by automake, edited by automake manager
include_HEADERS = TFMAXWT.h TFOAM_INTEGRAND.h TFOAM.h ROOT_DEF.h TFHST.h
lib_LTLIBRARIES = libFOAM.la
libFOAM_la_SOURCES = TFMAXWT.cxx TFOAM_INTEGRAND.cxx TFOAM.cxx

