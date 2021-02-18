SUBDIRS = farming

######INCLUDES = -I$(top_srcdir)/Foam2 -I$(top_srcdir)/MCdev $(all_includes) -I$(ROOTINCDIR)
all_includes :=  -I$(top_srcdir)/Foam2 -I$(top_srcdir)/MCdev $(all_includes)
AM_CPPFLAGS = $(all_includes) $(ROOTCFLAGS)   # used while compiling
AM_CXXFLAGS = $(ROOTAUXCFLAGS) -g -O2         # used while compiling and linking
#######CXXFLAGS = -g -O2 -std=c++11

# ROOT Dictionary file definition
dictC = MCdevDict.cxx
dictH = MCdevDict.h
nodist_libMCdev_la_SOURCES = $(dictC)
#Header files:
include_HEADERS = TMCgen.h BXFORMAT.h TSemaf.h TRobol.h THwtMon.h TKfig.h
#Additional flags:
CLEANFILES = $(dictC) $(dictH)
#---------------------------------------------------------------
# ROOT Directory generation:
ROOTCINTHDRS= $(include_HEADERS) $(noinst_HEADERS)
$(dictC) : $(ROOTCINTHDRS)
	$(ROOTCINT) -f $(dictC) -c $(all_includes) $(ROOTCINTHDRS)
	if [ -f "$(dictH:.h=_rdict.pcm)" -a -d ".libs/" ]; then ln -f -s "../$(dictH:.h=_rdict.pcm)" ".libs/"; fi # a brutal fix for ROOT 6 and libtool
#---------------------------------------------------------------
noinst_HEADERS  = LinkDef.h
lib_LTLIBRARIES = libMCdev.la
libMCdev_la_SOURCES = TMCgen.cxx  TSemaf.cxx TRobol.cxx THwtMon.cxx TKfig.cxx
libMCdev_la_LIBADD = \
	$(top_builddir)/Foam2/libFOAM.la \
	-L$(ROOTLIBDIR) $(ROOTGLIBS) -lFoam $(ROOTAUXLIBS)
#----- explicit target needed for ../ directory
$(top_builddir)/Foam2/libFOAM.la:
	(cd ../Foam2; make libFOAM.la)
#---------------------------------------------------------------
