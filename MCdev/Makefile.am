all_includes := -I$(top_srcdir)/Foam2 -I$(top_srcdir)/MCdev $(all_includes)
AM_CPPFLAGS = $(all_includes) $(ROOTCFLAGS) # used while compiling
AM_CXXFLAGS = $(ROOTAUXCFLAGS) # used while compiling and linking

SUBDIRS = farming

# ROOT Dictionary file definition
dictC = libMCdev.cxx
dictH = libMCdev.h
nodist_libMCdev_la_SOURCES = $(dictC)
#Header files:
include_HEADERS = TMCgen.h BXFORMAT.h TSemaf.h TRobol.h TKfig.h TWtMon.h
#Additional flags:
CLEANFILES = $(dictC) $(dictH) $(dictC:.cxx=_rdict.pcm)
#---------------------------------------------------------------
# ROOT Directory generation:
ROOTCINTHDRS= $(include_HEADERS) $(noinst_HEADERS)
$(dictC) : $(ROOTCINTHDRS)
	$(ROOTCINT) -f $(dictC) -c $(all_includes) $(ROOTCINTHDRS)
	if [ -f "$(dictC:.cxx=_rdict.pcm)" -a -d ".libs/" ]; then ln -f -s "../$(dictC:.cxx=_rdict.pcm)" ".libs/"; fi # a brutal fix for ROOT 6 and libtool
#---------------------------------------------------------------
noinst_HEADERS  = LinkDef.h
lib_LTLIBRARIES = libMCdev.la
libMCdev_la_SOURCES = TMCgen.cxx  TSemaf.cxx TRobol.cxx TKfig.cxx  TWtMon.cxx
libMCdev_la_LIBADD =	\
			$(top_builddir)/Foam2/libFOAM.la \
			-L$(ROOTLIBDIR) $(ROOTGLIBS) $(ROOTAUXLIBS)
libMCdev_la_LDFLAGS = #-Wl,-z,defs
#----- explicit target needed for ../ directory
$(top_builddir)/Foam2/libFOAM.la:
	(cd ../Foam2; make libFOAM.la)
#---------------------------------------------------------------
