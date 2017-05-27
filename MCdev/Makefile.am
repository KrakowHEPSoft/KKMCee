
SUBDIRS = farming

# ROOT Dictionary file definition
dictC = MCdevDict.cxx
dictH = MCdevDict.h
nodist_libMCdev_la_SOURCES = $(dictC)
#Header files:
include_HEADERS = TMCgen.h BXFORMAT.h TSemaf.h TRobol.h THwtMon.h TKfig.h
#Additional flags:
INCLUDES = -I$(top_srcdir)/Foam2 -I$(top_srcdir)/MCdev $(all_includes) \
	-I$(ROOTINCDIR)
CLEANFILES = $(dictC) $(dictH)
#---------------------------------------------------------------
# ROOT Directory generation:
ROOTCINTHDRS= $(include_HEADERS) $(noinst_HEADERS)
$(dictC) : $(ROOTCINTHDRS)
	$(ROOTCINT) -f $(dictC) -c $(INCLUDES) $(ROOTCINTHDRS)
#---------------------------------------------------------------
noinst_HEADERS  = LinkDef.h
lib_LTLIBRARIES = libMCdev.la
libMCdev_la_SOURCES = TMCgen.cxx  TSemaf.cxx TRobol.cxx THwtMon.cxx TKfig.cxx
libMCdev_la_LIBADD = $(top_builddir)/Foam2/libFOAM.la $(ROOTAUXLIBS) \
	$(ROOTGLIBS) -lFoam
libMCdev_la_LDFLAGS = -L$(ROOTLIBDIR)
#----- explicit target needed for ../ directory
$(top_builddir)/Foam2/libFOAM.la:
	(cd ../Foam2; make libFOAM.la)
#---------------------------------------------------------------
