AM_CPPFLAGS = $(all_includes) $(ROOTCFLAGS) # used while compiling
AM_CXXFLAGS = $(ROOTAUXCFLAGS) -g -O2       # used while compiling and linking
METASOURCES = AUTO

ROOTLDFLAGS = -L$(ROOTLIBDIR) $(ROOTLIBS)

AM_FFLAGS = @MY_F77_FLAGS@

###### libDZface.a is used only locally here #####
lib_LIBRARIES = libDZface.a

libDZface_a_SOURCES = hhDizet.f

################################################################
######## Program creating EW formfactor tables   ###############
bin_PROGRAMS = TabMainC

$(top_builddir)/dizet/libDizet2.a:
	(cd ../dizet; make)

#/////////////////////////////////////////////////////////
TabMainC_SOURCES = TabMainC.cxx
TabMainC_LDADD   = 	\
					libDZface.a \
					$(top_builddir)/dizet/libDizet2.a \
					$(ROOTAUXLIBS) $(ROOTGLIBS) -lFoam $(FLIBS)
TabMainC_LDFLAGS = $(ROOTLDFLAGS) -no-install -g
################################################################
