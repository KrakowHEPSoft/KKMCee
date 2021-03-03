AM_CPPFLAGS = $(all_includes) $(ROOTCFLAGS) # used while compiling
AM_CXXFLAGS = $(ROOTAUXCFLAGS) -g -O2       # used while compiling and linking
METASOURCES = AUTO

ROOTLDFLAGS = -L$(ROOTLIBDIR) $(ROOTLIBS)

AM_FFLAGS = @MY_F77_FLAGS@

lib_LIBRARIES = libBornv.a

libBornv_a_SOURCES = BornV.f hhDizet.f IO_aux.f

################################################################
######## Program creating EW formfactor tables   ###############
bin_PROGRAMS = TabMainC

$(top_builddir)/dizet/libDizet2.a:
	(cd ../dizet; make)

#/////////////////////////////////////////////////////////
TabMainC_SOURCES = TabMainC.cxx
TabMainC_LDADD   = 	\
					$(top_builddir)/DZface/libBornv.a \
					$(top_builddir)/dizet/libDizet2.a \
					$(ROOTAUXLIBS) $(ROOTGLIBS) -lgfortran -lFoam $(FLIBS)
TabMainC_LDFLAGS = $(ROOTLDFLAGS) -no-install -g
################################################################