SUBDIRS = workMini work1

all_includes :=  -I$(top_srcdir)/SRCee -I$(top_srcdir)/MCdev -I$(top_srcdir)/Foam2 $(all_includes)
AM_CPPFLAGS = $(all_includes) $(ROOTCFLAGS) # used while compiling
AM_CXXFLAGS = $(ROOTAUXCFLAGS) -g -O2       # used while compiling and linking

#/////////////////////////////////////////////////////////////
# Local executables with ROOT persistency
#-------------------------------------------------------------------------
bin_PROGRAMS = MainMini
#/////////////////////////////////////////////////////////
MainMini_SOURCES = MainMini.cxx 
MainMini_LDADD   = \
			$(top_builddir)/SRCee/libKKee.la \
			$(ROOTAUXLIBS) $(ROOTGLIBS) -lgfortran -lFoam $(FLIBS)
MainMini_LDFLAGS = -L$(ROOTLIBDIR) -rpath $(ROOTRPATH) -no-install -g

#/////////////////////////////////////////////////////////
