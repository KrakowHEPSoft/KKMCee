AM_CPPFLAGS = $(all_includes) $(ROOTCFLAGS) $(HEPMC_INCLUDES) # used while compiling
AM_CXXFLAGS = $(ROOTAUXCFLAGS) -g -O2 -std=c++11 -lstdc++     # used while compiling and linking
METASOURCES = AUTO

ROOTLDFLAGS = -L$(ROOTLIBDIR) $(ROOTLIBS)

AM_FFLAGS = @MY_F77_FLAGS@


################################################################
######## Program creating EW formfactor tables   ###############
bin_PROGRAMS = basic_tree

#/////////////////////////////////////////////////////////
basic_tree_SOURCES = basic_tree.cc
basic_tree_LDADD   = 	\
					$(ROOTAUXLIBS) $(ROOTGLIBS)
###basic_tree_LDFLAGS = $(ROOTLDFLAGS) -no-install -g
basic_tree_LDFLAGS = $(LDFLAGS_HEPMC) $(ROOTLDFLAGS) -no-install -g
################################################################