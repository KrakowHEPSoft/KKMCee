AUTOMAKE_OPTIONS = foreign subdir-objects
ACLOCAL_AMFLAGS = -I m4

# all_includes := $(all_includes)
# AM_CPPFLAGS = $(all_includes) -fPIC -Wall -Wno-write-strings # used while compiling
AM_CXXFLAGS = $(ROOTAUXCFLAGS) # used while compiling and linking

#F77=$(FC)
#AM_FFLAGS=$(AM_FCFLAGS) $(FCFLAGS)

lib_LTLIBRARIES = libTauolaPhotos.la

AM_FFLAGS = @MY_F77_FLAGS@

TAUOLA_SRC = \
    formf.f tauola.f curr_cleo.f f3pi.f pkorb.f \
    new-currents/other-currents/frho_pi_belle.f\
    new-currents/RChL-currents/rcht_3pi/f3pi_rcht.f \
    new-currents/RChL-currents/rcht_3pi/funct_3pi.f \
    new-currents/RChL-currents/rcht_common/FA1RCHL.f \
    new-currents/RChL-currents/rcht_common/ffwid3pi.f \
    new-currents/RChL-currents/rcht_common/funct_rpt.f \
    new-currents/RChL-currents/rcht_common/gaus_integr.f \
    new-currents/RChL-currents/rcht_common/gfact.f \
    new-currents/RChL-currents/rcht_common/initA1Tab.f \
    new-currents/RChL-currents/rcht_common/initA1TabKKpi.f \
    new-currents/RChL-currents/rcht_common/value_parameter.f \
    new-currents/RChL-currents/rcht_common/wid_a1_fit.f \
    new-currents/RChL-currents/rcht_common/wid_a1_fitKKpi.f \
    new-currents/other-currents/binp.f \
    new-currents/other-currents/curr_karl.f \
    tauola-c/ChannelForTauolaInterface.cxx \
    tauola-c/channels_wrappers.cxx

libTauolaPhotos_la_SOURCES = \
	$(TAUOLA_SRC) \
	Tauface.f \
	HepEvt.f \
	IO_aux.f
libTauolaPhotos_la_DEPENDENCIES = HepEvt.h

libTauolaPhotos_la_LIBADD = -lgfortran
##	$(top_builddir)/jetset/libJetset.a
libTauolaPhotos_la_LDFLAGS = #-Wl,-z,defs

DISTCLEANFILES = $(BUILT_SOURCES)
