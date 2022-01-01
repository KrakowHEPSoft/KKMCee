AUTOMAKE_OPTIONS = foreign
ACLOCAL_AMFLAGS = -I m4
#F77=$(FC)
#AM_FFLAGS=$(AM_FCFLAGS) $(FCFLAGS)

lib_LTLIBRARIES = libTauolaPhotos.la

AM_FFLAGS = @MY_F77_FLAGS@

####################
# KKMC doesn't have 'make.inc' and does not set C++ flags
# so we have to hard-code them here
CC      = g++
CFLAGS  = -fPIC -Wall -Wno-write-strings
.c.o:
	$(CC) $(CFLAGS) -c $< -o $@
####################

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
    tauola-c/ChannelForTauolaInterface.c \
    tauola-c/channels_wrappers.c

libTauolaPhotos_la_SOURCES = \
	$(TAUOLA_SRC) \
	photos.f \
	Tauface.f Taupair.f \
	KinLib.f \
	HepEvt.f 
libTauolaPhotos_la_DEPENDENCIES = Taupair.h HepEvt.h

libTauolaPhotos_la_LIBADD =  \
	$(top_builddir)/glibk/libGlibk.a \
	$(top_builddir)/jetset/libJetset.a \
	-lgfortran

DISTCLEANFILES = $(BUILT_SOURCES)
