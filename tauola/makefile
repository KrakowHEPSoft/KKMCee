include ../make.inc
#######################################################
LOCAL_LIB   =  glib.a 
LIB_OBJECTS =  formf.o tauola.o curr_cleo.o f3pi.o pkorb.o \
               new-currents/other-currents/frho_pi_belle.o \
               new-currents/RChL-currents/rcht_3pi/f3pi_rcht.o \
               new-currents/RChL-currents/rcht_3pi/funct_3pi.o \
               new-currents/RChL-currents/rcht_common/FA1RCHL.o \
               new-currents/RChL-currents/rcht_common/ffwid3pi.o \
               new-currents/RChL-currents/rcht_common/funct_rpt.o \
               new-currents/RChL-currents/rcht_common/gaus_integr.o \
               new-currents/RChL-currents/rcht_common/gfact.o \
               new-currents/RChL-currents/rcht_common/initA1Tab.o \
               new-currents/RChL-currents/rcht_common/initA1TabKKpi.o \
               new-currents/RChL-currents/rcht_common/value_parameter.o \
               new-currents/RChL-currents/rcht_common/wid_a1_fit.o \
               new-currents/RChL-currents/rcht_common/wid_a1_fitKKpi.o \
               new-currents/other-currents/binp.o \
               new-currents/other-currents/curr_karl.o \
               tauola-c/ChannelForTauolaInterface.o \
               tauola-c/channels_wrappers.o

##############   What has to be built
all:	$(LOCAL_LIB) 

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

%.o: %.f
	$(F77) $(FFLAGS) -c $< -o $@
$(LOCAL_LIB): $(LIB_OBJECTS)
	$(AR) rv  $(LOCAL_LIB) $(LIB_OBJECTS)
	$(RANLIB) $(LOCAL_LIB)
clean:
	rm -f *.o; rm -f *.a; rm -f *~ tauola-c/*.o tauola-no-c/*.o
	rm -f new-currents/other-currents/*.o new-currents/RChL-currents/*/*.o
#######################################################

