LOCAL_LIB   =  libfig.a 
LIB_OBJECTS =  semalib.o plotlib.o readat.o
##############   What has to be built
all:	$(LOCAL_LIB) 
###############
FFLAGS = -fno-second-underscore -fno-automatic -ffixed-line-length-132 -fbackslash -O
F77 = g77
AR = ar
RANLIB = ranlib
LD = $(F77)
#############
.f.o:
	$(F77) $(FFLAGS) -c $<
#############
semalib.f:	semalib.h
#############
$(LOCAL_LIB): $(LIB_OBJECTS)
	$(AR) rv  $(LOCAL_LIB) $(LIB_OBJECTS)
	$(RANLIB) $(LOCAL_LIB)
clean:
	rm -f *.o; rm -f *.a; rm -f *~ 
Clean:	clean
	rm -f *.o; rm -f *.a; rm -f *~ 
#############################

