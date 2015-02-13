FC	= gfortran
FFLAGS	= -O3  /usr/lib/crt1.o /usr/lib/libc.dylib
OBJS	= cosmo.o da.o growth.o linearpk.o spline.o hunt.o dvert.o dverk.o  rombint.o 
.SUFFIXES: .f90
.f90.o: 
	$(FC) $(FFLAGS) -c $<
default:  fisher_distance_bao_w0_wa_marginlized_over_Omegam_Omegak_h
fisher_distance_bao_w0_wa_marginlized_over_Omegam_Omegak_h: $(OBJS) fisher_distance_bao_w0_wa_marginlized_over_Omegam_Omegak_h.o
	$(FC) $(FFLAGS) -o $@ $(OBJS) $@.o $(LDFLAGS)
clean:
	-rm -f *.o *.mod
tidy: clean
	-rm -f fisher_distance_bao_w0_wa_marginlized_over_Omegam_Omegak_h
	-rm -f number.txt
