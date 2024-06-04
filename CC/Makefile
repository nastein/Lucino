xsec: mathtool.o event_module.o currents_opt_v2.o xsec_fact.o \
  nform_zexp.o integrated_diff_xsec.o ew_evgen.o
	gfortran -o $@ $^

%.o: %.f90
	gfortran -o $@ -c $^

.PHONY: clean

clean:
	$(RM) *.mod *.o xsec
