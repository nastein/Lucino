xsec: mathtool.o event_module.o currents_opt_EW_iso_sep.o xsec_fact_intf.o \
  nform_mod.o integrated_diff_xsec_intf.o ew_evgen.o
	gfortran -o $@ $^

%.o: %.f90
	gfortran -o $@ -c $^

.PHONY: clean

clean:
	$(RM) *.mod *.o xsec
