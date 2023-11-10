xsec: currents_opt_v1.o mathtool.o xsec_fact.o main_diff_xsec.o \
  nform_zexp.o
	gfortran -o $@ $^

%.o: %.f90
	gfortran -o $@ -c $^

.PHONY: clean

clean:
	$(RM) *.mod *.o xsec
