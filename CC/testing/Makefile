test: test_formfacts.o nform_zexp.o 
	gfortran -o $@ $^

#test_spinors: test_spinors.o currents_opt_v1.o 
#	gfortran -o $@ $^

%.o: %.f90
	gfortran -o $@ -c $^

.PHONY: clean

clean:
	$(RM) *.mod *.o xsec
