FC=gfortran
FCFLAGS=-std=f2008 -fcheck=all -O2

all : incomp_solver

linalg.o : linalg.f95
	$(FC) $(FCFLAGS) -c $<

incomp.o : incomp.f95 linalg.o
	$(FC) $(FCFLAGS) -c $<

solver.o : solver.f95 
	$(FC) $(FCFLAGS) -c $<

incomp_solver : incomp.o solver.o linalg.o
	$(FC) $(FCFLAGS) -o $@ $^

clean :
	rm -rf *.o *.tec *.dat *.out *.mod
