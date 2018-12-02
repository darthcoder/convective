# ~/Documents/Comprehensive/source/Makefile

FC = gfortran
FFLAGS = -O3 -Wall
LFLAGS = -lblas -llapack
.PHONY: clean test

%.o : %.f90
	$(FC) $(FFLAGS) $(LFLAGS) -c $<

compact.exe: compact.o 
	$(FC) compact.o $(LFLAGS) -o compact.exe

test: compact.exe
	./compact.exe
clean:
	rm -f *.o *.exe
