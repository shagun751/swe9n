SWITCH_ICC_RELEASE=-ipo -O3 -vec_report1
SWITCH_ICC_CUSTOM=-O3
SWITCH_ICC_DEBUG= -g -debug all -check all -implicitnone -fp-stack-check -heap-arrays -ftrapuv -check pointers -check bounds -traceback
SWITCH_GFORTRAN_RELEASE=-Ofast
SWITCH_GFORTRAN_O2=-O2
SWITCH_GFORTRAN_DEBUG=-O0 -g -fbounds-check -fbacktrace

PARA_FLAGS=-O3 -lstdc++ -fopenmp -L/usr/local/cuda/lib64

MAIN=swe9n

OBJ=$(MAIN).o
	
default: $(OBJ)
	mkdir -p Output Resume
	#gfortran -o swe $(OBJ) $(PARA_FLAGS) -Wall
	gfortran -o swe $(OBJ) $(PARA_FLAGS) 

$(MAIN).o: $(MAIN).f90
	gfortran $(PARA_FLAGS) -c $(MAIN).f90

clean:
	rm -rf *.log *.var *.mod *.dat *.o Output vol.png fort.20
	mkdir Output

rmout:
	rm -rf Output
	mkdir Output
	



