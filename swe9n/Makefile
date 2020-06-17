SWITCH_ICC_RELEASE=-ipo -O3 -vec_report1
SWITCH_ICC_CUSTOM=-O3
SWITCH_ICC_DEBUG= -g -debug all -check all -implicitnone -fp-stack-check -heap-arrays -ftrapuv -check pointers -check bounds -traceback
SWITCH_GFORTRAN_RELEASE=-Ofast
SWITCH_GFORTRAN_O2=-O2
SWITCH_GFORTRAN_DEBUG=-O0 -g -fbounds-check -fbacktrace

#COMPILER=ifort
#COMPILER=gfortran
COMPILER=pgfortran

ifeq ($(COMPILER),ifort)
	MODPATH=-module subroutines
	C_COMP=icc
	#PARA_FLAGS=-O3 -lstdc++ -fopenmp -L/usr/local/cuda/lib64	
	PARA_FLAGS=-lstdc++ -fopenmp -L/usr/local/cuda/lib64	
else ifeq ($(COMPILER),gfortran)
	MODPATH=-J subroutines
	C_COMP=g++
	#PARA_FLAGS=-O3 -lstdc++ -fopenmp -L/usr/local/cuda/lib64	
	PARA_FLAGS=-lstdc++ -fopenmp -L/usr/local/cuda/lib64	
else ifeq ($(COMPILER),pgfortran)
	MODPATH=-module subroutines
	C_COMP=pgfortran
	#PARA_FLAGS=-O3 -lstdc++ -fopenmp -L/usr/local/cuda/lib64	
	PARA_FLAGS=-fast -Minfo=accel
endif

MAIN=swe9n

OBJ=$(MAIN).f90

SUBS=\
subroutines/*.f90
	
default: $(OBJ) $(SUBS)
	mkdir -p Output
	$(COMPILER) -o swe $(OBJ) $(PARA_FLAGS) $(MODPATH)

cleanAll:
	rm -rf subroutines/*.mod *.o Output swe
	mkdir Output

clean:
	rm -rf Output
	mkdir Output
	