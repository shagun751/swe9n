

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
	CUDA_CHK = -Mcuda=ptxinfo
	GPU_FLAG = -ta=tesla:cc60,pinned
	PARA_FLAGS=-fast -acc -Minfo=accel $(CUDA_CHK) $(GPU_FLAG)
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

run:
	PGI_ACC_TIME=1 ./swe sqB

profile:
	nsys profile -t openacc --stats=true --force-overwrite true -o swe ./swe sqB
	