PETSC_DIR=/Users/victorsacek/Documents/petsc-3.7.5
PETSC_ARCH=arch-darwin-c-debug

INCFLAGS=-I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include

SOURCEC=main.cpp calc_drho.cpp calc_visc.cpp DM1.cpp DM1_v.cpp DMT.cpp DMv.cpp thermal_Ke.cpp reader.cpp
OBJECTS=$(SOURCEC:%.cpp=%.o)

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
 
all: ${OBJECTS} chkopts
	-${CLINKER} -o MD3D_3.5 ${OBJECTS} ${PETSC_LIB}
	
%.o: %.cpp
	mpicc -Wall -fdiagnostics-color -c $< -o $@ ${INCFLAGS}


