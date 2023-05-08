FC = gfortran
TARGET = main 
FFLAGS = -Wall -fmax-errors=1 -O3 -fcheck=all #-Wno-unused-variable
OBJ = mod_constant.o mod_spline.o mod_fileIO.o mod_neighbors.o mod_force.o main.o 

${TARGET}: ${OBJ}
	${FC} $^ -o $@ ${FFLAGS}

%.o: %.f90
	${FC} -c $< -o $@ ${FFLAGS}

.PHONY: clean

clean:
	rm -rf *.mod *.o main