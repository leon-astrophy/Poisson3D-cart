#include header#
include make.inc

#specify directories#
SOURCE := .

#source list#
SOURCE1 = $(shell cat $(SOURCE)/source.list)
SRC1 = $(notdir $(SOURCE1))
OBJ1 = $(patsubst %.f90,%.o,$(SRC1))

ALL: $(OBJ1)
	echo $(SRC1)
	$(F90) $(LDFLAGS) -o solver $(OBJ1)

$(OBJ1): %.o: $(SOURCE)/%.f90 
	$(F90) $(F90FLAGS) -c $< -o $@

clean:

	rm -rf solver
	rm -rf *.o
	rm -rf *.mod
	rm tmp.txt

cleanfile:

	rm -rf *.dat
	rm -rf *.hdf5