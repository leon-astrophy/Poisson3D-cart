# Define F90 compiler #
F90=/home/leon/HDF5/bin/h5fc

# Choose between CPU or GPU threading #
MULTICORE = -DCPU -mp
#MULTICORE = -DGPU -acc

# Define F90 flags #
F90FLAGS= -cpp $(MULTICORE)
LDFLAGS= -cpp $(MULTICORE)
