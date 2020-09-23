.PHONY: default help objects executable all clean

FC = gfortran

SOURCE_F90 = $(wildcard *.f90)
OBJECTS_F90 = $(patsubst %.f90, %.o, $(SOURCE_F90))
EXECUTABLE = main.out 

default: all
objects: $(OBJECTS_F90)
executable: $(EXECUTABLE)
all: objects executable

help:
	@echo "\
Options:\n\n\
  make objects:       compiler makes objects for every *.c\n\
  make executable:    compiler makes executable\n\
  make all:           build all previous\n\
  make clean:         delete output files\n\
  make help:          display this help"

%.o: %.f90 
	$(FC) -c -o $@ $^

$(EXECUTABLE): $(OBJECTS_F90)
	$(FC) -o $@ $^
clean:
	rm $(OBJECTS_F90) $(EXECUTABLE)
