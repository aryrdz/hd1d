#---------------------------------------------------
#   for Coral-1D
#---------------------------------------------------
#   Name of the executable
PROGRAM=hd1d
#
#COMPILER= mpif77
#COMPILER= ifort
COMPILER= gfortran
#
FLAGS=-fdefault-double-8
#
#####################################################
# There should be no need to modify below this line #
#####################################################
OBJECTS = constants.o \
           globals.o \
           hd1d.o
#
#---------------------------------------------------
# Compilation rules

$(PROGRAM)   :  prebuild ${OBJECTS}
	@echo Linking object files ...                                          
	@$(COMPILER) $(FLAGS) ${OBJECTS} -o $@  
	@echo Cleaning up ...
prebuild :
	@echo "hd1d build started `date`"

%.o:%.f90
	@echo Compiling $^ ...
	@$(COMPILER) $(FLAGS) -c $< -o $@
%.o:%.F90
	@echo Compiling $^ ...
	@$(COMPILER) $(FLAGS) -c $< -o $@
#
clean :
	rm -f *.o *.mod 
