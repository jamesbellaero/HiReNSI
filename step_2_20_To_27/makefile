# project : PML
# makefile parallel

F90    = ifort
FFLAGS = -O3

OBJS   = ShapeFunc.o \
		Copy.o \
		Main.o

inv   : $(OBJS)
	@echo Linking...
	$(F90) $(FFLAGS) $(OBJS) -o main
	@echo Done!

%.o   : %.F90
	@echo Compiling...
	$(F90) $(FFLAGS) -c $<



# clean rule for object files
cln  :
	@echo Clean Object files...
	rm -rf *.o main
	rm -rf *.mod main
	rm -rf fort.* main
