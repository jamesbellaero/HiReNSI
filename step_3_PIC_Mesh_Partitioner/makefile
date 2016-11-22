# =================================================================================================================
#
# Developer    : Babak Poursartip 
#                The University of Texas at Austin, TX, USA
# Code         : PIC (Partitioner and Input Creator)
# makefile     : PETSc Library
# Last Update  : 21 June 2012
# arash: 27 October
# =================================================================================================================

info ::
	@echo "do:"
	@echo "  make SOURCE"
	@echo "    where SOURCE is one of ${SOURCESF}"

FFLAGS	    = -O2 -I${METIS_DIR}/Lib -I${HDF5_DIR}/include
FLINKER	    = ${PETSC_DIR}/arch-icc-opt/bin/mpif90 
CPPFLAGS    = 
FPPFLAGS    =
LOCDIR      = /home/jamesbell/Research/Meshes/Mesh

#ALL_SRC	    :=$(wildcard *.F90)
#SRC         :=$(filter-out %GlobalVariables.F90,$(ALL_SRC))
SRC         :=  Parameters.o \
                PIC_Convert.o \
                Inf.o \
                PIC_Input.o \
                PIC_Output.o \
                PIC_Output_Binary.o \
                PIC_Output_HDF5.o \
                PIC_ReNumbering.o \
                PIC_Heterogeneous_Material.o \
                PIC_Inversion_Data_Structure.o \
                PIC_HDF5_Mat_Vis.o \
                PML_Parameters.o \
                ShapeFunc.o \
                PIC_MassDampStiffPMLNNZero.o \
                Assemble.o \
                PIC_ObjectSize.o \
                PIC_Main.o

OBJS        :=$(patsubst %.F90,%.o,$(SRC))

SOURCESH     = ${METIS_DIR}/Lib/metis.h

LIBS	    = -lmetis -lhdf5_fortran -lhdf5hl_fortran
LIBDIR	    = -L${METIS_DIR} -L${HDF5_DIR}

MANSEC      = KSP
#CLEANFILES = rhs.vtk solution.vtk

main2 :  $(SOURCESF)  
	@echo Linking...
	-$(FLINKER) $(FFLAGS) $(SOURCESF) -o main *.o
	@echo Done!

main :  $(OBJS)
	@echo Linking...
	@echo $(OBJS)
	ifort $(LIBDIR) $(LIBS) -o $@ $^
%.o:	%.F90
	ifort $(FFLAGS) -c $<
	@echo Done!




info ::
	@echo "end of make file "

# clean rule for object files
cln ::
	@echo "Clean Object files..."
	rm -rf *.o main
	rm -rf *.mod main
	rm -rf fort.* main



