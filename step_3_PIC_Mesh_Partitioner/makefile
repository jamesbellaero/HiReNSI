#================================================================================================================
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

FFLAGS	    = -O2 -fPIC -I${METIS_DIR}include -I${HDF5_DIR}/include -I${MPI_DIR}/include -g -check all -fpe0 -traceback -debug extended
FLINKER     = -ifort#/software/intel/impi/5.1.2.150/intel64/bin/mpiifort 
CPPFLAGS    = 
FPPFLAGS    =
#FINCLUDE    = -I${METIS_DIR}/Lib -I${METIS_DIR} -I${METIS_DIR}/Programs -I${HDF5_DIR}/include -I${MPI_DIR}/include
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
                #/home/jamesbell/Research/metis/metis-4.0.3/Lib/meshpart.o \
                PIC_Main.o

OBJS        :=$(patsubst %.F90,%.o,$(SRC))

SOURCESH     = ${METIS_DIR}/Lib/metis.h

LIBS	    =  -lmpifort -lfmpich -lmetis -lmpi -lmpigf -lhdf5_fortran -lhdf5hl_fortran
#LIBDIR	    = -L${METIS_DIR}/Lib -L${HDF5_DIR}/lib -L${MPI_DIR}/lib -L/software/intel/impi/5.1.2.150/intel64/lib/
LIBDIR      = -L${METIS_DIR}/build/Linux-x86_64/libmetis -L${HDF5_DIR}/lib -L${MPI_DIR}/lib -L/software/intel/impi/5.1.2.150/intel64/lib/


MANSEC      = KSP
#CLEANFILES = rhs.vtk solution.vtk

main :  $(OBJS)
	@echo Linking...
	$(FLINKER) $(FFLAGS) $(LIBDIR) $(LIBS) -o $@ $^
%.o:	%.F90
	@echo Compiling...
	$(FLINKER) $(FFLAGS) -c $<
	@echo Done!




info ::
	@echo "end of make file "

# clean rule for object files
cln ::
	@echo "Clean Object files..."
	rm -rf *.o main
	rm -rf *.mod main
	rm -rf fort.* main



