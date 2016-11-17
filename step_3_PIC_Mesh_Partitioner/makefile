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

FFLAGS	    = -O3
CPPFLAGS    = 
FPPFLAGS    =
LOCDIR      = home/jamesbell/Research/hirensi/step_3_PIC_Mesh_Partitioner

SOURCESF    =     Parameters.o \
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

SOURCEH     =

OBJSF       =  $(SOURCESF:%.F90=%.o)

MANSEC      = KSP
#CLEANFILES = rhs.vtk solution.vtk


include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

#main :  $(OBJSF) chkopts
#	@echo Linking...
#	-${FLINKER} $(FFLAGS) $(OBJSF) -o main ${PETSC_KSP_LIB}
#	@echo Done!

main :  $(SOURCESF) chkopts
	@echo Linking...
	@echo ${SOURCESF}
	@echo Relinking...
	@echo ${OBJSF}
	@echo Unlinking...
	-${FLINKER}  $(OBJSF) -o main *.o ${PETSC_KSP_LIB}
	@echo Done!

include ${PETSC_DIR}/conf/test


info ::
	@echo "end of make file "

# clean rule for object files
cln ::
	@echo "Clean Object files..."
	rm -rf *.o main
	rm -rf *.mod main




