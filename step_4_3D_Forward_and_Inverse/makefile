# project
# makefile PARLALLEL FOR PETSC
# IMPORTANT NOTICE: This make file works only for ".F90".
# If there exists a file with other extensions,
# the "OBJSF" MUST be modified.

info ::
	@echo "do:"
	@echo "  make SOURCE"
	@echo "    where SOURCE is one of ${SOURCESF}"

FFLAGS	         = -O2
CPPFLAGS         = 
FPPFLAGS         =
LOCDIR           =

SOURCESF         =         Parameters.o \
			PTC_Input.o \
			PTC_Input_HDF5.o \
			PTC_Results.o \
			ParaView.o \
			IO.o \
			External_Links.o \
			External_Solver.o \
			Element_Matrices.o \
			Element_Matrices_3D.o \
			Assemble_PETSc.o \
			Assemble_PETSc_3D.o \
			PML.o \
			Source1.o \
			Source2.o \
			Eigenvalues.o \
			Geometry.o \
			Shape_Functions.o \
			Verification.o \
			Visualization.o \
			Obsolete.o \
			Solver_PETSc_Implicit.o \
			Solver_PETSc_Explicit.o \
			Solver_PETSc_Inverse.o \
			Program_Main.o \

SOURCEH          =

OBJSF            =  $(SOURCESF:%.F90=%.o)

MANSEC           = KSP
#CLEANFILES       = rhs.vtk solution.vtk


include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

#main :  $(OBJSF) chkopts
#	@echo Linking...
#	-${FLINKER} $(FFLAGS) $(OBJSF) -o main ${PETSC_KSP_LIB}
#	@echo Done!


main :  $(SOURCESF) chkopts
	@echo Linking...
	-${FLINKER} $(OBJSF) -o main *.o ${PETSC_KSP_LIB}
	@echo Done!

include ${PETSC_DIR}/conf/test

# ----------------------------------------------------------------
# clean rule for object files
#clean:
#	@echo Clean Object files...
#	rm -rf *.o main












