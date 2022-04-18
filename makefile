#GFORTRAN OS X ###
FC = gfortran
FFLAGS = -m64 -fdefault-real-8 -fcheck=all -Wall \
			# -g -ffpe-trap=zero,invalid,overflow,underflow
TARGET_ARCH =
LDFLAGS = -m64
LIBS = 

EXE   = EulerSolver
VPATH = mod

.SUFFIXES: .f90 .o

GRAFIC_HEADER =

SRCMOD =    			\
	constants.f90	\
	variables.f90	\
	IO.f90		

SRCMAIN =      			\
	alloc.f90			\
	cellProp.f90		\
	init.f90			\
	calBCs.f90			\
	calDissipation.f90	\
	calFlux.f90			\
	calResidual.f90		\
	RK4.f90				\
	EulerSolver.f90

OBJMOD = ${SRCMOD:.f90=.o}

OBJMAIN = ${SRCMAIN:.f90=.o}

OBJ = $(OBJMOD) $(OBJMAIN)

$(EXE): $(OBJ)
	$(FC) $(LDFLAGS) $(OBJ) $(LIBS) -o $(EXE)

%.o  : %.f90
	$(FC) $(FFLAGS) -c $<

# Define dependencies for modules
$(OBJMAIN): $(OBJMOD)
IO.o : constants.o variables.o
variables.o : constants.o


# clean any generated files
clean:
	-rm -f *.o *~ core *.mod

solve:
	make
	./EulerSolver