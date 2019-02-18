
SHELL=/bin/sh

COMPILE.f77     = gfortran
COMPILE.f90     = gfortran
LINK.f90	= gfortran
AR = ar -rs

SINGLEPREC = -frecord-marker=4
DOUBLEPREC = -frecord-marker=4 -fdefault-real-8 -fdefault-double-8
PRECISION  = ${DOUBLEPREC}

SEARCH = 
DEBUGFLAG = -C -g -fbacktrace -ffpe-trap=invalid,zero,overflow -fcheck=all
DEBUG     = 

OPT0 = -O0
OPT1 = -O1
OPT2 = -O2
OPT3 = -O3
OPT4 = -O4

CFLAG = ${SEARCH} -c -w ${DEBUG}

Cflag0  = ${CFLAG} ${PRECISION} ${OPT0}
Cflag1  = ${CFLAG} ${PRECISION} ${OPT1}
Cflag2  = ${CFLAG} ${PRECISION} ${OPT2}
Cflag3  = ${CFLAG} ${PRECISION} ${OPT3}
Cflag4  = ${CFLAG} ${PRECISION} ${OPT4}

CFLAGS = ${CFLAG} -fno-automatic

Lflag1  = ${PRECISION} ${MPILIB} ${DEBUG}
Lflag2  = ${PRECISION} ${DEBUG}

.SUFFIXES:
.SUFFIXES: .f90 .F90 .f .for .ftn .o

.f90.o:
	${COMPILE.f90} ${Cflag3} $<

.f.o:
	${COMPILE.f77} ${Cflag3} -ffixed-line-length-132 $<

MODULES = \
	ModKind.o \
	ModNumConst.o \
	ModConst.o \
	ModAMIE_Interface.o \
	ModErrors.o \
	ModEIEConductance.o \
	ModEIEFiles.o \
	ModEIE_Interface.o \
	ModExtras.o \
	ModTimeConvert.o \
	ModWeimer05.o \
	ModWeimer.o


OBJECTS = \
	AMIE_Library.o \
	EIE_Initialize.o \
	EIE_Library.o \
	EIE_UaLibrary.o \
	EIE_set_inputs.o \
	ihp.o \
	library.o \
	merge_str.o \
	readAMIEoutput.o \
	apex_routines.o \
	hmr89.o \
	iz94.o \
	mh86.o

OBJECTS_EXE =  ${MODULES}  ${OBJECTS} get_currents.o

EXE = get_currents.exe

${EXE}:	${OBJECTS_EXE}
	${LINK.f90} -o ${EXE} ${OBJECTS_EXE}

clean:	
	rm -f *~ core *.o *.mod fort.* a.out *.exe *.a *.so *.protex


