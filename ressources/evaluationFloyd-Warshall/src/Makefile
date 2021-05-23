#=======================================================================
# Basic C: makefile example to use as a reminder or as a template
#-----------------------------------------------------------------------
# Anas Francis
#=======================================================================

#only ok for project with no templated classes


#compiler
COMPILER = mpicc
#linker
LINKER = mpicc

#options for linker and compiler
FLAGS =-Wall -std=gnu99 -lm -fopenmp

EXE_NAME= output.mp

#can have several ones separated by spaces, only cpp files
SOURCES = francis.c

#PATH to extra header used in the project (when using libs not installed in the OS)
INC_PATHS= -I

#extra libs used in the project (for instance -lpthread)
LIBS	= -L

LINKOBJ	= $(SOURCES:.c=.o)

$(EXE_NAME):
	$(LINKER) -g -o $(EXE_NAME) $(SOURCES) $(FLAGS)



# cleanup
clean:
	rm $(EXE_NAME)
