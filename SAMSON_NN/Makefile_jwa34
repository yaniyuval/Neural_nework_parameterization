# Makefile for various platforms
# Execute using Build csh-script only!
# Used together with Perl scripts in SRC/SCRIPT
# (C) 2005 Marat Khairoutdinov
#------------------------------------------------------------------
# uncomment to disable timers:
#
#NOTIMERS=-DDISABLE_TIMERS
#-----------------------------------------------------------------

SAM = SAM_$(ADV_DIR)_$(SGS_DIR)_$(RAD_DIR)_$(MICRO_DIR)

# Determine platform
PLATFORM := $(shell uname -s)

#----------------------------------------------------------------------
# Linux, Intel Compiler


ifeq ($(PLATFORM),Linux)

LIB_MPI = 
INC_MPI = 
INC_NETCDF = 
LIB_NETCDF = 


FF77 = mpif90 -c -fixed -extend_source
FF90 = mpif90 -c
CC = mpicc -c -DLINUX


FFLAGS = -O3
#FFLAGS = -g -ftrapuv -check all

# FFLAGS += -I${INC_MPI} -I${INC_NETCDF}
LD = mpif90
LDFLAGS = -lnetcdf -lnetcdff


endif

#----------------------------------------------
# you dont need to edit below this line


#compute the search path
dirs := . $(shell cat Filepath)
VPATH    := $(foreach dir,$(dirs),$(wildcard $(dir)))

.SUFFIXES:
.SUFFIXES: .f .f90 .c .o



all: $(SAM_DIR)/$(SAM)


SOURCES   := $(shell cat Srcfiles)

Depends: Srcfiles Filepath
				$(SAM_SRC)/SCRIPT/mkDepends Filepath Srcfiles > $@

Srcfiles: Filepath
				$(SAM_SRC)/SCRIPT/mkSrcfiles > $@

OBJS      := $(addsuffix .o, $(basename $(SOURCES)))

$(SAM_DIR)/$(SAM): $(OBJS)
				$(LD) -o $@ $(OBJS) $(LDFLAGS)


.f90.o:
				${FF90}  ${FFLAGS} $<
.f.o:
				${FF77}  ${FFLAGS} $<
.c.o:
				${CC}  ${CFLAGS} -I$(SAM_SRC)/TIMING $(NOTIMERS) $<



include Depends



clean:
				rm ./OBJ/*



