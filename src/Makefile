# ----- Makefile FOR mgpt CODE -----
#

OPTIMIZE =
OPT1    =

OPTIONS =  $(OPTIMIZE) \
	   $(OPT1)


EXEC = mgpt

all: $(EXEC)

# C compiler:
CC = gcc
# With icc (Intel):
#CC = icc

# Default CFLAGS:
CFLAGS = -g -O3 $(OPTIONS)
#
# With OpenMP:
#CFLAGS = -g -O3 -fopenmp $(OPTIONS)
# For icc use instead:
#CFLAGS = -g -O3 -qopenmp $(OPTIONS)


#
# Nothing to do below
#
H_PATH = general_libs

OBJS	= main.o mgpt_io.o mgpt.o startrun.o \
	models.o mgpt_diffeqs.o mgpt_quads.o \
	mgpt_postprocess.o \
	$(H_PATH)/diffeqs.o $(H_PATH)/mathutil.o \
	$(H_PATH)/clib.o $(H_PATH)/getparam.o \
	$(H_PATH)/mathfns.o $(H_PATH)/numrec.o \
	$(H_PATH)/inout.o $(H_PATH)/quads.o

INCL	= globaldefs.h cmdline_defs.h protodefs.h \
	models.h $(H_PATH)/diffeqs.h \
	$(H_PATH)/getparam.h $(H_PATH)/machines.h \
	$(H_PATH)/mathfns.h $(H_PATH)/stdinc.h \
	$(H_PATH)/strings.h $(H_PATH)/precision.h \
	$(H_PATH)/vectdefs.h $(H_PATH)/vectmath.h \
	$(H_PATH)/numrec.h \
	$(H_PATH)/inout.h $(H_PATH)/mathutil.h \
	$(H_PATH)/switchs.h $(H_PATH)/quads.h

$(EXEC): $(OBJS) 
	($(CC) $(OBJS) $(LIBS) $(CFLAGS) -o $@ -lm; cp $(EXEC) ../)

$(OBJS): $(INCL)

clean:
	(rm -f $(OBJS) $(EXEC); rm -fR mgpt.dSYM; rm ../$(EXEC))


.PHONY : all clean check install uninstall

