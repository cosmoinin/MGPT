# ----- Makefile FOR mgpt CODE -----
#

OPT1    =  -DNOGNU

OPTIONS =  $(OPTIMIZE) \
	   $(OPT1)


EXEC = mgpt

all: $(EXEC)

CC = gcc

CFLAGS = -g -O3 -fopenmp $(OPTIONS)

H_PATH = general_libs

OBJS	= main.o mglpt_io.o mglpt.o startrun.o \
	models.o mglpt_fns.o mglpt_diffeqs.o mglpt_quads.o \
	mglpt_postprocess.o \
	$(H_PATH)/diffeqs.o $(H_PATH)/mathutil.o \
	$(H_PATH)/clib.o $(H_PATH)/getparam.o \
	$(H_PATH)/mathfns.o $(H_PATH)/numrec.o \
	$(H_PATH)/inout.o $(H_PATH)/quads.o

INCL	= globaldefs.h data_struc_defs.h cmdline_defs.h protodefs.h \
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

check: all

install:
	install -d $(DESTDIR)$(bindir)
	install -m 0755 $(EXEC) $(DESTDIR)$(bindir)

uninstall:
	-rm $(DESTDIR)$(bindir)/$(EXEC)


.PHONY : all clean check install uninstall

