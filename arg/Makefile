# ***************************************************************************
#
#       Makefile template for programs. - JBA - Rev Nov 1993
#
#       This template works fine on the Sun.  There are two changes that
#       must be made for it to work on the DecStation:
#               - The conditional stuff for GRAPH & DEBUG must be commented out
#               - The libX.a must be changed to libX11.a
#
# ***************************************************************************

MACH = linux

OPENWINHOME = /usr/openwin

SHELL := /bin/sh

PROGRAM	      = demo

ARGC	      = argtool

CC            = gcc

#LABELS        = -DPE_LABEL=eos$(PROGRAM) -DPE_NAME=\"$(PROGRAM)\"

LABELS       =

CFLAGS        = -I${HOME}/incl -g $(LABELS)

FC            = g77

F77           = g77

FFLAGS        = -I. -g

DEST	      = /home/joea/bin/linux

MAKEFILE      = Makefile

EXTHDRS	      =

#EXTLIBS       =  ${EOS}/${MACH}/include/pe_main.c

HDRS	      = div_server.h

LDFLAGS	      =

SLIBS	      = -lm

ULIBS	      =

SUFFIX        = .arg:sC .arf:sF

LIBS	      = ${HOME}/lib/${MACH}/libjba.a \
		${HOME}/lib/${MACH}/libargtool.a
#		${EOS}/${MACH}/lib/libeos.a \

GLIBS        = 	${HOME}/lib/${MACH}/libgpac.a \
		$(OPENWINHOME)/lib/libxview.a \
		$(OPENWINHOME)/lib/libolgx.a \
		/usr//X11R6/lib/libX11.so

EXTLIBS       =

#EXTLIBS       =  ${EOS}/${MACH}/include/pe_main.c

LINKER	      = g++
#LINKER	      = g77

MAKEFILE      = Makefile

OBJS	      = demo.o \
		dserv.o

PRINT	      = lpr

SRCS	      = demo.arg \
		dserv.c

.SUFFIXES: .arg .arf

graph := EXTLIBS = $(GLIBS)

debug := FFLAGS += -g 

.arf.f:
	argtool $*.arf -ftn

.arg.c:
	$(ARGC) $*.arg -c
	
.arg.o:
	$(ARGC) $*.arg -c
	$(CC) $(CFLAGS) -c $*.c
#	rm -f $*.c

.arf.o:
	argtool $*.arf -ftn
	$(FC) $(FFLAGS) -c $*.f
#	rm -f $*.f

.y.o:	
	yacc -d $*.y
	mv y.tab.h $*_y.tab.h
	awk '/SED_TRANSLATE/ {print $$3}' $*.y > $*_y.awk
	echo 's/c/c/' >> $*_y.awk
	sed -f $*_y.awk y.tab.c > $*.c
	cc $(CFLAGS) -c $*.c
	rm $*.c
	rm y.tab.c
	rm $*_y.awk

.l.c:
	lex $*.l
	awk '/SED_TRANSLATE/ {print $$3}' $*.l > $*_l.awk
	echo 's/c/c/' >> $*_l.awk
	sed -f $*_l.awk lex.yy.c > $*.c
	rm lex.yy.c
	rm $*_l.awk

all:		$(PROGRAM)

$(PROGRAM):     $(OBJS) $(LIBS) $(ULIBS) $(GLIBS)
		@echo -n "Loading $(PROGRAM) ... "
		@if [ ! -f ../.gpac ]; \
		then \
			$(LINKER) $(LDFLAGS) $(OBJS) $(ULIBS) $(LIBS) $(SLIBS) $(EXTLIBS) -o $(PROGRAM); \
		else \
			$(LINKER) $(LDFLAGS) $(OBJS) $(ULIBS) $(LIBS) $(SLIBS) $(EXTLIBS) $(GLIBS) -o $(PROGRAM); \
		fi;
		@echo "done"

clean:;		@rm -f $(OBJS) $(PROGRAM)

depend:;	@mkmf -f $(MAKEFILE) PROGRAM=$(PROGRAM) DEST=$(DEST)

index:;		@ctags -wx $(HDRS) $(SRCS)

install:	$(PROGRAM)
		@echo Installing $(PROGRAM) in $(DEST)
		@install -s $(PROGRAM) $(DEST)
		@rm $(PROGRAM) *.o

print:;		@$(PRINT) $(HDRS) $(SRCS)

lprint:;	fold -w 132 $(HDRS) $(SRCS) | rsh adiron /usr/local/bin/imprint -J -Pip -L -h$(PROGRAM)

program graph debug:        $(PROGRAM)

tags:           $(HDRS) $(SRCS); @ctags $(HDRS) $(SRCS)

update:		$(DEST)/$(PROGRAM)

$(DEST)/$(PROGRAM): $(SRCS) $(LIBS) $(HDRS) $(EXTHDRS)
		@make -f $(MAKEFILE) DEST=$(DEST) install
###
dserv.o: div_server.h
