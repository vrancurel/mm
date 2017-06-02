SRCDIR = src
INCDIR = include
OBJDIR = objs
BINDIR = bin
LIBDIR = lib

AR = ar cq
RANLIB = ranlib

#OFLAGS = -O2
MFLAGS =

INCS = -I$(INCDIR)
WFLAGS = #-Wall -Wunused #-Werror 
CFLAGS = $(MFLAGS) $(INCS) $(OFLAGS) -g -D_REENTRANT -fPIC $(WFLAGS) -pg

CC  = gcc $(CFLAGS)
DO_PROTO = ../build-tools/do_proto

LIBMM=$(LIBDIR)/libmm.a

OBJS = \
$(OBJDIR)/memory.o \
$(OBJDIR)/err.o \
$(OBJDIR)/machine.o \
$(OBJDIR)/init.o \
$(OBJDIR)/pxop.o \
$(OBJDIR)/copy.o \
$(OBJDIR)/solve.o \
$(OBJDIR)/chfactor.o \
$(OBJDIR)/vecop.o \
$(OBJDIR)/matop.o \
$(OBJDIR)/submat.o \
$(OBJDIR)/lufactor.o \
$(OBJDIR)/mvmix.o \
$(OBJDIR)/mvutils.o \
$(OBJDIR)/mvmodnorm.o \
$(OBJDIR)/invnorm.o \
$(OBJDIR)/mvmoduni.o \
$(OBJDIR)/unirandp.o \
$(OBJDIR)/mvnsamp.o \
$(OBJDIR)/bvnsamp.o \
$(OBJDIR)/misc.o \
$(OBJDIR)/gfx.o

all: $(LIBMM)

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) -o $@ -c $<

$(LIBMM): $(OBJS)
	rm -f $(LIBMM)
	$(AR) $(LIBMM) $(OBJS)
	$(RANLIB) $(LIBMM)

clean:
	rm -f $(OBJS) $(LIBMM)

fproto:
	touch $(SRCDIR)/*.c
	$(DO_PROTO) $(INCS) -P$(SRCDIR) $(INCDIR)/*.h 
