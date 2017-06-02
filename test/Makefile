
RM = rm -f

#OFLAGS = -O2
MFLAGS =

INCS = -I../libmm/include -I..
LIBS = -L../libmm/lib -lmm -lm
DEPENDS = ../libmm/lib/libmm.a

CFLAGS = $(MFLAGS) -g $(OFLAGS) -Wall -Werror $(INCS) -D_GNU_SOURCE -pg
LDFLAGS = $(MFLAGS) $(LIBS) -pg

PROGS = \
prob_test \

PROB_TEST_OBJS = prob_test.o

all: $(PROGS)

prob_test: $(PROB_TEST_OBJS) $(DEPENDS)
	$(CC) -o prob_test $(PROB_TEST_OBJS) $(LDFLAGS)

install:

clean:
	$(RM) *.o $(PROGS)
