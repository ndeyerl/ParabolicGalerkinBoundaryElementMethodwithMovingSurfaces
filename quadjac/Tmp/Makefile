#implementation: 
#make //this builds the library libgaussq.a containing the objects
#make quadjac
#
# makefile for the pierreQP package to calculate Gauss quadrature points
#

# compile flags
CFLAGS = -g
#CFLAGS = -O2 

CC = gcc
AR = ar
RANLIB = ranlib

OBJS =  GaussCore.o Hermite.o Jacobi.o Laguerre.o

libgaussq.a: $(OBJS)
	$(AR) rv libgaussq.a $(OBJS)
	$(RANLIB) libgaussq.a

quadjac:  quadjac.o panelIA.o 
	gcc $(CFLAGS) -o quadjac quadjac.o panelIA.o libgaussq.a  -llapack -lm

quadjacNew:  quadjac.o panelIANew.o 
	gcc $(CFLAGS) -o quadjacNew quadjac.o panelIANew.o libgaussq.a  -llapack -lm

quadjactest:  quadjactest.o panelIA.o 
	gcc $(CFLAGS) -o quadjactest quadjactest.o panelIA.o libgaussq.a -lm

quadjactestNew:  quadjactest.o panelIANew.o 
	gcc $(CFLAGS) -o quadjactestNew quadjactest.o panelIANew.o libgaussq.a -lm




.c.o:
	gcc -c $(CFLAGS) $<

