# DO NOT CHANGE ANYTHING HERE!!! ##################################
# (unless you know *exactly* what you're doing...) 

# objects for automa
AOBJ = zygotic.o maternal.o integrate.o \
         ../lam/error.o solvers.o automa.o

#fly_sa objects
FOBJ = zygotic.o maternal.o integrate.o solvers.o \
       score.o translate.o ../lam/distributions.o ../lam/random.o ../lam/error.o
# serial code
FSOBJ = moves.o fly_sa.o ../lam/lsa.o savestate.o
# parallel code
FPOBJ = moves-mpi.o fly_sa-mpi.o ../lam/lsa-mpi.o savestate-mpi.o


#printscore objects
POBJ = zygotic.o maternal.o integrate.o \
       ../lam/error.o solvers.o score.o printscore.o

#unfold objects
UOBJ = zygotic.o maternal.o integrate.o \
	 ../lam/error.o solvers.o unfold.o score.o

#unstable-manifold objects
UMOBJ = zygotic.o maternal.o integrate.o \
	 ../lam/error.o solvers.o unstable-manifold.o score.o \

#calculate_volumes objects
COBJ = zygotic.o maternal.o integrate.o \
	 ../lam/error.o solvers.o calculate_volumes.o score.o \
	../linpack/dasum.o ../linpack/daxpy.o ../linpack/dcopy.o \
	../linpack/ddot.o ../linpack/dgebak.o ../linpack/dgebal.o \
	../linpack/dgeco.o ../linpack/dgeev.o ../linpack/dgefa.o \
	../linpack/dgehd2.o ../linpack/dgehrd.o ../linpack/dgemm.o \
	../linpack/dgemv.o ../linpack/dger.o ../linpack/dgesl.o \
	../linpack/dhseqr.o ../linpack/dlabad.o ../linpack/dlacpy.o \
	../linpack/dladiv.o ../linpack/dlaexc.o ../linpack/dlahqr.o \
	../linpack/dlahr2.o ../linpack/dlaln2.o ../linpack/dlamch.o \
	../linpack/dlange.o ../linpack/dlanv2.o ../linpack/dlapy2.o \
	../linpack/dlaqr0.o ../linpack/dlaqr1.o ../linpack/dlaqr2.o \
	../linpack/dlaqr3.o ../linpack/dlaqr4.o ../linpack/dlaqr5.o \
	../linpack/dlarf.o ../linpack/dlarfb.o ../linpack/dlarfg.o \
	../linpack/dlarft.o ../linpack/dlarfx.o ../linpack/dlartg.o \
	../linpack/dlascl.o ../linpack/dlaset.o ../linpack/dlassq.o \
	../linpack/dlasy2.o ../linpack/dnrm2.o ../linpack/dorg2r.o \
	../linpack/dorghr.o ../linpack/dorgqr.o ../linpack/drot.o \
	../linpack/dscal.o ../linpack/dswap.o ../linpack/dtrevc.o \
	../linpack/dtrexc.o ../linpack/dtrmm.o ../linpack/dtrmv.o \
	../linpack/idamax.o ../linpack/lsame.o ../linpack/xerbla.o \
	../linpack/ieeeck.o ../linpack/iparmq.o ../linpack/ilaenv.o \
	../linpack/dgedi.o 

#newtonraphson objects
NOBJ = zygotic.o maternal.o integrate.o \
	 ../lam/error.o solvers.o newtonraphson.o score.o \
	../linpack/dasum.o ../linpack/daxpy.o ../linpack/dcopy.o \
	../linpack/ddot.o ../linpack/dgebak.o ../linpack/dgebal.o \
	../linpack/dgeco.o ../linpack/dgeev.o ../linpack/dgefa.o \
	../linpack/dgehd2.o ../linpack/dgehrd.o ../linpack/dgemm.o \
	../linpack/dgemv.o ../linpack/dger.o ../linpack/dgesl.o \
	../linpack/dhseqr.o ../linpack/dlabad.o ../linpack/dlacpy.o \
	../linpack/dladiv.o ../linpack/dlaexc.o ../linpack/dlahqr.o \
	../linpack/dlahr2.o ../linpack/dlaln2.o ../linpack/dlamch.o \
	../linpack/dlange.o ../linpack/dlanv2.o ../linpack/dlapy2.o \
	../linpack/dlaqr0.o ../linpack/dlaqr1.o ../linpack/dlaqr2.o \
	../linpack/dlaqr3.o ../linpack/dlaqr4.o ../linpack/dlaqr5.o \
	../linpack/dlarf.o ../linpack/dlarfb.o ../linpack/dlarfg.o \
	../linpack/dlarft.o ../linpack/dlarfx.o ../linpack/dlartg.o \
	../linpack/dlascl.o ../linpack/dlaset.o ../linpack/dlassq.o \
	../linpack/dlasy2.o ../linpack/dnrm2.o ../linpack/dorg2r.o \
	../linpack/dorghr.o ../linpack/dorgqr.o ../linpack/drot.o \
	../linpack/dscal.o ../linpack/dswap.o ../linpack/dtrevc.o \
	../linpack/dtrexc.o ../linpack/dtrmm.o ../linpack/dtrmv.o \
	../linpack/idamax.o ../linpack/lsame.o ../linpack/xerbla.o \
	../linpack/ieeeck.o ../linpack/iparmq.o ../linpack/ilaenv.o

#scramble objects
SOBJ = zygotic.o maternal.o integrate.o \
	 ../lam/error.o score.o scramble.o translate.o solvers.o

SOURCES = `ls *.c`

#Below here are the rules for building things

all: $(FLYEXECS)

# special cases: dependencies and flags for individual .c files

automa.o: automa.c
	$(CC) -c $(CFLAGS) automa.c

fly_sa.o: fly_sa.c
	$(CC) -c $(CFLAGS) $(VFLAGS) fly_sa.c

printscore.o: printscore.c
	$(CC) -c $(CFLAGS) $(VFLAGS) printscore.c

scramble.o: scramble.c
	$(CC) -c $(CFLAGS) $(VFLAGS) scramble.c

unfold.o: unfold.c
	$(CC) -c $(CFLAGS) $(VFLAGS) unfold.c

unstable-manifold.o: unstable-manifold.c
	$(CC) -c $(CFLAGS) $(VFLAGS) unstable-manifold.c

calculate_volumes.o: calculate_volumes.c
	$(CC) -c $(CFLAGS) $(VFLAGS) calculate_volumes.c

newtonraphson.o: newtonraphson.c
	$(CC) -c $(CFLAGS) $(VFLAGS) newtonraphson.c

zygotic.o: zygotic.c
	$(CC) -c $(CFLAGS) $(KFLAGS) zygotic.c

# parallel stuff

fly_sa-mpi.o: fly_sa.c
	$(MPICC) -c -o fly_sa-mpi.o $(MPIFLAGS) $(CFLAGS) $(VFLAGS) fly_sa.c

lsa-mpi.o: ../lam/lsa.c
	$(MPICC) -c -o ../lam/lsa-mpi.o $(MPIFLAGS) $(CFLAGS) ../lam/lsa.c

moves-mpi.o: moves.c 
	$(MPICC) -c -o moves-mpi.o $(MPIFLAGS) $(CFLAGS) moves.c

savestate-mpi.o: savestate.c
	$(MPICC) -c -o savestate-mpi.o $(MPIFLAGS) $(CFLAGS) savestate.c

# executable targets: serial ...

automa: $(AOBJ)
	$(CC) -o automa $(CFLAGS) $(AOBJ) $(LIBS)

fly_sa: $(FOBJ) $(FSOBJ)
	$(CC) -o fly_sa $(CFLAGS) $(FOBJ) $(FSOBJ) $(FLIBS) 

printscore: $(POBJ)
	$(CC) -o printscore $(CFLAGS) $(POBJ) $(LIBS) 

unfold: $(UOBJ)
	$(CC) -o unfold $(CFLAGS) $(UOBJ) $(LIBS) 

unstable-manifold: $(UMOBJ)
	$(CC) -o unstable-manifold $(CFLAGS) $(UMOBJ) $(LIBS) $(FORTLIBS)

calculate_volumes: $(COBJ)
	$(CC) -o calculate_volumes $(CFLAGS) $(COBJ) $(LIBS) $(FORTLIBS)

newtonraphson: $(NOBJ)
	$(CC) -o newtonraphson $(CFLAGS) $(NOBJ) $(LIBS) $(FORTLIBS)

scramble: $(SOBJ)
	$(CC) -o scramble $(CFLAGS) $(SOBJ) $(LIBS) 

# ... and parallel

fly_sa.mpi: $(FOBJ) $(FPOBJ)
	$(MPICC) -o fly_sa.mpi $(CFLAGS) $(FOBJ) $(FPOBJ) $(FLIBS)

# special target to make the dynamical systems analysis targets

dynamicalsystems: newtonraphson unstable-manifold calculate_volumes

# ... and here are the cleanup and make deps rules

clean:
	rm -f *.o core*

Makefile: ${FRC}
	rm -f $@
	cp basic.mk $@
	echo "#Automatically generated dependencies list#" >> $@
	${CC} $(INCLUDES) -M ${SOURCES} >> $@
	chmod -w $@

