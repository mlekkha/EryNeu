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

#Automatically generated dependencies list#
automa.o: automa.c /usr/include/stdc-predef.h /usr/include/stdlib.h \
 /usr/include/features.h /usr/include/sys/cdefs.h \
 /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h \
 /usr/include/gnu/stubs-64.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/stddef.h \
 /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h \
 /usr/include/endian.h /usr/include/bits/endian.h \
 /usr/include/bits/byteswap.h /usr/include/bits/types.h \
 /usr/include/bits/typesizes.h /usr/include/bits/byteswap-16.h \
 /usr/include/sys/types.h /usr/include/time.h /usr/include/sys/select.h \
 /usr/include/bits/select.h /usr/include/bits/sigset.h \
 /usr/include/bits/time.h /usr/include/sys/sysmacros.h \
 /usr/include/bits/pthreadtypes.h /usr/include/alloca.h \
 /usr/include/bits/stdlib-float.h /usr/include/stdio.h \
 /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/stdarg.h \
 /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h \
 /usr/include/string.h /usr/include/xlocale.h /usr/include/unistd.h \
 /usr/include/bits/posix_opt.h /usr/include/bits/environments.h \
 /usr/include/bits/confname.h /usr/include/getopt.h /usr/include/ctype.h \
 ../lam/error.h ../lam/global.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/float.h integrate.h \
 maternal.h /usr/include/gsl/gsl_errno.h /usr/include/errno.h \
 /usr/include/bits/errno.h /usr/include/linux/errno.h \
 /usr/include/asm/errno.h /usr/include/asm-generic/errno.h \
 /usr/include/asm-generic/errno-base.h /usr/include/gsl/gsl_types.h \
 /usr/include/gsl/gsl_spline.h /usr/include/gsl/gsl_interp.h \
 /usr/include/gsl/gsl_inline.h maternal.h solvers.h zygotic.h \
 /usr/include/sys/resource.h /usr/include/bits/resource.h
calculate_volumes.o: calculate_volumes.c /usr/include/stdc-predef.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/float.h \
 /usr/include/stdlib.h /usr/include/features.h /usr/include/sys/cdefs.h \
 /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h \
 /usr/include/gnu/stubs-64.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/stddef.h \
 /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h \
 /usr/include/endian.h /usr/include/bits/endian.h \
 /usr/include/bits/byteswap.h /usr/include/bits/types.h \
 /usr/include/bits/typesizes.h /usr/include/bits/byteswap-16.h \
 /usr/include/sys/types.h /usr/include/time.h /usr/include/sys/select.h \
 /usr/include/bits/select.h /usr/include/bits/sigset.h \
 /usr/include/bits/time.h /usr/include/sys/sysmacros.h \
 /usr/include/bits/pthreadtypes.h /usr/include/alloca.h \
 /usr/include/bits/stdlib-float.h /usr/include/stdio.h \
 /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/stdarg.h \
 /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h \
 /usr/include/string.h /usr/include/xlocale.h /usr/include/unistd.h \
 /usr/include/bits/posix_opt.h /usr/include/bits/environments.h \
 /usr/include/bits/confname.h /usr/include/getopt.h /usr/include/math.h \
 /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h \
 /usr/include/bits/huge_vall.h /usr/include/bits/inf.h \
 /usr/include/bits/nan.h /usr/include/bits/mathdef.h \
 /usr/include/bits/mathcalls.h /usr/include/sys/resource.h \
 /usr/include/bits/resource.h ../lam/error.h ../lam/global.h integrate.h \
 maternal.h /usr/include/gsl/gsl_errno.h /usr/include/errno.h \
 /usr/include/bits/errno.h /usr/include/linux/errno.h \
 /usr/include/asm/errno.h /usr/include/asm-generic/errno.h \
 /usr/include/asm-generic/errno-base.h /usr/include/gsl/gsl_types.h \
 /usr/include/gsl/gsl_spline.h /usr/include/gsl/gsl_interp.h \
 /usr/include/gsl/gsl_inline.h solvers.h zygotic.h score.h
fly_sa.o: fly_sa.c /usr/include/stdc-predef.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/float.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/limits.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/syslimits.h \
 /usr/include/limits.h /usr/include/features.h /usr/include/sys/cdefs.h \
 /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h \
 /usr/include/gnu/stubs-64.h /usr/include/bits/posix1_lim.h \
 /usr/include/bits/local_lim.h /usr/include/linux/limits.h \
 /usr/include/bits/posix2_lim.h /usr/include/stdlib.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/stddef.h \
 /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h \
 /usr/include/endian.h /usr/include/bits/endian.h \
 /usr/include/bits/byteswap.h /usr/include/bits/types.h \
 /usr/include/bits/typesizes.h /usr/include/bits/byteswap-16.h \
 /usr/include/sys/types.h /usr/include/time.h /usr/include/sys/select.h \
 /usr/include/bits/select.h /usr/include/bits/sigset.h \
 /usr/include/bits/time.h /usr/include/sys/sysmacros.h \
 /usr/include/bits/pthreadtypes.h /usr/include/alloca.h \
 /usr/include/bits/stdlib-float.h /usr/include/stdio.h \
 /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/stdarg.h \
 /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h \
 /usr/include/string.h /usr/include/xlocale.h /usr/include/unistd.h \
 /usr/include/bits/posix_opt.h /usr/include/bits/environments.h \
 /usr/include/bits/confname.h /usr/include/getopt.h ../lam/error.h \
 ../lam/global.h ../lam/distributions.h integrate.h maternal.h \
 /usr/include/gsl/gsl_errno.h /usr/include/errno.h \
 /usr/include/bits/errno.h /usr/include/linux/errno.h \
 /usr/include/asm/errno.h /usr/include/asm-generic/errno.h \
 /usr/include/asm-generic/errno-base.h /usr/include/gsl/gsl_types.h \
 /usr/include/gsl/gsl_spline.h /usr/include/gsl/gsl_interp.h \
 /usr/include/gsl/gsl_inline.h moves.h ../lam/sa.h score.h \
 ../lam/random.h solvers.h zygotic.h /usr/include/sys/resource.h \
 /usr/include/bits/resource.h
integrate.o: integrate.c /usr/include/stdc-predef.h /usr/include/stdio.h \
 /usr/include/features.h /usr/include/sys/cdefs.h \
 /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h \
 /usr/include/gnu/stubs-64.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/stddef.h \
 /usr/include/bits/types.h /usr/include/bits/typesizes.h \
 /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/stdarg.h \
 /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h \
 /usr/include/stdlib.h /usr/include/bits/waitflags.h \
 /usr/include/bits/waitstatus.h /usr/include/endian.h \
 /usr/include/bits/endian.h /usr/include/bits/byteswap.h \
 /usr/include/bits/byteswap-16.h /usr/include/sys/types.h \
 /usr/include/time.h /usr/include/sys/select.h /usr/include/bits/select.h \
 /usr/include/bits/sigset.h /usr/include/bits/time.h \
 /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h \
 /usr/include/alloca.h /usr/include/bits/stdlib-float.h \
 /usr/include/math.h /usr/include/bits/huge_val.h \
 /usr/include/bits/huge_valf.h /usr/include/bits/huge_vall.h \
 /usr/include/bits/inf.h /usr/include/bits/nan.h \
 /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/float.h \
 /usr/include/string.h /usr/include/xlocale.h \
 /usr/include/gsl/gsl_errno.h /usr/include/errno.h \
 /usr/include/bits/errno.h /usr/include/linux/errno.h \
 /usr/include/asm/errno.h /usr/include/asm-generic/errno.h \
 /usr/include/asm-generic/errno-base.h /usr/include/gsl/gsl_types.h \
 /usr/include/gsl/gsl_spline.h /usr/include/gsl/gsl_interp.h \
 /usr/include/gsl/gsl_inline.h ../lam/error.h ../lam/global.h maternal.h \
 integrate.h solvers.h zygotic.h /usr/include/sys/resource.h \
 /usr/include/bits/resource.h score.h
maternal.o: maternal.c /usr/include/stdc-predef.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/float.h \
 /usr/include/stdlib.h /usr/include/features.h /usr/include/sys/cdefs.h \
 /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h \
 /usr/include/gnu/stubs-64.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/stddef.h \
 /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h \
 /usr/include/endian.h /usr/include/bits/endian.h \
 /usr/include/bits/byteswap.h /usr/include/bits/types.h \
 /usr/include/bits/typesizes.h /usr/include/bits/byteswap-16.h \
 /usr/include/sys/types.h /usr/include/time.h /usr/include/sys/select.h \
 /usr/include/bits/select.h /usr/include/bits/sigset.h \
 /usr/include/bits/time.h /usr/include/sys/sysmacros.h \
 /usr/include/bits/pthreadtypes.h /usr/include/alloca.h \
 /usr/include/bits/stdlib-float.h /usr/include/math.h \
 /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h \
 /usr/include/bits/huge_vall.h /usr/include/bits/inf.h \
 /usr/include/bits/nan.h /usr/include/bits/mathdef.h \
 /usr/include/bits/mathcalls.h /usr/include/ctype.h \
 /usr/include/xlocale.h /usr/include/string.h ../lam/error.h \
 ../lam/global.h maternal.h /usr/include/stdio.h /usr/include/libio.h \
 /usr/include/_G_config.h /usr/include/wchar.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/stdarg.h \
 /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h \
 integrate.h /usr/include/gsl/gsl_errno.h /usr/include/errno.h \
 /usr/include/bits/errno.h /usr/include/linux/errno.h \
 /usr/include/asm/errno.h /usr/include/asm-generic/errno.h \
 /usr/include/asm-generic/errno-base.h /usr/include/gsl/gsl_types.h \
 /usr/include/gsl/gsl_spline.h /usr/include/gsl/gsl_interp.h \
 /usr/include/gsl/gsl_inline.h zygotic.h /usr/include/sys/resource.h \
 /usr/include/bits/resource.h
moves.o: moves.c /usr/include/stdc-predef.h /usr/include/math.h \
 /usr/include/features.h /usr/include/sys/cdefs.h \
 /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h \
 /usr/include/gnu/stubs-64.h /usr/include/bits/huge_val.h \
 /usr/include/bits/huge_valf.h /usr/include/bits/huge_vall.h \
 /usr/include/bits/inf.h /usr/include/bits/nan.h \
 /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h \
 /usr/include/stdlib.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/stddef.h \
 /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h \
 /usr/include/endian.h /usr/include/bits/endian.h \
 /usr/include/bits/byteswap.h /usr/include/bits/types.h \
 /usr/include/bits/typesizes.h /usr/include/bits/byteswap-16.h \
 /usr/include/sys/types.h /usr/include/time.h /usr/include/sys/select.h \
 /usr/include/bits/select.h /usr/include/bits/sigset.h \
 /usr/include/bits/time.h /usr/include/sys/sysmacros.h \
 /usr/include/bits/pthreadtypes.h /usr/include/alloca.h \
 /usr/include/bits/stdlib-float.h /usr/include/stdio.h \
 /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/stdarg.h \
 /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h \
 /usr/include/string.h /usr/include/xlocale.h ../lam/error.h \
 ../lam/global.h /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/float.h \
 moves.h ../lam/sa.h score.h maternal.h ../lam/random.h \
 ../lam/distributions.h
newtonraphson.o: newtonraphson.c /usr/include/stdc-predef.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/float.h \
 /usr/include/stdlib.h /usr/include/features.h /usr/include/sys/cdefs.h \
 /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h \
 /usr/include/gnu/stubs-64.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/stddef.h \
 /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h \
 /usr/include/endian.h /usr/include/bits/endian.h \
 /usr/include/bits/byteswap.h /usr/include/bits/types.h \
 /usr/include/bits/typesizes.h /usr/include/bits/byteswap-16.h \
 /usr/include/sys/types.h /usr/include/time.h /usr/include/sys/select.h \
 /usr/include/bits/select.h /usr/include/bits/sigset.h \
 /usr/include/bits/time.h /usr/include/sys/sysmacros.h \
 /usr/include/bits/pthreadtypes.h /usr/include/alloca.h \
 /usr/include/bits/stdlib-float.h /usr/include/stdio.h \
 /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/stdarg.h \
 /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h \
 /usr/include/string.h /usr/include/xlocale.h /usr/include/unistd.h \
 /usr/include/bits/posix_opt.h /usr/include/bits/environments.h \
 /usr/include/bits/confname.h /usr/include/getopt.h /usr/include/math.h \
 /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h \
 /usr/include/bits/huge_vall.h /usr/include/bits/inf.h \
 /usr/include/bits/nan.h /usr/include/bits/mathdef.h \
 /usr/include/bits/mathcalls.h /usr/include/sys/resource.h \
 /usr/include/bits/resource.h ../lam/error.h ../lam/global.h integrate.h \
 maternal.h /usr/include/gsl/gsl_errno.h /usr/include/errno.h \
 /usr/include/bits/errno.h /usr/include/linux/errno.h \
 /usr/include/asm/errno.h /usr/include/asm-generic/errno.h \
 /usr/include/asm-generic/errno-base.h /usr/include/gsl/gsl_types.h \
 /usr/include/gsl/gsl_spline.h /usr/include/gsl/gsl_interp.h \
 /usr/include/gsl/gsl_inline.h solvers.h zygotic.h score.h
printscore.o: printscore.c /usr/include/stdc-predef.h /usr/include/math.h \
 /usr/include/features.h /usr/include/sys/cdefs.h \
 /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h \
 /usr/include/gnu/stubs-64.h /usr/include/bits/huge_val.h \
 /usr/include/bits/huge_valf.h /usr/include/bits/huge_vall.h \
 /usr/include/bits/inf.h /usr/include/bits/nan.h \
 /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h \
 /usr/include/stdio.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/stddef.h \
 /usr/include/bits/types.h /usr/include/bits/typesizes.h \
 /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/stdarg.h \
 /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h \
 /usr/include/stdlib.h /usr/include/bits/waitflags.h \
 /usr/include/bits/waitstatus.h /usr/include/endian.h \
 /usr/include/bits/endian.h /usr/include/bits/byteswap.h \
 /usr/include/bits/byteswap-16.h /usr/include/sys/types.h \
 /usr/include/time.h /usr/include/sys/select.h /usr/include/bits/select.h \
 /usr/include/bits/sigset.h /usr/include/bits/time.h \
 /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h \
 /usr/include/alloca.h /usr/include/bits/stdlib-float.h \
 /usr/include/string.h /usr/include/xlocale.h /usr/include/unistd.h \
 /usr/include/bits/posix_opt.h /usr/include/bits/environments.h \
 /usr/include/bits/confname.h /usr/include/getopt.h ../lam/error.h \
 ../lam/global.h /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/float.h \
 integrate.h maternal.h /usr/include/gsl/gsl_errno.h /usr/include/errno.h \
 /usr/include/bits/errno.h /usr/include/linux/errno.h \
 /usr/include/asm/errno.h /usr/include/asm-generic/errno.h \
 /usr/include/asm-generic/errno-base.h /usr/include/gsl/gsl_types.h \
 /usr/include/gsl/gsl_spline.h /usr/include/gsl/gsl_interp.h \
 /usr/include/gsl/gsl_inline.h score.h solvers.h zygotic.h \
 /usr/include/sys/resource.h /usr/include/bits/resource.h
savestate.o: savestate.c /usr/include/stdc-predef.h /usr/include/math.h \
 /usr/include/features.h /usr/include/sys/cdefs.h \
 /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h \
 /usr/include/gnu/stubs-64.h /usr/include/bits/huge_val.h \
 /usr/include/bits/huge_valf.h /usr/include/bits/huge_vall.h \
 /usr/include/bits/inf.h /usr/include/bits/nan.h \
 /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h \
 /usr/include/stdio.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/stddef.h \
 /usr/include/bits/types.h /usr/include/bits/typesizes.h \
 /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/stdarg.h \
 /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h \
 /usr/include/string.h /usr/include/xlocale.h /usr/include/stdlib.h \
 /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h \
 /usr/include/endian.h /usr/include/bits/endian.h \
 /usr/include/bits/byteswap.h /usr/include/bits/byteswap-16.h \
 /usr/include/sys/types.h /usr/include/time.h /usr/include/sys/select.h \
 /usr/include/bits/select.h /usr/include/bits/sigset.h \
 /usr/include/bits/time.h /usr/include/sys/sysmacros.h \
 /usr/include/bits/pthreadtypes.h /usr/include/alloca.h \
 /usr/include/bits/stdlib-float.h /usr/include/unistd.h \
 /usr/include/bits/posix_opt.h /usr/include/bits/environments.h \
 /usr/include/bits/confname.h /usr/include/getopt.h ../lam/error.h \
 ../lam/global.h /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/float.h \
 moves.h ../lam/sa.h score.h maternal.h ../lam/random.h
score.o: score.c /usr/include/stdc-predef.h /usr/include/math.h \
 /usr/include/features.h /usr/include/sys/cdefs.h \
 /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h \
 /usr/include/gnu/stubs-64.h /usr/include/bits/huge_val.h \
 /usr/include/bits/huge_valf.h /usr/include/bits/huge_vall.h \
 /usr/include/bits/inf.h /usr/include/bits/nan.h \
 /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h \
 /usr/include/stdlib.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/stddef.h \
 /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h \
 /usr/include/endian.h /usr/include/bits/endian.h \
 /usr/include/bits/byteswap.h /usr/include/bits/types.h \
 /usr/include/bits/typesizes.h /usr/include/bits/byteswap-16.h \
 /usr/include/sys/types.h /usr/include/time.h /usr/include/sys/select.h \
 /usr/include/bits/select.h /usr/include/bits/sigset.h \
 /usr/include/bits/time.h /usr/include/sys/sysmacros.h \
 /usr/include/bits/pthreadtypes.h /usr/include/alloca.h \
 /usr/include/bits/stdlib-float.h /usr/include/ctype.h \
 /usr/include/xlocale.h /usr/include/string.h /usr/include/sys/time.h \
 ../lam/error.h ../lam/global.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/float.h integrate.h \
 /usr/include/stdio.h /usr/include/libio.h /usr/include/_G_config.h \
 /usr/include/wchar.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/stdarg.h \
 /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h maternal.h \
 /usr/include/gsl/gsl_errno.h /usr/include/errno.h \
 /usr/include/bits/errno.h /usr/include/linux/errno.h \
 /usr/include/asm/errno.h /usr/include/asm-generic/errno.h \
 /usr/include/asm-generic/errno-base.h /usr/include/gsl/gsl_types.h \
 /usr/include/gsl/gsl_spline.h /usr/include/gsl/gsl_interp.h \
 /usr/include/gsl/gsl_inline.h score.h zygotic.h \
 /usr/include/sys/resource.h /usr/include/bits/resource.h solvers.h
scramble.o: scramble.c /usr/include/stdc-predef.h /usr/include/stdlib.h \
 /usr/include/features.h /usr/include/sys/cdefs.h \
 /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h \
 /usr/include/gnu/stubs-64.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/stddef.h \
 /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h \
 /usr/include/endian.h /usr/include/bits/endian.h \
 /usr/include/bits/byteswap.h /usr/include/bits/types.h \
 /usr/include/bits/typesizes.h /usr/include/bits/byteswap-16.h \
 /usr/include/sys/types.h /usr/include/time.h /usr/include/sys/select.h \
 /usr/include/bits/select.h /usr/include/bits/sigset.h \
 /usr/include/bits/time.h /usr/include/sys/sysmacros.h \
 /usr/include/bits/pthreadtypes.h /usr/include/alloca.h \
 /usr/include/bits/stdlib-float.h /usr/include/stdio.h \
 /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/stdarg.h \
 /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h \
 /usr/include/string.h /usr/include/xlocale.h /usr/include/unistd.h \
 /usr/include/bits/posix_opt.h /usr/include/bits/environments.h \
 /usr/include/bits/confname.h /usr/include/getopt.h moves.h \
 ../lam/global.h /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/float.h \
 ../lam/sa.h score.h maternal.h ../lam/error.h zygotic.h integrate.h \
 /usr/include/gsl/gsl_errno.h /usr/include/errno.h \
 /usr/include/bits/errno.h /usr/include/linux/errno.h \
 /usr/include/asm/errno.h /usr/include/asm-generic/errno.h \
 /usr/include/asm-generic/errno-base.h /usr/include/gsl/gsl_types.h \
 /usr/include/gsl/gsl_spline.h /usr/include/gsl/gsl_interp.h \
 /usr/include/gsl/gsl_inline.h /usr/include/sys/resource.h \
 /usr/include/bits/resource.h
solvers.o: solvers.c /usr/include/stdc-predef.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/float.h \
 /usr/include/math.h /usr/include/features.h /usr/include/sys/cdefs.h \
 /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h \
 /usr/include/gnu/stubs-64.h /usr/include/bits/huge_val.h \
 /usr/include/bits/huge_valf.h /usr/include/bits/huge_vall.h \
 /usr/include/bits/inf.h /usr/include/bits/nan.h \
 /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h \
 /usr/include/stdlib.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/stddef.h \
 /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h \
 /usr/include/endian.h /usr/include/bits/endian.h \
 /usr/include/bits/byteswap.h /usr/include/bits/types.h \
 /usr/include/bits/typesizes.h /usr/include/bits/byteswap-16.h \
 /usr/include/sys/types.h /usr/include/time.h /usr/include/sys/select.h \
 /usr/include/bits/select.h /usr/include/bits/sigset.h \
 /usr/include/bits/time.h /usr/include/sys/sysmacros.h \
 /usr/include/bits/pthreadtypes.h /usr/include/alloca.h \
 /usr/include/bits/stdlib-float.h /usr/include/stdio.h \
 /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/stdarg.h \
 /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h \
 /usr/include/string.h /usr/include/xlocale.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/limits.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/syslimits.h \
 /usr/include/limits.h /usr/include/bits/posix1_lim.h \
 /usr/include/bits/local_lim.h /usr/include/linux/limits.h \
 /usr/include/bits/posix2_lim.h ../lam/error.h ../lam/global.h solvers.h \
 maternal.h zygotic.h integrate.h /usr/include/gsl/gsl_errno.h \
 /usr/include/errno.h /usr/include/bits/errno.h \
 /usr/include/linux/errno.h /usr/include/asm/errno.h \
 /usr/include/asm-generic/errno.h /usr/include/asm-generic/errno-base.h \
 /usr/include/gsl/gsl_types.h /usr/include/gsl/gsl_spline.h \
 /usr/include/gsl/gsl_interp.h /usr/include/gsl/gsl_inline.h \
 /usr/include/sys/resource.h /usr/include/bits/resource.h
translate.o: translate.c /usr/include/stdc-predef.h /usr/include/ctype.h \
 /usr/include/features.h /usr/include/sys/cdefs.h \
 /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h \
 /usr/include/gnu/stubs-64.h /usr/include/bits/types.h \
 /usr/include/bits/typesizes.h /usr/include/endian.h \
 /usr/include/bits/endian.h /usr/include/bits/byteswap.h \
 /usr/include/bits/byteswap-16.h /usr/include/xlocale.h \
 /usr/include/stdlib.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/stddef.h \
 /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h \
 /usr/include/sys/types.h /usr/include/time.h /usr/include/sys/select.h \
 /usr/include/bits/select.h /usr/include/bits/sigset.h \
 /usr/include/bits/time.h /usr/include/sys/sysmacros.h \
 /usr/include/bits/pthreadtypes.h /usr/include/alloca.h \
 /usr/include/bits/stdlib-float.h /usr/include/string.h ../lam/error.h \
 ../lam/global.h /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/float.h \
 maternal.h /usr/include/stdio.h /usr/include/libio.h \
 /usr/include/_G_config.h /usr/include/wchar.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/stdarg.h \
 /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h moves.h \
 ../lam/sa.h score.h zygotic.h integrate.h /usr/include/gsl/gsl_errno.h \
 /usr/include/errno.h /usr/include/bits/errno.h \
 /usr/include/linux/errno.h /usr/include/asm/errno.h \
 /usr/include/asm-generic/errno.h /usr/include/asm-generic/errno-base.h \
 /usr/include/gsl/gsl_types.h /usr/include/gsl/gsl_spline.h \
 /usr/include/gsl/gsl_interp.h /usr/include/gsl/gsl_inline.h \
 /usr/include/sys/resource.h /usr/include/bits/resource.h
unfold.o: unfold.c /usr/include/stdc-predef.h /usr/include/stdlib.h \
 /usr/include/features.h /usr/include/sys/cdefs.h \
 /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h \
 /usr/include/gnu/stubs-64.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/stddef.h \
 /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h \
 /usr/include/endian.h /usr/include/bits/endian.h \
 /usr/include/bits/byteswap.h /usr/include/bits/types.h \
 /usr/include/bits/typesizes.h /usr/include/bits/byteswap-16.h \
 /usr/include/sys/types.h /usr/include/time.h /usr/include/sys/select.h \
 /usr/include/bits/select.h /usr/include/bits/sigset.h \
 /usr/include/bits/time.h /usr/include/sys/sysmacros.h \
 /usr/include/bits/pthreadtypes.h /usr/include/alloca.h \
 /usr/include/bits/stdlib-float.h /usr/include/stdio.h \
 /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/stdarg.h \
 /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h \
 /usr/include/string.h /usr/include/xlocale.h /usr/include/unistd.h \
 /usr/include/bits/posix_opt.h /usr/include/bits/environments.h \
 /usr/include/bits/confname.h /usr/include/getopt.h \
 /usr/include/sys/resource.h /usr/include/bits/resource.h ../lam/error.h \
 ../lam/global.h /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/float.h \
 integrate.h maternal.h /usr/include/gsl/gsl_errno.h /usr/include/errno.h \
 /usr/include/bits/errno.h /usr/include/linux/errno.h \
 /usr/include/asm/errno.h /usr/include/asm-generic/errno.h \
 /usr/include/asm-generic/errno-base.h /usr/include/gsl/gsl_types.h \
 /usr/include/gsl/gsl_spline.h /usr/include/gsl/gsl_interp.h \
 /usr/include/gsl/gsl_inline.h solvers.h zygotic.h score.h
unstable-manifold.o: unstable-manifold.c /usr/include/stdc-predef.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/float.h \
 /usr/include/stdlib.h /usr/include/features.h /usr/include/sys/cdefs.h \
 /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h \
 /usr/include/gnu/stubs-64.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/stddef.h \
 /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h \
 /usr/include/endian.h /usr/include/bits/endian.h \
 /usr/include/bits/byteswap.h /usr/include/bits/types.h \
 /usr/include/bits/typesizes.h /usr/include/bits/byteswap-16.h \
 /usr/include/sys/types.h /usr/include/time.h /usr/include/sys/select.h \
 /usr/include/bits/select.h /usr/include/bits/sigset.h \
 /usr/include/bits/time.h /usr/include/sys/sysmacros.h \
 /usr/include/bits/pthreadtypes.h /usr/include/alloca.h \
 /usr/include/bits/stdlib-float.h /usr/include/stdio.h \
 /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/stdarg.h \
 /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h \
 /usr/include/string.h /usr/include/xlocale.h /usr/include/unistd.h \
 /usr/include/bits/posix_opt.h /usr/include/bits/environments.h \
 /usr/include/bits/confname.h /usr/include/getopt.h /usr/include/math.h \
 /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h \
 /usr/include/bits/huge_vall.h /usr/include/bits/inf.h \
 /usr/include/bits/nan.h /usr/include/bits/mathdef.h \
 /usr/include/bits/mathcalls.h /usr/include/sys/resource.h \
 /usr/include/bits/resource.h ../lam/error.h ../lam/global.h integrate.h \
 maternal.h /usr/include/gsl/gsl_errno.h /usr/include/errno.h \
 /usr/include/bits/errno.h /usr/include/linux/errno.h \
 /usr/include/asm/errno.h /usr/include/asm-generic/errno.h \
 /usr/include/asm-generic/errno-base.h /usr/include/gsl/gsl_types.h \
 /usr/include/gsl/gsl_spline.h /usr/include/gsl/gsl_interp.h \
 /usr/include/gsl/gsl_inline.h solvers.h zygotic.h score.h
zygotic.o: zygotic.c /usr/include/stdc-predef.h /usr/include/ctype.h \
 /usr/include/features.h /usr/include/sys/cdefs.h \
 /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h \
 /usr/include/gnu/stubs-64.h /usr/include/bits/types.h \
 /usr/include/bits/typesizes.h /usr/include/endian.h \
 /usr/include/bits/endian.h /usr/include/bits/byteswap.h \
 /usr/include/bits/byteswap-16.h /usr/include/xlocale.h \
 /usr/include/stdio.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/stddef.h \
 /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h \
 /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/stdarg.h \
 /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h \
 /usr/include/stdlib.h /usr/include/bits/waitflags.h \
 /usr/include/bits/waitstatus.h /usr/include/sys/types.h \
 /usr/include/time.h /usr/include/sys/select.h /usr/include/bits/select.h \
 /usr/include/bits/sigset.h /usr/include/bits/time.h \
 /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h \
 /usr/include/alloca.h /usr/include/bits/stdlib-float.h \
 /usr/include/string.h /usr/include/math.h /usr/include/bits/huge_val.h \
 /usr/include/bits/huge_valf.h /usr/include/bits/huge_vall.h \
 /usr/include/bits/inf.h /usr/include/bits/nan.h \
 /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h ../lam/error.h \
 ../lam/global.h /usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/float.h \
 maternal.h solvers.h integrate.h /usr/include/gsl/gsl_errno.h \
 /usr/include/errno.h /usr/include/bits/errno.h \
 /usr/include/linux/errno.h /usr/include/asm/errno.h \
 /usr/include/asm-generic/errno.h /usr/include/asm-generic/errno-base.h \
 /usr/include/gsl/gsl_types.h /usr/include/gsl/gsl_spline.h \
 /usr/include/gsl/gsl_interp.h /usr/include/gsl/gsl_inline.h zygotic.h \
 /usr/include/sys/resource.h /usr/include/bits/resource.h
