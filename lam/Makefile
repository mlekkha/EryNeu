# gen_deviates
GDOBJ= deviates.o distributions.o error.o random.o

# LSA objectives
#

ifeq ($(MPI), on)
	LSAOBJ = lsa.o lsa-mpi.o
else
	LSAOBJ = lsa.o
endif	

# header files

LSA_HEADS = global.h sa.h MPI.h error.h
RND_HEADS = global.h random.h error.h
DIS_HEADS = global.h distributions.h error.h random.h

#targets

#all: gen_deviates deviates.o distributions.o error.o lsa.o random.o
all: gen_deviates $(GDOBJ) $(LSAOBJ)

gen_deviates: $(GDOBJ)
	$(CC) -o gen_deviates $(CFLAGS) $(GDOBJ) $(LIBS)


deviates.o:  $(DIS_HEADS) deviates.c 
	$(CC) $(CFLAGS) -c deviates.c -o deviates.o

distributions.o: $(DIS_HEADS) distributions.c
	$(CC) $(CFLAGS)  -c distributions.c 

error.o: error.c
	$(CC) $(CFLAGS) -c error.c -o error.o

lsa.o: $(LSA_HEADS) lsa.c
	$(CC) $(CFLAGS) -c lsa.c -o lsa.o

random.o: $(RND_HEADS) random.c
	$(CC) $(CFLAGS) -c random.c -o random.o

# parallel stuff

lsa-mpi.o: lsa.c
	$(MPICC) -c -o lsa-mpi.o $(MPIFLAGS) $(CFLAGS) lsa.c

# ... and here are the cleanup and make deps rules

clean:
	rm -f *.o core*

