CC=gcc
CFLAGS= -lgmp -lmpfr -I./include/

all: object_files tests
object_files: obj/big_float.o obj/bf_work_mem.o obj/bf_blas.o
tests: tests/blas_test 

tests/blas_test: test_src/blas_test.c
	$(CC) test_src/blas_test.c obj/* -o tests/blas_test $(CFLAGS)

obj/big_float.o: src/big_float.c include/big_float.h
	$(CC) src/big_float.c -c -o $@ $(CFLAGS)

obj/bf_work_mem.o: src/bf_work_mem.c include/bf_work_mem.h
	$(CC) src/bf_work_mem.c -c -o $@ $(CFLAGS) 

obj/bf_blas.o: src/bf_blas.c include/bf_blas.h
	$(CC) src/bf_blas.c -c -o $@ $(CFLAGS)



clean:
	rm obj/* tests/*
