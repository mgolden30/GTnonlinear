CC=gcc
CFLAGS= -lgmp -lmpfr -I./include/

all: object_files tests main
object_files: obj/big_float.o obj/bf_work_mem.o obj/bf_blas.o obj/bf_nonlinear.o
tests: tests/blas_test tests/bf_nonlinear_test 

main: 3body/main.c $(object_files)
	$(CC) 3body/main.c obj/* -o $@ $(CFLAGS)

tests/blas_test: test_src/blas_test.c $(object_files)
	$(CC) test_src/blas_test.c obj/* -o $@ $(CFLAGS)

tests/bf_nonlinear_test: test_src/bf_nonlinear_test.c $(object_files)
	$(CC) test_src/bf_nonlinear_test.c obj/* -o $@ $(CFLAGS)


obj/bf_nonlinear.o: src/bf_nonlinear.c include/bf_nonlinear.h 
	$(CC) src/bf_nonlinear.c -c -o $@ $(CFLAGS)

obj/big_float.o: src/big_float.c include/big_float.h
	$(CC) src/big_float.c -c -o $@ $(CFLAGS)

obj/bf_work_mem.o: src/bf_work_mem.c include/bf_work_mem.h
	$(CC) src/bf_work_mem.c -c -o $@ $(CFLAGS) 

obj/bf_blas.o: src/bf_blas.c include/bf_blas.h
	$(CC) src/bf_blas.c -c -o $@ $(CFLAGS)






clean:
	rm obj/* tests/* main
