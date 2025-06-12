#!/bin/bash
make clean
make
gcc -g -o test test.c -L. -lmceliece8192128_clean -I../../../common/
./test
