#! /bin/sh
gcc -O3 -o l1tf l1tf.c main.c -lblas -llapack -lm
