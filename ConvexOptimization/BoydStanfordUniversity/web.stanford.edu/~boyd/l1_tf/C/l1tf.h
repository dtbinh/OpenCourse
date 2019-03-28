/* l1tf.h
 *
 * Copyright (C) 2007 Kwangmoo Koh, Seung-Jean Kim and Stephen Boyd.
 *
 * See the file "COPYING.TXT" for copyright and warranty information.
 *
 * Author: Kwangmoo Koh (deneb1@stanford.edu)
 */
#ifndef L1TF_H
#define L1TF_H

#ifdef __cplusplus
extern "C" {
#endif

/* main routine for l1 trend filtering */
int l1tf(const int n, const double *y, const double lambda, double *x);

/* utilit to compte the maximum value of lambda */
double l1tf_lambdamax(const int n, double *y);

/* utility to print a vector */
void print_dvec(int n, const double *x);

#ifdef __cplusplus
}
#endif

#endif /* L1TF_H */
