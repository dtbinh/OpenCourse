l1tf: A Matlab solver for l1 trend filtering.

Version 0.7 , Apr 23 2007

--------------------------------------------------------------------------------

Please see the comments in the source code, for more information.

--------------------------------------------------------------------------------
Files of interest:

SOLVERS
    l1tf.m         : main solver
    l1tf_cvx.m     : solves l1 trend filtering problem using CVX package

EXAMPLES
    ex_basic.m     : example for the basic usage
    ex_snp500.m    : example for S&P500 data    (used in the paper)
    ex_synthetic.m : example for synthetic data (used in the paper)

UTILITY
    l1tf_lambdamax : utility file for finding maximum value of lambda
    pwl_fit        : piecewise linear fitting utility

ETC
    plot_signals   : utility for plotting original and filtered signals
                     (called from ex_snp500.m)
    snp500.dat     : Standard & Poor data file



COPYRIGHT (c) 2007  Kwangmoo Koh, Seung-Jean Kim and Stephen Boyd.

Permission to use, copy, modify, and distribute this software for
any purpose without fee is hereby granted, provided that this entire
notice is included in all copies of any software which is or includes
a copy or modification of this software and in all copies of the
supporting documentation for such software.
This software is being provided "as is", without any express or
implied warranty.  In particular, the authors do not make any
representation or warranty of any kind concerning the merchantability
of this software or its fitness for any particular purpose.

