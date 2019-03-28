l1tf: A C solver for l1 trend filtering

COPYRIGHT (c) 2007  Kwangmoo Koh, Seung-Jean Kim and Stephen Boyd.

See the file "COPYING.TXT" for copyright and warranty information.

--------------------------------------------------------------------------------

Files of interest:

    l1tf.c     : main solver for l1 trend filtering
    main.c     : command-line wrapper for l1tf
    snp500.txt : time-series data (daily closing values of the S&P 500)

--------------------------------------------------------------------------------

Compilation:
    To compile the solver, you need gcc and BLAS/LAPACK. Then, run make.sh or
    Makefile.

    $ ./make.sh
    or
    $ make

--------------------------------------------------------------------------------

Usage example:
    To perform l1 trend filtering on snp500.txt for lambda = 50 and store
    the result to result.txt, you may want to run

        $ l1tf snp500.txt 50 > result.txt

    To perform l1 trend filtering on snp500.txt for lambda = 0.01*lambda_max
    and store the result to result.txt, you may want to run

        $ l1tf -r snp500.txt 0.01 > result.txt

