function record_bfgs_iters(t, f, x, auxdata)
    global BFGS_ITERS;
    BFGS_ITERS = t;
end