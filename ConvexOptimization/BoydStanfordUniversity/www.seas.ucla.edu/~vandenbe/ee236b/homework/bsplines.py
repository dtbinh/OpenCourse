import numpy as np

def bsplines(x):

# [y, yp, ypp] = bsplines(x)
#
# Computes the values and the first and second derivatives of the 13 cubic
# B-splines with breakpoints 0, 1, 2, ..., 10.
# x must be in the interval [0, 10].
#
# y:    a numpy array with the values of the 13 cubic B-splines evaluated
#       at x.
# yp:   a numpy array with the first derivatives of the 13 cubic B-splines
#       B-splines evaluated at x.
# ypp:  a numpy array with the second derivatives of the 13 cubic
#       B-splines evaluated at x.

    if x < 0 or x > 10:
        raise ValueError('x must be in [0,10].')

    # knots
    a = np.array([0,0,0,0,1,2,3,4,5,6,7,8,9,10,10,10,10])

    # 0 <= x < 10 and a(i) <= x < a(i+1) or x = 10 and i is 14.
    i = np.floor(x) + 4

    # y1, y2, y are vectors of linear, quadratic and cubic B-splines at x.
    if i == 14:
        y1 = np.zeros(11); y1[-1] = 1
        y2 = np.zeros(12); y2[-1] = 1
        y  = np.zeros(13);  y[-1] = 1
    else:
        j = int(i - 1); x = 1.0 * x;
        w1 = (x - a[j]) / (a[j+1] - a[j])
        G1 = np.array([1-w1, w1]).reshape(2,1)
        w2 = np.array([(x - a[j-1]) / (a[j+1] - a[j-1]), (x - a[j]) / (a[j+2] - a[j])])
        G2 = np.vstack((np.diag(1-w2), np.zeros((1,2)))) + np.vstack((np.zeros((1,2)), np.diag(w2)))
        w3 = np.array([(x - a[j-2]) / (a[j+1] - a[j-2]), (x - a[j-1]) / (a[j+2] - a[j-1]), (x - a[j]) / (a[j+3] - a[j])])
        G3 = np.vstack((np.diag(1-w3), np.zeros((1,3)))) + np.vstack((np.zeros((1,3)), np.diag(w3)))
        y  = np.zeros(13)
        y[j-3:j+1] = np.matmul(np.matmul(G3, G2), G1).reshape(4,)
        y2 = np.zeros(12)
        y2[j-3:j] = np.matmul(G2,G1).reshape(3,)
        y1 = np.zeros(11)
        y1[j-3:j-1] = G1.reshape(2,)

    # width of support of linear, quadratic, and cubic B-splines.
    d1 = a[4:15] - a[2:13]
    d2 = a[4:16] - a[1:13]
    d3 = a[4:17] - a[0:13]

    # derivatives of quadratic B-splines.
    y2p = np.zeros(12)
    y2p[0] = -2 * y1[0] / d1[0]
    y2p[1:11] = 2 * ( np.divide(y1[0:10], d1[0:10]) - np.divide(y1[1:11], d1[1:11]) )
    y2p[11] = 2.0 * y1[10] / d1[10]

    # derivatives of cubic B-spines
    yp = np.zeros(13)
    yp[0] = -3.0 * y2[0] / d2[0]
    yp[1:12] = 3 * ( np.divide(y2[0:11], d2[0:11]) - np.divide(y2[1:12], d2[1:12]) )
    yp[12] = 3.0 * y2[11] / d2[11]

    # second derivatives of cubic B-spines.
    ypp = np.zeros(13)
    ypp[0] = -3.0 * y2p[0] / d2[0]
    ypp[1:12] = 3 * ( np.divide(y2p[0:11], d2[0:11]) - np.divide(y2p[1:12], d2[1:12]) )
    ypp[12] = 3.0 * y2p[11] / d2[11]

    return y, yp, ypp
