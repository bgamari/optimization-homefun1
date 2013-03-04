#/usr/bin/python

import numpy as np
import numpy.ma as ma

n = 100
M = np.random.random_sample((n,n))

def observations(nobs, ndim):
    a = ma.masked_all((n,n), dtype=float)
    for i in range(nobs):
        (i,j) = np.random.randint(0, n, 2)
        a[(i,j)] = np.random.random_sample()
        a[(j,i)] = a[(i,j)]
    return a

def proj_m_to_x(m):
    """ Project from symmetric matrix M to positive semi-definite matrix X """
    w,v = np.linalg.eigh(m)
    w[w < 0] = 0
    return np.dot(np.diag(w), np.dot(v, v.T))

def proj_x_to_m(obs, x):
    """ Project from positive semi-definite matrix X to completion matrix M """
    xx = x.copy()
    xx[~obs.mask] = obs[~obs.mask]
    return xx

obs = observations(10, 100)
m = obs.filled(0)
err = []
for i in range(100):
    x = proj_m_to_x(m)
    new_m = proj_x_to_m(obs, x)
    e = np.sum((m - new_m)**2)
    m = new_m

    err.append(e)
    print i, e

import matplotlib.pyplot as pl
pl.plot(err)
pl.xlabel('iteration')
pl.ylabel('Froebenius norm')
pl.savefig('q3-convergence.pdf')
