#!/usr/bin/python

import numpy as np

def generate_sys(neq, nvar):
    """ Generate a linear system """
    A = np.random.randn(neq, nvar)
    A = np.dot(A.T, A)
    assert(np.all(np.linalg.eigvals(A) > 0))
    b = np.random.randn(neq)
    return (A,b)

A,b = generate_sys(10,10)

def gradient_descent_step(A, b, x):
    """ Steepest descent step for linear system """
    r = b - np.dot(A, x)
    t = np.dot(r,r) / np.dot(r, np.dot(A, r))
    x = x + t*r
    return x

def conj_gradient(A, b, x0, niter=50, eps=1e-10):
    """ Conjugate gradient step for linear system """
    x = x0
    r = b - np.dot(A, x0)
    d = r
    resids = []
    for i in range(niter):
        a = np.dot(r, d) / np.dot(d, np.dot(A, d))
        x = x + a*d
        r_old = r
        r = r - np.dot(a, np.dot(A, d))
        b = np.dot(r, r) / np.dot(r_old, r_old)
        d = r + b*d

        resids.append(np.linalg.norm(np.dot(A, x) - b))
        if np.linalg.norm(r) < eps:
            print 'Converged after %d iterations' % i
            break

    return x, resids

def newton(A, b, x):
    """ Newton's method for linear system """
    pass
    #return x - 
    
x0 = np.random.randn(10)
x = x0
gd_resids = []
for i in range(50):
    x = gradient_descent_step(A, b, x)
    gd_resids.append(np.linalg.norm(np.dot(A, x) - b))

_,cg_resids = conj_gradient(A, b, x0, eps=0)

from matplotlib import pyplot as pl
pl.plot(gd_resids, label='steepest descent')
pl.plot(cg_resids, label='conj. gradient')
pl.legend()
pl.yscale('log')
pl.xlabel('iteration')
pl.ylabel('residual')
pl.savefig('q4-convergence.pdf')

