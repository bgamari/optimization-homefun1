#!/usr/bin/python

import numpy as np
from numpy import dot
from numpy.linalg import norm

def generate_sys(neq, nvar):
    """ Generate a linear system """
    A = np.random.randn(neq, nvar)
    A = dot(A.T, A)
    assert(np.all(np.linalg.eigvals(A) > 0))
    b = np.random.randn(neq)
    return (A,b)

n = 10
A,b = generate_sys(n,n)
print 'condition number', np.linalg.cond(A)
x0 = np.random.randn(n)
print 'start', x0

def gradient_descent_step(A, b, x):
    """ Steepest descent step for linear system """
    r = b - dot(A, x)
    t = dot(r,r) / dot(r, dot(A, r))
    x = x + t*r
    return x

def conj_gradient(A, b, x0, niter=50, eps=1e-12):
    """ Conjugate gradient for linear system """
    x = x0
    r = b - dot(A, x0)
    p = r
    resids = []
    for i in range(niter):
        a = dot(r,r) / dot(p, dot(A, p))
        x = x + a*p
        r_old = r
        r = r - a*dot(A, p)
        beta = dot(r,r) / dot(r_old,r_old)
        p = r + beta*p

        resids.append(norm(dot(A, x) - b))
        if norm(r) < eps:
            print 'Converged after %d iterations' % i
            break

    return x, resids

def newton(A, b, x):
    """ Newton's method for linear system """
    pass
    #return x -

x = x0.copy()
gd_resids = []
for i in range(50):
    x = gradient_descent_step(A, b, x)
    gd_resids.append(norm(dot(A, x) - b))
print 'gradient descent', x

x,cg_resids = conj_gradient(A, b, x0)
print 'conjugate gradient', x

from matplotlib import pyplot as pl
pl.plot(gd_resids, label='steepest descent')
pl.plot(cg_resids, label='conj. gradient')
pl.legend()
pl.yscale('log')
pl.xlabel('iteration')
pl.ylabel('residual')
pl.savefig('q4-convergence.pdf')
