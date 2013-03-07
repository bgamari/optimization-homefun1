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

def newton_step(A, b, x):
    """ Newton's method for linear system """
    return x - dot(np.linalg.inv(A), dot(A, x) - b)


from matplotlib import pyplot as pl
for i in range(100):
    n = 10
    A,b = generate_sys(n,n)
    print 'condition number', np.linalg.cond(A)
    x0 = np.random.randn(n)
    print 'start', x0

    x = x0.copy()
    gd_resids = []
    for i in range(50):
        gd_resids.append(norm(dot(A, x) - b))
        x = gradient_descent_step(A, b, x)
    print 'gradient descent', x

    x = x0.copy()
    newton_resids = []
    for i in range(50):
        newton_resids.append(norm(dot(A, x) - b))
        x = newton_step(A, b, x)
    print 'newton', x

    x,cg_resids = conj_gradient(A, b, x0, eps=0)
    print 'conjugate gradient', x

    pl.plot(gd_resids, '+', label='steepest descent', c='r')
    pl.plot(cg_resids, '+', label='conj. gradient', c='g')
    pl.plot(newton_resids, '+', label='newton', c='b')

pl.legend([])
pl.yscale('log')
pl.xlabel('iteration')
pl.ylabel('residual')
pl.savefig('q4-convergence.pdf')
