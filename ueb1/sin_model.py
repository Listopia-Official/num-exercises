"""
This is the solution of Florian Haas (3382958) and Pascal Bauer (3383821) for
the fourth exercise (first programming exercise) on the first worksheet of Numerik I.

This program solves the problem u''(x) = sin(2*pi*x) for x in (0,9) with the
Dirichlet boundary constraints u(0)=2 and u(9)=3.

Based on the provided template from Dominik Goeddeke, modified so that we can
use a non constant function; parts we didn't need were removed.

Copyright information from D. Goeddeke:
(c) 2019-20, dominik.goeddeke@mathematik.uni-stuttgart.de
"""

#######################################################
# imports
#######################################################

import matplotlib.pyplot as plt    # for plotting
import numpy as np                 # for linear algebra
import math                        # for ... math

#######################################################
# assembly of the linear system
#######################################################

""" 
n: number of inner unknowns in the finite difference discretisation

Note that is extremely inefficient to store a tridiagonal matrix as a dense one.
We set the datatype explicitly in preparation for later experiments
"""
def poissonmatrix(n) :
  A = np.zeros((n,n), dtype=np.float64)  
  for i in range(0,n,1) :
    A[i,i] = 2.
  for i in range(1,n,1) :
    A[i,i-1] = -1.
  for i in range(0,n-1,1) :
    A[i,i+1] = -1.
  return A


"""
n: number of inner unknowns in the finite difference discretisation
h: corresponding mesh width
f: function to evaluate (hasn't to be constant)
ga, gb: Dirichlet boundary values of the solution at the left and right interval boundary
"""
def rhs(n,h,f,ga,gb) :
  b = np.zeros(n, dtype=np.float64)
  for i in range(0,n,1) :
    b[i] = f(i*h) * h**2
  b[0] += ga
  b[n-1] += gb
  return b
  

"""
n: number of inner unknowns in the finite difference discretisation
a,b: interval ]a,b[
ga,gb: Dirichlet boundary values of the solution at the left and right interval boundary
f: function to evaluate (hasn't to be constant)
"""
def solve_problem(n,a,b,ga,gb,f) :
  h = (b-a)/(n+1)
  A = poissonmatrix(n)
  f = rhs(n,h,f,ga,gb)
  u = np.linalg.solve(A,f)
  # expand solution to include Dirichlet values
  u = np.concatenate(([ga],u))
  u = np.concatenate((u,[gb]))
  return u

#######################################################
# functions
#######################################################

 # evaluates u'' in x  
def uxx(x) :
  return math.sin(math.pi*2*x)

#######################################################
# actual scripting area
#######################################################

# test case for u''(x) = sin(2*pi*x) with the specified boundaries (see initial comment)
  
n = 500 # The amount of grid points (decrease when experiencing performance problems)
f = uxx

"""
Invokes the solver function to create the solution for the 
specified grid point count, grid width, interval, boundary values and function.
"""
u = solve_problem(n, 0,9, 2,3, f)

# Display the result
x = np.linspace(0, 9, n+2) # Create our x axis, match dimensions with y axis
plt.plot(x,u, color='blue', linestyle='-', linewidth=2, markersize=2)
plt.show()
