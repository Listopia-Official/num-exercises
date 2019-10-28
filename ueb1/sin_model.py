# Demo codes to perform all numerical experiments in chapter 1 of my lecture notes 
# "Numerische Mathematik 1". 
#
# Tested on the Linux console with Python 3.7.3, NumPy and MatPlotLib from Anaconda3-2019.03, 
# but any recent Python, NumPy and MatPlotLib should do the trick.
#
# Note that this file is optimised for "lecturer mode" and thus my console background to 
# programming. This programming style may be useful for scripted experiments while 
# writing a thesis, but is not didactically optimal for you in the sense of using existing
# demo codes right out of the box. See the preface of the brilliant textbook by Svein Linge
# and Hans Petter Langtangen, "Programming for Computations -- Python", Springer (2016)
# for an explanation of the underlying didactical background. Unfortunately, the code samples
# in this textbook have not been updated after Hans Petter passed away much too early, in
# particular related to the print() command. Some caution is advised. 
#
# Notes on diactics:
#  (1) Python allows many approaches, and mine is perfectly legitimate.
#  (2) You are forced to understand, rather than only use, the code when you insert it
#      into your favourite Python environment to reproduce the computer experiments.
#  (3) This extends to the programming excercises for this lecture, which are much easier to 
#      master once you understand the demos: Our programming excercises mostly comprise 
#      performing your own experiments with a given scheme, rather than implementing it
#      yourself from scratch.
# Assuming some fluency in Python and in mastering your programming environment, only minor
# additions and modifications to the provided code samples are necessary.
#
# Scroll to the end of the file to comment in/out selected functionality.
#
# Execution: python3 chap1.py
#
# (c) 2019-20, dominik.goeddeke@mathematik.uni-stuttgart.de


#######################################################
# imports
#######################################################

import matplotlib.pyplot as plt    # for plotting
import numpy as np                 # for linear algebra


#######################################################
# assembly of the linear system
#######################################################

def poissonmatrix(n) :
  # n: number of inner unknowns in the finite difference discretisation

  # Note that is extremely inefficient to store a tridiagonal matrix as a dense one.
  # We set the datatype explicitly in preparation for later experiments
  A = np.zeros((n,n), dtype=np.float64)  
  for i in range(0,n,1) :
    A[i,i] = 2.
  for i in range(1,n,1) :
    A[i,i-1] = -1.
  for i in range(0,n-1,1) :
    A[i,i+1] = -1.
  return A


def rhs(n,h,f,ga,gb) :
  # n: number of inner unknowns in the finite difference discretisation
  # h: corresponding mesh width
  # f: force function, assumed constant
  # ga, gb: Dirichlet boundary values of the solution at the left and right interval boundary
  b = np.zeros(n, dtype=np.float64)
  for i in range(0,n,1) :
    b[i] = f * h**2
  b[0] += ga
  b[n-1] += gb
  return b
  

def solve_problem(n,a,b,ga,gb,f) :
  # n: number of inner unknowns in the finite difference discretisation
  # a,b: interval ]a,b[
  # ga,gb: Dirichlet boundary values of the solution at the left and right interval boundary
  # f: force function, assumed constant
  h = (b-a)/(n+1)
  A = poissonmatrix(n)
  f = rhs(n,h,f,ga,gb)
  u = np.linalg.solve(A,f)
  # expand solution to include Dirichlet values
  u = np.concatenate(([ga],u))
  u = np.concatenate((u,[gb]))
  return u
  
  

#######################################################
# functions to generate all plots in the first part of chapter 1
#######################################################

def plot1():
  n = 100
  u = solve_problem(n, -3,3, 0,0, -1)
  x = np.linspace(-3, 3, n+2)   # this generates the discretisation points
  plt.plot(x,u, color='blue', linestyle='-', linewidth=2, markersize=2)
  u = solve_problem(n, -2,4, 0,0, -1)
  x = np.linspace(-2, 4, n+2)
  plt.plot(x,u, color='green', linestyle='-', linewidth=2, markersize=2)
  u = solve_problem(n, -1,5, 0,0, -1)
  x = np.linspace(-1, 5, n+2)
  plt.plot(x,u, color='red', linestyle='-', linewidth=2, markersize=2)
  u = solve_problem(n, 0,6, 0,0, -1)
  x = np.linspace(0, 6, n+2)
  plt.plot(x,u, color='cyan', linestyle='-', linewidth=2, markersize=2)
  plt.savefig('simple_examples_1.pdf')
  plt.close()

def plot2():
  n = 100
  u = solve_problem(n, 0,1, 0,0, -1)
  x = np.linspace(0, 1, n+2)
  plt.plot(x,u, color='blue', linestyle='-', linewidth=2, markersize=2)
  u = solve_problem(n, -1,2, 0,0, -1)
  x = np.linspace(-1, 2, n+2)
  plt.plot(x,u, color='green', linestyle='-', linewidth=2, markersize=2)
  u = solve_problem(n, -2,3, 0,0, -1)
  x = np.linspace(-2, 3, n+2)
  plt.plot(x,u, color='red', linestyle='-', linewidth=2, markersize=2)
  u = solve_problem(n, -3,4, 0,0, -1)
  x = np.linspace(-3, 4, n+2)
  plt.plot(x,u, color='cyan', linestyle='-', linewidth=2, markersize=2)
  plt.savefig('simple_examples_2.pdf')
  plt.close()

def plot3():
  n = 100
  a = -3
  b = 3
  u = solve_problem(n, a,b, 0,0, -1)
  x = np.linspace(a, b, n+2)
  plt.plot(x,u, color='blue', linestyle='-', linewidth=2, markersize=2)
  u = solve_problem(n, a,b, 0,0, -0.75)
  x = np.linspace(a, b, n+2)
  plt.plot(x,u, color='green', linestyle='-', linewidth=2, markersize=2)
  u = solve_problem(n, a,b, 0,0, -0.5)
  x = np.linspace(a, b, n+2)
  plt.plot(x,u, color='red', linestyle='-', linewidth=2, markersize=2)
  u = solve_problem(n, a,b, 0,0, -0.25)
  x = np.linspace(a, b, n+2)
  plt.plot(x,u, color='cyan', linestyle='-', linewidth=2, markersize=2)
  plt.savefig('simple_examples_3.pdf')
  plt.close()

def plot4():
  n = 100
  a = -3
  b = 3
  u = solve_problem(n, a,b, 0,0, -1)
  x = np.linspace(a, b, n+2)
  plt.plot(x,u, color='blue', linestyle='-', linewidth=2, markersize=2)
  u = solve_problem(n, a,b, 0,1, -1)
  x = np.linspace(a, b, n+2)
  plt.plot(x,u, color='green', linestyle='-', linewidth=2, markersize=2)
  u = solve_problem(n, a,b, 0,5, -1)
  x = np.linspace(a, b, n+2)
  plt.plot(x,u, color='red', linestyle='-', linewidth=2, markersize=2)
  u = solve_problem(n, a,b, 0,10, -1)
  x = np.linspace(a, b, n+2)
  plt.plot(x,u, color='cyan', linestyle='-', linewidth=2, markersize=2)
  plt.savefig('simple_examples_4.pdf')
  plt.close()

#######################################################
# functions to generate the table at the end of chapter 1
#######################################################

def u(x) :
  # evaluates u in x
  return (x*(1.0-x))**2
  
def uxx(x) :
  # evaluates u'' in x
  return 12.0*x**2 - 12.0*x + 2.0

def error_table() :
  err_norm = []
  for exp in range(1,11,1) :
    n = 2**exp
    a = 0
    b = 1
    h = (b-a)/(n+1)
    x = np.linspace(a+h, b-h, n)
    A = poissonmatrix(n)
    f = uxx(x) * (-1) * h**2
    exact_sol = u(x)
    approx_sol = np.linalg.solve(A,f)
    # no need to concatenate to full solution, as Dirichlet values are exact
    f = exact_sol - approx_sol
    err_norm.append(np.linalg.norm(f,np.inf))
    if exp == 1 :
      print(f'{n:4} {err_norm[exp-1]:6.4e}')
    else :
      print(f'{n:4} {err_norm[exp-1]:6.4e} {err_norm[exp-2]/err_norm[exp-1]:.2f}')
  
#######################################################
# toying with precisions
#######################################################  

def singleprec() :
  a64 = np.float64(1.0)
  b64 = np.float64(1e-8)
  c64 = np.float64(0.9998)
  d64 = np.float64(1.0002)
  e64 = np.float64(1e8)
  a32 = np.float32(1.00000000)
  b32 = np.float32(1e-8)
  c32 = np.float32(0.9998)
  d32 = np.float32(1.0002)
  e32 = np.float32(1e8)
  
  ab64 = a64+b64
  cd64 = c64*d64
  ab32 = a32+b32
  cd32 = c32*d32
  
  print(f'double: {a64:.8f} + {b64:.8f} = {ab64:.8f}')
  print(f'float : {a32:.8f} + {b32:.8f} = {ab32:.8f}')
  print(f'double: {c64:.4f} * {d64:.4f} = {cd64:.8f}')
  print(f'float : {c32:.4f} * {d32:.4f} = {cd32:.8f}')
  print(f'double: ({ab64:.8f} - {a64:.8f})*10^8 = {(ab64-a64)*e64:.8f}')
  print(f'float : ({ab32:.8f} - {a32:.8f})*10^8 = {(ab32-a32)*e32:.8f}')
  print(f'double: ({cd64:.8f} - {a64:.8f})*10^8 = {(cd64-a64)*e64:.8f}')
  print(f'float : ({cd32:.8f} - {a32:.8f})*10^8 = {(cd32-a32)*e32:.8f}')
  print(f'double: ({a64:.8f} - {a64:.8f} + {b64:.8f})*10^8 = {(a64-a64+b64)*e64:.8f}')
  print(f'float : ({a32:.8f} - {a32:.8f} + {b32:.8f})*10^8 = {(a32-a32+b32)*e32:.8f}')

#######################################################
# actual scripting area
#######################################################

#plot1()
#plot2()
#plot3()
#plot4()

## test case for -u''=0
u = solve_problem(100, 0,1, 10,12, lambda x: math.sin(x))
#x = np.linspace(0, 1, 100+2)
#plt.plot(x,u, color='blue', linestyle='-', linewidth=2, markersize=2)
#plt.savefig('test.pdf')

#singleprec()

error_table()



