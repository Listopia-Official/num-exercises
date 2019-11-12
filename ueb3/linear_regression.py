"""
This is the solution of Florian Haas (3382958) for
the fourth exercise (first programming exercise) on the third worksheet of Numerik I.

Firstly, it assembles the Matrix V as specified in the exercise instructions and computed it's condition number
dependent on a specified n, then it solves the problem of linear regression for pairs (x, y) in R^2.
"""

import numpy as np
import matplotlib.pyplot as plot

"""
Returns the matrix V as specified in the exercise instructions.

Input: A tuple x = (x0, ..., xn)^T
"""
def assembleVMatrix(x):
    return np.array([[1, x[i]] for i in range(len(x))])

# Computes conditional number of the supplied matrix
def cond(A):
    return np.linalg.cond(A)

# Implements the test as specified in a)
def testCond(n):
    V = assembleVMatrix(np.linspace(0,1,n)) # Create V with linspace(0,1,n)
    print("n=",n," with cond(V)=",cond(V)) # Print conditional number
    
"""
Implementation of the linear regression as specified in b)

Input: A list of pairs in R^2
Output: A tuple (a_1, a_0) containing the coefficients of the linear solution polynome
"""
def linearRegression(pairs):
    x = np.asarray([[pairs[i][0]] for i in range(len(pairs))]) # x in R^n
    y = np.asarray([[pairs[i][1]] for i in range(len(pairs))]) # y in R^n
    V = assembleVMatrix(np.matrix.transpose(x)[0]) # V in R^(n x 2), we use [0] because transpose returns an array in an array
    V_T = np.matrix.transpose(V) # V^T in R^(2 x n)
    a = np.linalg.solve(np.matmul(V_T,V), np.matmul(V_T,y)) # a in R^n, solves the LES as specified in the exercise
    return a

# Creates the regression and displays it in comparison to the data points
def plotSolution(pairs):
    solution = linearRegression(pairs)
    
    print("Solution for the pairs", pairs, ": p(x) = ",solution[1][0],"*x + ",solution[0][0])
    
    minX = min([pairs[i][0] for i in range(len(pairs))])
    maxX = max([pairs[i][0] for i in range(len(pairs))])
    
    for pair in pairs:
        plot.scatter(pair[0], pair[1])
        
    x = np.linspace(minX, maxX, 1000)
        
    plot.plot(x, [solution[1][0]*x_comp + solution[0][0] for x_comp in x])
    plot.show()
    
  
# Print cond(V) for different n
testCond(5)
testCond(10)
testCond(20)
testCond(50)

print()

# The solution for b)
plotSolution([(0,0), (1, 10), (2, 10), (3, 20)])

# Another test case
#plotSolution([(0,0), (1,2), (2,1), (3,3), (4,5), (5,4), (6,6), (7,7), (8,9)])
