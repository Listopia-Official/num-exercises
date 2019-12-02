"""
This is the solution of Florian Haas (3382958) and Pascal Bauer (3383821) for
the fifth exercise (first programming exercise) on the second worksheet of Numerik I.
"""

import numpy as np
import matplotlib.pyplot as plt

"""
Computes the divided differences for the points supplied to the function,
and returns the coefficients for the interpolation polynome with the Newton base.
For that, it creates and stores a scheme simliar to that specified in the lecture, but
mirrorwd and turned, so that the first row are the coefficients.
Input: An array containing the points (x_i, f_i)
Output: An array containing the coefficients c_0, ..., c_k
"""
def newtonCoefficients(points):
    pointCount = len(points) # Helper variable
    
    pointMap = {point[0]:point[1] for point in points} # A map containing the x_i as key and the f_i as value
    
    knownDiffs = {} # Empty map where the known computed, divided differences are stored
    
    c = [] # The coefficient list
    
    for i in range(pointCount):
        for j in range(pointCount-i): # The scheme is in triangular form, so we don't need to iterate through the whole square
            
            arg = tuple([points[k+j][0] for k in range(i+1)]) # Compute the argument [x_k, ..., x_j] of the divided difference
            
            dividedDifference = diff(pointMap, arg, knownDiffs) # Compute the divided difference
            knownDiffs[arg] = dividedDifference # Add the difference to the map of known ones
            
            if j is 0:
                c.append(dividedDifference) # If we're on the first column, we have a coefficient c_i - store that
            
    return c

"""
A helper function which computes the divided difference f[x_k, ..., x_j] like specified in the lecture.
Input: 
    pointMap: A map containing the x_i as key and the f_i as value
    arg: A list containing the x_k, ..., x_j (the argument of the divided difference)
    knownDiffs: A map containing the argument [x_i, ..., x_l] as key and the divided difference f[x_i, ..., x_l] as value
Output:
    The value of the divided difference (a scalar value).
"""
def diff(pointMap, arg, knownDiffs):
    if(len(arg) == 1): # One argument
        return pointMap[arg[0]] # f[x_i] = f_i
    else: # At least two arguments
        """
        Computes f[x_(k+1), ..., x_j] - f[x_k, ..., x_(j-1)] / (x_j - x_k)
        """
        return (knownDiffs[arg[1:len(arg)]] - knownDiffs[arg[:len(arg)-1]])/(arg[len(arg)-1] - arg[0])
    
"""
Computes the Newton base polynome functions as specified on the lecture.
Input:
    points: An array containing the points relevant for the base
    
Output:
    An array containing the Newton base functions as lambdas
"""
def newtonBase(points):
    base = []
    for j in range(len(points)):
        if j == 0:
            base.append(lambda x: 1) # The 0-th base is p(x) = 1
        else:
            """
                p(x) = product from k=0 to j-1 of (x-x_k), as specified in lecture
            """
            base.append(lambda x, j=j, points=points: np.prod(np.array([x-points[i][0] for i in range(j)])))
    
    return base
    
"""
Plots the interpolatio polynome and the specified points.
Input:
    points: The points to plot (x_i, f_i)
    base: The Newton base
    coefficients: The corresponding Newton coefficients
"""
def plotSolution(points, polynome, function):    
    minX = min([points[i][0] for i in range(len(points))]) # Compute plot area from points
    maxX = max([points[i][0] for i in range(len(points))])
    
    for point in points:
        plt.scatter(point[0], point[1]) # Plot the points
        
    x = np.linspace(minX, maxX, 1000) # x-Axis
        
    plt.plot(x, [polynome(x_comp) for x_comp in x]) # Display the interpolation polynome
    plt.plot(x, [function(x_comp) for x_comp in x]) # Display the interpolation polynome
    
    plt.show()
    
def createInterpolationPolynome(deg, function, lower, upper):
    points = []
    delta = (upper- lower) / (deg+1)
    for i in range(deg+1):
        x = lower + i * delta
        points.append((x,function(x)))    
    base = newtonBase(points)
    coefficients = newtonCoefficients(points)
    return (points, lambda x, base=base, coefficients=coefficients: np.sum([base[i](x) * coefficients[i] for i in range(len(coefficients))]))

def computeErrors(degs, function, lower, upper, spaceGrain=1000):
    errors = []
    
    space = np.linspace(lower, upper, spaceGrain)
    
    for deg in degs:
        tple = createInterpolationPolynome(deg,function,lower,upper)
        errors.append(np.max([np.abs(tple[1](x) - function(x)) for x in space]))
        plotSolution(tple[0], tple[1],function)
    
    return errors
    
degs = [0, 4, 8, 16, 32]

print("sin(2*pi*x):",computeErrors(degs, lambda x : np.sin(2*np.pi*x), 0, 1))
print("(1+x^2)^-1:",computeErrors(degs, lambda x : 1/(1+x**2), -5, 5))
print("|x|:",computeErrors(degs, lambda x : np.abs(x), -1, 1))
