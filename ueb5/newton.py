"""
This is the solution of Florian Haas (3382958) and Pascal Bauer (3383821) for
the fifth exercise (first programming exercise) on the second worksheet of Numerik I.
"""

import numpy as np
import matplotlib.pyplot as plt

def divDiff(points): # points: An array containing the f_i with key x_i
    pointCount = len(points)
    
    pointMap = {point[0]:point[1] for point in points}
    
    diffMatrix = np.zeros(shape=(pointCount,pointCount))
    
    knownDiffs = {}
    
    c = []
    
    for i in range(pointCount):
        for j in range(pointCount-i):
            arg = tuple([points[k+j][0] for k in range(i+1)])
            diff2 = diff(pointMap, arg, knownDiffs)
            knownDiffs[arg] = diff2
            diffMatrix[i][j] = diff2
            
            if j is 0:
                c.append(diff2)
            
    return c

def diff(pointMap, arg, knownDiffs):
    if(len(arg) == 1):
        return pointMap[arg[0]]
    else:
        return (knownDiffs[arg[1:len(arg)]] - knownDiffs[arg[:len(arg)-1]])/(arg[len(arg)-1] - arg[0])
    
def newtonBase(points):
    base = []
    for j in range(len(points)):
        if j == 0:
            base.append(lambda x: 1)
        else:
            base.append(lambda x, j=j, points=points: np.prod(np.array([x-points[i][0] for i in range(j)])))
    
    return base
    
def plotSolution(points, base, solution):    
    minX = min([points[i][0] for i in range(len(points))])
    maxX = max([points[i][0] for i in range(len(points))])
    
    for point in points:
        plt.scatter(point[0], point[1]) # Plot the points
        
    x = np.linspace(minX, maxX, 1000) # x-Axis
    
    polynome = lambda x : sum([base[i](x) * solution[i] for i in range(len(solution))])
        
    plt.plot(x, [polynome(x_comp) for x_comp in x]) # Display the interpolation
    plt.show()
    
points = [(-1, -1), (0,0), (1,3)] # Interpolation constraints
base = newtonBase(points) # The Newton base
coefficients = divDiff(points) # Coefficients of the Newton base

print("Points:",points)
print("Newton Base Coefficients:",coefficients)

plotSolution(points, base, coefficients)
