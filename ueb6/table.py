"""
This is the solution of Florian Haas (3382958) and Pascal Bauer (3383821) for
the sixth exercise (first programming exercise) on the second worksheet of Numerik I.
"""

import numpy as np

###
### Old implementation start
###

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

###
### Old implementation until here - now comes the new part
###
    
"""
Input:
    deg: The degree of the interpolation polynome
    function: The original function which should be interpolated
    lower: The lower value of the compact interval (for [A,B] it would be A)
    upper: The upper value of the compact interval (for [A,B] it would be B)
    
Output:
    A tuple containing the interpolation points (a list containing pairs of x-values and f(x)) in the first part
    and the interpolation polynome as a function in the second part.
"""
def createInterpolationPolynome(deg, function, lower, upper):
    points = [] # A list containing deg+1 tuples (x, f(x)) (Stützpunkte)
    delta = (upper- lower) / (deg+1) # Equal distance between the x-values (Äquidistante Stützstellen)
    
    # Populate the point list
    for i in range(deg+1):
        x = lower + i * delta # Compute the current x-value
        points.append((x,function(x)))    # Store the tuple in the list
    
    base = newtonBase(points) # Create the newton base for the points
    coefficients = newtonCoefficients(points) # Compute the newton coefficients
    
    # Assemble the interpolation polynome
    interpolationPolynome = lambda x, base=base, coefficients=coefficients: np.sum([base[i](x) * coefficients[i] for i in range(len(coefficients))])
    
    return (points, interpolationPolynome) # Return tuple as specified in head docs

"""
Input:
    degs: A list containing the degrees of the interpolation polynomes
    function: The original function to interpolate
    lower: The lower interval bound (for [A,B] it would be A)
    upper: The upper interval bound (for [A,B] it would be B)
    pointCount: The count of points to use for the x-space when computing the unit norm
    
Output: A list containing the computed error value for each supplied degree - the first error value corresponds to the first degree and so on.
"""
def computeErrors(degs, function, lower, upper, pointCount=1000):
    errors = []
    space = np.linspace(lower, upper, pointCount) # Create the space
    
    for deg in degs: # Compute the error for each degree
        tple = createInterpolationPolynome(deg,function,lower,upper) # Compute the points and the interpolation polynome - this is a tuple!
        
        """
        np.abs(tple[1](x) - function(x)) is the absolute error between the interpolation and the real function value for each x
        np.max(...) computed the maximum value for the created array (we can use max instead of sup, because the unit norm in a compact interval is equal to the max norm).
        """
        errors.append(np.max([np.abs(tple[1](x) - function(x)) for x in space]))
    
    return errors
    
degs = [0, 4, 8, 16, 32] # The degrees as specified in the table

# Compute the error values for each function
sine_errors = computeErrors(degs, lambda x : np.sin(2*np.pi*x), 0, 1) # sin(2*pi*x)
hyp_errors = computeErrors(degs, lambda x : 1/(1+x**2), -5, 5) # 1/(1+x^2)
abs_errors = computeErrors(degs, lambda x : np.abs(x), -1, 1) # |x|

# Print the table

# Header
print("-"*72)
print("| {:>2} | {:s} | {:s} | {:s} |".format("n","sin(2*pi*x) auf [0, 1]","1/(1+x^2) auf [-5, 5]", "|x| auf [-1,1]"))
print("-"*72)

# Body
for i in range(len(degs)):
    print("| {:>2} | {:^22.4e} | {:^21.4e} | {:^14.4e} |".format(degs[i], sine_errors[i], hyp_errors[i], abs_errors[i]))

print("-"*72)

"""
Our computed values probably diverge from the values in the table in the script, because of the floating-point arithmetic.
(or undetected implementation errors)
"""
