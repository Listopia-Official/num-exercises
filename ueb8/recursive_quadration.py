"""
This is the solution of Florian Haas (3382958) and Pascal Bauer (3383821) for
the fourth exercise (first programming exercise) on the eigth worksheet of Numerik I.
"""

def unity_weights():
    return (1/6, 4/6, 1/6) # Weights for the Simpson-Rule in the unity interval [0, 1] (from the script)

"""
This function applies theorem 6.5 (Transformationssatz) - it transforms the supplied
unity weights to the weigths in the supplied interval [lower, upper];
additionally it returns the unity grid points in that interval.

The function returns a tuple, where the first value are the weigths, and the second value are the grid points.
"""
def transform(unity_weights, lower, upper):
    length = upper-lower # |[a,b]| = b-a
    h = length / (len(unity_weights)-1) # Compute the distance between two adjacent grid points
    
    # Returns the tuple as specified in the description
    return ([length * unity_weight for unity_weight in unity_weights], [lower + i*h for i in range(len(unity_weights))])

"""
This function is an implementation of the adaptive, recursive quadrature.

Input:
    function: The function to integrate
    unity_weights: The unity weights in the unity interval  [0, 1] of the quadrature to use
    lower: The lower bound of the integration interval [lower, upper]
    upper: The upper bound in the integration interval (see above)
    tolerance: The tolerance to use for comparison, defaults to 1e-4
    no_recursion: Internal parameter, if True, only the coarse integral value will be computed and immediately returned

Output:
    The integral value computed with the adaptive, recursive quadrature.
"""
def adaptive_rekursive_quadrature(function, unity_weights, lower, upper, tolerance = 1e-4, no_recursion=False):
    # Create the weights and grid points for the current interval
    transformed_weights, transformed_grid_points = transform(unity_weights, lower, upper)
    
    # Compute the coarse_value - just use the quadrature formula as in the script
    coarse_value = sum([transformed_weights[i] * function(transformed_grid_points[i]) for i in range(len(transformed_grid_points))])
    
    # Returns the coarse_value if no_recursion is True - the sence of this is explained below
    if no_recursion:
        return coarse_value
    
    center_point = (upper + lower)/2 # The center point of the currentl interval
    
    """
    Compute the fine value without recursion - the interval is divided in the middle, and for each interval of those two 
    the COARSE integral value will be computed with the call below.
    The sum of those two is considered the fine value.
    """
    fine_value = adaptive_rekursive_quadrature(function, unity_weights, lower, center_point, tolerance, no_recursion=True) + adaptive_rekursive_quadrature(function, unity_weights, center_point, upper, tolerance, no_recursion=True)
    
    """
    Compute epsilon_rel as specified in the instructions. There did explicitely stand that epsilon_rel = |fine - coarse| / |fine|
    with coarse as the currently computed integral value and fine the value from the next iteration stage.
    """
    epsilon_rel = abs(fine_value - coarse_value) / abs(fine_value)
    
    if epsilon_rel <= tolerance: # If we fulfill the tolerance, coarse is good enough
        return coarse_value
    else:
        # Otherwise divide in two intervals and compute the fine value resursively - this time recursion is enabled again.
        return adaptive_rekursive_quadrature(function, unity_weights, lower, center_point, tolerance) + adaptive_rekursive_quadrature(function, unity_weights, center_point, upper, tolerance)    
    
function = lambda x: 10/ (1+x**2) # The function 10/(1+x2)

tolerances = [1e10, 1e5, 1e1, 1, 1e-1, 1e-2, 1e-5, 1e-7, 1e-10] # Some tolerances

print("Function: 10/(1+x^2)\n")
print("Value comparison based on tolerance:")

# Generate the table header
print("+" + "-"*41 + "+")
print("| Tolerance: | Integral value in [0, 10]: |")
print("+" + "-"*41 + "+")

# Generate the table body
for tolerance in tolerances:
    # Print the computed integral values for the supplied function in [0,10]
    print("| {:.4e} | {:^26.20f} |".format(tolerance, adaptive_rekursive_quadrature(function, unity_weights(), 0, 10, tolerance)))
    print("+" + "-"*41 + "+")
