"""
This is the solution of Florian Haas (3382958) and Pascal Bauer (3383821) for
the fourth exercise (first programming exercise) on the ninth worksheet of Numerik I.
"""

"""
Evaluates the composite Romberg trapezoidal quadrature

Input:
    function: The function to evaluate
    a: The lower bound
    b: The upper bound
    h: The step width
    
Return: The approximate integral value
"""
def romberg_trapezoidal(function, a, b, h):
    """
    The formula is implemented as specified in the lecture
    """
    n = (b-a) / h # from h = (b-a)/n
    
    integral = (function(a) + function(b)) / 2 # First part - the h is multiplied later
    
    for i in range(1,int(n)): # from 1 to n-1
        x_i = a + i*h # Supporting x value
        integral += function(x_i) # Evaluate the function at the supporting value
        
    integral *= h # Finally multiply h
    
    return integral
 
function = lambda x: 10/ (1+x**2) # The function 10/(1+x^2)

# The results from exercise 8 - saved here for a later comparison
comparison_data = {1e10:19.39578573241939096761,1e5:19.39578573241939096761  ,1e1:19.39578573241939096761 ,
                   1e0: 19.39578573241939096761  , 1e-1: 14.23681000488590875364 , 1e-2: 14.66139178571979861943  ,
                   1e-5:14.71129494610046428704 , 1e-7:  14.71127671063525532702  , 1e-10: 14.71127674323924239275  }

print("Function: 10/(1+x^2)\n")

"""
The maximum i (exclusive) for the sequence of h_i := 1 / 2^i.
i starts from 0
The default value for h is 1 (thus h/2^i is 1/2^i).
"""
max_support_iteration = 8 # You can adjust that

print("Resursive quadrature:")
print("+" + "-"*41 + "+")
print("| Tolerance: | Integral value in [0, 10]: |")
print("+" + "-"*41 + "+")

# Generate the table body
for entry in comparison_data.items():
    # Print the computed integral values for the supplied function in [0,10]
    print("| {:.4e} | {:^26.20f} |".format(entry[0], entry[1]))
    print("+" + "-"*41 + "+")

print("\nTrapezoidal Romberg integration:")
print("+" + "-"*41 + "+")
print("|     h:     | Integral value in [0, 10]: |")
print("+" + "-"*41 + "+")

# Compute the sequence of supporting values for the specified count of iterations
for i in range(max_support_iteration):
    h_i = 2**(-i) # h is 1
    integral = romberg_trapezoidal(function, 0, 10, h_i)
    print("| {:.4e} | {:^26.20f} |".format(h_i, integral))
    print("+" + "-"*41 + "+")
    
print("\nComparison: Look at the two tables and interpret it in any other meaningful way you can imagine.")
print("We can see that very small values for h are needed so that for example we are axact for the first three digits - with the resursive one a tolerance of ca. 1e-5 (eventually larger) is needed so that we do achieve that too. One can generally see that the values do converge really fast for Romberg.")
    
