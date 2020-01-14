"""
This is the solution of Florian Haas (3382958) and Pascal Bauer (3383821) for
the fourth exercise (first programming exercise) on the tenth worksheet of Numerik I.
"""

import numpy as np
import scipy as sp

"""
Returns the monomial base polynome of the specified index (p_0 = 1, p_1 = x, p_2 = x^2, ...)
"""
def monomial(index):
    return lambda x: x**index

# Just two lists in which computed values for the base polynomes are stores, so that they don't have to be computed again
legendre_cache = [lambda x : 1, lambda x: x]
normalized_legendre_cache = []

"""
Returns the Legendre polynomial of the specified index
"""
def legendre(index):
    # If in cache, return that
    if len(legendre_cache) > index:
        return legendre_cache[index]
    else:
        # Compute it with the recursion formula from the lecture
        p_k = legendre(index-1)
        p_k1 = legendre(index-2)
        
        function =  lambda x: x * p_k(x) - (index-1)**2/(4*(index-1)**2-1) * p_k1(x)
        legendre_cache.append(function) # Store in cache
        
        return function

"""
Returns the normalized Legendre polynomial of the specified index
"""    
def legendre_norm(index):
    # If in cache, return that
    if len(normalized_legendre_cache) > index:
        return normalized_legendre_cache[index]
    else:
        # Otherwise compute the Legendre polynomial and normalize it as specified in the lecture
        function = lambda x: np.math.factorial(2*index) / (2**index * np.math.factorial(index)**2) * legendre(index)(x)
        normalized_legendre_cache.append(function) # Store in cache
        return function

"""
Computes the dot-product with the L2-Norm in [-1,1]
"""
def l2_dot_product(function):
    return sp.integrate.quadrature(function,-1,1)[0]
    
"""
Computes the gram matrix for the specified base functions
"""
def gram(base_functions):
    n = len(base_functions) # The gram matrix is a nxn matrix
    return [[l2_dot_product(lambda x: base_functions[i](x) * base_functions[j](x)) for i in range(n)] for j in range(n)] # Compute it as specified in the lecture

"""
A helper function which computes the gram matrices for the specified dimensions for the specified base functions,
then their conditional numbers for several norm functions is computes and nicely displayed in a table.

basename is printed - it's the name of the base, surprise
"""
def print_conds(dimensions, base, basename):
    # Header with general information
    print("\nResults for the "+basename+" Base: ")
    print("| n  | cond_1(A)  | cond_2(A)  | cond_inf(A) |")

    # Compute the data for every dimension
    for dimension in dimensions:
        A = gram([base(i) for i in range(dimension)])
    
        # Conditional numbers
        cond_1 = np.linalg.cond(A, 1)
        cond_2 = np.linalg.cond(A, 2)
        cond_inf = np.linalg.cond(A, np.inf)
    
        # Print the data
        print("|{:^4}| {:.4e} | {:.4e} | {:.4e}  |".format(dimension, cond_1, cond_2, cond_inf))

dimensions = [5, 10, 15, 20] # The dimensions - FOR DIMENSION 20 IT TAKES A WHILE UNTIL THE RESULTS ARE DISPLAYED

# Test cases
print_conds(dimensions, monomial, "Monomial")
print_conds(dimensions, legendre, "Legendre")
print_conds(dimensions, legendre_norm, "Normalized Legendre")

"""
Vergleich mit den Daten aus der Vorlesung:
 - Für die Monombasis stimmen die Werte mit den dort angezeigten Daten modulo Gleitkommafehler überein
 - Die Gram-Matrix ist mit der Legendre-Basis am Anfang (n=5, 10) etwas besser konditioniert als mit der Monombasis - später aber wesentlich schlechter
 - Für die normierte Legendre-Basis ist die Gram-Matrix dagegen (vergleichweise) sehr gut konditioniert - für n=20 beträgt die Konditionszahl "nur" ca. 40
"""
