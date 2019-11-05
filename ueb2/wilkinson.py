"""
This is the solution of Florian Haas (3382958) and Pascal Bauer (3383821) for
the fourth exercise (first programming exercise) on the second worksheet of Numerik I.

This program computes the Wilkinson Matrix for every dimension n in the range of 0 and the maximum specified one,
and solves the specified linear equation system with the right side defined as in the exercise description.

Then the absolute error between the computed and expected solution vector is computed, it'll be compared with three
pre-defined error values, and the first dimension when the error exceeds a certain of one of those three thresholds
is printed, until the highes threshold was reached; then the program terminates.

We compute in single precision.
"""
import numpy as np 

n_max = 10000 # Maximum matrix evaluated

# The specified error values
err_1 = 1e-4
err_2 = 1e0
err_3 = 1e1

# Flags which state whether the assigned error value was reached or not
first_gt_err1 = False
first_gt_err2 = False
first_gt_err3 = False

# Dimensional iteration
for n in range(n_max):
    A = np.zeros(shape=(n,n)) # Matrix in right dimensions filled with zeros

    # Compute Wilkinson matrix
    for i in range(n):
        for j in range(n):
            if i == j or j == (n-1):
                A[i][j]=1
            elif j < i:
                A[i][j] = -1
            
    b = np.matmul(A,np.ones(shape=(n,1))) # Compute b as specified

    As = A.astype(np.float32)  # Ensure single precision
    bs = b.astype(np.float32) # Ensure single precision
    xs = np.linalg.solve(As,bs) # Solve linear equation system
    
    original_x = np.ones(shape=(n,1)).astype(np.float32) # The "real" solution, base of comparison

    delta_x = original_x - xs # Absolute error vector

    abs_error = np.linalg.norm(delta_x) # Absolute error scalar value
    
    """ 
    Compare the error with the pre-defines ones, notify the user when the threshold was reached for the first time.
    Yes, the three blocks are redundant, but you know, students are lazy...
    """
    
    if not first_gt_err1 and abs_error >= err_1:
        first_gt_err1 = True # Don't check this again
        print("The absolute error (",abs_error,") is greater than err_1=",err_1," for the first time with n=",n)
        
    if not first_gt_err2 and abs_error >= err_2:
        first_gt_err2 = True
        print("The absolute error (",abs_error,") is greater than err_2=",err_2," for the first time with n=",n)
        
    if not first_gt_err3 and abs_error >= err_3:
        first_gt_err3 = True
        print("The absolute error (",abs_error,") is greater than err_3=",err_3," for the first time with n=",n)
        
    # Terminate if all tresholds were reached
    if first_gt_err1 and first_gt_err2 and first_gt_err3:
        break
