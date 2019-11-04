# Demo codes to perform all numerical experiments in chapter 2 of my lecture notes 
# "Numerische Mathematik 1". 
#
# Tested on the Linux console with Python 3.7.3, NumPy and MatPlotLib from Anaconda3-2019.03,
# but any recent Python, NumPy and MatPlotLib should do the trick.
#
# Execution: python3 chap2.py
#
# (c) 2019-20, dominik.goeddeke@mathematik.uni-stuttgart.de

import numpy as np 

print(f'numpy version: {np.version.version}')
print(f'machine epsilon in double precision: {np.finfo(float).eps}') 
print(f'machine epsilon in single precision: {np.finfo(np.float32).eps}\n')

n_max = 10000

err_1 = 1e-4
err_2 = 1e0
err_3 = 1e1

first_gt_err1 = False
first_gt_err2 = False
first_gt_err3 = False

for n in range(n_max):
    A = np.zeros(shape=(n,n))

    for i in range(n):
        for j in range(n):
            if i == j or j == (n-1):
                A[i][j]=1
            elif j < i:
                A[i][j] = -1
            
    b = np.matmul(A,np.ones(shape=(n,1)))

    As = A.astype(np.float32) 
    bs = b.astype(np.float32)
    xs = np.linalg.solve(As,bs)
    
    original_x = np.ones(shape=(n,1)).astype(np.float32)

    delta_x = original_x - xs

    abs_error = np.linalg.norm(delta_x)
    
    if not first_gt_err1 and abs_error >= err_1:
        first_gt_err1 = True
        print("The absolute error (",abs_error,") is greater than err_1=",err_1," for the first time with n=",n)
        
    if not first_gt_err2 and abs_error >= err_2:
        first_gt_err2 = True
        print("The absolute error (",abs_error,") is greater than err_2=",err_2," for the first time with n=",n)
        
    if not first_gt_err3 and abs_error >= err_3:
        first_gt_err3 = True
        print("The absolute error (",abs_error,") is greater than err_3=",err_3," for the first time with n=",n)
        
    if first_gt_err1 and first_gt_err2 and first_gt_err3:
        break



