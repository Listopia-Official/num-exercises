"""
This is the solution of Florian Haas (3382958) for
the fourth exercise (first programming exercise) on the fourth worksheet of Numerik I.

A matrix, as specified in the instructions, will be computed for the specified dimensions.
Then we'll create the LR decomposition of this matrix step by step - every step will be visualized
and also saved at the filesystem as a video.

Finally we output the computed LR decomposition in one matrix.
"""

import numpy as np
import matplotlib.pyplot as plt

"""
Creates the matrix A as specified in the exercise instructions.

Input: n, the dimension of the quadratic matrix
Output: The computed matrix
"""
def createMatrix(n):
    mat = np.zeros(shape=(n,n))

    for i in range(n):
        for j in range(n):
            if i == 0:
                mat[i][j] = 1 # 1 on first row
            elif i == j:
                mat[i][j] = 2 # 2 at main diagonal (except at i=0)
            elif (i-1) == j or (i+1) == j: # Populate the diagonals above/below the main diagonal with -1
                mat[i][j] = -1

    return mat

"""
Creates the LR decomposition of the supplied matrix with the algorithm from the lecture.
Every step the decomposition matrix is visualized - a step is here interpreted as a operation on the matrix copy.
We use the indicies i, j, k in the ranges specified in the script, which we do for readibility purposes. For
array access, those will be shifted.

Input: The matrix which should be decomposed
Output: One matrix containing all relevant data about the LR decomposition
"""
def solve(mat):
    matrix = mat.copy() # Copy the input matrix
    
    step = 0 # Counter for every step
    
    i = 1 # Starts with 1, not with zero
    while i <= (n-1):
        i_array = i-1 # Shifted index for array access
        
        j = i+1
        while j <= n:
            j_array = j-1 # Shifted index for acces access
            
            matrix[j_array][i_array] = matrix[j_array][i_array] / matrix[i_array][i_array] # Don't visualize here, because this doesn't change the output graphic
             
            k = i+1
            while k <= n:
                k_array = k-1 # Shifted index
                
                matrix[j_array][k_array] = matrix[j_array][k_array] - matrix[j_array][i_array] * matrix[i_array][k_array] # Visualize after this step
                
                # Visualize here and increment step counter
                step+=1
                visualizeStep(matrix,step)
                
                k+=1
            j+=1
        i+=1

    return matrix

"""
Creates a graphic visualizing the supplied matrix and prints the step in the title.

Input: matrix: The matrix to visualize
       step:   The current step
"""
def visualizeStep(matrix, step):    
    plt.title("Step:"+ str(step) + "\n")
    plt.grid(linestyle='-', linewidth=1)
    plt.spy(matrix, marker=".")
    
    plt.savefig("step_"+str(step)+"_n="+str(matrix.shape[0])+".png") # Save on filesystem, with step and dimension in filename
    plt.show() # Show on console

n = 5 # Dimension of the matrix to create, change, if you want

matrix = createMatrix(n)
decomposition = solve(matrix)

print("Computed Matrix:\n",matrix)
print("\nLR decomposition matrix (in situ):\n",decomposition,"\n")
