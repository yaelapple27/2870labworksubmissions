import numpy as np
def lu_factorisation(A):
    """
    Compute the LU factorisation of a square matrix A.

    The function decomposes a square matrix ``A`` into the product of a lower
    triangular matrix ``L`` and an upper triangular matrix ``U`` such that:

    .. math::
        A = L U

    where ``L`` has unit diagonal elements and ``U`` is upper triangular.

    Parameters
    ----------
    A : numpy.ndarray
        A 2D NumPy array of shape ``(n, n)`` representing the square matrix to
        factorise.

    Returns
    -------
    L : numpy.ndarray
        A lower triangular matrix with shape ``(n, n)`` and unit diagonal.
    U : numpy.ndarray
        An upper triangular matrix with shape ``(n, n)``.
    """
    n, m = A.shape
    if n != m:
        raise ValueError(f"Matrix A is not square {A.shape=}")

    L, U = np.zeros_like(A), np.zeros_like(A)

    #diagonal line made
    for p in range(n):
        L[p,p] = 1.0

    for j in range(n):
        #computing all elements of U 
        for i in range(j+1):
            totalU = 0.0
            for k in range(i):
                totalU = totalU + (L[i,k] * U[k,j])
            U[i,j] = A[i,j] - totalU

        #computing all element of L
        #next line is j+1 to go below the diagonal line
        for i in range(j+1,n): 
            totalL = 0.0
            for k in range(j):
                totalL = totalL + (L[i,k] * U[k,j])
            L[i,j] = (A[i,j] - totalL)/U[j,j]


    return L, U

"""next question to display it"""
A = np.array([[4,2,0], 
              [2,3,1], 
              [0,1,2.5]])
L, U = lu_factorisation(A)
print("L =\n", L)
print("U =\n", U)


"""next question to see times"""
from time import perf_counter as now
#Ax = b, this means LUx = b
#this can be broken into Lz = b, Ux = z

#forwards substitution used to find z
def forwards(L, b):
    n = L.shape[0]
    z = np.zeros(n)
    for i in range(n):
        total = 0.0
        for k in range(i):
            total = total + (L[i, k] * z[k])
        z[i] = b[i] - total
    return z

#backwards substitution used to find x
def backwards(U, z):
    n = U.shape[0]
    x = np.zeros(n)
    i = n - 1
    while i >= 0:
        total = 0.0
        #above diaganol
        for k in range(i + 1, n):
            total = total + (U[i, k] * x[k])
        x[i] = (z[i] - total) / U[i, i]
        i = i - 1
    return x

def gaussian_elimination(A, b):
    # find shape of system
    n = A.shape[0]

    # perform forwards elimination
    for i in range(n - 1):
        # eliminate column i
        for j in range(i + 1, n):
            # row j
            factor = A[j, i] / A[i, i]
            

sizes = [2**j for j in range(1, 6)]

for n in sizes:
    # generate a random system of linear equations of size n
    A, b, x = generate_safe_system(n)
    b = b.flatten() 

    # do the solve

    #calculate using LU factorisation
    timer1 = now()
    L, U = lu_factorisation(A)
    timerEnd1 = now()
    solveTimeLUF = timerEnd1 - timer1
    
    #calculate using forwards an backwards substitution functions
    timer = now()
    z = forwards(L, b)
    x2 = backwards(U, z)
    timerEnd = now()
    solveTimeFB = timerEnd - timer

    #calculate using gaussian elimination
    timer2 = now()
    GE = gaussian_elimination(A, b)
    timerEnd2 = now()
    solveTimeGE = timerEnd2 - timer2

    print(f"N(sides):{n:6d}  Forwards/Backwards Substitution:{solveTimeFB:10.6f}  LU Factorisation:{solveTimeLUF:10.6f}  Gaussian Elimination:{solveTimeGE:10.6f}")