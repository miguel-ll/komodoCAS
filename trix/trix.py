import numpy as np
from math import comb,factorial
from scipy import linalg
from numpy import pi, exp,sqrt, prod
from vec import Vector,proj

#def print_mat(mat):
#    for x in mat:
#        print("[",end=' ')
#        for y in x:
#            print(y, end=' ')
#        print("]")

# The input matrix must be a (m x n) matrix such that m < n. The purpose of this is to solve a linear system of equations, but it can be used to calculate determinants via the gauss_det() function (although the bareiss_det() function is preferred).
def solve_sys(mat):
    num = len(mat)
    for i in range(0, num):
        # Searching the maximum value of a particular column
        max_el = abs(mat[i][i])
        # Row having the element of maximum value
        max_row = i
        for k in range(i + 1, num):
            if abs(mat[k][i]) > max_el:
                max_el = abs(mat[k][i])
                max_row = k
        # Swapping the maximum row with the current row
        for k in range(i, num + 1):
            temp = mat[max_row][k]
            mat[max_row][k] = mat[i][k]
            mat[i][k] = temp
        # Changing the value of the rows below the current row to 0
        for k in range(i + 1, num):
            curr = -mat[k][i] / mat[i][i]
            for j in range(i, num + 1):
                if i == j:
                    mat[k][j] = 0
                else:
                    mat[k][j] += curr * mat[i][j]
    # Solving the equation Ax = b for the created upper triangular matrix mat
    l = [0 for i in range(num)]
    for j in range(num - 1, -1, -1):
        l[j] = mat[j][num] / mat[j][j]
        for k in range(j - 1, -1, -1):
            mat[k][num] -= mat[k][j] * l[j]
    return l

# Computes the determinant, using Bareiss algorithm by default.
def det(A, triangular=False,circulant=False):
    n = len(A)
    if triangular:
        res = 1
        for i in range(n):
            res*=A[i][i]
        return res
    if circulant:
        b = A[0]
        evl = []
        for j in range(n):
	        evl.append(c_eval(b,n,j))
        return prod(evl)
    divider = 1
    for i in range(n):
        p = A[i][i]
        for j in range(n):
            if(i != j):
                for k in range(i+1,n):
                    if(i==0):
                        divider = 1;
                    else:
                        divider = A[i-1][i-1];
                    A[j][k] = ((p*A[j][k]) - (A[j][i]*A[i][k])) / divider;
    return A[n-1][n-1]

def getCofactor(A, temp, p, q):
    n = len(A)
    i = 0
    j = 0
    for row in range(n):
        for col in range(n):
            # Copying into a temporary matrix
            if (row != p and col != q):
                temp[i][j] = A[row][col]
                j += 1
                # Row is filled, so increase row index and reset col index
                if (j == n - 1):
                    j = 0
                    i += 1

# Recursive determinant
def determinant(A, n):
    D = 0
    if (n == 1):
        return A[0][0]
    temp = []
    for i in range(n):
        temp.append([None for _ in range(n)])
    sign = 1
    for f in range(n):
        # getting cofactor of A[0][f]
        getCofactor(A, temp, 0, f)
        D += sign * A[0][f] * determinant(temp, n - 1)
        sign = -sign
    return D

def adjoint(A, adj):
    N = len(A)
    if (N == 1):
        adj[0][0] = 1
        return
    sign = 1
    temp = []
    for i in range(N):
        temp.append([None for _ in range(N)])
    for i in range(N):
        for j in range(N):
            getCofactor(A, temp, i, j)
            sign = int((-1)**(i+j))
            adj[j][i] = (sign)*(determinant(temp, N-1))

# Calculate the inverse of a matrix
def inverse(A):
    N = len(A)
    inverse = [None for _ in range(N)]
    for i in range(N):
        inverse[i] = [None for _ in range(N)]
    det = bareiss_det(A)
    if (det == 0):
        print("Singular matrix, can't find its inverse")
        return 0
    adj = []
    for i in range(N):
        adj.append([None for _ in range(N)])
    adjoint(A, adj)
    for i in range(N):
        for j in range(N):
            inverse[i][j] = adj[i][j] / det
    return inverse

def lu_decomp(mat):
    n = len(mat)
    lower = [[0 for x in range(n)] for y in range(n)]
    upper = [[0 for x in range(n)] for y in range(n)]
    for i in range(n):
        # Upper triangular
        for k in range(i, n):
            # Summation of L(i, j) * U(j, k)
            sm = 0
            for j in range(i):
                sm += (lower[i][j] * upper[j][k])
            # Evaluating U(i, k)
            upper[i][k] = mat[i][k] - sm
        # Lower triangular
        for k in range(i, n):
            if (i == k):
                lower[i][i] = 1
            else:
                sm = 0
                for j in range(i):
                    sm += (lower[k][j] * upper[j][i])
                lower[k][i] = (mat[k][i] - sm) / upper[i][i]
    return lower, upper

def cholesky(matrix):
    n = len(matrix)
    lower = [[0 for x in range(n)] for y in range(n)];
    # Decomposing a matrix into lower triangular
    for i in range(n):
        for j in range(i + 1):
            sum1 = 0;
            if (j == i):
                for k in range(j):
                    sum1 += pow(lower[j][k], 2);
                lower[j][j] = sqrt(matrix[j][j] - sum1);
            else:
                # Evaluating L(i, j) using L(j, j)
                for k in range(j):
                    sum1 += (lower[i][k] *lower[j][k]);
                if(lower[j][j] > 0):
                    lower[i][j] = (matrix[i][j] - sum1) / lower[j][j];
    return lower

def transpose(B):
    A = np.copy(B)
    N = len(A)
    for i in range(N):
        for j in range(i+1,N):
            A[i][j], A[j][i] = A[j][i], A[i][j]
    return A

#matrix = [[4, 12, -16], [12, 37, -43], [-16, -43, 98]];
#print(transpose(matrix))
#print(cholesky(matrix));

def mrank(A):
    return np.linalg.matrix_rank(A)

# Back substitution for upper triangular matrices
def back_subs(A,b):
    n = len(A)
    for i in range(n-1,-1,-1):
        for j in range(i+1,n):
            b[i] = b[i] - A[i][j]*b[j]
        b[i] = b[i]/A[i][i]
    return b

# Forward substitution for lower triangular matrices
def forward_subs(mat,b):
    num = len(b)
    b[0] = b[0]/mat[0][0]
    for i in range(1,num):
        for j in range(0,i):
            b[i] = b[i] - mat[i][j] * b[j]
        b[i] = b[i]/mat[i][i]
    return b

# Solve the equation ax = b, assuming a is a triangular matrix.
def solve_tri(A,b,lower=False):
    if lower:
        return forward_subs(A,b)
    else:
        return back_subs(A,b)

# Solve a linear system using lu factorization
def lu_solve(A,b):
    l, u = lu_decomp(A)
    c = forward_subs(l,b)
    d = back_subs(u,c)
    return d

# Function to calculate the nth root of unity
def wn(n):
	if n == 0:
		return 1
	ii = complex(0,1)
	return exp((2*pi*ii)/n)

# Eigenvectors for circulant matrices
def c_evec(n,j):
	res = []
	for i in range(n):
		res.append(wn(n)**(i*j))
	return res

# Eigenvalues for circulant matrices
def c_eval(b,n,j):
	res = 0
	for i in range(n):
		res+= b[i]*wn(n)**(i*j)
	return res

# Eigenvalues and eigenvectors for circular matrices. The eigenvectors are stored horizontally; i.e P = evec.
def circ_eig(A):
	n = len(A)
	b = A[0]
	evec = []
	evl = []
	for j in range(n):
	   evec.append(c_evec(n,j))
	   evl.append(c_eval(b,n,j))
	return evec, evl

# Eigenvalues of triangular matrices
def tri_eig(A):
    z = []
    n = len(A)
    for i in range(n):
        z.append(A[i][i])
    return z

# Polar decomposition. A = UP, where U is a unitary matrix and P is a positive semi-definite hermitian matrix.
def polar_decomp(A):
    A = np.array(A)
    P = linalg.sqrtm( np.matmul(A.T,A) )
    U = np.matmul(A, np.linalg.inv(P))
    return U,P


def circulant(l):
    n = len(l)
    b = np.zeros((n,n))
    for j in range(n):
        for i in range(n):
            b[i][i-j] = l[j]
    return b

def tril(a,n=0):
    m = len(a)
    for j in range(m-1,0,-1):
        for i in range(j):
            if j+n < m:
                a[i][j+n] = 0;
            else:
                break

def pascal(n):
    a = np.array([[0 for _ in range(n)] for _ in range(n)], dtype=np.float32)
    z=1
    for i in range(n-1):
        a[i][z] = z
        z+=1
    for x in range(2,n):
        a+= (1/factorial(x))*np.linalg.matrix_power(a,x)
    a += np.eye(n)
    return np.matmul(a.T,a)

def hilbert(n):
    h = np.zeros((n,n))
    for i in range(1,n+1):
        for j in range(1,n+1):
            h[i-1][j-1] = 1/(i+j-1)
    return h

def invhilbert(n):
    h = np.zeros((n,n))
    for i in range(1,n+1):
        for j in range(1,n+1):
            h[i-1][j-1] = (-1)**(i+j) * (i+j-1) * comb(n+i-1,n-j) * comb(n+j-1,n-i) * (comb(i+j-2,i-1))**2
    return h

def hdet(n):
    cn = 1
    for i in range(1,n):
        cn*= factorial(i)
    cn2 = 1
    for j in range(1,2*n):
        cn2*=factorial(j)
    return (cn**4)/cn2

def col(a,m):
    l = []
    n = len(a)
    for i in range(n):
        l.append(a[i][m])
    return l

def gram_schmidt(A):
    n = len(A)
    res = []
    for i in range(n):
        q = Vector(col(A,i))
        for j in range(i):
            if i >= 2:
                z = Vector(res[j])
                y = q
            else:
                z = Vector(col(A,j))
                y = Vector(col(A,i))
            q = q - proj(z,y)
        res.append(list(q.unit()))
    return res

def make_householder(a):
    n = len(a)
    v = a / (a[0] + np.copysign(np.linalg.norm(a), a[0]))
    v[0] = 1
    H = np.eye(n)
    H -= (2 / np.dot(v, v)) * np.dot(v[:, None], v[None, :])
    return H

def qr_decomp(A):
    n = len(A)
    Q = np.eye(n)
    for i in range(n - 1):
        H = np.eye(n)
        H[i:, i:] = make_householder(A[i:, i])
        Q = np.dot(Q, H)
        A = np.dot(H, A)
    return Q, A

def qr_decomp2(A):
    Q = np.array(gram_schmidt(A))
    R = np.matmul(Q, A)
    return Q.T, R
