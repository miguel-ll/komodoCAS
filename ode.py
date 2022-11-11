from sympy import integrate,diff,Poly,sin,cos
from sympy.abc import x,y

def exact_solve(m,n):
    #integrate m in respect to x
    f = integrate(m,x)
    # subtract the derivative of f in respect to y from n and stores the value on g'(y).
    g = n-diff(f,y)
    # integrate g'(y) in respect to y.
    h = integrate(g,y)
    # return f + g(y), the solution to the exact diff eq.
    return str(f+h)

m = 2*x*y**4
n = 4*x**2 * y**3 - sin(y)
print( exact_solve(m,n) )
# sol. x**2 * y**4 + cos(y)
