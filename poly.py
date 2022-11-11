import numpy as np
from alg import brent, dkmethod
#from sympy import Poly,apart,symbols,gcd,div
#from sympy.abc import x,y,z

# Cauchy's bound to provide upper bounds on all complex roots.
def cauchy_bound(li):
    n = len(li)
    an = li[0]
    res = []
    for i in range(1,n):
	    res.append(abs(li[i]/an))
    return 1+max(res)

# Hoelder's inequality for upper bounds on all complex roots.
def hoelder_bound(li):
    le = []
    n = len(li)
    an = li[0]
    for j in range(n):
	    le.append(li[j]**2)
    return (1/an)*sum(le)**(0.5)

# Fujiwara's bound. Better for trickier cases, when the constant term is big.
def fbound(li):
    n = len(li)
    an = li[0]
    nn = []
    for i in range(1,n):
        if i == n-1:
            nn.append(abs(li[i]/(2*an))**(1/i))
        else:
	        nn.append(abs(li[i])**(1/i)/an)
    return 2*max(nn)

def sgn(x):
    if x<0:
        return -1
    if x>0:
        return 1
    else:
        return 0

def catchitvs(fx,a,b,ep):
    f = lambda x: eval(fx)
    itvs = []
    while a<b:
        if sgn(f(a)) != sgn(f(a+ep)):
            itvs.append(a)
            itvs.append(a+ep)
        a+=ep
    return itvs

# Durand-kerner method to find the complex roots
def poly_solve(eq, tol=10e-8):
    crts = dkmethod.dk_method(eq)
    return crts

# Brent's method to find only the real roots
def poly_solve_real(eq, a, b, ep=0.2,tol=10e-8):
    rts = []
    li = catchitvs(eq,a,b,ep)
    if li == []:
        print("No solutions have been found.")
        return
    for i in range(0, len(li), 2):
        root = brent.brent(eq, li[i], li[i+1], tolerance=tol)
        rts.append(root)
    return rts

"""
def normalize(poly):
    while poly and poly[-1] == 0:
        poly.pop()
    if poly == []:
        poly.append(0)

def poly_div(num, den):
    #Create normalized copies of the args
    num = num[:]
    normalize(num)
    den = den[:]
    normalize(den)

    if len(num) >= len(den):
        #Shift den towards right so it's the same degree as num
        shiftlen = len(num) - len(den)
        den = [0] * shiftlen + den
    else:
        return [0], num

    quot = []
    divisor = float(den[-1])
    for i in range(shiftlen + 1):
        #Get the next coefficient of the quotient.
        mult = num[-1] / divisor
        quot = [mult] + quot

        #Subtract mult * den from num, but don't bother if mult == 0
        #Note that when i==0, mult!=0; so quot is automatically normalized.
        if mult != 0:
            d = [mult * u for u in den]
            num = [u - v for u, v in zip(num, d)]

        num.pop()
        den.pop(0)

    normalize(num)
    return quot, num
"""

f = 2*x**3 - 5*x**2 - 8*x + 15
g = x-3

# Division of polynomials
q, _ = div(f, g, domain='QQ')
#print(q)

# GCD of polynomials
#print(gcd(f,g,domain='QQ'))

def mult(f,g):
    st = str(Poly(f)*Poly(g))
    t = st.split(',')[0]
    return t[5:]

# Multiply polynomials
#print(mult(f,g))

# Add and subtract polynomials
#print(f+g)
#print(f-g)

w = (x**2 + 2*x + 3)/(x**3 + 4*x**2 + 5*x + 2)

# partial fraction decomposition
#print(apart(w))

#eq = "x**4+2*x**2+6*x+1"
#a = []
#a = poly_solve(eq)
#print(a)
