def steff(eq,x0,tol=10e-7):
    h = 10e-6
    f = lambda x: eval(eq)
    g = lambda y: (f(y+h)-f(y-h))/(2*h)
    xn = x0
    for i in range(40):
        xn = xn - f(xn)/g(xn)
        if f(xn) > -tol and f(xn) < tol:
            return xn
    return xn
