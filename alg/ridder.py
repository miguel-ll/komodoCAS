from numpy import sqrt

def ridders(y, a, b, tol=1.0**(-9)):
    f = lambda x: eval(y)
    fa = f(a)
    fb = f(b)
    if fa == 0.0:
        return a
    if fb == 0.0:
        return b
    if fa*fb > 0:
        print("No root was found in the given interval.")
        return
    for i in range(100):
        c = 0.5*(a+b)
        fc = f(c)
        s = sqrt(fc**2 - fa*fb)
        if s == 0:
            return
        dx = (c-a)*(fc/s)
        if (fa - fb) == 0.0:
            dx = -dx
        x = c+dx
        fx = f(x)
        if i > 0:
            #if (abs(x-xOld) < tol*max(abs(x),1.0)):
            if f(x) == 0:
                return x
        xOld = x
        if fc*fx == 0.0:
            if fa*fx < 0.0:
                b = x
                fb = fx
            else:
                a = x
                fa = fx
        else:
            a = c
            b = x
            fa = fc
            fb = fx
    print("Iteration limit exceeded.")
    return
