from numpy import sin

def numdiff(fx,x,h):
    f = lambda x: eval(fx)
    y = (-f(x+2*h) + 8*f(x+h) - 8*f(x-h) + f(x-2*h))/(12*h)
    return y

#print(numdiff("sin(x)",2,0.00000001))
