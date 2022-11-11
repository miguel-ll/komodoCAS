from numpy import sqrt

def rmax(x,y):
	if x > y:
		return y
	else:
		return x

# f is the function, df the derivative, and ddf its second derivative.
def laguerre(f, df, ddf, x0,n):
	a = 0
	while (f(x0) != 0):
	    g = df(x0)/f(x0)
	    h = g**2 - ddf(x0)/f(x0)
	    p = g + sqrt((n-1)*(n*h - g**2))
	    q = g - sqrt((n-1)*(n*h - g**2))
	    a = n/(rmax(p,q))
	    x0 = x0 - a
	return x0
