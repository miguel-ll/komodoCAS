# This is not meant to be used yet. It's just a sketch of functions i pretend to use in the future if numpy ever become unnecessary for basic functions implementations.

import numpy as np

def sgn(x):
	if x > 0:
		return 1
	if x < 0:
		return -1
	else:
		return 0

def arctan2(y,x):
	if x == 0:
		return sgn(y)*np.pi/2
	if x < 0:
		return np.arctan(y/x) + sgn(y)*np.pi
	else:
		return np.arctan(y/x)

def clog(z):
	a = z.real
	b = z.imag
	c = arctan2(b,a)
	return complex(0.5*np.log(a*a + b*b), c)

def cexp(z):
    a = z.real
    b = z.imag
    ea = np.exp(a)
    m = np.sin(b)
    l = np.cos(b)
    if m > -0.000000000000001 and m < 0.000000000000001:
        m = 0
    if l > -0.000000000000001 and l < 0.000000000000001:
        l = 0
    v = complex(ea*l, ea*m)
    return v

def ccos(z):
	a = z.real
	b = z.imag
	return complex(np.cos(a)*np.cosh(b), -np.sin(a)*np.sinh(b))

def csin(z):
	a = z.real
	b = z.imag
	return complex(np.sin(a)*np.cosh(b), np.cos(a)*np.sinh(b))
	
def ctan(z):
	return csin(z)/ccos(z)
	
def csinh(z):
	a = z.real
	b = z.imag
	print(a)
	print(b)
	v = complex(np.cos(b)*np.sinh(a), np.sin(b)*np.cosh(a))
	return v

# these are lacking: complex cosh, tanh, asinh, acosh, atanh.

def carccos(z):
    i = 0 + 1j
    u = z + i*((1-(z*z))**(0.5))
    return -1*i*clog(u)

def carcsin(z):
    i = 0+1j
    u = z*i + (1-(z*z))**(0.5)
    return -1*i*clog(u)

def carctan(z):
    i = 0+1j
    return (0.5)*(i*clog((i+z)/(i-z)))

def clog2(base,z):
    return clog(z)/clog(base)
