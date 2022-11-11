#import numpy as np
from sympy.integrals import laplace_transform
from sympy.abc import t, s, a
from sympy.integrals.transforms import inverse_laplace_transform
from sympy import exp, sin, Symbol

#a = Symbol('a', positive = True)
# Using inverse_laplace_transform() method
#gff = inverse_laplace_transform(1/(6*s + 3), s, t)
gff = str(inverse_laplace_transform(1/(s**2 + 2*s), s, t))
gff = gff.replace("Heaviside(t)","1")

#print(gff)

# Using laplace_transform() method
gfg = laplace_transform(sin(t)/t, t, s)

#print(gfg[0])
