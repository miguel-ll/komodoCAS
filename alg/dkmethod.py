import random
import re

def prod_diff(value, roots):
    prod = 1
    for i in roots:
        if i != value:
            prod *= (value-i)
    return prod

def gen_cnums(num):
    vals = []
    for i in range(num):
        a = complex(random.random(),random.random())
        vals.append(a)
    return vals

def dk_method(fx, iterations=20):
    r = re.search("\*\*\d+", fx)
    degree = int(r.group()[2])
    func = lambda x: eval(fx)
    vals = gen_cnums(degree)
    cpy = list(vals)
    for _ in range(iterations):
        for i in range(len(vals)):
            cpy[i] = cpy[i] - func(cpy[i])/prod_diff(cpy[i],cpy)
    return cpy

def dk_method2(fx, iterations=20):
    r = re.search("\*\*\d+", fx)
    degree = int(r.group()[2])
    func = lambda x: eval(fx)
    res = []
    vals = gen_cnums(degree)
    cpy = list(vals)
    for _ in range(iterations):
        for i in range(len(vals)):
            cpy[i] = cpy[i] - func(cpy[i])/prod_diff(cpy[i],cpy)
    for j in range(len(cpy)):
        if not (cpy[j].imag > -0.0000000001 and cpy[j].imag < 0.0000000001):
            res.append(cpy[j])
    return res
