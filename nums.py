from fractions import Fraction as Fr
from math import factorial

# Euler–Mascheroni constant
gam=0.57721566490153286060651209008240243104215933593992

# returns the nth bernoulli number (where bernoulli(1) = -0.5)
def bernoulli(n,k=1):
    if n == 1:
        if k ==1:
            return -0.5
        else:
            return 0.5
    A = [0] * (n+1)
    for m in range(n+1):
        A[m] = Fr(1, m+1)
        for j in range(m, 0, -1):
          A[j-1] = j*(A[j-1] - A[j])
    return float(A[0])

# returns the nth bernoulli number as a fraction
def bernoullifrac(n):
    if n == 1:
        return Fr(-1,2)
    A = [0] * (n+1)
    for m in range(n+1):
        A[m] = Fr(1, m+1)
        for j in range(m, 0, -1):
          A[j-1] = j*(A[j-1] - A[j])
    return A[0]

# returns the nth harmonic number
def harmonic(n):
    sm = 0
    for i in range(n):
        sm += 1/(i+1)
    return sm

# returns the nth fibonacci number
def fibonacci(n):
    sqf = 5**(0.5)
    fn = 1/sqf*(((1+sqf)*0.5)**(n) - ((1-sqf)*0.5)**(n))
    return int(fn)

# prints the pascal triangle up to the nth row
def pascal(n):
    for i in range(n):
        for j in range(n-i+1):

            print(end=" ")  # for left spacing

        for j in range(i+1):

            # nCr = n!/((n-r)!*r!)
            print(factorial(i)//(factorial(j)*factorial(i-j)), end=" ")
        print()

# get a list of the nth power of divisors of a number s
def divisors(s,n=1):
    li = []
    i = 1
    while (i * i < s):
        if (s % i == 0):
            li.append(i**n)
        i += 1
    for i in range(int(sqrt(s)), 0, -1):
        if (s % i == 0):
            li.append((s//i)**n)
    return li

# returns the nth triangular number (or the sum from k=1 to n (k) = 1+2+3...n
def triangular(n):
    return int((n*(n+1))/2)
