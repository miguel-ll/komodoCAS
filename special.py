from numpy import sin, pi, exp, log,tan, sign
from math import sqrt,comb,ceil,factorial
from integrate import trapz
from nums import gam,bernoulli,harmonic
from scipy import special

p = [676.5203681218851
    ,-1259.1392167224028
    ,771.32342877765313
    ,-176.61502916214059
    ,12.507343278686905
    ,-0.13857109526572012
    ,9.9843695780195716e-6
    ,1.5056327351493116e-7]

# The gamma function, used as an extension of the factorial function to complex numbers, an approximation.
def gamma(z):
    EPSILON = 1e-10
    z = complex(z)
    if z.real < 0.5:
        y = pi / (sin(pi * z) * gamma(1 - z))  # reflection formula
    else:
        z -= 1
        x = 0.99999999999980993
        for (i, pval) in enumerate(p):
            x += pval / (z + i + 1)
        t = z + len(p) - 0.5
        y = ((2 * pi)**(0.5)) * t ** (z + 0.5) * exp(-t) * x
    if abs(y.imag) <= EPSILON:
        y = y.real
    return y

# A fairly good approximation for the log of the Gamma function.
def lgamma(z):
    cof = [76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5]
    y = x = z
    tmp = x+5.5
    tmp -= (x+0.5)*log(tmp)
    ser = 1.000000000190015
    for j in range(6):
        y+=1
        ser+=cof[j]/y
    return -tmp+log(2.5066282746310005*ser/x)

# The beta function is a function that is closely related to the Gamma function.
def beta(z1,z2):
	return (gamma(z1)*gamma(z2))/gamma(z1+z2)

# The default function for the evaluation of the Riemann Zeta function. A fast converging representation.
def rzeta(s,iter=30):
    res = 0
    for n in range(iter):
        for k in range(iter):
            res += 1/(2**(n+1))*comb(n,k)*(((-1)**k) / (k+1)**s )
    junk = 1/(1-2**(1-s))*res
#    if junk.real < 0.000000000000001 or junk.imag < 0.000000000000001:
#        return 0
    return junk

# Slowly convergent series for the Riemann Zeta function when Re(s) =/= 0. The functional equation is used for negative inputs.
def rzeta_2(s,iter=61):
    if s.real < 0:
         a = ((2*pi)**s)/pi*sin(pi*s*0.5)*gamma(1-s)*rzeta_2(1-s)
         if a.real < 0.000000000000001 or a.imag < 0.000000000000001:
             return 0
         return a
    res = 0
    for n in range(1,iter):
        res += ((-1)**(n+1))/n**s
    return res*(1/(1-2**(1-s)))

# Riemann Zeta function. Converges very slowly for Re(z) > 0. Don't even use this.
def r_zeta3(s,v=10000):
    res = 0
    for n in range(1,v):
        res += n/(n+1)**s - (n-s)/n**s
    return 1/(s-1)*res

# Dirichlet eta function. This Dirichlet series is the alternating sum corresponding to the Dirichlet series expansion of the Riemann zeta function.
# The fast convergence was acheived using an Euler transform.
def eta(s,iter=30):
    res = 0
    for n in range(iter):
        for k in range(iter):
            res += 1/(2**(n+1))*comb(n,k)*(1/(k+1)**s)*(-1)**k
    return res

# The log derivative of the Gamma function.
def digamma(z):
    # Definition in positive integers
    if z >= 0 and ceil(z)==z:
        # Undefined at zero
        if z == 0:
            return 0
        return harmonic(z-1) - gam
    # Particular values of the Digamma function
    if z == 0.5:
        return -2*log(2) - gam
    if z == 0.25:
        return -pi/2 - 3*log(2) - gam
    # Integral formula for the Digamma function
    if z > 0 and z < 4:
        return -gam - 1/z + trapz(f"(1-x**{z})/(1-x)",0,0.999999999999,40)
    # Reflection formula for the Digamma function
    if z < 0:
        return digamma(1-z) - pi/tan(pi*z)
    # Asymptotic expansion using the Bernoulli numbers
    return log(z) - 1/(2*z) - 1/(12*z**2) + 1/(120*z**4) - 1/(252*z**6) + 1/(240*z**8) - 1/(132*z**10) + 691/(32760*z**12)

# The derivative of the Digamma function. Bad when z < 3.
# TODO. Fix values for z < 3.
def trigamma(z):
    # Reflection formula for the Trigamma function
    if z < 0:
        return (pi/(sin(pi*z)))**2 - trigamma(1-z)
    # Particular values of the Trigamma function
    if z == 1:
        return pi**2/6
    if z == 0.5:
        return pi**2/2
    if z == 2:
        return pi**2/6 - 1
    if z == 1.5:
        return pi**2/2 - 4
    # Asymptotic expansion using the Bernoulli numbers of the second kind
    sm = 0
    for k in range(24):
        sm += bernoulli(k,k=2)/z**(k+1)
    return sm

def polygamma(m,z):
    sm = 0
    n = 1
    # The first derivative and second derivatives of the Gamma function are the Digamma and Trigamma, respectively.
    if m == 0:
        return digamma(z)
    if m == 1:
        return trigamma(z)
    # If m is big enough, the number of iterations in the for loop can be reduced.
    if m == 9:
        n = 10
    if m == 16:
        n = 50
    for k in range(500//n):
        sm+= 1/(z+k)**(m+1)
    return sm*factorial(m)*(-1)**(m+1)

# Error function.
def erf(x):
    sm = 0
    if x.imag == 0:
        if x >= 5.2:
            return 1.0
        for n in range(100):
            sm+= (((-1)**n)*x**(2*n+1))/(factorial(n)*(2*n+1))
        return 2/sqrt(pi) * sm
    else:
        a = [0.254829592, -0.284496736, 1.421413741, -1.453152027, 1.061405429]
        p = 0.3275911
        t = 1/(1+p*x)
        for i in range(5):
            sm+= a[i]*t**(i+1)
        return 1 - sm*exp(-x**2)

# Inverse error function. The approximation gets worse as z approaches 1.
def inverf(z):
    res = z + (pi/12)*z**3 + ((7*pi**2)/480)*z**5 + ((127*pi**3)/40320)*z**7 + ((4369*pi**4)/5806080)*z**9 + ((34807*pi**5)/182476800)*z**11
    return sqrt(pi)/2 * res

# Inverse of the complementary error function.
def inverfc(z):
	return -1*inverf(z - 1)

# Complementary error function.
def erfc(x):
    return 1-erf(x)

# Imaginary error function.
def erfi(z):
    sm=0
    for n in range(50):
        sm+= (z**(2*n+1))/(factorial(n)*(2*n+1))
    return 2/sqrt(pi)*sm

# Dawson function, where tp denotes the type. (Type 1 is the D+ Dawson function, and type 2 is D-).
def dawson(x,tp=1):
    if tp == 2:
        return sqrt(pi)/2 * exp(x**2) * erf(x)
    else:
        return sqrt(pi)/2 * exp(-x**2) * erfi(x)

# the sine integral Si(x).
def ssi(x):
    if x.real >= 25:
        return pi/2 - cos(x)/x *(1- 2/x**2 + 24/x**4 -720/x**6) - sin(x)/x * (1/x - 6/x**3 + 120/x**5 - 5040/x**7)
    sm = 0
    for n in range(40):
        sm+= ((-1)**n * x**(2*n+1))/((2*n+1)*factorial(2*n+1))
    return sm

# the cosine integral Ci(x).
def ci(x):
    if x.real >= 25:
        return sin(x)/x *(1- 2/x**2 + 24/x**4 -720/x**6) - cos(x)/x * (1/x - 6/x**3 + 120/x**5 - 5040/x**7)
    sm = 0
    for n in range(1,50):
        sm+=((-1)**n * x**(2*n))/(2*n*factorial(2*n))
    return gam + log(x) + sm

# the cosine integral cin(x).
def cin(x):
    return -1*ci(x) + gam + log(x)

# the sine integral si(x).
def si(x):
    return ssi(x) - pi/2

# the hyperbolic sine integal shi(x)
def shi(x):
    z = complex(0,1)
    r= ssi(z*x)/z
    if r.imag == 0:
        return r.real
    return r

# the hyperbolic cosine integral chi(x). The approximation is very bad for x > 4.
def chi(x):
    return gam + log(x) + x**2/4 + x**4/96 + x**6/4320 + x**8/322560 + x**10/36288000

# the logarithmic integral function
def li(x,iter=50):
    sm = 0
    for k in range(1,iter):
        sm+= (log(x)**k)/(factorial(k)*k)
    return gam + log(log(x)) + sm

# a bad approximation for the logarithmic integral function
def li2(x):
    return x/log(x) + x/log(x)**2 + (2*x)/log(x)**3 + (6*x)/log(x)**4 + (24*x)/log(x)**5 + (120*x)/log(x)**6
    #return x/log(x) + x/log(x)**2 + (2*x)/log(x)**3 + (6*x)/log(x)**4 + (24*x)/log(x)**5

# asymptotic expansion of the logarithmic integral function
def li3(x):
    z = 6
    if x > 1100:
        z = 7
    if x >= 3500:
        z = 8
    if x >= 10000:
        z = 9
    if x>= 18500:
        z = 10
    if x>= 55000:
        z = 11
    if x>=3000000:
        z = 15
    sm = 0
    for k in range(z):
        sm += factorial(k)/log(x)**k
    return x/log(x) * sm

# the fresnel integral S(x).
def fresnel_sin(x):
    sm = 0
    for n in range(35):
        sm+= (-1)**n * x**(4*n+3)/( factorial(2*n+1)*(4*n+3) )
    return sm

# the fresnel integral C(x).
def fresnel_cos(x):
    sm = 0
    for n in range(35):
        sm+= (-1)**n * x**(4*n+1)/( factorial(2*n)*(4*n+1) )
    return sm

# lower and upper incomplete gamma functions.
def lower_igamma(s,x):
    # the function is not defined for s=0.
    if s == 0:
        return 0
    sm = 0
    for k in range(60):
        sm+= ((-x)**k)/(factorial(k) * (s+k))
    return x**s * sm

def upper_igamma(s,x):
    if s==0:
        sm = 0
        for k in range(1,60):
            sm+= ((-x)**k) / (k*factorial(k))
        return -gam - log(x) - sm
    return gamma(s) - lower(s,x)

# Exponential integral E_n
def e_n(n,x):
    return x**(n-1)*upper_igamma(1-n,x)

# Exponential integral Ei
def ei(x):
    sm = 0
    for k in range(1,60):
        sm+= ((-1)**(k+1) * (-x)**k) / (factorial(k)*k)
    return gam + log(abs(x)) - sm

def W(z):
	x = 0
	for _ in range(9):
		x = x - (x*exp(x) - z)/(exp(x) + x*exp(x))
	return x

#print(ei(1.8))
