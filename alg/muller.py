from math import sqrt, floor;
 
 # a, b and c are three initial guesses of the root.
def muller(fx, a, b, c):
    MAX_ITERATIONS=100;
    f = lambda x: eval(fx)
    res = 0;
    i = 0;
    while (True):
        # Calculating various constants
        # required to calculate x3
        f1 = f(a); f2 = f(b); f3 = f(c);
        d1 = f1 - f3;
        d2 = f2 - f3;
        h1 = a - c;
        h2 = b - c;
        a0 = f3;
        a1 = (((d2 * pow(h1, 2)) -
               (d1 * pow(h2, 2))) /
              ((h1 * h2) * (h1 - h2)));
        a2 = (((d1 * h2) - (d2 * h1)) /
              ((h1 * h2) * (h1 - h2)));
        x = ((-2 * a0) / (a1 +
             abs(sqrt(a1 * a1 - 4 * a0 * a2))));
        y = ((-2 * a0) / (a1 -
            abs(sqrt(a1 * a1 - 4 * a0 * a2))));
 
        # Taking the root which is
        # closer to x2
        if (x >= y):
            res = x + c;
        else:
            res = y + c;
 
        # checking for resemblance of x3
        # with x2 till two decimal places
        m = res * 100;
        n = c * 100;
        m = floor(m);
        n = floor(n);
        if (m == n):
            break;
        a = b;
        b = c;
        c = res;
        if (i > MAX_ITERATIONS):
            print("Root couldn't have been found using Muller's method");
            break;
        i += 1;
    if (i <= MAX_ITERATIONS):
        return res;

#a = 0; b=1; c=2;
#print(muller("x**5-5*x +3",a, b, c));
