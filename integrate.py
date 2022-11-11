from numpy import *

def gquad(f,a,b):
    u=(b-a)/2.*(5./9*f((b-a)/2.*-1.*sqrt(3./5)+(b+a)/2.)+8./9*f((b+a)/2.)+5./9*f((b-a)/2.*sqrt(3./5)+(b+a)/2.))
    return u

def sub_division(f,a,b,tol,entire):
    # tol is the tolerance.
    # This function splits an integral into the right and left half and compares it to the integral on the whole function.
    # If tolerance is not satisfied, the function splits the left and right into smaller intervals and repeats the function. 
    # When the integral interval values satisfy the tolerance, they are added to the sum.
    results = []
    a_z=a+(b-a)/2.0  #sub-divides intervals
    b_k=a+(b-a)/2.0
    entire=gquad(f,a,b)
    right=gquad(f,a_z,b)
    left=gquad(f,a,b_k)
    if abs(entire-(left+right)) < tol * max(abs(entire), (abs(left)+abs(right))):
        results.append(entire)
        return entire
    x=sub_division(f,a_z,b,tol,right)+sub_division(f,a,b_k,tol,left)
    results.append(x)
    return x

def quad(fx,a,b,tol):
    f = lambda x: eval(fx)
    return sub_division(f,a,b,tol,gquad(f,a,b))

def mcarlo(fx,a,b,N):
    f = lambda x: eval(fx)
    res = 0
    xrand = zeros(N)
    for i in range(len(xrand)):
        xrand[i] = random.uniform(a,b)
    for j in range(N):
        res += f(xrand[j])
    return (b-a)/float(N)*res

def trapz(fx,a,b,n):
    f = lambda x: eval(fx)
    h = (b-a)/n
    xn = a+h
    sum2 = 0
    summ = 0.5*f(a+0.0000000000000001)+0.5*f(b+0.0000000000000001)
    for i in range(n-1):
        summ+= f(xn)
        xn += h
    z = 2*n
    for j in range(1,z):
        sum2 += f(j/z)
    sum2 += 0.5*(f(a+0.0000000000000001)+f(b+0.0000000000000001))
    sum2 *= (b-a)/(2*n)
    ehalf = (sum2 - (summ*h))/3
    return sum2+ehalf

def simps(fx,x0,xn,n):
    f = lambda x: eval(fx)
    h = (xn - x0)/n
    integration = f(x0+0.0000000000000000001) + f(xn+ 0.0000000000000000001)
    for i in range(1,n):
        k = x0 + i*h
        if i%2 == 0:
            integration = integration + 2 * f(k)
        else:
            integration = integration + 4 * f(k)
    integration = integration * h/3
    return integration

def weedle(fx, a, b):
    y = lambda x: eval(fx)
    h = (b - a) / 6;
    sum = 0;
    a+=0.00000000000000000000000000001
    sum = sum + (((3 * h) / 10) * (y(a) + y(a + 2 * h) + 5 * y(a + h) + 6 * y(a + 3 * h) + y(a + 4 * h) + 5 * y(a + 5 * h) + y(a + 6 * h)));
    return sum;

def boolerule(fx,a, b):
    y = lambda x: eval(fx)
    n = 4
    h = ((b - a) / n)
    sum = 0
    bl = (7 * y(a) + 32 * y(a + h) + 12 *
        y(a + 2 * h)+32 * y(a + 3 * h)+7 *
        y(a + 4 * h))* 2 * h / 45
    sum = sum + bl
    return sum

#f1 = "sin(x)/x"
#tol = 10e-15
#print(quad(f1,0,2, tol))
#print(mcarlo(f1,0,2,1000))
#print(trapz(f1,0,2,200))
#print(simps(f1,0,2,100))
#print(weedle(f1,0,2))
