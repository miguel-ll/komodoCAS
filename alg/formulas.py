from numpy import sin,cos,arccos

# Quadratic formula
def solvequad(a,b,c):
    x1 = (-b + (b*b-4*a*c)**(0.5))/(2*a)
    x2 = (-b - (b*b-4*a*c)**(0.5))/(2*a)
    return [x1,x2]

# Cardano's formula and trigonometric solution
def solvecubic(p,q):
    dt = -(4*p*p*p + 27*q*q)
    # if delta > 0, the cubic has three distinct real roots. if delta < 0, the cubic has one real root and two complex roots. if delta == 0, the cubic has a multiple root.
    if dt < 0 or dt == 0:
        n = (3**(0.5))/2
        unity1 = complex(-0.5, n)
        unity2 = complex(-0.5, -n)
        m = ((q*q)/4 + (p*p*p)/27)**(0.5)
        c = ((-q/2) + m)**(1/3)
        w = -p/(3*c)
        r1 = c+w
        r2 = c*unity1 + w*unity2
        r3 = c*unity2 + w*unity1
        return [r1,r2,r3]
    else:
        tk = []
        for k in range(3):
            a = ((3*q)/(2*p))*((-3/p)**(0.5))
            b = (1/3)*np.arccos(a) - (2*np.pi*k)/3
            tk.append(2*((-p/3)**(0.5))*np.cos(b))
        return tk
