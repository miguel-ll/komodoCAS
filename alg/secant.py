def secant(fx, x0,x1,e=0.000001,N=20):
    f = lambda x: eval(fx)
    step = 1
    condition = True
    while condition:
        if f(x0) == f(x1):
            print('ERROR: Dividing by zero.')
            return

        x2 = x0 - (x1-x0)*f(x0)/( f(x1) - f(x0) )
        x0 = x1
        x1 = x2
        step = step + 1

        if step > N:
            print('ERROR: The value did not converge. Perhaps more iterations are necessary.')
            return

        condition = abs(f(x2)) > e
    return x2
